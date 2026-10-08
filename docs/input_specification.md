# Input arrays

Inputs stay on your computer. The command-line interface reads NumPy `.npz` files
with `allow_pickle=False`; it does not read image headers or copy source metadata.
Do not commit your input or output arrays to this repository.

## Calibration

Each file contains:

| Key | Shape | Meaning |
| --- | --- | --- |
| `b0` | X × Y × Z | Mean b=0 signal |
| `d1000` | X × Y × Z | Mean b≈1000 s/mm² signal |
| `fat` | X × Y × Z | Dixon fat signal on the DWI grid |
| `body` | X × Y × Z | Boolean body mask |
| `ff` | X × Y × Z | Dixon fat fraction, 0–1 within the body |
| `scale` | scalar | Product of water-to-DWI scaling and water/fat SPGR ratio |

Use aligned, finite arrays. Zero-fill outside the field of view; define the body
mask only where valid Dixon and DWI coverage overlap. Array axis 2 is the slice
axis. Set PE to axis 0 or 1 explicitly; these are array indices, not anatomical
directions. Reorientation changes the sign convention.

```sh
dixon-fat calibrate --input calibration_a.npz --input calibration_b.npz --pe-axis 0 --slices 2 3 4 5 --output kernel.json
```

Repeat `--input` for additional independent calibration cases. No case names or
individual results are written to the kernel file. At least 500 valid voxels per
selected slice are required, matching the study's calibration rule. This is not a
minimum calibration sample-size recommendation.

## Correction

A correction file contains `dwi` (X × Y × Z × volumes), `fat` (X × Y × Z), and
`scale` (scalar). Use the same preparation and intensity convention as calibration.

```sh
dixon-fat correct --input evaluation.npz --kernel kernel.json --output corrected.npz
```

The output contains the corrected signal array only. Every volume, including b=0,
receives the same subtraction. Negative signals are retained at this stage. The
study's fitting routine applies its signal floor and b=0 normalization afterward.
The output is still individual imaging data and must remain local.

## Preparation used in the study

Dixon water/fat images were resampled to the DWI grid with trilinear interpolation
and 27 subvoxel samples. Slice-wise in-plane translations were chosen by water/b=0
correlation within ±8 pixels. The primary analysis used slice indices 2–17.
Water-to-DWI scaling was the median b=0 / aligned-water ratio in lean voxels
(FF≤0.12). The additional SPGR ratio was water response / fat response, using
T1 values of 1420 and 370 ms and the Dixon flip angle and TR.

Registration, water/fat assignment, masking, and intensity scaling remain explicit
preparation steps. They are not inferred from filenames or scanner metadata by
this package. The Python helpers in `dixon_dwi.preparation` provide the numerical
resampling and scaling operations; users supply their own arrays and transforms.
