# Method

Let F be the aligned Dixon fat image and s its conversion to DWI signal units.
The residual estimate is s times the sum of kappa[k] × F shifted by k PE pixels.
Shifts use zero padding, never circular wrapping. A positive shift moves signal
toward increasing array indices. The estimate is subtracted from every DWI volume.

## Calibration

Interior lean voxels have FF≤0.12 and no FF≥0.80 voxel within 12 PE pixels in either
direction. Their median b≈1000 / b=0 ratio estimates water attenuation, a. The
non-attenuating component is (D1000 − a × D0) / (1 − a).

Calibration regresses this component, divided by s, on 81 shifted fat images
(−40…+40 pixels) and an intercept. Kernel coefficients are nonnegative; the
intercept is unconstrained. Normal equations are pooled across calibration cases.
The solver retains the study's ridge term, 0.001 times the mean design diagonal,
and its split-variable representation of the intercept. The intercept is not
included in the subtracted fat estimate.

Do not fit the kernel on evaluation cases. Minimum calibration sample size and
cross-scanner transfer were not established. Kernel mass is the sum of weights;
its centroid is their weighted mean displacement in reconstructed pixels.

## Comparators and fitting

B0 leaves signal unchanged. B2 uses the same kernel mass at zero displacement.
B3 uses olefinic and glycerol peaks with linear interpolation onto integer shifts.
Its polarity and mass were calibration-derived; it is not calibration-free.
`physics_kernel` requires water–fat shift, acquisition-to-reconstruction scaling,
polarity, and mass as explicit arguments.

`Protocol` accepts b-values in s/mm² and b-vectors with shape 3 × volumes. It retains
the study's weighted tensor fit and dictionary-based NODDI implementation,
including 6000 integration directions, diffusivities 1.7/3.0 µm²/ms, and the
original parameter grids. `fit_methods` takes voxel × volume signals and named
per-voxel fat estimates (`None` for B0). It floors corrected signals to 0.001 times
the absolute corrected b=0 mean before normalization.

The B4 comparator adds a nonnegative constant column to the same dictionary.
That offset is excluded from tissue fractions and removed before the tensor fit.
It was a post hoc comparator in the study. See `comparators.fit_offset`.

These modules contain the study's numerical fitting operations, not a replacement
for independent clinical or acquisition validation. Directional boundary metrics,
cohort selection, and participant-derived phantom experiments are outside this
release's runnable workflow.
