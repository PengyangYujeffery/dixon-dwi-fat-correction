# Dixon-guided fat correction for diffusion MRI

Residual fat can appear away from its anatomical location in diffusion MRI. This
package uses a Dixon fat image and a calibrated displacement kernel to estimate
that signal before diffusion-model fitting.

It accompanies *Displaced Fat Signal in Thigh Diffusion MRI: A Dixon-Guided
Correction*, prepared for Investigative Radiology.

## Try it

Use Python 3.11 or newer. From this folder:

```sh
python -m venv .venv
```

Activate the environment (`.venv\Scripts\activate` on Windows, or
`source .venv/bin/activate` on Linux/macOS), then run:

```sh
python -m pip install -r requirements.txt
python -m pip install --no-deps -e .
python examples/synthetic_demo.py
python -m unittest discover -s tests
```

The example creates artificial signals in memory and reports correction errors.
It reads no images, downloads nothing, and saves no dataset. It checks the software
workflow; it is not an additional study experiment.

## What is included

- Kernel calibration and signal correction.
- The study's diffusion-model fitting and physics/offset comparators.
- Selected cohort-level summaries and the protocol-specific fitted kernel.
- Editable TikZ sources for the method schematic and aggregate kernel figure.

No participant-level data, medical images, anatomical phantoms, or study logs are
included. The cohort analyses cannot be recomputed without the original data.

For your own prepared arrays, see [the input guide](docs/input_specification.md).
The package does not convert DICOM or automatically register Dixon images.

## Scope

The correction was studied with one acquisition family. Recalibrate for a new
protocol, keep calibration and evaluation data separate, and verify the PE axis
and displacement sign. The supplied kernel is a study reference, not a universal
correction. Reduced boundary asymmetry does not establish unbiased microstructure
estimates or improved repeatability.

See [method details](docs/method.md) and [release contents](docs/release_scope.md).

## License and citation

Code: [BSD-3-Clause](LICENSE), copyright Pengyang Yu. This software license does
not grant access to, or rights over, the collaborators' study data.

Version 1.0.0 is archived on Zenodo: [10.5281/zenodo.23242353](https://doi.org/10.5281/zenodo.23242353).
For the software citation, see [CITATION.cff](CITATION.cff).
