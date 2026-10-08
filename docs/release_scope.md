# Release contents

This release contains software, documentation, and selected cohort-level outputs.
There is no dataset. The example creates artificial signals at runtime, without
participant anatomy or acquisition files.

`protocol_kernel.json` contains the primary calibration kernel: integer shifts,
weights, total mass, centroid, and the study array PE axis. It is protocol-specific.
The CSV files contain one row per method/metric, with cohort estimates, pointwise
95% bootstrap intervals, and counts. MD is in µm²/ms; slopes are per fat-fraction
percentage point. Boundary asymmetry is receiving minus mirror muscle. These
summaries preserve the recorded results; they cannot regenerate individual points.

`paper/figures` contains only the method schematic and aggregate kernel figure
from the manuscript, as editable TikZ sources. No image or participant scatter
figure is included. The kernel figure also shows pooled phantom-calibration
kernels; it contains no anatomical geometry from those phantoms.

Boundary and interior summaries use the common B3 comparison run for B0/B2/B3/M
and the B4 run for B4. Both runs reproduced the primary methods. No raw result
JSON, per-case table, sufficient-statistic matrix, bootstrap sample, image, mask,
cache, individual coefficient, or patient-derived phantom is distributed.

The study used its own preparation and cohort orchestration. This package provides
the algorithm on explicitly prepared arrays. It does not reproduce the complete
cohort workflow without those private inputs.

Before making a release, run `python tools/check_release.py`, review the complete
file list and contents, and archive only the allowlisted files. The scanner is a
backstop; it does not anonymize files or certify arbitrary future additions.
