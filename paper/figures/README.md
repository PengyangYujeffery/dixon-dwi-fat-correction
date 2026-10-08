# Figure sources

These are the manuscript's method schematic and aggregate kernel figure.
All elements are editable TikZ or PGFPlots. No image files are needed.

With a TeX installation containing `standalone`, `tikz`, and `pgfplots`, run from
the repository root after creating an `outputs` folder:

```sh
pdflatex -output-directory=outputs paper/figures/method.tex
pdflatex -output-directory=outputs paper/figures/kernels.tex
```

Rendered files stay in `outputs` and are excluded from the release.
