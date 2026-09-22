# figures/

Structured layout (see `conventions/visualization.md`, *Where figures live*). Every figure here is a regenerable artifact of a script — never a hand-made image.

- `metadata/<family>/` — Stage-1 cohort / metadata figures.
- `qc/<family>/` — Stage-3 QC figures.
- `analysis/<family>/<label>/` — Stage-4+ results figures.

Rules:
- `<family>` is the plot/analysis template name (e.g. `volcano`, `pca`, `missingness`).
- Processing state (e.g. `log2`, `median-norm`, `imputed`) goes in the **file stem**, never as a directory.
- A figure attached to a finding has a stem that **starts with the finding id** (e.g. `F0007_volcano_log2-median-norm.svg`).
- The tree never goes deeper than three levels below `figures/`.
