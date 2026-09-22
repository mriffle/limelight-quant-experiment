---
schema_version: 1
generated: 2026-09-22
---

# Script registry

| Task | Path | Kind | Status | Provides | Uses | Seeded from | Description |
|------|------|------|--------|----------|------|-------------|-------------|
| metadata-characterize | scripts/promoted/metadata_characterize.py | analysis | promoted | — | common.hashing, common.design | batch-correct-combat@0.2 | Stage-1 metadata validity checks, design-invariant hypothesis tests (batch/run-order/pairing), confounding statistics (Cramer's V raw + bias-corrected), and cohort-characterization tables for the raloxifene/HLM LFQ metadata. Reads header lines only; makes no figures. |
| metadata-figures | scripts/promoted/metadata_figures.py | analysis | promoted | — | common.hashing, common.design, common.figures.colors, common.figures.figure_io | — | Stage-1 cohort figures (level counts, condition x batch, condition x run half, within-batch run layout) rendered from the precomputed results/metadata tables; registry colors, dual export + separate legend, per-figure provenance JSON. Refuses to render if the upstream run failed (FAILED.json present). |
| — | scripts/promoted/common/hashing.py | module | promoted | sha256_of_file, compute_file_hashes, compute_data_version, duplicates | — | — | sha256 file hashing + role-keyed data-version stamp, shared by metadata_characterize.py and metadata_figures.py so there is exactly one hashing implementation in the project. |
| — | scripts/promoted/common/design.py | module | promoted | run_half_label, FAILURE_MARKER_NAME | — | — | The single run-half (ceil(n/2), middle sample -> early) rule and the upstream-failure-marker filename shared by metadata_characterize.py and metadata_figures.py. |
| — | scripts/promoted/common/figures/colors.py | module | promoted | BACKGROUND_COLOR, Palette, CategoricalPaletteExceededError, DEFAULT_REGISTRY_PATH, load_registry, load_palette, assign_colors, get_color | — | okabe-ito-colors@0.2 | Project color registry: reads/extends state/color_registry.json so a given (category, value) gets the same Okabe-Ito color in every figure. Deterministic next-unused-color assignment, persisted; gray background labels outside the palette; raises CategoricalPaletteExceededError past 8 categories (the >8-category guard). Study-agnostic; fail-loud. |
| — | scripts/promoted/common/figures/figure_io.py | module | promoted | PUBLICATION_RCPARAMS, FigureArtifacts, publication_style, save_figure | — | figure-io@0.3 | Figure save helpers enforcing the visualization conventions: dual export (SVG vector + 300-DPI PNG) plus an optional companion legend figure exported as <base>.legend.{svg,png} (kept out of the plot so it cannot overlap the data), and a shared publication matplotlib style. Deviation from figure-io@0.3: SVG output is byte-deterministic (fixed svg.hashsalt, no Date metadata), and save_figure must be called inside publication_style(). Study-agnostic; fail-loud. |
