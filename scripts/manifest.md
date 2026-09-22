---
schema_version: 1
generated: 2026-09-22
---

# Script registry

| Task | Path | Kind | Status | Provides | Uses | Seeded from | Description |
|------|------|------|--------|----------|------|-------------|-------------|
| metadata-characterize | scripts/promoted/metadata_characterize.py | analysis | scratch | — | — | batch-correct-combat@0.2 | Stage-1 metadata validity checks, design-invariant hypothesis tests (batch/run-order/pairing), confounding statistics (Cramer's V raw + bias-corrected), and cohort-characterization tables for the raloxifene/HLM LFQ metadata. Reads header lines only; makes no figures. |
