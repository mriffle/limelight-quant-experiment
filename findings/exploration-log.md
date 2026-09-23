# Exploration log

Append-only record of what was looked at and discarded — the multiplicity context that informs each finding's caveats (doc 03.6). One dated entry per exploratory thread.

## 2026-09-23 — Differential abundance, raloxifene-d0 vs control (Stage 4, first pass)

- **What was run:** a moderated (limma-style) linear model, BH-corrected per quantity × design. Runner `scripts/scratch/de_raloxifene_vs_control.py`; outputs in `results/de/raloxifene-vs-control/`.
- **Quantities (4):**
  - LFQ protein (1,801) and LFQ peptide (10,229), both median-normalized log2;
  - NSAF (log2, as-is; 2,158 after dropping 10 constant);
  - log2 PSM counts (unnormalized; 2,167 after dropping 1 constant).
- **Designs (3):**
  - paired, condition + candidate_pair: primary, chosen by the scientist from QC (finding 0003);
  - condition + batch;
  - condition only.
- **That makes 12 analyses in total, plus a common-set run** (1,640 groups) comparing LFQ, NSAF and PSM head to head.
- **Unrequested diagnostic:** the statistician added a within-pair relabelling diagnostic (8 labellings).
- **Results:**
  - Nothing reaches q < 0.05 anywhere.
  - The only q < 0.10 results are 10 LFQ peptides in the paired design (q 0.052–0.077), and none of their proteins is a protein-level hit.
  - The batch and unadjusted designs are conservative (π0 = 1.00).
- **Post-hoc look, not a claim:** CYP3A4 protein log2FC −0.24 [−0.53, +0.06], p 0.096. This is context only (raloxifene inactivates CYP3A4) and was not pre-specified.
- **Discarded / not pursued:** none yet. The designs other than paired are kept as sensitivity analyses.
