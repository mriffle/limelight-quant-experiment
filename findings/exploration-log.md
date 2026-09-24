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

## 2026-09-24 — limma-trend sensitivity + ATPK case study

- **Model choice was post hoc.** limma-trend was chosen *after* ATPK was seen near the top of the no-trend ranking. That adds multiplicity on top of the 12 quantity × design fits.
- **What was run:** trend=TRUE on all 4 quantities × 3 designs (`scripts/scratch/de_trend_sensitivity.py`), matching R limma 3.58.1 to about 1e-12.
- **Results:**
  - LFQ protein, paired: 1 hit (ATPK, q 0.040).
  - LFQ peptide, paired: 10 hits at q < 0.05, the same 10 peptides as no-trend q < 0.10.
  - Every other quantity × design: 0 hits.
- **Within-pair relabelling:** the observed labelling ranks first, but the floor is p = 1/8, and it is the only labelling aligned with run order.
- **ATPK case study** (`scripts/scratch/atpk_case_study.py`):
  - The protein was examined because it was the top hit, so its effect size is selection-optimistic.
  - The LFQ-vs-spectral disagreement comes from a single Met-oxidized peptide whose PSM count and MS1 intensity move in opposite directions.
- **Running in parallel:** the abundance-agreement analysis (LFQ vs NSAF / PSM) and the peptide/modification analysis (raloxifene adducts, presence/absence, modification classes).

## 2026-09-24 — Peptide/modification analysis and LFQ-vs-spectral abundance agreement

- **Peptide/modification analysis** (`scripts/scratch/peptide_mod_analysis.py`). Mass decomposition, the adduct class, presence/absence under 3 detection definitions × 5 classes, intensity by class, protein-adjusted peptide changes, modified-minus-reference forms, and the oxidation index.
  - **Stats review:** passed after wording fixes.
    - The detection asymmetry must not be attributed to drift.
    - The oxidation shift must be reported together with its non-significant related tests: protein-adjusted p = 0.072, global index p = 0.55/0.10, 0/64 proteins.
  - **Multiplicity:** many pair-level class tests. BH was applied only within small families; a BH over all of them is still owed.
- **Abundance-agreement analysis** (`scripts/scratch/abundance_agreement.py`). Correlations, slopes, the count floor, residual drivers, detection, and within-protein tracking.
  - **Stats review:** failed on interpretation only.
    - The "~13% compression" slope claim was **withdrawn**: the slope is not identified, and it is 1.33 on a per-peptide LFQ scale.
    - The "68 proteins lose LFQ to contaminant copies" count was **corrected to 18** (11 share exclusively with contaminant copies). 68 is only "share ≥1 peptide".
    - LFQ − log2 L vs NSAF agreement was **dropped**: it is circular.
    - Per-protein permutation p-values over 8 runs assume run exchangeability, which the pairs violate.
- **Scientist decision:** handle contaminant peptide sharing as a caveat now, with the fix made upstream later (a contaminant FASTA without human proteins).
