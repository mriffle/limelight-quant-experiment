# QC report — raloxifene vs control HLM, FlashLFQ LFQ

*Stage 3 (integrity gate) — **signed off 2026-09-23 by Michael Riffle.** Data version `sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74`. These are descriptive sanity checks, not findings. Every figure is regenerable from `scripts/promoted/qc_fig_*.py` and the processing states in `results/qc_states/` (built by `scripts/promoted/qc_prep.py`).*

## 1. Data read — verified

Loaders: `scripts/promoted/loaders/`. Verification: `scripts/promoted/verify_loaders.py` → `results/stage3/verification.md`. **42/42 checks pass:**

- **Counts reconcile with the raw files:** 4,344 proteins and 33,358 peptides, each × 8 samples. Orientation is confirmed: features are rows, runs are columns.
- **Spot reconciliation:** 250 random cells per matrix, re-read through an independent `csv`-module parse, 0 mismatches.
- **Two-way derivation:** per-sample totals and medians agree between the loaders and the independent parse.
- **Identifiers round-trip exactly as text.** Dtypes are explicit, and all values are positive when present.
- **Missing values:**
  - Protein: 13,084 `0` + 14 literal `NaN` → missing, with the two tokens tracked separately.
  - Peptide: 89,181 `0`.
  - The Detection-Type invariant (positive ⇔ MSMS/MBR) holds for every cell.
- **Contaminants (revised rule):** 33 protein entries (no-accession *and* not grouped with a real accession) and 449 peptides. Another 19 no-accession rows stand in for real proteins in indistinguishable groups (e.g. `sp|CATD_HUMAN|,sp|P07339|CATD_HUMAN`) and are kept, flagged `contaminant_grouped_with_real`.
- **Complete non-contaminant analysis set:** 1,801 proteins and 10,229 peptides.
- **Sample ↔ metadata pairing is exact:** 8 ↔ 8, no orphans or duplicates. Condition 4/4, batch 6/2. **All 8 samples are experimental; there are no control or QC samples.**
- **Limelight dump labels are resolved by value correspondence in code.** Each label matches exactly one run (100% vs ≤ 44% for the runner-up) and reproduces `results/stage2/limelight_label_map.tsv`.
- **PSM counts and NSAF are fully verified:** all 34,752 cells of each match the raw dump. Dump ↔ quant rows are 1:1, one per protein group.
- **Pairing re-derived independently** from `metadata.tsv` + the mapping file and checked row by row.
- **Data version:** the canonical six-file stamp, reproduced by every script.
- **Scale:** `linear`, recorded on the loaded `Dataset`.
- **Loader tests:** fixtures, properties, planted truth (a 2× effect recovered after median normalization) and edge cases.

Processing states (`results/qc_states/<level>/<state>`, `qc_prep.py`):
- **Built on the complete non-contaminant set:** raw log2, median-normalized (linear and log2), and a batch-corrected **preview**. The preview is ComBat given the **batch label only**; it ran without warnings.
- **Batch confounding (`assess_batch_confounding`):**
  - batch × condition: Cramér's V = 0.00, not confounded.
  - batch × candidate pair: V = 0.82. This is expected by construction, since each pair lies within one batch.

## 2. Identification depth

![ID depth](../figures/qc/id-depth/id-depth-raw-linear.png)
![ID depth legend](../figures/qc/id-depth/id-depth-raw-linear.legend.png)

Detected features per run (raw, linear), in acquisition order, colored by condition.
- **The six 2021 runs:** 2,682–2,972 protein groups and 21.7k–25.2k peptides.
- **Both 2022 runs (AZ941, AZ942):** about 20–25% lower, at about 2,250 protein groups and 17.7k peptides.
- The drop follows **batch, not condition.** AZ908 is the shallowest 2021 run.

## 3. Missingness

![Missingness, protein](../figures/qc/missingness/missingness-protein-raw-linear.png)
![Missingness, protein legend](../figures/qc/missingness/missingness-protein-raw-linear.legend.png)
![Missingness, peptide](../figures/qc/missingness/missingness-peptide-raw-linear.png)
![Missingness, peptide legend](../figures/qc/missingness/missingness-peptide-raw-linear.legend.png)

- **Completeness (left panels):**
  - Requiring detection in all 8 runs keeps 1,805 of the 3,423 detected protein groups and 10,392 of the 31,635 detected peptides.
  - The 2022 batch curve sits well below the 2021 curve.
- **Missingness depends on abundance (right panels):** detection rate rises with abundance (Pearson r = 0.68 for proteins, 0.51 for peptides). That is the pattern of left-censoring, where low-abundance features drop out, i.e. missing not at random.
- **Bearing on your Stage-2 missing-value decision (complete features only):** the complete-feature set is biased toward abundant proteins. Any imputation used later should be left-censored, not KNN or mean imputation.

## 4. Dynamic range

![Dynamic range](../figures/qc/dynamic-range/dynamic-range-protein-raw-linear.png)
![Dynamic range legend](../figures/qc/dynamic-range/dynamic-range-protein-raw-linear.legend.png)

- **Range:** protein median intensities span about 4.7 orders of magnitude, and 3.9 orders among proteins detected in every run. The head of the curve rolls off smoothly, with no single dominant protein.
- **Top of the ranking:** ER and liver proteins (BIP, PDIA1/3/4, ACSL1, NCPR, ATPB, …) plus some serum proteins (HPT, CO3), consistent with microsomes.
- **Contaminants (orange; 33 pure entries, 9 of them detected)** sit mid-curve or lower.
- Real proteins that share a group with a contaminant-list entry (e.g. ALBU, CATD, CRP) are now treated as ordinary proteins.

## 5. Abundance distributions — raw → normalized → batch-corrected

![Abundance boxplots, protein](../figures/qc/abundance-boxplot/abundance-boxplot-protein-raw-normalized-batchcorrected-log2.png)
![Abundance boxplots, protein legend](../figures/qc/abundance-boxplot/abundance-boxplot-protein-raw-normalized-batchcorrected-log2.legend.png)
![Abundance boxplots, peptide](../figures/qc/abundance-boxplot/abundance-boxplot-peptide-raw-normalized-batchcorrected-log2.png)
![Abundance boxplots, peptide legend](../figures/qc/abundance-boxplot/abundance-boxplot-peptide-raw-normalized-batchcorrected-log2.legend.png)

- **Raw:** per-sample medians spread over 1.17 log2 for proteins (about 2.3-fold) and 1.34 log2 for peptides. AZ905 and AZ906 are highest and AZ942 lowest.
- **Median normalization** aligns them (range 0.00).
- **The ComBat preview** leaves them essentially flat (0.03 / 0.04).
- **Box widths are similar across samples,** so median scaling alone is enough to align the distributions.

## 6. CV by processing state

![CV, protein](../figures/qc/cv/cv-experimental-protein-raw-normalized-batchcorrected-linear.png)
![CV, protein legend](../figures/qc/cv/cv-experimental-protein-raw-normalized-batchcorrected-linear.legend.png)
![CV, peptide](../figures/qc/cv/cv-experimental-peptide-raw-normalized-batchcorrected-linear.png)
![CV, peptide legend](../figures/qc/cv/cv-experimental-peptide-raw-normalized-batchcorrected-linear.legend.png)

| Level | Raw | Median-normalized | ComBat preview (de-logged) |
|---|---|---|---|
| Protein (1,801) | 0.397 | 0.268 | 0.171 |
| Peptide (10,229) | 0.452 | 0.343 | 0.209 |

- Normalization removes a large share of the between-run variation.
- These CVs are computed across all 8 samples, both conditions, so they combine biological and technical variation. **No pooled-QC runs exist, so technical CV can't be separated out.**
- **The ComBat drop is partly mechanical:** it sets a 2-sample batch's mean to the grand mean. Read it as "batch is a real variance component", not as improved quantification.

## 7. Sample correlation

![Sample correlation, protein](../figures/qc/sample-correlation/sample-correlation-protein-normalized-log2.png)
![Sample correlation, protein legend](../figures/qc/sample-correlation/sample-correlation-protein-normalized-log2.legend.png)
![Sample correlation, peptide](../figures/qc/sample-correlation/sample-correlation-peptide-normalized-log2.png)
![Sample correlation, peptide legend](../figures/qc/sample-correlation/sample-correlation-peptide-normalized-log2.legend.png)

**The main structural result of QC.** On the median-normalized log2 data:
- **Every candidate pair clusters first:**
  - At the protein level each control sits next to its raloxifene partner: AZ905–AZ906 r = 0.992, AZ909–AZ910 0.992, AZ907–AZ908 0.987, AZ941–AZ942 0.986.
  - At the peptide level the pair correlations are 0.96–0.98.
- **Batch** is the top-level split: cross-batch r is 0.90–0.945 for proteins.
- **Condition plays no part in the clustering.**

## 8. PCA

![PCA by batch, protein](../figures/qc/pca/pca-by-batch-protein-raw-normalized-batchcorrected-log2.png)
![PCA by batch, protein legend](../figures/qc/pca/pca-by-batch-protein-raw-normalized-batchcorrected-log2.legend.png)
![PCA by condition, protein](../figures/qc/pca/pca-by-condition-protein-raw-normalized-batchcorrected-log2.png)
![PCA by condition, protein legend](../figures/qc/pca/pca-by-condition-protein-raw-normalized-batchcorrected-log2.legend.png)
![PCA by batch, peptide](../figures/qc/pca/pca-by-batch-peptide-raw-normalized-batchcorrected-log2.png)
![PCA by batch, peptide legend](../figures/qc/pca/pca-by-batch-peptide-raw-normalized-batchcorrected-log2.legend.png)
![PCA by condition, peptide](../figures/qc/pca/pca-by-condition-peptide-raw-normalized-batchcorrected-log2.png)
![PCA by condition, peptide legend](../figures/qc/pca/pca-by-condition-peptide-raw-normalized-batchcorrected-log2.legend.png)

- **Raw PC1** (55% protein, 49% peptide) is driven by the high-loading AZ905/AZ906 pair and disappears after median normalization.
- **After normalization, PC1 is batch** (about 43%): the 2022 pair against all 2021 runs. The remaining components separate the 2021 pairs from each other.
- **The ComBat preview** moves the 2022 pair to the center, which is expected from centering a 2-sample batch. After that, the 2021 pair-versus-pair structure dominates.
- **Condition does not separate on PC1 or PC2 in any state, at either level.**

## 9. Spectral-count views — NSAF (as-is) and PSM counts (raw)

These are **comparators only**; the LFQ intensities above are the analysis quantity. NSAF is used as delivered: already normalized, not batch-corrected, log2 for box/correlation/PCA. PSM counts are used **raw**, with no normalization, batch correction or log transform, as they are commonly used. The PSM plots are *expected* to look poor. The box plot shows log2 PSM counts; CV, correlation and PCA use raw linear counts. Box, CV, correlation and PCA panels use the 2,168 non-contaminant protein groups with ≥ 1 PSM in all 8 runs.

**Source-precision caveat:** the Limelight dump prints NSAF ≥ 0.001 to 3 decimal places. The top-abundance NSAF values are therefore coarse: 24 distinct values across 1,629 cells, and 10 groups are identical in all runs. This shows as steps at the head of the NSAF dynamic-range curve and a small spike at CV ≈ 0. The NSAF PCA drops those 10 constant groups (the PSM PCA drops 1).

![ID depth, PSM ≥ 1](../figures/qc/id-depth/id-depth-nsaf-psm-raw-linear.png)
![ID depth, PSM ≥ 1 legend](../figures/qc/id-depth/id-depth-nsaf-psm-raw-linear.legend.png)

- **Protein groups identified (≥ 1 PSM) per run:** 2,901–3,600. The 2022 runs are about 400 lower.
- **Depth follows batch.** More groups are identified by PSMs than are quantified by LFQ (≈ 3,300 vs ≈ 2,840 median).

![Missingness, NSAF](../figures/qc/missingness/missingness-nsaf-raw-linear.png)
![Missingness, NSAF legend](../figures/qc/missingness/missingness-nsaf-raw-linear.legend.png)
![Missingness, PSM](../figures/qc/missingness/missingness-psm-raw-linear.png)
![Missingness, PSM legend](../figures/qc/missingness/missingness-psm-raw-linear.legend.png)

- **Missingness again depends on abundance** (MNAR r = 0.66 for NSAF, 0.67 for PSM).
- **2,190 groups are identified in all 8 runs.**

![Dynamic range, NSAF](../figures/qc/dynamic-range/dynamic-range-nsaf-raw-linear.png)
![Dynamic range, NSAF legend](../figures/qc/dynamic-range/dynamic-range-nsaf-raw-linear.legend.png)
![Dynamic range, PSM](../figures/qc/dynamic-range/dynamic-range-psm-raw-linear.png)
![Dynamic range, PSM legend](../figures/qc/dynamic-range/dynamic-range-psm-raw-linear.legend.png)

- **NSAF spans 5.0 orders of magnitude and PSM counts 3.3,** compared with 4.7 for LFQ.
- **The PSM tail is stepped** because counts are whole numbers.
- **Contaminants rank higher here:** CYB5 at NSAF rank 6 and pig trypsin at 13. Small proteins gain rank under NSAF's length correction.

![Box plot, NSAF](../figures/qc/abundance-boxplot/abundance-boxplot-nsaf-raw-log2.png)
![Box plot, NSAF legend](../figures/qc/abundance-boxplot/abundance-boxplot-nsaf-raw-log2.legend.png)
![Box plot, PSM](../figures/qc/abundance-boxplot/abundance-boxplot-psm-raw-log2.png)
![Box plot, PSM legend](../figures/qc/abundance-boxplot/abundance-boxplot-psm-raw-log2.legend.png)

- **NSAF per-sample medians** are already aligned, within 0.22 log2.
- **PSM counts on log2** (no normalization; no pseudocount, since every complete count is ≥ 1; revised at your request from the original linear view): per-sample medians spread over 0.53 log2. AZ905/AZ906 are highest, and AZ908 and the 2022 runs have lower medians. This is the unnormalized depth/loading offset that NSAF's total-count normalization removes (NSAF spread 0.22).

![CV, NSAF](../figures/qc/cv/cv-experimental-nsaf-raw-linear.png)
![CV, NSAF legend](../figures/qc/cv/cv-experimental-nsaf-raw-linear.legend.png)
![CV, PSM](../figures/qc/cv/cv-experimental-psm-raw-linear.png)
![CV, PSM legend](../figures/qc/cv/cv-experimental-psm-raw-linear.legend.png)

| Measure | Features | Median CV |
|---|---|---|
| LFQ intensity, raw | 1,801 | 0.397 |
| LFQ intensity, median-normalized | 1,801 | 0.268 |
| NSAF (as-is) | 2,168 | 0.345 |
| PSM counts (raw) | 2,168 | 0.343 |

![Correlation, NSAF](../figures/qc/sample-correlation/sample-correlation-nsaf-raw-log2.png)
![Correlation, NSAF legend](../figures/qc/sample-correlation/sample-correlation-nsaf-raw-log2.legend.png)
![Correlation, PSM](../figures/qc/sample-correlation/sample-correlation-psm-raw-linear.png)
![Correlation, PSM legend](../figures/qc/sample-correlation/sample-correlation-psm-raw-linear.legend.png)

- **Pair structure:** in both, every candidate pair is its own tightest cluster (NSAF r 0.958–0.967; PSM 0.990–0.995) and the 2022 pair splits off first. This is the same structure as the LFQ data.
- **The raw-count Pearson r values are inflated** by a few very high-count proteins, which is why they are higher than NSAF's.

![PCA by batch, NSAF](../figures/qc/pca/pca-by-batch-nsaf-raw-log2.png)
![PCA by batch, NSAF legend](../figures/qc/pca/pca-by-batch-nsaf-raw-log2.legend.png)
![PCA by condition, NSAF](../figures/qc/pca/pca-by-condition-nsaf-raw-log2.png)
![PCA by condition, NSAF legend](../figures/qc/pca/pca-by-condition-nsaf-raw-log2.legend.png)
![PCA by batch, PSM](../figures/qc/pca/pca-by-batch-psm-raw-linear.png)
![PCA by batch, PSM legend](../figures/qc/pca/pca-by-batch-psm-raw-linear.legend.png)
![PCA by condition, PSM](../figures/qc/pca/pca-by-condition-psm-raw-linear.png)
![PCA by condition, PSM legend](../figures/qc/pca/pca-by-condition-psm-raw-linear.legend.png)

- **PC1 is batch** in both measures (NSAF 36.5%, PSM 34.8%).
- **PC2 spreads the 2021 pairs.**
- **Condition does not separate.**

## 10. What QC says for the analysis

1. **The data read is sound.** Pairing is exact and the Limelight dumps are resolved.
2. **The candidate pairs are real structure.** Each control/raloxifene pair (905/906, 907/908, 909/910, 941/942) correlates far more closely with itself than with anything else. They most likely share a microsome source or prep. In Stage 4, **pair should be treated as a blocking factor**: a paired design, or pair as a covariate. An unpaired test would add between-pair variance to the error term and lose power. *(Recorded as a caveat finding — see below.)*
3. **Batch is a real variance component,** carried by the two 2022 runs, which are also shallower. This supports your Stage-2 decision to use batch as a covariate. Note that with pair blocking, batch is nested within pair: the 2022 pair *is* its own batch. A pair term therefore absorbs batch.
4. **No condition signal is visible in the global structure,** consistent with your near-null expectation.
5. **Missingness is abundance-dependent.** The complete-feature analysis set is biased toward abundant proteins.
6. **Run order:** caveat [0001](../findings/0001-run-order-aliased-with-condition.md) still stands. Nothing here separates drift from treatment, since there are no pooled-QC injections.
7. **The three quantities agree on the sample structure** (pair > batch > condition = nothing). NSAF precision in the export is coarse for abundant proteins, which the Stage-4 LFQ-vs-NSAF comparisons must account for.

Caveats: [0001 run order aliased with condition](../findings/0001-run-order-aliased-with-condition.md) · [0002 6-vs-2 batch structure](../findings/0002-two-acquisition-batches-6-vs-2.md) · 0003 pair structure (recorded at this gate).
