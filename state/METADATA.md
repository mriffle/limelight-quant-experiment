# Metadata — verified description

*Stage 1, 2026-09-22. Generated from `scripts/promoted/metadata_characterize.py` → `results/metadata/` and confirmed by the scientist at the Stage-1 checkpoint. Regenerate from those; never edit into inconsistency with the data.*

**Data version:** `sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74`
(combined sha256 over the six `data/` files, keyed by role; per-file hashes in `results/metadata/data_version.json`).

## Source

`data/metadata.tsv` — tab-separated, 8 data rows, 2 columns. All 17 validity checks pass (`results/metadata/validity_checks.tsv`).

## Columns

| Column | Meaning | Type / domain (validated) |
|---|---|---|
| `Replicate` | The mzML run file name for one LC-MS/MS acquisition. It is also the join key (see below). | text; unique; all 8 match `UWPRExp480_<YYYY>_<MMDD>_AZ_<NNN>_<AZnnn>_AZ_complex.mzML` |
| `condition` | Treatment arm (variable of interest). | {`control`, `raloxifene-d0`}; `raloxifene-d0` = **unlabeled (non-deuterated) raloxifene** (scientist-confirmed) |

### Fields derived from the file name (not separate metadata columns)

| Derived field | Source | Status |
|---|---|---|
| `acq_date` / `batch` | `<YYYY>_<MMDD>` → 2021-05-06 or 2022-03-18; batch label `B2021-05-06` / `B2022-03-18` | **Proxy batch.** Unknown whether the later pair was also prepared separately or only acquired later (scientist: don't know). |
| `seq_number`, `run_position_within_batch`, `run_half` | `<NNN>` (017, 019, 034 … 049) | **Presumed within-day injection order** (scientist: don't know). Not independently confirmed; mzML start timestamps could confirm it later. |
| `sample_id`, `sample_number` | `<AZnnn>` (AZ905 … AZ942) | Unique sample identifier. |
| `candidate_pair` | consecutive IDs (905/906, 907/908, 909/910, 941/942) | **Candidate matched pairs.** Each pair holds one control (odd ID) and one raloxifene (even ID). Whether they share an aliquot or incubation is unknown (scientist: don't know). |
| Constant tokens | `UWPRExp480`, `AZ`, `AZ_complex` | Same for every run, so they carry no information. |

## Experimental design

- **Human liver microsomes**, raloxifene-d0-treated vs. untreated control. **4 vs. 4 biological replicates**, one injection each (no technical replicates), **8 runs**.
- Two acquisition dates: **B2021-05-06 has 6 runs (3 control + 3 raloxifene)** and **B2022-03-18 has 2 runs (1 + 1)**.
- No other covariates are recorded (no donor, sex, age, dose or incubation time).

## Experimental vs. control samples

- **Rule:** scientist-confirmed. All 8 runs are biological samples of the two treatment arms, and neither metadata column nor file names contain a pool, QC, reference, standard or blank marker.
- **Counts:** experimental 8 (control arm 4, raloxifene-d0 arm 4); control/QC samples **0**.
- **The analysis set is all 8 runs.** No control samples exist to exclude. (Here "control" is a biological arm name, not a QC-sample role.)
- **Limitation:** there are no pooled-QC or reference runs, so technical reproducibility and run-order drift can't be measured separately from biological variation.

## Tested relationships (`results/metadata/hypotheses.tsv`, `associations.tsv`)

| # | Hypothesis | Result |
|---|---|---|
| H1 | Arms balanced 4 : 4 | holds |
| H2 | Every batch contains both conditions | holds (3/3 and 1/1). Cramér's V condition × batch = 0.00 (raw and bias-corrected) |
| H3 | Each candidate pair holds exactly one control and one raloxifene; control ↔ odd ID | holds for all 4 pairs / 8 samples. V condition × pair = 0.00 |
| H4 | Run order is independent of condition | **fails.** In every batch, every control was run before every raloxifene sample. B2021-05-06: U = 9/9, exact permutation p = 0.10 two-sided (2/20). Stratified over both batches: U = 10/10, rank-biserial = 1.00 (control earlier), **exact p = 0.05 two-sided** (2/40 permutations). This is the most extreme arrangement possible. The direction was not specified in advance, so the two-sided p is the one reported; the one-sided value in the observed direction is 0.025. V condition × run-half = 1.00. |
| H5 | Samples run in sample-ID order | no (B2021-05-06 ran 905, 907, 909, 906, 908, 910) |

With n = 8, Cramér's V is unstable in both forms. The exact permutation result (H4) is the main evidence on run order.

Run layout (`results/metadata/run_layout.tsv`):

| Batch | pos 1 | pos 2 | pos 3 | pos 4 | pos 5 | pos 6 |
|---|---|---|---|---|---|---|
| B2021-05-06 | ctrl AZ905 (034) | ctrl AZ907 (039) | ctrl AZ909 (041) | ralox AZ906 (045) | ralox AZ908 (047) | ralox AZ910 (049) |
| B2022-03-18 | ctrl AZ941 (017) | ralox AZ942 (019) | | | | |

## Imbalances, skews, confounds (caveat findings)

- **Run order is aliased with condition** (H4). If the sequence numbers are injection order, any drift over the run (instrument sensitivity, LC column, carry-over) coincides exactly with the treatment contrast and can't be separated from it. There are no pooled QC injections to measure drift. → caveat finding **0001** (to be recorded).
- **Batch structure** (H2). Condition is balanced within each batch, so batch doesn't bias the contrast. But batch sizes are unequal (6 vs. 2), the batch meaning is unknown (prep + acquisition, or acquisition only), and a ten-month gap between acquisitions can bring a large batch effect. That effect would add variance and could dominate PCA and clustering. The later batch contains a single pair, so a batch effect can't be estimated apart from that one pair's own variation. Batch should enter differential models as a covariate and be checked in QC. → caveat finding **0002** (to be recorded).
- **Candidate pairing** (H3). Not a confound, but a possible blocking factor. Stage 4 should treat a paired analysis (pair as a covariate) as a sensitivity analysis next to the default unpaired one.
- **Small n.** With 4 vs. 4, all differential results are exploratory. Stage 0 records this dataset as a near-null test (little true difference expected).

## Join key to the data matrices (`results/metadata/join_key.tsv`)

`Replicate` (mzML file name) ↔ `search_scan_file_id` via `data/scan-file-search-id-mapping.txt` (a verified bijection, 8 ↔ 8) ↔ quant-file column names:

| Sample | Condition | Batch | search_scan_file_id | Quant columns |
|---|---|---|---|---|
| AZ905 | control | B2021-05-06 | 18461 | `Intensity_search_scan_file_id_18461`, `Detection Type_search_scan_file_id_18461` |
| AZ907 | control | B2021-05-06 | 18459 | `…_18459` |
| AZ909 | control | B2021-05-06 | 18454 | `…_18454` |
| AZ906 | raloxifene-d0 | B2021-05-06 | 18460 | `…_18460` |
| AZ908 | raloxifene-d0 | B2021-05-06 | 18456 | `…_18456` |
| AZ910 | raloxifene-d0 | B2021-05-06 | 18457 | `…_18457` |
| AZ941 | control | B2022-03-18 | 18455 | `…_18455` |
| AZ942 | raloxifene-d0 | B2022-03-18 | 18458 | `…_18458` |

- `data/protein-quants.tsv` (`Intensity_…`) and `data/peptide-quants.tsv` (`Intensity_…`, `Detection Type_…`) each carry exactly these 8 ids (verified from the headers).
- **Unresolved:** the Limelight table dumps (`protein-limelight-table-dump.txt`, `peptide-limelight-table-dump.txt`) label their per-run columns with 8 truncated, ambiguous suffixes: `1_0`, `1_05`, `1_050`, `1_0506`, `1_0506_`, `1_0506_A`, `2_0318_`, `2_0318_A`. The labels do show which batch a column belongs to (6 labels start `1_`, 2 start `2_`), but not which sample within the batch. The mapping to runs must be established by matching values, e.g. dump `Quant (FlashLFQ)` against the quant-file intensities. That happens in Stage 2 and is verified at the Stage 3 integrity gate. Until then, PSM counts and NSAF can't be joined to samples.

## Figures

- `figures/metadata/distributions/cohort-counts` — counts per condition / batch / candidate pair.
- `figures/metadata/crosstabs/0002-condition-by-batch`, `figures/metadata/crosstabs/0001-condition-by-run-half`.
- `figures/metadata/run-layout/0001-run-layout` — run order by batch, colored by condition.
- Table 1: `results/metadata/table1.md`.
