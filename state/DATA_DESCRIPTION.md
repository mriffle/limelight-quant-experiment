# Data description — verified

*Stage 2, 2026-09-22. Derived from `scripts/scratch/stage2_explore.py` → `results/stage2/` plus the checks recorded below. The preprocessing decisions were confirmed by the scientist. Regenerate from those; never edit into inconsistency with the data.*

**Data version:** `sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74` (same as `state/METADATA.md`).

## Files

| File | Role | Orientation / shape |
|---|---|---|
| `data/protein-quants.tsv` | FlashLFQ **protein** intensities (primary protein matrix) | rows = proteins (4,344), columns = 3 annotation + 8 sample columns |
| `data/peptide-quants.tsv` | FlashLFQ **peptide** intensities + per-run detection type (primary peptide matrix) | rows = peptides (33,358), columns = 5 annotation + 8 `Intensity_…` + 8 `Detection Type_…` |
| `data/protein-limelight-table-dump.txt` | Limelight protein table: **PSMs, NSAF**, FlashLFQ quant (display-rounded), Limelight quant | rows = protein groups (4,344), 8 runs × 4 measures + group number |
| `data/peptide-limelight-table-dump.txt` | Limelight peptide table: PSMs, quant, MBR flag, shared-group flag per run | rows = modified peptides (33,486) |

All four files are features × samples (wide). They are tab-separated ASCII with single-row headers, no embedded metadata and no locale or decimal-comma issues. Every quant-file sample column parses as numeric; the one text token is `NaN` (below).

## Sample identifiers → metadata

- **Quant files:** `Intensity_search_scan_file_id_<id>` and `Detection Type_search_scan_file_id_<id>`. Each `<id>` maps to a metadata `Replicate` through `scan-file-search-id-mapping.txt`, a verified 8 ↔ 8 bijection (`state/METADATA.md`).
- **Limelight dumps: resolved.** The per-run suffixes are truncated labels that follow neither run order nor ID order. They were mapped by value correspondence. Each dump `Quant (FlashLFQ)` column equals exactly one protein-quants run once that run is rounded to 3 significant figures (100% of 4,317 matched proteins). No other run exceeds 44% agreement. The peptide dump independently confirms the same mapping (log-r = 1.000, median relative difference 0.07% ≈ rounding). Map: `results/stage2/limelight_label_map.tsv`.

| Dump label | Sample | Condition | Batch |
|---|---|---|---|
| `1_0506_A` | AZ905 | control | B2021-05-06 |
| `1_050` | AZ907 | control | B2021-05-06 |
| `1_05` | AZ909 | control | B2021-05-06 |
| `1_0506_` | AZ906 | raloxifene-d0 | B2021-05-06 |
| `1_0506` | AZ908 | raloxifene-d0 | B2021-05-06 |
| `1_0` | AZ910 | raloxifene-d0 | B2021-05-06 |
| `2_0318_` | AZ941 | control | B2022-03-18 |
| `2_0318_A` | AZ942 | raloxifene-d0 | B2022-03-18 |

- **Control/QC samples:** none in the data either. Every sample column maps to one of the 8 experimental runs, and no pool, QC, reference or blank columns exist. This matches Stage 1.

## Feature identifiers

- **Protein (quants):** `psvid_<n>_sp|<accession>|<ENTRY>_HUMAN` (`psvid` = Limelight protein-sequence-version id). All 4,344 are unique, all use the `sp` database, and all are human. `Gene Name` and `Organism` are empty for every row. No row contains `;`: each row is a single protein, and the quant file splits indistinguishable groups into separate rows.
- **Protein (dump):** `sp|<accession>|<ENTRY>` without the `psvid` prefix. 27 dump rows are **comma-joined indistinguishable groups** (e.g. `sp|P0DPH7|TBA3C_HUMAN,sp|P0DPH8|TBA3D_HUMAN`). For the PSM/NSAF join, map a group row to its member quant rows, one group → many rows. 4,317 of 4,344 rows match by exact accession.
- **Peptide `Sequence`** = base sequence + **one bracketed total modification mass appended at the C-terminal end** (e.g. `…LCR[+114.042928]` = 2 × carbamidomethyl; `[+0.0]` = unmodified). This means modification *positions* are not encoded: positional isoforms are collapsed into one row. Masses present:

  | Mass | Rows | Composition |
  |---|---|---|
  | +0.0 | 26,267 | unmodified |
  | +57.02 | 4,946 | carbamidomethyl |
  | +15.99 | 802 | oxidation |
  | +471.15 | 48 | presumed raloxifene-derived adduct |
  | +528.17 | 5 | +471.15 + carbamidomethyl |

  Other combinations also occur. 33,358 unique sequences, 32,008 base sequences.
- **Peptide `Protein Groups`:** `;`-joined `psvid_…` members for shared peptides (2,711 rows). The protein-quants rows are a subset of the peptide protein-group members (4,344/4,344 present).
- **No spreadsheet corruption:** identifiers are accessions and entry names (no gene-symbol-as-date risk; `Gene Name` is empty), read as text, and no scientific-notation or leading-zero issues were found.

## Contaminants and decoys — decision

- **Decoys:** none. There are no `DECOY`/`REV`/`rev_` entries (Percolator-filtered output).
- **Contaminants:** **52 protein entries** use a no-accession form, `sp|<ENTRY>_HUMAN|` (keratins K1H2/K2C1/KRHB4, HBA, ALBU, TRFE, SODC, ANXA5, CO5, …). This looks like an appended common-contaminants FASTA, and some duplicate real accessions (`sp|ALBU_HUMAN|` alongside `sp|P02768|ALBU_HUMAN`). They account for 1.1–1.7% of summed protein intensity per run. The one plain `CON` substring hit, ACON_HUMAN (aconitase), is **not** a contaminant.
  - **Decision (scientist-confirmed):** **flag the no-accession form as contaminant and exclude it** from normalization and all differential analysis. Keep it visible in QC as contaminant share per run. At peptide level, exclude any peptide whose `Protein Groups` includes a contaminant entry (675 peptides).
  - Detection rule (regex on a protein id / group member): `^psvid_\d+_sp\|[^|]+\|$`.

## Transformation / normalization state

- **Scale: `linear`.** Values are raw MS1 intensities: non-integer, strictly positive when quantified, and spanning about 3e4 to 5e9 (protein; log10 range 4.0–5.1 per run). Median log2 is about 23 (protein) and about 21.5 (peptide). No negatives. The loader records `scale = "linear"`.
- **Not normalized.** On features quantified in all runs, per-sample median log2 offsets span **−0.39 … +0.66 (protein)** and **−0.51 … +0.66 (peptide)**, about 1.5-fold; column sums vary from 0.66 to 1.53 × the mean. If FlashLFQ's own normalization had been on, these offsets would be near zero.
  - *Observation for QC:* the two highest-offset runs are **AZ905 (+0.66) and AZ906 (+0.53)**, which form one candidate pair (P905_906). This hints at a pair or loading effect, relevant to the pairing question in [caveat 0002](../findings/0002-two-acquisition-batches-6-vs-2.md). To examine in Stage 3.
- **Match-between-runs is on.** 25,105 peptide values (14% of positive peptide cells) are `MBR`, and FlashLFQ protein intensities include them. They are kept as delivered.

## Missing-value semantics

| File | Token | Count | Meaning |
|---|---|---|---|
| protein-quants | `0` | 13,084 cells (30–48% per run) | not quantified in that run (FlashLFQ convention), **not** a true zero |
| protein-quants | `NaN` (literal text) | 14 cells / 7 proteins | FlashLFQ "not quantifiable across runs"; the dump shows `FlashLFQ: NaN — not quantifiable across runs` |
| peptide-quants | `0` | 89,181 cells | not quantified; the reason is in `Detection Type` (below) |
| dumps | `''` (blank) | e.g. 8,811 Limelight-quant cells | no PSMs in that run |
| dumps | `overlapping signal` | 3,424 `Quant (Limelight)` cells | Limelight's own quant not reported |

- `NaN` and `0` are both treated as **missing (not quantified)**. They are recorded as distinct tokens, and neither is a measured zero.
- **Peptide `Detection Type` fully determines quantification:**

  | Detection type | Cells | Intensity |
  |---|---|---|
  | `MSMS` | 152,578 | always > 0 |
  | `MBR` | 25,105 | always > 0 |
  | `NotDetected` | 75,363 | always 0 |
  | `MSMSIdentifiedButNotQuantified` | 10,135 | always 0 (identified, no MS1 peak) |
  | `MSMSAmbiguousPeakfinding` | 3,683 | always 0 |

- **Missingness structure:**
  - The 2022 runs (AZ941/AZ942) have about 48% missing vs. 31–38% for 2021, and their minimum quantified intensity is about 3–7× higher. That points to a **shallower batch with a higher detection floor**, consistent with left-censoring (MNAR); the Stage-3 MNAR diagnostic will confirm.
  - Proteins by number of runs quantified: 0: 921 · 1: 157 · 2: 219 · 3: 160 · 4: 212 · 5: 209 · 6: 381 · 7: 280 · 8: 1,805.
  - **921 proteins (21%) have no intensity in any run** even though they were identified: median 8 PSMs in total, maximum 2,719. They are present in the dump with PSMs/NSAF but have no LFQ value, which matters for the LFQ-vs-spectral-count comparison.
  - Peptides: 1,723 all-zero rows; 10,392 quantified in all 8 runs.

## Preprocessing decisions (scientist-confirmed, 2026-09-22)

1. **Normalization: median.** Scale each run on the linear scale so its median matches, computed over non-contaminant features quantified in all runs. Then log2 for analysis. `lib/common/normalize`, method `median`.
2. **Missing values: complete features only.** Normalization and differential analysis use features quantified in all 8 runs, with no imputation:
   - **1,789 proteins** and **10,121 peptides** after contaminant exclusion (`handle_missing` with `max_missing_fraction = 0`).
   - Revisit against the Stage-3 `missingness` MNAR diagnostic. A filter + left-censored imputation run is a candidate Stage-4 sensitivity analysis.
   - QC figures (ID depth, missingness, dynamic range) use the full raw matrices.
3. **Batch axis: `batch`** (the acquisition-date proxy, B2021-05-06 vs B2022-03-18).
   - Condition is balanced within batch (V = 0), so there is no batch ↔ condition confound.
   - **Batch enters the differential model as a covariate** (condition + batch). Results are reported with and without batch adjustment.
   - Stage-3 QC previews **batch-label-only** ComBat for PCA/CV visuals only.
   - Standing caveats: [0001 run order aliased with condition](../findings/0001-run-order-aliased-with-condition.md); [0002 6-vs-2 batch structure](../findings/0002-two-acquisition-batches-6-vs-2.md).
4. **Contaminants:** excluded as above. **Decoys:** none present.
5. **Match-between-runs values:** kept as delivered, not separately filtered.

## Known data-quality issues

- The Limelight dump numbers are **display-rounded** (3 significant figures; thousands separators in PSMs). Use the dumps only for PSMs/NSAF (exact integer PSMs; NSAF to 3 significant figures). Never use them as an intensity source: `protein-quants.tsv`/`peptide-quants.tsv` are the full-precision intensities.
- NSAF is shown to 3 significant figures (`0.018`, `2.34e-5`; `0` = no PSMs, 8,811 cells), which is adequate for rank and correlation comparisons. Column sums are 0.982–0.993 rather than 1, probably from rounding and/or rows outside the export.
- The two indistinguishable-group representations differ: the dump uses comma-joined groups, and quants use split rows.
- Peptide modification positions are not encoded (positional isoforms are collapsed).
- The 2022 batch is shallower (fewer IDs, higher intensity floor).
