# Project framing

*Stage 0 — stated by the scientist, 2026-09-22. This is a living document: it will be revised as understanding of the data deepens (e.g. after Stages 1–2), with changes noted below.*

## Domain

Bottom-up, data-dependent acquisition (DDA) mass-spectrometry proteomics with **label-free quantification (LFQ)**. Peptide spectrum matches come from an existing **Comet** database search post-processed with **Percolator**; label-free quantities were computed with **FlashLFQ** on top of those identifications. The data are exported from **Limelight** (the proteomics results platform), both as quant tables and as Limelight table dumps.

The biological system is **human liver microsomes (HLM)**, treated with **raloxifene** or left untreated (controls).

## What is being examined and why

This is a **methods-development study**, not a biological investigation. The scientist is developing methods to provide **QC plots and high-level exploration and analysis** of LFQ proteomics results — feature finding, differential analysis, and related views. The raloxifene/HLM dataset serves as a realistic example on which to build and assess those methods.

The biology (raloxifene is a known mechanism-based inactivator of CYP3A4) is context, not the target of inquiry. Biological interpretation of any differences is out of scope unless it bears on judging whether a method behaves sensibly.

## Experimental design

- **Two conditions:** raloxifene-treated (`raloxifene-d0` — unlabeled, non-deuterated raloxifene) vs. control (untreated) human liver microsomes.
- **Replication:** 4 biological replicates per condition; **no technical replicates** → 8 runs total (4 vs. 4). *(Corrected at Stage 1: Stage 0 stated 3 vs. 3; the scientist confirmed 4 vs. 4 against the metadata.)*
- Manipulated: raloxifene treatment. Measured: peptide- and protein-level LFQ intensities (FlashLFQ), plus identification-derived measures (PSM counts, NSAF).
- Batch / run order / pairing: not recorded by the scientist. The run file names encode two acquisition dates (six runs 2021-05-06, two runs 2022-03-18), a sequence number, and consecutive control/raloxifene sample-ID pairs — see `state/METADATA.md`. Treatment details (dose, incubation time) not stated.

**Power note:** n = 4 per group is a small design. Differential tests will have limited power; effect-size estimates will be imprecise, and variance moderation / multiplicity handling will matter a great deal.

## Scientific goals

Exploratory/descriptive and methodological, at **both the protein level and the peptide level**:

1. **Standard QC** — e.g. identification counts per run, missing-value structure, intensity distributions, replicate agreement/correlation.
2. **PCA** (and related sample-level structure views) — do samples group by condition, or does run-level variation dominate?
3. **Differential abundance with volcano plots** — raloxifene vs. control.
4. **Compare the LFQ quants against spectral-counting measures** — PSM counts and NSAF values (present in some of the data files) — to assess how well the MS1-intensity–based quantities agree with count-based abundance estimates.

**Primary quantity (scientist, 2026-09-22):** the **FlashLFQ MS1-intensity (LFQ) quants are the quantity for all analysis** (differential abundance, PCA, etc.). NSAF and PSM counts are **comparators only**. They appear in QC as parallel views (NSAF as-is; PSM counts raw, as commonly used). In Stage 4 they get head-to-head comparisons against the LFQ-derived quants, e.g. per-protein correlation scatter plots.

**Prior expectation:** the scientist expects **little or no true difference** between treated and control. An in vitro microsome incubation should not substantially change protein abundance. This dataset therefore acts approximately as a **near-null test** of the differential pipeline: a well-calibrated method should report few or no confident hits, and the p-value distribution should look roughly uniform. A large number of "significant" proteins would more likely point to a technical artifact (normalization, missing-value handling, run effects) than to biology.

## Guarding against motivated reasoning

Stating goals and expectations up front helps keep the work relevant, but it can also bias the work. Two ways that could happen here:

- **"No difference" is itself a stated expectation.** It must not become a reason to dismiss real signal or to tune methods until the volcano plot looks empty. If differences show up, they are recorded and investigated on their merits (technical cause or real effect), not explained away.
- Because the purpose is methods development, a method looking good on this dataset is **not** evidence that it is correct. Where possible, method behavior should be judged against checkable properties (null calibration, agreement with spectral counts, replicate consistency) rather than against how the plots look.

Exploration is generous and promotion is ruthless: any finding still has to pass the skepticism gates and independent validation, whatever it says relative to these hopes.

## Revision log

- 2026-09-22 — Initial framing (Stage 0).
- 2026-09-22 — Stage 1: replication corrected from 3 vs. 3 to 4 vs. 4 (scientist-confirmed); "d0" = unlabeled raloxifene; batch/run order/pairing unknown.
- 2026-09-22 — Stage 3: LFQ intensities fixed as the analysis quantity; NSAF/PSM counts are comparators (QC views + Stage-4 head-to-head).
