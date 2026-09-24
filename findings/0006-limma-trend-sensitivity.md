---
id: 6
title: "Adding a mean-variance trend to the moderated prior yields at most one protein hit (ATPK) and 10 peptide hits, none distinguishable from within-pair relabelling"
status: candidate
phase: exploratory
kind: discovery
created: 2026-09-24
updated: 2026-09-24

summary: "A limma-trend prior (mean-variance trend in the empirical-Bayes variance) changes the near-null raloxifene-d0 vs control result only for the pair-blocked LFQ model: 1 protein hit (ATPK, q 0.040) and 10 peptide hits at q < 0.05 — the same 10 that were q < 0.10 without a trend — while NSAF, PSM, and all batch/unadjusted fits stay at 0 hits. None of the hits survive within-pair relabelling above the 1/8 permutation floor."
verdict: "The trend prior is a defensible model (residual SD falls with abundance, Spearman −1.0), but it was adopted after ATPK was seen near the top, and choosing it is what produces the hits — a post hoc model choice layered on 12 fits. The result remains effectively null: the observed labelling cannot beat the p = 1/8 floor of within-pair relabelling, and every hit is confounded with run order."

entities:
  - { db: uniprot, id: "P56134", label: "ATPK" }

relationships:
  - { type: refines, target: 4, note: "Sharpens the near-null verdict of 0004: a trend prior surfaces one protein and ten peptides, but they do not clear the relabelling floor, so the near-null conclusion holds." }
  - { type: relates_to, target: 1, note: "Run order is aliased with condition, so any hit here is condition OR run order." }
  - { type: relates_to, target: 3, note: "The hits exist only in the pair-blocked design, the structure documented in 0003." }

provenance:
  data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
  script: { path: "scripts/scratch/de_trend_sensitivity.py", commit: "fc32cbc" }
  params:
    module: "scripts/scratch/analysis/differential_abundance.py (method='moderated', trend=True)"
    sample_set: "experimental"
    n_experimental: 8
    n_controls_excluded: 0
    contrast: "condition: raloxifene-d0 vs control (positive log2FC = higher in raloxifene-d0)"
    method: "limma-style moderated t, eBayes(trend=TRUE): prior variance = natural cubic spline (4 df incl. intercept) of mean log2 abundance fitted to log residual variances (limma fitFDist with covariate); one d0; non-robust"
    designs: { paired: "condition + candidate_pair (primary)", batch: "condition + batch (sensitivity)", unadjusted: "condition only (sensitivity)" }
    quantities: "protein (LFQ, 1801 complete), peptide (LFQ), nsaf, psm_log2 — all 4 x 3 designs"
    residual_df_paired: 3
    correction: "Benjamini-Hochberg, per quantity x design over the features tested (contrast term)"
    ci: "95% t interval on the moderated SE (df = residual_df + d0)"
    pi0_estimator: "Storey, fixed lambda = 0.5"
    relabelling_diagnostic: "8 within-pair labellings (observed + 7 swaps); min attainable p = 1/8"
    validation: "matches R limma 3.58.1 lmFit + eBayes(trend=TRUE) on all 12 fits to ~1e-12 (r_agreement.json)"
    outputs: "results/de/raloxifene-vs-control/trend/{summary.json, <quantity>_<design>.tsv, results/}"
  environment: "pyproject.toml + uv.lock (Python 3.12.3)"
  seeded_from: { template: "differential-abundance", version: "0.1" }
  seed: null
  result_id: null

evidence:
  - metric: "log2FC"
    value: -0.516
    ci: [-0.6488, -0.3833]
    p_value: 2.19e-05
    p_adjusted: 0.0395
    correction: "BH"
    test: "limma-trend moderated t, paired (condition + candidate_pair), df 3"
    n: 8
    note: "ATPK (P56134), LFQ protein, paired, trend=TRUE — the single protein hit; q 0.21 without a trend. residual SD 0.0166, trend prior SD 0.1025 at mean log2 abundance 26.7."
  - metric: "n_hits_q<0.05 (protein, LFQ, paired)"
    value: 1
    ci: null
    p_value: null
    p_adjusted: null
    correction: "BH"
    test: "limma-trend, trend=TRUE"
    n: 1801
    note: "1 down, 0 up. Without a trend: 0 hits at q<0.05."
  - metric: "n_hits_q<0.05 (peptide, LFQ, paired)"
    value: 10
    ci: null
    p_value: null
    p_adjusted: null
    correction: "BH"
    test: "limma-trend, trend=TRUE"
    n: null
    note: "3 up, 7 down. These are exactly the 10 peptides at q<0.10 without a trend (no-trend q 0.052-0.077 -> trend q 0.028-0.042). Without a trend: 0 hits at q<0.05."
  - metric: "n_hits_q<0.05 (NSAF / PSM, all designs; protein & peptide batch/unadjusted)"
    value: 0
    ci: null
    p_value: null
    p_adjusted: null
    correction: "BH"
    test: "limma-trend, trend=TRUE"
    n: null
    note: "NSAF and PSM-count at 0 hits under every design; LFQ protein & peptide at 0 hits under batch and unadjusted. The trend changes the answer only for the pair-blocked LFQ fits."
  - metric: "trend prior SD range vs abundance (protein, paired)"
    value: null
    ci: [0.0787, 0.3686]
    p_value: null
    p_adjusted: null
    correction: null
    test: "eBayes trend spline; Spearman(prior SD, mean abundance) = -1.0"
    n: 1801
    note: "Trend prior SD falls monotonically 0.369 -> 0.079 with abundance (d0 4.79) vs the constant no-trend prior SD 0.137 (d0 3.01). It inflates low-abundance SEs and shrinks high-abundance ones, so the hits are high-abundance features."
  - metric: "within-pair relabelling rank (observed labelling)"
    value: null
    ci: null
    p_value: 0.125
    p_adjusted: null
    correction: null
    test: "8 within-pair labellings (observed + 7 swaps); min attainable p = 1/8"
    n: 8
    note: "Observed labelling ranks 1/8 on hits for both protein (1 hit; all 7 relabellings 0) and peptide (10 hits; all 7 relabellings 0-1) but a within-pair permutation cannot go below 1/8. Storey pi0 ranks 1/8 (protein), 2/8 (peptide)."

figures:
  - png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-prior-sd-trend-protein-paired.png"
    svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-prior-sd-trend-protein-paired.svg"
    legend_png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-prior-sd-trend-protein-paired.legend.png"
    legend_svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-prior-sd-trend-protein-paired.legend.svg"
    caption: "Per-protein residual SD vs mean log2 abundance (LFQ, paired, n=1801), with the trend prior SD curve (falls 0.369->0.079) against the constant no-trend prior SD (0.137); ATPK highlighted at residual SD 0.017, trend prior SD 0.10."
    script: { path: "scripts/scratch/fig_trend_atpk.py", commit: "fc32cbc" }
    data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
    result_id: null
    params: { design: "paired", quantity: "protein" }
  - png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-protein-paired.png"
    svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-protein-paired.svg"
    legend_png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-protein-paired.legend.png"
    legend_svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-protein-paired.legend.svg"
    caption: "BH q with a trend (x) vs without a trend (y) per protein (LFQ, paired), points colored by mean log2 abundance; q=0.05/0.10 guide lines. Only ATPK crosses q<0.05 under the trend."
    script: { path: "scripts/scratch/fig_trend_atpk.py", commit: "fc32cbc" }
    data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
    result_id: null
    params: { design: "paired", quantity: "protein" }
  - png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-peptide-paired.png"
    svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-peptide-paired.svg"
    legend_png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-peptide-paired.legend.png"
    legend_svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-peptide-paired.legend.svg"
    caption: "BH q with a trend (x) vs without a trend (y) per peptide (LFQ, paired), colored by mean log2 abundance; q=0.05/0.10 guide lines. The 10 peptides that were q<0.10 without a trend move to q<0.05 with it."
    script: { path: "scripts/scratch/fig_trend_atpk.py", commit: "fc32cbc" }
    data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
    result_id: null
    params: { design: "paired", quantity: "peptide" }
  - png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-protein-trend-paired.png"
    svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-protein-trend-paired.svg"
    legend_png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-protein-trend-paired.legend.png"
    legend_svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-protein-trend-paired.legend.svg"
    caption: "Volcano of raloxifene-d0 vs control, LFQ protein, paired, trend=TRUE — log2FC (x) vs -log10(BH q) (y), n=1801; q=0.05/0.10 guide lines; the single hit (ATPK, log2FC -0.52) labelled."
    script: { path: "scripts/scratch/fig_trend_atpk.py", commit: "fc32cbc" }
    data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
    result_id: null
    params: { design: "paired", quantity: "protein" }
  - png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-peptide-trend-paired.png"
    svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-peptide-trend-paired.svg"
    legend_png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-peptide-trend-paired.legend.png"
    legend_svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-peptide-trend-paired.legend.svg"
    caption: "Volcano of raloxifene-d0 vs control, LFQ peptide, paired, trend=TRUE — log2FC (x) vs -log10(BH q) (y); q=0.05/0.10 guide lines; the 10 q<0.05 hits (3 up, 7 down) labelled."
    script: { path: "scripts/scratch/fig_trend_atpk.py", commit: "fc32cbc" }
    data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
    result_id: null
    params: { design: "paired", quantity: "peptide" }
  - png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-relabel-diagnostic-trend.png"
    svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-relabel-diagnostic-trend.svg"
    legend_png: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-relabel-diagnostic-trend.legend.png"
    legend_svg: "figures/analysis/differential-abundance/raloxifene-vs-control/0006-relabel-diagnostic-trend.legend.svg"
    caption: "Within-pair relabelling diagnostic (trend=TRUE): hit count per quantity (protein/peptide/NSAF/PSM) across the 8 within-pair labellings, observed labelling highlighted vs the 7 swaps. Observed ranks 1/8 for protein and peptide; NSAF/PSM are 0 for all labellings."
    script: { path: "scripts/scratch/fig_trend_atpk.py", commit: "fc32cbc" }
    data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
    result_id: null
    params: {}

references: []

validation:
  computational_reproduction: { status: not_attempted }
  analytic_replication:       { status: not_attempted }
  data_replication:           { status: not_attempted }

integrity_signoff: true
---

# Adding a mean-variance trend to the moderated prior yields at most one protein hit (ATPK) and 10 peptide hits, none distinguishable from within-pair relabelling

## Summary
The primary raloxifene-d0 vs control result (near-null; see [finding 0004](0004-no-differential-abundance-raloxifene-vs-control.md)) used a moderated-t model with a constant empirical-Bayes prior. Refitting all 4 quantities × 3 designs with a **mean-variance trend** in the prior (limma-trend, matched to R limma 3.58.1 to ~1e-12 on all 12 fits) changes the answer only for the pair-blocked LFQ fits: **1 protein hit** (ATPK, P56134, q 0.040) and **10 peptide hits** at q < 0.05 — exactly the 10 peptides that were q < 0.10 without a trend. NSAF, PSM counts, and every batch/unadjusted fit stay at 0 hits. None of these hits can beat the p = 1/8 floor of within-pair relabelling.

## Verdict
The trend prior is a **defensible model** — residual SD falls monotonically with abundance (Spearman −1.0), which is exactly the condition a trend is designed for — but it was chosen **after** ATPK was seen near the top of the no-trend ranking, and choosing it is what produces the hits. This is a post hoc model choice layered on the 12 quantity × design fits already run. The result stays effectively null: the observed labelling ranks best possible but cannot go below 1/8, and every hit is aliased with run order ([finding 0001](0001-run-order-aliased-with-condition.md)).

## Evidence

**A trend prior falls with abundance, and that is what lets high-abundance features become hits.** The empirical-Bayes prior SD, held constant at 0.137 in the primary model, becomes a curve falling 0.369 → 0.079 across the abundance range under a trend (d0 4.79 vs 3.01). It therefore inflates the SE of low-abundance proteins and shrinks it for high-abundance ones. ATPK — a high-abundance protein with a tiny residual SD (0.017) — is moderated toward a *smaller* prior SD (0.10) than the constant model would use, so its already-tight fold change becomes significant.

![Per-protein residual SD vs mean log2 abundance (LFQ, paired, n=1801) with the trend prior SD curve falling 0.369 to 0.079 against the constant no-trend prior SD 0.137; ATPK highlighted.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-prior-sd-trend-protein-paired.png)

![Legend for Figure 1 — residual SD points, the trend prior curve (blue), the constant no-trend prior (orange), and the highlighted ATPK marker.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-prior-sd-trend-protein-paired.legend.png)

*Figure 1. Residual SD vs abundance with both priors. Produced by `scripts/scratch/fig_trend_atpk.py` (fc32cbc) from data `sha256:bc6b73d3…`.*

Each point is one protein at its mean abundance (x) and residual SD (y). The trend prior curve slopes down to the right while the no-trend prior is a flat line; ATPK sits far right with a residual SD well below both priors, so under the trend it is moderated toward a smaller variance and clears q < 0.05 — the mechanism, not a biological signal, is what makes it a hit.

**Only ATPK crosses the protein significance line under the trend.** Plotting the BH q with a trend against the q without one, per protein, shows the whole proteome hugging the diagonal except ATPK, which drops from q 0.21 (no trend) to q 0.040 (trend).

![BH q with a trend (x) vs without a trend (y) per protein (LFQ, paired), colored by mean log2 abundance, with q=0.05/0.10 guide lines; only ATPK crosses q<0.05.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-protein-paired.png)

![Legend for Figure 2 — points colored by mean log2 abundance (viridis), q=0.05 and q=0.10 guide lines.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-protein-paired.legend.png)

*Figure 2. Protein q, trend vs no-trend. Produced by `scripts/scratch/fig_trend_atpk.py` (fc32cbc) from data `sha256:bc6b73d3…`.*

Read the lower-left corner: a single high-abundance point (ATPK) sits left of the vertical q=0.05 line but above the horizontal one — significant with a trend, not without — while every other protein stays on the near-diagonal cloud above q=0.05 in both. One protein moving is the entire protein-level effect of the trend.

![Volcano of raloxifene-d0 vs control, LFQ protein, paired, trend=TRUE — log2FC (x) vs -log10(BH q) (y), n=1801, ATPK labelled.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-protein-trend-paired.png)

![Legend for Figure 3 — significance classes (up / down / not significant) and the q=0.05/0.10 guide lines.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-protein-trend-paired.legend.png)

*Figure 3. Protein volcano, trend. Produced by `scripts/scratch/fig_trend_atpk.py` (fc32cbc) from data `sha256:bc6b73d3…`.*

The volcano confirms the count: a single labelled point (ATPK, log2FC −0.52) rises above the q=0.05 line, with the rest of the proteome in the gray cloud below — a one-protein result, not a coherent signature.

**The peptide hits are the no-trend borderline cases, promoted.** The 10 peptides that reach q < 0.05 with a trend are precisely the 10 that sat at q < 0.10 without one (no-trend q 0.052–0.077 → trend q 0.028–0.042); no-trend had 0 hits at q < 0.05.

![BH q with a trend (x) vs without a trend (y) per peptide (LFQ, paired), colored by mean log2 abundance, q=0.05/0.10 guide lines.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-peptide-paired.png)

![Legend for Figure 4 — points colored by mean log2 abundance (viridis), q=0.05 and q=0.10 guide lines.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-q-trend-vs-notrend-peptide-paired.legend.png)

*Figure 4. Peptide q, trend vs no-trend. Produced by `scripts/scratch/fig_trend_atpk.py` (fc32cbc) from data `sha256:bc6b73d3…`.*

The points between the horizontal q=0.10 and vertical q=0.05 guide lines are the promoted set: they were just under 0.10 without a trend and just under 0.05 with it — a threshold nudge, not a new population of effects.

![Volcano of raloxifene-d0 vs control, LFQ peptide, paired, trend=TRUE — log2FC (x) vs -log10(BH q) (y); the 10 q<0.05 hits labelled (3 up, 7 down).](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-peptide-trend-paired.png)

![Legend for Figure 5 — significance classes (up / down / not significant) and the q=0.05/0.10 guide lines.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-volcano-peptide-trend-paired.legend.png)

*Figure 5. Peptide volcano, trend. Produced by `scripts/scratch/fig_trend_atpk.py` (fc32cbc) from data `sha256:bc6b73d3…`.*

Ten labelled points (3 up, 7 down) clear the q=0.05 line; all are modest fold changes just past threshold, consistent with borderline cases being pushed over rather than a strong response emerging.

**None of the hits survives within-pair relabelling above the permutation floor.** With four matched pairs, there are 8 within-pair labellings (observed + 7 swaps), so the smallest attainable p is 1/8. The observed labelling produces the most hits for both protein (1; all 7 swaps give 0) and peptide (10; the 7 swaps give 0 or 1) and the lowest Storey π0 — but "ranks first of 8" is the best a permutation *can* do and still means p = 1/8. NSAF and PSM counts give 0 hits under every labelling.

![Within-pair relabelling diagnostic (trend=TRUE): hit count per quantity across the 8 within-pair labellings, observed labelling vs the 7 swaps.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-relabel-diagnostic-trend.png)

![Legend for Figure 6 — observed labelling (orange) vs within-pair relabellings (gray), by quantity.](../figures/analysis/differential-abundance/raloxifene-vs-control/0006-relabel-diagnostic-trend.legend.png)

*Figure 6. Relabelling diagnostic, trend. Produced by `scripts/scratch/fig_trend_atpk.py` (fc32cbc) from data `sha256:bc6b73d3…`.*

For protein and peptide the orange (observed) bar is the tallest and every gray (relabelled) bar is at or near zero — the observed labelling is special, but with only 8 labellings the strongest statement possible is p = 1/8, above conventional α; for NSAF and PSM even the observed bar is zero.

## Methods / how to produce
Run `scripts/scratch/de_trend_sensitivity.py` (commit fc32cbc) with module `scripts/scratch/analysis/differential_abundance.py` (`method='moderated', trend=True`, seeded from `differential-abundance@0.1` with limma-trend added), on data version `sha256:bc6b73d3…`, environment `pyproject.toml + uv.lock` (Python 3.12.3). All 4 quantities × 3 designs are refit with `eBayes(trend=TRUE)` (prior variance = natural cubic spline, 4 df, of mean log2 abundance fitted to log residual variances); BH per quantity × design over the contrast term; Storey π0 at λ=0.5; the 8 within-pair labelling diagnostic. Fits match R limma 3.58.1 `lmFit + eBayes(trend=TRUE)` to ~1e-12 across all 12 fits (`r_agreement.json`). Figures from `scripts/scratch/fig_trend_atpk.py` (module `scripts/scratch/analysis_figures/trend_atpk.py`). Outputs under `results/de/raloxifene-vs-control/trend/`.

## Discussion
The finding is a near-null framing. A mean-variance trend is a legitimate *a priori* modelling choice — residual variance does fall with abundance here, which is the assumption the trend encodes — but its adoption in this analysis was outcome-driven: it was tried after ATPK was near the top under the plain model. That the "significant" set is one high-abundance protein and ten borderline peptides, all confined to the single pair-blocked LFQ design and none able to clear the relabelling floor, is what keeps the honest conclusion at "no differential abundance." The lesson is about model degrees of freedom under a small n, not about biology: with four pairs the choice between two defensible priors flips the headline count between 0 and 11.

## Caveats
- **Exploratory.** Four pairs, residual df 3. Effects are derived on the same data they were observed in; the forking-paths risk is not retired.
- **The trend model was chosen post hoc**, after ATPK was seen near the top of the no-trend ranking. This is multiplicity on top of the 12 quantity × design fits already run (see the exploration log, 2026-09-24).
- **Run order is aliased with condition** ([finding 0001](0001-run-order-aliased-with-condition.md)): every control ran before its raloxifene-d0 partner, so any hit is condition **or** run order and the two cannot be separated.
- **The 1/8 permutation floor** means no hit here is distinguishable from within-pair relabelling at conventional α; "ranks 1/8" is the best a four-pair permutation can produce.
- **The trend prior is defensible but decisive**: residual SD falls with abundance (Spearman −1.0), so a trend is reasonable — but adopting it, rather than any signal in the data, changes the answer from 0 hits to 11.
- The calibration is unchanged in kind from the primary model: the paired fits still show a small-p excess concentrated at low abundance (see [finding 0004](0004-no-differential-abundance-raloxifene-vs-control.md)).

## Follow-ups
- The ATPK protein hit's LFQ-vs-spectral direction disagreement was pursued as a case study (see Related findings).
- If the trend model is to be reported at all, pre-specify it before analysis on a held-out or orthogonal dataset (confirmatory phase) rather than adopting it after seeing the ranking.

## Related findings
- Refines [finding 0004](0004-no-differential-abundance-raloxifene-vs-control.md): 0004 established the near-null result under a constant prior; this finding shows a trend prior surfaces a handful of hits that still fail the relabelling floor, so the near-null verdict stands.
- Relates to [finding 0001](0001-run-order-aliased-with-condition.md): run order is aliased with condition, so every hit here is condition or run order.
- Relates to [finding 0003](0003-samples-structured-by-matched-pairs.md): the hits appear only in the pair-blocked design, the matched-pair structure documented there.

## References
None.
