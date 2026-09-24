---
id: 9
title: "Small, consistent modification-class shifts with treatment (CAM/Cys peptides lower, oxidized forms higher), all confounded with run order; oxidation evidence is mixed"
status: candidate
phase: exploratory
kind: discovery
created: 2026-09-24
updated: 2026-09-24

summary: "In the complete-case paired peptide set (10,229 peptides), carbamidomethyl-only (Cys) peptides shift slightly lower with treatment (log2FC −0.066 [−0.101, −0.032], 4/4 pairs, q 0.009) and oxidized peptides slightly higher (+0.073 [0.040, 0.106], 4/4, q 0.009), with matching detection-asymmetry ORs (CAM 0.73, oxidized 1.42). All effects are small (<0.12 log2, <9%) and every one is fully confounded with run order (control was run first in every pair). The CAM-only shift survives protein adjustment (−0.056, p 0.0077); the oxidation shift does NOT (+0.076, p 0.072) and its global index (p 0.55 complete / 0.10 missing-as-zero) and per-protein tests (0/64 at q<0.10) are non-significant."
verdict: "Small, coherent, treatment-associated shifts in modification-class abundance and detection, all of them confounded with run order — mechanism undetermined, and explicitly NOT to be called 'drift'. The CAM-only (Cys-peptide) decrease is the firmer of the two: it survives protein-level adjustment. The oxidation-increase signal is mixed and must be reported as such — the +0.073 class shift does not survive protein adjustment (p 0.072), the global oxidation index is non-significant (p 0.55 / 0.10), and no protein passes per-protein testing (0/64 at q<0.10). Exploratory (4 pairs); a global BH across the many pair-level class tests is still owed."

entities: []

relationships:
  - { type: relates_to, target: 1, note: "Every class/detection shift is fully confounded with run order (control run first in every pair); the run-order aliasing of 0001 is why none of these can be attributed to treatment mechanism." }
  - { type: relates_to, target: 3, note: "The analysis is paired on candidate pair (0003); the class shifts are per-pair paired differences." }
  - { type: relates_to, target: 4, note: "Same complete-case paired peptide set (10,229) and paired design as the near-null differential-abundance analysis; these are small within-peptide class effects beneath that null." }

provenance:
  data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
  script: { path: "scripts/scratch/peptide_mod_analysis.py", commit: "fc32cbc" }
  params:
    module: "scripts/scratch/analysis/peptide_mods.py"
    sample_set: "experimental (all 8 runs; no QC/pool controls exist)"
    n_experimental: 8
    n_controls_excluded: 0
    pairs: ["P905_906", "P907_908", "P909_910", "P941_942"]
    contrast: "raloxifene-d0 - control (positive = higher in raloxifene-d0)"
    quantity: "LFQ peptide, median-normalized log2, 10,229 complete features"
    classes: ["unmodified", "CAM-only", "oxidized", "raloxifene-adduct"]
    effects:
      class_log2fc: "paired log2FC (ralox - control) per class; class − unmodified pair-level t (df 3)"
      protein_adjusted: "peptide − protein paired log2FC, pair-median centred"
      modified_vs_reference: "mean over forms of (modified − reference) paired log2FC"
      oxidation_index: "per-pair log2(global oxidized/unoxidized index treated / control); complete-case and missing-as-zero; per-protein"
    correction: "BH within small families (per test); a global BH across all pair-level class tests is owed"
    outputs: "results/peptide-mods/{summary.json, class_log2fc.tsv, class_log2fc_tests.tsv, protein_adjusted.tsv, protein_adjusted_class_tests.tsv, modified_vs_reference*.tsv, oxidation_index_sample.tsv, oxidation_index_protein.tsv, discordance_class_tests.tsv, discordance_by_pair.tsv}"
  environment: "pyproject.toml + uv.lock (Python 3.12.3; numpy 2.5.3, pandas 3.0.6)"
  seeded_from: null
  seed: null
  result_id: null

evidence:
  - metric: "CAM-only − unmodified paired log2FC"
    value: -0.066
    ci: [-0.101, -0.032]
    p_value: 0.0087
    p_adjusted: 0.0087
    correction: "BH within class family"
    test: "pair-level paired t on (class − unmodified) mean shift, df 3"
    n: 8
    note: "4/4 pairs negative (−0.045 / −0.068 / −0.057 / −0.096). Firmer of the two class shifts: survives protein adjustment."
  - metric: "oxidized − unmodified paired log2FC"
    value: 0.073
    ci: [0.040, 0.106]
    p_value: 0.0058
    p_adjusted: 0.0087
    correction: "BH within class family"
    test: "pair-level paired t, df 3"
    n: 8
    note: "4/4 pairs positive (0.066 / 0.048 / 0.083 / 0.096). MUST be reported alongside its non-significant related tests (below) — oxidation evidence is mixed."
  - metric: "oxidized-form − unoxidized-form paired log2FC (matched forms)"
    value: 0.121
    ci: [0.024, 0.218]
    p_value: 0.029
    p_adjusted: 0.029
    correction: "BH"
    test: "pair-level paired t over 123 matched form pairs, df 3"
    n: 8
    note: "123 base sequences with both an oxidized and an unoxidized form quantified; 4/4 pairs positive."
  - metric: "detection asymmetry by class — pair-level OR vs unmodified"
    value: "CAM-only 0.73; oxidized 1.42"
    ci: null
    p_value: "CAM 0.027; oxidized 0.012"
    p_adjusted: null
    correction: "BH within class family"
    test: "pair-stratified OR on discordant detection cells; pair-level sign"
    n: 8
    note: "CAM-only 4/4 pairs OR<1 (fewer CAM detections in treated); oxidized 4/4 pairs OR>1 (more oxidized detections in treated). Directionally matches the intensity shifts."
  - metric: "protein-adjusted CAM-only shift (survives)"
    value: -0.056
    ci: [-0.069, -0.042]
    p_value: 0.0077
    p_adjusted: 0.015
    correction: "BH"
    test: "peptide − protein paired log2FC, pair-median centred; pair-level t df 3"
    n: 8
    note: "The CAM-only decrease persists after removing each peptide's own protein-level fold change → not merely a protein-abundance artifact."
  - metric: "protein-adjusted oxidation shift (does NOT survive)"
    value: 0.076
    ci: [-0.012, 0.164]
    p_value: 0.072
    p_adjusted: 0.072
    correction: "BH"
    test: "peptide − protein paired log2FC, pair-median centred; pair-level t df 3"
    n: 8
    note: "Non-significant after protein adjustment (pair-level CI now spans 0). This is one of the three mixed-evidence results the oxidation shift must be reported with."
  - metric: "global oxidation index, treated/control"
    value: "+0.038 (complete-case); +0.096 (missing-as-zero)"
    ci: null
    p_value: "0.55 (complete); 0.10 (missing-as-zero)"
    p_adjusted: null
    correction: null
    test: "paired t on per-pair log2(oxidized/unoxidized index treated/control), df 3"
    n: 8
    note: "Both non-significant. Second mixed-evidence result for oxidation."
  - metric: "per-protein oxidation index tests"
    value: "0 of 64 proteins at q < 0.10"
    ci: null
    p_value: null
    p_adjusted: null
    correction: "BH over 64 proteins"
    test: "per-protein paired oxidation index"
    n: 8
    note: "min q 0.447. Third mixed-evidence result: no individual protein carries the oxidation signal."
  - metric: "baseline unmodified-peptide detection, treated-only share by pair"
    value: "0.41 / 0.29 / 0.36 (2021 pairs); 0.52 (2022 pair)"
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "treated-only share of discordant unmodified-peptide detection cells, per pair"
    n: 8
    note: "The three 2021-05-06 pairs detect FEWER unmodified peptides in the treated (later-run) member than in the paired control; the single 2022-03-18 pair is near 0.5. Consistent with a run-position detection effect, and confounded with it."

figures:
  - png: "figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-by-pair.png"
    svg: "figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-by-pair.svg"
    legend_png: "figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-by-pair.legend.png"
    legend_svg: "figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-by-pair.legend.svg"
    caption: "Per-pair paired class shifts (ralox − control), mean ± 95% CI (t, df 3): CAM-only − unmodified −0.066 [−0.101, −0.032] (4/4 pairs negative), oxidized − unmodified +0.073 [0.040, 0.106] (4/4 positive), oxidized-form − unoxidized-form +0.121 [0.024, 0.218] (4/4, 123 forms); per-pair points colored by candidate pair; dashed line at 0. n=8."
    script: { path: "scripts/scratch/fig_peptide_mods.py", commit: "fc32cbc" }
    data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
    result_id: null
    params: { module: "scripts/scratch/analysis_figures/peptide_mods_figs.py", input: "results/peptide-mods/{class_log2fc_tests.tsv, modified_vs_reference_class_tests.tsv}", provenance: "results/peptide-mods/figure_provenance.json", dpi: 300 }
  - png: "figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-distributions.png"
    svg: "figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-distributions.svg"
    legend_png: "figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-distributions.legend.png"
    legend_svg: "figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-distributions.legend.svg"
    caption: "Kernel-density distributions of per-peptide paired log2FC (ralox − control) by modification class: unmodified (n=8,284, median 0.017), CAM-only (n=1,801, median −0.049), oxidized (n=144, median 0.092); view ±1.25 log2. The class medians differ by <0.12 log2 while each distribution's spread is ~0.25 SD — the shifts are small relative to between-peptide variability. n=8."
    script: { path: "scripts/scratch/fig_peptide_mods.py", commit: "fc32cbc" }
    data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
    result_id: null
    params: { module: "scripts/scratch/analysis_figures/peptide_mods_figs.py", input: "results/peptide-mods/class_log2fc.tsv", kde: "gaussian_kde Scott bandwidth", provenance: "results/peptide-mods/figure_provenance.json", dpi: 300 }
  - png: "figures/analysis/peptide-mods/class-shifts/0009-detection-by-run-position.png"
    svg: "figures/analysis/peptide-mods/class-shifts/0009-detection-by-run-position.svg"
    legend_png: "figures/analysis/peptide-mods/class-shifts/0009-detection-by-run-position.legend.png"
    legend_svg: "figures/analysis/peptide-mods/class-shifts/0009-detection-by-run-position.legend.svg"
    caption: "Treated-only share of discordant unmodified-peptide detection cells, per pair, with 95% CI: P905_906 0.41, P907_908 0.29, P909_910 0.36 (all three 2021-05-06 pairs below 0.5 → fewer detections in the later-run treated member), P941_942 0.52 (the 2022-03-18 pair); dashed line at 0.5. n=8."
    script: { path: "scripts/scratch/fig_peptide_mods.py", commit: "fc32cbc" }
    data_version: "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
    result_id: null
    params: { module: "scripts/scratch/analysis_figures/peptide_mods_figs.py", input: "results/peptide-mods/discordance_by_pair.tsv", provenance: "results/peptide-mods/figure_provenance.json", dpi: 300 }

references: []

validation:
  computational_reproduction: { status: not_attempted }
  analytic_replication:       { status: not_attempted }
  data_replication:           { status: not_attempted }

integrity_signoff: true
---

# Small, consistent modification-class shifts with treatment (CAM/Cys peptides lower, oxidized forms higher), all confounded with run order; oxidation evidence is mixed

## Summary
In the complete-case paired peptide set (10,229 peptides), two small modification-class shifts move consistently with treatment: carbamidomethyl-only (Cys) peptides are slightly lower in the raloxifene member of each pair (paired log2FC −0.066, 95% CI −0.101 to −0.032, 4/4 pairs, q 0.009) and oxidized peptides slightly higher (+0.073, CI 0.040 to 0.106, 4/4 pairs, q 0.009). The class-level detection asymmetry runs the same way (pair-level OR vs unmodified 0.73 for CAM-only, 1.42 for oxidized). Every effect is small (<0.12 log2, <9%) and every one is fully confounded with run order: control was run first in every pair. The CAM-only decrease survives protein-level adjustment (−0.056, p 0.0077); the oxidation increase does not (+0.076, p 0.072), and the global oxidation index and per-protein tests are non-significant.

## Verdict
These are small, coherent, treatment-associated shifts in both modification-class abundance and detection — but all of them are confounded with run order, so the mechanism is undetermined and must *not* be called "drift". The CAM-only (Cys-peptide) decrease is the firmer of the two, because it persists after removing each peptide's own protein-level fold change. The oxidation-increase signal is genuinely mixed and has to be reported that way: the +0.073 class shift does not survive protein adjustment (p 0.072), the global oxidation index is non-significant (p 0.55 complete-case, 0.10 missing-as-zero), and 0 of 64 proteins pass per-protein testing (q<0.10). Recorded as an exploratory candidate (4 pairs); a global BH across the many pair-level class tests is still owed. The CAM-only shift is not attributable to the observed Cys–raloxifene adducts of [finding 0008](0008-raloxifene-cys-adducts-treatment-specific.md): there are only 39 adduct features at a ≈2.5% fraction.

## Evidence

**Two small class shifts move consistently across all four pairs.** The CAM-only − unmodified paired shift is −0.066 (CI −0.101 to −0.032), negative in all four pairs (−0.045 / −0.068 / −0.057 / −0.096); the oxidized − unmodified shift is +0.073 (CI 0.040 to 0.106), positive in all four (0.066 / 0.048 / 0.083 / 0.096); and among the 123 base sequences with both an oxidized and an unoxidized form quantified, the oxidized-form − unoxidized-form shift is +0.121 (CI 0.024 to 0.218), again 4/4 pairs.

![Per-pair paired class shifts, mean with 95% CI: CAM-only minus unmodified near −0.066 (4/4 pairs below zero), oxidized minus unmodified near +0.073 (4/4 above zero), oxidized-form minus unoxidized-form near +0.121 (4/4 above zero); per-pair points colored by pair; dashed line at zero; n=8.](../figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-by-pair.png)

![Legend for Figure 1 — series: CAM-only − unmodified, oxidized − unmodified, oxidized-form − unoxidized-form; large marker = pair-level mean, whiskers = 95% t CI (df 3), small markers = the four candidate pairs; dashed line at 0.](../figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-by-pair.legend.png)

*Figure 1. Per-pair class shifts. Produced by `scripts/scratch/fig_peptide_mods.py` (module `scripts/scratch/analysis_figures/peptide_mods_figs.py`, commit fc32cbc) from `results/peptide-mods/class_log2fc_tests.tsv` and `modified_vs_reference_class_tests.tsv`, data `sha256:bc6b73d3…1ba74`.*

Each row is one class contrast; the large marker is the pair-level mean and the four small markers are the individual pairs. The CAM-only mean sits left of the dashed zero line with all four pairs on the same (negative) side; the two oxidation contrasts sit right of zero with all four pairs positive. The consistency across pairs is what makes these effects real as *associations* — but note that all four pairs share the same control-first run order, which Figure 3 and the caveats return to.

**The shifts are small relative to between-peptide spread.** The class medians differ by less than 0.12 log2 (under 9%), while each class's per-peptide log2FC distribution has a spread of roughly 0.25 SD.

![Kernel-density distributions of per-peptide paired log2FC by class — unmodified centered near 0.017, CAM-only shifted slightly negative (median −0.049), oxidized shifted slightly positive (median 0.092); heavily overlapping; view plus or minus 1.25 log2; n=8.](../figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-distributions.png)

![Legend for Figure 2 — density curves per modification class (unmodified, CAM-only, oxidized) with class medians marked; x = paired log2FC (ralox − control).](../figures/analysis/peptide-mods/class-shifts/0009-class-log2fc-distributions.legend.png)

*Figure 2. Class log2FC distributions. Produced by `scripts/scratch/fig_peptide_mods.py` (module `scripts/scratch/analysis_figures/peptide_mods_figs.py`, commit fc32cbc) from `results/peptide-mods/class_log2fc.tsv`, data `sha256:bc6b73d3…1ba74`.*

The three density curves sit almost on top of one another, centered near zero; the CAM-only curve's median is nudged left and the oxidized curve's median nudged right, but the offsets are tiny next to the width of each curve. The message is that these are small mean displacements of broad, overlapping distributions — detectable because they are consistent across pairs, not because any individual peptide moves a lot.

**The detection asymmetry runs the same direction as the intensity shifts, and the unmodified baseline itself depends on run position.** By class, the pair-level detection OR versus unmodified is 0.73 for CAM-only (4/4 pairs OR<1, p 0.027) and 1.42 for oxidized (4/4 pairs OR>1, p 0.012): fewer CAM detections and more oxidized detections in the treated member. Underneath that, the unmodified-peptide detection baseline is itself run-position-dependent: in the three 2021-05-06 pairs the treated-only share of discordant unmodified cells is 0.41 / 0.29 / 0.36 (the later-run treated member detects *fewer* peptides than its control), while the single 2022-03-18 pair is 0.52.

![Treated-only share of discordant unmodified-peptide detection cells per pair, with 95% CI — the three 2021-05-06 pairs at 0.41, 0.29 and 0.36 all below 0.5, the 2022-03-18 pair at 0.52; dashed line at 0.5; n=8.](../figures/analysis/peptide-mods/class-shifts/0009-detection-by-run-position.png)

![Legend for Figure 3 — per-pair treated-only detection share (unmodified peptides) with Clopper-Pearson 95% CI; pairs colored/labeled by acquisition batch (2021-05-06 vs 2022-03-18); dashed reference at 0.5.](../figures/analysis/peptide-mods/class-shifts/0009-detection-by-run-position.legend.png)

*Figure 3. Baseline detection by run position. Produced by `scripts/scratch/fig_peptide_mods.py` (module `scripts/scratch/analysis_figures/peptide_mods_figs.py`, commit fc32cbc) from `results/peptide-mods/discordance_by_pair.tsv`, data `sha256:bc6b73d3…1ba74`.*

Three of the four pairs sit below the dashed 0.5 line, meaning the treated (later-acquired) run identifies fewer unmodified peptides than its paired control. Because condition and run order are aliased ([finding 0001](0001-run-order-aliased-with-condition.md)), this baseline shift cannot be separated from treatment, and it is the backdrop against which the small class shifts must be read.

**The CAM-only shift survives protein adjustment; the oxidation shift does not, and its other tests are non-significant.** Removing each peptide's own protein-level fold change (pair-median centred), the CAM-only decrease persists (−0.056, p 0.0077), so it is not merely a protein-abundance artifact. The oxidation increase, by contrast, does *not* survive (+0.076, p 0.072); the global oxidized/unoxidized index is non-significant both complete-case (+0.038, p 0.55) and missing-as-zero (+0.096, p 0.10); and no protein passes per-protein testing (0 of 64 at q<0.10, min q 0.447). The oxidation evidence is therefore mixed and is reported here as the full set, not the significant subset alone.

## Methods / how to produce
Run `scripts/scratch/peptide_mod_analysis.py` (module `scripts/scratch/analysis/peptide_mods.py`, written from scratch) at commit fc32cbc on data `sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74`, under `pyproject.toml` + `uv.lock` (Python 3.12.3, numpy 2.5.3, pandas 3.0.6). The quantity is LFQ peptide (median-normalized log2, 10,229 complete-case features), paired within the four candidate pairs (P905_906, P907_908, P909_910, P941_942). Per class, the per-peptide paired log2FC (ralox − control) is summarized and the (class − unmodified) mean shift is tested at the pair level with a paired t (df 3), BH-corrected within the class family. Protein adjustment subtracts each peptide's protein-level paired log2FC (pair-median centred) before the same test. The matched-form contrast (oxidized − unoxidized) uses the 123 base sequences with both forms quantified. The global oxidation index is the per-pair log2 ratio of summed oxidized to unoxidized intensity, tested complete-case and missing-as-zero, and per protein (BH over 64). Detection asymmetry by class uses a pair-stratified OR on discordant detection cells. Outputs are under `results/peptide-mods/` (`class_log2fc*.tsv`, `protein_adjusted*.tsv`, `modified_vs_reference*.tsv`, `oxidation_index_*.tsv`, `discordance_*`). Figures are from `scripts/scratch/fig_peptide_mods.py` (module `scripts/scratch/analysis_figures/peptide_mods_figs.py`); provenance sidecar `results/peptide-mods/figure_provenance.json`. The statistics review passed with mandatory wording: the shifts must not be called "drift"; the oxidation shift must be reported together with its non-significant protein-adjusted, global-index and per-protein results.

## Discussion
These are the kind of small, systematic effects that a modification-aware peptide analysis surfaces beneath a protein-abundance null ([finding 0004](0004-no-differential-abundance-raloxifene-vs-control.md)): the class means move only a few percent, but they move the same way in every pair, which is why they clear a pair-level test. The direction is internally coherent — CAM-only (Cys) peptides both decrease in intensity and are detected less often in the treated member, while oxidized peptides both increase and are detected more often.

What cannot be said is *why*. Because control was run first in every pair, run order is perfectly aliased with condition ([finding 0001](0001-run-order-aliased-with-condition.md)), and the run-position dependence of the unmodified detection baseline (Figure 3) shows the acquisition sequence alone moves these readouts. The mechanism is therefore undetermined, and this must not be described as "drift" — that names a specific cause the design cannot support. It is tempting to link the CAM-only decrease to the raloxifene Cys adducts of [finding 0008](0008-raloxifene-cys-adducts-treatment-specific.md) (an adducted cysteine cannot also be carbamidomethylated), but the arithmetic does not support that: there are only 39 adduct features at a ≈2.5% median fraction, far too few to move a class of 1,801 CAM-only peptides by 0.06 log2. The two observations are directionally suggestive but not causally connected by this evidence.

## Caveats
1. **Exploratory, minimal design.** 4 matched pairs, residual df 3; all pair-level tests rest on 4 values.
2. **All shifts are confounded with run order; mechanism undetermined — do NOT call it "drift".** Control was run first in every pair ([finding 0001](0001-run-order-aliased-with-condition.md)), so condition and run position are aliased. The unmodified detection baseline is itself run-position-dependent (Figure 3). No mechanism (drift, column state, sample age, treatment) can be singled out.
3. **Oxidation evidence is mixed and must be reported in full.** The +0.073 class shift is significant, but the protein-adjusted version is not (p 0.072), the global oxidation index is non-significant (p 0.55 complete-case, 0.10 missing-as-zero), and 0 of 64 proteins pass per-protein testing (q<0.10). Do not headline only the significant subset.
4. **The CAM-only shift is firmer but not attributable to the observed Cys adducts.** It survives protein adjustment (−0.056, p 0.0077), unlike the oxidation shift. It is not explained by the 39 raloxifene Cys adducts of [finding 0008](0008-raloxifene-cys-adducts-treatment-specific.md) (≈2.5% fraction is far too small).
5. **Effects are small relative to spread.** All class means differ by <0.12 log2 (<9%), against a per-peptide log2FC spread of ~0.25 SD (Figure 2).
6. **Multiplicity — a global BH is owed.** BH was applied only within small class families (per test). A global BH across all the pair-level class tests run here has not been applied and is still owed.
7. **The 2022 pair shows the largest shifts.** P941_942 (the single 2022-03-18 pair, run back-to-back) has the largest per-pair CAM-only (−0.096) and oxidation (0.096) shifts. That weakly argues against a mechanism proportional to the inter-run time gap, but it is n=1 and cannot be leaned on.

## Follow-ups
- Apply a global BH across all pair-level class tests and report which survive.
- Ask whether a run-order-only (non-treatment) model can reproduce the class shifts and the baseline detection pattern of Figure 3.
- Have a blinded verifier re-derive the CAM-only shift and its protein-adjusted version under a pre-specified concordance criterion.
- Promote the analysis and figure scripts to `scripts/promoted/` before any validation attempt.

## Related findings
- [Finding 0001 (run order aliased with condition)](0001-run-order-aliased-with-condition.md) — `relates_to`. Control was run first in every pair, so all class and detection shifts are confounded with run position; this is why the mechanism is undetermined.
- [Finding 0003 (samples structured by matched pairs)](0003-samples-structured-by-matched-pairs.md) — `relates_to`. The analysis is paired on candidate pair, and the shifts are per-pair paired differences.
- [Finding 0004 (no differential abundance)](0004-no-differential-abundance-raloxifene-vs-control.md) — `relates_to`. Same complete-case paired peptide set (10,229) and design; these are small within-peptide class effects sitting beneath that protein-level null.
- [Finding 0008 (raloxifene Cys adducts, treatment-specific)](0008-raloxifene-cys-adducts-treatment-specific.md) — mentioned as a candidate but rejected explanation for the CAM-only shift (too few adduct features). No edge asserted; surfaced for the scientist to type if wanted.

## References
None. All statements are about this dataset or are general statistical background (blocking, BH, run-order aliasing); no external-knowledge claim is made here.
</content>
