"""Stage 4: differential abundance, raloxifene-d0 vs control (moderated linear model).

Runs the project copy of the ``differential-abundance`` template
(``analysis/differential_abundance.py``, ``method="moderated"``) on four quantities:

  * ``protein``  -- FlashLFQ protein LFQ, median-normalized log2 (1,801 complete
                    non-contaminant proteins; the analysis quantity),
  * ``peptide``  -- FlashLFQ peptide LFQ, median-normalized log2 (10,229 peptides),
  * ``nsaf``     -- Limelight NSAF, log2(NSAF + 1e-9) as-is (comparator),
  * ``psm_log2`` -- Limelight PSM counts, log2(count) as-is, no pseudocount (every
                    count >= 1 in the complete set) and no normalization (comparator),

each under three designs (contrast ``condition``, reference ``control``, so a positive
log2FC means higher in raloxifene-d0):

  * ``paired``     -- covariate ``candidate_pair`` (PRIMARY; pair absorbs batch),
  * ``batch``      -- covariate ``batch`` (sensitivity),
  * ``unadjusted`` -- no covariate (sensitivity).

Features constant across all 8 runs (NSAF/PSM rounding) are untestable and are dropped
before testing, with the count reported. It also runs the three protein-level
quantities on their **common feature set** for a like-for-like precision comparison.

Outputs under ``--out-dir`` (default ``results/de/raloxifene-vs-control``):
``<quantity>_<design>.tsv`` (per-feature contrast table),
``results/<quantity>_<design>/`` (the full :class:`DifferentialAbundanceResult`, incl.
covariate-term tests, via ``analysis.result_io.save_result``),
``common_set_comparison.tsv`` and ``summary.json``.
Deterministic (closed-form); no seed consumed.

Run:
    ./.venv/bin/python scripts/scratch/de_raloxifene_vs_control.py
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import sys
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis.differential_abundance import (  # noqa: E402
    DifferentialAbundanceResult,
    Method,
    ZeroResidualVarianceWarning,
    differential_abundance,
)
from analysis.result_io import save_result  # noqa: E402
from common.hashing import sha256_of_file  # noqa: E402
from loaders.data_loading import Dataset  # noqa: E402
from loaders.dataset_io import load_dataset  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "de-raloxifene-vs-control",
    "kind": "analysis",
    "provides": [],
    "uses": [
        "analysis.differential_abundance",
        "analysis.result_io",
        "loaders.dataset_io",
        "common.hashing",
    ],
    "seeded_from": "differential-abundance@0.1 (via analysis/differential_abundance)",
    "description": (
        "Moderated (limma-style) differential abundance raloxifene-d0 vs control on "
        "protein/peptide LFQ (median-normalized log2) and the NSAF / log2 PSM-count "
        "comparators, under paired (candidate_pair, primary), batch and unadjusted "
        "designs; per-feature tables, saved result objects, a common-feature-set "
        "LFQ-vs-spectral comparison and a summary JSON (hits, Storey pi0, p-value "
        "calibration read, design concordance)."
    ),
}

LOG = logging.getLogger("de_raloxifene_vs_control")

CONTRAST = "condition"
REFERENCE = "control"
TREATED = "raloxifene-d0"
METHOD: Method = "moderated"
DESIGNS: dict[str, tuple[str, ...]] = {
    "paired": ("candidate_pair",),
    "batch": ("batch",),
    "unadjusted": (),
}
PRIMARY_DESIGN = "paired"
Q_THRESHOLDS: tuple[float, ...] = (0.05, 0.10)
CONCORDANCE_Q = 0.10
PI0_LAMBDA = 0.5
N_HIST_BINS = 20
# A bin-density deviation beyond this many binomial SDs from uniform is flagged.
CALIBRATION_Z = 3.0
TOP_N = 15


@dataclass(frozen=True)
class QuantitySpec:
    """One quantity to test: its prep-once state and how to present it."""

    state: str  # path under the qc_states root
    annotation_columns: tuple[str, ...]  # feature_metadata columns copied to tables
    log2_counts: bool = False  # True: raw counts -> log2 here (no pseudocount)


QUANTITIES: dict[str, QuantitySpec] = {
    "protein": QuantitySpec("protein/normalized_log", ("accession", "entry")),
    "peptide": QuantitySpec(
        "peptide/normalized_log", ("base_sequence", "protein_groups")
    ),
    "nsaf": QuantitySpec("nsaf/raw_log", ("first_member_id",)),
    "psm_log2": QuantitySpec(
        "psm/raw_linear_complete", ("first_member_id",), log2_counts=True
    ),
}
COMMON_SET_QUANTITIES: tuple[str, ...] = ("protein", "nsaf", "psm_log2")


# --------------------------------------------------------------------------- #
# Loading + validation
# --------------------------------------------------------------------------- #
def validate_design(metadata: pd.DataFrame) -> None:
    """Fail loud unless the sample design is the verified 4-vs-4 paired layout."""
    if not (metadata["sample_role"] == "experimental").all():
        raise ValueError("Non-experimental samples present; DE uses experimental only.")
    counts = metadata[CONTRAST].value_counts().to_dict()
    if counts != {REFERENCE: 4, TREATED: 4}:
        raise ValueError(f"Expected 4 {REFERENCE} + 4 {TREATED}; got {counts}.")
    for pair, rows in metadata.groupby("candidate_pair"):
        if sorted(rows[CONTRAST]) != [REFERENCE, TREATED]:
            raise ValueError(f"Pair {pair!r} is not one {REFERENCE} + one {TREATED}.")
        if len(set(rows["batch"])) != 1:
            raise ValueError(f"Pair {pair!r} spans batches; pair would not nest batch.")


def drop_constant_features(dataset: Dataset) -> tuple[Dataset, list[str]]:
    """Remove features whose value is identical in every sample (untestable)."""
    constant = np.ptp(dataset.abundances, axis=0) == 0
    dropped = [str(x) for x in dataset.feature_names[constant]]
    keep = ~constant
    kept = Dataset(
        abundances=dataset.abundances[:, keep],
        feature_names=dataset.feature_names[keep],
        feature_metadata=dataset.feature_metadata.loc[keep].reset_index(drop=True),
        metadata=dataset.metadata,
        scale=dataset.scale,
    )
    return kept, dropped


def load_quantity(qc_root: Path, spec: QuantitySpec) -> tuple[Dataset, list[str]]:
    """Load one quantity's state, apply its transform, drop constant features."""
    dataset = load_dataset(qc_root / spec.state)
    abundances = np.asarray(dataset.abundances, dtype=float)
    if not np.all(np.isfinite(abundances)):
        raise ValueError(f"{spec.state}: non-finite abundances in a complete state.")
    if spec.log2_counts:
        if dataset.scale != "linear":
            raise ValueError(
                f"{spec.state}: expected linear counts, got {dataset.scale}"
            )
        if abundances.min() < 1 or not np.all(abundances == np.round(abundances)):
            raise ValueError(f"{spec.state}: expected integer counts >= 1.")
        abundances = np.log2(abundances)
        dataset = Dataset(
            abundances=abundances,
            feature_names=dataset.feature_names,
            feature_metadata=dataset.feature_metadata,
            metadata=dataset.metadata,
            scale="log2",
        )
    if dataset.scale != "log2":
        raise ValueError(f"{spec.state}: expected log2 scale, got {dataset.scale}.")
    if "is_contaminant" in dataset.feature_metadata and bool(
        dataset.feature_metadata["is_contaminant"].any()
    ):
        raise ValueError(f"{spec.state}: contaminant features present.")
    validate_design(dataset.metadata)
    return drop_constant_features(dataset)


def subset_features(dataset: Dataset, keep_ids: Sequence[str], key: str) -> Dataset:
    """Restrict ``dataset`` to the features whose ``key`` id is in ``keep_ids``."""
    ids = (
        dataset.feature_names
        if key == "feature"
        else dataset.feature_metadata[key].to_numpy()
    )
    position = {str(x): i for i, x in enumerate(ids)}
    if len(position) != len(ids):
        raise ValueError(f"Duplicate ids in join key {key!r}.")
    idx = np.array([position[k] for k in keep_ids], dtype=int)
    return Dataset(
        abundances=dataset.abundances[:, idx],
        feature_names=np.asarray(list(keep_ids), dtype=object),
        feature_metadata=dataset.feature_metadata.iloc[idx].reset_index(drop=True),
        metadata=dataset.metadata,
        scale=dataset.scale,
    )


# --------------------------------------------------------------------------- #
# Testing + tables
# --------------------------------------------------------------------------- #
def run_design(
    dataset: Dataset, covariates: tuple[str, ...]
) -> DifferentialAbundanceResult:
    """Moderated DE of condition (raloxifene-d0 vs control) with ``covariates``."""
    with warnings.catch_warnings():
        # Recorded on the result (n_prior_floored) and reported in the summary.
        warnings.simplefilter("ignore", ZeroResidualVarianceWarning)
        return differential_abundance(
            dataset,
            CONTRAST,
            covariates=covariates,
            reference={CONTRAST: REFERENCE},
            method=METHOD,
        )


def contrast_table(
    result: DifferentialAbundanceResult, dataset: Dataset, spec: QuantitySpec
) -> pd.DataFrame:
    """Per-feature contrast table (feature, annotations, log2FC, CI, SE, t, p, q...)."""
    ct = result.contrast_table
    if len(result.contrast_terms) != 1:
        raise ValueError(f"Expected one contrast term; got {result.contrast_terms}.")
    annotations = dataset.feature_metadata.loc[:, list(spec.annotation_columns)]
    annotations.insert(0, "feature", dataset.feature_names.astype(str))
    table = ct.rename(
        columns={
            "effect": "log2fc",
            "statistic": "t",
            "sigma": "residual_sd",
            "mean_abundance": "mean_log2_abundance",
        }
    ).loc[
        :,
        [
            "feature",
            "log2fc",
            "ci_low",
            "ci_high",
            "se",
            "t",
            "p",
            "q",
            "residual_sd",
            "mean_log2_abundance",
        ],
    ]
    merged = annotations.merge(table, on="feature", how="inner", validate="1:1")
    if len(merged) != len(table):
        raise ValueError("Annotation join lost features.")
    return merged.sort_values(["p", "feature"], kind="stable").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Diagnostics
# --------------------------------------------------------------------------- #
def storey_pi0(pvalues: np.ndarray, lam: float = PI0_LAMBDA) -> float:
    """Storey's fixed-lambda pi0 = #{p > lam} / (m (1 - lam)), clamped to <= 1.

    Same estimator as ``lib/figures/pvalue_hist.estimate_pi0`` (lambda 0.5).
    """
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        raise ValueError("No finite p-values.")
    return float(min(np.sum(p > lam) / (p.size * (1.0 - lam)), 1.0))


def calibration_read(
    pvalues: np.ndarray, n_bins: int = N_HIST_BINS
) -> dict[str, object]:
    """Qualitative read of a raw-p histogram against Uniform(0, 1).

    Bin densities are scaled so uniform = 1. A region is flagged when its mean density
    departs from 1 by more than ``CALIBRATION_Z`` binomial SDs (the sampling noise of
    a uniform histogram with this many p-values, scaled by the region's bin count):
    the first bin (p < 0.05: spike / deficit), the low quarter (p < 0.25: broad
    excess / deficit) and the top two bins (p > 0.9: hump / dip). The KS distance to
    uniform is descriptive (at these m it rejects almost any departure).
    """
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p)]
    m = p.size
    if m == 0:
        raise ValueError("No finite p-values.")
    counts, _ = np.histogram(p, bins=n_bins, range=(0.0, 1.0))
    density = counts / (m / n_bins)
    sd = float(np.sqrt((1.0 - 1.0 / n_bins) * n_bins / m))
    n_low = n_bins // 4
    regions = (
        ("first bin", float(density[0]), 1, "spike at 0", "deficit near 0"),
        (
            "low quarter",
            float(density[:n_low].mean()),
            n_low,
            "broad excess of small p (p < 0.25)",
            "broad deficit of small p (p < 0.25)",
        ),
        ("top", float(density[-2:].mean()), 2, "hump near 1", "dip near 1"),
    )
    flags: list[str] = []
    for _name, value, n_region, high, low in regions:
        tol = CALIBRATION_Z * sd / np.sqrt(n_region)
        if value > 1.0 + tol:
            flags.append(high)
        elif value < 1.0 - tol:
            flags.append(low)
    read = "approximately uniform" if not flags else "; ".join(flags)
    return {
        "read": read,
        "density_first_bin_p_lt_0.05": round(regions[0][1], 3),
        "density_p_lt_0.25": round(regions[1][1], 3),
        "density_p_gt_0.9": round(regions[2][1], 3),
        "uniform_bin_sd": round(sd, 3),
        "bin_densities": [round(float(d), 3) for d in density],
        "ks_distance_to_uniform": round(float(stats.kstest(p, "uniform").statistic), 4),
    }


def hit_counts(table: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Hits at each q threshold, split by direction (up = higher in raloxifene)."""
    out: dict[str, dict[str, int]] = {}
    for thr in Q_THRESHOLDS:
        hit = table["q"] < thr
        out[f"q<{thr:.2f}"] = {
            "total": int(hit.sum()),
            "up": int((hit & (table["log2fc"] > 0)).sum()),
            "down": int((hit & (table["log2fc"] < 0)).sum()),
        }
    return out


def design_summary(
    result: DifferentialAbundanceResult, table: pd.DataFrame
) -> dict[str, object]:
    """The per-quantity-per-design summary block."""
    p = table["p"].to_numpy(dtype=float)
    ci_half = (table["ci_high"] - table["ci_low"]) / 2.0
    return {
        "covariates": list(result.covariates),
        "n_features_tested": int(np.isfinite(p).sum()),
        "n_samples": result.n_samples,
        "residual_df": result.residual_df,
        "prior_df_d0": _finite_or_str(result.prior_df),
        "prior_variance_s0sq": result.prior_variance,
        "moderated_df_total": _finite_or_str(
            None
            if result.prior_df is None or result.residual_df is None
            else result.residual_df + result.prior_df
        ),
        "n_zero_residual_variance_floored": result.n_prior_floored,
        "hits": hit_counts(table),
        # A global shift (e.g. an unnormalized loading difference) shows up here.
        "median_log2fc": round(float(table["log2fc"].median()), 4),
        "fraction_log2fc_positive": round(float((table["log2fc"] > 0).mean()), 4),
        "storey_pi0_lambda0.5": round(storey_pi0(p), 4),
        "calibration": calibration_read(p),
        "median_residual_sd": round(float(table["residual_sd"].median()), 4),
        "median_moderated_se": round(float(table["se"].median()), 4),
        "median_ci95_half_width": round(float(ci_half.median()), 4),
        "covariate_terms_tested": sorted(
            set(result.table.loc[~result.table["is_contrast"], "term"])
        ),
    }


def _finite_or_str(x: float | None) -> float | str | None:
    """JSON-safe: ``inf`` becomes the string ``"inf"``."""
    if x is None:
        return None
    return float(x) if np.isfinite(x) else "inf"


def hit_set(table: pd.DataFrame, thr: float) -> set[str]:
    return set(table.loc[table["q"] < thr, "feature"].astype(str))


def pairwise_concordance(
    tables: dict[str, pd.DataFrame], key: str = "feature"
) -> list[dict[str, object]]:
    """Spearman of log2FC / t / p and q<0.10 hit-set overlap for every pair."""
    out: list[dict[str, object]] = []
    for a, b in itertools.combinations(tables, 2):
        ta = tables[a].set_index(key)
        tb = tables[b].set_index(key).reindex(ta.index)
        if tb["log2fc"].isna().any():
            raise ValueError(f"{a} and {b} do not share a feature set.")
        sa, sb = hit_set(tables[a], CONCORDANCE_Q), hit_set(tables[b], CONCORDANCE_Q)
        union = sa | sb
        out.append(
            {
                "a": a,
                "b": b,
                "n_features": len(ta),
                "spearman_log2fc": _spearman(ta["log2fc"], tb["log2fc"]),
                "max_abs_log2fc_difference": float(
                    np.max(np.abs(ta["log2fc"] - tb["log2fc"]))
                ),
                "spearman_t": _spearman(ta["t"], tb["t"]),
                "spearman_p": _spearman(ta["p"], tb["p"]),
                f"hits_q<{CONCORDANCE_Q:.2f}": {
                    "a": len(sa),
                    "b": len(sb),
                    "both": len(sa & sb),
                    "jaccard": (len(sa & sb) / len(union)) if union else None,
                },
            }
        )
    return out


def _spearman(x: pd.Series, y: pd.Series) -> float:
    return round(float(stats.spearmanr(x.to_numpy(), y.to_numpy()).statistic), 4)


def pair_relabel_diagnostic(dataset: Dataset) -> dict[str, object]:
    """Paired design under every within-pair relabelling (the restricted null).

    Swapping control/raloxifene inside a subset of pairs is a valid null relabelling
    for a paired design; fixing the first pair's orientation enumerates each distinct
    |t| configuration once (a global swap only flips signs): 2^(n_pairs-1)
    labellings, the observed one included. If the moderated paired test is calibrated
    and there is no condition effect, the observed labelling's p-value distribution
    should look like the others'. Non-observed labellings also break the run-order
    aliasing (control ran first in every pair; finding 0001), so an observed
    labelling that stands out reflects condition *or* run order, not the model.
    """
    md = dataset.metadata
    pairs = sorted(md["candidate_pair"].unique())
    swap = {REFERENCE: TREATED, TREATED: REFERENCE}
    rows: list[dict[str, object]] = []
    for flips in itertools.product((False, True), repeat=len(pairs) - 1):
        flipped = [pr for pr, f in zip(pairs[1:], flips, strict=True) if f]
        in_flip = md["candidate_pair"].isin(flipped)
        relabelled = md[CONTRAST].where(~in_flip, md[CONTRAST].map(swap))
        ds = Dataset(
            abundances=dataset.abundances,
            feature_names=dataset.feature_names,
            feature_metadata=dataset.feature_metadata,
            metadata=md.assign(**{CONTRAST: relabelled}),
            scale=dataset.scale,
        )
        result = run_design(ds, DESIGNS[PRIMARY_DESIGN])
        p = result.contrast_table["p"].to_numpy(dtype=float)
        q = result.contrast_table["q"].to_numpy(dtype=float)
        cal = calibration_read(p)
        rows.append(
            {
                "flipped_pairs": flipped,
                "observed": not flipped,
                "pi0": round(storey_pi0(p), 4),
                "density_p_lt_0.25": cal["density_p_lt_0.25"],
                "density_first_bin_p_lt_0.05": cal["density_first_bin_p_lt_0.05"],
                "hits_q<0.10": int(np.sum(q < 0.10)),
                "prior_df_d0": _finite_or_str(result.prior_df),
            }
        )
    observed = next(r for r in rows if r["observed"])
    pi0s = [float(str(r["pi0"])) for r in rows]
    return {
        "n_labellings": len(rows),
        "observed_pi0_rank_lowest_first": 1
        + sum(v < float(str(observed["pi0"])) for v in pi0s),
        "note": "rank 1 of n = observed has the most small-p excess; the smallest "
        "attainable restricted-permutation p is 1/n_labellings",
        "labellings": rows,
    }


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run_quantity(
    name: str,
    dataset: Dataset,
    spec: QuantitySpec,
    out_dir: Path,
    prefix: str = "",
) -> tuple[dict[str, pd.DataFrame], dict[str, object]]:
    """All designs for one quantity: write tables + result objects; return tables."""
    tables: dict[str, pd.DataFrame] = {}
    summaries: dict[str, object] = {}
    for design, covariates in DESIGNS.items():
        result = run_design(dataset, covariates)
        table = contrast_table(result, dataset, spec)
        stem = f"{prefix}{name}_{design}"
        if not prefix:
            table.to_csv(out_dir / f"{stem}.tsv", sep="\t", index=False)
            save_result(result, out_dir / "results" / stem)
        tables[design] = table
        summaries[design] = design_summary(result, table)
        LOG.info(
            "%s %s: n=%d df=%s d0=%s hits q<0.05=%d q<0.10=%d",
            stem,
            covariates,
            len(table),
            result.residual_df,
            result.prior_df,
            int((table["q"] < 0.05).sum()),
            int((table["q"] < 0.10).sum()),
        )
    return tables, summaries


def top_proteins(table: pd.DataFrame, n: int = TOP_N) -> list[dict[str, object]]:
    rows = table.sort_values(["q", "p"], kind="stable").head(n)
    return [
        {
            "accession": str(r.accession),
            "entry": str(r.entry),
            "log2fc": round(float(r.log2fc), 4),
            "ci95": [round(float(r.ci_low), 4), round(float(r.ci_high), 4)],
            "p": float(r.p),
            "q": round(float(r.q), 4),
        }
        for r in rows.itertuples()
    ]


def common_set(
    datasets: dict[str, Dataset],
) -> tuple[list[str], dict[str, Dataset]]:
    """Protein ids present (and non-constant) in every common-set quantity."""
    lfq_ids = [str(x) for x in datasets["protein"].feature_names]
    spectral = [
        {str(x) for x in datasets[q].feature_metadata["first_member_id"]}
        for q in COMMON_SET_QUANTITIES
        if q != "protein"
    ]
    common = [i for i in lfq_ids if all(i in s for s in spectral)]
    subsets = {
        q: subset_features(
            datasets[q], common, "feature" if q == "protein" else "first_member_id"
        )
        for q in COMMON_SET_QUANTITIES
    }
    # Carry LFQ protein annotations onto every subset so tables share keys.
    lfq_meta = subsets["protein"].feature_metadata.loc[:, ["accession", "entry"]]
    for q in COMMON_SET_QUANTITIES:
        if q != "protein":
            subsets[q].feature_metadata = subsets[q].feature_metadata.assign(
                accession=lfq_meta["accession"].to_numpy(),
                entry=lfq_meta["entry"].to_numpy(),
            )
    return common, subsets


def common_set_table(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Wide per-feature table: one block of columns per quantity (primary design)."""
    base = tables["protein"].loc[:, ["feature", "accession", "entry"]]
    cols = ["log2fc", "se", "ci_low", "ci_high", "t", "p", "q", "residual_sd"]
    for q, t in tables.items():
        block = t.loc[:, ["feature", *cols]].rename(
            columns={c: f"{q}_{c}" for c in cols}
        )
        base = base.merge(block, on="feature", how="inner", validate="1:1")
    return base.sort_values("protein_p", kind="stable").reset_index(drop=True)


def main(argv: Sequence[str] | None = None) -> int:
    root = _SCRATCH.parent.parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--qc-root", type=Path, default=root / "results" / "qc_states")
    parser.add_argument(
        "--out-dir", type=Path, default=root / "results/de/raloxifene-vs-control"
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.qc_root / "manifest.json").read_text(encoding="utf-8"))
    LOG.info("data_version %s; qc_root %s", manifest["data_version"], args.qc_root)

    datasets: dict[str, Dataset] = {}
    dropped: dict[str, list[str]] = {}
    per_quantity: dict[str, object] = {}
    concordance: dict[str, object] = {}
    all_tables: dict[str, dict[str, pd.DataFrame]] = {}
    for name, spec in QUANTITIES.items():
        dataset, dropped[name] = load_quantity(args.qc_root, spec)
        if datasets and not dataset.metadata.index.equals(
            next(iter(datasets.values())).metadata.index
        ):
            raise ValueError(f"{name}: sample order differs from other quantities.")
        datasets[name] = dataset
        LOG.info(
            "%s: %d features (%d constant dropped)",
            name,
            dataset.abundances.shape[1],
            len(dropped[name]),
        )
        tables, summaries = run_quantity(name, dataset, spec, out_dir)
        all_tables[name] = tables
        per_quantity[name] = {
            "state": spec.state,
            "transform": "log2(count), no pseudocount (all counts >= 1)"
            if spec.log2_counts
            else "as stored (log2)",
            "n_constant_features_dropped": len(dropped[name]),
            "constant_features_dropped": dropped[name],
            "designs": summaries,
            "pair_relabel_diagnostic": pair_relabel_diagnostic(dataset),
        }
        concordance[name] = pairwise_concordance(tables)

    # Common feature set: LFQ protein vs NSAF vs log2 PSM, like for like.
    common_ids, subsets = common_set(datasets)
    common_tables: dict[str, dict[str, pd.DataFrame]] = {}
    common_summary: dict[str, object] = {}
    for q in COMMON_SET_QUANTITIES:
        spec = QuantitySpec(QUANTITIES[q].state, ("accession", "entry"))
        tables, summaries = run_quantity(q, subsets[q], spec, out_dir, prefix="common_")
        common_tables[q] = tables
        common_summary[q] = summaries
    by_design = {
        d: pairwise_concordance({q: common_tables[q][d] for q in COMMON_SET_QUANTITIES})
        for d in DESIGNS
    }
    common_set_table(
        {q: common_tables[q][PRIMARY_DESIGN] for q in COMMON_SET_QUANTITIES}
    ).to_csv(out_dir / "common_set_comparison.tsv", sep="\t", index=False)

    summary = {
        "question": "Differential abundance raloxifene-d0 vs control (positive log2FC "
        "= higher in raloxifene-d0).",
        "data_version": manifest["data_version"],
        "script": "scripts/scratch/de_raloxifene_vs_control.py",
        "script_sha256": sha256_of_file(Path(__file__)),
        "seeded_from": "differential-abundance@0.1 (lib/analysis/"
        "differential_abundance.py), project copy scripts/scratch/analysis/",
        "method": METHOD,
        "contrast": {"column": CONTRAST, "reference": REFERENCE, "level": TREATED},
        "designs": {d: list(c) for d, c in DESIGNS.items()},
        "primary_design": PRIMARY_DESIGN,
        "correction": "Benjamini-Hochberg FDR, per quantity x design over the features "
        "tested (contrast term); covariate terms BH-adjusted separately per term",
        "ci": "95% t interval on the moderated SE (df = residual_df + d0)",
        "sample_set": "experimental (8 of 8; no QC/pool controls exist)",
        "samples": list(datasets["protein"].metadata.index.astype(str)),
        "pi0_estimator": f"Storey, fixed lambda = {PI0_LAMBDA}",
        "calibration_rule": f"{N_HIST_BINS} bins; flag if the first bin or the p>0.9 "
        f"mean departs from uniform by > {CALIBRATION_Z} binomial SDs",
        "quantities": per_quantity,
        "design_concordance": concordance,
        "top_proteins_primary": top_proteins(all_tables["protein"][PRIMARY_DESIGN]),
        "common_set": {
            "definition": "LFQ-complete proteins (1,801) whose protein-quants id "
            "equals a spectral-complete group's first_member_id, excluding groups "
            "constant in NSAF or PSM counts",
            "n_features": len(common_ids),
            "quantities": common_summary,
            "quantity_concordance": by_design,
            "table": "common_set_comparison.tsv (primary design)",
        },
        "result_objects": "results/<quantity>_<design>/ -- load with "
        "analysis.result_io.load_result(path, DifferentialAbundanceResult)",
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    LOG.info("wrote %s", out_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
