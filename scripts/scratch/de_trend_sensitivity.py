"""Stage 4: limma-trend sensitivity of the raloxifene-d0 vs control DE (all quantities).

Re-runs the moderated DE of ``de_raloxifene_vs_control.py`` with the empirical-Bayes
prior variance trending with mean log-abundance (``differential_abundance(...,
trend=True)``, limma ``eBayes(trend=TRUE)``; validated against R limma 3.58.1 by
``validate_limma_trend.py``) on the **same inputs**: the four quantities (LFQ protein,
LFQ peptide, NSAF, log2 PSM), the same states/transforms/constant-feature drops (the
base runner's ``load_quantity``) and the same three designs (paired = primary, batch,
unadjusted). Each trend fit is set beside the no-trend fit of the same data; the
no-trend fit is recomputed here and checked against the on-disk base tables.

Outputs under ``--out-dir`` (default ``results/de/raloxifene-vs-control/trend``):

  * ``<quantity>_<design>.tsv`` -- per-feature trend table (log2FC, CI, SE, t, p, q,
    residual SD, mean abundance, per-feature prior variance/SD) joined to the no-trend
    SE/t/p/q and hit flags at q < 0.05 / 0.10 under each model (volcano + q-scatter
    inputs),
  * ``results/<quantity>_<design>/`` -- the ``DifferentialAbundanceTrendResult``
    (``analysis.result_io.save_result``; load with that class),
  * ``summary.json`` -- per quantity x design: hits (up/down), Storey pi0, calibration
    read, d0 and prior range, and the trend-vs-no-trend comparison (hit-set changes,
    Spearman of p and t, features whose q crosses 0.05 / 0.10, SE ratio by abundance
    tertile); plus, per quantity, the within-pair relabelling diagnostic under both
    models (pi0 + hit counts per labelling, observed rank).

Deterministic (closed-form); no seed consumed.

Run:
    ./.venv/bin/python scripts/scratch/de_trend_sensitivity.py
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import sys
import warnings
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import de_raloxifene_vs_control as base  # noqa: E402
from analysis.differential_abundance import (  # noqa: E402
    DifferentialAbundanceResult,
    DifferentialAbundanceTrendResult,
    ZeroResidualVarianceWarning,
    differential_abundance,
)
from analysis.result_io import save_result  # noqa: E402
from common.hashing import sha256_of_file  # noqa: E402
from loaders.data_loading import Dataset  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "de-raloxifene-vs-control-trend",
    "kind": "analysis",
    "provides": [],
    "uses": [
        "de_raloxifene_vs_control",
        "analysis.differential_abundance",
        "analysis.result_io",
        "common.hashing",
    ],
    "seeded_from": "differential-abundance@0.1 (via analysis/differential_abundance, "
    "project adaptation 5: trend=True)",
    "description": (
        "limma-trend sensitivity analysis of the raloxifene-d0 vs control moderated "
        "DE: trend vs no-trend for 4 quantities x 3 designs, with hit-set changes, "
        "p-value calibration and the within-pair relabelling diagnostic under both "
        "models."
    ),
}

LOG = logging.getLogger("de_trend_sensitivity")

Q_THRESHOLDS = base.Q_THRESHOLDS
TOP_N = 12
N_ABUNDANCE_TERTILES = 3


# --------------------------------------------------------------------------- #
# Fitting
# --------------------------------------------------------------------------- #
def run_design(
    dataset: Dataset, covariates: tuple[str, ...], *, trend: bool
) -> DifferentialAbundanceResult:
    """Moderated DE of condition with ``covariates``; limma-trend prior if ``trend``."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ZeroResidualVarianceWarning)
        return differential_abundance(
            dataset,
            base.CONTRAST,
            covariates=covariates,
            reference={base.CONTRAST: base.REFERENCE},
            method="moderated",
            trend=trend,
        )


def _thr_tag(thr: float) -> str:
    return f"q{round(thr * 100):02d}"


def trend_table(
    trend_result: DifferentialAbundanceResult,
    notrend_table: pd.DataFrame,
    dataset: Dataset,
    spec: base.QuantitySpec,
) -> pd.DataFrame:
    """Per-feature trend table joined to the no-trend SE/t/p/q and hit flags."""
    table = base.contrast_table(trend_result, dataset, spec)
    prior = trend_result.contrast_table.loc[:, ["feature", "prior_variance"]]
    table = table.merge(prior, on="feature", how="left", validate="1:1")
    table["prior_sd"] = np.sqrt(table["prior_variance"])
    nt = notrend_table.loc[:, ["feature", "log2fc", "se", "t", "p", "q"]].rename(
        columns={c: f"notrend_{c}" for c in ("log2fc", "se", "t", "p", "q")}
    )
    table = table.merge(nt, on="feature", how="inner", validate="1:1")
    if len(table) != len(notrend_table):
        raise ValueError("Trend and no-trend tables do not share a feature set.")
    if not np.allclose(table["log2fc"], table["notrend_log2fc"], rtol=0, atol=1e-12):
        raise ValueError("log2FC differs between trend and no-trend fits.")
    table = table.drop(columns="notrend_log2fc")
    for thr in Q_THRESHOLDS:
        tag = _thr_tag(thr)
        table[f"hit_{tag}_trend"] = table["q"] < thr
        table[f"hit_{tag}_notrend"] = table["notrend_q"] < thr
    return table


def check_against_base(table: pd.DataFrame, base_path: Path) -> None:
    """Fail loud unless the recomputed no-trend fit equals the on-disk base table."""
    disk = pd.read_csv(base_path, sep="\t", dtype={"feature": str})
    merged = disk.merge(
        table.loc[:, ["feature", "notrend_p", "notrend_q"]],
        on="feature",
        how="inner",
        validate="1:1",
    )
    if len(merged) != len(disk) or len(disk) != len(table):
        raise ValueError(f"{base_path}: feature set differs from the recomputed fit.")
    for col in ("p", "q"):
        if not np.allclose(merged[col], merged[f"notrend_{col}"], rtol=1e-12, atol=0):
            raise ValueError(
                f"{base_path}: no-trend {col} differs from the base table."
            )


# --------------------------------------------------------------------------- #
# Summaries
# --------------------------------------------------------------------------- #
def _feature_label(row: pd.Series, spec: base.QuantitySpec) -> str:
    return " | ".join(str(row[c]) for c in spec.annotation_columns)


def hit_changes(
    table: pd.DataFrame, spec: base.QuantitySpec, max_listed: int = 40
) -> dict[str, object]:
    """Features entering / leaving the hit set at each q threshold under trend."""
    out: dict[str, object] = {}
    for thr in Q_THRESHOLDS:
        tag = _thr_tag(thr)
        t_hit, n_hit = table[f"hit_{tag}_trend"], table[f"hit_{tag}_notrend"]
        gained = table[t_hit & ~n_hit]
        lost = table[~t_hit & n_hit]
        out[f"q<{thr:.2f}"] = {
            "n_trend": int(t_hit.sum()),
            "n_notrend": int(n_hit.sum()),
            "n_both": int((t_hit & n_hit).sum()),
            "n_crossing": int((t_hit != n_hit).sum()),
            "n_gained_under_trend": len(gained),
            "n_lost_under_trend": len(lost),
            "gained": [
                {
                    "feature": str(r["feature"]),
                    "label": _feature_label(r, spec),
                    "log2fc": round(float(r["log2fc"]), 4),
                    "q_trend": round(float(r["q"]), 4),
                    "q_notrend": round(float(r["notrend_q"]), 4),
                    "mean_log2_abundance": round(float(r["mean_log2_abundance"]), 3),
                }
                for _, r in gained.head(max_listed).iterrows()
            ],
            "lost": [
                {
                    "feature": str(r["feature"]),
                    "label": _feature_label(r, spec),
                    "q_trend": round(float(r["q"]), 4),
                    "q_notrend": round(float(r["notrend_q"]), 4),
                }
                for _, r in lost.head(max_listed).iterrows()
            ],
        }
    return out


def se_ratio_by_abundance(table: pd.DataFrame) -> list[dict[str, object]]:
    """Median trend/no-trend SE ratio in equal-count mean-abundance tertiles."""
    tertile = pd.qcut(
        table["mean_log2_abundance"], N_ABUNDANCE_TERTILES, labels=False
    ).to_numpy()
    ratio = (table["se"] / table["notrend_se"]).to_numpy()
    out: list[dict[str, object]] = []
    for k in range(N_ABUNDANCE_TERTILES):
        sel = tertile == k
        sub = table.loc[sel, "mean_log2_abundance"]
        out.append(
            {
                "tertile": k + 1,
                "mean_log2_abundance_range": [
                    round(float(sub.min()), 3),
                    round(float(sub.max()), 3),
                ],
                "median_se_ratio_trend_over_notrend": round(
                    float(np.median(ratio[sel])), 4
                ),
                "n_q<0.10_trend": int((table.loc[sel, "q"] < 0.10).sum()),
                "n_q<0.10_notrend": int((table.loc[sel, "notrend_q"] < 0.10).sum()),
            }
        )
    return out


def top_features(
    table: pd.DataFrame, spec: base.QuantitySpec, n: int = TOP_N
) -> list[dict[str, object]]:
    rows = table.sort_values(["q", "p"], kind="stable").head(n)
    return [
        {
            "feature": str(r["feature"]),
            "label": _feature_label(r, spec),
            "log2fc": round(float(r["log2fc"]), 4),
            "ci95": [round(float(r["ci_low"]), 4), round(float(r["ci_high"]), 4)],
            "p": float(r["p"]),
            "q": round(float(r["q"]), 4),
            "q_notrend": round(float(r["notrend_q"]), 4),
            "mean_log2_abundance": round(float(r["mean_log2_abundance"]), 3),
            "residual_sd": round(float(r["residual_sd"]), 4),
            "prior_sd": round(float(r["prior_sd"]), 4),
        }
        for _, r in rows.iterrows()
    ]


def design_block(
    trend_result: DifferentialAbundanceResult,
    notrend_result: DifferentialAbundanceResult,
    table: pd.DataFrame,
    spec: base.QuantitySpec,
) -> dict[str, object]:
    """Trend summary + no-trend reference + comparison for one quantity x design."""
    if not isinstance(trend_result, DifferentialAbundanceTrendResult):
        raise TypeError("Expected a DifferentialAbundanceTrendResult.")
    p_t = table["p"].to_numpy(dtype=float)
    p_n = table["notrend_p"].to_numpy(dtype=float)
    notrend_view = table.loc[:, ["notrend_q", "log2fc"]].rename(
        columns={"notrend_q": "q"}
    )
    prior_sd = table["prior_sd"]
    return {
        "covariates": list(trend_result.covariates),
        "n_features_tested": int(np.isfinite(p_t).sum()),
        "residual_df": trend_result.residual_df,
        "trend": {
            "prior_df_d0": base._finite_or_str(trend_result.prior_df),
            "prior_sd_range": [
                round(float(prior_sd.min()), 4),
                round(float(prior_sd.max()), 4),
            ],
            "prior_sd_median": round(float(prior_sd.median()), 4),
            "spearman_prior_sd_vs_mean_abundance": base._spearman(
                prior_sd, table["mean_log2_abundance"]
            ),
            "spline_df": trend_result.trend_spline_df,
            "knots_mean_log2_abundance": [
                round(k, 4) for k in trend_result.trend_knots
            ],
            "n_zero_residual_variance_floored": trend_result.n_prior_floored,
            "hits": base.hit_counts(table),
            "storey_pi0_lambda0.5": round(base.storey_pi0(p_t), 4),
            "calibration": base.calibration_read(p_t),
            "median_moderated_se": round(float(table["se"].median()), 4),
        },
        "notrend": {
            "prior_df_d0": base._finite_or_str(notrend_result.prior_df),
            "prior_sd": None
            if notrend_result.prior_variance is None
            else round(float(np.sqrt(notrend_result.prior_variance)), 4),
            "hits": base.hit_counts(notrend_view),
            "storey_pi0_lambda0.5": round(base.storey_pi0(p_n), 4),
            "calibration_read": base.calibration_read(p_n)["read"],
            "median_moderated_se": round(float(table["notrend_se"].median()), 4),
        },
        "trend_vs_notrend": {
            "spearman_p": base._spearman(table["p"], table["notrend_p"]),
            "spearman_t": base._spearman(table["t"], table["notrend_t"]),
            "hit_set_changes": hit_changes(table, spec),
            "se_ratio_by_abundance_tertile": se_ratio_by_abundance(table),
        },
        "top_features_trend": top_features(table, spec),
    }


# --------------------------------------------------------------------------- #
# Within-pair relabelling diagnostic (paired design, both models)
# --------------------------------------------------------------------------- #
def _relabel_row(result: DifferentialAbundanceResult) -> dict[str, object]:
    ct = result.contrast_table
    p = ct["p"].to_numpy(dtype=float)
    q = ct["q"].to_numpy(dtype=float)
    cal = base.calibration_read(p)
    return {
        "pi0": round(base.storey_pi0(p), 4),
        "density_p_lt_0.25": cal["density_p_lt_0.25"],
        "density_first_bin_p_lt_0.05": cal["density_first_bin_p_lt_0.05"],
        "hits_q<0.05": int(np.sum(q < 0.05)),
        "hits_q<0.10": int(np.sum(q < 0.10)),
        "min_q": round(float(np.nanmin(q)), 4),
        "prior_df_d0": base._finite_or_str(result.prior_df),
    }


def _rank_block(rows: list[dict[str, object]], model: str) -> dict[str, object]:
    """Observed labelling's standing among all within-pair labellings, one model."""
    obs = next(r for r in rows if r["observed"])
    vals = [r[model] for r in rows]
    o = obs[model]

    def _f(d: object, key: str) -> float:
        return float(d[key])  # type: ignore[index]

    return {
        "observed_pi0_rank_lowest_first": 1
        + sum(_f(v, "pi0") < _f(o, "pi0") for v in vals),
        "observed_min_q_rank_lowest_first": 1
        + sum(_f(v, "min_q") < _f(o, "min_q") for v in vals),
        "n_labellings_with_hits_q<0.05_ge_observed": sum(
            _f(v, "hits_q<0.05") >= _f(o, "hits_q<0.05") for v in vals
        ),
        "n_labellings_with_hits_q<0.10_ge_observed": sum(
            _f(v, "hits_q<0.10") >= _f(o, "hits_q<0.10") for v in vals
        ),
        "max_hits_q<0.05_non_observed": max(
            int(_f(r[model], "hits_q<0.05")) for r in rows if not r["observed"]
        ),
        "max_hits_q<0.10_non_observed": max(
            int(_f(r[model], "hits_q<0.10")) for r in rows if not r["observed"]
        ),
    }


def pair_relabel_diagnostic(dataset: Dataset) -> dict[str, object]:
    """Paired design under every within-pair relabelling, trend and no-trend.

    Same enumeration as the base runner (first pair's orientation fixed: 2^(n_pairs-1)
    labellings, observed included). With 8 labellings the smallest attainable
    restricted-permutation p is 1/8 = 0.125: the observed labelling can at best be
    *consistent with* a signal, never significant on this check alone. Non-observed
    labellings also break the run-order aliasing (finding 0001).
    """
    md = dataset.metadata
    pairs = sorted(md["candidate_pair"].unique())
    swap = {base.REFERENCE: base.TREATED, base.TREATED: base.REFERENCE}
    rows: list[dict[str, object]] = []
    for flips in itertools.product((False, True), repeat=len(pairs) - 1):
        flipped = [pr for pr, f in zip(pairs[1:], flips, strict=True) if f]
        in_flip = md["candidate_pair"].isin(flipped)
        relabelled = md[base.CONTRAST].where(~in_flip, md[base.CONTRAST].map(swap))
        ds = Dataset(
            abundances=dataset.abundances,
            feature_names=dataset.feature_names,
            feature_metadata=dataset.feature_metadata,
            metadata=md.assign(**{base.CONTRAST: relabelled}),
            scale=dataset.scale,
        )
        covariates = base.DESIGNS[base.PRIMARY_DESIGN]
        rows.append(
            {
                "flipped_pairs": flipped,
                "observed": not flipped,
                "trend": _relabel_row(run_design(ds, covariates, trend=True)),
                "notrend": _relabel_row(run_design(ds, covariates, trend=False)),
            }
        )
    return {
        "design": base.PRIMARY_DESIGN,
        "n_labellings": len(rows),
        "smallest_attainable_restricted_permutation_p": 1.0 / len(rows),
        "observed_standing": {
            "trend": _rank_block(rows, "trend"),
            "notrend": _rank_block(rows, "notrend"),
        },
        "labellings": rows,
    }


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def main(argv: Sequence[str] | None = None) -> int:
    root = _SCRATCH.parent.parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--qc-root", type=Path, default=root / "results" / "qc_states")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=root / "results/de/raloxifene-vs-control",
        help="no-trend base-runner outputs (checked against the recomputed fits)",
    )
    parser.add_argument(
        "--out-dir", type=Path, default=root / "results/de/raloxifene-vs-control/trend"
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.qc_root / "manifest.json").read_text(encoding="utf-8"))

    per_quantity: dict[str, object] = {}
    d0_trend: list[float] = []
    for name, spec in base.QUANTITIES.items():
        dataset, dropped = base.load_quantity(args.qc_root, spec)
        designs: dict[str, object] = {}
        for design, covariates in base.DESIGNS.items():
            trend_result = run_design(dataset, covariates, trend=True)
            notrend_result = run_design(dataset, covariates, trend=False)
            notrend_table = base.contrast_table(notrend_result, dataset, spec)
            table = trend_table(trend_result, notrend_table, dataset, spec)
            check_against_base(table, args.base_dir / f"{name}_{design}.tsv")
            stem = f"{name}_{design}"
            table.to_csv(out_dir / f"{stem}.tsv", sep="\t", index=False)
            save_result(trend_result, out_dir / "results" / stem)
            designs[design] = design_block(trend_result, notrend_result, table, spec)
            if trend_result.prior_df is not None:
                d0_trend.append(float(trend_result.prior_df))
            LOG.info(
                "%s: d0 %.3f (no-trend %.3f); q<0.05 %d (no-trend %d); "
                "q<0.10 %d (no-trend %d)",
                stem,
                trend_result.prior_df,
                notrend_result.prior_df,
                int((table["q"] < 0.05).sum()),
                int((table["notrend_q"] < 0.05).sum()),
                int((table["q"] < 0.10).sum()),
                int((table["notrend_q"] < 0.10).sum()),
            )
        per_quantity[name] = {
            "state": spec.state,
            "n_constant_features_dropped": len(dropped),
            "designs": designs,
            "pair_relabel_diagnostic": pair_relabel_diagnostic(dataset),
        }

    summary = {
        "question": "Does the raloxifene-d0 vs control DE depend on a mean-variance "
        "trend in the empirical-Bayes prior (limma-trend vs the primary no-trend "
        "model)? Positive log2FC = higher in raloxifene-d0.",
        "data_version": manifest["data_version"],
        "script": "scripts/scratch/de_trend_sensitivity.py",
        "script_sha256": sha256_of_file(Path(__file__)),
        "module": "scripts/scratch/analysis/differential_abundance.py "
        "(method='moderated', trend=True: project adaptation 5)",
        "module_sha256": sha256_of_file(
            _SCRATCH / "analysis" / "differential_abundance.py"
        ),
        "seeded_from": "differential-abundance@0.1",
        "method": "limma-style moderated t, eBayes(trend=TRUE): prior variance = "
        "natural cubic spline (4 df incl. intercept) of mean log2 abundance fitted to "
        "log residual variances (limma fitFDist with covariate); one d0; non-robust",
        "validation": "matches R limma 3.58.1 lmFit + eBayes(trend=TRUE) on all 12 "
        "fits to ~1e-12 (r_agreement.json)",
        "inputs": "identical to the no-trend runner (states, log2 transforms, "
        "constant-feature drops, designs); no-trend fits recomputed and checked equal "
        "to the base tables",
        "contrast": {
            "column": base.CONTRAST,
            "reference": base.REFERENCE,
            "level": base.TREATED,
        },
        "designs": {d: list(c) for d, c in base.DESIGNS.items()},
        "primary_design": base.PRIMARY_DESIGN,
        "correction": "Benjamini-Hochberg FDR per quantity x design over the features "
        "tested (contrast term)",
        "ci": "95% t interval on the moderated SE (df = residual_df + d0)",
        "sample_set": "experimental (8 of 8)",
        "pi0_estimator": f"Storey, fixed lambda = {base.PI0_LAMBDA}",
        "d0_range_trend_all_fits": [round(min(d0_trend), 4), round(max(d0_trend), 4)],
        "quantities": per_quantity,
        "table_columns": "feature, annotations, log2fc, ci_low, ci_high, se, t, p, q "
        "(trend); residual_sd, mean_log2_abundance, prior_variance, prior_sd (trend "
        "prior at the feature's mean abundance); notrend_se, notrend_t, notrend_p, "
        "notrend_q; hit_q05_/hit_q10_ trend/notrend flags",
        "result_objects": "results/<quantity>_<design>/ -- load with analysis."
        "result_io.load_result(path, DifferentialAbundanceTrendResult)",
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    LOG.info("wrote %s", out_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
