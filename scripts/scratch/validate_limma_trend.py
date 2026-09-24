"""Validate the project limma-trend option against R limma ``eBayes(trend=TRUE)``.

Exports the exact matrices the trend runner tests (same states, transforms and constant
drops, via ``de_raloxifene_vs_control.load_quantity``) to a temporary directory, runs
``limma_trend_check.R`` (``lmFit`` + ``eBayes(trend=TRUE)`` and ``trend=FALSE``) for all
4 quantities x 3 designs, and compares per feature: t, p, BH q, the per-feature prior
variance and the prior df ``d0`` (plus the no-trend t/p/d0 as a regression check).
Writes ``<out-dir>/r_agreement.json``. Needs ``Rscript`` with limma on PATH.

Run:
    ./.venv/bin/python scripts/scratch/validate_limma_trend.py
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import tempfile
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
    ZeroResidualVarianceWarning,
    differential_abundance,
)
from loaders.data_loading import Dataset  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "de-raloxifene-vs-control-trend",
    "kind": "validation",
    "provides": [],
    "uses": ["de_raloxifene_vs_control", "analysis.differential_abundance"],
    "seeded_from": None,
    "description": (
        "Numerical agreement of the project limma-trend implementation with R limma "
        "(lmFit + eBayes(trend=TRUE)) on all 12 quantity x design fits."
    ),
}

LOG = logging.getLogger("validate_limma_trend")
R_SCRIPT = _SCRATCH / "limma_trend_check.R"


def python_fit(
    dataset: Dataset, covariates: tuple[str, ...], trend: bool
) -> DifferentialAbundanceResult:
    """The project fit (contrast table + prior) for one design."""
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


def _max_abs(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b)))


def _max_rel(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))


def compare(
    dataset: Dataset, covariates: tuple[str, ...], r_table: pd.DataFrame
) -> dict[str, object]:
    """Per-feature agreement of the Python trend / no-trend fits with R."""
    tr = python_fit(dataset, covariates, trend=True)
    nt = python_fit(dataset, covariates, trend=False)
    ct = tr.contrast_table.set_index("feature").reindex(r_table["feature"])
    cn = nt.contrast_table.set_index("feature").reindex(r_table["feature"])
    if ct["p"].isna().any() or cn["p"].isna().any():
        raise ValueError("Feature sets differ between Python and R.")
    if tr.prior_df is None or nt.prior_df is None:
        raise ValueError("Moderated fit returned no prior df.")
    r_d0 = float(r_table["df_prior"].iloc[0])
    r_d0_nt = float(r_table["notrend_df_prior"].iloc[0])
    return {
        "n_features": len(r_table),
        "d0": {"python": tr.prior_df, "r": r_d0, "abs_diff": abs(tr.prior_df - r_d0)},
        "max_abs_diff_log2fc": _max_abs(ct["effect"].to_numpy(), r_table["coef"]),
        "max_abs_diff_t": _max_abs(ct["statistic"].to_numpy(), r_table["t"]),
        "max_rel_diff_p": _max_rel(ct["p"].to_numpy(), r_table["p"].to_numpy()),
        "max_abs_diff_q": _max_abs(ct["q"].to_numpy(), r_table["q"]),
        "max_rel_diff_prior_variance": _max_rel(
            ct["prior_variance"].to_numpy(), r_table["s2_prior"].to_numpy()
        ),
        "max_abs_diff_mean_abundance": _max_abs(
            ct["mean_abundance"].to_numpy(), r_table["amean"]
        ),
        "hits_q<0.05": {
            "python": int((ct["q"] < 0.05).sum()),
            "r": int((r_table["q"] < 0.05).sum()),
        },
        "hits_q<0.10": {
            "python": int((ct["q"] < 0.10).sum()),
            "r": int((r_table["q"] < 0.10).sum()),
        },
        "notrend": {
            "d0_abs_diff": abs(nt.prior_df - r_d0_nt),
            "max_abs_diff_t": _max_abs(
                cn["statistic"].to_numpy(), r_table["notrend_t"]
            ),
            "max_rel_diff_p": _max_rel(
                cn["p"].to_numpy(), r_table["notrend_p"].to_numpy()
            ),
        },
    }


def main(argv: Sequence[str] | None = None) -> int:
    root = _SCRATCH.parent.parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--qc-root", type=Path, default=root / "results" / "qc_states")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "results/de/raloxifene-vs-control/trend",
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    datasets = {
        name: base.load_quantity(args.qc_root, spec)[0]
        for name, spec in base.QUANTITIES.items()
    }
    with tempfile.TemporaryDirectory() as tmp:
        in_dir, r_out = Path(tmp) / "in", Path(tmp) / "out"
        in_dir.mkdir()
        md = next(iter(datasets.values())).metadata
        md.loc[:, ["sample_id", base.CONTRAST, "candidate_pair", "batch"]].to_csv(
            in_dir / "samples.tsv", sep="\t", index=False
        )
        for name, ds in datasets.items():
            if not ds.metadata.index.equals(md.index):
                raise ValueError(f"{name}: sample order differs.")
            frame = pd.DataFrame(
                ds.abundances.T, columns=list(md["sample_id"].astype(str))
            )
            frame.insert(0, "feature", ds.feature_names.astype(str))
            frame.to_csv(
                in_dir / f"{name}.tsv", sep="\t", index=False, float_format="%.17g"
            )
        proc = subprocess.run(
            ["Rscript", str(R_SCRIPT), str(in_dir), str(r_out)],
            check=True,
            capture_output=True,
            text=True,
        )
        limma_version = proc.stdout.strip()
        LOG.info("R: %s", limma_version)
        agreement: dict[str, object] = {}
        for name, ds in datasets.items():
            for design, covariates in base.DESIGNS.items():
                r_table = pd.read_csv(
                    r_out / f"{name}_{design}.tsv",
                    sep="\t",
                    quoting=3,
                    keep_default_na=False,
                    dtype={"feature": str},
                )
                res = compare(ds, covariates, r_table)
                agreement[f"{name}_{design}"] = res
                LOG.info(
                    "%s %s: d0 diff %.2e, max|dt| %.2e, max rel dp %.2e, max|dq| %.2e",
                    name,
                    design,
                    res["d0"]["abs_diff"],  # type: ignore[index]
                    res["max_abs_diff_t"],
                    res["max_rel_diff_p"],
                    res["max_abs_diff_q"],
                )
    out = {
        "reference": f"R {limma_version}: lmFit + eBayes(trend=TRUE) (and "
        "trend=FALSE); BH via p.adjust on the condition coefficient",
        "script": "scripts/scratch/validate_limma_trend.py + limma_trend_check.R",
        "fits": agreement,
    }
    path = args.out_dir / "r_agreement.json"
    path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    LOG.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
