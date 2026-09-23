"""Stage-3 QC runner: protein dynamic-range (rank-abundance) figure, raw linear.

Loads the prepared ``raw_linear`` protein state
(``results/qc_states/protein/raw_linear``, written once by ``qc_prep.py``; all 4344
protein groups x 8 experimental samples, NaN = not detected, contaminants flagged in
``feature_metadata.is_contaminant``) and renders the whole-cohort rank-abundance curve
(median of detected values + IQR band across samples) with every detected contaminant
protein marked (registry group ``highlight_group`` / ``contaminant``) and only the
``--label-top-n`` most-abundant contaminants text-labelled.

Writes ``<output-dir>/<base-name>.{svg,png}`` + ``.legend.{svg,png}`` and a small JSON
summary (key numbers + provenance: script sha256, data version, params) to
``--summary-file``. Deterministic (no stochastic step).

Run:
    ./.venv/bin/python scripts/promoted/qc_fig_dynamic_range.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from loaders.dataset_io import load_dataset  # noqa: E402
from qc_figures.dynamic_range import save_dynamic_range  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "stage3-qc-dynamic-range",
    "kind": "script",
    "provides": [],
    "uses": ["loaders.dataset_io", "qc_figures.dynamic_range"],
    "seeded_from": None,
    "description": (
        "Render the Stage-3 protein dynamic-range figure (raw, linear) with "
        "contaminant proteins highlighted; write figure + legend + JSON summary."
    ),
}

LOG = logging.getLogger(__name__)
CONTAMINANT_GROUP = "contaminant"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--state-dir", type=Path, default=Path("results/qc_states/protein/raw_linear")
    )
    parser.add_argument(
        "--manifest", type=Path, default=Path("results/qc_states/manifest.json")
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("figures/qc/dynamic-range")
    )
    parser.add_argument(
        "--base-name", type=str, default="dynamic-range-protein-raw-linear"
    )
    parser.add_argument(
        "--summary-file",
        type=Path,
        default=Path("results/stage3/dynamic_range/protein_raw_linear.summary.json"),
    )
    parser.add_argument("--label-top-n", type=int, default=5)
    parser.add_argument("--min-intensity", type=float, default=0.0)
    parser.add_argument(
        "--registry", type=Path, default=Path("state/color_registry.json")
    )
    parser.add_argument("--log-level", type=str.upper, default="INFO")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(level=args.log_level, format="%(levelname)s %(message)s")

    dataset = load_dataset(args.state_dir)
    if dataset.scale != "linear":
        raise ValueError(f"expected a linear state; got {dataset.scale!r}")
    fm = dataset.feature_metadata
    for col in ("is_contaminant", "entry"):
        if col not in fm.columns:
            raise ValueError(f"feature_metadata lacks required column {col!r}.")
    names = np.asarray(dataset.feature_names).astype(str)
    if len(fm) != names.size:
        raise ValueError("feature_metadata is not aligned to feature_names.")
    is_contam = fm["is_contaminant"].to_numpy(dtype=bool)
    entries = fm["entry"].astype(str).to_numpy()

    abund = np.asarray(dataset.abundances, dtype=float)
    detected_any = (np.isfinite(abund) & (abund > args.min_intensity)).any(axis=0)
    contam_ids = names[is_contam]
    undetected_contam = names[is_contam & ~detected_any]
    if undetected_contam.size:
        LOG.warning(
            "%d contaminant(s) never detected; not marked: %s",
            undetected_contam.size,
            list(undetected_contam),
        )
    mark = is_contam & detected_any
    highlight = {
        str(fid): str(entry)
        for fid, entry in zip(names[mark], entries[mark], strict=True)
    }
    groups = dict.fromkeys(highlight, CONTAMINANT_GROUP)

    params: dict[str, Any] = {
        "state_dir": str(args.state_dir),
        "label_top_n": args.label_top_n,
        "min_intensity": args.min_intensity,
        "highlight": "feature_metadata.is_contaminant (detected in >=1 sample)",
        "mode": "whole-cohort median + IQR band",
    }
    artifacts, plot = save_dynamic_range(
        dataset,
        args.output_dir,
        args.base_name,
        highlight_features=highlight,
        highlight_groups=groups,
        label_top_n=args.label_top_n,
        min_intensity=args.min_intensity,
        title="Protein dynamic range (raw, linear; all 8 samples)",
        registry_path=args.registry,
    )

    res = plot.result
    ranked = res.feature_names_ranked.astype(str)
    rank_of = {fid: i + 1 for i, fid in enumerate(ranked)}
    entry_of = dict(zip(names, entries, strict=True))
    contam_set = set(contam_ids)
    contam_ranks = sorted(rank_of[f] for f in highlight)
    k = res.n_features_detected

    def _row(i: int) -> dict[str, Any]:
        fid = ranked[i]
        return {
            "rank": i + 1,
            "entry": entry_of[fid],
            "feature": fid,
            "contaminant": fid in contam_set,
            "log2_median": round(float(res.log2_median[i]), 3),
            "n_detected": int(res.n_detected[i]),
        }

    manifest = json.loads(args.manifest.read_text()) if args.manifest.exists() else {}
    summary: dict[str, Any] = {
        "figure": {
            "svg": str(artifacts.svg),
            "png": str(artifacts.png),
            "legend_svg": str(artifacts.legend_svg),
            "legend_png": str(artifacts.legend_png),
        },
        "color_map": plot.color_map,
        "n_samples": res.n_samples,
        "n_features_total": res.n_features_total,
        "n_features_detected": k,
        "dynamic_range_orders": round(res.dynamic_range_orders, 3),
        "log2_median_max": round(float(res.log2_median[0]), 3),
        "log2_median_min": round(float(res.log2_median[-1]), 3),
        "n_detected_all_samples": int((res.n_detected == res.n_samples).sum()),
        "n_detected_single_sample": int((res.n_detected == 1).sum()),
        "top20": [_row(i) for i in range(min(20, k))],
        "contaminants": {
            "n_flagged": int(is_contam.sum()),
            "n_marked": len(highlight),
            "n_in_top10": sum(r <= 10 for r in contam_ranks),
            "n_in_top50": sum(r <= 50 for r in contam_ranks),
            "n_in_top100": sum(r <= 100 for r in contam_ranks),
            "n_in_top_decile": sum(r <= k / 10 for r in contam_ranks),
            "median_rank": float(np.median(contam_ranks)) if contam_ranks else None,
            "ranks": contam_ranks,
            "labelled": [_row(r - 1) for r in contam_ranks[: args.label_top_n]],
        },
        "provenance": {
            "script": "scripts/promoted/qc_fig_dynamic_range.py",
            "script_sha256": _sha256(Path(__file__).resolve()),
            "module": "scripts/promoted/qc_figures/dynamic_range.py",
            "module_sha256": _sha256(_SCRATCH / "qc_figures" / "dynamic_range.py"),
            "git_commit": None,
            "data_version": manifest.get("data_version"),
            "params": params,
        },
    }
    args.summary_file.parent.mkdir(parents=True, exist_ok=True)
    args.summary_file.write_text(json.dumps(summary, indent=2) + "\n")
    LOG.info(
        "wrote %s, %s (+ legend) and %s",
        artifacts.svg,
        artifacts.png,
        args.summary_file,
    )


if __name__ == "__main__":
    main()
