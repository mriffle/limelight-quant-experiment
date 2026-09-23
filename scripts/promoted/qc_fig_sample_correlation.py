"""Stage-3 QC runner: sample-correlation heatmaps (protein + peptide, normalized log2).

Loads the prepared processing state ``results/qc_states/<level>/normalized_log/``
(median-normalized, log2, complete non-contaminant features — written once by
``qc_prep.py``), renders a clustered Pearson sample-correlation heatmap with
condition / batch / candidate-pair stripes via :mod:`qc_figures.correlation`, and
writes ``<out>/sample-correlation-<level>-normalized-log2.{svg,png}`` + legend.

Also writes a regenerable numeric companion under ``results/qc_figures/
sample-correlation/``: the clustered correlation matrix per level (TSV) and a
provenance JSON (script sha256, data version, params, off-diagonal min/max).

Usage::

    uv run python scripts/promoted/qc_fig_sample_correlation.py \
        [--out-dir figures/qc/sample-correlation] [--levels protein peptide]
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

import matplotlib

matplotlib.use("Agg")

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402
from loaders.dataset_io import load_dataset  # noqa: E402
from qc_figures.correlation import save_sample_correlation  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "stage3-qc-sample-correlation",
    "kind": "script",
    "provides": ["main"],
    "uses": ["loaders.dataset_io", "qc_figures.correlation"],
    "seeded_from": {"template": "sample-correlation", "version": "0.1"},
    "description": (
        "Renders the Stage-3 QC sample-correlation family (Pearson, average-linkage "
        "clustered, condition/batch/candidate_pair stripes) for protein and peptide "
        "from the prepared normalized_log state; writes figures + a numeric/provenance "
        "companion under results/qc_figures/sample-correlation/."
    ),
}

LOG = logging.getLogger("qc_fig_sample_correlation")

PROJECT_ROOT = _SCRATCH.parent.parent
STATE = "normalized_log"
METHOD = "pearson"
DEFAULT_OUT_DIR = PROJECT_ROOT / "figures" / "qc" / "sample-correlation"
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "results" / "qc_figures" / "sample-correlation"
STATES_ROOT = PROJECT_ROOT / "results" / "qc_states"
REGISTRY = PROJECT_ROOT / "state" / "color_registry.json"
ANNOTATIONS: dict[str, str] = {
    "condition": "condition",
    "batch": "batch",
    "candidate_pair": "candidate pair",
}
FEATURE_NOUN: dict[str, str] = {"protein": "proteins", "peptide": "peptides"}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def render_level(
    level: str, out_dir: Path, results_dir: Path, stem_prefix: str = ""
) -> dict[str, Any]:
    """Render one level's heatmap and write its numeric companion; return a summary.

    ``stem_prefix`` (default empty) is prepended to the figure stem and to the
    companion TSV name, e.g. ``"0003-"`` for a finding-id-prefixed copy.
    """
    dataset = load_dataset(STATES_ROOT / level / STATE)
    if dataset.scale != "log2":
        raise ValueError(f"{level}/{STATE} is on scale {dataset.scale!r}, not log2.")
    n_samples, n_features = dataset.abundances.shape
    title = (
        f"Sample correlation, {level} ({n_features:,} {FEATURE_NOUN[level]}, "
        f"n = {n_samples})\nmedian-normalized, log2 · Pearson · average linkage"
    )
    base = f"{stem_prefix}sample-correlation-{level}-normalized-log2"
    artifacts, plot = save_sample_correlation(
        dataset,
        out_dir,
        base,
        annotations=ANNOTATIONS,
        method=METHOD,
        title=title,
        registry_path=REGISTRY,
    )
    result = plot.result
    frame = result.ordered_frame()
    results_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(results_dir / f"{stem_prefix}{level}_pearson_clustered.tsv", sep="\t")

    off = result.off_diagonal()
    iu = np.triu_indices(n_samples, k=1)
    pairs = [
        (
            str(result.sample_ids[i]),
            str(result.sample_ids[j]),
            float(result.matrix[i, j]),
        )
        for i, j in zip(iu[0], iu[1], strict=True)
    ]
    lo = min(pairs, key=lambda p: p[2])
    hi = max(pairs, key=lambda p: p[2])
    summary: dict[str, Any] = {
        "level": level,
        "n_samples": int(n_samples),
        "n_features": int(n_features),
        "clustered_order": [str(s) for s in result.sample_ids[result.order]],
        "off_diagonal_min": {"pair": lo[:2], "r": lo[2]},
        "off_diagonal_max": {"pair": hi[:2], "r": hi[2]},
        "off_diagonal_median": float(np.median(off)),
        "color_maps": plot.color_maps,
        "artifacts": {
            k: str(Path(v).resolve().relative_to(PROJECT_ROOT)) if v else None
            for k, v in {
                "svg": artifacts.svg,
                "png": artifacts.png,
                "legend_svg": artifacts.legend_svg,
                "legend_png": artifacts.legend_png,
            }.items()
        },
    }
    LOG.info("%s: %s", level, json.dumps(summary, indent=1))
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument(
        "--levels",
        nargs="+",
        default=["protein", "peptide"],
        choices=["protein", "peptide"],
    )
    parser.add_argument(
        "--stem-prefix",
        default="",
        help="Prefix for figure stems + companion files (e.g. '0003-'); default none.",
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    manifest = json.loads((STATES_ROOT / "manifest.json").read_text())
    summaries = [
        render_level(lv, args.out_dir, args.results_dir, args.stem_prefix)
        for lv in args.levels
    ]
    module_path = _SCRATCH / "qc_figures" / "correlation.py"
    provenance = {
        "scripts": {
            str(Path(__file__).resolve().relative_to(PROJECT_ROOT)): _sha256(
                Path(__file__).resolve()
            ),
            str(module_path.relative_to(PROJECT_ROOT)): _sha256(module_path),
        },
        "commit": None,  # project is not a git repository; script sha256 pins it.
        "data_version": manifest["data_version"],
        "input_state": STATE,
        "stem_prefix": args.stem_prefix,
        "params": {
            "method": METHOD,
            "linkage": "average",
            "linkage_metric": "euclidean (over correlation-matrix rows)",
            "annotations": list(ANNOTATIONS),
            "diagonal": "masked (gray); color scale = off-diagonal min..max",
            "cell_text_decimals": 3,
            "dpi": 300,
        },
        "levels": summaries,
    }
    args.results_dir.mkdir(parents=True, exist_ok=True)
    prov_path = args.results_dir / f"{args.stem_prefix}provenance.json"
    prov_path.write_text(json.dumps(provenance, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
