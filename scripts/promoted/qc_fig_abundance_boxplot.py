"""Stage-3 QC figure family ``abundance-boxplot``: per-sample log2 abundance boxes.

For each level in {protein, peptide}, loads the three prepared log2 processing states
written by ``qc_prep.py`` (``results/qc_states/<level>/{raw_log,normalized_log,
batch_corrected_log}/``; complete non-contaminant features, the 8 experimental samples
in acquisition order) and renders one figure of three stacked box-plot panels
(raw -> median-normalized -> batch-corrected ComBat preview), boxes labeled by
``sample_id`` with condition + batch annotation stripes (colors from
``state/color_registry.json``). Load -> check -> plot only: no normalization or ComBat
is re-run here.

Writes, per level, to ``--output-dir`` (default ``figures/qc/abundance-boxplot``):
``abundance-boxplot-<level>-raw-normalized-batchcorrected-log2.{svg,png}`` and its
``.legend.{svg,png}``; plus a provenance/key-numbers JSON (``--provenance-out``).

No stochastic step, so no seed.

Run:
    ./.venv/bin/python scripts/promoted/qc_fig_abundance_boxplot.py \\
        --qc-states-dir results/qc_states \\
        --output-dir figures/qc/abundance-boxplot \\
        --provenance-out results/stage3/figure_provenance/abundance-boxplot.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from loaders.data_loading import Dataset  # noqa: E402
from loaders.dataset_io import load_dataset  # noqa: E402
from qc_figures.abundance_boxplot import (  # noqa: E402
    median_range,
    save_abundance_boxplots,
)

__script_meta__: dict[str, object] = {
    "task": "qc-figures-abundance-boxplot",
    "kind": "figure",
    "provides": [],
    "uses": [
        "loaders.dataset_io",
        "loaders.data_loading",
        "qc_figures.abundance_boxplot",
    ],
    "seeded_from": None,
    "description": (
        "Stage-3 QC abundance-boxplot family runner: loads the prepared raw/"
        "normalized/batch-corrected log2 states per level (protein, peptide), "
        "verifies sample alignment + acquisition order + experimental-only, and "
        "renders the three-panel per-sample box plot with sample_id ticks and "
        "condition/batch stripes; records per-sample medians and provenance."
    ),
}

LOGGER = logging.getLogger("qc_fig_abundance_boxplot")

LEVELS: tuple[str, ...] = ("protein", "peptide")

# (state directory, panel label: processing state + scale) in top-to-bottom order.
STATES: tuple[tuple[str, str], ...] = (
    ("raw_log", "raw  ·  log2"),
    ("normalized_log", "median-normalized  ·  log2"),
    (
        "batch_corrected_log",
        "batch-corrected (ComBat preview, batch only)  ·  log2",
    ),
)

CATEGORICAL_ANNOTATIONS: tuple[str, ...] = ("condition", "batch")
SAMPLE_LABEL_COLUMN = "sample_id"
ORDER_COLUMNS: tuple[str, ...] = ("acq_date", "seq_number")
BASE_NAME_TEMPLATE = "abundance-boxplot-{level}-raw-normalized-batchcorrected-log2"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--qc-states-dir", type=Path, default=Path("results/qc_states"))
    parser.add_argument(
        "--output-dir", type=Path, default=Path("figures/qc/abundance-boxplot")
    )
    parser.add_argument(
        "--provenance-out",
        type=Path,
        default=Path("results/stage3/figure_provenance/abundance-boxplot.json"),
    )
    parser.add_argument(
        "--registry", type=Path, default=Path("state/color_registry.json")
    )
    parser.add_argument("--levels", nargs="+", default=list(LEVELS), choices=LEVELS)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args(argv)


def check_acquisition_order(dataset: Dataset, level: str) -> None:
    """Refuse samples that are not experimental-only and in acquisition order."""
    meta = dataset.metadata
    for column in (*ORDER_COLUMNS, "sample_role", SAMPLE_LABEL_COLUMN):
        if column not in meta.columns:
            raise ValueError(f"{level}: metadata lacks required column {column!r}.")
    roles = set(meta["sample_role"].astype(str))
    if roles != {"experimental"}:
        raise ValueError(
            f"{level}: expected experimental samples only (controls render in their "
            f"own call); got roles {sorted(roles)}."
        )
    order_key = pd.DataFrame(
        {
            "date": pd.to_datetime(meta["acq_date"].astype(str), errors="raise"),
            "seq": pd.to_numeric(meta["seq_number"], errors="raise"),
        }
    )
    sorted_idx = order_key.sort_values(["date", "seq"], kind="stable").index
    if not sorted_idx.equals(meta.index):
        raise ValueError(
            f"{level}: samples are not in acquisition order (acq_date, seq_number): "
            f"{meta.index.tolist()} vs expected {sorted_idx.tolist()}."
        )


def load_level(qc_states_dir: Path, level: str) -> dict[str, Dataset]:
    """Load the three log2 states for ``level`` keyed by panel label (panel order)."""
    datasets: dict[str, Dataset] = {}
    for state_dir, label in STATES:
        path = qc_states_dir / level / state_dir
        dataset = load_dataset(path)
        if dataset.scale != "log2":
            raise ValueError(f"{path}: expected scale 'log2', got {dataset.scale!r}.")
        if not bool(np.isfinite(dataset.abundances).all()):
            raise ValueError(
                f"{path}: expected a complete (all-finite) matrix; found "
                f"{int((~np.isfinite(dataset.abundances)).sum())} non-finite cells."
            )
        check_acquisition_order(dataset, level)
        LOGGER.info(
            "%s/%s: %d samples x %d features (%s)",
            level,
            state_dir,
            dataset.abundances.shape[0],
            dataset.abundances.shape[1],
            dataset.scale,
        )
        datasets[label] = dataset
    n_features = {d.abundances.shape[1] for d in datasets.values()}
    if len(n_features) != 1:
        raise ValueError(f"{level}: states differ in feature count {n_features}.")
    return datasets


def _sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def _git_commit(repo_dir: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo_dir), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.stdout.strip() or None


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    logging.basicConfig(
        level=args.log_level, format="%(asctime)s %(levelname)s %(message)s"
    )
    LOGGER.info("params: %s", args)

    manifest_path = args.qc_states_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    data_version = manifest["data_version"]

    module_path = _SCRATCH / "qc_figures" / "abundance_boxplot.py"
    provenance: dict[str, Any] = {
        "family": "abundance-boxplot",
        "scripts": {
            "runner": {
                "path": "scripts/promoted/qc_fig_abundance_boxplot.py",
                "sha256": _sha256(Path(__file__).resolve()),
            },
            "module": {
                "path": "scripts/promoted/qc_figures/abundance_boxplot.py",
                "sha256": _sha256(module_path),
            },
            "git_head": _git_commit(_SCRATCH),
        },
        "data_version": data_version,
        "inputs_manifest": str(manifest_path),
        "params": {
            "qc_states_dir": str(args.qc_states_dir),
            "output_dir": str(args.output_dir),
            "states": [s for s, _ in STATES],
            "categorical_annotations": list(CATEGORICAL_ANNOTATIONS),
            "sample_label_column": SAMPLE_LABEL_COLUMN,
            "box_colormap": "cool",
            "show_outliers": True,
            "dpi": args.dpi,
        },
        "levels": {},
    }

    for level in args.levels:
        datasets = load_level(args.qc_states_dir, level)
        reference = next(iter(datasets.values()))
        n_features = reference.abundances.shape[1]
        base_name = BASE_NAME_TEMPLATE.format(level=level)
        title = (
            f"{level.capitalize()} abundance per sample  ·  {n_features:,} complete "
            f"non-contaminant {level}s"
        )
        artifacts, result = save_abundance_boxplots(
            datasets,
            args.output_dir,
            base_name,
            categorical_annotations=CATEGORICAL_ANNOTATIONS,
            feature_type=level,
            title=title,
            sample_label_column=SAMPLE_LABEL_COLUMN,
            annotate_median_range=True,
            registry_path=args.registry,
            persist_colors=False,
            dpi=args.dpi,
        )
        panels: dict[str, Any] = {}
        for (state_dir, label), medians in zip(
            STATES, result.medians.values(), strict=True
        ):
            spread = median_range(medians)
            panels[state_dir] = {
                "panel_label": label,
                "per_sample_median": {
                    str(sid): round(float(m), 6)
                    for sid, m in zip(result.sample_ids, medians, strict=True)
                },
                "median_range": round(spread, 6),
            }
            LOGGER.info("%s %s: per-sample median range %.4f", level, state_dir, spread)
        provenance["levels"][level] = {
            "n_features": n_features,
            "n_samples": int(result.sample_ids.size),
            "sample_order": [str(s) for s in result.sample_ids],
            "panels": panels,
            "artifacts": {
                "svg": str(artifacts.svg),
                "png": str(artifacts.png),
                "legend_svg": str(artifacts.legend_svg),
                "legend_png": str(artifacts.legend_png),
            },
        }
        LOGGER.info("%s: wrote %s", level, artifacts)

    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(
        json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8"
    )
    LOGGER.info("wrote provenance to %s", args.provenance_out)


if __name__ == "__main__":
    main()
