"""Stage-3 QC figure family ``cv``: per-feature CV across processing states.

For each level in {protein, peptide}, loads three prepared processing states from
``results/qc_states/<level>/`` (written once by ``scripts/promoted/qc_prep.py``) --

  * ``raw_linear_complete``    -> "Raw"
  * ``normalized_linear``      -> "Median-normalized"
  * ``batch_corrected_linear`` -> "ComBat batch-corrected" (ComBat on log2 data,
    de-logged to linear by qc_prep; an approximation of a linear-scale correction)

-- and overlays their per-feature CV (std/mean, ddof=1) distributions computed across
ALL 8 samples, writing ``figures/qc/cv/cv-experimental-<level>-raw-normalized-
batchcorrected-linear.{svg,png}`` plus its ``.legend.{svg,png}``.

Every sample in this project is experimental (``sample_role == "experimental"`` for all
8; there are no control/QC pools), so there is no per-control-type CV figure. The
script asserts this (fail loud) rather than assuming it.

It is load -> plot only: no normalization / ComBat is re-run here. Per-figure provenance
(script + module sha256, input file sha256s, the qc_states manifest data_version,
parameters, medians, color map, artifact paths) is written to
``results/stage3/qc_figures/cv_provenance.json``.

No stochastic step (gaussian_kde is deterministic), so there is no seed.

Run (from the project root):
    ./.venv/bin/python scripts/promoted/qc_fig_cv.py
"""

from __future__ import annotations

import argparse
import json
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

from common.hashing import sha256_of_file  # noqa: E402
from loaders.data_loading import Dataset  # noqa: E402
from loaders.dataset_io import load_dataset  # noqa: E402
from qc_figures import cv as cv_module  # noqa: E402

__script_meta__: dict[str, object] = {
    "kind": "script",
    "uses": ["qc_figures.cv", "loaders.dataset_io", "common.hashing"],
    "description": (
        "Stage-3 QC CV family: overlay per-feature CV distributions of the raw, "
        "median-normalized and ComBat batch-corrected (de-logged) linear states at "
        "protein and peptide level, all 8 (experimental) samples."
    ),
}

LEVELS: tuple[str, ...] = ("protein", "peptide")

# (state directory, display label) in draw / legend order.
STATES: tuple[tuple[str, str], ...] = (
    ("raw_linear_complete", "Raw"),
    ("normalized_linear", "Median-normalized"),
    ("batch_corrected_linear", "ComBat batch-corrected"),
)

EXPECTED_N_SAMPLES = 8
LEGEND_TITLE = "Processing state (linear)"


def base_name(level: str) -> str:
    """File stem for one level's CV figure (processing states + scale in the stem)."""
    return f"cv-experimental-{level}-raw-normalized-batchcorrected-linear"


def _check_experimental(level: str, datasets: dict[str, Dataset]) -> None:
    """Assert every state holds the same 8 experimental samples (fail loud)."""
    reference: list[str] | None = None
    for label, ds in datasets.items():
        meta = ds.metadata
        if "sample_role" not in meta.columns:
            raise ValueError(f"{level}/{label}: metadata lacks a 'sample_role' column.")
        roles = sorted({str(r) for r in meta["sample_role"]})
        if roles != ["experimental"]:
            raise ValueError(
                f"{level}/{label}: expected only experimental samples, got roles "
                f"{roles}; control samples must be plotted separately."
            )
        ids = [str(s) for s in meta.index]
        if len(ids) != EXPECTED_N_SAMPLES or ds.abundances.shape[0] != len(ids):
            raise ValueError(
                f"{level}/{label}: expected {EXPECTED_N_SAMPLES} samples, got "
                f"{len(ids)} metadata rows / {ds.abundances.shape[0]} abundance rows."
            )
        if reference is None:
            reference = ids
        elif ids != reference:
            raise ValueError(
                f"{level}/{label}: sample order {ids} differs from {reference}."
            )
        n_nonfinite = int(np.count_nonzero(~np.isfinite(ds.abundances)))
        n_nonpos = int(np.count_nonzero(ds.abundances <= 0))
        if n_nonfinite or n_nonpos:
            raise ValueError(
                f"{level}/{label}: complete linear state has {n_nonfinite} non-finite "
                f"and {n_nonpos} non-positive values; expected none."
            )


def _state_files(state_dir: Path) -> dict[str, dict[str, str]]:
    """sha256 of every file in one prepared state directory."""
    return {
        p.name: {"path": str(p), "sha256": sha256_of_file(p)}
        for p in sorted(state_dir.iterdir())
        if p.is_file()
    }


def run(
    qc_states_dir: Path,
    output_dir: Path,
    registry_path: Path,
    provenance_path: Path,
    n_bins: int,
) -> dict[str, Any]:
    """Render both levels; return (and write) the provenance record."""
    manifest_path = qc_states_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    record: dict[str, Any] = {
        "script": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256_of_file(Path(__file__).resolve()),
        },
        "module": {
            "path": str(Path(cv_module.__file__).resolve()),
            "sha256": sha256_of_file(Path(cv_module.__file__).resolve()),
            "seeded_from": cv_module.__script_meta__["seeded_from"],
        },
        "git_commit": None,  # project is not a git repository
        "data_version": manifest.get("data_version"),
        "qc_states_manifest": {
            "path": str(manifest_path),
            "sha256": sha256_of_file(manifest_path),
        },
        "color_registry": str(registry_path),
        "params": {
            "n_bins": n_bins,
            "max_cv": None,
            "x_upper": "99.5th percentile of pooled CVs",
            "cv": "std(ddof=1)/mean per feature across all 8 samples, linear scale",
            "category": cv_module.DEFAULT_CV_CATEGORY,
            "states": [s for s, _ in STATES],
        },
        "sample_subset": "all 8 samples (all experimental; no control/QC samples)",
        "figures": {},
    }
    for level in LEVELS:
        datasets: dict[str, Dataset] = {}
        inputs: dict[str, Any] = {}
        for state, label in STATES:
            state_dir = qc_states_dir / level / state
            datasets[label] = load_dataset(state_dir)
            inputs[state] = _state_files(state_dir)
        _check_experimental(level, datasets)
        n_features = int(next(iter(datasets.values())).abundances.shape[1])
        noun = "Protein" if level == "protein" else "Peptide"
        artifacts, plot = cv_module.save_cv(
            datasets,
            output_dir,
            base_name(level),
            feature_type=level,
            xlabel=f"{noun} CV (std / mean), linear scale",
            n_bins=n_bins,
            title=(
                f"{noun} CV across 8 experimental samples · "
                f"n = {n_features:,} complete {level}s"
            ),
            legend_title=LEGEND_TITLE,
            registry_path=registry_path,
            persist_colors=True,
        )
        finite_counts = {
            k: int(np.count_nonzero(np.isfinite(v))) for k, v in plot.result.cvs.items()
        }
        record["figures"][base_name(level)] = {
            "level": level,
            "n_features": n_features,
            "n_samples": EXPECTED_N_SAMPLES,
            "median_cv": plot.result.medians,
            "n_finite_cv": finite_counts,
            "color_map": plot.color_map,
            "inputs": inputs,
            "artifacts": {
                "svg": str(artifacts.svg),
                "png": str(artifacts.png),
                "legend_svg": str(artifacts.legend_svg),
                "legend_png": str(artifacts.legend_png),
            },
        }
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--qc-states-dir", type=Path, default=Path("results/qc_states"))
    parser.add_argument("--output-dir", type=Path, default=Path("figures/qc/cv"))
    parser.add_argument(
        "--registry-path", type=Path, default=Path("state/color_registry.json")
    )
    parser.add_argument(
        "--provenance-path",
        type=Path,
        default=Path("results/stage3/qc_figures/cv_provenance.json"),
    )
    parser.add_argument("--n-bins", type=int, default=60)
    args = parser.parse_args(argv)
    record = run(
        args.qc_states_dir,
        args.output_dir,
        args.registry_path,
        args.provenance_path,
        args.n_bins,
    )
    for name, fig in record["figures"].items():
        meds = ", ".join(f"{k}={v:.4f}" for k, v in fig["median_cv"].items())
        print(f"{name}: n={fig['n_features']} | {meds}")
        print(f"  colors: {fig['color_map']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
