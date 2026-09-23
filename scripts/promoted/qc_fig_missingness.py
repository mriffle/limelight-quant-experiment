"""Stage-3 QC figure family *missingness*: completeness-by-batch + MNAR diagnostic.

Loads the prep-once ``raw_linear`` state (``results/qc_states/<level>/raw_linear/``;
NaN = not quantified; all features incl. flagged contaminants) for each requested
level via ``loaders.dataset_io.load_dataset`` and renders
``<output-dir>/missingness-<level>-raw-linear.{svg,png}`` + ``.legend.{svg,png}``
through ``qc_figures.missingness``. Writes ``missingness.provenance.json`` beside the
figures (script hashes, data_version from the prep manifest, params, key numbers).
No stochastic step.

Run (from the project root):
    ./.venv/bin/python scripts/promoted/qc_fig_missingness.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from loaders.dataset_io import load_dataset  # noqa: E402
from qc_figures.missingness import (  # noqa: E402
    save_missingness,
    summarize_completeness,
)

__script_meta__: dict[str, object] = {
    "task": "qc-fig-missingness",
    "kind": "figure",
    "provides": [],
    "uses": ["loaders.dataset_io", "qc_figures.missingness"],
    "seeded_from": None,
    "description": (
        "Stage-3 QC missingness family: for protein and peptide raw_linear states, "
        "renders a completeness curve overlaid by batch (+ all runs) and an MNAR "
        "panel (detection rate vs mean log2 intensity, hexbin + binned median + "
        "Pearson r) to figures/qc/missingness/, with a provenance JSON."
    ),
}

_LEVEL_LABEL = {"protein": "Protein", "peptide": "Peptide"}
_STATE = "raw_linear"


def _sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--qc-states-dir", type=Path, default=Path("results/qc_states"))
    p.add_argument("--output-dir", type=Path, default=Path("figures/qc/missingness"))
    p.add_argument("--levels", nargs="+", default=["protein", "peptide"])
    p.add_argument("--color-by", default="batch")
    p.add_argument("--min-intensity", type=float, default=0.0)
    p.add_argument("--registry", type=Path, default=Path("state/color_registry.json"))
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    manifest_path = args.qc_states_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    data_version = manifest.get("data_version")
    if not data_version:
        raise ValueError(f"{manifest_path} has no data_version.")

    figures: dict[str, Any] = {}
    for level in args.levels:
        if level not in _LEVEL_LABEL:
            raise ValueError(f"unknown level {level!r}; expected {list(_LEVEL_LABEL)}")
        state_dir = args.qc_states_dir / level / _STATE
        dataset = load_dataset(state_dir)
        n_samples, n_features = dataset.abundances.shape
        base = f"missingness-{level}-raw-linear"
        title = (
            f"{_LEVEL_LABEL[level]} missingness — raw, linear "
            f"({n_samples} runs, {n_features:,} features)"
        )
        artifacts, plot = save_missingness(
            dataset,
            args.output_dir,
            base,
            color_by=args.color_by,
            min_intensity=args.min_intensity,
            title=title,
            legend_title=args.color_by,
            registry_path=args.registry,
            dpi=args.dpi,
        )
        summary = summarize_completeness(
            dataset,
            fractions=(1.0, 0.5),
            group_by=args.color_by,
            min_intensity=args.min_intensity,
        )
        n_contam = int(dataset.feature_metadata["is_contaminant"].sum())
        figures[base] = {
            "input": str(state_dir),
            "n_samples": int(n_samples),
            "n_features": int(n_features),
            "n_contaminant_features_included": n_contam,
            "n_never_detected": int((plot.result.feature_detection_rate == 0).sum()),
            "mnar_pearson_r": plot.result.mnar_correlation,
            "n_mnar_features": plot.result.n_mnar_features,
            "features_retained": summary,
            "color_map": plot.color_map,
            "artifacts": {
                k: (str(v) if v is not None else None)
                for k, v in vars(artifacts).items()
            },
        }
        print(json.dumps({base: figures[base]}, indent=2))

    module_path = _SCRATCH / "qc_figures" / "missingness.py"
    provenance = {
        "scripts": {
            "runner": {
                "path": "scripts/promoted/qc_fig_missingness.py",
                "sha256": _sha256(Path(__file__).resolve()),
            },
            "module": {
                "path": "scripts/promoted/qc_figures/missingness.py",
                "sha256": _sha256(module_path),
                "seeded_from": {"template": "missingness", "version": "0.1"},
            },
        },
        "commit": None,
        "commit_note": "project is not a git repository; script sha256 pinned instead",
        "data_version": data_version,
        "params": {
            "state": _STATE,
            "levels": args.levels,
            "color_by": args.color_by,
            "min_intensity": args.min_intensity,
            "registry": str(args.registry),
            "dpi": args.dpi,
            "detection_rule": "finite and > min_intensity",
        },
        "figures": figures,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "missingness.provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
