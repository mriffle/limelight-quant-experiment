"""Stage-3 QC: PCA processing-state series (family ``pca``), protein + peptide.

Loads the prep-once states ``results/qc_states/<level>/{raw_log, normalized_log,
batch_corrected_log}`` (log2; non-contaminant features complete in all 8 runs) and,
per level, renders two small-multiple figures (one PC1/PC2 panel per state):

  * ``pca-by-batch-<level>-raw-normalized-batchcorrected-log2``     -- batch coloring
  * ``pca-by-condition-<level>-raw-normalized-batchcorrected-log2`` -- condition
    coloring (this study has no control samples, so the sample-class coloring the
    convention asks for is replaced by the condition coloring, per the dispatch)

each with a separate legend image, under ``figures/qc/pca/``. Writes a per-figure
provenance JSON (script sha256 + git commit, data_version, input hashes, params,
per-state % variance and PC1/PC2 scores).

Deterministic (full-SVD PCA; greedy label placement); no seed is consumed.

Run:
    ./.venv/bin/python scripts/promoted/qc_fig_pca.py
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from common.hashing import sha256_of_file  # noqa: E402
from loaders.dataset_io import load_dataset  # noqa: E402
from qc_figures.pca import StatePanel, save_pca_state_series  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "qc-fig-pca",
    "kind": "figure",
    "provides": [],
    "uses": ["loaders.dataset_io", "qc_figures.pca", "common.hashing"],
    "seeded_from": None,
    "description": (
        "Stage-3 QC PCA family: per level (protein, peptide) a batch-colored and a "
        "condition-colored PC1/PC2 small-multiple series over raw -> median-"
        "normalized -> ComBat batch-only preview (log2), sample-id labeled; writes "
        "figures/qc/pca/ + provenance JSON."
    ),
}

LEVELS: tuple[str, ...] = ("protein", "peptide")
# (state directory, panel label) in left-to-right order.
STATES: tuple[tuple[str, str], ...] = (
    ("raw_log", "Raw (log2)"),
    ("normalized_log", "Median-normalized (log2)"),
    ("batch_corrected_log", "Batch-corrected, ComBat batch-only (log2)"),
)
COLORINGS: tuple[tuple[str, str], ...] = (
    ("batch", "Batch"),
    ("condition", "Condition"),
)
LABEL_COLUMN = "sample_id"
STEM = "pca-by-{coloring}-{level}-raw-normalized-batchcorrected-log2"


def _git_commit(root: Path) -> str | None:
    """HEAD commit of the project repo, or ``None`` if not a git repository."""
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def _rel(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _state_hashes(state_dir: Path, root: Path) -> dict[str, str]:
    files = ("abundances.npy", "metadata.parquet", "dataset.json")
    return {_rel(state_dir / f, root): sha256_of_file(state_dir / f) for f in files}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--qc-states-dir", type=Path, default=Path("results/qc_states"))
    parser.add_argument("--output-dir", type=Path, default=Path("figures/qc/pca"))
    parser.add_argument(
        "--registry", type=Path, default=Path("state/color_registry.json")
    )
    parser.add_argument(
        "--provenance-file",
        type=Path,
        default=None,
        help="Default: results/stage3/figure_provenance/<stem-prefix>pca.json.",
    )
    parser.add_argument(
        "--stem-prefix",
        default="",
        help="Prefix for figure stems + provenance file (e.g. '0003-'); default none.",
    )
    parser.add_argument(
        "--levels", nargs="+", default=list(LEVELS), choices=list(LEVELS)
    )
    parser.add_argument(
        "--colorings",
        nargs="+",
        default=[c for c, _ in COLORINGS],
        choices=[c for c, _ in COLORINGS],
    )
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.provenance_file is None:
        args.provenance_file = Path(
            f"results/stage3/figure_provenance/{args.stem_prefix}pca.json"
        )
    root: Path = args.project_root
    manifest_path: Path = args.qc_states_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    data_version = manifest["data_version"]
    script = Path(__file__).resolve()
    module = _SCRATCH / "qc_figures" / "pca.py"
    code = {
        "script": {
            "path": _rel(script, root),
            "sha256": sha256_of_file(script),
            "module": _rel(module, root),
            "module_sha256": sha256_of_file(module),
            "git_commit": _git_commit(root),
            "seeded_from": "lib/figures/pca.py (pca-plot v0.3)",
        },
        "registry_sha256": sha256_of_file(args.registry),
    }

    records: dict[str, object] = {}
    for level in args.levels:
        panels = []
        inputs: dict[str, str] = {}
        for state_dir_name, label in STATES:
            state_dir = args.qc_states_dir / level / state_dir_name
            panels.append(StatePanel(label=label, dataset=load_dataset(state_dir)))
            inputs.update(_state_hashes(state_dir, root))
        n_features = panels[0].dataset.abundances.shape[1]
        n_samples = panels[0].dataset.abundances.shape[0]
        for coloring, legend_title in COLORINGS:
            if coloring not in args.colorings:
                continue
            stem = args.stem_prefix + STEM.format(coloring=coloring, level=level)
            title = (
                f"PCA of {level}s by {coloring} · {n_features:,} complete "
                f"{level}s · n = {n_samples}"
            )
            artifacts, plot = save_pca_state_series(
                panels,
                coloring,
                args.output_dir,
                stem,
                label_by=LABEL_COLUMN,
                title=title,
                legend_title=legend_title,
                registry_path=args.registry,
                dpi=args.dpi,
            )
            sample_ids = panels[0].dataset.metadata[LABEL_COLUMN].astype(str).tolist()
            summary = {
                label: {
                    "pc1_pct": round(float(res.explained_variance_ratio[0]) * 100, 2),
                    "pc2_pct": round(float(res.explained_variance_ratio[1]) * 100, 2),
                    "n_features": res.n_features,
                    "scores": {
                        sid: [round(float(v), 3) for v in res.scores[i, :2]]
                        for i, sid in enumerate(sample_ids)
                    },
                }
                for label, res in plot.results.items()
            }
            records[stem] = {
                "svg": _rel(artifacts.svg, root),
                "png": _rel(artifacts.png, root),
                "legend_svg": _rel(artifacts.legend_svg, root)
                if artifacts.legend_svg
                else None,
                "legend_png": _rel(artifacts.legend_png, root)
                if artifacts.legend_png
                else None,
                **code,
                "data_version": data_version,
                "inputs": inputs,
                "params": {
                    "level": level,
                    "color_by": coloring,
                    "label_by": LABEL_COLUMN,
                    "states": [s for s, _ in STATES],
                    "standardize": True,
                    "n_components": 2,
                    "svd_solver": "full",
                    "dpi": args.dpi,
                },
                "color_map": plot.color_map,
                "processing_state": "raw_log | normalized_log | batch_corrected_log",
                "summary": summary,
            }
            print(f"wrote {_rel(artifacts.png, root)}")

    args.provenance_file.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_file.write_text(json.dumps(records, indent=2) + "\n")
    print(f"wrote {_rel(args.provenance_file, root)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
