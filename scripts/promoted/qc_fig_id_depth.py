"""Stage-3 QC runner: identification depth (detected features per run), raw linear.

Loads the prep-once processing state ``raw_linear`` for the protein and peptide levels
(``results/qc_states/<level>/raw_linear/``, written by ``qc_prep.py``) and renders ONE
stacked figure -- protein-group depth over peptide depth -- via
:mod:`qc_figures.id_depth`, bars in acquisition order (batch, then run position),
colored by ``condition`` through ``state/color_registry.json``, with a dotted divider
between the two acquisition batches and a dashed reference median over the
experimental samples (all 8 runs here: the study has no pooled/control-reference
samples). ``raw_linear`` is the full matrix incl. flagged contaminants, NaN = missing,
so a feature counts as detected when finite and > 0.

Writes ``<out-dir>/id-depth-raw-linear.{svg,png}`` + ``.legend.{svg,png}`` (default
``figures/qc/id-depth``) and a provenance JSON (script sha256 + git commit, data
version, inputs, params, the per-sample counts).

Usage::

    .venv/bin/python scripts/promoted/qc_fig_id_depth.py [--out-dir DIR]
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
_PROJECT_ROOT = _SCRATCH.parent.parent
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from common.hashing import sha256_of_file  # noqa: E402
from loaders.data_loading import Dataset  # noqa: E402
from loaders.dataset_io import load_dataset  # noqa: E402
from qc_figures.id_depth import save_id_depth  # noqa: E402

LOGGER = logging.getLogger("qc_fig_id_depth")

STATE = "raw_linear"
BASE_NAME = "id-depth-raw-linear"
LEVELS: dict[str, str] = {"Protein groups": "protein", "Peptides": "peptide"}
YLABELS: dict[str, str] = {
    "Protein groups": "Detected protein groups",
    "Peptides": "Detected peptides",
}
COLOR_BY = "condition"
DIVIDER_BY = "batch"
REFERENCE_ROLE = "experimental"
TITLE = "Identification depth per run (raw, linear)"
ORDER_COLUMNS = ("batch", "run_position_within_batch")


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--states-root",
        type=Path,
        default=_PROJECT_ROOT / "results" / "qc_states",
        help="Root of the prep-once processing states.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=_PROJECT_ROOT / "figures" / "qc" / "id-depth",
        help="Target figure directory (structured figures/ layout).",
    )
    parser.add_argument(
        "--provenance",
        type=Path,
        default=_PROJECT_ROOT
        / "results"
        / "stage3"
        / "qc_figures"
        / f"{BASE_NAME}.provenance.json",
        help="Where to write this figure's provenance JSON.",
    )
    parser.add_argument(
        "--registry",
        type=Path,
        default=_PROJECT_ROOT / "state" / "color_registry.json",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args(argv)


def _verify_acquisition_order(ds: Dataset, level: str) -> None:
    """Fail loud unless rows are already sorted by batch then run position."""
    meta = ds.metadata
    for col in ORDER_COLUMNS:
        if col not in meta.columns:
            raise ValueError(f"{level}: metadata lacks order column {col!r}.")
    batch = meta["batch"].astype(str).to_numpy()
    pos = meta["run_position_within_batch"].astype(int).to_numpy()
    keys = list(zip(batch, pos, strict=True))
    if keys != sorted(keys):
        raise ValueError(
            f"{level}: samples are not in acquisition order (batch, run position): "
            f"{keys}."
        )


def _git_commit(root: Path) -> str | None:
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


def _rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(_PROJECT_ROOT.resolve()))
    except ValueError:
        return str(resolved)


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args(argv)

    datasets: dict[str, Dataset] = {}
    input_dirs: dict[str, Path] = {}
    for label, level in LEVELS.items():
        path = args.states_root / level / STATE
        ds = load_dataset(path)
        _verify_acquisition_order(ds, level)
        datasets[label] = ds
        input_dirs[level] = path
        LOGGER.info("Loaded %s %s: %s", level, STATE, ds.abundances.shape)

    first = next(iter(datasets.values()))
    roles = first.metadata["sample_role"].astype(str).to_numpy()
    reference_mask = roles == REFERENCE_ROLE
    if not reference_mask.any():
        raise ValueError(f"No {REFERENCE_ROLE!r} samples for the reference median.")
    LOGGER.info(
        "Reference median over %d/%d %s samples",
        int(reference_mask.sum()),
        reference_mask.size,
        REFERENCE_ROLE,
    )

    artifacts, result, color_map = save_id_depth(
        datasets,
        args.out_dir,
        BASE_NAME,
        color_by=COLOR_BY,
        reference_mask=reference_mask,
        ylabel=YLABELS,
        title=TITLE,
        legend_title="Condition",
        annotate_counts=True,
        divider_by=DIVIDER_BY,
        registry_path=args.registry,
        persist_colors=False,
        dpi=args.dpi,
    )

    # Non-contaminant counts are recorded alongside (not plotted) so the finding can
    # state how much of the depth is contaminant features.
    noncontam: dict[str, list[int]] = {}
    for label, ds in datasets.items():
        keep = ~ds.feature_metadata["is_contaminant"].astype(bool).to_numpy()
        a = np.asarray(ds.abundances, dtype=float)[:, keep]
        noncontam[label] = [int(v) for v in (np.isfinite(a) & (a > 0)).sum(axis=1)]

    manifest_path = args.states_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    script = Path(__file__).resolve()
    module = _SCRATCH / "qc_figures" / "id_depth.py"
    record: dict[str, object] = {
        "figure": BASE_NAME,
        "svg": _rel(artifacts.svg),
        "png": _rel(artifacts.png),
        "legend_svg": _rel(artifacts.legend_svg) if artifacts.legend_svg else None,
        "legend_png": _rel(artifacts.legend_png) if artifacts.legend_png else None,
        "script": {
            "runner": _rel(script),
            "runner_sha256": sha256_of_file(script),
            "module": _rel(module),
            "module_sha256": sha256_of_file(module),
            "git_commit": _git_commit(_PROJECT_ROOT),
            "seeded_from": "lib/figures/id_depth.py (id-depth v0.1)",
        },
        "data_version": manifest.get("data_version"),
        "processing_state": STATE,
        "inputs": {lvl: _rel(p) for lvl, p in input_dirs.items()},
        "color_registry_sha256": sha256_of_file(args.registry),
        "params": {
            "color_by": COLOR_BY,
            "divider_by": DIVIDER_BY,
            "reference_subset": f"sample_role == {REFERENCE_ROLE!r}",
            "min_intensity": 0.0,
            "dpi": args.dpi,
            "order": list(ORDER_COLUMNS),
        },
        "color_map": color_map,
        "sample_ids": [str(s) for s in result.sample_ids],
        "counts": {k: [int(v) for v in c] for k, c in result.counts.items()},
        "reference_median": result.reference_median,
        "counts_non_contaminant": noncontam,
    }
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    LOGGER.info(
        "Wrote %s, %s; provenance %s", artifacts.svg, artifacts.png, args.provenance
    )
    for label, c in result.counts.items():
        LOGGER.info(
            "%s: %s (ref median %.0f)",
            label,
            dict(zip(record["sample_ids"], c.tolist(), strict=True)),  # type: ignore[call-overload]
            result.reference_median[label],
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
