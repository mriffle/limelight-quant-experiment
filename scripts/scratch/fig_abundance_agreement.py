"""Stage-4 figures for findings 0010 (LFQ vs spectral abundance agreement) and 0011
(contaminant peptide sharing), family ``quant-comparison``.

Renders nine figures from the tables in ``results/abundance-agreement/`` (see
``analysis_figures.abundance_agreement_figs``); every on-canvas number that
``summary.json`` records is cross-checked against it (fail loud). Each figure is
dual-exported (SVG + 300-DPI PNG) with a separate ``.legend.{svg,png}``. Provenance
(script hashes, data_version, input hashes, params, key numbers, artifact paths)
goes to ``--provenance``. Deterministic (the only randomness is a fixed-seed strip
jitter in the no-LFQ figure).

Run (from the project root):
    ./.venv/bin/python scripts/scratch/fig_abundance_agreement.py
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis_figures import abundance_agreement_figs as aaf  # noqa: E402
from common.figures.figure_io import (  # noqa: E402
    FigureArtifacts,
    publication_style,
    save_figure,
)
from common.hashing import sha256_of_file  # noqa: E402
from loaders.peptide_loader import load_peptide_dataset  # noqa: E402
from loaders.protein_loader import load_protein_dataset  # noqa: E402
from loaders.samples import read_samples  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "fig-abundance-agreement",
    "kind": "figure",
    "provides": [],
    "uses": [
        "analysis_figures.abundance_agreement_figs",
        "common.figures.figure_io",
        "common.hashing",
        "loaders.samples",
        "loaders.protein_loader",
        "loaders.peptide_loader",
    ],
    "seeded_from": None,
    "description": (
        "Findings 0010/0011 figures (LFQ vs PSM/NSAF agreement; contaminant peptide "
        "sharing) from results/abundance-agreement tables, numbers cross-checked "
        "against summary.json, provenance JSON beside the results."
    ),
}

LOG = logging.getLogger("fig_abundance_agreement")
EXPECTED_DATA_VERSION_PREFIX = "sha256:bc6b73d3"


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    p.add_argument(
        "--result-dir", type=Path, default=Path("results/abundance-agreement")
    )
    p.add_argument("--samples", type=Path, default=Path("results/metadata/samples.tsv"))
    p.add_argument(
        "--output-dir-0010",
        type=Path,
        default=Path("figures/analysis/quant-comparison/abundance-agreement"),
    )
    p.add_argument(
        "--output-dir-0011",
        type=Path,
        default=Path("figures/analysis/quant-comparison/contaminant-sharing"),
    )
    p.add_argument(
        "--provenance",
        type=Path,
        default=Path("results/abundance-agreement/figure_provenance.json"),
    )
    p.add_argument(
        "--protein-quants", type=Path, default=Path("data/protein-quants.tsv")
    )
    p.add_argument(
        "--peptide-quants", type=Path, default=Path("data/peptide-quants.tsv")
    )
    p.add_argument(
        "--protein-dump",
        type=Path,
        default=Path("data/protein-limelight-table-dump.txt"),
    )
    p.add_argument(
        "--contaminant-table",
        type=Path,
        default=Path("results/abundance-agreement/contaminant_sharing_corrected.tsv"),
    )
    p.add_argument("--registry", type=Path, default=Path("state/color_registry.json"))
    p.add_argument("--top-no-lfq", type=int, default=20)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args(argv)


def _git_head() -> str | None:
    """HEAD commit of the project repository (None outside git)."""
    out = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False
    )
    return out.stdout.strip() if out.returncode == 0 else None


def _committed(paths: list[str]) -> dict[str, bool]:
    """Whether each path is tracked and unmodified at render time."""
    out: dict[str, bool] = {}
    for path in paths:
        st = subprocess.run(
            ["git", "status", "--porcelain", "--", path],
            capture_output=True,
            text=True,
            check=False,
        )
        out[path] = st.returncode == 0 and st.stdout.strip() == ""
    return out


def _paths(art: FigureArtifacts) -> dict[str, str | None]:
    return {k: (str(v) if v is not None else None) for k, v in vars(art).items()}


def _contaminant_table(args: argparse.Namespace, data: aaf.AgreementData) -> Any:
    """Recompute contaminant sharing from peptide Protein Groups (verified loaders)."""
    prot = load_protein_dataset(
        args.protein_quants, args.samples, protein_limelight_file=args.protein_dump
    )
    fm = prot.dataset.feature_metadata
    ids = fm["protein_group"].astype(str).tolist()
    is_cont = fm["is_contaminant"].to_numpy(dtype=bool)
    cont_ids = frozenset(i for i, c in zip(ids, is_cont, strict=True) if c)
    if len(cont_ids) != data.summary["inputs"]["n_contaminants_excluded"]:
        raise ValueError(f"{len(cont_ids)} contaminants != analysis count")
    ab = np.asarray(prot.dataset.abundances, dtype=np.float64)
    if ab.shape != (len(data.samples), len(ids)):
        raise ValueError(
            f"protein matrix shape {ab.shape} (expected samples x features)"
        )
    lfq_runs = {
        i: int(np.isfinite(ab[:, j]).sum()) for j, i in enumerate(ids) if i in cont_ids
    }
    ref = data.t("contaminants").set_index("protein_id")["n_runs_lfq"]
    for cid, n in lfq_runs.items():
        if int(ref.loc[cid]) != n:
            raise ValueError(f"{cid}: LFQ runs {n} != contaminants.tsv {ref.loc[cid]}")
    pep = load_peptide_dataset(
        args.peptide_quants, args.samples, contaminant_ids=cont_ids
    )
    pab = np.asarray(pep.dataset.abundances, dtype=np.float64)
    quantified = np.isfinite(pab).any(axis=0)
    groups = pep.dataset.feature_metadata["protein_groups"].astype(str).tolist()
    return aaf.contaminant_sharing_table(
        groups, quantified, cont_ids, lfq_runs, data.t("no_lfq_proteins")
    )


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args(argv)
    LOG.info("params %s", {k: str(v) for k, v in vars(args).items()})
    order = list(read_samples(args.samples).sample_id)
    data = aaf.load_agreement(args.result_dir, order)
    if not str(data.summary["data_version"]).startswith(EXPECTED_DATA_VERSION_PREFIX):
        raise ValueError(f"unexpected data_version {data.summary['data_version']}")
    qcol = aaf.registry_colors("quantity", ["LFQ", "NSAF", "PSM"], args.registry)
    bcol = aaf.registry_colors(
        "batch", sorted(set(data.batch_of.values())), args.registry
    )

    contam = _contaminant_table(args, data)
    args.contaminant_table.parent.mkdir(parents=True, exist_ok=True)
    contam.to_csv(args.contaminant_table, sep="\t", index=False)
    LOG.info("wrote %s (%d proteins)", args.contaminant_table, len(contam))

    jobs: list[tuple[str, Path, Callable[[], aaf.Result]]] = [
        (
            "0010-density-lfq-vs-psm-nsaf",
            args.output_dir_0010,
            lambda: aaf.plot_density(data, qcol),
        ),
        (
            "0010-per-sample-correlations",
            args.output_dir_0010,
            lambda: aaf.plot_per_sample(data, bcol),
        ),
        (
            "0010-slope-estimators",
            args.output_dir_0010,
            lambda: aaf.plot_slopes(data, qcol),
        ),
        (
            "0010-agreement-by-psm-count",
            args.output_dir_0010,
            lambda: aaf.plot_by_psm_count(data, qcol),
        ),
        (
            "0010-residual-drivers",
            args.output_dir_0010,
            lambda: aaf.plot_residual_drivers(data, qcol),
        ),
        (
            "0010-detection-classes",
            args.output_dir_0010,
            lambda: aaf.plot_detection_classes(data, args.registry),
        ),
        (
            "0010-no-lfq-reasons",
            args.output_dir_0010,
            lambda: aaf.plot_no_lfq_reasons(data, args.registry, args.top_no_lfq),
        ),
        (
            "0010-within-protein-tracking",
            args.output_dir_0010,
            lambda: aaf.plot_within_protein(data, qcol),
        ),
        (
            "0011-contaminant-sharing",
            args.output_dir_0011,
            lambda: aaf.plot_contaminant_sharing(data, qcol, contam),
        ),
    ]
    figures: dict[str, Any] = {}
    numbers: dict[str, Any] = {}
    with publication_style():
        for stem, out, make in jobs:
            fig, legend, nums = make()
            art = save_figure(fig, out, stem, legend_fig=legend, dpi=args.dpi)
            figures[stem] = _paths(art)
            numbers[stem] = nums
            LOG.info("wrote %s", art.png)

    module = _SCRATCH / "analysis_figures" / "abundance_agreement_figs.py"
    module_rel = Path("scripts/scratch/analysis_figures/abundance_agreement_figs.py")
    provenance = {
        "findings": ["0010", "0011"],
        "scripts": {
            "runner": {
                "path": "scripts/scratch/fig_abundance_agreement.py",
                "sha256": sha256_of_file(Path(__file__).resolve()),
            },
            "module": {
                "path": "scripts/scratch/analysis_figures/abundance_agreement_figs.py",
                "sha256": sha256_of_file(module),
            },
            "figure_io": {
                "path": "scripts/promoted/common/figures/figure_io.py",
                "sha256": sha256_of_file(
                    _PROMOTED / "common" / "figures" / "figure_io.py"
                ),
            },
            "colors": {
                "path": "scripts/promoted/common/figures/colors.py",
                "sha256": sha256_of_file(
                    _PROMOTED / "common" / "figures" / "colors.py"
                ),
            },
        },
        "commit": _git_head(),
        "scripts_committed": _committed(
            ["scripts/scratch/fig_abundance_agreement.py", str(module_rel)]
        ),
        "data_version": data.summary["data_version"],
        "analysis_script": {
            "path": data.summary["script"],
            "sha256": data.summary["script_sha256"],
            "module": data.summary["module"],
            "module_sha256": data.summary["module_sha256"],
        },
        "inputs": {
            name: {
                "path": str(args.result_dir / name),
                "sha256": sha256_of_file(args.result_dir / name),
            }
            for name in ["summary.json", *data.summary["tables"]]
        }
        | {
            str(path): {"path": str(path), "sha256": sha256_of_file(path)}
            for path in (args.protein_quants, args.peptide_quants, args.protein_dump)
        },
        "outputs": {
            "contaminant_table": {
                "path": str(args.contaminant_table),
                "sha256": sha256_of_file(args.contaminant_table),
            }
        },
        "params": {
            "result_dir": str(args.result_dir),
            "samples": str(args.samples),
            "registry": str(args.registry),
            "top_no_lfq": args.top_no_lfq,
            "dpi": args.dpi,
            "contaminant_rule": "refined 33-entry set (protein_loader."
            "refine_contaminants); a peptide row becomes unique "
            "to real protein P if P is its only non-contaminant "
            "member and it has >= 1 contaminant member",
            "registry_read": ["quantity", "batch"],
            "figure_local_colors_not_persisted": ["detection_class", "no_lfq_reason"],
            "fit_line_colors": {
                "LOWESS": aaf.FIT_LOWESS,
                "SMA": aaf.FIT_LINE,
                "partial_fit": aaf.FIT_BLACK,
            },
            "density": "hexbin gridsize (60, 42), log color, shared norm across panels",
            "residual_drivers": "component-plus-residual: residual minus the other 3 "
            "terms' fitted contributions (centered); 4-term OLS "
            "refit here and matched to residual_regression.tsv",
            "acquisition_order": order,
        },
        "color_map": {"quantity": qcol, "batch": bcol},
        "key_numbers": numbers,
        "figures": figures,
    }
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(json.dumps(provenance, indent=2, default=str) + "\n")
    LOG.info("wrote %s", args.provenance)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
