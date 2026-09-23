"""Stage-4 supplementary DE figures for findings 0004 / 0005 (raloxifene-vs-control).

Renders three figures for claims that no existing figure shows:

  * ``0004-pvalue-hist-protein-paired-by-abundance-tercile`` -- LFQ protein, paired
    design: raw-p histograms by equal-count mean-log2-abundance tercile, Storey pi0
    (lambda 0.5) per panel;
  * ``0004-residual-sd-by-design-protein`` -- LFQ protein: per-protein residual SD
    for the paired / batch / unadjusted designs, medians marked;
  * ``0005-mean-abundance-lfq-vs-spectral`` -- mean log2 abundance, LFQ vs PSM and
    LFQ vs NSAF, on the 1,640 common protein groups (Spearman rho, LOWESS trend).

Each is dual-exported (SVG + 300-DPI PNG) with a separate ``.legend.{svg,png}``. The
two 0004 figures go to ``--da-output-dir``, the 0005 figure to ``--quant-output-dir``.
Every number that also exists in ``summary.json`` (full-set pi0, n, per-design median
residual SD, common-set n) is recomputed and checked against it (fail loud). Provenance
(script hashes, data_version, inputs, params, key numbers, artifact paths) goes to
``--provenance``. Deterministic; no seed consumed.

Run (from the project root):
    ./.venv/bin/python scripts/scratch/fig_de_supplementary.py
"""

from __future__ import annotations

import argparse
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

from analysis_figures.de_supplementary import (  # noqa: E402
    DESIGNS,
    PVALUE_CATEGORY,
    AgreementPanel,
    design_colors,
    plot_abundance_agreement,
    plot_pvalue_by_tercile,
    plot_residual_sd_by_design,
    read_table,
    storey_pi0,
)
from analysis_figures.quant_comparison import load_comparison  # noqa: E402
from common.figures.figure_io import (  # noqa: E402
    FigureArtifacts,
    publication_style,
    save_figure,
)
from common.hashing import sha256_of_file  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "fig-de-supplementary",
    "kind": "figure",
    "provides": [],
    "uses": [
        "analysis_figures.de_supplementary",
        "analysis_figures.quant_comparison",
        "common.figures.figure_io",
        "common.hashing",
    ],
    "seeded_from": None,
    "description": (
        "Supplementary figures for findings 0004/0005: LFQ protein p-value histograms "
        "by abundance tercile, residual SD by design, LFQ vs spectral mean-abundance "
        "agreement; numbers cross-checked against summary.json; provenance JSON."
    ),
}

LOG = logging.getLogger("fig_de_supplementary")

PRIMARY = "paired"
EXPECTED_COMMON_N = 1640
STEMS = {
    "tercile": "0004-pvalue-hist-protein-paired-by-abundance-tercile",
    "design_sd": "0004-residual-sd-by-design-protein",
    "abundance": "0005-mean-abundance-lfq-vs-spectral",
}
PANELS = (
    AgreementPanel("psm_log2", "PSM", "mean log2 PSM count"),
    AgreementPanel("nsaf", "NSAF", "mean log2 NSAF"),
)
PROCESSING_0005 = "LFQ median-normalized log2 · NSAF as-is log2 · PSM log2 unnormalized"
_SUMMARY_ATOL = 6e-5  # summary.json rounds to 4 dp


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    p.add_argument(
        "--de-dir", type=Path, default=Path("results/de/raloxifene-vs-control")
    )
    p.add_argument(
        "--da-output-dir",
        type=Path,
        default=Path("figures/analysis/differential-abundance/raloxifene-vs-control"),
    )
    p.add_argument(
        "--quant-output-dir",
        type=Path,
        default=Path("figures/analysis/quant-comparison/raloxifene-vs-control"),
    )
    p.add_argument(
        "--provenance",
        type=Path,
        default=Path(
            "results/de/raloxifene-vs-control/supplementary_figure_provenance.json"
        ),
    )
    p.add_argument("--registry", type=Path, default=Path("state/color_registry.json"))
    p.add_argument("--n-bins", type=int, default=20)
    p.add_argument("--pi0-lambda", type=float, default=0.5)
    p.add_argument("--bins-per-decade", type=int, default=20)
    p.add_argument("--lowess-frac", type=float, default=0.3)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args(argv)


def _check(label: str, computed: float, reported: object) -> None:
    if not isinstance(reported, int | float):
        raise ValueError(f"summary.json {label}: not a number ({reported!r}).")
    if abs(computed - float(reported)) > _SUMMARY_ATOL:
        raise ValueError(
            f"{label}: computed {computed:.6f} != summary.json {reported}; the figure "
            "would disagree with the recorded analysis."
        )
    LOG.info("check %s: computed %.4f == summary %s", label, computed, reported)


def _artifact_paths(artifacts: FigureArtifacts) -> dict[str, str | None]:
    return {k: (str(v) if v is not None else None) for k, v in vars(artifacts).items()}


def _round(stats: dict[str, dict[str, float | int]]) -> dict[str, dict[str, Any]]:
    return {
        k: {n: (round(v, 4) if isinstance(v, float) else v) for n, v in s.items()}
        for k, s in stats.items()
    }


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args(argv)
    LOG.info("params %s", vars(args))
    summary: dict[str, Any] = json.loads((args.de_dir / "summary.json").read_text())
    if summary.get("primary_design") != PRIMARY:
        raise ValueError(f"primary_design is {summary.get('primary_design')!r}.")
    protein_designs = summary["quantities"]["protein"]["designs"]

    # ---- inputs: LFQ protein, three designs (same feature set) -------------------
    tables = {
        d.key: read_table(
            args.de_dir / f"protein_{d.key}.tsv",
            ["p", "residual_sd", "mean_log2_abundance"],
            "feature",
        )
        for d in DESIGNS
    }
    ref_ids = set(tables[PRIMARY]["feature"])
    for key, frame in tables.items():
        if set(frame["feature"]) != ref_ids:
            raise ValueError(f"protein_{key}.tsv feature set differs from paired.")
        if len(frame) != protein_designs[key]["n_features_tested"]:
            raise ValueError(f"protein_{key}: n {len(frame)} != summary.")
    paired = tables[PRIMARY]
    p = paired["p"].to_numpy(dtype=np.float64)
    abundance = paired["mean_log2_abundance"].to_numpy(dtype=np.float64)
    _check(
        "protein paired pi0 (all)",
        storey_pi0(p, args.pi0_lambda),
        protein_designs[PRIMARY]["storey_pi0_lambda0.5"],
    )
    colors = design_colors(args.registry)
    primary_label = next(d.label for d in DESIGNS if d.key == PRIMARY)
    n_prot = len(paired)

    # ---- inputs: common set (mean abundance joined from <q>_paired.tsv) ----------
    common = load_comparison(args.de_dir, design=PRIMARY, expected_n=EXPECTED_COMMON_N)
    if summary["common_set"]["n_features"] != common.n_features:
        raise ValueError("common-set n disagrees with summary.json.")

    with publication_style():
        fig, legend, tercile_stats = plot_pvalue_by_tercile(
            p,
            abundance,
            color=colors[primary_label],
            title="LFQ protein raw p-value by abundance tercile, raloxifene-d0 vs "
            "control",
            subtitle=(
                f"paired design (condition + pair), moderated t · median-normalized "
                f"log2 · n = {n_prot:,} proteins · equal-count terciles of mean log2 "
                "intensity"
            ),
            abundance_name="mean log2 intensity",
            n_bins=args.n_bins,
            lam=args.pi0_lambda,
        )
        art_tercile = save_figure(
            fig, args.da_output_dir, STEMS["tercile"], legend_fig=legend, dpi=args.dpi
        )

        sds = {
            k: f["residual_sd"].to_numpy(dtype=np.float64) for k, f in tables.items()
        }
        fig, legend, sd_stats = plot_residual_sd_by_design(
            sds,
            colors,
            title="LFQ protein residual SD by design",
            subtitle=(
                f"per-protein unmoderated OLS residual SD · median-normalized log2 · "
                f"n = {n_prot:,} proteins per design"
            ),
            bins_per_decade=args.bins_per_decade,
        )
        art_sd = save_figure(
            fig, args.da_output_dir, STEMS["design_sd"], legend_fig=legend, dpi=args.dpi
        )

        fig, legend, agree_stats = plot_abundance_agreement(
            common,
            PANELS,
            x_key="protein",
            xlabel="mean log2 LFQ intensity (a.u.)",
            title=(
                f"Abundance level (mean over 8 runs), {common.n_features:,} common "
                "protein groups"
            ),
            subtitle=PROCESSING_0005,
            lowess_frac=args.lowess_frac,
        )
        art_abund = save_figure(
            fig,
            args.quant_output_dir,
            STEMS["abundance"],
            legend_fig=legend,
            dpi=args.dpi,
        )

    for d in DESIGNS:
        _check(
            f"protein {d.key} median residual SD",
            float(sd_stats[d.key]["median"]),
            protein_designs[d.key]["median_residual_sd"],
        )
    if sum(int(s["n"]) for s in tercile_stats.values()) != n_prot:
        raise ValueError("tercile sizes do not sum to n.")

    def _script(rel: str) -> dict[str, str]:
        return {"path": rel, "sha256": sha256_of_file(_SCRATCH.parent.parent / rel)}

    input_names = [
        "summary.json",
        "common_set_comparison.tsv",
        *(f"protein_{d.key}.tsv" for d in DESIGNS),
        "nsaf_paired.tsv",
        "psm_log2_paired.tsv",
        "nsaf_batch.tsv",
        "psm_log2_batch.tsv",
    ]
    provenance = {
        "findings": ["0004", "0005"],
        "scripts": {
            "runner": _script("scripts/scratch/fig_de_supplementary.py"),
            "module": _script("scripts/scratch/analysis_figures/de_supplementary.py"),
            "quant_comparison_module": _script(
                "scripts/scratch/analysis_figures/quant_comparison.py"
            ),
            "figure_io": _script("scripts/promoted/common/figures/figure_io.py"),
            "colors": _script("scripts/promoted/common/figures/colors.py"),
        },
        "commit": None,
        "commit_note": "scripts uncommitted at render time (HEAD 416cb10 predates "
        "them); script sha256 pinned instead",
        "data_version": summary["data_version"],
        "de_script": {"path": summary["script"], "sha256": summary["script_sha256"]},
        "inputs": {
            n: {
                "path": str(args.de_dir / n),
                "sha256": sha256_of_file(args.de_dir / n),
            }
            for n in input_names
        },
        "params": {
            "n_bins": args.n_bins,
            "pi0_estimator": f"Storey fixed lambda = {args.pi0_lambda}, clamped <= 1",
            "terciles": "equal-count split of protein_paired mean_log2_abundance "
            "(stable argsort, np.array_split; larger group first)",
            "p_values": "protein_paired.tsv p (moderated t, contrast term)",
            "residual_sd": "per-feature unmoderated OLS residual SD from "
            "protein_<design>.tsv; histogram on log10(SD), "
            f"{args.bins_per_decade} bins per decade",
            "mean_abundance": "mean over 8 runs of the log2 value, "
            "<q>_paired.tsv mean_log2_abundance joined on feature / first_member_id "
            "(analysis_figures.quant_comparison.load_comparison, design-verified)",
            "trend": f"statsmodels LOWESS, frac = {args.lowess_frac}, it = 3",
            "color_category": f"{PVALUE_CATEGORY} (figure-local, persist=False; same "
            "mapping as the 0004 p-value figures)",
            "registry": str(args.registry),
            "dpi": args.dpi,
        },
        "color_map": colors,
        "key_numbers": {
            "protein_paired_pi0_all": round(storey_pi0(p, args.pi0_lambda), 4),
            "tercile": _round(tercile_stats),
            "residual_sd_by_design": _round(sd_stats),
            "mean_abundance_agreement": _round(agree_stats),
        },
        "figures": {
            STEMS["tercile"]: _artifact_paths(art_tercile),
            STEMS["design_sd"]: _artifact_paths(art_sd),
            STEMS["abundance"]: _artifact_paths(art_abund),
        },
    }
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(json.dumps(provenance, indent=2) + "\n")
    LOG.info("key numbers %s", json.dumps(provenance["key_numbers"]))
    LOG.info("wrote %s", args.provenance)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
