"""Stage-4 figure family *quant-comparison*, label ``raloxifene-vs-control``.

Renders the three finding-0005 figures comparing LFQ / NSAF / PSM on the 1,640
common protein groups under the PRIMARY paired design (condition + candidate_pair):

  * ``0005-log2fc-scatter-lfq-nsaf-psm``   -- pairwise per-protein log2FC scatter;
  * ``0005-residual-sd-by-quantity``       -- residual SD + 95% CI half-width
                                              distributions, medians marked;
  * ``0005-residual-sd-vs-abundance``      -- residual SD vs mean log2 abundance.

Each is dual-exported (SVG + 300-DPI PNG) with a separate ``.legend.{svg,png}`` to
``--output-dir``. Every on-canvas statistic is recomputed from the tables and checked
against ``summary.json`` (fail loud on disagreement). Provenance (script hashes,
data_version, params, key numbers, artifact paths) goes to ``--provenance``.
Deterministic; no seed consumed.

Run (from the project root):
    ./.venv/bin/python scripts/scratch/fig_quant_comparison.py
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

from analysis_figures.quant_comparison import (  # noqa: E402
    PAIRS,
    QUANTITIES,
    ComparisonData,
    load_comparison,
    plot_log2fc_scatter,
    plot_precision_distributions,
    plot_sd_vs_abundance,
    quantity_colors,
)
from common.figures.figure_io import (  # noqa: E402
    FigureArtifacts,
    publication_style,
    save_figure,
)
from common.hashing import sha256_of_file  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "fig-quant-comparison",
    "kind": "figure",
    "provides": [],
    "uses": [
        "analysis_figures.quant_comparison",
        "common.figures.figure_io",
        "common.hashing",
    ],
    "seeded_from": None,
    "description": (
        "Finding-0005 quant-comparison figures (LFQ vs NSAF vs PSM, 1,640 common "
        "groups, paired design) to figures/analysis/quant-comparison/"
        "raloxifene-vs-control/, numbers cross-checked against summary.json, "
        "provenance JSON beside the DE results."
    ),
}

LOG = logging.getLogger("fig_quant_comparison")

FINDING_ID = "0005"
DESIGN = "paired"
EXPECTED_N = 1640
STEMS = {
    "scatter": f"{FINDING_ID}-log2fc-scatter-lfq-nsaf-psm",
    "distributions": f"{FINDING_ID}-residual-sd-by-quantity",
    "abundance": f"{FINDING_ID}-residual-sd-vs-abundance",
}
PROCESSING = "LFQ median-normalized log2 · NSAF as-is log2 · PSM log2 unnormalized"
# summary.json stores rounded values; allow the rounding plus float noise.
_SUMMARY_ATOL = 6e-5


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    p.add_argument(
        "--de-dir", type=Path, default=Path("results/de/raloxifene-vs-control")
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figures/analysis/quant-comparison/raloxifene-vs-control"),
    )
    p.add_argument(
        "--provenance",
        type=Path,
        default=Path(
            "results/de/raloxifene-vs-control/quant_comparison_figure_provenance.json"
        ),
    )
    p.add_argument("--registry", type=Path, default=Path("state/color_registry.json"))
    p.add_argument("--n-bins", type=int, default=10)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args(argv)


def _check(label: str, computed: float, reported: object) -> None:
    if not isinstance(reported, int | float):
        raise ValueError(f"summary.json {label}: not a number ({reported!r}).")
    if abs(computed - float(reported)) > _SUMMARY_ATOL:
        raise ValueError(
            f"{label}: computed {computed:.6f} != summary.json {reported} -- the "
            f"figure would disagree with the recorded analysis."
        )
    LOG.info("check %s: computed %.4f == summary %s", label, computed, reported)


def _verify_summary(
    summary: dict[str, Any],
    data: ComparisonData,
    rhos: dict[str, float],
    medians: dict[str, dict[str, float]],
) -> None:
    """Fail loud unless design, n and every on-canvas number match summary.json."""
    if summary.get("primary_design") != DESIGN:
        raise ValueError(f"primary_design is {summary.get('primary_design')!r}.")
    if summary["designs"][DESIGN] != ["candidate_pair"]:
        raise ValueError(f"{DESIGN} covariates: {summary['designs'][DESIGN]!r}.")
    common = summary["common_set"]
    if common["n_features"] != data.n_features or data.n_features != EXPECTED_N:
        raise ValueError(
            f"common-set n: summary {common['n_features']}, table "
            f"{data.n_features}, expected {EXPECTED_N}."
        )
    if "primary design" not in str(common.get("table", "")):
        raise ValueError("summary.json does not label the common table as primary.")
    labels = {q.key: q.label for q in QUANTITIES}
    for row in common["quantity_concordance"][DESIGN]:
        if row["n_features"] != data.n_features:
            raise ValueError(f"concordance n mismatch: {row}")
        _check(
            f"spearman {row['a']}-{row['b']}",
            rhos[f"{labels[row['a']]}-{labels[row['b']]}"],
            row["spearman_log2fc"],
        )
    for q in QUANTITIES:
        block = common["quantities"][q.key][DESIGN]
        if block["n_features_tested"] != data.n_features:
            raise ValueError(f"{q.key}: n_features_tested {block['n_features_tested']}")
        _check(
            f"{q.key} median residual SD",
            medians["residual_sd"][q.label],
            block["median_residual_sd"],
        )
        _check(
            f"{q.key} median CI half-width",
            medians["ci_half_width"][q.label],
            block["median_ci95_half_width"],
        )


def _artifact_paths(artifacts: FigureArtifacts) -> dict[str, str | None]:
    return {k: (str(v) if v is not None else None) for k, v in vars(artifacts).items()}


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args(argv)
    LOG.info("params %s", vars(args))
    summary_path = args.de_dir / "summary.json"
    summary: dict[str, Any] = json.loads(summary_path.read_text())

    data = load_comparison(args.de_dir, design=DESIGN, expected_n=EXPECTED_N)
    LOG.info(
        "loaded %d common-set groups x %d quantities (%s design verified)",
        data.n_features,
        len(QUANTITIES),
        data.design,
    )
    colors = quantity_colors(args.registry)
    subtitle = (
        f"paired design (condition + pair) · {data.n_features:,} common protein "
        f"groups · {PROCESSING}"
    )

    with publication_style():
        fig, legend, rhos = plot_log2fc_scatter(
            data,
            title="Per-protein log2 fold change, raloxifene \u2212 control",
            subtitle=subtitle,
        )
        art_scatter = save_figure(
            fig, args.output_dir, STEMS["scatter"], legend_fig=legend, dpi=args.dpi
        )
        fig, legend, medians = plot_precision_distributions(
            data,
            colors,
            title="Per-protein precision by quantity",
            subtitle=subtitle,
        )
        art_dist = save_figure(
            fig,
            args.output_dir,
            STEMS["distributions"],
            legend_fig=legend,
            dpi=args.dpi,
        )
        fig, legend, trend = plot_sd_vs_abundance(
            data,
            colors,
            title="Residual SD vs mean abundance, by quantity",
            subtitle=subtitle,
            n_bins=args.n_bins,
        )
        art_abund = save_figure(
            fig, args.output_dir, STEMS["abundance"], legend_fig=legend, dpi=args.dpi
        )

    _verify_summary(summary, data, rhos, medians)

    zero_sd = {
        q.label: int(np.sum(data.column(q.key, "residual_sd") < 1e-10))
        for q in QUANTITIES
    }
    module = _SCRATCH / "analysis_figures" / "quant_comparison.py"
    provenance = {
        "finding": FINDING_ID,
        "scripts": {
            "runner": {
                "path": "scripts/scratch/fig_quant_comparison.py",
                "sha256": sha256_of_file(Path(__file__).resolve()),
            },
            "module": {
                "path": "scripts/scratch/analysis_figures/quant_comparison.py",
                "sha256": sha256_of_file(module),
            },
            "figure_io": {
                "path": "scripts/promoted/common/figures/figure_io.py",
                "sha256": sha256_of_file(
                    _PROMOTED / "common" / "figures" / "figure_io.py"
                ),
            },
        },
        "commit": None,
        "commit_note": "project is not a git repository; script sha256 pinned instead",
        "data_version": summary["data_version"],
        "inputs": {
            name: {
                "path": str(args.de_dir / name),
                "sha256": sha256_of_file(args.de_dir / name),
            }
            for name in [
                "common_set_comparison.tsv",
                "summary.json",
                *(f"{q.key}_{d}.tsv" for q in QUANTITIES for d in (DESIGN, "batch")),
            ]
        },
        "de_script": {"path": summary["script"], "sha256": summary["script_sha256"]},
        "params": {
            "design": DESIGN,
            "design_check": "common-set log2fc + residual_sd == <q>_paired.tsv "
            "exactly and residual_sd != <q>_batch.tsv",
            "mean_log2_abundance_source": "<q>_paired.tsv mean_log2_abundance "
            "(joined on feature / first_member_id)",
            "ci_half_width": "(ci_high - ci_low) / 2; 95% t on moderated SE",
            "residual_sd": "per-feature unmoderated OLS residual SD, residual df 3",
            "n_bins": args.n_bins,
            "binning": "equal-count bins of mean log2 abundance; median x, median y",
            "histogram_bin_width": 0.025,
            "registry": str(args.registry),
            "registry_category": "quantity",
            "dpi": args.dpi,
        },
        "color_map": colors,
        "key_numbers": {
            "n_features": data.n_features,
            "spearman_log2fc": {k: round(v, 4) for k, v in rhos.items()},
            "median": {
                k: {q: round(v, 4) for q, v in m.items()} for k, m in medians.items()
            },
            "n_residual_sd_exactly_zero": zero_sd,
            "abundance_trend": trend,
            "pairs": [list(p) for p in PAIRS],
        },
        "figures": {
            STEMS["scatter"]: _artifact_paths(art_scatter),
            STEMS["distributions"]: _artifact_paths(art_dist),
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
