"""Stage-4 figures for findings 0006 (limma-trend sensitivity) and 0007 (ATPK case).

Loads the **saved** outputs of ``scripts/scratch/de_trend_sensitivity.py`` and
``scripts/scratch/atpk_case_study.py`` under ``results/de/raloxifene-vs-control/trend/``
(never re-fits a model) and renders:

into ``figures/analysis/differential-abundance/raloxifene-vs-control/`` (0006)

1. ``0006-prior-sd-trend-protein-paired`` — LFQ protein, paired design: per-protein
   residual SD vs mean log2 abundance with the constant no-trend prior SD and the
   limma-trend prior-SD curve.
2. ``0006-q-trend-vs-notrend-{protein,peptide}-paired`` — ``-log10 q`` limma-trend vs
   no-trend, colored by mean log2 abundance, trend hits (q < 0.05) labelled.
3. ``0006-volcano-{protein,peptide}-trend-paired`` — volcanoes of the limma-trend fit
   (project volcano module), hits labelled.
4. ``0006-relabel-diagnostic-trend`` — hit count and pi0 of the 8 within-pair
   labellings per quantity (observed highlighted).

into ``figures/analysis/quant-comparison/atpk-case/`` (0007)

5. ``0007-atpk-per-pair-by-quantity`` — ATPK per pair, control -> raloxifene-d0, for
   LFQ protein (normalized log2), NSAF (log2) and PSM count (log2 axis).
6. ``0007-atpk-peptides`` — per ATPK peptide: LFQ peptide intensity and Limelight PSM
   count per pair.

Every input is cross-checked before plotting (result objects vs TSVs; no-trend columns
vs the base DE tables; TSV-derived hits / pi0 vs summary.json; per-pair differences and
PSM totals vs atpk_case.json) — fail loud on any mismatch. Provenance (script + module
sha256, data_version, input sha256s, params, colors, drawn numbers, artifact paths) is
written to ``results/de/raloxifene-vs-control/trend/figure_provenance.json``.

Deterministic. Run from the project root:
    ./.venv/bin/python scripts/scratch/fig_trend_atpk.py
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis.differential_abundance import (  # noqa: E402
    DifferentialAbundanceTrendResult,
)
from analysis.result_io import load_result  # noqa: E402
from analysis_figures import pvalue_hist as ph  # noqa: E402
from analysis_figures import trend_atpk as ta  # noqa: E402
from analysis_figures import volcano as vo  # noqa: E402
from common.figures.colors import assign_colors  # noqa: E402
from common.figures.figure_io import (  # noqa: E402
    FigureArtifacts,
    publication_style,
    save_figure,
)
from common.hashing import sha256_of_file  # noqa: E402
from fig_de_raloxifene_vs_control import label_map  # noqa: E402

__script_meta__: dict[str, object] = {
    "kind": "script",
    "uses": [
        "analysis_figures.trend_atpk",
        "analysis_figures.volcano",
        "analysis_figures.pvalue_hist",
        "analysis.result_io",
        "analysis.differential_abundance",
        "common.figures.colors",
        "common.figures.figure_io",
        "common.hashing",
        "fig_de_raloxifene_vs_control.label_map",
    ],
    "description": (
        "Stage-4 figures for findings 0006 (limma-trend sensitivity of the "
        "raloxifene-d0 vs control DE) and 0007 (ATPK LFQ-vs-spectral case study). "
        "Loads cached results; no re-fitting."
    ),
}

LOG = logging.getLogger("fig_trend_atpk")

FDR = 0.05
FDR_LINES = (0.05, 0.10)
PI0_LAMBDA = 0.5
DESIGN = "paired"
CONTRAST_TERM = "condition[raloxifene-d0 vs control]"
EFFECT_LABEL = "log2 fold change (raloxifene-d0 / control)"
DIRECTION_LABELS = ("higher in raloxifene-d0", "higher in control")
DESIGN_TEXT = "paired design (condition + pair)"
PROCESSING = "LFQ median-normalized log2"
ABUNDANCE_LABEL = "mean log2 abundance (LFQ, median-normalized)"
CONDITIONS = ("control", "raloxifene-d0")
PAIRS = ("P905_906", "P907_908", "P909_910", "P941_942")
ATPK_GROUP = "psvid_86283_sp|P56134|ATPK_HUMAN"
_TOL = 1e-9
VOLCANO_SIZE = (12.0, 6.5)

RELABEL_QUANTITIES: tuple[tuple[str, str], ...] = (
    ("protein", "Protein\n(LFQ)"),
    ("peptide", "Peptide\n(LFQ)"),
    ("nsaf", "NSAF"),
    ("psm_log2", "PSM count\n(log2)"),
)
NOUN = {
    "protein": ("Protein", "proteins", "protein"),
    "peptide": ("Peptide", "peptides", "peptide"),
}

# ATPK peptide columns: (display, LFQ feature id, Limelight dump sequence).
ATPK_PEPTIDES: tuple[tuple[str, str, str], ...] = (
    ("DFSPSGIFGAFQR", "DFSPSGIFGAFQR[+0.0]", "DFSPSGIFGAFQR"),
    ("LGELPSWILM(ox)R", "LGELPSWILMR[+15.9949]", "LGELPSWILM[15.99]R"),
    ("LGELPSWILMR", "LGELPSWILMR[+0.0]", "LGELPSWILMR"),
)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _hits_word(n: int) -> str:
    return "hit" if n == 1 else "hits"


def _read_tsv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    if frame.empty:
        raise ValueError(f"{path}: empty table.")
    return frame


def _git(args: list[str], cwd: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", *args], cwd=cwd, capture_output=True, text=True, check=True
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.stdout.strip()


def _artifacts(a: FigureArtifacts) -> dict[str, str | None]:
    return {
        "svg": str(a.svg),
        "png": str(a.png),
        "legend_svg": str(a.legend_svg) if a.legend_svg else None,
        "legend_png": str(a.legend_png) if a.legend_png else None,
    }


def _check_result_vs_tsv(
    name: str, result: DifferentialAbundanceTrendResult, tsv: pd.DataFrame
) -> None:
    """Contrast rows of the cached trend result must equal the exported TSV."""
    if tuple(result.contrast_terms) != (CONTRAST_TERM,):
        raise ValueError(f"{name}: contrast terms {result.contrast_terms}.")
    if result.trend is not True:
        raise ValueError(f"{name}: cached result is not a limma-trend fit.")
    rows = result.table[result.table["term"] == CONTRAST_TERM].set_index("feature")
    t = tsv.set_index("feature")
    if set(rows.index) != set(t.index) or len(rows) != len(t):
        raise ValueError(f"{name}: result features differ from the TSV features.")
    t = t.loc[rows.index]
    for col_r, col_t in (
        ("effect", "log2fc"),
        ("p", "p"),
        ("q", "q"),
        ("mean_abundance", "mean_log2_abundance"),
        ("prior_variance", "prior_variance"),
        ("sigma", "residual_sd"),
    ):
        diff = np.abs(rows[col_r].to_numpy(float) - t[col_t].to_numpy(float))
        if not np.all(diff <= _TOL):
            raise ValueError(f"{name}: {col_r} differs from TSV by {diff.max():.3g}.")


def _check_notrend_vs_base(name: str, tsv: pd.DataFrame, base: pd.DataFrame) -> None:
    """The trend table's no-trend columns must equal the primary (no-trend) DE table."""
    merged = tsv.merge(base, on="feature", how="outer", suffixes=("", "_base"))
    if len(merged) != len(tsv) or len(merged) != len(base):
        raise ValueError(f"{name}: feature sets differ from the base DE table.")
    for col in ("p", "q"):
        diff = np.abs(
            merged[f"notrend_{col}"].to_numpy(float)
            - merged[f"{col}_base"].to_numpy(float)
        )
        if not np.all(diff <= _TOL):
            raise ValueError(f"{name}: notrend_{col} differs from base {col}.")


def _pairs_block(pairs_tsv: pd.DataFrame, quantity: str, feature: str) -> pd.DataFrame:
    block = pairs_tsv[
        (pairs_tsv["quantity"] == quantity) & (pairs_tsv["feature"] == feature)
    ].set_index("pair")
    if sorted(block.index) != sorted(PAIRS) or block.index.duplicated().any():
        raise ValueError(f"{quantity}/{feature}: pairs {list(block.index)}.")
    return block.loc[list(PAIRS)]


def _n_up_down(diff: np.ndarray) -> tuple[int, int, int]:
    d = diff[np.isfinite(diff)]
    return int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum())


def _direction_text(diff: np.ndarray) -> str:
    up, down, tie = _n_up_down(diff)
    n = up + down + tie
    if down == n:
        return f"{n}/{n} pairs down"
    if up == n:
        return f"{n}/{n} pairs up"
    parts = [f"{up}/{n} up"]
    if down:
        parts.append(f"{down} down")
    if tie:
        parts.append(f"{tie} tied")
    return ", ".join(parts)


# --------------------------------------------------------------------------- #
# Run
# --------------------------------------------------------------------------- #
def run(
    trend_dir: Path,
    base_de_dir: Path,
    da_output_dir: Path,
    atpk_output_dir: Path,
    registry_path: Path,
    qc_manifest_path: Path,
    provenance_path: Path,
) -> dict[str, Any]:
    """Render both families; return (and write) the provenance record."""
    summary_path = trend_dir / "summary.json"
    case_path = trend_dir / "atpk_case.json"
    pairs_path = trend_dir / "atpk_case_pairs.tsv"
    samples_path = trend_dir / "atpk_case_samples.tsv"
    summary = json.loads(summary_path.read_text())
    case = json.loads(case_path.read_text())
    qc_manifest = json.loads(qc_manifest_path.read_text())
    data_version = summary["data_version"]
    for label, dv in (
        ("atpk_case.json", case["data_version"]),
        ("qc manifest", qc_manifest.get("data_version")),
    ):
        if dv != data_version:
            raise ValueError(f"data_version mismatch: {label} {dv} vs {data_version}.")
    if summary["primary_design"] != DESIGN:
        raise ValueError(f"primary design is {summary['primary_design']!r}.")
    if case["protein_group"] != ATPK_GROUP:
        raise ValueError(f"atpk_case protein group {case['protein_group']!r}.")
    LOG.info("data_version %s", data_version)

    project_root = Path.cwd()
    script = Path(__file__).resolve()
    modules = {
        "trend_atpk": Path(ta.__file__).resolve(),
        "volcano": Path(vo.__file__).resolve(),
        "pvalue_hist": Path(ph.__file__).resolve(),
        "fig_de_raloxifene_vs_control(label_map)": (
            _SCRATCH / "fig_de_raloxifene_vs_control.py"
        ),
    }
    record: dict[str, Any] = {
        "script": {
            "path": str(script.relative_to(project_root)),
            "sha256": sha256_of_file(script),
        },
        "modules": {
            k: {"path": str(p.relative_to(project_root)), "sha256": sha256_of_file(p)}
            for k, p in modules.items()
        },
        "git_head": _git(["rev-parse", "HEAD"], project_root),
        "data_version": data_version,
        "inputs": {
            p.name: {"path": str(p), "sha256": sha256_of_file(p)}
            for p in (summary_path, case_path, pairs_path, samples_path)
        },
        "color_registry": str(registry_path),
        "params": {
            "fdr": FDR,
            "fdr_guide_lines": list(FDR_LINES),
            "hit_rule": "q < fdr (strict)",
            "pi0_lambda": PI0_LAMBDA,
            "design": DESIGN,
            "contrast_term": CONTRAST_TERM,
            "volcano_labels": (
                "all q < 0.05 hits; deterministic side-column layout "
                "(trend_atpk.label_hits_in_columns), boxes above the q line; "
                "figure resized to VOLCANO_SIZE"
            ),
            "paired_panels_equal_log2_span": True,
            "paired_panels_x_dodge": 0.05,
        },
        "figures": {},
    }

    # ---------------- load + cross-check trend tables / results ----------------
    tsvs: dict[str, pd.DataFrame] = {}
    results: dict[str, DifferentialAbundanceTrendResult] = {}
    for qty, _ in RELABEL_QUANTITIES:
        name = f"{qty}_{DESIGN}"
        tsv_path = trend_dir / f"{name}.tsv"
        base_path = base_de_dir / f"{name}.tsv"
        rdir = trend_dir / "results" / name
        tsv = _read_tsv(tsv_path)
        _check_notrend_vs_base(name, tsv, _read_tsv(base_path))
        res = load_result(rdir, DifferentialAbundanceTrendResult)
        _check_result_vs_tsv(name, res, tsv)
        tsvs[qty] = tsv
        results[qty] = res
        record["inputs"][f"{name}.tsv"] = {
            "path": str(tsv_path),
            "sha256": sha256_of_file(tsv_path),
        }
        record["inputs"][f"base/{name}.tsv"] = {
            "path": str(base_path),
            "sha256": sha256_of_file(base_path),
        }
        record["inputs"][f"result/{name}"] = {
            "result_dir": str(rdir),
            "files_sha256": {
                p.name: sha256_of_file(p) for p in sorted(rdir.iterdir()) if p.is_file()
            },
        }
        dsum = summary["quantities"][qty]["designs"][DESIGN]
        n_trend = int((tsv["q"] < FDR).sum())
        if n_trend != int(dsum["trend"]["hits"]["q<0.05"]["total"]):
            raise ValueError(f"{name}: trend hits {n_trend} != summary.")
        if int(tsv["hit_q05_trend"].sum()) != n_trend:
            raise ValueError(f"{name}: hit flag column disagrees with q < {FDR}.")
        pi0 = ph.estimate_pi0(tsv["p"].to_numpy(float), PI0_LAMBDA)
        if abs(pi0 - float(dsum["trend"]["storey_pi0_lambda0.5"])) > 5e-5:
            raise ValueError(f"{name}: pi0 {pi0:.4f} != summary.")
        LOG.info(
            "%s: %d features, %d trend hits, pi0 %.4f", name, len(tsv), n_trend, pi0
        )

    # ---------------- colors (registry) ----------------
    prior_colors = assign_colors(
        "eb_prior", ["no-trend", "limma-trend"], registry_path=registry_path
    )
    labelling_colors = assign_colors(
        "pair_labelling",
        ["observed"],
        background_values=["within-pair relabelling"],
        registry_path=registry_path,
    )
    pair_colors = assign_colors(
        "candidate_pair", list(PAIRS), registry_path=registry_path
    )
    record["colors"] = {
        "eb_prior": prior_colors,
        "pair_labelling": labelling_colors,
        "candidate_pair": pair_colors,
        "significance": "registry category 'Significance' (volcano module)",
        "q_scatter": f"sequential {ta.SEQUENTIAL_CMAP} by mean log2 abundance",
    }

    def save(
        stem: str, out_dir: Path, fig: ta.RenderedFigure, meta: dict[str, Any]
    ) -> None:
        arts = save_figure(fig.figure, out_dir, stem, legend_fig=fig.legend_figure)
        record["figures"][stem] = {
            **meta,
            "stats": fig.stats,
            "artifacts": _artifacts(arts),
        }
        LOG.info("wrote %s", arts.png)

    with publication_style():
        # ---------------- 1. prior SD trend (protein) ----------------
        prot = tsvs["protein"]
        psum = summary["quantities"]["protein"]["designs"][DESIGN]
        atpk_pos = np.flatnonzero(prot["feature"].to_numpy() == ATPK_GROUP)
        if atpk_pos.size != 1:
            raise ValueError("ATPK not found exactly once in the protein table.")
        notrend_sd = float(psum["notrend"]["prior_sd"])
        fig1 = ta.plot_prior_sd_trend(
            prot["mean_log2_abundance"],
            prot["residual_sd"],
            prot["prior_sd"],
            notrend_prior_sd=notrend_sd,
            series_colors=prior_colors,
            highlight={"ATPK_HUMAN": int(atpk_pos[0])},
            residual_df=int(psum["residual_df"]),
            d0_trend=float(psum["trend"]["prior_df_d0"]),
            d0_notrend=float(psum["notrend"]["prior_df_d0"]),
            feature_noun="protein",
            abundance_label=ABUNDANCE_LABEL,
            title=(
                "Empirical-Bayes prior SD: no-trend vs limma-trend\n"
                f"Protein · {PROCESSING} · {DESIGN_TEXT} · n = {len(prot):,} proteins"
            ),
        )
        lo_hi = psum["trend"]["prior_sd_range"]
        if (
            abs(float(prot["prior_sd"].min()) - lo_hi[0]) > 5e-4
            or abs(float(prot["prior_sd"].max()) - lo_hi[1]) > 5e-4
        ):
            raise ValueError("trend prior SD range != summary.")
        save(
            "0006-prior-sd-trend-protein-paired",
            da_output_dir,
            fig1,
            {"kind": "residual SD vs abundance + priors", "quantity": "protein"},
        )

        # ---------------- 2 + 3. q-q scatter and trend volcano ----------------
        for qty in ("protein", "peptide"):
            tsv = tsvs[qty]
            noun, plural, singular = NOUN[qty]
            lmap = label_map(qty, tsv)
            lmap_full = label_map(qty, tsv, abbreviate=False)
            hits = tsv[tsv["q"] < FDR].sort_values(["q", "p"], kind="stable")
            hit_pos = [int(tsv.index.get_loc(i)) for i in hits.index]
            n_notrend = int((tsv["notrend_q"] < FDR).sum())
            qq = ta.plot_q_trend_vs_notrend(
                tsv["q"],
                tsv["notrend_q"],
                tsv["mean_log2_abundance"],
                label_index=hit_pos,
                labels=[lmap[str(f)] for f in hits["feature"]],
                fdr_lines=FDR_LINES,
                abundance_label=ABUNDANCE_LABEL,
                feature_noun=singular,
                title=(
                    f"{noun} · BH q: limma-trend vs no-trend prior\n"
                    f"{PROCESSING} · {DESIGN_TEXT}\n"
                    f"n = {len(tsv):,} {plural} · q < {FDR:g}: {len(hits)} trend / "
                    f"{n_notrend} no-trend"
                ),
            )
            save(
                f"0006-q-trend-vs-notrend-{qty}-{DESIGN}",
                da_output_dir,
                qq,
                {
                    "kind": "q trend vs q no-trend",
                    "quantity": qty,
                    "labelled": [
                        {
                            "feature": str(r.feature),
                            "label": lmap[str(r.feature)],
                            "label_full": lmap_full[str(r.feature)],
                            "log2fc": float(r.log2fc),
                            "q_trend": float(r.q),
                            "q_notrend": float(r.notrend_q),
                            "mean_log2_abundance": float(r.mean_log2_abundance),
                        }
                        for r in hits.itertuples()
                    ],
                },
            )

            res = results[qty]
            dsum = summary["quantities"][qty]["designs"][DESIGN]["trend"]["hits"][
                "q<0.05"
            ]
            n_hits = int(dsum["total"])
            vplot, features = vo.volcano_from_result(
                res,
                term=CONTRAST_TERM,
                fdr=FDR,
                effect_label=EFFECT_LABEL,
                annotate_top=0,
                direction_labels=DIRECTION_LABELS,
                title=(
                    f"{noun}: raloxifene-d0 vs control · {DESIGN_TEXT} · "
                    f"limma-trend\n{PROCESSING} · n = {len(tsv):,} {plural} · "
                    f"{n_hits} {_hits_word(n_hits)} at q < {FDR:g} "
                    f"({dsum['up']} up / {dsum['down']} down)"
                ),
                legend_title=f"Significance (BH, n = {len(tsv):,})",
                registry_path=registry_path,
            )
            if (vplot.counts.up, vplot.counts.down) != (
                int(dsum["up"]),
                int(dsum["down"]),
            ):
                raise ValueError(f"{qty}: volcano up/down != summary.")
            rows = res.table[res.table["term"] == CONTRAST_TERM].reset_index(drop=True)
            if [str(f) for f in rows["feature"]] != features:
                raise ValueError(f"{qty}: volcano feature order changed.")
            hit_rows = rows[rows["q"] < FDR].sort_values(["q", "p"], kind="stable")
            hit_idx = [int(i) for i in hit_rows.index]
            # Wider canvas than the module's 8 x 6.5 in so the side label columns fit
            # without stretching the symmetric x-axis far past the data.
            vplot.figure.set_size_inches(VOLCANO_SIZE)
            cols = ta.label_hits_in_columns(
                vplot.figure.axes[0],
                rows["effect"].to_numpy(float),
                -np.log10(rows["q"].to_numpy(float)),
                hit_idx,
                [lmap[features[i]] for i in hit_idx],
                y_floor=-math.log10(FDR),
                fontsize=8.0,
            )
            # Rebuild the legend so it carries the labelled-ring entry (all hits).
            plt.close(vplot.legend_figure)
            legend_fig = vo._legend_figure(
                vplot.color_map,
                vplot.counts,
                DIRECTION_LABELS,
                f"Significance (BH, n = {len(tsv):,})",
                fdr=FDR,
                n_labelled=len(hit_idx),
                n_labelled_hits=len(hit_idx),
            )
            vol = ta.RenderedFigure(
                vplot.figure,
                legend_fig,
                {
                    "counts": {
                        "up": vplot.counts.up,
                        "down": vplot.counts.down,
                        "ns": vplot.counts.ns,
                    },
                    "min_q": float(rows["q"].min()),
                    "label_groups": cols.texts,
                    "leader_near_misses_lt_12px": cols.n_leader_near_misses,
                    "lowest_label_box_bottom_neg_log10_q": cols.lowest_box_bottom,
                    "q_line_neg_log10": -math.log10(FDR),
                    "x_limits": list(vplot.figure.axes[0].get_xlim()),
                    "effect_range": [
                        float(rows["effect"].min()),
                        float(rows["effect"].max()),
                    ],
                },
            )
            save(
                f"0006-volcano-{qty}-trend-{DESIGN}",
                da_output_dir,
                vol,
                {
                    "kind": "volcano (limma-trend)",
                    "quantity": qty,
                    "color_map": vplot.color_map,
                    "labelled": [
                        {
                            "feature": features[i],
                            "label": text,
                            "log2fc": float(rows.loc[i, "effect"]),
                            "q": float(rows.loc[i, "q"]),
                        }
                        for i, text in ((i, lmap[features[i]]) for i in hit_idx)
                    ],
                },
            )

        # ---------------- 4. relabel diagnostic ----------------
        relabel: list[ta.RelabelQuantity] = []
        relabel_meta: dict[str, Any] = {}
        min_p = set()
        for qty, display in RELABEL_QUANTITIES:
            diag = summary["quantities"][qty]["pair_relabel_diagnostic"]
            if diag["design"] != DESIGN:
                raise ValueError(f"{qty}: relabel design {diag['design']!r}.")
            labs = diag["labellings"]
            if len(labs) != int(diag["n_labellings"]):
                raise ValueError(f"{qty}: labelling count mismatch.")
            obs = [i for i, lab in enumerate(labs) if lab["observed"]]
            if len(obs) != 1 or labs[obs[0]]["flipped_pairs"]:
                raise ValueError(
                    f"{qty}: expected exactly one unflipped observed labelling."
                )
            hits = [int(lab["trend"]["hits_q<0.05"]) for lab in labs]
            pi0s = [float(lab["trend"]["pi0"]) for lab in labs]
            o = obs[0]
            n_ge = sum(h >= hits[o] for h in hits)
            rank = 1 + sum(p < pi0s[o] for i, p in enumerate(pi0s) if i != o)
            standing = diag["observed_standing"]["trend"]
            if n_ge != int(standing["n_labellings_with_hits_q<0.05_ge_observed"]):
                raise ValueError(f"{qty}: hits >= observed {n_ge} != summary.")
            if rank != int(standing["observed_pi0_rank_lowest_first"]):
                raise ValueError(f"{qty}: pi0 rank {rank} != summary.")
            tsv = tsvs[qty]
            if hits[o] != int((tsv["q"] < FDR).sum()):
                raise ValueError(f"{qty}: observed labelling hits != TSV.")
            if (
                abs(pi0s[o] - ph.estimate_pi0(tsv["p"].to_numpy(float), PI0_LAMBDA))
                > 5e-5
            ):
                raise ValueError(f"{qty}: observed labelling pi0 != TSV.")
            min_p.add(float(diag["smallest_attainable_restricted_permutation_p"]))
            relabel.append(ta.RelabelQuantity(display, hits, pi0s, o, n_ge, rank))
            relabel_meta[qty] = [
                {"flipped_pairs": lab["flipped_pairs"], "hits": h, "pi0": p}
                for lab, h, p in zip(labs, hits, pi0s, strict=True)
            ]
        if min_p != {1 / len(relabel[0].hits)}:
            raise ValueError(f"min attainable p {min_p} != 1/{len(relabel[0].hits)}.")
        n_lab = len(relabel[0].hits)
        fig4 = ta.plot_relabel_diagnostic(
            relabel,
            observed_color=labelling_colors["observed"],
            other_color=labelling_colors["within-pair relabelling"],
            fdr=FDR,
            title=(
                f"Within-pair relabelling · limma-trend · {DESIGN_TEXT}\n"
                f"{n_lab} labellings (observed + {n_lab - 1} pair-flip sets) · "
                f"min attainable permutation p = 1/{n_lab}"
            ),
        )
        save(
            "0006-relabel-diagnostic-trend",
            da_output_dir,
            fig4,
            {"kind": "pair-relabel diagnostic", "labellings": relabel_meta},
        )

        # ---------------- 5. ATPK per pair by quantity ----------------
        pairs_tsv = _read_tsv(pairs_path)
        samples_tsv = _read_tsv(samples_path)
        batch_of = (
            pairs_tsv.drop_duplicates("pair").set_index("pair")["batch"].to_dict()
        )
        pair_labels = {p: f"{p} · {batch_of[p]}" for p in PAIRS}
        de = case["de_results"]
        q_specs = (
            (
                "lfq_protein_norm_log2",
                "LFQ protein",
                "log2 intensity (median-norm.)",
                False,
                de["lfq_protein"]["paired"],
            ),
            ("nsaf_log2", "NSAF", "log2 NSAF (as-is)", False, de["nsaf_paired"]),
            (
                "psm_count",
                "PSM count",
                "PSM count (log2 axis)",
                True,
                de["psm_log2_paired"],
            ),
        )
        panels5: list[ta.PairedPanel] = []
        per_pair5: dict[str, Any] = {}
        for quantity, title, ylabel, count_axis, de_row in q_specs:
            block = _pairs_block(pairs_tsv, quantity, ATPK_GROUP)
            c = block["control_value"].to_numpy(float)
            t = block["raloxifene_value"].to_numpy(float)
            diff = np.log2(t / c) if count_axis else t - c
            ref = case["protein_within_pair"][quantity]["per_pair"]
            for p, d in zip(PAIRS, diff, strict=True):
                if abs(float(d) - float(ref[p])) > 5e-4:
                    raise ValueError(f"{quantity}/{p}: diff {d} != atpk_case {ref[p]}.")
            if abs(float(diff.mean()) - float(de_row["log2fc"])) > 5e-4:
                raise ValueError(f"{quantity}: mean diff != DE log2FC.")
            note = (
                f"log2FC {ta.fmt_signed(float(de_row['log2fc']))} · "
                f"q = {ta.fmt_q(float(de_row['q']))} · {_direction_text(diff)}"
            )
            panels5.append(
                ta.PairedPanel(title, ylabel, list(c), list(t), count_axis, note)
            )
            per_pair5[quantity] = {
                "control": dict(zip(PAIRS, c.tolist(), strict=True)),
                "raloxifene": dict(zip(PAIRS, t.tolist(), strict=True)),
                "within_pair_log2_difference": dict(
                    zip(PAIRS, diff.tolist(), strict=True)
                ),
                "de_paired_limma_trend": de_row,
            }
        lims5 = ta.equal_span_limits([[p] for p in panels5])
        fig5 = ta.plot_paired_panels(
            [panels5],
            [lims5],
            pairs=PAIRS,
            pair_colors=pair_colors,
            pair_legend_labels=pair_labels,
            condition_labels=CONDITIONS,
            title=(
                "ATPK (P56134) within-pair change by quantity, "
                "control → raloxifene-d0\n"
                f"{DESIGN_TEXT} · limma-trend q · n = 4 pairs (8 runs) · "
                "equal log2 span per panel"
            ),
            panel_size=(3.7, 4.1),
        )
        save(
            "0007-atpk-per-pair-by-quantity",
            atpk_output_dir,
            fig5,
            {"kind": "per-pair segments by quantity", "per_pair": per_pair5},
        )

        # ---------------- 6. ATPK peptides: LFQ vs PSM count ----------------
        lfq_row: list[ta.PairedPanel] = []
        psm_row: list[ta.PairedPanel] = []
        pep_meta: dict[str, Any] = {}
        dump = {d["dump_sequence"]: d for d in case["peptide_dump_psms"]}
        lfq_info = {d["peptide"]: d for d in case["lfq_peptides"]}
        for col, (display, lfq_id, dump_seq) in enumerate(ATPK_PEPTIDES):
            b = _pairs_block(pairs_tsv, "lfq_peptide_norm_log2", lfq_id)
            c = b["control_value"].to_numpy(float)
            t = b["raloxifene_value"].to_numpy(float)
            det = samples_tsv[
                (samples_tsv["quantity"] == "lfq_peptide_norm_log2")
                & (samples_tsv["feature"] == lfq_id)
            ]
            if len(det) != 8:
                raise ValueError(f"{lfq_id}: {len(det)} sample rows, expected 8.")
            n_quant = int(det["value"].notna().sum())
            if n_quant != int(np.isfinite(c).sum() + np.isfinite(t).sum()):
                raise ValueError(f"{lfq_id}: sample/pair quantified counts disagree.")
            det_counts = det["detection_type"].value_counts().to_dict()
            if det_counts != lfq_info[lfq_id]["detection_type_counts"]:
                raise ValueError(f"{lfq_id}: detection types != atpk_case.")
            complete = np.isfinite(c) & np.isfinite(t)
            if complete.all():
                diff = t - c
                ref = lfq_info[lfq_id]["within_pair_norm_log2_difference"]["per_pair"]
                for p, d in zip(PAIRS, diff, strict=True):
                    if abs(float(d) - float(ref[p])) > 5e-4:
                        raise ValueError(f"{lfq_id}/{p}: diff != atpk_case.")
                lfq_note: str | None = (
                    f"mean Δ {ta.fmt_signed(float(diff.mean()))} · "
                    f"{_direction_text(diff)}"
                )
                missing_note: str | None = None
            else:
                if complete.any():
                    raise ValueError(f"{lfq_id}: partially complete pairs not handled.")
                n_amb = int(det_counts.get("MSMSAmbiguousPeakfinding", 0))
                lfq_note = "no complete pair"
                missing_note = (
                    f"LFQ in {n_quant}/8 runs\n{n_amb}/8 MSMSAmbiguousPeakfinding"
                )
            lfq_row.append(
                ta.PairedPanel(
                    display,
                    "LFQ log2 intensity (median-norm.)" if col == 0 else "",
                    list(c),
                    list(t),
                    False,
                    lfq_note,
                    allow_missing=not complete.all(),
                    missing_note=missing_note,
                )
            )
            bp = _pairs_block(pairs_tsv, "psm_count_peptide_dump", dump_seq)
            pc = bp["control_value"].to_numpy(float)
            pt = bp["raloxifene_value"].to_numpy(float)
            dref = dump[dump_seq]
            if (int(pc.sum()), int(pt.sum())) != (
                int(dref["psms_total_control"]),
                int(dref["psms_total_raloxifene"]),
            ):
                raise ValueError(f"{dump_seq}: PSM totals != atpk_case.")
            pdiff = np.log2(pt / pc)
            psm_row.append(
                ta.PairedPanel(
                    display,
                    "PSM count (log2 axis)" if col == 0 else "",
                    list(pc),
                    list(pt),
                    True,
                    (
                        f"Σ PSMs {int(pc.sum())} → {int(pt.sum())} · "
                        f"{_direction_text(pdiff)}"
                    ),
                )
            )
            pep_meta[display] = {
                "lfq_feature": lfq_id,
                "dump_sequence": dump_seq,
                "lfq_control": dict(
                    zip(PAIRS, [None if math.isnan(v) else v for v in c], strict=True)
                ),
                "lfq_raloxifene": dict(
                    zip(PAIRS, [None if math.isnan(v) else v for v in t], strict=True)
                ),
                "lfq_detection_types": det_counts,
                "psm_control": dict(zip(PAIRS, pc.astype(int).tolist(), strict=True)),
                "psm_raloxifene": dict(
                    zip(PAIRS, pt.astype(int).tolist(), strict=True)
                ),
                "psm_totals": [int(pc.sum()), int(pt.sum())],
            }
        lims6 = ta.equal_span_limits([lfq_row, psm_row])
        fig6 = ta.plot_paired_panels(
            [lfq_row, psm_row],
            [[lims6[0]] * len(lfq_row), [lims6[1]] * len(psm_row)],
            pairs=PAIRS,
            pair_colors=pair_colors,
            pair_legend_labels=pair_labels,
            condition_labels=CONDITIONS,
            title=(
                "ATPK peptides: LFQ MS1 intensity (top) vs PSM count (bottom), "
                "control → raloxifene-d0\n"
                "FlashLFQ peptide, median-normalized log2 · Limelight PSMs · "
                "n = 4 pairs (8 runs) · equal log2 span per row"
            ),
            panel_size=(3.6, 3.7),
        )
        save(
            "0007-atpk-peptides",
            atpk_output_dir,
            fig6,
            {"kind": "per-peptide LFQ vs PSM per pair", "peptides": pep_meta},
        )

    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    LOG.info("provenance -> %s", provenance_path)
    return record


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    parser.add_argument(
        "--trend-dir", type=Path, default=Path("results/de/raloxifene-vs-control/trend")
    )
    parser.add_argument(
        "--base-de-dir", type=Path, default=Path("results/de/raloxifene-vs-control")
    )
    parser.add_argument(
        "--da-output-dir",
        type=Path,
        default=Path("figures/analysis/differential-abundance/raloxifene-vs-control"),
    )
    parser.add_argument(
        "--atpk-output-dir",
        type=Path,
        default=Path("figures/analysis/quant-comparison/atpk-case"),
    )
    parser.add_argument(
        "--registry-path", type=Path, default=Path("state/color_registry.json")
    )
    parser.add_argument(
        "--qc-manifest-path", type=Path, default=Path("results/qc_states/manifest.json")
    )
    parser.add_argument(
        "--provenance-path",
        type=Path,
        default=Path("results/de/raloxifene-vs-control/trend/figure_provenance.json"),
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run(
        args.trend_dir,
        args.base_de_dir,
        args.da_output_dir,
        args.atpk_output_dir,
        args.registry_path,
        args.qc_manifest_path,
        args.provenance_path,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
