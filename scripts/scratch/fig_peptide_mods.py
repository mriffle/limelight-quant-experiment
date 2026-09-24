"""Stage-4 peptide-modification figures for findings 0008 / 0009.

Finding 0008 (raloxifene Cys adducts), to ``--adduct-output-dir``:
  * ``0008-adduct-detection-heatmap`` -- 60 adduct peptides x 8 runs (pair order),
    FlashLFQ detection type per cell, Limelight PSM count in the cell;
  * ``0008-adduct-discordance-by-residue`` -- treated-only share of discordant pair
    cells (Clopper-Pearson 95 % CI) for Cys / Tyr-Trp adducts and unmodified peptides,
    per detection definition (quantified / MS/MS only / PSM);
  * ``0008-cyp3a4-adduct-peptide`` -- the CYP3A4 GFCMFDMECHK adduct: per-run log2
    intensity and PSMs.
Finding 0009 (class-level shifts), to ``--class-output-dir``:
  * ``0009-class-log2fc-by-pair`` -- per-pair class shifts vs unmodified (CAM-only,
    oxidized) and oxidized modified-minus-reference, mean +/- t CI across 4 pairs;
  * ``0009-class-log2fc-distributions`` -- paired peptide log2FC densities by class;
  * ``0009-detection-by-run-position`` -- unmodified peptides quantified per run in
    acquisition order + per-pair treated-only share.

Inputs: the plot-ready tables of ``scripts/scratch/peptide_mod_analysis.py`` in
``--in-dir`` (default ``results/peptide-mods``); nothing is recomputed from raw data.
Every plotted number that the analysis also recorded is cross-checked against it (fail
loud). Provenance (script hashes, data_version, input hashes, params, key numbers,
artifact paths, colors) goes to ``--provenance``. Deterministic; no seed consumed.
Registry: adds project categories ``detection_type`` and ``mod_class`` (deterministic
Okabe-Ito assignment) if absent; ``condition`` / ``candidate_pair`` are read.

Run (from the project root):
    ./.venv/bin/python scripts/scratch/fig_peptide_mods.py
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

import numpy as np
import pandas as pd

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
_ROOT = _SCRATCH.parent.parent
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis_figures.peptide_mods_figs import (  # noqa: E402
    CONTROL,
    DETECTION_LABELS,
    TREATED,
    DiscordancePanel,
    DiscordancePoint,
    PairShare,
    SampleInfo,
    ShiftGroup,
    clopper_pearson,
    plot_adduct_detection_heatmap,
    plot_class_log2fc_distributions,
    plot_class_shift_by_pair,
    plot_detection_by_run,
    plot_discordance_by_residue,
    plot_single_adduct_peptide,
)
from common.figures.colors import assign_colors  # noqa: E402
from common.figures.figure_io import (  # noqa: E402
    FigureArtifacts,
    publication_style,
    save_figure,
)
from common.hashing import sha256_of_file  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "fig-peptide-mods",
    "kind": "figure",
    "provides": [],
    "uses": [
        "analysis_figures.peptide_mods_figs",
        "common.figures.colors",
        "common.figures.figure_io",
        "common.hashing",
    ],
    "seeded_from": None,
    "description": (
        "Figures for findings 0008 (raloxifene adducts: detection heatmap, "
        "discordance by residue, CYP3A4 peptide) and 0009 (class log2FC shifts by "
        "pair, class distributions, detection by run position); numbers cross-checked "
        "against the peptide-mods analysis outputs; provenance JSON."
    ),
}

LOG = logging.getLogger("fig_peptide_mods")

STEMS = {
    "heatmap": "0008-adduct-detection-heatmap",
    "discordance": "0008-adduct-discordance-by-residue",
    "cyp3a4": "0008-cyp3a4-adduct-peptide",
    "shift": "0009-class-log2fc-by-pair",
    "dist": "0009-class-log2fc-distributions",
    "runs": "0009-detection-by-run-position",
}
CYP3A4_FEATURE = "GFCMFDMECHK[+528.171864]"
DEFINITIONS = (
    ("quantified", "quantified (MS/MS or MBR)"),
    ("msms", "MS/MS only"),
    ("psm", "Limelight PSM > 0"),
)
PLOT_CLASSES = ("unmodified", "CAM-only", "oxidized")
FLASHLFQ_TYPES = ("MSMS", "MBR", "MSMSIdentifiedButNotQuantified")
GROUP_CYS = "Cys adducts"
GROUP_TYRTRP = "Tyr/Trp adducts"
GROUP_UNMOD = "unmodified peptides"
INPUTS = (
    "summary.json",
    "adduct_peptides_long.tsv",
    "adduct_peptides.tsv",
    "adduct_pair_summary.tsv",
    "adduct_proteins.tsv",
    "discordance_class_tests.tsv",
    "discordance_by_pair.tsv",
    "class_log2fc.tsv",
    "class_log2fc_tests.tsv",
    "modified_vs_reference.tsv",
    "modified_vs_reference_class_tests.tsv",
)
_ATOL = 1e-6
SEP = " \u00b7 "  # middle-dot separator
TIMES = "\u00d7"
MINUS = "\u2212"
_ATOL_4DP = 6e-5  # per-pair strings are rounded to 4 dp


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    p.add_argument("--in-dir", type=Path, default=Path("results/peptide-mods"))
    p.add_argument(
        "--adduct-output-dir",
        type=Path,
        default=Path("figures/analysis/peptide-mods/raloxifene-adducts"),
    )
    p.add_argument(
        "--class-output-dir",
        type=Path,
        default=Path("figures/analysis/peptide-mods/class-shifts"),
    )
    p.add_argument(
        "--provenance",
        type=Path,
        default=Path("results/peptide-mods/figure_provenance.json"),
    )
    p.add_argument("--registry", type=Path, default=Path("state/color_registry.json"))
    p.add_argument("--dist-xlim", type=float, default=1.25)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args(argv)


def _read(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required input missing: {path}")
    return pd.read_csv(path, sep="\t")


def _close(label: str, a: float, b: float, atol: float = _ATOL) -> None:
    if not (math.isfinite(a) and math.isfinite(b)) or abs(a - b) > atol:
        raise ValueError(f"cross-check {label}: {a!r} != {b!r} (atol {atol}).")


def _per_pair(text: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for item in text.split(";"):
        key, val = item.split(":")
        out[key] = float(val)
    return out


def _samples(long: pd.DataFrame) -> list[SampleInfo]:
    meta = long[
        [
            "sample_id",
            "condition",
            "candidate_pair",
            "batch",
            "run_position_within_batch",
        ]
    ].drop_duplicates()
    if not meta["sample_id"].is_unique or len(meta) != 8:
        raise ValueError("expected 8 samples with one metadata row each.")
    if set(meta["condition"]) != {CONTROL, TREATED}:
        raise ValueError(f"unexpected conditions {set(meta['condition'])}.")
    out = [
        SampleInfo(str(r[0]), str(r[1]), str(r[2]), str(r[3]), int(r[4]))
        for r in meta.itertuples(index=False)
    ]
    for pair in {s.pair for s in out}:
        conds = sorted(s.condition for s in out if s.pair == pair)
        if conds != [CONTROL, TREATED]:
            raise ValueError(f"pair {pair} is not one control + one treated.")
    return out


def _artifact_paths(a: FigureArtifacts) -> dict[str, str | None]:
    return {k: (str(v) if v is not None else None) for k, v in vars(a).items()}


def _git_head() -> str | None:
    try:
        res = subprocess.run(
            ["git", "-C", str(_ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return res.stdout.strip() or None


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args(argv)
    LOG.info("params %s", vars(args))
    d = args.in_dir
    summary: dict[str, Any] = json.loads((d / "summary.json").read_text())
    long = _read(d / "adduct_peptides_long.tsv")
    adducts = _read(d / "adduct_peptides.tsv")
    pair_summary = _read(d / "adduct_pair_summary.tsv")
    proteins = _read(d / "adduct_proteins.tsv")
    disc_tests = _read(d / "discordance_class_tests.tsv")
    disc_pair = _read(d / "discordance_by_pair.tsv")
    fc = _read(d / "class_log2fc.tsv")
    fc_tests = _read(d / "class_log2fc_tests.tsv").set_index("mod_class")
    mvr = _read(d / "modified_vs_reference.tsv")
    mvr_tests = _read(d / "modified_vs_reference_class_tests.tsv").set_index(
        "mod_class"
    )

    samples = _samples(long)
    pair_order = sorted({s.pair for s in samples})
    by_pair = [
        s
        for p in pair_order
        for s in sorted(
            (x for x in samples if x.pair == p), key=lambda x: x.condition != CONTROL
        )
    ]
    run_order = sorted(samples, key=lambda s: (s.batch, s.run_position))

    # ---------------- colors (registry) ----------------
    reg = args.registry
    cond_colors = assign_colors("condition", [CONTROL, TREATED], registry_path=reg)
    pair_colors = assign_colors("candidate_pair", pair_order, registry_path=reg)
    det_raw = assign_colors("detection_type", list(FLASHLFQ_TYPES), registry_path=reg)
    det_colors = {DETECTION_LABELS[k]: v for k, v in det_raw.items()}
    class_colors = assign_colors("mod_class", list(PLOT_CLASSES), registry_path=reg)
    group_colors = {
        GROUP_CYS: "#000000",
        GROUP_TYRTRP: "#000000",
        GROUP_UNMOD: class_colors["unmodified"],
    }
    pair_colors = {p: pair_colors[p] for p in pair_order}

    # ---------------- cross-checks: adduct grid ----------------
    n_add = int(summary["adducts"]["n_adduct_features"])
    if len(adducts) != n_add:
        raise ValueError("adduct_peptides.tsv row count != summary n_adduct_features.")
    quant = long["detection_type"].isin(["MSMS", "MBR"])
    per_feat = long.assign(q=quant).pivot_table(
        index="feature", columns="condition", values="q", aggfunc="sum"
    )
    chk = adducts.set_index("feature")
    if not (per_feat[TREATED] == chk.loc[per_feat.index, "n_treated_quantified"]).all():
        raise ValueError("per-feature treated quantified counts disagree.")
    if not (per_feat[CONTROL] == chk.loc[per_feat.index, "n_control_quantified"]).all():
        raise ValueError("per-feature control quantified counts disagree.")
    ps = pair_summary[pair_summary["definition"] == "quantified"].set_index(
        "candidate_pair"
    )
    per_sample_q = long.assign(q=quant).groupby("sample_id")["q"].sum()
    for pair, r in ps.iterrows():
        if per_sample_q[r["control_sample"]] != r["n_control_detected"]:
            raise ValueError(f"{pair}: control adduct count disagrees.")
        if per_sample_q[r["treated_sample"]] != r["n_treated_detected"]:
            raise ValueError(f"{pair}: treated adduct count disagrees.")
    cyp_entries = set(
        proteins.loc[proteins["is_cyp"].astype(bool), "entry"].astype(str)
    )
    if cyp_entries != set(summary["adducts"]["cyp_proteins"]):
        raise ValueError("CYP entry set disagrees with summary.json.")

    # ---------------- discordance panels ----------------
    tests = summary["adducts"]["tests"]
    panels: list[DiscordancePanel] = []
    disc_numbers: dict[str, Any] = {}
    for key, name in DEFINITIONS:
        pts: list[DiscordancePoint] = []
        for grp, src in ((GROUP_CYS, "Cys"), (GROUP_TYRTRP, "non-Cys (Tyr/Trp)")):
            r = tests[key]["by_residue"][src]
            k, n = int(r["k"]), int(r["n"])
            lo, hi = clopper_pearson(k, n)
            _close(f"{key} {src} ci_low", lo, float(r["ci_low"]))
            _close(f"{key} {src} ci_high", hi, float(r["ci_high"]))
            _close(f"{key} {src} prop", k / n, float(r["prop"]))
            pts.append(DiscordancePoint(grp, k, n, k / n, lo, hi))
        row = disc_tests[
            (disc_tests["definition"] == key)
            & (disc_tests["mod_class"] == "unmodified")
        ].iloc[0]
        k = int(row["treated_only_cells"])
        n = k + int(row["control_only_cells"])
        lo, hi = clopper_pearson(k, n)
        _close(f"{key} unmodified ci_low", lo, float(row["prop_ci_low"]))
        _close(f"{key} unmodified ci_high", hi, float(row["prop_ci_high"]))
        pts.append(DiscordancePoint(GROUP_UNMOD, k, n, k / n, lo, hi))
        adr = disc_tests[
            (disc_tests["definition"] == key)
            & (disc_tests["mod_class"] == "raloxifene-adduct")
        ].iloc[0]
        orv, olo, ohi = (
            float(adr["pair_level_or"]),
            float(adr["pair_level_ci_low"]),
            float(adr["pair_level_ci_high"]),
        )
        n_gt1 = int(adr["pair_level_n_pairs_or_gt1"])
        note = (
            "all adducts vs unmodified:\n"
            f"pair-level OR {orv:.0f} [{olo:.0f}, {ohi:.0f}], {n_gt1}/4 pairs > 1"
        )
        panels.append(DiscordancePanel(name, tuple(pts), note))
        disc_numbers[key] = {
            "points": {
                p.group: {
                    "k": p.k,
                    "n": p.n,
                    "prop": round(p.prop, 4),
                    "ci": [round(p.ci_low, 4), round(p.ci_high, 4)],
                }
                for p in pts
            },
            "adduct_vs_unmodified_pair_level_or": [
                round(orv, 3),
                round(olo, 3),
                round(ohi, 3),
            ],
            "pairs_or_gt1": n_gt1,
        }

    # ---------------- class shifts ----------------
    shift_groups: list[ShiftGroup] = []
    shift_numbers: dict[str, Any] = {}
    for cls, label in (
        ("CAM-only", f"CAM-only\n{MINUS} unmodified"),
        ("oxidized", f"oxidized\n{MINUS} unmodified"),
    ):
        r = fc_tests.loc[cls]
        pp = _per_pair(str(r["pair_shift_per_pair"]))
        mean = float(r["pair_shift_mean"])
        _close(
            f"{cls} mean of per-pair",
            float(np.mean(list(pp.values()))),
            mean,
            _ATOL_4DP,
        )
        n_pos = sum(v > 0 for v in pp.values())
        if n_pos != int(r["pair_shift_n_pairs_positive"]):
            raise ValueError(f"{cls}: n pairs positive disagrees.")
        sign = f"{n_pos}/4 pairs > 0" if n_pos >= 2 else f"{4 - n_pos}/4 pairs < 0"
        note = f"{sign}\nBH q {float(r['pair_shift_q_BH']):.3f}"
        shift_groups.append(
            ShiftGroup(
                label,
                pp,
                mean,
                float(r["pair_shift_ci_low"]),
                float(r["pair_shift_ci_high"]),
                note,
            )
        )
        shift_numbers[cls] = {
            "per_pair": pp,
            "mean": round(mean, 4),
            "ci": [
                round(float(r["pair_shift_ci_low"]), 4),
                round(float(r["pair_shift_ci_high"]), 4),
            ],
            "q_BH": round(float(r["pair_shift_q_BH"]), 4),
        }
    ox = mvr[mvr["mod_class"] == "oxidized"]
    rt = mvr_tests.loc["oxidized"]
    if len(ox) != int(rt["n_forms"]):
        raise ValueError("oxidized modified-vs-reference form count disagrees.")
    pp_mvr = {p: float(ox[f"diff_{p}"].mean()) for p in pair_order}
    for p, v in _per_pair(str(rt["pair_level_per_pair"])).items():
        _close(f"mvr per-pair {p}", pp_mvr[p], v, _ATOL_4DP)
    mean_mvr = float(rt["pair_level_mean"])
    _close("mvr mean", float(np.mean(list(pp_mvr.values()))), mean_mvr)
    n_pos = sum(v > 0 for v in pp_mvr.values())
    shift_groups.append(
        ShiftGroup(
            f"oxidized form\n{MINUS} unoxidized form",
            pp_mvr,
            mean_mvr,
            float(rt["pair_level_ci_low"]),
            float(rt["pair_level_ci_high"]),
            f"{n_pos}/4 pairs > 0\nBH q {float(rt['pair_level_q_BH']):.3f}",
        )
    )
    shift_numbers["oxidized_modified_minus_reference"] = {
        "per_pair": {k: round(v, 4) for k, v in pp_mvr.items()},
        "mean": round(mean_mvr, 4),
        "ci": [
            round(float(rt["pair_level_ci_low"]), 4),
            round(float(rt["pair_level_ci_high"]), 4),
        ],
        "n_forms": len(ox),
        "q_BH": round(float(rt["pair_level_q_BH"]), 4),
    }

    # ---------------- distributions ----------------
    if len(fc) != int(summary["intensity_by_class"]["n_complete_peptides"]):
        raise ValueError("class_log2fc row count != n_complete_peptides.")
    dist_vals = {
        c: fc.loc[fc["mod_class"] == c, "log2fc"].to_numpy(dtype=np.float64)
        for c in PLOT_CLASSES
    }
    for c, v in dist_vals.items():
        if v.size != int(fc_tests.loc[c, "n"]):
            raise ValueError(f"{c}: n disagrees with class_log2fc_tests.")
        _close(f"{c} median", float(np.median(v)), float(fc_tests.loc[c, "median"]))

    # ---------------- detection by run ----------------
    dq = disc_pair[
        (disc_pair["definition"] == "quantified")
        & (disc_pair["mod_class"] == "unmodified")
    ].set_index("candidate_pair")
    counts: dict[str, int] = {}
    shares: list[PairShare] = []
    for p in pair_order:
        r = dq.loc[p]
        ctrl = next(s for s in samples if s.pair == p and s.condition == CONTROL)
        trt = next(s for s in samples if s.pair == p and s.condition == TREATED)
        counts[ctrl.sample_id] = int(r["n_detected_control"])
        counts[trt.sample_id] = int(r["n_detected_treated"])
        k = int(r["treated_only_cells"])
        n = k + int(r["control_only_cells"])
        lo, hi = clopper_pearson(k, n)
        _close(f"{p} unmodified share", k / n, float(r["prop_treated_only"]))
        _close(f"{p} unmodified ci_low", lo, float(r["prop_ci_low"]))
        _close(f"{p} unmodified ci_high", hi, float(r["prop_ci_high"]))
        shares.append(PairShare(p, str(r["batch"]), k, n, k / n, lo, hi))

    n_cys = int(tests["quantified"]["by_residue"]["Cys"]["n_features"])
    unmod_rows = disc_tests[disc_tests["mod_class"] == "unmodified"]
    n_unmod = int(unmod_rows["n_peptides"].iloc[0])
    with publication_style():
        fig, leg, heat_stats = plot_adduct_detection_heatmap(
            long,
            adducts,
            by_pair,
            cyp_entries,
            det_colors,
            cond_colors,
            title=(
                f"Raloxifene-adduct peptides: detection per run "
                f"({n_add} peptides {TIMES} 8 runs)"
            ),
            subtitle=SEP.join(
                [
                    "FlashLFQ detection type (raw, per run)",
                    "runs grouped by pair",
                    "rows sorted by raloxifene runs quantified",
                ]
            ),
        )
        if heat_stats["n_cys"] != n_cys:
            raise ValueError("heatmap Cys row count disagrees with summary.")
        art_heat = save_figure(
            fig, args.adduct_output_dir, STEMS["heatmap"], legend_fig=leg, dpi=args.dpi
        )

        fig, leg = plot_discordance_by_residue(
            panels,
            group_colors,
            title=(
                "Treated-only share of discordant pair cells, "
                "adducts by residue vs unmodified"
            ),
            subtitle=SEP.join(
                [
                    f"per pair {TIMES} peptide cell detected in one arm only",
                    "4 pairs",
                    f"{n_cys} Cys + {n_add - n_cys} Tyr/Trp adduct peptides",
                    f"{n_unmod:,} unmodified",
                ]
            ),
        )
        art_disc = save_figure(
            fig,
            args.adduct_output_dir,
            STEMS["discordance"],
            legend_fig=leg,
            dpi=args.dpi,
        )

        fig, leg, cyp_stats = plot_single_adduct_peptide(
            long,
            CYP3A4_FEATURE,
            by_pair,
            cond_colors,
            title="CYP3A4 adduct peptide GFC[+471.15]MFDMEC[+57.02]HK",
            subtitle=SEP.join(
                [
                    "C3 raloxifene adduct, C9 carbamidomethyl",
                    "per run, 4 pairs",
                    "raw FlashLFQ intensity (not normalized)",
                ]
            ),
        )
        art_cyp = save_figure(
            fig, args.adduct_output_dir, STEMS["cyp3a4"], legend_fig=leg, dpi=args.dpi
        )

        fig, leg = plot_class_shift_by_pair(
            shift_groups,
            pair_colors,
            title="Class-level shifts in paired peptide log2FC, per pair",
            subtitle=SEP.join(
                [
                    f"raloxifene-d0 {MINUS} control",
                    "median-normalized log2",
                    f"{len(fc):,} complete peptides; "
                    f"{len(ox)} oxidized/unoxidized form pairs",
                ]
            ),
            ylabel="class shift in paired log2FC (log2 units)",
        )
        art_shift = save_figure(
            fig, args.class_output_dir, STEMS["shift"], legend_fig=leg, dpi=args.dpi
        )

        fig, leg, dist_stats = plot_class_log2fc_distributions(
            dist_vals,
            class_colors,
            title="Paired peptide log2FC by modification class",
            subtitle=SEP.join(
                [
                    f"raloxifene-d0 {MINUS} control, paired design (condition + pair)",
                    "median-normalized log2",
                    f"{len(fc):,} peptides quantified in all 8 runs",
                ]
            ),
            xlabel=f"paired log2FC, raloxifene-d0 {MINUS} control (log2 units)",
            xlim=(-args.dist_xlim, args.dist_xlim),
        )
        art_dist = save_figure(
            fig, args.class_output_dir, STEMS["dist"], legend_fig=leg, dpi=args.dpi
        )

        fig, leg = plot_detection_by_run(
            counts,
            run_order,
            shares,
            cond_colors,
            pair_colors,
            title="Unmodified-peptide detection by run and treated-only share by pair",
            subtitle=SEP.join(
                [
                    "FlashLFQ quantified (MS/MS or MBR), raw detection calls",
                    f"{len(shares)} pairs",
                    "run order aliased with condition",
                ]
            ),
            count_label="unmodified peptides quantified",
            share_label="treated-only share of discordant cells\n(unmodified peptides)",
        )
        art_runs = save_figure(
            fig, args.class_output_dir, STEMS["runs"], legend_fig=leg, dpi=args.dpi
        )

    def _script(rel: str) -> dict[str, str]:
        return {"path": rel, "sha256": sha256_of_file(_ROOT / rel)}

    provenance = {
        "findings": ["0008", "0009"],
        "scripts": {
            "runner": _script("scripts/scratch/fig_peptide_mods.py"),
            "module": _script("scripts/scratch/analysis_figures/peptide_mods_figs.py"),
            "figure_io": _script("scripts/promoted/common/figures/figure_io.py"),
            "colors": _script("scripts/promoted/common/figures/colors.py"),
        },
        "commit": _git_head(),
        "commit_note": (
            "figure scripts uncommitted at render time; script sha256 pinned"
        ),
        "data_version": summary["data_version"],
        "analysis_script": {
            "path": summary["script"],
            "sha256": summary["script_sha256"],
        },
        "inputs": {
            n: {"path": str(d / n), "sha256": sha256_of_file(d / n)} for n in INPUTS
        },
        "params": {
            "detection_categories": (
                "FlashLFQ MSMS; MBR; MSMSIdentifiedButNotQuantified + "
                "MSMSAmbiguousPeakfinding merged as 'MS/MS ID, not quantified'; "
                "NotDetected blank"
            ),
            "heatmap_row_order": "Cys first, then Tyr/Trp; n_treated_quantified desc, "
            "n_control_quantified asc, treated_psms_total desc, feature asc",
            "discordance_ci": "Clopper-Pearson 95% (recomputed; equals analysis CI)",
            "pair_level_or": (
                "discordance_class_tests pair_level_or "
                "(raloxifene-adduct vs unmodified, all residues)"
            ),
            "class_shift": "pair_shift_* of class_log2fc_tests (mean(class) - "
            "mean(unmodified) paired diff per pair; t CI df 3; BH over classes); "
            "per-pair values parsed from 4-dp strings",
            "modified_minus_reference": "per-pair mean over 123 oxidized forms of "
            "diff_<pair> in modified_vs_reference.tsv; CI/q from class tests",
            "distribution": "scipy gaussian_kde (Scott bandwidth) per class, x view "
            f"+/-{args.dist_xlim}",
            "detection_counts": "discordance_by_pair quantified/unmodified "
            "n_detected_control/treated",
            "registry": str(reg),
            "dpi": args.dpi,
        },
        "color_map": {
            "condition": cond_colors,
            "candidate_pair": pair_colors,
            "detection_type": det_raw,
            "mod_class": class_colors,
            "discordance_groups": group_colors,
        },
        "key_numbers": {
            "heatmap": {k: v for k, v in heat_stats.items() if k != "row_order"},
            "heatmap_row_order": heat_stats["row_order"],
            "discordance": disc_numbers,
            "cyp3a4": cyp_stats,
            "class_shift": shift_numbers,
            "distributions": {
                c: {
                    k: (round(v, 4) if isinstance(v, float) else v)
                    for k, v in s.items()
                }
                for c, s in dist_stats.items()
            },
            "detection_by_run": {
                "counts": counts,
                "shares": {
                    s.pair: {
                        "k": s.k,
                        "n": s.n,
                        "prop": round(s.prop, 4),
                        "ci": [round(s.ci_low, 4), round(s.ci_high, 4)],
                    }
                    for s in shares
                },
            },
        },
        "figures": {
            STEMS["heatmap"]: _artifact_paths(art_heat),
            STEMS["discordance"]: _artifact_paths(art_disc),
            STEMS["cyp3a4"]: _artifact_paths(art_cyp),
            STEMS["shift"]: _artifact_paths(art_shift),
            STEMS["dist"]: _artifact_paths(art_dist),
            STEMS["runs"]: _artifact_paths(art_runs),
        },
    }
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(json.dumps(provenance, indent=2) + "\n")
    LOG.info(
        "key numbers %s", json.dumps(provenance["key_numbers"], default=str)[:3000]
    )
    LOG.info("wrote %s", args.provenance)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
