"""Peptide-modification figures: findings 0008 (raloxifene adducts), 0009 (classes).

From-scratch Stage-4 family ``peptide-mods`` (no ``lib/figures`` template exists for
these views). Reads the plot-ready tables written by
``scripts/scratch/peptide_mod_analysis.py`` under ``results/peptide-mods/`` and renders:

  * :func:`plot_adduct_detection_heatmap` -- the 60 raloxifene-adduct peptides x 8 runs,
    cell = FlashLFQ detection type, PSM count printed in the cell;
  * :func:`plot_discordance_by_residue` -- share of discordant pair cells that are
    treated-only (Clopper-Pearson 95 % CI) for Cys adducts, Tyr/Trp adducts and the
    unmodified baseline, one panel per detection definition;
  * :func:`plot_single_adduct_peptide` -- one adduct peptide (the CYP3A4 one): per-run
    log2 intensity and PSMs, control vs raloxifene within each pair;
  * :func:`plot_class_shift_by_pair` -- per-pair class shifts (points per pair + mean
    and t CI across pairs);
  * :func:`plot_class_log2fc_distributions` -- paired peptide log2FC densities by class;
  * :func:`plot_detection_by_run` -- unmodified peptides quantified per run in
    acquisition order + per-pair treated-only share.

Every function only *plots* numbers the analysis already computed (or recomputes simple
sums that the runner cross-checks against the analysis tables); nothing is re-derived
from raw data. Colors come from the registry via the ``color_map`` arguments the runner
builds with :func:`common.figures.colors.assign_colors`. Legends are separate figures.
All inputs are validated fail-loud.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.artist import Artist
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter
from numpy.typing import NDArray
from scipy import stats

__script_meta__: dict[str, object] = {
    "task": "fig-peptide-mods",
    "kind": "module",
    "provides": [
        "SampleInfo",
        "DETECTION_LABELS",
        "detection_category",
        "row_label",
        "plot_adduct_detection_heatmap",
        "plot_discordance_by_residue",
        "plot_single_adduct_peptide",
        "plot_class_shift_by_pair",
        "plot_class_log2fc_distributions",
        "plot_detection_by_run",
    ],
    "uses": [],
    "seeded_from": None,
    "description": (
        "Figures for findings 0008/0009 (peptide modifications): adduct detection "
        "heatmap, adduct discordance by residue, CYP3A4 adduct peptide, class log2FC "
        "shifts by pair, class log2FC distributions, detection by run position. "
        "Plots precomputed results/peptide-mods tables; fail-loud input checks."
    ),
}

FloatArray = NDArray[np.float64]

MINUS = "\u2212"
EMDASH = "\u2014"
MIDDOT = "\u00b7"
CONTROL = "control"
TREATED = "raloxifene-d0"
NEUTRAL = "#000000"
GRID_GRAY = "#bdbdbd"
EDGE_GRAY = "#4d4d4d"

# FlashLFQ detection type -> the 3 displayed evidence categories (NotDetected = blank).
# MSMSAmbiguousPeakfinding (1 of 480 adduct cells) has an MS2 ID but no intensity, so it
# joins the identified-not-quantified category.
DETECTION_LABELS: dict[str, str] = {
    "MSMS": "MS/MS (quantified)",
    "MBR": "MBR transfer (quantified)",
    "MSMSIdentifiedButNotQuantified": "MS/MS ID, not quantified",
    "MSMSAmbiguousPeakfinding": "MS/MS ID, not quantified",
}
NOT_DETECTED = "NotDetected"


@dataclass(frozen=True)
class SampleInfo:
    """One run: id, arm, pair, batch, 1-based run position within its batch."""

    sample_id: str
    condition: str
    pair: str
    batch: str
    run_position: int


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def detection_category(detection_type: str) -> str | None:
    """Displayed category for a FlashLFQ detection type (``None`` = not detected)."""
    if detection_type == NOT_DETECTED:
        return None
    if detection_type not in DETECTION_LABELS:
        raise ValueError(f"Unknown FlashLFQ detection type {detection_type!r}.")
    return DETECTION_LABELS[detection_type]


def _short_entry(entry: str) -> str:
    return entry.removesuffix("_HUMAN")


def _abbrev_peptide(seq: str, max_len: int) -> str:
    if len(seq) <= max_len:
        return seq
    head = max_len // 2
    tail = max_len - head - 1
    return f"{seq[:head]}\u2026{seq[-tail:]}"


def row_label(entries: str, sites: str, base_sequence: str, max_pep: int = 15) -> str:
    """``ENTRY[/ENTRY] · site(s) · peptide``; _HUMAN stripped, peptide abbreviated."""
    ents = "/".join(sorted(_short_entry(e) for e in entries.split(";")))
    site = sites.replace(":Ralox", "").replace(" | ", "|").replace(";", "+")
    return f"{ents} {MIDDOT} {site} {MIDDOT} {_abbrev_peptide(base_sequence, max_pep)}"


def _require(frame: pd.DataFrame, cols: Sequence[str], name: str) -> None:
    missing = [c for c in cols if c not in frame.columns]
    if missing:
        raise ValueError(f"{name}: missing columns {missing}.")


def _category_legend(
    handles: Sequence[Artist], labels: Sequence[str], title: str | None = None
) -> Figure:
    fig = plt.figure(figsize=(3.2, 0.4 + 0.28 * len(handles)))
    fig.legend(
        list(handles),
        list(labels),
        loc="center left",
        frameon=False,
        title=title,
        alignment="left",
    )
    return fig


def _stacked_legend(
    sections: Sequence[tuple[str, Sequence[Artist], Sequence[str]]],
) -> Figure:
    """One legend with a bold header row per section (evenly spaced)."""
    handles: list[Artist] = []
    labels: list[str] = []
    header_rows: list[int] = []
    for i, (title, hs, ls) in enumerate(sections):
        if i:
            handles.append(Patch(facecolor="none", edgecolor="none"))
            labels.append(" ")
        header_rows.append(len(labels))
        handles.append(Patch(facecolor="none", edgecolor="none"))
        labels.append(title)
        handles.extend(hs)
        labels.extend(ls)
    fig = plt.figure(figsize=(3.6, 0.3 + 0.26 * len(handles)))
    leg = fig.legend(handles, labels, loc="center left", frameon=False)
    texts = leg.get_texts()
    for row in header_rows:
        texts[row].set_fontweight("bold")
    return fig


def _header(fig: Figure, title: str, subtitle: str, top: float, pad_in: float) -> None:
    """Place title + subtitle a fixed physical distance above the top axes edge."""
    height = fig.get_figheight()
    y_sub = top + pad_in / height
    fig.text(0.01, y_sub, subtitle, fontsize=8.5, ha="left", va="bottom")
    fig.text(
        0.01,
        y_sub + 0.2 / height,
        title,
        fontsize=11.5,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def _marker(
    color: str, marker: str = "o", size: float = 8.0, edge: str = NEUTRAL
) -> Line2D:
    return Line2D(
        [],
        [],
        ls="none",
        marker=marker,
        markersize=size,
        markerfacecolor=color,
        markeredgecolor=edge,
    )


# --------------------------------------------------------------------------- #
# 1. Adduct detection heatmap
# --------------------------------------------------------------------------- #
def plot_adduct_detection_heatmap(
    long: pd.DataFrame,
    adducts: pd.DataFrame,
    samples: Sequence[SampleInfo],
    cyp_entries: set[str],
    detection_colors: Mapping[str, str],
    condition_colors: Mapping[str, str],
    *,
    title: str,
    subtitle: str,
    row_height: float = 0.17,
) -> tuple[Figure, Figure, dict[str, object]]:
    """Adduct peptides (rows) x runs (columns, pair order), cell = detection category.

    Rows: Cys adducts first, then Tyr/Trp; within a group sorted by the number of
    treated runs quantified (desc), then control runs quantified (asc), then treated
    PSMs (desc), then feature. The in-cell number is the Limelight PSM count (> 0 only).
    CYP rows are bold.
    """
    _require(
        long,
        ["feature", "sample_id", "detection_type", "dump_psms"],
        "adduct_peptides_long",
    )
    _require(
        adducts,
        [
            "feature",
            "base_sequence",
            "entries",
            "dump_ralox_sites",
            "dump_ralox_residues",
            "n_treated_quantified",
            "n_control_quantified",
            "treated_psms_total",
        ],
        "adduct_peptides",
    )
    order_ids = [s.sample_id for s in samples]
    if not adducts["feature"].is_unique:
        raise ValueError("adduct_peptides: duplicate features.")
    feats = set(adducts["feature"])
    if set(long["feature"]) != feats:
        raise ValueError("long/wide adduct feature sets differ.")
    if len(long) != len(feats) * len(order_ids):
        raise ValueError(
            "adduct_peptides_long is not a complete feature x sample grid."
        )

    ad = adducts.copy()
    ad["is_cys"] = ad["dump_ralox_residues"].astype(str) == "C"
    ad = ad.sort_values(
        [
            "is_cys",
            "n_treated_quantified",
            "n_control_quantified",
            "treated_psms_total",
            "feature",
        ],
        ascending=[False, False, True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    n_cys = int(ad["is_cys"].sum())
    n_rows = len(ad)

    cats = list(dict.fromkeys(DETECTION_LABELS.values()))
    cat_index = {c: i + 1 for i, c in enumerate(cats)}  # 0 = not detected
    grid = np.zeros((n_rows, len(order_ids)), dtype=np.int64)
    psm = np.zeros((n_rows, len(order_ids)), dtype=np.int64)
    lookup = long.set_index(["feature", "sample_id"])
    if not lookup.index.is_unique:
        raise ValueError("adduct_peptides_long: duplicate (feature, sample) cells.")
    for i, feat in enumerate(ad["feature"]):
        for j, sid in enumerate(order_ids):
            row = lookup.loc[(feat, sid)]
            cat = detection_category(str(row["detection_type"]))
            grid[i, j] = 0 if cat is None else cat_index[cat]
            val = float(row["dump_psms"])
            if not math.isfinite(val) or val < 0 or val != int(val):
                raise ValueError(f"bad PSM count {val} for {feat} {sid}.")
            psm[i, j] = int(val)

    cmap = ListedColormap(["#ffffff", *(detection_colors[c] for c in cats)])
    height = 1.3 + row_height * n_rows
    fig, ax = plt.subplots(figsize=(6.4, height))
    ax.imshow(
        grid,
        cmap=cmap,
        vmin=-0.5,
        vmax=len(cats) + 0.5,
        aspect="auto",
        interpolation="nearest",
    )
    # cell grid
    for x in np.arange(-0.5, len(order_ids), 1.0):
        ax.axvline(x, color=GRID_GRAY, lw=0.4)
    for y in np.arange(-0.5, n_rows, 1.0):
        ax.axhline(y, color=GRID_GRAY, lw=0.4)
    # pair separators (every 2 columns) and batch separator
    for k in range(2, len(order_ids), 2):
        ax.axvline(k - 0.5, color=NEUTRAL, lw=1.0)
    ax.axhline(n_cys - 0.5, color=NEUTRAL, lw=1.4)
    for i in range(n_rows):
        for j in range(len(order_ids)):
            if psm[i, j] > 0:
                ax.text(
                    j,
                    i,
                    str(psm[i, j]),
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color=NEUTRAL,
                )

    labels = [
        row_label(str(r.entries), str(r.dump_ralox_sites), str(r.base_sequence))
        for r in ad.itertuples(index=False)
    ]
    is_cyp = [any(e in cyp_entries for e in str(x).split(";")) for x in ad["entries"]]
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(labels, fontsize=7, family="monospace")
    for tick, cyp in zip(ax.get_yticklabels(), is_cyp, strict=True):
        if cyp:
            tick.set_fontweight("bold")
    ax.set_xticks(range(len(order_ids)))
    ax.set_xticklabels(order_ids, rotation=90, fontsize=8)
    ax.tick_params(axis="both", length=0)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)

    # condition stripe + pair labels above the heatmap
    for j, s in enumerate(samples):
        ax.add_patch(
            Rectangle(
                (j - 0.5, -2.1),
                1.0,
                1.2,
                facecolor=condition_colors[s.condition],
                edgecolor="white",
                lw=0.6,
                clip_on=False,
            )
        )
    for k in range(0, len(order_ids), 2):
        pair = samples[k].pair
        ax.text(k + 0.5, -2.4, pair, ha="center", va="bottom", fontsize=7.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    # residue-group labels on the right
    for y0, y1, name in (
        (0, n_cys, f"Cys adducts (n = {n_cys})"),
        (n_cys, n_rows, f"Tyr/Trp adducts (n = {n_rows - n_cys})"),
    ):
        ax.annotate(
            "",
            xy=(len(order_ids) - 0.2, y0 - 0.4),
            xytext=(len(order_ids) - 0.2, y1 - 0.6),
            arrowprops={"arrowstyle": "-", "lw": 1.0, "color": NEUTRAL},
            annotation_clip=False,
        )
        ax.text(
            len(order_ids) + 0.05,
            (y0 + y1 - 1) / 2,
            name,
            rotation=270,
            ha="left",
            va="center",
            fontsize=9,
        )
    for text, off, kw in (
        (subtitle, 44, {"fontsize": 8.5}),
        (title, 58, {"fontsize": 11, "fontweight": "bold"}),
    ):
        ax.annotate(
            text,
            xy=(0.01, 1.0),
            xycoords=("figure fraction", "axes fraction"),
            xytext=(0, off),
            textcoords="offset points",
            ha="left",
            va="bottom",
            **kw,
        )

    handles: list[Artist] = [
        Patch(facecolor=detection_colors[c], edgecolor=EDGE_GRAY) for c in cats
    ]
    handles.append(Patch(facecolor="#ffffff", edgecolor=EDGE_GRAY))
    det_labels = [*cats, "not detected"]
    cond_handles: list[Artist] = [
        Patch(facecolor=condition_colors[c]) for c in (CONTROL, TREATED)
    ]
    legend = _stacked_legend(
        [
            ("Cell: FlashLFQ detection", handles, det_labels),
            ("Top stripe: condition", cond_handles, [CONTROL, TREATED]),
            (
                "Cell number",
                [Patch(facecolor="none", edgecolor="none")],
                ["Limelight PSMs in run (if > 0)"],
            ),
            (
                "Row label",
                [Patch(facecolor="none", edgecolor="none")],
                ["entry \u00b7 adducted site \u00b7 peptide (bold = CYP)"],
            ),
        ]
    )

    counts = {c: int((grid == cat_index[c]).sum()) for c in cats}
    counts["not detected"] = int((grid == 0).sum())
    stats_out: dict[str, object] = {
        "n_rows": n_rows,
        "n_cys": n_cys,
        "n_tyr_trp": n_rows - n_cys,
        "cell_counts": counts,
        "n_cyp_rows": int(sum(is_cyp)),
        "row_order": list(ad["feature"]),
    }
    return fig, legend, stats_out


# --------------------------------------------------------------------------- #
# 2. Discordance by residue
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class DiscordancePoint:
    """One group's treated-only share among discordant cells (k of n)."""

    group: str
    k: int
    n: int
    prop: float
    ci_low: float
    ci_high: float


@dataclass(frozen=True)
class DiscordancePanel:
    """One detection definition: the groups + the pair-level OR annotation."""

    title: str
    points: tuple[DiscordancePoint, ...]
    note: str


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Exact two-sided Clopper-Pearson interval for k successes of n."""
    if n <= 0 or not 0 <= k <= n:
        raise ValueError(f"invalid k/n = {k}/{n}.")
    lo = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def plot_discordance_by_residue(
    panels: Sequence[DiscordancePanel],
    group_colors: Mapping[str, str],
    *,
    title: str,
    subtitle: str,
) -> tuple[Figure, Figure]:
    """Per definition panel: treated-only share (CI) per group; dashed 0.5 line."""
    if not panels:
        raise ValueError("no panels.")
    groups = [p.group for p in panels[0].points]
    for panel in panels:
        if [p.group for p in panel.points] != groups:
            raise ValueError("every panel must show the same groups in the same order.")
    fig, axes = plt.subplots(
        1, len(panels), figsize=(3.0 * len(panels) + 0.6, 4.3), sharey=True
    )
    ax_list = list(np.atleast_1d(axes))
    for ax, panel in zip(ax_list, panels, strict=True):
        ax.axhline(0.5, color=EDGE_GRAY, ls="--", lw=1.0, zorder=1)
        for x, pt in enumerate(panel.points):
            ax.errorbar(
                x,
                pt.prop,
                yerr=[[pt.prop - pt.ci_low], [pt.ci_high - pt.prop]],
                fmt="o",
                color=group_colors[pt.group],
                markeredgecolor=NEUTRAL,
                ecolor=NEUTRAL,
                elinewidth=1.3,
                capsize=4,
                markersize=8,
                zorder=3,
            )
            ax.text(
                x + 0.16,
                pt.prop,
                f"{pt.k:,}/{pt.n:,}",
                ha="left",
                va="center",
                fontsize=8.5,
            )
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([g.replace(" ", "\n", 1) for g in groups], fontsize=9)
        ax.set_xlim(-0.5, len(groups) - 0.2)
        ax.set_ylim(0.0, 1.05)
        ax.set_title(panel.title, fontsize=10.5)
        ax.text(
            0.02,
            0.03,
            panel.note,
            transform=ax.transAxes,
            fontsize=8,
            ha="left",
            va="bottom",
        )
    ax_list[0].set_ylabel("treated-only share of discordant\npair \u00d7 peptide cells")
    top = 0.8
    fig.subplots_adjust(top=top, bottom=0.17, left=0.1, right=0.98, wspace=0.22)
    _header(fig, title, subtitle, top, 0.35)

    handles: list[Artist] = [_marker(group_colors[g]) for g in groups]
    handles.append(Line2D([], [], color=NEUTRAL, lw=1.3, marker="|", markersize=8))
    handles.append(Line2D([], [], color=EDGE_GRAY, ls="--", lw=1.0))
    labels = [
        *groups,
        "95% Clopper\u2013Pearson CI",
        "0.5 = no treated/control imbalance",
    ]
    return fig, _category_legend(handles, labels)


# --------------------------------------------------------------------------- #
# 3. Single adduct peptide (CYP3A4)
# --------------------------------------------------------------------------- #
def plot_single_adduct_peptide(
    long: pd.DataFrame,
    feature: str,
    samples: Sequence[SampleInfo],
    condition_colors: Mapping[str, str],
    *,
    title: str,
    subtitle: str,
) -> tuple[Figure, Figure, dict[str, object]]:
    """Two panels sharing x (pairs; control then raloxifene): log2 intensity, PSMs."""
    _require(
        long,
        ["feature", "sample_id", "detection_type", "log2_intensity", "dump_psms"],
        "adduct_peptides_long",
    )
    sub = long[long["feature"] == feature].set_index("sample_id")
    if len(sub) != len(samples) or set(sub.index) != {s.sample_id for s in samples}:
        raise ValueError(f"{feature}: expected exactly one row per sample.")
    pairs = list(dict.fromkeys(s.pair for s in samples))
    offset = {CONTROL: -0.18, TREATED: 0.18}
    fig, (ax_i, ax_p) = plt.subplots(
        2,
        1,
        figsize=(5.2, 4.6),
        sharex=True,
        gridspec_kw={"height_ratios": [1.6, 1.0], "hspace": 0.12},
    )
    detected_log2: list[float] = []
    for s in samples:
        row = sub.loc[s.sample_id]
        x = pairs.index(s.pair) + offset[s.condition]
        color = condition_colors[s.condition]
        det = str(row["detection_type"])
        cat = detection_category(det)
        l2 = float(row["log2_intensity"])
        if cat is not None and "quantified)" in cat:
            if not math.isfinite(l2):
                raise ValueError(f"{s.sample_id}: quantified but no intensity.")
            detected_log2.append(l2)
            ax_i.plot(
                x, l2, "o", color=color, markeredgecolor=NEUTRAL, markersize=9, zorder=3
            )
            ax_i.text(
                x,
                l2 + 0.35,
                "MS/MS" if det == "MSMS" else "MBR",
                ha="center",
                va="bottom",
                fontsize=7.5,
            )
        else:
            if math.isfinite(l2):
                raise ValueError(
                    f"{s.sample_id}: intensity present but not quantified."
                )
        psms = float(row["dump_psms"])
        ax_p.bar(x, psms, width=0.32, color=color, edgecolor=NEUTRAL, lw=0.6)
        if psms == 0:
            ax_p.text(x, 0.05, "0", ha="center", va="bottom", fontsize=8)
    if not detected_log2:
        raise ValueError(f"{feature}: never quantified.")
    lo = math.floor(min(detected_log2)) - 2.0
    hi = math.ceil(max(detected_log2)) + 1.0
    ax_i.set_ylim(lo, hi)
    # not-quantified runs marked at the floor
    for s in samples:
        row = sub.loc[s.sample_id]
        if not math.isfinite(float(row["log2_intensity"])):
            x = pairs.index(s.pair) + offset[s.condition]
            ax_i.plot(
                x,
                lo + 0.45,
                marker="x",
                color=condition_colors[s.condition],
                markersize=8,
                mew=2.0,
                ls="none",
            )
    ax_i.axhline(lo + 0.9, color=GRID_GRAY, lw=0.6)
    ax_i.text(-0.62, lo + 0.45, "not\nquant.", fontsize=7, ha="left", va="center")
    ax_i.set_ylabel("log2 intensity\n(FlashLFQ, raw, a.u.)")
    ax_p.set_ylabel("Limelight PSMs")
    ax_p.yaxis.get_major_locator().set_params(integer=True)
    batch_of = {s.pair: s.batch.removeprefix("B")[:4] for s in samples}
    ax_p.set_xticks(range(len(pairs)))
    ax_p.set_xticklabels([f"{p}\n({batch_of[p]})" for p in pairs], fontsize=9)
    ax_p.set_xlim(-0.65, len(pairs) - 0.35)
    ax_p.set_xlabel("candidate pair (control left, raloxifene-d0 right)")
    top = 0.88
    fig.subplots_adjust(top=top)
    _header(fig, title, subtitle, top, 0.12)

    handles: list[Artist] = [
        Patch(facecolor=condition_colors[c], edgecolor=NEUTRAL)
        for c in (CONTROL, TREATED)
    ]
    handles.append(_marker("#ffffff", "o", edge=NEUTRAL))
    handles.append(
        Line2D([], [], ls="none", marker="x", color=NEUTRAL, markersize=8, mew=2.0)
    )
    labels = [
        CONTROL,
        TREATED,
        "quantified (label = detection type)",
        "not quantified in run",
    ]
    stats_out: dict[str, object] = {
        "feature": feature,
        "per_sample": {
            s.sample_id: {
                "detection_type": str(sub.loc[s.sample_id, "detection_type"]),
                "log2_intensity": (
                    None
                    if not math.isfinite(float(sub.loc[s.sample_id, "log2_intensity"]))
                    else round(float(sub.loc[s.sample_id, "log2_intensity"]), 4)
                ),
                "psms": int(float(sub.loc[s.sample_id, "dump_psms"])),
            }
            for s in samples
        },
    }
    return fig, _category_legend(handles, labels), stats_out


# --------------------------------------------------------------------------- #
# 4. Class shift by pair
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ShiftGroup:
    """One contrast: per-pair values + mean and its CI across pairs."""

    label: str
    per_pair: Mapping[str, float]
    mean: float
    ci_low: float
    ci_high: float
    note: str


def plot_class_shift_by_pair(
    groups: Sequence[ShiftGroup],
    pair_colors: Mapping[str, str],
    *,
    title: str,
    subtitle: str,
    ylabel: str,
) -> tuple[Figure, Figure]:
    """Points per pair (pair color) + mean +/- CI (black) per contrast; zero line.

    All contrasts are drawn with identical styling (no contrast is emphasized).
    """
    pairs = list(pair_colors)
    fig, ax = plt.subplots(figsize=(1.9 * len(groups) + 1.6, 4.4))
    ax.axhline(0.0, color=EDGE_GRAY, lw=1.0, zorder=1)
    jitter = np.linspace(-0.15, 0.15, len(pairs))
    for x, g in enumerate(groups):
        if set(g.per_pair) != set(pairs):
            raise ValueError(f"{g.label}: per-pair values must cover {pairs}.")
        for dx, p in zip(jitter, pairs, strict=True):
            ax.plot(
                x - 0.12 + dx * 0.6,
                g.per_pair[p],
                "o",
                color=pair_colors[p],
                markeredgecolor=NEUTRAL,
                markersize=7,
                zorder=3,
            )
        ax.errorbar(
            x + 0.2,
            g.mean,
            yerr=[[g.mean - g.ci_low], [g.ci_high - g.mean]],
            fmt="D",
            color=NEUTRAL,
            markersize=7,
            capsize=4,
            elinewidth=1.4,
            zorder=4,
        )
        ax.text(
            x + 0.3,
            g.mean,
            f"{g.mean:+.3f}".replace("-", MINUS),
            ha="left",
            va="center",
            fontsize=8.5,
        )
        ax.text(
            x,
            1.0,
            g.note,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g.label for g in groups], fontsize=9)
    ax.set_xlim(-0.6, len(groups) - 0.25)
    ax.set_ylabel(ylabel)
    top = 0.78
    fig.subplots_adjust(top=top, bottom=0.16, left=0.14, right=0.98)
    _header(fig, title, subtitle, top, 0.5)

    handles: list[Artist] = [_marker(pair_colors[p], size=7) for p in pairs]
    handles.append(Line2D([], [], ls="-", marker="D", color=NEUTRAL, markersize=6))
    labels = [*pairs, "mean over 4 pairs \u00b1 95% CI (t, df 3)"]
    return fig, _category_legend(handles, labels, title="candidate pair")


# --------------------------------------------------------------------------- #
# 5. Class log2FC distributions
# --------------------------------------------------------------------------- #
def plot_class_log2fc_distributions(
    values: Mapping[str, FloatArray],
    class_colors: Mapping[str, str],
    *,
    title: str,
    subtitle: str,
    xlabel: str,
    xlim: tuple[float, float] = (-1.25, 1.25),
) -> tuple[Figure, Figure, dict[str, dict[str, float | int]]]:
    """Gaussian-KDE (Scott) density per class, normalized per class; median lines.

    Also returns per-class n / median / IQR. Values outside ``xlim`` still enter the
    KDE and the statistics; they are only clipped from view (count reported).
    """
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    grid = np.linspace(xlim[0], xlim[1], 600)
    stats_out: dict[str, dict[str, float | int]] = {}
    ymax = 0.0
    for cls, v in values.items():
        arr = np.asarray(v, dtype=np.float64)
        if arr.size < 10 or not np.all(np.isfinite(arr)):
            raise ValueError(f"{cls}: need >= 10 finite values.")
        dens = stats.gaussian_kde(arr)(grid)
        ymax = max(ymax, float(dens.max()))
        color = class_colors[cls]
        ax.plot(grid, dens, color=color, lw=2.0)
        ax.fill_between(grid, dens, color=color, alpha=0.12, lw=0)
        med = float(np.median(arr))
        ax.axvline(med, color=color, ls="--", lw=1.4)
        stats_out[cls] = {
            "n": int(arr.size),
            "median": med,
            "q25": float(np.quantile(arr, 0.25)),
            "q75": float(np.quantile(arr, 0.75)),
            "n_outside_view": int(((arr < xlim[0]) | (arr > xlim[1])).sum()),
        }
    ax.axvline(0.0, color=EDGE_GRAY, lw=0.8)
    ax.set_xlim(*xlim)
    ax.set_ylim(0, ymax * 1.08)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density (per class)")
    top = 0.86
    fig.subplots_adjust(top=top, bottom=0.15, left=0.11, right=0.98)
    _header(fig, title, subtitle, top, 0.12)

    handles: list[Artist] = []
    labels: list[str] = []
    for cls, s in stats_out.items():
        handles.append(Line2D([], [], color=class_colors[cls], lw=2.0))
        med_txt = f"{s['median']:+.3f}".replace("-", MINUS)
        labels.append(f"{cls} (n = {int(s['n']):,}; median {med_txt})")
    handles.append(Line2D([], [], color=NEUTRAL, ls="--", lw=1.4))
    labels.append("dashed = class median")
    return fig, _category_legend(handles, labels, title="peptide class"), stats_out


# --------------------------------------------------------------------------- #
# 6. Detection by run position
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PairShare:
    """Treated-only share of discordant cells for one pair."""

    pair: str
    batch: str
    k: int
    n: int
    prop: float
    ci_low: float
    ci_high: float


def plot_detection_by_run(
    counts: Mapping[str, int],
    samples_in_run_order: Sequence[SampleInfo],
    shares: Sequence[PairShare],
    condition_colors: Mapping[str, str],
    pair_colors: Mapping[str, str],
    *,
    title: str,
    subtitle: str,
    count_label: str,
    share_label: str,
) -> tuple[Figure, Figure]:
    """(a) peptides quantified per run in acquisition order, bars by condition;
    (b) per-pair treated-only share (CI) in pair colors; dashed 0.5."""
    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(10.0, 4.4),
        gridspec_kw={"width_ratios": [1.9, 1.0], "wspace": 0.32},
    )
    xs: list[float] = []
    x = 0.0
    prev_batch: str | None = None
    batch_spans: dict[str, list[float]] = {}
    for s in samples_in_run_order:
        if prev_batch is not None and s.batch != prev_batch:
            x += 0.8
        xs.append(x)
        batch_spans.setdefault(s.batch, []).append(x)
        prev_batch = s.batch
        x += 1.0
    for xi, s in zip(xs, samples_in_run_order, strict=True):
        n = counts[s.sample_id]
        ax_a.bar(
            xi,
            n,
            width=0.75,
            color=condition_colors[s.condition],
            edgecolor=NEUTRAL,
            lw=0.6,
        )
        ax_a.text(xi, n, f"{n:,}", ha="center", va="bottom", fontsize=7.5)
    ax_a.set_xticks(xs)
    ax_a.set_xticklabels(
        [f"{s.run_position}\n{s.sample_id}" for s in samples_in_run_order],
        fontsize=8.5,
    )
    ymax = max(counts.values())
    ax_a.set_ylim(0, ymax * 1.12)
    for batch, bx in batch_spans.items():
        ax_a.text(
            (min(bx) + max(bx)) / 2,
            ymax * 1.10,
            batch.removeprefix("B"),
            ha="center",
            va="top",
            fontsize=9,
        )
    ax_a.set_ylabel(count_label)
    ax_a.set_xlabel("run position within batch / sample")
    ax_a.yaxis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{v:,.0f}"))
    ax_a.set_title("a  per run, acquisition order", loc="left", fontsize=10.5)

    for xb, sh in enumerate(shares):
        ax_b.errorbar(
            xb,
            sh.prop,
            yerr=[[sh.prop - sh.ci_low], [sh.ci_high - sh.prop]],
            fmt="o",
            color=pair_colors[sh.pair],
            markeredgecolor=NEUTRAL,
            ecolor=NEUTRAL,
            capsize=4,
            markersize=8,
        )
        ax_b.text(
            xb + 0.2, sh.prop, f"{sh.prop:.2f}", ha="left", va="center", fontsize=8.5
        )
    ax_b.axhline(0.5, color=EDGE_GRAY, ls="--", lw=1.0)
    ax_b.set_xticks(range(len(shares)))
    ax_b.set_xticklabels(
        [
            f"{s.pair.removeprefix('P')}\n({s.batch.removeprefix('B')[:4]})"
            for s in shares
        ],
        fontsize=8,
    )
    ax_b.set_xlim(-0.5, len(shares) - 0.3)
    ax_b.set_ylim(0.0, 1.0)
    ax_b.set_ylabel(share_label)
    ax_b.set_title("b  per pair", loc="left", fontsize=10.5)
    ax_b.set_xlabel("candidate pair (batch year)")
    top = 0.82
    fig.subplots_adjust(top=top, bottom=0.2, left=0.08, right=0.98)
    _header(fig, title, subtitle, top, 0.35)

    legend = _stacked_legend(
        [
            (
                "a  bar color: condition",
                [
                    Patch(facecolor=condition_colors[c], edgecolor=NEUTRAL)
                    for c in (CONTROL, TREATED)
                ],
                [CONTROL, TREATED],
            ),
            (
                "b  point color: candidate pair",
                [
                    *(_marker(pair_colors[p], size=7) for p in pair_colors),
                    Line2D([], [], color=NEUTRAL, lw=1.3, marker="|", markersize=8),
                    Line2D([], [], color=EDGE_GRAY, ls="--"),
                ],
                [
                    *pair_colors,
                    "95% Clopper\u2013Pearson CI",
                    "0.5 = no treated/control imbalance",
                ],
            ),
        ]
    )
    return fig, legend
