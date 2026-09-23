"""Volcano figure for a differential-abundance contrast (project copy of volcano@0.2).

Seeded from the plugin template ``lib/figures/volcano.py`` (v0.2). What it draws (one
figure): the per-feature **effect size** (log2 fold change) on x vs ``-log10(BH q)`` on
y, with three-way significance coloring — not significant (gray background), significant
& up, significant & down — and the **hit counts carried in the separate legend**.

Project deviations from the template (each deliberate, each documented):

* **Imports** point at the project packages (``common.figures.colors``,
  ``common.figures.figure_io``; ``analysis.differential_abundance`` from
  ``scripts/scratch``).
* **``annotate_scope``** — the template labels only *significant* hits, so a zero-hit
  contrast gets no labels at all. ``annotate_scope="all"`` instead labels the
  ``annotate_top`` smallest-q features among **all** testable ones (ties on a BH
  plateau broken by raw p, then by input order). Labelled points get a thin black ring
  so the eye finds the dot each label belongs to, and the legend gains a
  ``labelled: k smallest q`` entry that states how many of them are hits, so a labelled
  non-significant feature is never mistaken for a hit.
* **Threshold visibility** — the y-axis always extends above the ``q = fdr`` guide (the
  template relied on autoscaling); the guide's value and the smallest q
  (``min q = …``) are printed in the right margin, so a zero-hit volcano shows where
  the bar is and how near the best feature came.
* **Labels of non-hits stay below the bar** — when no labelled feature is a hit, all
  label boxes are confined below the threshold line (verified after drawing, fail
  loud), so a non-significant label can never read as significant.
* **Symmetric x-axis** about 0 (``symmetric_x=True``) so up/down asymmetry reads
  honestly.
* **Coincident labels are merged** — labelled points whose rings overlap (closer than
  ``_MERGE_WITHIN`` of the axes, single linkage) share one stacked label (members in
  left-to-right order) with one leader to their centroid, because separate leaders
  to overlapping rings cannot say which dot is which.
* :func:`volcano_from_result` takes ``label_map`` (feature -> display label; e.g. entry
  names) and passes the raw p for label tie-breaking.

* **Strict hit rule** ``q < fdr`` (template: ``q <= fdr``), matching the project's DE
  summary (``hits["q<0.05"]``) so the legend counts and summary.json agree by
  definition.

Significance: a feature is a hit when ``q < fdr`` (and, with ``effect_threshold``,
``|effect| >= effect_threshold``). A ``q`` that underflowed to 0 is floored to the
smallest positive double so ``-log10`` stays finite.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from functools import partial
from itertools import permutations
from pathlib import Path
from typing import Literal, cast

import matplotlib.pyplot as plt
import numpy as np
import textalloc as ta
from analysis.differential_abundance import DifferentialAbundanceResult
from common.figures.colors import DEFAULT_REGISTRY_PATH, assign_colors
from common.figures.figure_io import publication_style
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.text import Text

__script_meta__: dict[str, object] = {
    "template": {"name": "volcano", "version": "0.2"},
    "kind": "module",
    "provides": [
        "VolcanoCounts",
        "VolcanoPlot",
        "plot_volcano",
        "volcano_from_result",
    ],
    "uses": [
        "analysis.differential_abundance",
        "common.figures.colors",
        "common.figures.figure_io",
    ],
    "seeded_from": "volcano@0.2",
    "description": (
        "Volcano figure for a differential-abundance contrast: log2 fold change vs "
        "-log10(BH q), NS gray / up / down registry colors, hit counts in a separate "
        "legend. Project additions: annotate_scope='all' labels the k smallest-q "
        "features even when none is a hit (ringed points + a legend entry stating how "
        "many are hits), the q-threshold guide is always in view and value-labelled, "
        "symmetric x-axis, label_map for display names. Collision-free labels via "
        "textalloc. Fail-loud."
    ),
}

DEFAULT_SIGNIFICANCE_CATEGORY = "Significance"
_UP_KEY = "up"
_DOWN_KEY = "down"
_NS_KEY = "NS"

AnnotateScope = Literal["hits", "all"]

# textalloc tuning (template values).
_LABEL_FONTSIZE = 8
_LEADER_COLOR = "0.45"
_LEADER_WIDTH = 0.6
_LABEL_CANDIDATES = 900
_LABEL_MARGIN = 0.012
_LABEL_MIN_DISTANCE = 0.018
_LABEL_MAX_DISTANCE = 0.28
# Keep label boxes a hair inside the top AND bottom of the axes (off the spines).
_LABEL_TOP_PAD_FRAC = 0.02
# Labelled points closer than this (axes fraction, ~one ring diameter) share one
# stacked label + one leader: separate leaders to overlapping rings are ambiguous.
_MERGE_WITHIN = 0.012
# Leader lengths (axes fraction) tried after the requested one until a layout has no
# leader-through-label crossing and no overlapping label boxes.
_LABEL_DISTANCE_FALLBACKS: tuple[float, ...] = (0.2, 0.28, 0.35, 0.45, 0.6)
# x-axis widening factors for the column fallback. Widening does not create room for
# text (label widths in data units scale with the range), so only 1.0 is used; kept as
# a knob for layouts dominated by leader geometry rather than label width.
_COLUMN_X_EXPANSIONS: tuple[float, ...] = (1.0,)
# Headroom above the highest point / threshold guide (fraction of that height).
_Y_HEADROOM = 0.12
# When no labelled feature is a hit, label boxes are confined below the q-threshold
# line by this margin (axes fraction) so a non-significant label never reads as a hit.
_LABEL_THRESHOLD_GAP_FRAC = 0.03
# Ring around labelled points (points^2); small so a ring just under the threshold
# line does not touch it.
_RING_SIZE = 14
_X_PAD = 1.08


@dataclass(frozen=True)
class VolcanoCounts:
    """The three significance bucket counts drawn (and shown in the legend)."""

    up: int
    down: int
    ns: int


@dataclass
class VolcanoPlot:
    """A rendered volcano figure plus its companion legend figure.

    Attributes
    ----------
    figure, legend_figure:
        The main figure (no baked legend) and the standalone legend figure.
    counts:
        The :class:`VolcanoCounts` actually drawn.
    color_map:
        ``{"up": hex, "down": hex, "NS": hex}`` used.
    labelled_index:
        Positions (into the input arrays) of the labelled features, in rank order.
    labelled_text:
        The label strings drawn, aligned with ``labelled_index``.
    n_labelled_hits:
        How many labelled features are significant hits.
    label_groups:
        Positions (into ``labelled_index``) sharing one stacked label because their
        points coincide at print scale (see :func:`_group_coincident`).
    label_box_top:
        Highest label-box top in data units (``-log10 q``); ``None`` if no labels.
    label_conflicts:
        Leader-through-label crossings + overlapping label pairs in the kept layout.
    """

    figure: Figure
    legend_figure: Figure
    counts: VolcanoCounts
    color_map: dict[str, str]
    labelled_index: list[int]
    labelled_text: list[str]
    n_labelled_hits: int
    label_groups: list[list[int]]
    label_box_top: float | None
    label_conflicts: int


def _safe_neg_log10(qvalues: np.ndarray) -> np.ndarray:
    """``-log10(q)`` with zeros floored to the tiniest positive double; ``NaN`` kept."""
    q = np.asarray(qvalues, dtype=float)
    floor = float(np.finfo(float).tiny)
    safe = np.where(np.isnan(q), np.nan, np.maximum(q, floor))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.asarray(-np.log10(safe), dtype=float)


def _rank_for_labels(
    q: np.ndarray,
    p: np.ndarray | None,
    candidates: np.ndarray,
    k: int,
) -> np.ndarray:
    """Indices of the ``k`` smallest-q candidates; ties by p, then input order."""
    idx = np.flatnonzero(candidates)
    if idx.size == 0 or k == 0:
        return np.asarray([], dtype=int)
    tiebreak = p[idx] if p is not None else np.zeros(idx.size)
    # lexsort: last key is primary -> (q, then p, then position).
    order = np.lexsort((idx, tiebreak, q[idx]))
    return np.asarray(idx[order][:k], dtype=int)


# --------------------------------------------------------------------------- #
# Core (decoupled — plain arrays)
# --------------------------------------------------------------------------- #
def plot_volcano(
    effect: np.ndarray,
    qvalues: np.ndarray,
    *,
    pvalues: np.ndarray | None = None,
    fdr: float = 0.05,
    effect_threshold: float | None = None,
    effect_label: str = "log2 fold change",
    direction_labels: tuple[str, str] = ("up", "down"),
    labels: np.ndarray | None = None,
    annotate_top: int = 0,
    annotate_scope: AnnotateScope = "hits",
    symmetric_x: bool = True,
    label_max_distance: float = _LABEL_MAX_DISTANCE,
    category: str = DEFAULT_SIGNIFICANCE_CATEGORY,
    title: str | None = None,
    legend_title: str = "Significance",
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
) -> VolcanoPlot:
    """Render a volcano from per-feature effects and BH q-values.

    Parameters
    ----------
    effect, qvalues:
        ``(n_features,)`` effect sizes (x) and BH q-values. ``NaN`` in either drops the
        feature from the plot and the counts.
    pvalues:
        Optional raw p-values (same length) — used only to break ties between equal
        q-values when ranking features for labels. Must lie in ``[0, 1]``.
    fdr:
        q threshold (hit iff ``q < fdr``) and the horizontal guide (default 0.05).
    effect_threshold:
        Optional ``|effect|`` gate (vertical guides drawn). ``None`` = q only.
    effect_label:
        x-axis label.
    direction_labels:
        ``(positive, negative)`` legend nouns for the effect sign.
    labels:
        Optional ``(n_features,)`` display names, required when ``annotate_top > 0``.
    annotate_top:
        Number of features to label (default ``0``).
    annotate_scope:
        ``"hits"`` (template behaviour: label only significant features) or ``"all"``
        (label the smallest-q features whether or not they are hits).
    symmetric_x:
        Make the x-axis symmetric about 0 (default ``True``).
    label_max_distance:
        textalloc's longest allowed leader (axes fraction; template value 0.28). Raise
        it for a crowded cluster of near-coincident labelled points.
    category, title, legend_title, registry_path, persist_colors:
        As in the template.
    """
    eff = np.asarray(effect, dtype=float)
    q = np.asarray(qvalues, dtype=float)
    if eff.shape != q.shape:
        raise ValueError(
            f"effect and qvalues must have the same shape; got {eff.shape} and "
            f"{q.shape}."
        )
    if eff.ndim != 1:
        raise ValueError(f"effect/qvalues must be 1D; got {eff.ndim}D.")
    finite_q = q[np.isfinite(q)]
    if finite_q.size and (finite_q.min() < 0.0 or finite_q.max() > 1.0):
        raise ValueError("qvalues must lie in [0, 1].")
    p_arr: np.ndarray | None = None
    if pvalues is not None:
        p_arr = np.asarray(pvalues, dtype=float)
        if p_arr.shape != eff.shape:
            raise ValueError(
                f"pvalues must match effect length; got {p_arr.shape} vs {eff.shape}."
            )
        fp = p_arr[np.isfinite(p_arr)]
        if fp.size and (fp.min() < 0.0 or fp.max() > 1.0):
            raise ValueError("pvalues must lie in [0, 1].")
    if not 0.0 < fdr < 1.0:
        raise ValueError(f"fdr must be in (0, 1); got {fdr}.")
    if effect_threshold is not None and effect_threshold < 0:
        raise ValueError(
            f"effect_threshold must be non-negative when given; got {effect_threshold}."
        )
    if annotate_top < 0:
        raise ValueError(f"annotate_top must be >= 0; got {annotate_top}.")
    if annotate_scope not in ("hits", "all"):
        raise ValueError(
            f"annotate_scope must be 'hits' or 'all'; got {annotate_scope}"
        )
    if annotate_top > 0 and labels is None:
        raise ValueError("annotate_top > 0 requires labels (feature names).")
    name_arr = np.asarray(labels) if labels is not None else None
    if name_arr is not None and name_arr.shape != eff.shape:
        raise ValueError(
            f"labels must match effect length; got {name_arr.shape} vs {eff.shape}."
        )

    testable = np.isfinite(q) & np.isfinite(eff)
    if not testable.any():
        raise ValueError("No testable feature (all effect/q are NaN).")
    gate = np.abs(eff) >= effect_threshold if effect_threshold is not None else True
    sig = testable & (q < fdr) & gate
    up = sig & (eff > 0)
    down = sig & (eff < 0)
    ns = testable & ~sig
    counts = VolcanoCounts(up=int(up.sum()), down=int(down.sum()), ns=int(ns.sum()))
    if counts.up + counts.down + counts.ns != int(testable.sum()):
        raise AssertionError("bucket counts do not partition the testable features.")

    neg_log_q = _safe_neg_log10(q)
    candidates = sig if annotate_scope == "hits" else testable
    label_idx = (
        _rank_for_labels(q, p_arr, candidates, annotate_top)
        if annotate_top > 0
        else np.asarray([], dtype=int)
    )
    label_text = [str(name_arr[i]) for i in label_idx] if name_arr is not None else []
    n_labelled_hits = int(sig[label_idx].sum()) if label_idx.size else 0

    with publication_style():
        fig, ax = plt.subplots(figsize=(8, 6.5))
        try:
            color_map = assign_colors(
                category,
                [_UP_KEY, _DOWN_KEY],
                registry_path=registry_path,
                background_values=[_NS_KEY],
                persist=persist_colors,
            )
            _draw_volcano(
                ax,
                eff=eff,
                neg_log_q=neg_log_q,
                testable=testable,
                up=up,
                down=down,
                ns=ns,
                color_map=color_map,
                fdr=fdr,
                effect_threshold=effect_threshold,
                effect_label=effect_label,
                symmetric_x=symmetric_x,
                min_q=float(np.nanmin(q[testable])),
            )
            label_groups: list[list[int]] = []
            label_box_top: float | None = None
            label_conflicts = 0
            if label_idx.size:
                label_groups, label_box_top, label_conflicts = _annotate(
                    ax,
                    eff,
                    neg_log_q,
                    testable,
                    label_idx,
                    label_text,
                    label_max_distance,
                    threshold_y=float(-np.log10(fdr)),
                    cap_below_threshold=n_labelled_hits == 0,
                )
            if title is not None:
                fig.suptitle(title, fontsize=12, weight="bold")
            legend_figure = _legend_figure(
                color_map,
                counts,
                direction_labels,
                legend_title,
                fdr=fdr,
                n_labelled=int(label_idx.size),
                n_labelled_hits=n_labelled_hits,
            )
        except BaseException:
            plt.close(fig)
            raise

    return VolcanoPlot(
        figure=fig,
        legend_figure=legend_figure,
        counts=counts,
        color_map=color_map,
        labelled_index=[int(i) for i in label_idx],
        labelled_text=label_text,
        n_labelled_hits=n_labelled_hits,
        label_groups=label_groups,
        label_box_top=label_box_top,
        label_conflicts=label_conflicts,
    )


# --------------------------------------------------------------------------- #
# Family bridge — read a DifferentialAbundanceResult
# --------------------------------------------------------------------------- #
def volcano_from_result(
    result: DifferentialAbundanceResult,
    *,
    term: str | None = None,
    fdr: float = 0.05,
    effect_threshold: float | None = None,
    effect_label: str | None = None,
    annotate_top: int = 0,
    annotate_scope: AnnotateScope = "hits",
    label_map: Mapping[str, str] | None = None,
    label_max_distance: float = _LABEL_MAX_DISTANCE,
    direction_labels: tuple[str, str] | None = None,
    title: str | None = None,
    legend_title: str = "Significance",
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
) -> tuple[VolcanoPlot, list[str]]:
    """Volcano for one contrast term of a :class:`DifferentialAbundanceResult`.

    Returns the plot and the result-table feature ids in plotted (input) order, so the
    caller can map ``labelled_index`` back to features. ``label_map`` maps feature id ->
    display label; every feature must be present when labels are requested.
    """
    term = _resolve_term(result, term)
    rows = result.table[result.table["term"] == term]
    if rows["feature"].duplicated().any():
        raise ValueError(f"term {term!r}: duplicated feature ids in the result table.")
    features = [str(f) for f in rows["feature"]]
    eff = rows["effect"].to_numpy(dtype=float)
    q = rows["q"].to_numpy(dtype=float)
    p = rows["p"].to_numpy(dtype=float)
    labels_arg: np.ndarray | None = None
    if annotate_top > 0:
        if label_map is None:
            labels_arg = np.asarray(features)
        else:
            missing = [f for f in features if f not in label_map]
            if missing:
                raise KeyError(
                    f"label_map lacks {len(missing)} feature(s), e.g. {missing[:3]}."
                )
            labels_arg = np.asarray([label_map[f] for f in features])
    dlabels = direction_labels if direction_labels is not None else ("up", "down")
    plot = plot_volcano(
        eff,
        q,
        pvalues=p,
        fdr=fdr,
        effect_threshold=effect_threshold,
        effect_label=effect_label if effect_label is not None else result.effect_label,
        direction_labels=dlabels,
        labels=labels_arg,
        annotate_top=annotate_top,
        annotate_scope=annotate_scope,
        label_max_distance=label_max_distance,
        title=title if title is not None else term,
        legend_title=legend_title,
        registry_path=registry_path,
        persist_colors=persist_colors,
    )
    return plot, features


def _resolve_term(result: DifferentialAbundanceResult, term: str | None) -> str:
    """Pick the contrast term to plot; fail loud on an ambiguous or unknown choice."""
    terms = result.contrast_terms
    if term is None:
        if len(terms) != 1:
            raise ValueError(
                f"This result has {len(terms)} contrast terms {list(terms)}; pass "
                f"term= to choose which to plot."
            )
        return terms[0]
    if term not in terms:
        raise ValueError(
            f"term {term!r} is not a contrast term; choose one of {list(terms)}."
        )
    return term


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #
def _draw_volcano(
    ax: Axes,
    *,
    eff: np.ndarray,
    neg_log_q: np.ndarray,
    testable: np.ndarray,
    up: np.ndarray,
    down: np.ndarray,
    ns: np.ndarray,
    color_map: dict[str, str],
    fdr: float,
    effect_threshold: float | None,
    effect_label: str,
    symmetric_x: bool,
    min_q: float,
) -> None:
    """Scatter the three buckets, the threshold guides, and fix the axis limits."""
    ax.scatter(
        eff[ns],
        neg_log_q[ns],
        c=color_map[_NS_KEY],
        s=10,
        alpha=0.6,
        edgecolors="none",
        zorder=1,
        rasterized=bool(ns.sum() > 5000),
    )
    for mask, key in ((down, _DOWN_KEY), (up, _UP_KEY)):
        ax.scatter(
            eff[mask],
            neg_log_q[mask],
            c=color_map[key],
            s=22,
            alpha=0.85,
            edgecolors="white",
            linewidths=0.3,
            zorder=3,
        )
    threshold_y = float(-np.log10(fdr))
    ax.axvline(0.0, color="gray", linestyle=":", linewidth=1, zorder=0)
    ax.axhline(threshold_y, color="0.3", linestyle="--", linewidth=1, zorder=0)
    if effect_threshold is not None and effect_threshold > 0:
        for x in (-effect_threshold, effect_threshold):
            ax.axvline(x, color="gray", linestyle="--", linewidth=0.8, zorder=0)

    y_max = max(float(np.nanmax(neg_log_q[testable])), threshold_y)
    ax.set_ylim(0.0, y_max * (1.0 + _Y_HEADROOM))
    if symmetric_x:
        x_abs = float(np.nanmax(np.abs(eff[testable])))
        if effect_threshold is not None:
            x_abs = max(x_abs, effect_threshold)
        ax.set_xlim(-x_abs * _X_PAD, x_abs * _X_PAD)
    # Value label in the right margin, outside the axes, so no point label can
    # collide with it.
    ax.text(
        1.01,
        threshold_y,
        f"q = {fdr:g}",
        transform=ax.get_yaxis_transform(),
        ha="left",
        va="center",
        clip_on=False,
        fontsize=9,
        color="0.3",
    )
    # Smallest q beneath it: makes a near miss (or the distance to the bar) explicit.
    ax.annotate(
        f"min q = {min_q:.2g}",
        xy=(1.01, threshold_y),
        xycoords=ax.get_yaxis_transform(),
        xytext=(0, -11),
        textcoords="offset points",
        ha="left",
        va="center",
        annotation_clip=False,
        fontsize=9,
        color="0.3",
    )
    ax.set_xlabel(effect_label)
    ax.set_ylabel(r"$-\log_{10}$(BH $q$)")
    ax.grid(True, alpha=0.25)


def _layout_conflicts(
    ax: Axes,
    renderer: object,
    text_objs: list[Text | None],
    line_objs: list[Line2D | None],
    points_data: np.ndarray,
) -> int:
    """Count label-layout defects, in display coordinates.

    Defects: a leader line passing through *another* label's box; two leaders
    crossing; a pair of overlapping label boxes; a data point (any plotted feature)
    inside a label box (``points_data``, data coords — the caller passes the labelled
    points; transformed here so a changed axis range is honoured).
    Boxes are shrunk by 1 px so mere contact is not counted. ``text_objs[i]`` and
    ``line_objs[i]`` belong to the same label.
    """
    points_display = ax.transData.transform(points_data)
    boxes: list[tuple[float, float, float, float] | None] = []
    for t in text_objs:
        if t is None:
            boxes.append(None)
            continue
        bb = t.get_window_extent(renderer)  # type: ignore[arg-type]
        boxes.append((bb.x0 + 1.0, bb.y0 + 1.0, bb.x1 - 1.0, bb.y1 - 1.0))

    def inside(pts: np.ndarray, box: tuple[float, float, float, float]) -> np.ndarray:
        bx0, by0, bx1, by1 = box
        return np.asarray(
            (pts[:, 0] >= bx0)
            & (pts[:, 0] <= bx1)
            & (pts[:, 1] >= by0)
            & (pts[:, 1] <= by1)
        )

    conflicts = 0
    seg = np.linspace(0.0, 1.0, 200)[:, None]
    polylines = [
        None if line is None else ax.transData.transform(line.get_xydata())
        for line in line_objs
    ]
    for i, pts in enumerate(polylines):
        if pts is None:
            continue
        samples = np.vstack(
            [pts[m] + seg * (pts[m + 1] - pts[m]) for m in range(len(pts) - 1)]
        )
        for j, box in enumerate(boxes):
            if j != i and box is not None and bool(inside(samples, box).any()):
                conflicts += 1
    # leader-leader crossings
    for i in range(len(polylines)):
        for j in range(i + 1, len(polylines)):
            pi, pj = polylines[i], polylines[j]
            if pi is None or pj is None:
                continue
            for m in range(len(pi) - 1):
                for n in range(len(pj) - 1):
                    if _segments_cross(pi[m], pi[m + 1], pj[n], pj[n + 1]):
                        conflicts += 1
    for i in range(len(boxes)):
        bi = boxes[i]
        if bi is None:
            continue
        conflicts += int(inside(points_display, bi).sum())
        for j in range(i + 1, len(boxes)):
            bj = boxes[j]
            if bj is None:
                continue
            if bi[0] < bj[2] and bj[0] < bi[2] and bi[1] < bj[3] and bj[1] < bi[3]:
                conflicts += 1
    return conflicts


def _segments_cross(
    p1: np.ndarray, p2: np.ndarray, q1: np.ndarray, q2: np.ndarray
) -> bool:
    """Proper intersection of segments p1-p2 and q1-q2 (shared endpoints excluded)."""

    def orient(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
        return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))

    d1, d2 = orient(q1, q2, p1), orient(q1, q2, p2)
    d3, d4 = orient(p1, p2, q1), orient(p1, p2, q2)
    return d1 * d2 < 0 and d3 * d4 < 0


def _column_layout(
    ax: Axes,
    renderer: object,
    anchors_x: list[float],
    anchors_y: list[float],
    label_strs: list[str],
    points_data: np.ndarray,
    *,
    y_top: float,
    y_floor: float,
) -> tuple[list[Text | None], list[Line2D | None]]:
    """Deterministic fallback: elbow-leader label columns, one row per label.

    Labels of points left of x = 0 form a column just left or just right of those
    points; likewise for points right of 0. All labels share ONE top-down row sequence
    (from ``y_top``), so labels of the two sides never share a row. Each point is
    joined to its label by a vertical-then-horizontal leader. Every combination of
    column sides, row order (all permutations for <= 7 labels, else x-sorted orders)
    and start depth is scored analytically — leader-leader crossings, labels covering
    a data point, labels outside the axes, rows below ``y_floor`` — and the best is
    drawn. The column block starts just above the highest labelled point (capped at
    ``y_top``) so leaders stay short. The caller re-verifies the drawn result with
    :func:`_layout_conflicts`.
    """
    base_x0, base_x1 = ax.get_xlim()
    y_lo, y_hi = ax.get_ylim()
    yr = y_hi - y_lo
    gap_y = 0.012 * yr
    ax_w_px = ax.get_window_extent(renderer).width  # type: ignore[arg-type]
    ax_h_px = ax.get_window_extent(renderer).height  # type: ignore[arg-type]
    n = len(label_strs)
    px_sizes: list[tuple[float, float]] = []
    for text in label_strs:
        probe = ax.text(0.0, 0.0, text, fontsize=_LABEL_FONTSIZE)
        bb = probe.get_window_extent(renderer)  # type: ignore[arg-type]
        probe.remove()
        px_sizes.append((bb.width, bb.height))
    heights = [h / ax_h_px * yr for _, h in px_sizes]
    block = sum(heights) + gap_y * (n - 1)
    y_start = min(y_top, max(anchors_y) + 0.03 * yr + block)
    left = [k for k in range(n) if anchors_x[k] < 0.0]
    right = [k for k in range(n) if anchors_x[k] >= 0.0]

    def column(members: list[int], side: str) -> tuple[float, str] | None:
        if not members:
            return None
        if side == "left":
            return min(anchors_x[k] for k in members) - gap_x, "right"
        return max(anchors_x[k] for k in members) + gap_x, "left"

    if n <= 7:
        orders: list[tuple[int, ...]] = list(permutations(range(n)))
    else:
        by_x = sorted(range(n), key=lambda k: anchors_x[k])
        orders = [tuple(by_x), tuple(reversed(by_x))]

    Seg = tuple[np.ndarray, np.ndarray]
    best_score = float("inf")
    best_rows: list[tuple[int, float, float, str, list[Seg]]] = []
    best_lim = (base_x0, base_x1)
    side_opts = [(a, b) for a in ("left", "right") for b in ("left", "right")]
    centre = 0.5 * (base_x0 + base_x1)
    half = 0.5 * (base_x1 - base_x0)
    for factor in _COLUMN_X_EXPANSIONS:
        x0, x1 = centre - half * factor, centre + half * factor
        xr = x1 - x0
        gap_x = 0.02 * xr
        # label widths in data units at this axis range (axes pixel width is fixed)
        sizes = [(w / ax_w_px * xr, h / ax_h_px * yr) for (w, h) in px_sizes]
        for lside, rside in side_opts:
            cols: dict[int, tuple[float, str]] = {}
            for members, side in ((left, lside), (right, rside)):
                col = column(members, side)
                if col is not None:
                    for k in members:
                        cols[k] = col
            for order in orders:
                for m in range(0, 30):
                    top = y_start - m * 0.01 * yr
                    rows: list[tuple[int, float, float, str, list[Seg]]] = []
                    score = 0.0
                    n_covered = 0
                    for k in order:
                        edge, ha = cols[k]
                        w, h = sizes[k]
                        lx0, lx1 = (
                            (edge - w, edge) if ha == "right" else (edge, edge + w)
                        )
                        if lx0 < x0 or lx1 > x1 or top - h < y_floor:
                            score = float("inf")
                            break
                        covered = (
                            (points_data[:, 0] >= lx0)
                            & (points_data[:, 0] <= lx1)
                            & (points_data[:, 1] >= top - h)
                            & (points_data[:, 1] <= top)
                        ).sum()
                        ymid = top - h / 2.0
                        a = np.array([anchors_x[k], anchors_y[k]])
                        elbow = np.array([anchors_x[k], ymid])
                        end = np.array([edge, ymid])
                        # leader through this label's own box region is fine; others are
                        # checked below. Mild tie-breaker: total leader length.
                        score += 100.0 * float(covered)
                        n_covered += int(covered)
                        score += (
                            abs(anchors_y[k] - ymid) / yr
                            + 0.1 * abs(anchors_x[k] - edge) / xr
                        )
                        rows.append((k, edge, top, ha, [(a, elbow), (elbow, end)]))
                        top = top - h - gap_y
                    if score == float("inf"):
                        break
                    # leader-leader crossings and leaders through other labels' boxes
                    for r1 in range(len(rows)):
                        k1, e1, t1, h1, segs1 = rows[r1]
                        w1, hh1 = sizes[k1]
                        bx0, bx1 = (e1 - w1, e1) if h1 == "right" else (e1, e1 + w1)
                        for r2 in range(len(rows)):
                            if r2 == r1:
                                continue
                            segs2 = rows[r2][4]
                            if r2 > r1:
                                for sa in segs1:
                                    for sb in segs2:
                                        if _segments_cross(sa[0], sa[1], sb[0], sb[1]):
                                            score += 100.0
                            for sb in segs2:
                                lo = np.minimum(sb[0], sb[1])
                                hi = np.maximum(sb[0], sb[1])
                                if (
                                    lo[0] < bx1
                                    and hi[0] > bx0
                                    and lo[1] < t1
                                    and hi[1] > t1 - hh1
                                ):
                                    score += 100.0
                    if score < best_score:
                        best_score, best_rows, best_lim = score, rows, (x0, x1)
                    if n_covered == 0:
                        break  # deeper rows only lengthen leaders
                if best_score < 100.0:
                    break
            if best_score < 100.0:
                break
        if best_score < 100.0:
            break

    texts: list[Text | None] = [None] * n
    lines: list[Line2D | None] = [None] * n
    ax.set_xlim(best_lim)
    for k, edge, top, ha, segs in best_rows:
        texts[k] = ax.text(
            edge,
            top,
            label_strs[k],
            ha=ha,
            va="top",
            fontsize=_LABEL_FONTSIZE,
            color="black",
            zorder=5,
        )
        (a, elbow), (_, end) = segs
        (lines[k],) = ax.plot(
            [a[0], elbow[0], end[0]],
            [a[1], elbow[1], end[1]],
            color=_LEADER_COLOR,
            linewidth=_LEADER_WIDTH,
            zorder=3,
        )
    return texts, lines


def _group_coincident(xf: np.ndarray, yf: np.ndarray, within: float) -> list[list[int]]:
    """Single-linkage groups of points closer than ``within`` (axes-fraction coords).

    Each group is sorted left-to-right (then bottom-to-top); groups are ordered by
    their first member's position in the input.
    """
    n = xf.size
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            if float(np.hypot(xf[i] - xf[j], yf[i] - yf[j])) < within:
                parent[find(i)] = find(j)
    groups: dict[int, list[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    out = [
        sorted(g, key=lambda k: (float(xf[k]), float(yf[k]))) for g in groups.values()
    ]
    return sorted(out, key=min)


def _annotate(
    ax: Axes,
    eff: np.ndarray,
    neg_log_q: np.ndarray,
    testable: np.ndarray,
    order: np.ndarray,
    texts: list[str],
    max_distance: float,
    *,
    threshold_y: float,
    cap_below_threshold: bool,
) -> tuple[list[list[int]], float, int]:
    """Ring the labelled points and place their labels collision-free (textalloc).

    Points that coincide at print scale share one stacked label anchored at their
    centroid. With ``cap_below_threshold`` (no labelled feature is a hit) every label
    box is confined below the q-threshold line, and this is verified after drawing
    (fail loud). Several leader lengths are tried (:data:`_LABEL_DISTANCE_FALLBACKS`)
    and the first layout free of leader-through-label crossings and overlapping label
    boxes is kept, falling back to a deterministic elbow-leader column
    layout (:func:`_column_layout`) when textalloc cannot find one. Returns the groups
    (positions into ``order``), the highest label-box top in data units, and the kept
    layout's defect count.
    """
    ax.scatter(
        eff[order],
        neg_log_q[order],
        s=_RING_SIZE,
        facecolors="none",
        edgecolors="black",
        linewidths=0.7,
        zorder=4,
    )
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xs = eff[order]
    ys = neg_log_q[order]
    groups = _group_coincident(
        (xs - x0) / (x1 - x0), (ys - y0) / (y1 - y0), _MERGE_WITHIN
    )
    pad = _LABEL_TOP_PAD_FRAC * (y1 - y0)
    y_top = (
        threshold_y - _LABEL_THRESHOLD_GAP_FRAC * (y1 - y0)
        if cap_below_threshold
        else y1 - pad
    )
    anchors_x = [float(np.mean(xs[g])) for g in groups]
    anchors_y = [float(np.mean(ys[g])) for g in groups]
    label_strs = ["\n".join(texts[k] for k in g) for g in groups]
    renderer = ax.figure.canvas.get_renderer()  # type: ignore[attr-defined]

    finite = testable & np.isfinite(neg_log_q)
    # Final verification counts only the ringed (labelled) points as "covered": a label
    # over a gray background point is tolerable, over a labelled point it is not.
    ringed_data = np.column_stack([xs, ys])

    def run_textalloc(dist: float) -> tuple[list[Text | None], list[Line2D | None]]:
        _, _, raw_texts, raw_lines = ta.allocate(
            ax,
            anchors_x,
            anchors_y,
            label_strs,
            x_scatter=eff[testable],
            y_scatter=neg_log_q[testable],
            textsize=_LABEL_FONTSIZE,
            draw_lines=True,
            linecolor=_LEADER_COLOR,
            linewidth=_LEADER_WIDTH,
            textcolor="black",
            margin=_LABEL_MARGIN,
            min_distance=_LABEL_MIN_DISTANCE,
            max_distance=dist,
            nbr_candidates=_LABEL_CANDIDATES,
            avoid_label_lines_overlap=True,
            avoid_crossing_label_lines=True,
            xlims=(x0, x1),
            ylims=(y0 + pad, y_top),
        )
        return cast(list[Text | None], raw_texts), cast(list[Line2D | None], raw_lines)

    points_data = np.column_stack([eff[finite], neg_log_q[finite]])

    def run_column() -> tuple[list[Text | None], list[Line2D | None]]:
        return _column_layout(
            ax,
            renderer,
            anchors_x,
            anchors_y,
            label_strs,
            points_data,
            y_top=y_top,
            y_floor=y0 + pad,
        )

    # Candidate layouts in preference order: textalloc at the requested leader length,
    # then the fallback lengths, then the deterministic elbow-leader column layout.
    # The first defect-free one (see _layout_conflicts; a label box above the cap also
    # counts) is kept; otherwise the least-defective, whose count is returned.
    tried = [max_distance, *(d for d in _LABEL_DISTANCE_FALLBACKS if d != max_distance)]
    Layout = tuple[list[Text | None], list[Line2D | None]]
    candidates: Sequence[Callable[[], Layout]] = [
        partial(run_textalloc, d) for d in tried
    ] + [run_column]
    to_data_sel = ax.transData.inverted()
    best: tuple[int, list[Text], list[Line2D]] | None = None
    xlim_orig = ax.get_xlim()
    best_xlim = xlim_orig
    for make in candidates:
        ax.set_xlim(xlim_orig)
        text_objs, line_objs = make()
        t_objs: list[Text] = [t for t in text_objs if t is not None]
        l_objs: list[Line2D] = [ln for ln in line_objs if ln is not None]
        if len(t_objs) != len(groups):
            raise AssertionError(
                f"layout drew {len(t_objs)} labels for {len(groups)} label groups."
            )
        n_conf = _layout_conflicts(ax, renderer, text_objs, line_objs, ringed_data)
        top_now = max(
            float(to_data_sel.transform(t.get_window_extent(renderer)).max(axis=0)[1])
            for t in t_objs
        )
        if cap_below_threshold and top_now >= threshold_y:
            n_conf += 1000
        if best is None or n_conf < best[0]:
            if best is not None:
                for artist in (*best[1], *best[2]):
                    artist.remove()
            best = (n_conf, t_objs, l_objs)
            best_xlim = ax.get_xlim()
        else:
            for artist in (*t_objs, *l_objs):
                artist.remove()
        if n_conf == 0:
            break
    if best is None:
        raise AssertionError("no label layout was produced.")
    n_conflicts, new_texts, _ = best
    ax.set_xlim(best_xlim)

    to_data = ax.transData.inverted()
    tops = [
        float(to_data.transform(t.get_window_extent(renderer)).max(axis=0)[1])
        for t in new_texts
    ]
    box_top = max(tops)
    if cap_below_threshold and box_top >= threshold_y:
        raise AssertionError(
            f"a non-significant label box reaches y={box_top:.3f} >= the q-threshold "
            f"line at {threshold_y:.3f}."
        )
    return groups, box_top, n_conflicts


def _legend_figure(
    color_map: dict[str, str],
    counts: VolcanoCounts,
    direction_labels: tuple[str, str],
    legend_title: str,
    *,
    fdr: float,
    n_labelled: int,
    n_labelled_hits: int,
) -> Figure:
    """Standalone legend: up / down / NS swatches with the hit counts (+ label key)."""
    up_label, down_label = direction_labels
    entries = [
        (_UP_KEY, f"q < {fdr:g}, {up_label} (n={counts.up:,})"),
        (_DOWN_KEY, f"q < {fdr:g}, {down_label} (n={counts.down:,})"),
        (_NS_KEY, f"not significant (n={counts.ns:,})"),
    ]
    handles: list[Line2D] = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=color_map[key],
            markeredgecolor="white",
            markersize=9,
        )
        for key, _ in entries
    ]
    texts = [text for _, text in entries]
    handles.append(
        Line2D([0], [0], color="0.3", linestyle="--", linewidth=1.2),
    )
    texts.append(f"BH q = {fdr:g} threshold")
    if n_labelled:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor="none",
                markeredgecolor="black",
                markersize=8,
            )
        )
        texts.append(
            f"labelled: {n_labelled} smallest q ({n_labelled_hits} of {n_labelled} "
            f"significant)"
        )
    height = 0.3 * len(texts) + 0.45
    fig, ax = plt.subplots(figsize=(5.2, height))
    ax.axis("off")
    ax.legend(
        handles,
        texts,
        title=legend_title,
        loc="center",
        frameon=True,
        fontsize=11,
        title_fontsize=12,
    )
    return fig
