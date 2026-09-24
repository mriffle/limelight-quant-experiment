"""Plot primitives for the limma-trend sensitivity (0006) and ATPK case (0007) figures.

From-scratch project module (no ``lib/`` template covers these plot types); it reuses
the project's shared figure machinery: colors come from the registry via the caller
(``common.figures.colors``), and label placement for the q-vs-q scatter reuses the
volcano module's collision-free annotator (``analysis_figures.volcano._annotate``:
coincident points share one stacked label, textalloc placement with a deterministic
fallback) so labelled hits look the same as on the volcanoes.

Plot primitives (all decoupled from file I/O — plain arrays / small dataclasses in,
``RenderedFigure`` out; the caller saves via ``common.figures.figure_io.save_figure``
inside ``publication_style``):

* :func:`plot_prior_sd_trend` — per-feature residual SD vs mean log2 abundance with
  the constant no-trend prior SD and the limma-trend prior-SD curve overlaid.
* :func:`plot_q_trend_vs_notrend` — ``-log10 q`` (limma-trend) vs ``-log10 q``
  (no-trend) with ``y = x``, q-threshold guides on both axes, points colored by mean
  log2 abundance (sequential; the colorbar is the separate legend image).
* :func:`plot_relabel_diagnostic` — per quantity, the hit count and Storey pi0 of all
  within-pair labellings as dot strips, the observed labelling highlighted.
* :func:`plot_paired_panels` — small multiples of per-pair control -> treated
  segments (one line per pair, colored by pair), on log2-scaled axes that share one
  log2 span so slopes are comparable across panels.

Fail loud on shape / range mismatches; no silent NaN coercion (a NaN in a paired
panel is drawn as a missing point and must be declared by the caller via
``allow_missing``).
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from itertools import permutations

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from matplotlib import colors as mcolors
from matplotlib import patheffects
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

from analysis_figures import volcano as vo

__script_meta__: dict[str, object] = {
    "kind": "module",
    "provides": [
        "RenderedFigure",
        "RelabelQuantity",
        "PairedPanel",
        "plot_prior_sd_trend",
        "plot_q_trend_vs_notrend",
        "plot_relabel_diagnostic",
        "plot_paired_panels",
        "equal_span_limits",
    ],
    "uses": ["analysis_figures.volcano", "common.figures.figure_io"],
    "seeded_from": None,
    "description": (
        "Figure primitives for finding 0006 (limma-trend sensitivity: prior-SD "
        "trend, q trend-vs-no-trend, pair-relabel diagnostic) and 0007 (ATPK "
        "per-pair LFQ vs spectral segments). Colors are supplied by the caller from "
        "the registry; the colorbar / swatch keys are separate legend figures."
    ),
}

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]

MINUS = chr(0x2212)  # typographic minus sign
POINT_GRAY = "#505050"  # neutral background points (not a category)
REFERENCE_GRAY = "0.45"  # y = x / guide lines
GUIDE_STYLES: tuple[str, ...] = ("--", ":")  # q = 0.05, q = 0.10
SEQUENTIAL_CMAP = "viridis"
COUNT_TICKS: tuple[float, ...] = (1, 2, 3, 4, 6, 8, 10, 12, 16, 24, 32, 48, 64)


def fmt_signed(value: float, digits: int = 2) -> str:
    """Signed fixed-point text with a typographic minus: -0.516 -> (minus)0.52."""
    text = f"{value:+.{digits}f}"
    return text.replace("-", MINUS)


def fmt_q(q: float) -> str:
    """BH q for on-canvas notes: 3 decimals below 0.1, else 2 (0.039 / 0.70)."""
    if not 0.0 <= q <= 1.0:
        raise ValueError(f"q must be in [0, 1]; got {q}.")
    return f"{q:.3f}" if q < 0.1 else f"{q:.2f}"


@dataclass
class RenderedFigure:
    """A figure, its separate legend figure, and the numbers drawn (for provenance)."""

    figure: Figure
    legend_figure: Figure
    stats: dict[str, object] = field(default_factory=dict)


def _as_float(name: str, values: object) -> FloatArray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1-D; got shape {arr.shape}.")
    return arr


def _legend_figure(
    handles: Sequence[Line2D], labels: Sequence[str], title: str | None = None
) -> Figure:
    if len(handles) != len(labels) or not handles:
        raise ValueError("legend needs one label per handle (and >= 1 handle).")
    width = max(3.0, 0.085 * max(len(s) for s in labels) + 1.0)
    height = 0.32 * len(handles) + (0.55 if title else 0.25)
    fig = plt.figure(figsize=(width, height))
    fig.legend(
        handles,
        labels,
        loc="center",
        frameon=False,
        handlelength=2.4,
        title=title,
        title_fontsize=11,
        fontsize=10,
    )
    return fig


def _suptitle(fig: Figure, title: str, y: float = 0.995) -> None:
    fig.suptitle(title, fontsize=12, fontweight="bold", y=y, va="top")


# --------------------------------------------------------------------------- #
# 1. Residual SD vs abundance with the two priors
# --------------------------------------------------------------------------- #
def plot_prior_sd_trend(
    mean_abundance: object,
    residual_sd: object,
    trend_prior_sd: object,
    *,
    notrend_prior_sd: float,
    series_colors: Mapping[str, str],
    highlight: Mapping[str, int] | None = None,
    residual_df: int,
    d0_trend: float,
    d0_notrend: float,
    feature_noun: str,
    abundance_label: str,
    title: str,
) -> RenderedFigure:
    """Per-feature residual SD (log axis) vs mean abundance + both priors.

    ``series_colors`` must hold ``"no-trend"`` and ``"limma-trend"`` (registry
    colors). ``highlight`` maps a display label to a feature position to ring + label.
    """
    x = _as_float("mean_abundance", mean_abundance)
    y = _as_float("residual_sd", residual_sd)
    prior = _as_float("trend_prior_sd", trend_prior_sd)
    if not (x.shape == y.shape == prior.shape):
        raise ValueError("mean_abundance, residual_sd, trend_prior_sd lengths differ.")
    if not (np.isfinite(x).all() and np.isfinite(y).all() and np.isfinite(prior).all()):
        raise ValueError("non-finite abundance / residual SD / prior SD.")
    if (y <= 0).any():
        raise ValueError("residual SD must be > 0 for a log axis.")
    if not notrend_prior_sd > 0:
        raise ValueError(f"notrend_prior_sd must be > 0; got {notrend_prior_sd}.")
    for key in ("no-trend", "limma-trend"):
        if key not in series_colors:
            raise KeyError(f"series_colors lacks {key!r}.")
    order = np.argsort(x, kind="stable")
    curve_x = x[order]
    curve_y = prior[order]
    # The trend prior is a smooth function of abundance: sorted by x it must be a
    # function (equal x -> equal prior) — check rather than assume.
    if np.any(np.diff(curve_x) == 0) and np.any(
        np.abs(np.diff(curve_y))[np.diff(curve_x) == 0] > 1e-12
    ):
        raise ValueError("trend prior SD is not a function of mean abundance.")

    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    fig.subplots_adjust(left=0.12, right=0.88, bottom=0.12, top=0.86)
    ax.scatter(
        x, y, s=7, color=POINT_GRAY, alpha=0.4, linewidths=0, zorder=2, rasterized=False
    )
    ax.axhline(
        notrend_prior_sd,
        color=series_colors["no-trend"],
        lw=2.2,
        ls="--",
        zorder=4,
    )
    ax.plot(curve_x, curve_y, color=series_colors["limma-trend"], lw=2.6, zorder=5)
    ax.set_yscale("log")
    lo = 10 ** (math.floor(math.log10(float(y.min())) * 4) / 4)
    hi = 10 ** (math.ceil(math.log10(float(y.max())) * 4) / 4)
    ax.set_ylim(lo, hi)
    ticks = [t for t in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0) if lo <= t <= hi]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{v:g}"))
    ax.set_xlabel(abundance_label)
    ax.set_ylabel(f"residual SD (log2 units, df = {residual_df}; log axis)")
    ax.grid(True, which="major", alpha=0.25)

    # Right-margin value labels for the two priors (outside the data area).
    ax.text(
        1.01,
        notrend_prior_sd,
        f"{notrend_prior_sd:.3f}",
        transform=ax.get_yaxis_transform(),
        ha="left",
        va="center",
        fontsize=9.5,
        color=series_colors["no-trend"],
        fontweight="bold",
    )
    ax.text(
        1.01,
        float(curve_y[-1]),
        f"{float(curve_y[-1]):.3f}",
        transform=ax.get_yaxis_transform(),
        ha="left",
        va="center",
        fontsize=9.5,
        color=series_colors["limma-trend"],
        fontweight="bold",
    )
    ax.text(
        -0.01,
        float(curve_y[0]),
        f"{float(curve_y[0]):.3f}",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=9.5,
        color=series_colors["limma-trend"],
        fontweight="bold",
    )

    highlighted: dict[str, dict[str, float]] = {}
    for label, idx in (highlight or {}).items():
        if not 0 <= idx < x.size:
            raise IndexError(f"highlight {label!r}: index {idx} out of range.")
        ax.scatter(
            [x[idx]],
            [y[idx]],
            s=46,
            facecolors="none",
            edgecolors="black",
            linewidths=1.1,
            zorder=6,
        )
        ax.annotate(
            label,
            xy=(float(x[idx]), float(y[idx])),
            xytext=(22, -14),
            textcoords="offset points",
            fontsize=10,
            fontweight="bold",
            ha="left",
            va="center",
            arrowprops={"arrowstyle": "-", "color": "0.3", "lw": 0.8},
            zorder=7,
        )
        highlighted[label] = {
            "mean_log2_abundance": float(x[idx]),
            "residual_sd": float(y[idx]),
            "trend_prior_sd": float(prior[idx]),
        }
    _suptitle(fig, title)

    handles = [
        Line2D([], [], ls="none", marker="o", ms=5, color=POINT_GRAY, alpha=0.6),
        Line2D([], [], color=series_colors["no-trend"], lw=2.2, ls="--"),
        Line2D([], [], color=series_colors["limma-trend"], lw=2.6),
    ]
    labels = [
        f"per-{feature_noun} residual SD (n = {x.size:,})",
        f"no-trend prior SD = {notrend_prior_sd:.3f} (d0 = {d0_notrend:.1f})",
        (
            f"limma-trend prior SD {float(curve_y[0]):.3f} → "
            f"{float(curve_y[-1]):.3f} (d0 = {d0_trend:.1f})"
        ),
    ]
    if highlight:
        handles.append(
            Line2D(
                [],
                [],
                ls="none",
                marker="o",
                ms=8,
                mfc="none",
                mec="black",
                mew=1.1,
            )
        )
        labels.append("labelled: limma-trend hit (q < 0.05)")
    legend = _legend_figure(handles, labels, title="Key")
    stats: dict[str, object] = {
        "n_features": int(x.size),
        "notrend_prior_sd": notrend_prior_sd,
        "trend_prior_sd_at_min_abundance": float(curve_y[0]),
        "trend_prior_sd_at_max_abundance": float(curve_y[-1]),
        "trend_prior_sd_range": [float(prior.min()), float(prior.max())],
        "residual_sd_range": [float(y.min()), float(y.max())],
        "y_limits": [lo, hi],
        "highlighted": highlighted,
    }
    return RenderedFigure(fig, legend, stats)


# --------------------------------------------------------------------------- #
# 2. q (trend) vs q (no-trend)
# --------------------------------------------------------------------------- #
def _neg_log10(q: FloatArray) -> FloatArray:
    floor = float(np.finfo(float).tiny)
    return np.asarray(-np.log10(np.maximum(q, floor)), dtype=np.float64)


_QQ_LABEL_FONTSIZE = 8.5
_QQ_MERGE_WITHIN = 0.012  # axes fraction; as the volcano's coincident-label merge


def _place_qq_labels(
    ax: Axes,
    fig: Figure,
    x: FloatArray,
    y: FloatArray,
    idx: IntArray,
    labels: list[str],
) -> tuple[list[list[int]], int]:
    """Ring labelled points; deterministic labels: stacks above, singletons left.

    Points closer than :data:`_QQ_MERGE_WITHIN` (axes fraction; the volcano's
    single-linkage grouping) share one stacked label placed above the group (the
    axis carries headroom for it); a lone point gets its label to its left at the same
    height. Leaders join each label to its point / group centroid. After drawing,
    label boxes must not overlap each other or any ringed point (fail loud); the
    number of background points under a label box is returned as the conflict count.
    """
    xs, ys = x[idx], y[idx]
    ax.scatter(
        xs, ys, s=22, facecolors="none", edgecolors="black", linewidths=0.9, zorder=4
    )
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    groups = vo._group_coincident(  # reuse the volcano's coincident-label grouping
        (xs - x0) / (x1 - x0), (ys - y0) / (y1 - y0), _QQ_MERGE_WITHIN
    )
    dx = 0.03 * (x1 - x0)
    texts = []
    for g in groups:
        gx, gy = float(np.mean(xs[g])), float(np.mean(ys[g]))
        if len(g) > 1:
            text = "\n".join(labels[k] for k in g)
            t = ax.annotate(
                text,
                xy=(gx, gy),
                xytext=(0, 14),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=_QQ_LABEL_FONTSIZE,
                bbox={"fc": "white", "ec": "none", "alpha": 0.85, "pad": 1.0},
                arrowprops={"arrowstyle": "-", "color": "0.35", "lw": 0.7},
                zorder=6,
            )
        else:
            t = ax.annotate(
                labels[g[0]],
                xy=(gx, gy),
                xytext=(gx - dx, gy),
                textcoords="data",
                ha="right",
                va="center",
                fontsize=_QQ_LABEL_FONTSIZE,
                bbox={"fc": "white", "ec": "none", "alpha": 0.85, "pad": 1.0},
                arrowprops={"arrowstyle": "-", "color": "0.35", "lw": 0.7},
                zorder=6,
            )
        texts.append(t)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()  # type: ignore[attr-defined]
    # Text-only extents (Annotation.get_window_extent also unions the leader arrow).
    boxes = [Text.get_window_extent(t, renderer) for t in texts]
    axbox = ax.get_window_extent(renderer)
    for i, b in enumerate(boxes):
        if b.x0 < axbox.x0 or b.x1 > axbox.x1 or b.y0 < axbox.y0 or b.y1 > axbox.y1:
            raise AssertionError(f"label {i} extends outside the axes.")
        for j in range(i + 1, len(boxes)):
            if b.overlaps(boxes[j]):
                raise AssertionError(f"labels {i} and {j} overlap.")
    ringed = ax.transData.transform(np.column_stack([xs, ys]))
    allpts = ax.transData.transform(np.column_stack([x, y]))
    conflicts = 0
    for b in boxes:
        pad = 2.0
        if any(
            b.x0 - pad <= px <= b.x1 + pad and b.y0 - pad <= py <= b.y1 + pad
            for px, py in ringed
        ):
            raise AssertionError(
                f"label box {b} covers a labelled point; ringed {ringed.tolist()}."
            )
        conflicts += int(
            np.sum(
                (allpts[:, 0] >= b.x0)
                & (allpts[:, 0] <= b.x1)
                & (allpts[:, 1] >= b.y0)
                & (allpts[:, 1] <= b.y1)
            )
        )
    return groups, conflicts


def plot_q_trend_vs_notrend(
    q_trend: object,
    q_notrend: object,
    mean_abundance: object,
    *,
    label_index: Sequence[int],
    labels: Sequence[str],
    fdr_lines: tuple[float, float] = (0.05, 0.10),
    abundance_label: str,
    feature_noun: str,
    title: str,
) -> RenderedFigure:
    """``-log10 q`` trend (y) vs no-trend (x), colored by mean abundance (viridis)."""
    qt = _as_float("q_trend", q_trend)
    qn = _as_float("q_notrend", q_notrend)
    ab = _as_float("mean_abundance", mean_abundance)
    if not (qt.shape == qn.shape == ab.shape):
        raise ValueError("q_trend, q_notrend, mean_abundance lengths differ.")
    for name, arr in (("q_trend", qt), ("q_notrend", qn)):
        if not np.isfinite(arr).all() or arr.min() < 0 or arr.max() > 1:
            raise ValueError(f"{name} must be finite and in [0, 1].")
    if not np.isfinite(ab).all():
        raise ValueError("mean_abundance has non-finite values.")
    if len(label_index) != len(labels):
        raise ValueError("label_index and labels lengths differ.")
    idx = np.asarray(label_index, dtype=np.int64)
    if idx.size and (idx.min() < 0 or idx.max() >= qt.size):
        raise IndexError("label_index out of range.")
    lo_fdr, hi_fdr = sorted(fdr_lines)
    if not 0 < lo_fdr < hi_fdr < 1:
        raise ValueError(f"fdr_lines must be two values in (0, 1); got {fdr_lines}.")

    x = _neg_log10(qn)
    y = _neg_log10(qt)
    top = max(float(x.max()), float(y.max()), -math.log10(lo_fdr))
    lim = math.ceil((top * 1.3) * 10) / 10  # headroom for stacked labels
    norm = mcolors.Normalize(vmin=float(ab.min()), vmax=float(ab.max()))

    fig, ax = plt.subplots(figsize=(7.0, 7.0))
    fig.subplots_adjust(left=0.12, right=0.86, bottom=0.1, top=0.86)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.plot([0, lim], [0, lim], color=REFERENCE_GRAY, lw=1.0, zorder=1)
    for fdr, ls in zip((lo_fdr, hi_fdr), GUIDE_STYLES, strict=True):
        v = -math.log10(fdr)
        ax.axhline(v, color="0.3", ls=ls, lw=1.0, zorder=1)
        ax.axvline(v, color="0.3", ls=ls, lw=1.0, zorder=1)
        ax.text(
            1.01,
            v,
            f"q = {fdr:.2f}",
            transform=ax.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=9,
            color="0.3",
        )
        ax.text(
            v - 0.008 * lim,
            0.02,
            f"q = {fdr:.2f}",
            transform=ax.get_xaxis_transform(),
            ha="right",
            va="bottom",
            fontsize=9,
            color="0.3",
            rotation=90,
        )
    order = np.argsort(ab, kind="stable")
    ax.scatter(
        x[order],
        y[order],
        c=ab[order],
        cmap=SEQUENTIAL_CMAP,
        norm=norm,
        s=12,
        alpha=0.85,
        linewidths=0,
        zorder=2,
    )
    ax.set_xlabel(r"$-\log_{10}$(BH $q$), no-trend prior")
    ax.set_ylabel(r"$-\log_{10}$(BH $q$), limma-trend prior")
    ax.grid(True, alpha=0.2)

    groups: list[list[int]] = []
    conflicts = 0
    if idx.size:
        groups, conflicts = _place_qq_labels(ax, fig, x, y, idx, list(labels))
    _suptitle(fig, title)

    # Legend image: the abundance colorbar + the line/ring key.
    n_lines = 3 + (1 if idx.size else 0)
    legend = plt.figure(figsize=(4.6, 1.1 + 0.3 * n_lines))
    cax = legend.add_axes((0.08, 0.78, 0.84, 0.08))
    cbar = legend.colorbar(
        ScalarMappable(norm=norm, cmap=SEQUENTIAL_CMAP),
        cax=cax,
        orientation="horizontal",
    )
    cbar.set_label(abundance_label, fontsize=10)
    handles = [
        Line2D([], [], color=REFERENCE_GRAY, lw=1.0),
        Line2D([], [], color="0.3", lw=1.0, ls=GUIDE_STYLES[0]),
        Line2D([], [], color="0.3", lw=1.0, ls=GUIDE_STYLES[1]),
    ]
    texts = [
        "y = x (models agree)",
        f"BH q = {lo_fdr:.2f}",
        f"BH q = {hi_fdr:.2f}",
    ]
    if idx.size:
        handles.append(
            Line2D([], [], ls="none", marker="o", ms=7, mfc="none", mec="black")
        )
        texts.append(f"labelled: limma-trend hits, q < {lo_fdr:g} (n = {idx.size})")
    legend.legend(
        handles,
        texts,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        frameon=False,
        fontsize=10,
        title=f"Each point = one {feature_noun}",
        title_fontsize=10,
    )
    stats: dict[str, object] = {
        "n_features": int(qt.size),
        "axis_limit_neg_log10": lim,
        "n_trend_q_lt": {f"{f:g}": int((qt < f).sum()) for f in (lo_fdr, hi_fdr)},
        "n_notrend_q_lt": {f"{f:g}": int((qn < f).sum()) for f in (lo_fdr, hi_fdr)},
        "n_above_diagonal": int((y > x + 1e-12).sum()),
        "n_below_diagonal": int((y < x - 1e-12).sum()),
        "abundance_color_range": [float(ab.min()), float(ab.max())],
        "label_groups": [[str(labels[k]) for k in g] for g in groups],
        "label_layout_conflicts": conflicts,
    }
    return RenderedFigure(fig, legend, stats)


# --------------------------------------------------------------------------- #
# 3. Pair-relabel diagnostic
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class RelabelQuantity:
    """One quantity's labellings (observed included) for the strip plot."""

    name: str
    hits: Sequence[int]
    pi0: Sequence[float]
    observed: int  # position of the observed labelling within hits / pi0
    n_hits_ge_observed: int
    pi0_rank_lowest_first: int

    def __post_init__(self) -> None:
        if len(self.hits) != len(self.pi0) or not self.hits:
            raise ValueError(f"{self.name}: hits / pi0 lengths differ or are empty.")
        if not 0 <= self.observed < len(self.hits):
            raise IndexError(f"{self.name}: observed index out of range.")
        if any(h < 0 for h in self.hits):
            raise ValueError(f"{self.name}: negative hit count.")
        if any(not 0 <= p <= 1 for p in self.pi0):
            raise ValueError(f"{self.name}: pi0 outside [0, 1].")


def _dodge(values: FloatArray, tol: float, step: float) -> FloatArray:
    """Horizontal offsets that spread (near-)tied values symmetrically about 0."""
    offsets = np.zeros(values.size, dtype=np.float64)
    order = np.argsort(values, kind="stable")
    clusters: list[list[int]] = []
    for i in order:
        if clusters and abs(values[i] - values[clusters[-1][-1]]) <= tol:
            clusters[-1].append(int(i))
        else:
            clusters.append([int(i)])
    for cluster in clusters:
        members = sorted(cluster)
        m = len(members)
        for k, i in enumerate(members):
            offsets[i] = (k - (m - 1) / 2) * step
    return offsets


def plot_relabel_diagnostic(
    quantities: Sequence[RelabelQuantity],
    *,
    observed_color: str,
    other_color: str,
    fdr: float,
    title: str,
) -> RenderedFigure:
    """Two panels (hit count, pi0); x = quantity; 8 dots each, observed highlighted."""
    if not quantities:
        raise ValueError("no quantities to plot.")
    n_lab = {len(q.hits) for q in quantities}
    if len(n_lab) != 1:
        raise ValueError(f"quantities have differing labelling counts {n_lab}.")
    n_labellings = n_lab.pop()

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.9))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.22, top=0.8, wspace=0.22)
    ax_h, ax_p = axes
    positions = np.arange(len(quantities), dtype=np.float64)
    max_hits = max(max(q.hits) for q in quantities)
    step = 0.1
    for pos, q in zip(positions, quantities, strict=True):
        for ax, vals, tol in (
            (ax_h, np.asarray(q.hits, dtype=np.float64), 0.0),
            (ax_p, np.asarray(q.pi0, dtype=np.float64), 0.015),
        ):
            off = _dodge(vals, tol, step)
            other = np.ones(vals.size, dtype=bool)
            other[q.observed] = False
            ax.scatter(
                pos + off[other],
                vals[other],
                s=34,
                color=other_color,
                edgecolors="0.35",
                linewidths=0.7,
                zorder=3,
            )
            ax.scatter(
                [pos + off[q.observed]],
                [vals[q.observed]],
                s=95,
                marker="D",
                color=observed_color,
                edgecolors="black",
                linewidths=0.9,
                zorder=4,
            )
        ax_h.text(
            pos,
            -0.17,
            f"perm p = {q.n_hits_ge_observed}/{n_labellings}",
            transform=ax_h.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=9.5,
        )
        ax_p.text(
            pos,
            -0.17,
            f"obs rank {q.pi0_rank_lowest_first}/{n_labellings}",
            transform=ax_p.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=9.5,
        )
    for ax in axes:
        ax.set_xticks(positions)
        ax.set_xticklabels([q.name for q in quantities])
        ax.set_xlim(-0.5, len(quantities) - 0.5)
        ax.grid(True, axis="y", alpha=0.25)
    ax_h.set_ylim(-0.6, max(max_hits, 1) * 1.12 + 0.5)
    ax_h.yaxis.set_major_locator(
        FixedLocator(list(range(0, int(max_hits) + 2, 2 if max_hits > 6 else 1)))
    )
    ax_h.set_ylabel(f"hits at BH q < {fdr:g} (count)")
    ax_h.set_title("Hit count per labelling", fontsize=12)
    pmin = min(min(q.pi0) for q in quantities)
    ax_p.set_ylim(math.floor(pmin * 20) / 20 - 0.02, 1.02)
    ax_p.set_ylabel(r"Storey $\pi_0$ ($\lambda$ = 0.5)")
    ax_p.set_title(r"$\pi_0$ per labelling (lower = more signal)", fontsize=12)
    _suptitle(fig, title)

    handles = [
        Line2D(
            [],
            [],
            ls="none",
            marker="D",
            ms=9,
            mfc=observed_color,
            mec="black",
        ),
        Line2D([], [], ls="none", marker="o", ms=7, mfc=other_color, mec="0.35"),
    ]
    labels = [
        "observed labelling",
        f"within-pair relabellings ({n_labellings - 1})",
    ]
    legend = _legend_figure(handles, labels, title="Labelling")
    stats: dict[str, object] = {
        q.name: {
            "hits": list(q.hits),
            "pi0": list(q.pi0),
            "observed": q.observed,
            "n_hits_ge_observed": q.n_hits_ge_observed,
            "pi0_rank_lowest_first": q.pi0_rank_lowest_first,
        }
        for q in quantities
    }
    return RenderedFigure(fig, legend, stats)


# --------------------------------------------------------------------------- #
# 4. Per-pair control -> treated segments
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PairedPanel:
    """One panel of per-pair control -> treated values.

    ``control`` / ``treated`` are aligned with the figure's ``pairs``. With
    ``count_axis`` the values are raw counts drawn on a log2 axis (ticks in counts);
    otherwise they are already log2 values on a linear axis. NaN = not quantified,
    allowed only when ``allow_missing``.
    """

    title: str
    ylabel: str
    control: Sequence[float]
    treated: Sequence[float]
    count_axis: bool = False
    note: str | None = None
    allow_missing: bool = False
    missing_note: str | None = None

    def log2_values(self) -> tuple[FloatArray, FloatArray]:
        """Control / treated on the log2 scale the axis displays."""
        c = np.asarray(self.control, dtype=np.float64)
        t = np.asarray(self.treated, dtype=np.float64)
        if c.shape != t.shape or c.ndim != 1:
            raise ValueError(f"{self.title}: control/treated shapes differ.")
        if not self.allow_missing and not (
            np.isfinite(c).all() and np.isfinite(t).all()
        ):
            raise ValueError(f"{self.title}: missing values but allow_missing=False.")
        if self.count_axis:
            finite = np.concatenate([c[np.isfinite(c)], t[np.isfinite(t)]])
            if (finite <= 0).any():
                raise ValueError(f"{self.title}: counts must be > 0 on a log2 axis.")
            with np.errstate(invalid="ignore"):
                return np.log2(c), np.log2(t)
        return c, t


def equal_span_limits(
    groups: Sequence[Sequence[PairedPanel]], pad: float = 0.25
) -> list[tuple[float, float]]:
    """Per-group log2 y-limits (centered on each group's data) sharing ONE span.

    Every group gets the same span (the widest group's data range + ``2 * pad``), so
    a within-pair slope reads identically in every panel. Returned in log2 units.
    """
    ranges: list[tuple[float, float]] = []
    for group in groups:
        vals = np.concatenate([np.concatenate(p.log2_values()) for p in group])
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            raise ValueError("a panel group has no finite values.")
        ranges.append((float(vals.min()), float(vals.max())))
    span = max(hi - lo for lo, hi in ranges) + 2 * pad
    out: list[tuple[float, float]] = []
    for lo, hi in ranges:
        mid = (lo + hi) / 2
        out.append((mid - span / 2, mid + span / 2))
    return out


def plot_paired_panels(
    grid: Sequence[Sequence[PairedPanel]],
    ylims_log2: Sequence[Sequence[tuple[float, float]]],
    *,
    pairs: Sequence[str],
    pair_colors: Mapping[str, str],
    pair_legend_labels: Mapping[str, str],
    condition_labels: tuple[str, str],
    title: str,
    panel_size: tuple[float, float] = (3.5, 3.9),
    dodge: float = 0.05,
) -> RenderedFigure:
    """Rows x columns of per-pair segments (one colored line per pair), dodged in x."""
    n_rows = len(grid)
    if n_rows == 0 or len({len(r) for r in grid}) != 1:
        raise ValueError("grid must be a non-empty rectangular list of rows.")
    n_cols = len(grid[0])
    if len(ylims_log2) != n_rows or any(len(r) != n_cols for r in ylims_log2):
        raise ValueError("ylims_log2 must match the grid shape.")
    missing = [p for p in pairs if p not in pair_colors or p not in pair_legend_labels]
    if missing:
        raise KeyError(f"pairs without a color / legend label: {missing}.")
    width = panel_size[0] * n_cols + 0.9
    height = panel_size[1] * n_rows + 1.25
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(width, height), squeeze=False)
    fig.subplots_adjust(
        left=0.9 / width + 0.04,
        right=0.985,
        bottom=0.55 / height,
        top=1 - 1.25 / height,
        wspace=0.42,
        hspace=0.62 if n_rows > 1 else 0.3,
    )
    n = len(pairs)
    offsets = [(i - (n - 1) / 2) * dodge for i in range(n)]
    drawn: dict[str, object] = {}
    for r, row in enumerate(grid):
        for c, panel in enumerate(row):
            ax: Axes = axes[r][c]
            lc, lt = panel.log2_values()
            if lc.size != n:
                raise ValueError(f"{panel.title}: {lc.size} values for {n} pairs.")
            for i, pair in enumerate(pairs):
                color = pair_colors[pair]
                xs = np.array([0.0, 1.0]) + offsets[i]
                ys = np.array([lc[i], lt[i]])
                ok = np.isfinite(ys)
                stroke = [
                    patheffects.Stroke(linewidth=3.6, foreground="0.25"),
                    patheffects.Normal(),
                ]
                if ok.all():
                    ax.plot(
                        xs,
                        ys,
                        color=color,
                        lw=2.2,
                        zorder=3,
                        path_effects=stroke,
                        solid_capstyle="round",
                    )
                ax.scatter(
                    xs[ok],
                    ys[ok],
                    s=40,
                    color=color,
                    edgecolors="0.15",
                    linewidths=0.9,
                    zorder=4,
                )
            lo, hi = ylims_log2[r][c]
            ax.set_ylim(lo, hi)
            if panel.count_axis:
                ticks = [math.log2(t) for t in COUNT_TICKS if lo <= math.log2(t) <= hi]
                ax.yaxis.set_major_locator(FixedLocator(ticks))
                ax.yaxis.set_major_formatter(
                    FuncFormatter(lambda v, _pos: f"{2**v:.0f}")
                )
            ax.set_xlim(-0.45, 1.45)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(list(condition_labels))
            ax.set_ylabel(panel.ylabel)
            ax.grid(True, axis="y", alpha=0.25)
            ax.set_title(panel.title, fontsize=12, fontweight="bold", pad=24)
            if panel.note:
                ax.text(
                    0.5,
                    1.02,
                    panel.note,
                    transform=ax.transAxes,
                    ha="center",
                    va="bottom",
                    fontsize=9.5,
                )
            if panel.missing_note:
                ax.text(
                    0.5,
                    0.04,
                    panel.missing_note,
                    transform=ax.transAxes,
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    color="0.25",
                )
            drawn[f"r{r}c{c}:{panel.title}"] = {
                "ylim_log2": [lo, hi],
                "control_log2": [None if not math.isfinite(v) else v for v in lc],
                "treated_log2": [None if not math.isfinite(v) else v for v in lt],
            }
    _suptitle(fig, title)
    handles = [
        Line2D(
            [],
            [],
            color=pair_colors[p],
            lw=2.2,
            marker="o",
            ms=6,
            mec="0.15",
            path_effects=[
                patheffects.Stroke(linewidth=3.6, foreground="0.25"),
                patheffects.Normal(),
            ],
        )
        for p in pairs
    ]
    legend = _legend_figure(
        handles, [pair_legend_labels[p] for p in pairs], title="Candidate pair"
    )
    return RenderedFigure(fig, legend, {"panels": drawn, "x_dodge": offsets})


# --------------------------------------------------------------------------- #
# 5. Deterministic side-column hit labels (for the volcano)
# --------------------------------------------------------------------------- #
@dataclass
class ColumnLabels:
    """What :func:`label_hits_in_columns` drew."""

    groups: list[list[int]]
    texts: list[str]
    lowest_box_bottom: float
    n_leader_crossings: int
    n_leader_near_misses: int = 0


def _seg_cross(p1: FloatArray, p2: FloatArray, q1: FloatArray, q2: FloatArray) -> bool:
    def orient(a: FloatArray, b: FloatArray, c: FloatArray) -> float:
        return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))

    d1, d2 = orient(q1, q2, p1), orient(q1, q2, p2)
    d3, d4 = orient(p1, p2, q1), orient(p1, p2, q2)
    return d1 * d2 < 0 and d3 * d4 < 0


def _point_seg_dist(p: FloatArray, a: FloatArray, b: FloatArray) -> float:
    ab = b - a
    t = float(np.clip(np.dot(p - a, ab) / max(float(np.dot(ab, ab)), 1e-12), 0, 1))
    return float(np.hypot(*(a + t * ab - p)))


def label_hits_in_columns(
    ax: Axes,
    x: FloatArray,
    y: FloatArray,
    idx: Sequence[int],
    labels: Sequence[str],
    *,
    y_floor: float,
    fontsize: float = 8.5,
    ring_clearance_px: float = 12.0,
) -> ColumnLabels:
    """Label ``idx`` points in two side columns kept strictly above ``y_floor``.

    Points left of x = 0 get a right-aligned column left of all of them; points at or
    right of 0 a left-aligned column right of all of them. Coincident points (the
    volcano's single-linkage grouping) share one stacked label. Rows are stacked
    around the side's mean height, lifted to clear ``y_floor``; the row order is the
    permutation (exhaustive for <= 8 labels per side) with no leader crossings and no
    leader passing within ``ring_clearance_px`` of another labelled point, then the
    shortest total leader length. The x / y limits are widened (x symmetrically) to fit
    the columns. Verified after drawing (fail loud).
    """
    fig = ax.figure
    idx_arr = np.asarray(idx, dtype=np.int64)
    if idx_arr.size != len(labels):
        raise ValueError("idx and labels lengths differ.")
    if idx_arr.size == 0:
        return ColumnLabels([], [], math.inf, 0)
    xs, ys = x[idx_arr], y[idx_arr]
    ax.scatter(
        xs, ys, s=22, facecolors="none", edgecolors="black", linewidths=0.9, zorder=4
    )
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    groups = vo._group_coincident(  # the volcano's coincident-label grouping
        (xs - x0) / (x1 - x0), (ys - y0) / (y1 - y0), _QQ_MERGE_WITHIN
    )
    group_text = ["\n".join(labels[k] for k in g) for g in groups]
    gx = np.array([float(np.mean(xs[g])) for g in groups])
    gy = np.array([float(np.mean(ys[g])) for g in groups])

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()  # type: ignore[attr-defined]
    sizes_px: list[tuple[float, float]] = []
    for text in group_text:
        probe = ax.text(0, 0, text, fontsize=fontsize)
        bb = Text.get_window_extent(probe, renderer)
        probe.remove()
        sizes_px.append((bb.width, bb.height))
    axbox = ax.get_window_extent(renderer)
    gap_px = 6.0

    def to_data_h(px: float) -> float:
        return px / axbox.height * (ax.get_ylim()[1] - ax.get_ylim()[0])

    sides = {
        "left": [k for k in range(len(groups)) if gx[k] < 0],
        "right": [k for k in range(len(groups)) if gx[k] >= 0],
    }
    # Column placement per side: "outer" (beyond the side's points, away from 0) or
    # "inner" (between the side's points and the other side's). Outer columns set a
    # LOWER bound on the symmetric half-width, inner columns an UPPER bound (they must
    # clear the other side's labelled points). Pick the feasible combination (at most
    # one inner side) with the smallest half-width.
    width_px = axbox.width
    base_half = max(abs(x0), abs(x1))

    def bounds(side: str, mode: str) -> tuple[float, float]:
        members = sides[side]
        if not members:
            return 0.0, math.inf
        frac = (max(sizes_px[k][0] for k in members) + 4 * gap_px) * 1.05 / width_px
        if mode == "outer":
            far = abs(
                min(gx[k] for k in members)
                if side == "left"
                else max(gx[k] for k in members)
            )
            if frac >= 0.5:
                return math.inf, math.inf
            # |far| / (2 half) <= 1/2 - f
            return far / (1 - 2 * frac), math.inf
        other = sides["right" if side == "left" else "left"]
        near = (
            max(gx[k] for k in members)
            if side == "left"
            else min(gx[k] for k in members)
        )
        if not other:
            limit = 0.0
        else:
            limit = (
                min(gx[k] for k in other)
                if side == "left"
                else max(gx[k] for k in other)
            )
        span = abs(limit - near)
        # span * W / (2 half) >= w  ->  half <= span / (2 f)
        return 0.0, span / (2 * frac) if span > 0 else 0.0

    best_combo: tuple[float, dict[str, str]] | None = None
    for lmode in ("outer", "inner"):
        for rmode in ("outer", "inner"):
            if lmode == rmode == "inner":
                continue
            lo_l, hi_l = bounds("left", lmode)
            lo_r, hi_r = bounds("right", rmode)
            h = max(base_half, lo_l, lo_r)
            if h <= min(hi_l, hi_r) and (best_combo is None or h < best_combo[0]):
                best_combo = (h, {"left": lmode, "right": rmode})
    if best_combo is None:
        raise AssertionError("no column placement fits the labels.")
    half, modes = best_combo
    ax.set_xlim(-half, half)
    stack_px = max(
        (sum(sizes_px[k][1] + gap_px for k in m) for m in sides.values() if m),
        default=0.0,
    )
    top_needed = max(y_floor + to_data_h(stack_px + 2 * gap_px), float(gy.max()))
    if top_needed > y1:
        ax.set_ylim(y0, top_needed + to_data_h(3 * gap_px))
    fig.canvas.draw()
    to_disp = ax.transData.transform
    others_disp = to_disp(np.column_stack([xs, ys]))

    placed: dict[int, tuple[float, float, str]] = {}  # group -> (edge x, ymid, ha)
    for side, members in sides.items():
        if not members:
            continue
        gpx = to_disp(np.column_stack([gx[members], gy[members]]))
        outer_left = (side == "left") == (modes[side] == "outer")
        if outer_left:  # column to the LEFT of the side's points, right-aligned
            edge_px = float(gpx[:, 0].min()) - 3 * gap_px
        else:  # column to the RIGHT of the side's points, left-aligned
            edge_px = float(gpx[:, 0].max()) + 3 * gap_px
        heights = [sizes_px[k][1] for k in members]
        if len(members) > 8:
            raise ValueError("more than 8 label groups on one side.")
        floor_px = float(to_disp([[0.0, y_floor]])[0, 1]) + gap_px
        top_px = float(axbox.y1) - gap_px
        centre = float(gpx[:, 1].mean())
        # Spread the rows (larger row gap) until some order clears every other
        # labelled point by ring_clearance_px; keep the least-spread such layout.
        chosen: tuple[float, tuple[int, ...], float, float] | None = None
        for spread in (1.0, 2.0, 3.0, 4.0):
            row_gap = gap_px * spread
            total = sum(heights) + row_gap * (len(members) - 1)
            start = min(top_px, max(centre + total / 2, floor_px + total))
            if start - total < floor_px - 1e-6:
                break  # no room for a wider spread
            best: tuple[float, tuple[int, ...]] | None = None
            for order in permutations(range(len(members))):
                top = start
                ends: dict[int, FloatArray] = {}
                for o in order:
                    h = heights[o]
                    ends[o] = np.array([edge_px, top - h / 2])
                    top -= h + row_gap
                cost = 0.0
                segs = {o: (gpx[o], ends[o]) for o in range(len(members))}
                for a_ in range(len(members)):
                    for b_ in range(a_ + 1, len(members)):
                        if _seg_cross(*segs[a_], *segs[b_]):
                            cost += 1e6
                    own = set(groups[members[a_]])
                    for pt_i, pt in enumerate(others_disp):
                        if pt_i in own:
                            continue
                        d = _point_seg_dist(pt, *segs[a_])
                        if d < ring_clearance_px:
                            # 1e4 per near miss + more the closer it passes
                            cost += 1e4 * (
                                1.0 + (ring_clearance_px - d) / ring_clearance_px
                            )
                    cost += float(np.hypot(*(segs[a_][1] - segs[a_][0])))
                if best is None or cost < best[0]:
                    best = (cost, order)
            assert best is not None
            if chosen is None or best[0] < chosen[0]:
                chosen = (best[0], best[1], start, row_gap)
            if best[0] < 1e4:
                break
        if chosen is None:
            raise AssertionError(f"{side} label column does not fit above the floor.")
        if chosen[0] >= 1e6:
            raise AssertionError(f"{side}: every leader layout has a crossing.")
        _, order_best, start, row_gap = chosen
        top = start
        for o in order_best:
            h = heights[o]
            ymid_px = top - h / 2
            top -= h + row_gap
            edge_d, ymid_d = ax.transData.inverted().transform([[edge_px, ymid_px]])[0]
            placed[members[o]] = (
                float(edge_d),
                float(ymid_d),
                "right" if outer_left else "left",
            )

    text_objs: list[Text] = []
    leader_pts: list[tuple[FloatArray, FloatArray]] = []
    for k, (edge, ymid, ha) in placed.items():
        text_objs.append(
            ax.text(
                edge,
                ymid,
                group_text[k],
                ha=ha,
                va="center",
                fontsize=fontsize,
                zorder=6,
                bbox={"fc": "white", "ec": "none", "alpha": 0.85, "pad": 1.0},
            )
        )
        ax.plot([gx[k], edge], [gy[k], ymid], color="0.35", lw=0.7, zorder=3)
        leader_pts.append((np.array([gx[k], gy[k]]), np.array([edge, ymid])))

    # Verify: boxes above floor, no overlaps, no crossings, inside the axes.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()  # type: ignore[attr-defined]
    axbox = ax.get_window_extent(renderer)
    boxes = [Text.get_window_extent(t, renderer) for t in text_objs]
    inv = ax.transData.inverted()
    bottoms = [float(inv.transform([[b.x0, b.y0]])[0, 1]) for b in boxes]
    if min(bottoms) <= y_floor:
        raise AssertionError("a label box reaches below the floor.")
    for i, b in enumerate(boxes):
        if b.x0 < axbox.x0 or b.x1 > axbox.x1 or b.y1 > axbox.y1:
            raise AssertionError(f"label {i} extends outside the axes.")
        for j in range(i + 1, len(boxes)):
            if b.overlaps(boxes[j]):
                raise AssertionError(f"labels {i} and {j} overlap.")
    ringed_px = ax.transData.transform(np.column_stack([xs, ys]))
    for b in boxes:
        inside = (
            (ringed_px[:, 0] >= b.x0 - 3)
            & (ringed_px[:, 0] <= b.x1 + 3)
            & (ringed_px[:, 1] >= b.y0 - 3)
            & (ringed_px[:, 1] <= b.y1 + 3)
        )
        if inside.any():
            raise AssertionError("a label box covers a labelled point.")
    disp = [(to_disp(a), to_disp(b)) for a, b in leader_pts]
    crossings = sum(
        _seg_cross(*disp[i], *disp[j])
        for i in range(len(disp))
        for j in range(i + 1, len(disp))
    )
    if crossings:
        raise AssertionError(f"{crossings} leader crossing(s).")
    near_misses = 0
    for (k, _), (a_px, b_px) in zip(placed.items(), disp, strict=True):
        own = set(groups[k])
        for pt_i, pt in enumerate(ringed_px):
            if pt_i not in own and _point_seg_dist(pt, a_px, b_px) < ring_clearance_px:
                near_misses += 1
    ordered = list(placed)
    return ColumnLabels(
        [groups[k] for k in ordered],
        [group_text[k] for k in ordered],
        min(bottoms),
        crossings,
        near_misses,
    )
