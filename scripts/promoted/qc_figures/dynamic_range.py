"""Dynamic-range / rank-abundance QC figures for a :class:`Dataset`.

Project module seeded from the plugin ``lib/figures`` template ``dynamic-range`` v0.1.
Held to the correctness charter (conventions/correctness.md): **assume nothing, verify
everything, fail loud.**

What it draws (one figure): features ranked by abundance (most -> least) on the x-axis
vs their log2 abundance on the y-axis -- the rank-abundance curve. It quantifies the
**dynamic range** (orders of magnitude) and shows whether a handful of hyper-abundant
features dominate. Two modes:

* **Whole-cohort** (default): the per-feature **median** of detected values, with a
  shaded **IQR band** (per-feature 25-75th percentile across samples) behind it.
  Optional ``highlight_features`` mark named features at their (rank, abundance);
  ``label_top_n`` limits text labels to the N most-abundant highlights (the rest are
  marked without a label) so a large highlight set (e.g. all contaminants) stays
  legible.
* **Per-class** (``class_by`` set): one independently-ranked median curve per sample
  class, registry-colored. Highlights are unavailable in this mode.

Scale (HARD refuse): the dynamic range is read on the **raw linear matrix** --
:func:`compute_dynamic_range` raises :class:`DynamicRangeScaleError` on a non-linear
``Dataset.scale``. Never-detected features are left off the curve.

DEVIATIONS FROM THE TEMPLATE (beyond header / ``__script_meta__``):

* Imports: ``loaders.data_loading.Dataset`` (project loader package) and
  ``common.figures.{colors,figure_io}`` (the promoted figure helpers) instead of
  ``common.data_loading`` / ``figures.*``.
* :func:`save_dynamic_range` builds AND saves inside ``publication_style()`` -- the
  project ``figure_io`` requires ``save_figure`` be called inside the style context.
* The whole-cohort key (median line, IQR band, highlight marker/groups) is ALWAYS
  rendered as the separate legend image; the template's on-axes ``ax.legend`` for the
  median/IQR is removed (no baked-in legend, conventions/visualization.md).
* New ``label_top_n`` (label only the top-N highlights by rank) and ``highlight_size``
  (marker area) parameters, and an ``n_samples`` count in the axes title.

This module makes NO study decisions: the detection threshold, class column, and
highlight set are the caller's.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from common.figures.colors import DEFAULT_REGISTRY_PATH, assign_colors
from common.figures.figure_io import FigureArtifacts, publication_style, save_figure
from loaders.data_loading import Dataset
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

__script_meta__: dict[str, object] = {
    "task": "stage3-qc-dynamic-range",
    "kind": "module",
    "provides": [
        "DynamicRangeScaleError",
        "DynamicRangeResult",
        "DynamicRangePlot",
        "compute_dynamic_range",
        "plot_dynamic_range",
        "save_dynamic_range",
    ],
    "uses": [
        "loaders.data_loading",
        "common.figures.colors",
        "common.figures.figure_io",
    ],
    "seeded_from": {"template": "dynamic-range", "version": "0.1"},
    "description": (
        "Dynamic-range / rank-abundance QC figures from a Dataset: features ranked by "
        "median detected abundance vs log2 abundance, whole-cohort median + IQR band "
        "(default) with optional highlight_features (label_top_n limits text labels), "
        "or a per-class overlay (class_by). Hard-refuses a non-linear scale. "
        "Dual-export plus a separate legend image (always, in whole-cohort mode). "
        "Study-agnostic; fail-loud."
    ),
}

# Default color for highlighted features when no ``highlight_groups`` are given.
DEFAULT_HIGHLIGHT_COLOR = "#0072B2"
_BAND_COLOR = "#bdbdbd"
_CURVE_COLOR = "#333333"
_BAND_ALPHA = 0.55


class DynamicRangeScaleError(ValueError):
    """Raised when the dynamic range is requested on a non-linear abundance scale."""


@dataclass(frozen=True)
class DynamicRangeResult:
    """The ranked per-feature abundance summary underlying a dynamic-range figure.

    Only features detected in at least one sample appear; arrays are ordered most ->
    least abundant.

    Attributes
    ----------
    feature_names_ranked:
        ``(k,)`` feature ids, most -> least abundant.
    log2_median, log2_q25, log2_q75:
        ``(k,)`` log2 of the per-feature median / 25th / 75th percentile of detected
        values across samples.
    n_detected:
        ``(k,)`` number of samples in which each ranked feature is detected.
    dynamic_range_orders:
        ``log10(max median / min median)`` -- orders of magnitude spanned.
    n_features_total:
        Features in the input Dataset.
    n_features_detected:
        Features detected in >= 1 sample (``k``).
    n_samples:
        Samples in the input Dataset.
    """

    feature_names_ranked: np.ndarray
    log2_median: np.ndarray
    log2_q25: np.ndarray
    log2_q75: np.ndarray
    n_detected: np.ndarray
    dynamic_range_orders: float
    n_features_total: int
    n_features_detected: int
    n_samples: int

    @property
    def ranks(self) -> np.ndarray:
        """``(k,)`` ranks ``1..k`` aligned to the ranked arrays."""
        return np.arange(1, self.n_features_detected + 1)


@dataclass
class DynamicRangePlot:
    """A rendered dynamic-range figure plus its companion legend figure."""

    figure: Figure
    legend_figure: Figure | None
    result: DynamicRangeResult
    color_map: dict[str, str]


# --------------------------------------------------------------------------- #
# Compute
# --------------------------------------------------------------------------- #


def compute_dynamic_range(
    dataset: Dataset, *, min_intensity: float = 0.0
) -> DynamicRangeResult:
    """Rank features by median detected abundance; return the ranked log2 summary.

    Raises :class:`DynamicRangeScaleError` unless ``dataset.scale == "linear"``, or
    ``ValueError`` on a non-2D / empty Dataset or a feature-name mismatch.
    """
    if dataset.scale != "linear":
        raise DynamicRangeScaleError(
            f"Dynamic range requires linear-scale abundances but the Dataset is on "
            f"scale {dataset.scale!r}. 'detected' (> {min_intensity}) and the log2 "
            f"y-axis are raw-linear properties. Run it on the unnormalized data."
        )
    abundances = np.asarray(dataset.abundances, dtype=float)
    if abundances.ndim != 2:
        raise ValueError(
            f"abundances must be 2D (n_samples, n_features); got {abundances.shape}."
        )
    names = np.asarray(dataset.feature_names)
    n_total = abundances.shape[1]
    if names.shape[0] != n_total:
        raise ValueError(
            f"feature_names has {names.shape[0]} entries but abundances has {n_total} "
            f"features; they are not aligned."
        )
    if abundances.shape[0] < 1 or n_total < 1:
        raise ValueError(
            f"Dataset is empty ({abundances.shape[0]} samples x {n_total} features)."
        )

    median, q25, q75, detected_count = _detected_quartiles(abundances, min_intensity)
    keep = detected_count >= 1
    if not bool(keep.any()):
        raise ValueError(
            "no feature is detected in any sample (all <= min_intensity); none to rank."
        )

    med_keep = median[keep]
    # Stable sort on the negated median: ties keep input order, deterministic.
    order = np.argsort(-med_keep, kind="stable")
    med_ranked = med_keep[order]
    orders = (
        float(np.log10(med_ranked[0] / med_ranked[-1]))
        if med_ranked[-1] > 0
        else float("nan")
    )

    return DynamicRangeResult(
        feature_names_ranked=names[keep][order],
        log2_median=np.log2(med_ranked),
        log2_q25=np.log2(q25[keep][order]),
        log2_q75=np.log2(q75[keep][order]),
        n_detected=detected_count[keep][order],
        dynamic_range_orders=orders,
        n_features_total=int(n_total),
        n_features_detected=int(keep.sum()),
        n_samples=int(abundances.shape[0]),
    )


def _detected_quartiles(
    abundances: np.ndarray, min_intensity: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-feature median / q25 / q75 of detected (> min_intensity) values + count."""
    mask = np.isfinite(abundances) & (abundances > min_intensity)
    positive = np.where(mask, abundances, np.nan)
    with warnings.catch_warnings():
        # all-NaN columns -> NaN (then dropped as never-detected)
        warnings.simplefilter("ignore", RuntimeWarning)
        median = np.nanmedian(positive, axis=0)
        q25 = np.nanpercentile(positive, 25, axis=0)
        q75 = np.nanpercentile(positive, 75, axis=0)
    return median, q25, q75, mask.sum(axis=0)


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #


def plot_dynamic_range(
    dataset: Dataset,
    *,
    class_by: str | None = None,
    highlight_features: Mapping[str, str] | None = None,
    highlight_groups: Mapping[str, str] | None = None,
    highlight_color: str = DEFAULT_HIGHLIGHT_COLOR,
    highlight_label: str = "highlighted",
    highlight_size: float = 30.0,
    label_top_n: int | None = None,
    min_intensity: float = 0.0,
    log_rank: bool = False,
    title: str | None = None,
    legend_title: str | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
) -> DynamicRangePlot:
    """Draw the rank-abundance curve for one Dataset (whole-cohort or per-class).

    Parameters
    ----------
    dataset:
        A single linear-scale :class:`Dataset` (raised otherwise).
    class_by:
        Metadata column to split the curve by. Mutually exclusive with highlights.
    highlight_features:
        ``{feature_id: display label}`` to mark at their (rank, log2 median). Unknown
        or never-detected ids are refused.
    highlight_groups:
        Optional ``{feature_id: group}`` coloring highlights by registry group (the >8
        guard applies). Every highlighted id must have a group if given.
    highlight_color:
        Color for highlights when ``highlight_groups`` is not given.
    highlight_label:
        Legend label for the single-color highlight marker (ignored with groups).
    highlight_size:
        Highlight marker area (points^2).
    label_top_n:
        Text-label only the N most-abundant highlights (others are marked only).
        ``None`` labels every highlight; ``0`` labels none.
    min_intensity:
        Detection threshold (default ``0.0``).
    log_rank:
        Log-scaled rank x-axis.
    title:
        Optional figure suptitle.
    legend_title:
        Title for the companion legend figure.
    registry_path:
        Color registry JSON. Default ``state/color_registry.json``.
    persist_colors:
        Write newly assigned colors back to the registry (default ``True``).

    Raises
    ------
    DynamicRangeScaleError
        If the Dataset is not on the ``"linear"`` scale.
    ValueError
        On an empty Dataset, ``class_by`` and highlights both given, a missing
        ``class_by`` column, an unknown / never-detected highlight, an incomplete
        ``highlight_groups``, or a negative ``label_top_n``.
    CategoricalPaletteExceededError
        If ``class_by`` / ``highlight_groups`` exceed the palette capacity.
    """
    if class_by is not None and highlight_features:
        raise ValueError(
            "class_by and highlight_features are mutually exclusive: per-class curves "
            "are ranked independently, so a highlight has no single rank. Pass one."
        )
    if label_top_n is not None and label_top_n < 0:
        raise ValueError(f"label_top_n must be >= 0 or None; got {label_top_n}.")
    result = compute_dynamic_range(dataset, min_intensity=min_intensity)

    color_map: dict[str, str] = {}
    legend_figure: Figure | None = None
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    try:
        if class_by is not None:
            color_map = _draw_per_class(
                ax,
                dataset,
                class_by,
                min_intensity=min_intensity,
                registry_path=registry_path,
                persist_colors=persist_colors,
            )
            handles: list[Artist] = [
                Line2D([0], [0], color=color_map[k], linewidth=2.5) for k in color_map
            ]
            legend_figure = _legend_figure(
                handles,
                list(color_map),
                legend_title if legend_title is not None else class_by,
            )
        else:
            _draw_curve(ax, result)
            handles = [
                Line2D([0], [0], color=_CURVE_COLOR, linewidth=2.0),
                Patch(facecolor=_BAND_COLOR, alpha=_BAND_ALPHA, edgecolor="none"),
            ]
            labels = [
                "median of detected values",
                "IQR across samples (25th-75th pct)",
            ]
            if highlight_features:
                color_map = _draw_highlights(
                    ax,
                    dataset,
                    result,
                    highlight_features,
                    highlight_groups,
                    highlight_color,
                    highlight_size=highlight_size,
                    label_top_n=label_top_n,
                    registry_path=registry_path,
                    persist_colors=persist_colors,
                )
                marker_colors = (
                    color_map
                    if highlight_groups
                    else {highlight_label: highlight_color}
                )
                for name, color in marker_colors.items():
                    handles.append(_marker_handle(color))
                    labels.append(name)
            legend_figure = _legend_figure(handles, labels, legend_title or "")

        if log_rank:
            ax.set_xscale("log")
        ax.set_xlabel("abundance rank (1 = most abundant)")
        ax.set_ylabel("log2 intensity (median of detected, a.u.)")
        ax.grid(True, linestyle=":", alpha=0.4)
        if title is not None:
            fig.suptitle(title, fontsize=15, weight="bold")
    except BaseException:
        plt.close(fig)
        if legend_figure is not None:
            plt.close(legend_figure)
        raise

    return DynamicRangePlot(
        figure=fig, legend_figure=legend_figure, result=result, color_map=color_map
    )


def save_dynamic_range(
    dataset: Dataset,
    output_dir: str | Path,
    base_name: str,
    *,
    class_by: str | None = None,
    highlight_features: Mapping[str, str] | None = None,
    highlight_groups: Mapping[str, str] | None = None,
    highlight_color: str = DEFAULT_HIGHLIGHT_COLOR,
    highlight_label: str = "highlighted",
    highlight_size: float = 30.0,
    label_top_n: int | None = None,
    min_intensity: float = 0.0,
    log_rank: bool = False,
    title: str | None = None,
    legend_title: str | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
    dpi: int = 300,
) -> tuple[FigureArtifacts, DynamicRangePlot]:
    """Render (:func:`plot_dynamic_range`) and save inside ``publication_style()``.

    Writes ``<base>.{svg,png}`` and ``<base>.legend.{svg,png}``; figures are closed
    even if saving fails. Returns the artifacts and the (closed) plot, whose
    ``result`` / ``color_map`` the caller may report.
    """
    with publication_style():
        plot = plot_dynamic_range(
            dataset,
            class_by=class_by,
            highlight_features=highlight_features,
            highlight_groups=highlight_groups,
            highlight_color=highlight_color,
            highlight_label=highlight_label,
            highlight_size=highlight_size,
            label_top_n=label_top_n,
            min_intensity=min_intensity,
            log_rank=log_rank,
            title=title,
            legend_title=legend_title,
            registry_path=registry_path,
            persist_colors=persist_colors,
        )
        artifacts = save_figure(
            plot.figure, output_dir, base_name, legend_fig=plot.legend_figure, dpi=dpi
        )
    return artifacts, plot


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #


def _draw_curve(ax: Axes, result: DynamicRangeResult) -> None:
    """The whole-cohort median curve + IQR band (key goes in the legend image)."""
    ranks = result.ranks
    ax.fill_between(
        ranks,
        result.log2_q25,
        result.log2_q75,
        color=_BAND_COLOR,
        alpha=_BAND_ALPHA,
        linewidth=0,
    )
    ax.plot(ranks, result.log2_median, color=_CURVE_COLOR, linewidth=1.7)
    ax.set_title(
        f"{result.dynamic_range_orders:.1f} orders of magnitude  "
        f"({result.n_features_detected} of {result.n_features_total} features "
        f"detected; {result.n_samples} samples)",
        loc="left",
        fontsize=12,
    )


def _draw_per_class(
    ax: Axes,
    dataset: Dataset,
    class_by: str,
    *,
    min_intensity: float,
    registry_path: str | Path,
    persist_colors: bool,
) -> dict[str, str]:
    """One independently-ranked median curve per class; returns the color map."""
    if class_by not in dataset.metadata.columns:
        raise ValueError(
            f"class_by {class_by!r} is not a metadata column "
            f"{list(dataset.metadata.columns)}."
        )
    series = dataset.metadata[class_by]
    if bool(series.isna().any()):
        raise ValueError(
            f"class_by {class_by!r} has {int(series.isna().sum())} missing value(s); "
            f"resolve or relabel them first."
        )
    classes = series.to_numpy().astype(str)
    abundances = np.asarray(dataset.abundances, dtype=float)
    ordered = _ordered_unique(classes)
    color_map = assign_colors(
        class_by, ordered, registry_path=registry_path, persist=persist_colors
    )
    for cls in ordered:
        members = classes == cls
        median, _, _, count = _detected_quartiles(abundances[members], min_intensity)
        med = np.sort(median[count >= 1])[::-1]
        ax.plot(
            np.arange(1, med.size + 1),
            np.log2(med),
            color=color_map[cls],
            linewidth=2.0,
        )
    ax.set_title("Dynamic range by sample class", loc="left", fontsize=12)
    result: dict[str, str] = dict(color_map)
    return result


@dataclass(frozen=True)
class _Highlight:
    """One resolved highlight: its rank position, log2 abundance, label, and color."""

    x: int
    y: float
    text: str
    color: str


def _draw_highlights(
    ax: Axes,
    dataset: Dataset,
    result: DynamicRangeResult,
    highlight_features: Mapping[str, str],
    highlight_groups: Mapping[str, str] | None,
    highlight_color: str,
    *,
    highlight_size: float,
    label_top_n: int | None,
    registry_path: str | Path,
    persist_colors: bool,
) -> dict[str, str]:
    """Mark the highlighted features, label the top-N; return the group color map."""
    rank_of = {str(name): i + 1 for i, name in enumerate(result.feature_names_ranked)}
    known = set(np.asarray(dataset.feature_names).astype(str))
    color_map: dict[str, str] = {}
    if highlight_groups is not None:
        missing = [f for f in highlight_features if f not in highlight_groups]
        if missing:
            raise ValueError(
                f"highlight_groups is missing group(s) for {missing[:5]}; every "
                f"highlighted feature needs a group when highlight_groups is given."
            )
        color_map = assign_colors(
            "highlight_group",
            _ordered_unique(np.asarray(list(highlight_groups.values()))),
            registry_path=registry_path,
            persist=persist_colors,
        )

    items: list[_Highlight] = []
    for fid, label in highlight_features.items():
        if fid not in known:
            raise ValueError(f"highlighted feature {fid!r} is not in the Dataset.")
        if fid not in rank_of:
            raise ValueError(
                f"highlighted feature {fid!r} is never detected, so it has no "
                f"abundance to mark; drop it from highlight_features."
            )
        idx = rank_of[fid] - 1
        color = (
            color_map[highlight_groups[fid]] if highlight_groups else highlight_color
        )
        ax.scatter(
            [rank_of[fid]],
            [result.log2_median[idx]],
            s=highlight_size,
            color=color,
            edgecolor="black",
            linewidth=0.5,
            zorder=3,
        )
        items.append(
            _Highlight(
                x=rank_of[fid],
                y=float(result.log2_median[idx]),
                text=label,
                color=color,
            )
        )

    items.sort(key=lambda h: h.x)
    to_label = items if label_top_n is None else items[:label_top_n]
    finite = result.log2_median[np.isfinite(result.log2_median)]
    _place_labels(
        ax,
        to_label,
        result.n_features_detected,
        float(finite.max()),
        float(finite.min()),
    )
    return color_map


def _place_labels(
    ax: Axes, items: list[_Highlight], n: int, y_top: float, y_bot: float
) -> None:
    """Leader-line labels, de-collided and clamped inside the axes.

    Head points (left half) are labelled in an evenly spaced column to the right, one
    leader per point (deviation: the template clamped each label to its point's y);
    tail points (right half) are labelled below their points. Both stacks are
    de-collided and the y-range padded so no label leaves the axes.
    """
    ax.set_ylim(y_bot - 1.5, y_top + 1.0)
    gap = (y_top - y_bot) * 0.058

    head = sorted((h for h in items if h.x < n * 0.5), key=lambda h: h.y, reverse=True)
    col_x = n * 0.18
    # Evenly spaced column starting at the top-most highlight: labels sit in the open
    # space above the (monotone-decreasing) curve, never beside an unlabelled point.
    cur = min(y_top, head[0].y) if head else y_top
    for h in head:
        _annotate(ax, h, col_x, cur, "left")
        cur -= gap

    tail = sorted((h for h in items if h.x >= n * 0.5), key=lambda h: h.y, reverse=True)
    cur = min((h.y for h in tail), default=y_bot) - gap * 1.6
    for h in tail:
        ly = max(cur, y_bot + gap)
        cur = ly - gap
        _annotate(ax, h, float(h.x - n * 0.01), ly, "right")


def _annotate(
    ax: Axes, item: _Highlight, label_x: float, label_y: float, ha: str
) -> None:
    """One leader-line annotation from its label position back to the point."""
    ax.annotate(
        item.text,
        xy=(item.x, item.y),
        xytext=(label_x, label_y),
        fontsize=10,
        color="black",
        va="center",
        ha=ha,
        arrowprops={"arrowstyle": "-", "color": item.color, "lw": 0.8},
    )


def _ordered_unique(values: np.ndarray) -> list[str]:
    """Distinct values in first-seen order (stable, study-agnostic)."""
    seen: list[str] = []
    for value in values.astype(str):
        if value not in seen:
            seen.append(value)
    return seen


def _marker_handle(color: str) -> Line2D:
    """A legend handle matching the highlight scatter marker."""
    return Line2D(
        [0],
        [0],
        linestyle="none",
        marker="o",
        markersize=7,
        markerfacecolor=color,
        markeredgecolor="black",
        markeredgewidth=0.5,
    )


# --------------------------------------------------------------------------- #
# Legend figure (rendered separately so it never overlaps the plot)
# --------------------------------------------------------------------------- #


def _legend_figure(
    handles: list[Artist], labels: list[str], legend_title: str
) -> Figure:
    """Standalone key: one entry per handle."""
    height = max(1.2, 0.38 * len(labels) + 0.6)
    width = max(3.4, 0.095 * max(len(s) for s in labels) + 1.2)
    fig, ax = plt.subplots(figsize=(width, height))
    ax.axis("off")
    ax.legend(
        handles,
        labels,
        title=legend_title or None,
        loc="center",
        frameon=False,
        fontsize=11,
        title_fontsize=12,
    )
    return fig
