"""Identification-depth (detected-features-per-run) QC bar charts for a `Dataset`.

PROJECT COPY seeded from the plugin template ``lib/figures/id_depth.py`` (``id-depth``
v0.1). Held to the correctness charter: **assume nothing, verify everything, fail
loud.** Project adaptations (all study-agnostic in behaviour):

* package layout -- ``loaders.data_loading`` for :class:`Dataset` and the promoted
  ``common.figures.{colors,figure_io}`` helpers instead of ``common.data_loading`` /
  ``figures.*``;
* the layout scales down for small studies (the template's 12-inch minimum width and
  5-pt tick labels were sized for 50+ runs; this study has 8), see :func:`_figsize`;
* ``ylabel`` may be a per-panel mapping (protein groups vs peptides count different
  things, so each panel names its own unit);
* optional ``annotate_counts`` -- the terse per-bar count printed above each bar
  (legible for small n; off by default because it crowds at large n);
* optional ``divider_by`` -- a thin dotted divider (plus a short text label) between
  *contiguous* runs of a metadata column, e.g. the acquisition batch. It is a
  position cue, not a color channel, so the single categorical color channel stays
  the bar color (the template's "no annotation stripes" rule is kept). It refuses a
  column whose values are not contiguous in the sample order (a divider would lie).

What it draws (one figure): a stack of per-sample bar charts, one panel per labelled
:class:`Dataset`, sharing the same samples. Each bar is the number of *detected*
features in that sample: finite and strictly greater than ``min_intensity`` (default
``0.0``), counted per sample over the features (``axis=1``). Bars are drawn in the
Dataset's row order -- order by acquisition upstream. Detection is a property of the
raw linear matrix, so a non-linear ``Dataset.scale`` is a HARD refuse
(:class:`IdDepthScaleError`). Each panel carries a dashed reference median, over the
``reference_mask`` subset when given, else all samples. The legend is rendered as its
own figure and saved beside the plot as ``<base>.legend.{svg,png}``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from common.figures.colors import DEFAULT_REGISTRY_PATH, assign_colors
from common.figures.figure_io import FigureArtifacts, publication_style, save_figure
from loaders.data_loading import Dataset
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

__script_meta__: dict[str, object] = {
    "template": {"name": "id-depth", "version": "0.1"},
    "kind": "module",
    "provides": [
        "IdDepthScaleError",
        "IdDepthResult",
        "IdDepthPlot",
        "compute_detection_counts",
        "plot_id_depth",
        "save_id_depth",
    ],
    "uses": [
        "loaders.data_loading",
        "common.figures.colors",
        "common.figures.figure_io",
    ],
    "seeded_from": "lib/figures/id_depth.py (id-depth v0.1)",
    "description": (
        "Identification-depth QC bar charts from one or more Datasets: detected "
        "features per sample (finite and > min_intensity), one stacked panel per "
        "feature level, bars in acquisition order. Bars colored by a categorical "
        "metadata column via the project registry (>8 guard); per-panel reference "
        "median; optional per-bar counts and contiguous-group dividers; hard-refuses a "
        "non-linear scale; dual-export plus a separate legend image."
    ),
}

# Bar color when no ``color_by`` is given: a neutral slate, outside the categorical
# palette so it consumes no registry slot and reads as "uncategorized".
UNCATEGORIZED_BAR_COLOR = "#6e7f8d"
_INK = "#333333"


class IdDepthScaleError(ValueError):
    """Raised when a detection count is requested on a non-linear abundance scale.

    On a log or centered scale a ``0`` is an ordinary value and ``> 0`` no longer means
    "detected", so the count is meaningless. Count depth on the unnormalized data.
    """


@dataclass(frozen=True)
class IdDepthResult:
    """The per-sample detection counts underlying an identification-depth figure.

    Attributes
    ----------
    sample_ids:
        ``(n_samples,)`` the sample identifiers (the shared order across panels).
    counts:
        ``{label: (n_samples,) detected-feature count}`` in panel order.
    reference_median:
        ``{label: median detected count}`` over the ``reference_mask`` subset when
        given, else all samples. ``NaN`` only if the subset was empty.
    reference_is_subset:
        ``True`` when the median was taken over a caller-supplied subset.
    """

    sample_ids: np.ndarray
    counts: dict[str, np.ndarray]
    reference_median: dict[str, float]
    reference_is_subset: bool


@dataclass
class IdDepthPlot:
    """A rendered identification-depth figure plus its companion legend figure."""

    figure: Figure
    legend_figure: Figure | None
    result: IdDepthResult
    color_map: dict[str, str]


# --------------------------------------------------------------------------- #
# Compute
# --------------------------------------------------------------------------- #


def compute_detection_counts(
    dataset: Dataset, *, min_intensity: float = 0.0
) -> np.ndarray:
    """Per-sample count of detected features (finite and ``> min_intensity``).

    Raises :class:`IdDepthScaleError` unless ``dataset.scale == "linear"``, or
    ``ValueError`` if abundances are not 2D. ``NaN`` and values ``<= min_intensity``
    both count as not detected.
    """
    if dataset.scale != "linear":
        raise IdDepthScaleError(
            f"Identification depth (a detected-feature count) requires linear-scale "
            f"abundances but the Dataset is on scale {dataset.scale!r}. On log/centered"
            f" scales a 0 is an ordinary value, so '> {min_intensity}' no longer means "
            f"'detected'. Count depth on the unnormalized data, before any transform."
        )
    abundances = np.asarray(dataset.abundances, dtype=float)
    if abundances.ndim != 2:
        raise ValueError(
            f"abundances must be 2D (n_samples, n_features); got {abundances.shape}."
        )
    detected = np.isfinite(abundances) & (abundances > min_intensity)
    return np.asarray(detected.sum(axis=1), dtype=int)


def _compute_result(
    datasets: Mapping[str, Dataset],
    sample_ids: np.ndarray,
    *,
    min_intensity: float,
    reference_mask: np.ndarray | None,
) -> IdDepthResult:
    """Per-sample counts + the per-level reference median (scale guard fires here)."""
    counts = {
        label: compute_detection_counts(ds, min_intensity=min_intensity)
        for label, ds in datasets.items()
    }
    reference: dict[str, float] = {}
    for label, count in counts.items():
        subset = count[reference_mask] if reference_mask is not None else count
        reference[label] = float(np.median(subset)) if subset.size else float("nan")
    return IdDepthResult(
        sample_ids=sample_ids,
        counts=counts,
        reference_median=reference,
        reference_is_subset=reference_mask is not None,
    )


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #


def plot_id_depth(
    datasets: Mapping[str, Dataset],
    *,
    color_by: str | None = None,
    min_intensity: float = 0.0,
    reference_mask: np.ndarray | None = None,
    show_reference_line: bool = True,
    ylabel: str | Mapping[str, str] = "Detected features",
    title: str | None = None,
    legend_title: str | None = None,
    annotate_counts: bool = False,
    divider_by: str | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
) -> IdDepthPlot:
    """Stack a per-sample identification-depth bar chart for each labelled level.

    Parameters
    ----------
    datasets:
        Ordered ``{label: Dataset}`` (insertion order = top-to-bottom panel order).
        Every Dataset must be ``"linear"`` and share the same samples in the same order.
    color_by:
        Metadata column (read from the first level) coloring each bar through the
        project color registry (capped at eight). ``None`` draws uniform neutral bars
        and no legend. Missing values are refused.
    min_intensity:
        Detection threshold (detected when finite and ``> min_intensity``).
    reference_mask:
        Optional boolean over the samples selecting the reference-median subset.
    show_reference_line:
        Draw + label the per-panel reference median (default ``True``).
    ylabel:
        Y-axis label for every panel, or a ``{label: ylabel}`` mapping covering every
        panel.
    title:
        Optional figure suptitle (keep it short: the comparison + processing state).
    legend_title:
        Title for the companion legend figure; defaults to ``color_by``.
    annotate_counts:
        Print each bar's count above it (default ``False``).
    divider_by:
        Metadata column whose *contiguous* groups are separated by a dotted divider and
        labelled at the top of the first panel; refused if a value is non-contiguous.
    registry_path, persist_colors:
        The color registry JSON and whether newly assigned colors are written back.

    Raises
    ------
    IdDepthScaleError
        If any Dataset is not on the ``"linear"`` scale.
    ValueError
        On empty/misaligned inputs, a bad ``reference_mask``, a missing or incomplete
        ``color_by`` / ``divider_by`` column, a non-contiguous ``divider_by``, or a
        ``ylabel`` mapping that does not cover every panel.
    CategoricalPaletteExceededError
        If ``color_by`` exceeds the palette's categorical capacity (the >8 guard).
    """
    if len(datasets) == 0:
        raise ValueError("datasets is empty; pass >= 1 labelled Dataset to plot.")

    # Validate + compute BEFORE creating the figure so the likely failures leak nothing.
    sample_ids = _validate_sample_alignment(datasets)
    n_samples = sample_ids.size
    reference_mask = _validate_reference_mask(reference_mask, n_samples)
    first = next(iter(datasets.values()))
    bar_values = _resolve_column(first, color_by)
    divider_values = _resolve_column(first, divider_by)
    groups = _contiguous_groups(divider_values, divider_by)
    labels = list(datasets.keys())
    ylabels = _resolve_ylabels(ylabel, labels)
    result = _compute_result(
        datasets, sample_ids, min_intensity=min_intensity, reference_mask=reference_mask
    )
    title_for_legend = legend_title if legend_title is not None else (color_by or "")

    color_map: dict[str, str] = {}
    legend_figure: Figure | None = None
    with publication_style():
        fig, raw_axes = plt.subplots(
            len(labels),
            1,
            figsize=_figsize(len(labels), n_samples),
            sharex=True,
            constrained_layout=True,
        )
        axes = list(np.atleast_1d(raw_axes))
        try:
            if bar_values is not None:
                color_map = assign_colors(
                    str(color_by),
                    list(bar_values),
                    registry_path=registry_path,
                    persist=persist_colors,
                )
                bar_colors = [color_map[v] for v in bar_values]
            else:
                bar_colors = [UNCATEGORIZED_BAR_COLOR] * n_samples

            for i, (ax, label) in enumerate(zip(axes, labels, strict=True)):
                _draw_bars(
                    ax,
                    result.counts[label],
                    bar_colors,
                    reference_median=result.reference_median[label],
                    reference_is_subset=result.reference_is_subset,
                    show_reference_line=show_reference_line,
                    annotate_counts=annotate_counts,
                    panel_label=label,
                    ylabel=ylabels[label],
                )
                _draw_dividers(ax, groups, label_groups=i == 0)

            tick_size = 10 if n_samples <= 30 else 5
            axes[-1].set_xticks(np.arange(n_samples))
            axes[-1].set_xticklabels(
                sample_ids,
                rotation=45 if n_samples <= 30 else 90,
                ha="right" if n_samples <= 30 else "center",
                rotation_mode="anchor" if n_samples <= 30 else "default",
                fontsize=tick_size,
            )
            axes[-1].set_xlabel(
                f"Sample (acquisition order, n={n_samples})", fontsize=11
            )
            if title is not None:
                fig.suptitle(title, fontsize=13, weight="bold")

            if bar_values is not None:
                legend_figure = _legend_figure(bar_values, color_map, title_for_legend)
        except BaseException:
            plt.close(fig)
            raise

    return IdDepthPlot(
        figure=fig, legend_figure=legend_figure, result=result, color_map=color_map
    )


def save_id_depth(
    datasets: Mapping[str, Dataset],
    output_dir: str | Path,
    base_name: str,
    *,
    color_by: str | None = None,
    min_intensity: float = 0.0,
    reference_mask: np.ndarray | None = None,
    show_reference_line: bool = True,
    ylabel: str | Mapping[str, str] = "Detected features",
    title: str | None = None,
    legend_title: str | None = None,
    annotate_counts: bool = False,
    divider_by: str | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
    dpi: int = 300,
) -> tuple[FigureArtifacts, IdDepthResult, dict[str, str]]:
    """Render (:func:`plot_id_depth`) and dual-export; return artifacts + result.

    Writes ``<base>.{svg,png}`` and, when ``color_by`` is set,
    ``<base>.legend.{svg,png}``. Both figures are closed even if saving fails. The
    save runs inside :func:`publication_style` so SVG settings apply at save time.
    """
    plot = plot_id_depth(
        datasets,
        color_by=color_by,
        min_intensity=min_intensity,
        reference_mask=reference_mask,
        show_reference_line=show_reference_line,
        ylabel=ylabel,
        title=title,
        legend_title=legend_title,
        annotate_counts=annotate_counts,
        divider_by=divider_by,
        registry_path=registry_path,
        persist_colors=persist_colors,
    )
    with publication_style():
        artifacts = save_figure(
            plot.figure, output_dir, base_name, legend_fig=plot.legend_figure, dpi=dpi
        )
    return artifacts, plot.result, plot.color_map


# --------------------------------------------------------------------------- #
# Validation (fail loud)
# --------------------------------------------------------------------------- #


def _validate_sample_alignment(datasets: Mapping[str, Dataset]) -> np.ndarray:
    """Refuse stacking panels whose samples differ; return the shared sample ids."""
    reference_ids: np.ndarray | None = None
    reference_label = ""
    for label, ds in datasets.items():
        abundances = np.asarray(ds.abundances, dtype=float)
        if abundances.ndim != 2:
            raise ValueError(
                f"Dataset {label!r} abundances must be 2D (n_samples, n_features); "
                f"got {abundances.shape}."
            )
        n_samples, n_features = abundances.shape
        if n_samples < 1 or n_features < 1:
            raise ValueError(
                f"Dataset {label!r} is empty ({n_samples} samples x {n_features} "
                f"features); need >= 1 of each to draw a depth bar chart."
            )
        ids = ds.metadata.index.to_numpy().astype(str)
        if ids.shape[0] != n_samples:
            raise ValueError(
                f"Dataset {label!r} has {ids.shape[0]} metadata rows but {n_samples} "
                f"abundance samples; it is not row-aligned."
            )
        if reference_ids is None:
            reference_ids = ids
            reference_label = label
            continue
        if ids.shape != reference_ids.shape or not np.array_equal(ids, reference_ids):
            raise ValueError(
                f"Dataset {label!r} has different samples (or a different order) than "
                f"{reference_label!r}; the stacked levels must align bar-for-bar."
            )
    if reference_ids is None:  # unreachable: datasets checked non-empty by caller
        raise ValueError("datasets is empty.")
    return reference_ids


def _validate_reference_mask(
    reference_mask: np.ndarray | None, n_samples: int
) -> np.ndarray | None:
    """Coerce + length-check the reference mask (a boolean over the samples)."""
    if reference_mask is None:
        return None
    mask = np.asarray(reference_mask, dtype=bool)
    if mask.ndim != 1 or mask.shape[0] != n_samples:
        raise ValueError(
            f"reference_mask must be a 1D boolean array of length n_samples="
            f"{n_samples}; got shape {mask.shape}."
        )
    return mask


def _resolve_column(dataset: Dataset, column: str | None) -> np.ndarray | None:
    """Validate + extract a metadata column as per-sample strings (or ``None``)."""
    if column is None:
        return None
    if column not in dataset.metadata.columns:
        raise ValueError(
            f"{column!r} is not a metadata column {list(dataset.metadata.columns)}."
        )
    series = dataset.metadata[column]
    if bool(series.isna().any()):
        n_bad = int(series.isna().sum())
        raise ValueError(
            f"Metadata column {column!r} has {n_bad} missing value(s); resolve or "
            f"relabel them first."
        )
    return np.asarray(series.to_numpy(), dtype=str)


def _contiguous_groups(
    values: np.ndarray | None, column: str | None
) -> list[tuple[str, int, int]]:
    """``[(value, start, stop)]`` contiguous blocks; refuse a non-contiguous value."""
    if values is None:
        return []
    groups: list[tuple[str, int, int]] = []
    start = 0
    for i in range(1, values.size + 1):
        if i == values.size or values[i] != values[start]:
            groups.append((str(values[start]), start, i))
            start = i
    names = [g[0] for g in groups]
    if len(set(names)) != len(names):
        raise ValueError(
            f"divider_by {column!r} is not contiguous in the sample order "
            f"(blocks: {names}); a divider would misrepresent the grouping."
        )
    return groups


def _resolve_ylabels(
    ylabel: str | Mapping[str, str], labels: list[str]
) -> dict[str, str]:
    """Expand a single y-label to every panel, or check a mapping covers them all."""
    if isinstance(ylabel, str):
        return dict.fromkeys(labels, ylabel)
    missing = [lab for lab in labels if lab not in ylabel]
    if missing:
        raise ValueError(f"ylabel mapping is missing panel(s) {missing}.")
    return {lab: ylabel[lab] for lab in labels}


# --------------------------------------------------------------------------- #
# Layout & drawing
# --------------------------------------------------------------------------- #


def _figsize(n_levels: int, n_samples: int) -> tuple[float, float]:
    """Width scales with sample count (6.5-inch floor); height with the panel stack."""
    width = max(6.5, 2.5 + n_samples * 0.45) if n_samples <= 30 else n_samples * 0.18
    height = 2.6 * n_levels + 1.3
    return (max(width, 6.5), height)


def _draw_bars(
    ax: Axes,
    counts: np.ndarray,
    bar_colors: list[str],
    *,
    reference_median: float,
    reference_is_subset: bool,
    show_reference_line: bool,
    annotate_counts: bool,
    panel_label: str,
    ylabel: str,
) -> None:
    """Draw one bar per sample (colored), the reference median, optional counts."""
    x = np.arange(counts.size)
    ax.bar(x, counts, color=bar_colors, edgecolor="white", linewidth=0.3, width=0.8)

    if annotate_counts:
        for xi, c in zip(x, counts, strict=True):
            # Lift the label a few points off the bar so its white box (which masks
            # the reference line behind the text) never notches the bar top.
            ax.annotate(
                f"{int(c):,}",
                xy=(float(xi), float(c)),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
                color=_INK,
                zorder=4,  # above the reference line, which passes behind the label
                bbox={"boxstyle": "square,pad=0.1", "fc": "white", "ec": "none"},
            )

    if show_reference_line and np.isfinite(reference_median):
        ax.axhline(
            reference_median, color=_INK, linestyle="--", linewidth=1.0, zorder=3
        )
        tag = "ref. median" if reference_is_subset else "median"
        ax.text(
            1.01,
            reference_median,
            f"{tag}\n{reference_median:,.0f}",
            transform=ax.get_yaxis_transform(),
            va="center",
            ha="left",
            fontsize=8,
            color=_INK,
            clip_on=False,
        )

    ymax = float(np.max(counts)) if counts.size else 1.0
    ax.set_ylim(0, ymax * (1.18 if annotate_counts else 1.08))
    ax.set_xlim(-0.5, counts.size - 0.5)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(panel_label, loc="left", fontsize=11, weight="bold")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{int(v):,}"))
    ax.grid(axis="y", linestyle=":", alpha=0.4)


def _draw_dividers(
    ax: Axes, groups: list[tuple[str, int, int]], *, label_groups: bool
) -> None:
    """Dotted divider between contiguous groups; group names on the top panel."""
    if len(groups) < 2:
        return
    for _, _, stop in groups[:-1]:
        ax.axvline(stop - 0.5, color=_INK, linestyle=":", linewidth=1.0, zorder=0)
    if label_groups:
        for name, start, stop in groups:
            ax.text(
                (start + stop - 1) / 2,
                1.0,
                name,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=8,
                color=_INK,
            )


# --------------------------------------------------------------------------- #
# Legend figure (rendered separately so it never overlaps the plot)
# --------------------------------------------------------------------------- #


def _legend_figure(
    bar_values: np.ndarray, color_map: dict[str, str], legend_title: str
) -> Figure:
    """Standalone swatch legend: one swatch per ``color_by`` value, first-seen order."""
    seen: list[str] = []
    for value in bar_values.astype(str):
        if value not in seen:
            seen.append(value)
    height = max(1.2, 0.3 * len(seen) + 0.7)
    fig, ax = plt.subplots(figsize=(2.6, height))
    ax.axis("off")
    handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markersize=12,
            markerfacecolor=color_map[value],
            markeredgecolor="none",
        )
        for value in seen
    ]
    ax.legend(
        handles,
        seen,
        title=legend_title or None,
        loc="center",
        frameon=False,
        fontsize=11,
        title_fontsize=12,
        ncol=1 if len(seen) <= 6 else 2,
    )
    return fig
