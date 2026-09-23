"""Coefficient-of-variation (CV) distribution figures for a :class:`Dataset`.

Project copy of the ``cv-plot`` template (lib/figures/cv.py v0.1), adapted for the
Stage-3 QC CV family of this project. Held to the correctness charter: assume nothing,
verify everything, fail loud.

What it draws (one figure): overlaid per-feature CV distributions -- one curve per
labelled :class:`Dataset` -- as translucent histograms plus Gaussian-KDE lines, with
each distribution's median CV marked by a dashed vertical line and annotated. Used here
to compare processing states of one dataset (raw vs median-normalized vs ComBat
batch-corrected).

CV is ``std / mean`` computed **per feature across the samples** on a **linear**
intensity scale (``ddof=1``). :func:`compute_cv` **raises** :class:`CVScaleError` on a
non-linear ``Dataset.scale`` (on log/centered scales ``std/mean`` is ill-defined).

Deviations from the template (all project wiring, no change to the computation):

* imports point at the project packages (``loaders.data_loading``,
  ``common.figures.colors``, ``common.figures.figure_io``);
* the project ``figure_io.save_figure`` must run inside ``publication_style()``, so
  :func:`save_cv` builds AND saves inside one style block;
* the median annotation stays terse (``median = 0.123``) and is colored like its
  state; the separate legend image names the colors (annotation budget);
* the title is drawn as the axes title (not a figure suptitle);
* :func:`plot_cv_distribution` accepts an optional ``xlabel`` so the caller can state
  the scale on the axis.

Control samples would be rendered separately from experimental samples; the caller
subsets the :class:`Dataset` before plotting.
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
from scipy import stats

__script_meta__: dict[str, object] = {
    "template": {"name": "cv-plot", "version": "0.1"},
    "kind": "module",
    "provides": [
        "CVScaleError",
        "CVResult",
        "CVPlot",
        "compute_cv",
        "plot_cv_distribution",
        "save_cv",
    ],
    "uses": [
        "loaders.data_loading",
        "common.figures.colors",
        "common.figures.figure_io",
    ],
    "seeded_from": {"template": "cv-plot", "version": "0.1"},
    "description": (
        "CV distribution figures from one or more Datasets: per-feature CV (std/mean) "
        "computed across samples on a linear scale, overlaid as histogram + "
        "Gaussian-KDE curves with annotated median lines. Hard-refuses a non-linear "
        "scale; state colors from the project registry with the >8-category guard; "
        "dual-export plus a separate legend image. Deviation from cv-plot@0.1: project "
        "import paths, save inside publication_style(), optional xlabel. Fail-loud."
    ),
}

# Color-registry namespace for methodological processing states, shared across every CV
# figure so e.g. "Median-normalized" keeps its color at protein and peptide level.
DEFAULT_CV_CATEGORY = "NormalizationState"

# Default x-axis upper bound: the 99.5th percentile of the pooled finite CVs.
_UPPER_PERCENTILE = 99.5


class CVScaleError(ValueError):
    """Raised when CV is requested on a non-linear abundance scale."""


@dataclass(frozen=True)
class CVResult:
    """The per-feature CVs underlying a CV figure.

    Attributes
    ----------
    cvs:
        ``{label: (n_features,) CV array}`` in input order; ``NaN`` marks a feature
        whose mean is <= 0 (undefined CV).
    medians:
        ``{label: median CV}`` over finite CVs (``NaN`` if none).
    """

    cvs: dict[str, np.ndarray]
    medians: dict[str, float]


@dataclass
class CVPlot:
    """A rendered CV figure plus its companion legend figure."""

    figure: Figure
    legend_figure: Figure
    result: CVResult
    color_map: dict[str, str]


# --------------------------------------------------------------------------- #
# Compute
# --------------------------------------------------------------------------- #


def compute_cv(dataset: Dataset) -> np.ndarray:
    """Per-feature coefficient of variation (``std / mean``, ``ddof=1``) across samples.

    A feature with mean <= 0 yields ``NaN``. Raises :class:`CVScaleError` if
    ``dataset.scale`` is not ``"linear"``; ``ValueError`` if abundances are not 2D or
    there are fewer than two samples.
    """
    if dataset.scale != "linear":
        raise CVScaleError(
            f"CV (std/mean) requires linear-scale abundances but the Dataset is on "
            f"scale {dataset.scale!r}. On log/centered scales the per-feature mean "
            f"sits near 0, so std/mean is ill-defined. Compute CV on linear data."
        )
    abundances = np.asarray(dataset.abundances, dtype=float)
    if abundances.ndim != 2:
        raise ValueError(
            f"abundances must be 2D (n_samples, n_features); got {abundances.shape}."
        )
    n_samples = abundances.shape[0]
    if n_samples < 2:
        raise ValueError(
            f"CV needs >= 2 samples to have a spread to measure; got {n_samples}."
        )

    means = np.nanmean(abundances, axis=0)
    stds = np.nanstd(abundances, axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cvs = stds / means
    cvs[means <= 0] = np.nan
    return np.asarray(cvs, dtype=float)


def _compute_result(datasets: Mapping[str, Dataset]) -> CVResult:
    """Compute each state's CV array + median (the scale guard fires per Dataset)."""
    cvs = {label: compute_cv(ds) for label, ds in datasets.items()}
    medians = {label: _nanmedian(cv) for label, cv in cvs.items()}
    return CVResult(cvs=cvs, medians=medians)


def _nanmedian(values: np.ndarray) -> float:
    """Median of the finite values (``NaN`` if none) without a runtime warning."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    return float(np.median(finite))


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #


def plot_cv_distribution(
    datasets: Mapping[str, Dataset],
    *,
    category: str = DEFAULT_CV_CATEGORY,
    feature_type: str = "feature",
    xlabel: str | None = None,
    max_cv: float | None = None,
    n_bins: int = 60,
    title: str | None = None,
    legend_title: str | None = None,
    show_medians: bool = True,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
) -> CVPlot:
    """Overlay the CV distributions of one or more labelled linear-scale Datasets.

    Must be called inside ``publication_style()`` if the style is to apply (see
    :func:`save_cv`, which does so). ``datasets`` insertion order is the draw/legend
    order; all must be linear and share identical ``feature_names`` (raised
    otherwise). ``xlabel`` overrides the default ``"<feature_type> CV (std / mean)"``.

    Raises
    ------
    CVScaleError
        If any Dataset is not on the ``"linear"`` scale.
    ValueError
        On an empty ``datasets``, invalid ``n_bins``/``max_cv``, or mismatched
        ``feature_names``.
    CategoricalPaletteExceededError
        If the labels exceed the palette's categorical capacity (the >8 guard).
    """
    if len(datasets) == 0:
        raise ValueError("datasets is empty; pass >= 1 labelled Dataset to plot.")
    if n_bins < 1:
        raise ValueError(f"n_bins must be >= 1; got {n_bins}.")
    if max_cv is not None and not max_cv > 0:
        raise ValueError(f"max_cv must be positive when given; got {max_cv}.")

    _validate_feature_alignment(datasets)
    result = _compute_result(datasets)
    labels = list(datasets.keys())
    title_for_legend = legend_title if legend_title is not None else category
    axis_label = xlabel if xlabel is not None else f"{feature_type} CV (std / mean)"

    fig, ax = plt.subplots(figsize=(10, 5.5))
    try:
        color_map = assign_colors(
            category,
            labels,
            registry_path=registry_path,
            persist=persist_colors,
        )
        bin_upper = _upper_bound(result, labels, max_cv)
        _draw_distributions(
            ax,
            result,
            labels,
            color_map,
            bin_upper=bin_upper,
            n_bins=n_bins,
            xlabel=axis_label,
        )
        if show_medians:
            _annotate_medians(ax, result, labels, color_map)
        if title is not None:
            ax.set_title(title, fontsize=14, weight="bold")
        legend_figure = _legend_figure(labels, color_map, title_for_legend)
    except BaseException:
        plt.close(fig)
        raise

    return CVPlot(
        figure=fig, legend_figure=legend_figure, result=result, color_map=color_map
    )


def save_cv(
    datasets: Mapping[str, Dataset],
    output_dir: str | Path,
    base_name: str,
    *,
    category: str = DEFAULT_CV_CATEGORY,
    feature_type: str = "feature",
    xlabel: str | None = None,
    max_cv: float | None = None,
    n_bins: int = 60,
    title: str | None = None,
    legend_title: str | None = None,
    show_medians: bool = True,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
    dpi: int = 300,
) -> tuple[FigureArtifacts, CVPlot]:
    """Render a CV figure and save it dual-export (``<base>.{svg,png}`` + legend).

    Build and save both happen inside one ``publication_style()`` block (required by
    the project ``save_figure``). Returns the artifact paths and the (now closed)
    :class:`CVPlot`, whose ``result``/``color_map`` the caller records.
    """
    with publication_style():
        plot = plot_cv_distribution(
            datasets,
            category=category,
            feature_type=feature_type,
            xlabel=xlabel,
            max_cv=max_cv,
            n_bins=n_bins,
            title=title,
            legend_title=legend_title,
            show_medians=show_medians,
            registry_path=registry_path,
            persist_colors=persist_colors,
        )
        artifacts = save_figure(
            plot.figure,
            output_dir,
            base_name,
            legend_fig=plot.legend_figure,
            dpi=dpi,
        )
    return artifacts, plot


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #


def _validate_feature_alignment(datasets: Mapping[str, Dataset]) -> None:
    """Refuse overlaying CVs computed over different feature sets (fail loud)."""
    reference_names: np.ndarray | None = None
    reference_label = ""
    for label, ds in datasets.items():
        names = np.asarray(ds.feature_names)
        if reference_names is None:
            reference_names = names
            reference_label = label
            continue
        if names.shape != reference_names.shape or not np.array_equal(
            names, reference_names
        ):
            raise ValueError(
                f"Dataset {label!r} has different feature_names than "
                f"{reference_label!r}; CV distributions are only comparable across the "
                f"same features in the same order. Align them before plotting."
            )


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #


def _upper_bound(result: CVResult, labels: list[str], max_cv: float | None) -> float:
    """The x-axis / bin / KDE upper bound: ``max_cv`` if given, else the 99.5th pct."""
    if max_cv is not None:
        return float(max_cv)
    pooled = _pooled_finite(result, labels)
    if pooled.size == 0:
        return 1.0
    upper = float(np.percentile(pooled, _UPPER_PERCENTILE))
    return upper if upper > 0 else 1.0


def _pooled_finite(result: CVResult, labels: list[str]) -> np.ndarray:
    """All finite CVs across every state, concatenated (for shared bins/limits)."""
    parts = [result.cvs[label][np.isfinite(result.cvs[label])] for label in labels]
    return np.concatenate(parts) if parts else np.asarray([], dtype=float)


def _draw_distributions(
    ax: Axes,
    result: CVResult,
    labels: list[str],
    color_map: dict[str, str],
    *,
    bin_upper: float,
    n_bins: int,
    xlabel: str,
) -> None:
    """Histogram + KDE for each state on shared bins; set the axis labels and limits."""
    bins = [float(edge) for edge in np.linspace(0.0, bin_upper, n_bins + 1)]
    for label in labels:
        finite = result.cvs[label][np.isfinite(result.cvs[label])]
        color = color_map[label]
        ax.hist(
            finite,
            bins=bins,
            density=True,
            alpha=0.30,
            color=color,
            label=label,
            edgecolor="none",
        )
        _overlay_kde(ax, finite, color=color, x_upper=bin_upper)

    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Density", fontsize=14)
    ax.tick_params(labelsize=12)
    ax.set_xlim(0.0, bin_upper)


def _overlay_kde(ax: Axes, finite: np.ndarray, *, color: str, x_upper: float) -> None:
    """Overlay a Gaussian-KDE curve for the finite CVs (skip if degenerate)."""
    if finite.size < 2 or float(finite.std()) <= 0.0:
        return
    lower = float(finite.min())
    if not lower < x_upper:
        return
    kde = stats.gaussian_kde(finite)
    x_grid = np.linspace(lower, x_upper, 300)
    density = np.asarray(kde(x_grid), dtype=float)
    ax.plot(x_grid, density, color=color, linewidth=2.25)


def _annotate_medians(
    ax: Axes, result: CVResult, labels: list[str], color_map: dict[str, str]
) -> None:
    """Dashed median line + a terse text label for each state (finite medians only)."""
    for label in labels:
        med = result.medians[label]
        if np.isfinite(med):
            ax.axvline(med, color=color_map[label], linestyle="--", linewidth=1.4)

    # Reserve headroom above the data so the stacked labels sit clear of the curves,
    # and back each label with white so neighbouring dashed lines pass behind it.
    y_data_top = ax.get_ylim()[1]
    n_labels = len(labels)
    y_top = y_data_top * (1.0 + 0.09 * n_labels + 0.04)
    ax.set_ylim(0.0, y_top)
    for i, label in enumerate(labels):
        med = result.medians[label]
        if not np.isfinite(med):
            continue
        ax.text(
            med,
            y_top * (0.985 - 0.075 * i),
            f"median = {med:.3f}",
            color=color_map[label],
            fontsize=12,
            weight="bold",
            va="top",
            ha="left",
            zorder=5,
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.9,
                "pad": 1.5,
            },
        )


# --------------------------------------------------------------------------- #
# Legend figure (rendered separately so it never overlaps the plot)
# --------------------------------------------------------------------------- #


def _legend_figure(
    labels: list[str], color_map: dict[str, str], legend_title: str
) -> Figure:
    """Standalone swatch legend: one colored line per state, in ``labels`` order."""
    height = max(1.4, 0.35 * len(labels) + 0.8)
    fig, ax = plt.subplots(figsize=(4.2, height))
    ax.axis("off")
    handles = [
        Line2D([0], [0], color=color_map[label], linewidth=2.5) for label in labels
    ]
    ax.legend(
        handles,
        labels,
        title=legend_title,
        loc="center",
        frameon=True,
        fontsize=12,
        title_fontsize=13,
    )
    return fig
