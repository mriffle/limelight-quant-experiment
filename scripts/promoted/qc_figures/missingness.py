"""Data-completeness / missingness diagnostic figures for a :class:`Dataset`.

PROJECT COPY seeded from the plugin template ``lib/figures/missingness.py``
(``missingness`` v0.1). Held to the correctness charter: assume nothing, verify
everything, fail loud.

What it draws (one figure, two panels):

* **Completeness curve** (left) -- features retained vs the required *detection
  fraction* t: how many features are detected in at least a fraction t of the runs in
  a group. Overlaid per ``color_by`` group (here: acquisition batch) plus a neutral
  whole-cohort curve (all runs), which is the feature-filter-relevant one.
* **MNAR diagnostic** (right) -- per-feature detection rate (over all runs) vs mean
  log2 abundance of its detected values, as a hexbin with a binned-median trend and
  the Pearson r annotated.

Detection: a feature is *detected* in a sample when its value is finite and strictly
greater than ``min_intensity`` (default ``0.0``). Hard-refuses a non-linear scale.

Deviations from the template (all deliberate):

1. Package layout: ``loaders.data_loading`` / ``common.figures.*`` imports.
2. Colors are read **read-only** from the registry (``persist`` never writes): a
   ``color_by`` value missing from the registry raises instead of being assigned a new
   color (sibling figure jobs run concurrently against the same registry file, and the
   registry has no cross-process lock).
3. Completeness step alignment fixed: the template drew ``where="post"``, which shows
   the count for ">= k of m" over thresholds ``[k/m, (k+1)/m)`` -- one step too
   permissive (e.g. at t = 0.75 with m = 2 it showed the ">= 1 of 2" count). Here the
   curve is drawn ``where="pre"`` from t = 0, so for every t in ``((k-1)/m, k/m]`` the
   plotted value is the number of features detected in ``>= ceil(t*m)`` runs, with a
   marker at each attainable threshold k/m.
4. A neutral whole-cohort ("all runs") curve is drawn in addition to the per-group
   curves, outside the categorical palette (no registry slot consumed).
5. Annotation budget: the on-canvas interpretation text ("low abundance -> more
   missing = MNAR") is removed; the MNAR panel carries only ``Pearson r`` and the
   feature count. The binned-median key moved from the axes into the separate legend
   image, so the only on-figure key is the hexbin density colorbar (the data's own
   count scale, as in the template).
6. ``summarize_completeness`` returns the load-bearing numbers (features retained at
   100% / >= 50% detection) for provenance and the finding.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from common.figures.colors import DEFAULT_REGISTRY_PATH, load_registry
from common.figures.figure_io import FigureArtifacts, publication_style, save_figure
from loaders.data_loading import Dataset
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

__script_meta__: dict[str, object] = {
    "template": {"name": "missingness", "version": "0.1"},
    "kind": "module",
    "provides": [
        "MissingnessScaleError",
        "CompletenessResult",
        "MissingnessPlot",
        "compute_completeness",
        "summarize_completeness",
        "plot_missingness",
        "save_missingness",
    ],
    "uses": [
        "loaders.data_loading",
        "common.figures.colors",
        "common.figures.figure_io",
    ],
    "seeded_from": {"template": "missingness", "version": "0.1"},
    "description": (
        "Data-completeness / missingness diagnostics from a linear-scale Dataset: a "
        "two-panel figure of (1) a feature completeness curve (features retained vs "
        "required detection fraction) overlaid per group (read-only registry colors, "
        ">8 guard) plus an all-runs curve and (2) an MNAR diagnostic (per-feature "
        "detection rate vs mean log2 abundance, hexbin + binned median + Pearson r). "
        "Detected = finite and > min_intensity; hard-refuses a non-linear scale. "
        "Dual-export plus a separate legend image."
    ),
}

# Whole-cohort curve color: a neutral dark gray outside the Okabe-Ito palette, so it
# consumes no registry slot and cannot collide with a group color.
ALL_RUNS_COLOR = "#4D4D4D"

# Hexbin colormap for the MNAR panel (perceptually uniform; density on a log scale).
_HEXBIN_CMAP = "viridis"

# Trend-line color for the MNAR binned median (Okabe-Ito vermillion, off the hexbin).
_TREND_COLOR = "#D55E00"

# Number of abundance bins for the MNAR binned-median trend.
_TREND_BINS = 20

# Maximum categorical colors (the >8 rule); re-checked against the registry palette.
_MAX_CATEGORICAL = 8


class MissingnessScaleError(ValueError):
    """Raised when completeness/missingness is requested on a non-linear scale."""


class RegistryColorMissingError(KeyError):
    """Raised when a ``color_by`` value has no color in the registry (read-only use)."""


@dataclass(frozen=True)
class CompletenessResult:
    """Per-feature / per-sample detection summaries underlying a missingness figure.

    Attributes
    ----------
    sample_ids:
        ``(n_samples,)`` sample identifiers (row order of the Dataset).
    feature_detection_rate:
        ``(n_features,)`` fraction of all samples in which each feature is detected.
    sample_detection_rate:
        ``(n_samples,)`` fraction of features each sample detects (not drawn).
    feature_mean_log_abundance:
        ``(n_features,)`` mean log2 of each feature's positive detected values
        (``NaN`` for a feature never detected).
    mnar_correlation:
        Pearson r between mean log2 abundance and detection rate over features where
        both are finite (``NaN`` if undefined).
    n_mnar_features:
        Number of features entering the Pearson r / hexbin.
    """

    sample_ids: np.ndarray
    feature_detection_rate: np.ndarray
    sample_detection_rate: np.ndarray
    feature_mean_log_abundance: np.ndarray
    mnar_correlation: float
    n_mnar_features: int


@dataclass
class MissingnessPlot:
    """A rendered missingness figure plus its companion legend figure."""

    figure: Figure
    legend_figure: Figure
    result: CompletenessResult
    color_map: dict[str, str]


# --------------------------------------------------------------------------- #
# Compute
# --------------------------------------------------------------------------- #


def compute_completeness(
    dataset: Dataset, *, min_intensity: float = 0.0
) -> CompletenessResult:
    """Per-feature / per-sample detection rates + the MNAR correlation.

    Raises :class:`MissingnessScaleError` unless ``dataset.scale == "linear"``, or
    ``ValueError`` on a non-2D / empty / non-row-aligned Dataset.
    """
    mask, abundances = _detection_mask(dataset, min_intensity)
    sample_ids = _sample_ids(dataset, mask.shape[0])
    return _build_result(mask, abundances, sample_ids)


def features_retained(mask: np.ndarray, fraction: float) -> int:
    """Number of features detected in ``>= ceil(fraction * n_samples)`` samples.

    ``fraction`` must lie in ``(0, 1]``. Exposed so the reported key numbers are
    computed by the same rule as the drawn curve.
    """
    if not 0.0 < fraction <= 1.0:
        raise ValueError(f"fraction must be in (0, 1]; got {fraction}.")
    n_samples = mask.shape[0]
    k = math.ceil(fraction * n_samples - 1e-12)
    return int((mask.sum(axis=0) >= k).sum())


def summarize_completeness(
    dataset: Dataset,
    *,
    fractions: tuple[float, ...] = (1.0, 0.5),
    group_by: str | None = None,
    min_intensity: float = 0.0,
) -> dict[str, dict[str, int]]:
    """Features retained at each ``fraction`` for all runs and (optionally) per group.

    Returns ``{"all runs": {"1.0": n, "0.5": n}, "<group>": {...}, ...}``.
    """
    mask, _ = _detection_mask(dataset, min_intensity)
    out: dict[str, dict[str, int]] = {
        "all runs": {str(f): features_retained(mask, f) for f in fractions}
    }
    classes = _resolve_color_by(dataset, group_by)
    if classes is not None:
        for cls in _ordered_unique(classes):
            sub = mask[classes == cls]
            out[cls] = {str(f): features_retained(sub, f) for f in fractions}
    return out


def _detection_mask(
    dataset: Dataset, min_intensity: float
) -> tuple[np.ndarray, np.ndarray]:
    """The boolean (n_samples, n_features) detected mask + the float abundances."""
    if dataset.scale != "linear":
        raise MissingnessScaleError(
            f"Completeness/missingness requires linear-scale abundances but the "
            f"Dataset is on scale {dataset.scale!r}. On log/centered scales a 0 is an "
            f"ordinary value, so '> {min_intensity}' no longer means 'detected'."
        )
    abundances = np.asarray(dataset.abundances, dtype=float)
    if abundances.ndim != 2:
        raise ValueError(
            f"abundances must be 2D (n_samples, n_features); got {abundances.shape}."
        )
    n_samples, n_features = abundances.shape
    if n_samples < 1 or n_features < 1:
        raise ValueError(
            f"Dataset is empty ({n_samples} samples x {n_features} features)."
        )
    mask = np.isfinite(abundances) & (abundances > min_intensity)
    return mask, abundances


def _sample_ids(dataset: Dataset, n_samples: int) -> np.ndarray:
    """The row-aligned sample ids (fail loud if metadata is not row-aligned)."""
    ids = np.asarray(dataset.metadata.index.to_numpy(), dtype=str)
    if ids.shape[0] != n_samples:
        raise ValueError(
            f"metadata has {ids.shape[0]} rows but {n_samples} abundance samples; it "
            f"is not row-aligned."
        )
    return ids


def _build_result(
    mask: np.ndarray, abundances: np.ndarray, sample_ids: np.ndarray
) -> CompletenessResult:
    """Assemble the detection rates + MNAR correlation from the detected mask."""
    feature_detection_rate = np.asarray(mask.mean(axis=0), dtype=float)
    sample_detection_rate = np.asarray(mask.mean(axis=1), dtype=float)
    mean_log_ab = _feature_mean_log_abundance(mask, abundances)
    ok = np.isfinite(mean_log_ab) & np.isfinite(feature_detection_rate)
    mnar = _mnar_correlation(mean_log_ab[ok], feature_detection_rate[ok])
    return CompletenessResult(
        sample_ids=sample_ids,
        feature_detection_rate=feature_detection_rate,
        sample_detection_rate=sample_detection_rate,
        feature_mean_log_abundance=mean_log_ab,
        mnar_correlation=mnar,
        n_mnar_features=int(ok.sum()),
    )


def _feature_mean_log_abundance(mask: np.ndarray, abundances: np.ndarray) -> np.ndarray:
    """Mean log2 of each feature's positive detected values (``NaN`` if it has none)."""
    log_mask = mask & (abundances > 0.0)
    counts = log_mask.sum(axis=0).astype(float)
    safe = np.where(log_mask, abundances, 1.0)  # placeholder -> log2 == 0, then zeroed
    log_sum = np.where(log_mask, np.log2(safe), 0.0).sum(axis=0)
    out = np.full(abundances.shape[1], np.nan, dtype=float)
    np.divide(log_sum, counts, out=out, where=counts > 0)
    return out


def _mnar_correlation(a: np.ndarray, b: np.ndarray) -> float:
    """Pearson r between two finite vectors (NaN if undefined)."""
    if a.size < 2 or float(a.std()) == 0.0 or float(b.std()) == 0.0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #


def plot_missingness(
    dataset: Dataset,
    *,
    color_by: str,
    min_intensity: float = 0.0,
    title: str | None = None,
    legend_title: str | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
) -> MissingnessPlot:
    """Draw the two-panel completeness + MNAR diagnostic for one Dataset.

    Parameters
    ----------
    dataset:
        A single linear-scale :class:`Dataset` (raised otherwise).
    color_by:
        Metadata column splitting the completeness curve into one registry-colored
        line per group (read-only registry lookup; <= 8 groups).
    min_intensity:
        Detection threshold (detected == finite and ``> min_intensity``).
    title:
        Optional figure suptitle.
    legend_title:
        Title of the companion legend figure; defaults to ``color_by``.
    registry_path:
        Color registry JSON (read only).

    Raises
    ------
    MissingnessScaleError, ValueError, RegistryColorMissingError
    """
    mask, abundances = _detection_mask(dataset, min_intensity)
    sample_ids = _sample_ids(dataset, mask.shape[0])
    result = _build_result(mask, abundances, sample_ids)
    classes = _resolve_color_by(dataset, color_by)
    if classes is None:  # pragma: no cover - color_by is required (typed str)
        raise ValueError("color_by is required.")
    ordered = _ordered_unique(classes)
    color_map = _registry_colors(color_by, ordered, registry_path)

    fig, (ax_curve, ax_mnar) = plt.subplots(
        1, 2, figsize=(13, 5.2), constrained_layout=True
    )
    try:
        _draw_completeness(ax_curve, mask, classes, ordered, color_map)
        _draw_mnar(ax_mnar, fig, result)
        if title is not None:
            fig.suptitle(title, fontsize=14, weight="bold")
        legend_figure = _legend_figure(
            ordered,
            classes,
            color_map,
            legend_title if legend_title is not None else color_by,
        )
    except BaseException:
        plt.close(fig)
        raise
    return MissingnessPlot(
        figure=fig, legend_figure=legend_figure, result=result, color_map=color_map
    )


def save_missingness(
    dataset: Dataset,
    output_dir: str | Path,
    base_name: str,
    *,
    color_by: str,
    min_intensity: float = 0.0,
    title: str | None = None,
    legend_title: str | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    dpi: int = 300,
) -> tuple[FigureArtifacts, MissingnessPlot]:
    """Render (:func:`plot_missingness`) and dual-export the figure + legend image.

    Built and saved inside :func:`publication_style` (the style must be active while
    artists are created and at save time). Figures are closed by ``save_figure``.
    """
    with publication_style():
        plot = plot_missingness(
            dataset,
            color_by=color_by,
            min_intensity=min_intensity,
            title=title,
            legend_title=legend_title,
            registry_path=registry_path,
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
# Validation / colors
# --------------------------------------------------------------------------- #


def _resolve_color_by(dataset: Dataset, color_by: str | None) -> np.ndarray | None:
    """Validate + extract ``color_by`` as per-sample strings (or ``None``)."""
    if color_by is None:
        return None
    if color_by not in dataset.metadata.columns:
        raise ValueError(
            f"color_by {color_by!r} is not a metadata column "
            f"{list(dataset.metadata.columns)}."
        )
    series = dataset.metadata[color_by]
    if bool(series.isna().any()):
        raise ValueError(
            f"color_by {color_by!r} has {int(series.isna().sum())} missing value(s)."
        )
    return np.asarray(series.to_numpy(), dtype=str)


def _ordered_unique(values: np.ndarray) -> list[str]:
    """Distinct values in first-seen order."""
    seen: list[str] = []
    for value in values.astype(str):
        if value not in seen:
            seen.append(value)
    return seen


def _registry_colors(
    category: str, values: list[str], registry_path: str | Path
) -> dict[str, str]:
    """Read-only registry lookup; raise on a missing value or > 8 categories."""
    if len(values) > _MAX_CATEGORICAL:
        raise ValueError(
            f"{category!r} has {len(values)} groups; color encodes at most "
            f"{_MAX_CATEGORICAL}. Facet or use a second channel instead."
        )
    registry = load_registry(registry_path)
    entry = registry.get(category)
    if not isinstance(entry, dict) or not isinstance(entry.get("values"), dict):
        raise RegistryColorMissingError(
            f"Registry {registry_path} has no category {category!r}."
        )
    known: dict[str, str] = {str(k): str(v) for k, v in entry["values"].items()}
    missing = [v for v in values if v not in known]
    if missing:
        raise RegistryColorMissingError(
            f"Registry category {category!r} has no color for {missing}; add them to "
            f"the registry first (this module never invents colors)."
        )
    return {v: known[v] for v in values}


# --------------------------------------------------------------------------- #
# Drawing -- completeness curve
# --------------------------------------------------------------------------- #


def _completeness_curve(counts: np.ndarray, m: int) -> tuple[np.ndarray, np.ndarray]:
    """(t, retained) for t = 0, 1/m, ..., 1: features detected in >= ceil(t*m) of m.

    The t = 0 point repeats the k = 1 count so a ``where="pre"`` step covers
    ``(0, 1/m]`` with the ">= 1 of m" value.
    """
    ks = np.arange(1, m + 1)
    retained = np.array([int((counts >= k).sum()) for k in ks], dtype=int)
    t = np.concatenate([[0.0], ks / m])
    y = np.concatenate([[retained[0]], retained])
    return t, y


def _draw_curve(
    ax: Axes, counts: np.ndarray, m: int, color: str, *, linestyle: str = "-"
) -> int:
    """One step curve + markers at attainable thresholds; returns its minimum."""
    t, y = _completeness_curve(counts, m)
    ax.step(t, y, where="pre", color=color, linewidth=2.2, linestyle=linestyle)
    ax.plot(t[1:], y[1:], linestyle="none", marker="o", markersize=5, color=color)
    return int(y.min())


def _draw_completeness(
    ax: Axes,
    mask: np.ndarray,
    classes: np.ndarray,
    ordered: list[str],
    color_map: dict[str, str],
) -> None:
    """Per-group curves (registry colors) + the neutral all-runs curve."""
    n_samples, n_features = mask.shape
    minima: list[int] = []
    for cls in ordered:
        members = classes == cls
        m = int(members.sum())
        minima.append(_draw_curve(ax, mask[members].sum(axis=0), m, color_map[cls]))
    minima.append(
        _draw_curve(ax, mask.sum(axis=0), n_samples, ALL_RUNS_COLOR, linestyle="--")
    )
    n_never = int((mask.sum(axis=0) == 0).sum())
    ax.axhline(n_features, color="#999999", linestyle=":", linewidth=0.9)
    ax.text(
        0.01,
        n_features,
        f"all features: {n_features:,} ({n_never:,} in 0/{n_samples} runs)",
        va="bottom",
        ha="left",
        fontsize=9,
        color="#555555",
    )
    low = min(minima)
    pad = 0.04 * (n_features - low) if n_features > low else 1.0
    ax.set_ylim(max(0.0, low - pad), n_features + 0.06 * (n_features - low + 1))
    ax.set_xlim(0.0, 1.02)
    ax.set_xlabel("required detection fraction t (≥ t of runs in group)")
    ax.set_ylabel("features retained (count)")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{v:,.0f}"))
    ax.set_title("Completeness by batch", loc="left", fontsize=12, weight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)


# --------------------------------------------------------------------------- #
# Drawing -- MNAR diagnostic
# --------------------------------------------------------------------------- #


def _draw_mnar(ax: Axes, fig: Figure, result: CompletenessResult) -> None:
    """Hexbin of detection rate vs mean log2 abundance + binned-median trend + r."""
    mean_log_ab = result.feature_mean_log_abundance
    rate = result.feature_detection_rate
    ok = np.isfinite(mean_log_ab) & np.isfinite(rate)
    x = mean_log_ab[ok]
    y = rate[ok]
    n_samples = result.sample_ids.shape[0]

    # Detection rates are discrete (k / n_samples). Pin the hex lattice so its primary
    # rows sit exactly on k / n_samples (extent y = [0, 1], ny = n_samples); every
    # point then bins into a primary-row hex (never the half-offset secondary rows).
    hexbin = ax.hexbin(
        x,
        y,
        gridsize=(40, n_samples),
        extent=(float(x.min()), float(x.max()), 0.0, 1.0),
        cmap=_HEXBIN_CMAP,
        bins="log",
        mincnt=1,
    )
    cbar = fig.colorbar(hexbin, ax=ax)
    cbar.set_label("features per bin (log scale)")

    centers, medians = _binned_median(x, y, _TREND_BINS)
    if centers.size:
        ax.plot(
            centers,
            medians,
            color=_TREND_COLOR,
            linewidth=2.4,
            marker="o",
            markersize=4,
        )

    annotation = "Pearson r undefined (no spread)"
    if np.isfinite(result.mnar_correlation):
        annotation = (
            f"Pearson r = {result.mnar_correlation:.2f}\n"
            f"n = {result.n_mnar_features:,} features"
        )
    ax.text(
        0.97,
        0.04,
        annotation,
        transform=ax.transAxes,
        fontsize=11,
        va="bottom",
        ha="right",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "#cccccc"},
    )
    ax.set_ylim(-0.02, 1.06)
    ax.set_yticks(np.arange(0, n_samples + 1) / n_samples)
    ax.set_yticklabels([f"{k}/{n_samples}" for k in range(n_samples + 1)])
    ax.set_xlabel("mean log2 intensity of detected values")
    ax.set_ylabel("detection rate (runs detected / all runs)")
    ax.set_title("MNAR diagnostic (all runs)", loc="left", fontsize=12, weight="bold")


def _binned_median(
    x: np.ndarray, y: np.ndarray, n_bins: int
) -> tuple[np.ndarray, np.ndarray]:
    """Median ``y`` within ``n_bins`` percentile bins of ``x`` (skips empty bins).

    Bins are half-open ``[lo, hi)`` except the last, which is closed, so a value on
    an interior edge is counted once.
    """
    if x.size < 2:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    edges = np.unique(np.percentile(x, np.linspace(0.0, 100.0, n_bins + 1)))
    centers: list[float] = []
    medians: list[float] = []
    last = len(edges) - 2
    for i, (lo, hi) in enumerate(pairwise(edges)):
        in_bin = (x >= lo) & ((x <= hi) if i == last else (x < hi))
        if bool(in_bin.any()):
            centers.append(float(np.median(x[in_bin])))
            medians.append(float(np.median(y[in_bin])))
    return np.asarray(centers, dtype=float), np.asarray(medians, dtype=float)


# --------------------------------------------------------------------------- #
# Legend figure (rendered separately so it never overlaps the plot)
# --------------------------------------------------------------------------- #


def _legend_figure(
    ordered: list[str],
    classes: np.ndarray,
    color_map: dict[str, str],
    legend_title: str,
) -> Figure:
    """Standalone key: group curves (with n), the all-runs curve, the MNAR trend."""
    counts = {cls: int((classes == cls).sum()) for cls in ordered}
    handles: list[Line2D] = [
        Line2D([0], [0], color=color_map[cls], linewidth=2.5, marker="o")
        for cls in ordered
    ]
    labels = [f"{cls} (n={counts[cls]} runs)" for cls in ordered]
    handles.append(
        Line2D(
            [0], [0], color=ALL_RUNS_COLOR, linewidth=2.5, linestyle="--", marker="o"
        )
    )
    labels.append(f"all runs (n={classes.shape[0]})")
    handles.append(
        Line2D([0], [0], color=_TREND_COLOR, linewidth=2.5, marker="o", markersize=4)
    )
    labels.append("binned median (MNAR panel)")
    fig, ax = plt.subplots(figsize=(3.4, 0.3 * len(labels) + 0.5))
    ax.axis("off")
    ax.legend(
        handles,
        labels,
        title=legend_title or None,
        loc="center",
        frameon=True,
        fontsize=10,
        title_fontsize=11,
        handlelength=3.2,
    )
    return fig
