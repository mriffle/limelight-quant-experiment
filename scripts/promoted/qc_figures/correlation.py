"""Sample-vs-sample correlation heatmap for a :class:`Dataset` (a QC sanity check).

PROJECT COPY seeded from the plugin template ``lib/figures/correlation.py``
(``sample-correlation`` v0.1). Held to the correctness charter: assume nothing,
verify everything, fail loud.

What it draws (one figure): a hierarchically-clustered **sample x sample** Pearson /
Spearman correlation heatmap over whole feature profiles, a top dendrogram, and one
registry-colored **annotation stripe** per categorical design column, so one can read
off whether the clustering lines up with the design. A value-scale colorbar sits
beside the heatmap (the documented on-axes exception); the stripe keys are a
**separate legend figure** (``<base>.legend.{svg,png}``).

Project adaptations vs the template (all deliberate, all small-n driven):

* **Sized for a handful of samples** (this study: 8). The template's 16 x 14 in canvas
  and 5 pt tick labels target hundreds of samples; here the canvas is ~7.5 x 7.5 in
  with legible ticks, and every cell is **annotated with its r** (``annotate_cells``)
  when n is small enough for the text to fit (``_MAX_ANNOTATED``).
* **Diagonal masked.** Self-correlation is trivially 1 and carries no information; it
  is drawn in neutral gray (unannotated) and the color scale spans the **off-diagonal
  min..max**, so the whole colormap resolves the informative comparisons. (The
  template anchored vmax at 1.0 with the diagonal colored.)
* **Registry is read-only here.** Every annotation value must already be registered
  in ``state/color_registry.json``; an unregistered value raises rather than being
  assigned a color at plot time (concurrent figure jobs must not race to mutate the
  shared registry, and the design colors were fixed at Stage 1).
* Imports follow the project layout: ``loaders.data_loading`` and
  ``common.figures.{colors,figure_io}``.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common.figures.colors import DEFAULT_REGISTRY_PATH, load_registry
from common.figures.figure_io import FigureArtifacts, publication_style, save_figure
from loaders.data_loading import LOG_SCALES, Dataset
from matplotlib.axes import Axes
from matplotlib.colors import Normalize, to_rgba
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram, linkage

__script_meta__: dict[str, object] = {
    "task": "stage3-qc-sample-correlation",
    "kind": "module",
    "provides": [
        "CorrelationScaleWarning",
        "CorrelationResult",
        "CorrelationPlot",
        "compute_correlation",
        "plot_sample_correlation",
        "save_sample_correlation",
    ],
    "uses": [
        "loaders.data_loading",
        "common.figures.colors",
        "common.figures.figure_io",
    ],
    "seeded_from": {"template": "sample-correlation", "version": "0.1"},
    "description": (
        "Small-n sample-vs-sample correlation heatmap: average-linkage clustered "
        "Pearson/Spearman correlation of whole sample profiles with a dendrogram, "
        "per-cell r annotation, masked diagonal (color scale = off-diagonal range), "
        "read-only registry-colored categorical annotation stripes, an on-axes "
        "r-colorbar, and a separate stripe-key legend image; dual export. Warns on "
        "Pearson over a non-log scale. Fail-loud on non-finite / constant samples, "
        "duplicate ids, missing annotation values, or unregistered colors."
    ),
}

_METHODS: frozenset[str] = frozenset({"pearson", "spearman"})
_METHOD_SYMBOL: dict[str, str] = {"pearson": "Pearson r", "spearman": "Spearman rho"}

_HEATMAP_CMAP = "RdYlBu_r"
_DIAGONAL_COLOR = "#d9d9d9"

_LINKAGE_METRIC = "euclidean"
_LINKAGE_METHOD = "average"

# Above this many samples the per-cell r text no longer fits legibly.
_MAX_ANNOTATED = 16

# Max categorical values per stripe (the registry's eight-colour budget).
_MAX_CATEGORIES = 8


class CorrelationScaleWarning(UserWarning):
    """Pearson correlation is running on a non-log (e.g. linear) scale."""


@dataclass(frozen=True)
class _Annotation:
    """One resolved categorical annotation dimension (a stripe above the heatmap)."""

    column: str
    label: str
    per_sample: np.ndarray  # str values, input sample order
    colors: dict[str, str]  # value -> hex, from the registry


@dataclass(frozen=True)
class CorrelationResult:
    """The sample-vs-sample correlation underlying the figure.

    Attributes
    ----------
    matrix:
        ``(n, n)`` correlation matrix in the **input** sample order. Symmetric, unit
        diagonal.
    sample_ids:
        ``(n,)`` unique sample identifiers (row/column order of ``matrix``).
    method:
        ``"pearson"`` or ``"spearman"``.
    order:
        ``(n,)`` clustered display order (dendrogram leaves).
    linkage_z:
        The scipy linkage matrix the order came from.
    n_features:
        Number of features correlated over.
    """

    matrix: np.ndarray
    sample_ids: np.ndarray
    method: str
    order: np.ndarray
    linkage_z: np.ndarray
    n_features: int

    def off_diagonal(self) -> np.ndarray:
        """The strictly-upper-triangle correlations (each pair once)."""
        iu = np.triu_indices(len(self.matrix), k=1)
        values: np.ndarray = self.matrix[iu]
        return values

    def ordered_frame(self) -> pd.DataFrame:
        """The matrix reordered by clustering, as a labelled frame."""
        ids = self.sample_ids[self.order]
        return pd.DataFrame(
            self.matrix[np.ix_(self.order, self.order)],
            index=pd.Index(ids, name="sample_id"),
            columns=pd.Index(ids),
        )


@dataclass
class CorrelationPlot:
    """A rendered correlation figure plus its companion stripe-key legend figure."""

    figure: Figure
    legend_figure: Figure | None
    result: CorrelationResult
    color_maps: dict[str, dict[str, str]]


# --------------------------------------------------------------------------- #
# Compute
# --------------------------------------------------------------------------- #


def compute_correlation(
    dataset: Dataset,
    *,
    method: str = "pearson",
    linkage_metric: str = _LINKAGE_METRIC,
    linkage_method: str = _LINKAGE_METHOD,
) -> CorrelationResult:
    """Correlate every pair of samples over their full feature profiles and cluster.

    Raises
    ------
    ValueError
        Unknown ``method``; non-2D / non-finite abundances; < 3 samples (clustering
        needs at least three leaves to be informative); non-unique or misaligned
        sample ids; a zero-variance sample.
    """
    if method not in _METHODS:
        raise ValueError(f"method {method!r} not supported; use {sorted(_METHODS)}.")
    abundances = np.asarray(dataset.abundances, dtype=float)
    if abundances.ndim != 2:
        raise ValueError(f"abundances must be 2D; got shape {abundances.shape}.")
    n_samples, n_features = abundances.shape
    if n_samples < 3:
        raise ValueError(f"need >= 3 samples to cluster; got {n_samples}.")
    if n_features < 3:
        raise ValueError(f"need >= 3 features to correlate; got {n_features}.")
    if not np.isfinite(abundances).all():
        n_bad = int((~np.isfinite(abundances)).sum())
        raise ValueError(
            f"compute_correlation requires finite abundances; found {n_bad} "
            f"non-finite value(s). Resolve missingness first."
        )
    sample_ids = dataset.metadata.index.to_numpy().astype(str)
    if len(sample_ids) != n_samples:
        raise ValueError(
            f"metadata has {len(sample_ids)} rows but abundances has {n_samples} "
            f"samples; the Dataset is not row-aligned."
        )
    dup = pd.Index(sample_ids).duplicated(keep=False)
    if bool(dup.any()):
        raise ValueError(f"duplicate sample ids: {sorted(set(sample_ids[dup]))}.")
    constant = np.std(abundances, axis=1) == 0.0
    if bool(constant.any()):
        raise ValueError(
            f"sample(s) constant across all features (correlation undefined): "
            f"{sample_ids[constant].tolist()}."
        )

    frame = pd.DataFrame(abundances.T, columns=pd.Index(sample_ids))
    corr_method = cast(Literal["pearson", "kendall", "spearman"], method)
    matrix = frame.corr(method=corr_method).to_numpy(dtype=float)
    if not np.isfinite(matrix).all():
        raise ValueError("correlation matrix contains non-finite values.")
    if not np.allclose(matrix, matrix.T) or not np.allclose(np.diag(matrix), 1.0):
        raise ValueError("correlation matrix is not symmetric with a unit diagonal.")

    linkage_z = np.asarray(
        linkage(matrix, metric=linkage_metric, method=linkage_method), dtype=float
    )
    leaves = dendrogram(linkage_z, no_plot=True)["leaves"]
    order = np.asarray(leaves, dtype=int)
    if sorted(order.tolist()) != list(range(n_samples)):
        raise ValueError("dendrogram leaf order is not a permutation of the samples.")
    return CorrelationResult(
        matrix=matrix,
        sample_ids=sample_ids,
        method=method,
        order=order,
        linkage_z=linkage_z,
        n_features=int(n_features),
    )


# --------------------------------------------------------------------------- #
# Annotation resolution (read-only registry, fail loud)
# --------------------------------------------------------------------------- #


def _resolve_annotations(
    dataset: Dataset,
    annotations: Mapping[str, str],
    registry_path: str | Path,
) -> list[_Annotation]:
    """Resolve ``{metadata column: display label}`` to registry-colored stripes."""
    registry = load_registry(registry_path)
    resolved: list[_Annotation] = []
    for column, label in annotations.items():
        if column not in dataset.metadata.columns:
            raise ValueError(
                f"annotation {column!r} is not a metadata column "
                f"{list(dataset.metadata.columns)}."
            )
        series = dataset.metadata[column]
        if bool(series.isna().any()):
            raise ValueError(f"annotation {column!r} has missing values.")
        values = series.to_numpy().astype(str)
        distinct = list(dict.fromkeys(values.tolist()))
        if len(distinct) > _MAX_CATEGORIES:
            raise ValueError(
                f"annotation {column!r} has {len(distinct)} values (> "
                f"{_MAX_CATEGORIES}); change the encoding strategy, do not add colors."
            )
        entry = registry.get(column)
        if not isinstance(entry, dict) or not isinstance(entry.get("values"), dict):
            raise ValueError(f"category {column!r} is not in the color registry.")
        registered: dict[str, str] = {
            str(k): str(v) for k, v in entry["values"].items()
        }
        missing = [v for v in distinct if v not in registered]
        if missing:
            raise ValueError(
                f"value(s) {missing} of {column!r} have no registry color; register "
                f"them in {registry_path} first."
            )
        resolved.append(
            _Annotation(
                column=column,
                label=label,
                per_sample=values,
                colors={v: registered[v] for v in distinct},
            )
        )
    return resolved


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class _Axes:
    heatmap: Axes
    colorbar: Axes
    stripes: list[Axes]
    dendrogram: Axes


def _build_layout(fig: Figure, n_annotations: int) -> _Axes:
    """Stack dendrogram / stripes / heatmap in one column + a thin colorbar column."""
    height_ratios = [1.6] + [0.32] * n_annotations + [8.0]
    gs = GridSpec(
        len(height_ratios),
        2,
        figure=fig,
        width_ratios=[24, 1],
        height_ratios=height_ratios,
        hspace=0.08,
        wspace=0.04,
        left=0.24,
        right=0.90,
        top=0.88,
        bottom=0.13,
    )
    dendro_ax = fig.add_subplot(gs[0, 0])
    stripes = [fig.add_subplot(gs[1 + i, 0]) for i in range(n_annotations)]
    heatmap_ax = fig.add_subplot(gs[-1, 0])
    colorbar_ax = fig.add_subplot(gs[-1, 1])
    return _Axes(
        heatmap=heatmap_ax, colorbar=colorbar_ax, stripes=stripes, dendrogram=dendro_ax
    )


def _draw_dendrogram(ax: Axes, linkage_z: np.ndarray, n_samples: int) -> None:
    """Top dendrogram mapped onto heatmap column coordinates (scipy leaf i at 10i+5)."""
    dendro = dendrogram(linkage_z, no_plot=True)
    for xs, ys in zip(dendro["icoord"], dendro["dcoord"], strict=True):
        ax.plot([(x - 5.0) / 10.0 for x in xs], ys, color="#444444", linewidth=1.0)
    ax.set_xlim(-0.5, n_samples - 0.5)
    ax.set_ylim(bottom=0.0)
    ax.axis("off")


def _draw_stripes(
    axes: list[Axes], annotations: list[_Annotation], order: np.ndarray
) -> None:
    for ax, annotation in zip(axes, annotations, strict=True):
        ordered = annotation.per_sample[order]
        rgba = np.array([to_rgba(annotation.colors[str(v)]) for v in ordered])
        n = len(ordered)
        ax.imshow(
            rgba.reshape(1, -1, 4),
            aspect="auto",
            interpolation="none",
            extent=(-0.5, n - 0.5, -0.5, 0.5),
        )
        # Thin white separators so adjacent same-colored cells stay countable.
        for x in np.arange(0.5, n - 0.5, 1.0):
            ax.axvline(x, color="white", linewidth=0.8)
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_ylabel(
            annotation.label,
            rotation=0,
            ha="right",
            va="center",
            fontsize=11,
            labelpad=8,
        )


def _text_color(rgba: tuple[float, float, float, float]) -> str:
    """Black or white text, whichever contrasts with the cell (relative luminance)."""
    r, g, b = rgba[:3]
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "black" if luminance > 0.5 else "white"


def _draw_heatmap(
    ax: Axes,
    colorbar_ax: Axes,
    matrix: np.ndarray,
    *,
    labels: np.ndarray,
    method: str,
    annotate_cells: bool,
    decimals: int,
) -> None:
    """Heatmap with a masked diagonal; color scale spans the off-diagonal range."""
    n = len(matrix)
    diag = np.eye(n, dtype=bool)
    off = matrix[~diag]
    vmin = float(off.min())
    vmax = float(off.max())
    if vmax <= vmin:
        vmax = vmin + 1e-6
    masked = np.ma.masked_array(matrix, mask=diag)
    cmap = plt.get_cmap(_HEATMAP_CMAP).with_extremes(bad=_DIAGONAL_COLOR)
    norm = Normalize(vmin=vmin, vmax=vmax)
    image = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto", interpolation="none")

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=10)
    ax.set_yticklabels(labels, fontsize=10)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if annotate_cells:
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                value = float(matrix[i, j])
                color = _text_color(cmap(norm(value)))
                ax.text(
                    j,
                    i,
                    f"{value:.{decimals}f}",
                    ha="center",
                    va="center",
                    fontsize=8.5,
                    color=color,
                )

    cbar = ax.figure.colorbar(image, cax=colorbar_ax)
    cbar.set_label(_METHOD_SYMBOL[method], fontsize=11)
    cbar.ax.tick_params(labelsize=9)


def _legend_figure(annotations: list[_Annotation]) -> Figure:
    """Standalone swatch key: one titled block per annotation stripe."""
    # Each block's height is proportional to its rows (title + one row per value),
    # so blocks stack compactly without large gaps.
    rows = [1 + len(a.colors) for a in annotations]
    fig, raw_axes = plt.subplots(
        len(annotations),
        1,
        figsize=(2.4, 0.27 * sum(rows) + 0.25 * len(rows)),
        gridspec_kw={"height_ratios": rows, "hspace": 0.0},
    )
    for ax, annotation in zip(list(np.atleast_1d(raw_axes)), annotations, strict=True):
        ax.axis("off")
        handles = [
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="",
                markersize=11,
                markerfacecolor=color,
                markeredgecolor="none",
            )
            for color in annotation.colors.values()
        ]
        ax.legend(
            handles,
            list(annotation.colors.keys()),
            title=annotation.label,
            loc="center left",
            frameon=False,
            fontsize=10,
            title_fontsize=11,
            alignment="left",
            borderaxespad=0.0,
            handletextpad=0.4,
            labelspacing=0.3,
        )
    fig.tight_layout()
    return fig


def plot_sample_correlation(
    dataset: Dataset,
    *,
    annotations: Mapping[str, str],
    method: str = "pearson",
    title: str | None = None,
    annotate_cells: bool | None = None,
    decimals: int = 3,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    linkage_metric: str = _LINKAGE_METRIC,
    linkage_method: str = _LINKAGE_METHOD,
) -> CorrelationPlot:
    """Render a clustered sample-correlation heatmap with annotation stripes.

    Parameters
    ----------
    annotations:
        Ordered ``{metadata column: display label}`` of categorical stripes; each
        column is a registry category and every value must already be registered.
    annotate_cells:
        Write r in each off-diagonal cell. ``None`` (default) = only when
        n <= ``_MAX_ANNOTATED``.
    decimals:
        Decimal places for the cell text.
    """
    resolved = _resolve_annotations(dataset, annotations, registry_path)
    result = compute_correlation(
        dataset,
        method=method,
        linkage_metric=linkage_metric,
        linkage_method=linkage_method,
    )
    if method == "pearson" and dataset.scale not in LOG_SCALES:
        warnings.warn(
            f"Pearson correlation on non-log scale {dataset.scale!r}; abundant "
            f"features dominate. Log-transform first or use method='spearman'.",
            CorrelationScaleWarning,
            stacklevel=2,
        )
    n = len(result.sample_ids)
    do_annotate = n <= _MAX_ANNOTATED if annotate_cells is None else annotate_cells
    order = result.order

    fig = plt.figure(figsize=(7.6, 7.8))
    try:
        axes = _build_layout(fig, len(resolved))
        if title is not None:
            fig.suptitle(title, fontsize=13, weight="bold", y=0.975)
        _draw_dendrogram(axes.dendrogram, result.linkage_z, n)
        _draw_stripes(axes.stripes, resolved, order)
        _draw_heatmap(
            axes.heatmap,
            axes.colorbar,
            result.matrix[np.ix_(order, order)],
            labels=result.sample_ids[order],
            method=method,
            annotate_cells=do_annotate,
            decimals=decimals,
        )
        legend_fig = _legend_figure(resolved) if resolved else None
    except BaseException:
        plt.close(fig)
        raise
    return CorrelationPlot(
        figure=fig,
        legend_figure=legend_fig,
        result=result,
        color_maps={a.column: dict(a.colors) for a in resolved},
    )


def save_sample_correlation(
    dataset: Dataset,
    output_dir: str | Path,
    base_name: str,
    *,
    annotations: Mapping[str, str],
    method: str = "pearson",
    title: str | None = None,
    annotate_cells: bool | None = None,
    decimals: int = 3,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    dpi: int = 300,
) -> tuple[FigureArtifacts, CorrelationPlot]:
    """Render inside the publication style and dual-export figure + legend."""
    with publication_style():
        plot = plot_sample_correlation(
            dataset,
            annotations=annotations,
            method=method,
            title=title,
            annotate_cells=annotate_cells,
            decimals=decimals,
            registry_path=registry_path,
        )
        artifacts = save_figure(
            plot.figure, output_dir, base_name, legend_fig=plot.legend_figure, dpi=dpi
        )
    return artifacts, plot
