"""PCA processing-state series (small multiples) for Stage-3 QC.

Project module seeded from the plugin ``lib/figures`` template ``pca-plot`` v0.3.
Held to the correctness charter (conventions/correctness.md): **assume nothing,
verify everything, fail loud.**

What it draws (one figure): one **PC1 vs PC2** scatter panel per processing state
(e.g. raw -> median-normalized -> batch-corrected), side by side, points colored by
one categorical sample-metadata column from the project color registry and labeled
directly with a sample identifier. PCA is computed **independently per state** on
per-feature standardized abundances (correlation-matrix PCA), so each panel's axis
labels carry that state's own % variance. PC signs are arbitrary per panel.

DEVIATIONS FROM THE TEMPLATE (pca-plot v0.3), and why:

  * **Small multiples of states instead of PC1/PC2 + PC3/PC4.** The QC convention
    (conventions/visualization.md, ``pca-plot``) renders the family as a
    processing-state series, one panel per state; the Stage-3 dispatch asked for
    PC1/PC2 only.
  * **No marginal KDEs / per-PC group tests.** This study has n = 8 samples with a
    6-vs-2 batch split: a 2-point KDE is meaningless and the smallest attainable
    two-sided Mann-Whitney p for 6 vs 2 is 2/28 = 0.071, so the annotation would be
    noise dressed as a test. Direct sample labels carry more information here.
  * **Direct point labels** (``label_by``), placed by a small greedy
    collision-avoiding search (``textalloc`` is not installed in this project).
  * **Registry colors are read-only.** Every value colored must already be in the
    registry (``state/color_registry.json``); an unknown value raises rather than
    being assigned a new color, so concurrent QC families never race on the registry.
  * Imports resolve against this project's layout (``loaders.*``,
    ``common.figures.*``).

Scale: PCA is meaningful on a log-ish scale; a non-log ``Dataset`` triggers
:class:`PCAScaleWarning` (as in the template), it does not refuse.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from common.figures.colors import DEFAULT_REGISTRY_PATH, load_registry
from common.figures.figure_io import FigureArtifacts, publication_style, save_figure
from loaders.data_loading import LOG_SCALES, Dataset
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

__script_meta__: dict[str, object] = {
    "template": {"name": "pca-plot", "version": "0.3"},
    "kind": "module",
    "provides": [
        "PCAScaleWarning",
        "PCAResult",
        "StatePanel",
        "PCASeriesPlot",
        "compute_pca",
        "registry_colors",
        "plot_pca_state_series",
        "save_pca_state_series",
    ],
    "uses": [
        "loaders.data_loading",
        "common.figures.colors",
        "common.figures.figure_io",
    ],
    "seeded_from": "lib/figures/pca.py (pca-plot v0.3)",
    "description": (
        "PCA processing-state series: one PC1/PC2 panel per state (raw -> normalized "
        "-> batch-corrected), per-state standardized PCA, points colored by a "
        "registry category (read-only, fail on unknown value, >8 guard) and directly "
        "labeled by sample id; dual export + separate swatch legend."
    ),
}

_N_COMPONENTS = 2
_MAX_CATEGORICAL = 8
_MARKER_SIZE = 90.0  # scatter `s` (points^2)
_LABEL_FONTSIZE = 10.0
_LABEL_OFFSET_PT = 7.0


class PCAScaleWarning(UserWarning):
    """Warning that PCA is running on a non-log-ish (e.g. linear) abundance scale."""


@dataclass(frozen=True)
class PCAResult:
    """Outcome of :func:`compute_pca`.

    Attributes
    ----------
    scores:
        ``(n_samples, n_components)`` PCA scores.
    explained_variance_ratio:
        ``(n_components,)`` fraction of total variance per component, in ``[0, 1]``.
    standardized:
        Whether features were z-scored before PCA.
    n_features:
        Number of features the PCA ran on.
    """

    scores: np.ndarray
    explained_variance_ratio: np.ndarray
    standardized: bool
    n_features: int


@dataclass(frozen=True)
class StatePanel:
    """One panel of the series: a processing state's label and its data."""

    label: str
    dataset: Dataset


@dataclass
class PCASeriesPlot:
    """A rendered state-series figure, its legend figure, and per-panel results."""

    figure: Figure
    legend_figure: Figure
    results: dict[str, PCAResult]
    color_map: dict[str, str]


# --------------------------------------------------------------------------- #
# Compute
# --------------------------------------------------------------------------- #


def compute_pca(
    dataset: Dataset,
    *,
    n_components: int = _N_COMPONENTS,
    standardize: bool = True,
    random_state: int = 0,
) -> PCAResult:
    """Run PCA on a :class:`Dataset`'s abundances (samples x features).

    Features are z-scored first when ``standardize`` (correlation-matrix PCA). The
    full SVD is deterministic; ``random_state`` is plumbed through only so it can be
    recorded. Raises on non-finite abundances or an impossible ``n_components``.
    """
    abundances = np.asarray(dataset.abundances, dtype=float)
    if abundances.ndim != 2:
        raise ValueError(
            f"abundances must be 2D (n_samples, n_features); got {abundances.shape}."
        )
    n_samples, n_features = abundances.shape
    if not np.isfinite(abundances).all():
        n_bad = int((~np.isfinite(abundances)).sum())
        raise ValueError(
            f"compute_pca requires finite abundances but found {n_bad} non-finite "
            f"value(s). Resolve missingness before PCA."
        )
    if n_samples < 2:
        raise ValueError(f"PCA needs >= 2 samples; got {n_samples}.")
    max_components = min(n_samples, n_features)
    if not 1 <= n_components <= max_components:
        raise ValueError(
            f"n_components={n_components} must be between 1 and "
            f"min(n_samples, n_features)={max_components}."
        )
    if standardize:
        constant = int((abundances.std(axis=0) == 0).sum())
        if constant:
            raise ValueError(
                f"{constant} feature(s) are constant across samples; standardization "
                f"would divide by zero. Drop them before PCA."
            )
    matrix = StandardScaler().fit_transform(abundances) if standardize else abundances
    pca = PCA(n_components=n_components, svd_solver="full", random_state=random_state)
    scores = np.asarray(pca.fit_transform(matrix), dtype=float)
    ratio = np.asarray(pca.explained_variance_ratio_, dtype=float)
    return PCAResult(
        scores=scores,
        explained_variance_ratio=ratio,
        standardized=standardize,
        n_features=n_features,
    )


# --------------------------------------------------------------------------- #
# Colors (read-only registry lookup)
# --------------------------------------------------------------------------- #


def registry_colors(
    category: str,
    values: Sequence[str],
    *,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
) -> dict[str, str]:
    """Return ``{value: hex}`` for ``values`` from the registry, read-only.

    Raises if the category or any value is absent (no silent new color), or if more
    than eight distinct values would be colored (the >8-category rule).
    """
    unique = list(dict.fromkeys(str(v) for v in values))
    if len(unique) > _MAX_CATEGORICAL:
        raise ValueError(
            f"{len(unique)} categories for {category!r} exceeds the "
            f"{_MAX_CATEGORICAL}-color budget; facet or use a second channel instead."
        )
    registry = load_registry(registry_path)
    entry = registry.get(category)
    if not isinstance(entry, dict) or not isinstance(entry.get("values"), dict):
        raise ValueError(
            f"Color registry {registry_path} has no category {category!r} with a "
            f"'values' map."
        )
    known: dict[str, object] = entry["values"]
    missing = [v for v in unique if v not in known]
    if missing:
        raise ValueError(
            f"Values {missing} of category {category!r} are not in the color registry "
            f"{registry_path}; add them there first (colors are read-only here)."
        )
    out: dict[str, str] = {}
    for value in unique:
        color = known[value]
        if not isinstance(color, str):
            raise ValueError(f"Registry color for {category}/{value} is not a string.")
        out[value] = color
    return out


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #


def _validate_panels(panels: Sequence[StatePanel], columns: Sequence[str]) -> None:
    """Fail loud unless every panel has the same samples, in the same order."""
    if not panels:
        raise ValueError("At least one state panel is required.")
    ref = panels[0].dataset.metadata
    for col in columns:
        if col not in ref.columns:
            raise ValueError(
                f"Column {col!r} is not in the sample metadata {list(ref.columns)}."
            )
        if ref[col].isna().any():
            raise ValueError(f"Metadata column {col!r} has missing values.")
    ref_ids = ref.index.astype(str).tolist()
    for panel in panels:
        meta = panel.dataset.metadata
        if meta.index.astype(str).tolist() != ref_ids:
            raise ValueError(
                f"Panel {panel.label!r} samples/order differ from panel "
                f"{panels[0].label!r}; all states must share one sample order."
            )
        for col in columns:
            if not meta[col].astype(str).equals(ref[col].astype(str)):
                raise ValueError(
                    f"Panel {panel.label!r} metadata column {col!r} differs from "
                    f"panel {panels[0].label!r}."
                )
        if panel.dataset.abundances.shape[0] != len(ref_ids):
            raise ValueError(
                f"Panel {panel.label!r} has {panel.dataset.abundances.shape[0]} rows "
                f"but {len(ref_ids)} metadata rows."
            )
        if panel.dataset.scale not in LOG_SCALES:
            warnings.warn(
                f"Panel {panel.label!r}: PCA on scale {panel.dataset.scale!r}, which "
                f"is not log-ish {sorted(LOG_SCALES)}.",
                PCAScaleWarning,
                stacklevel=3,
            )


def plot_pca_state_series(
    panels: Sequence[StatePanel],
    color_by: str,
    *,
    label_by: str | None = None,
    category: str | None = None,
    title: str | None = None,
    legend_title: str | None = None,
    standardize: bool = True,
    random_state: int = 0,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    panel_size: float = 4.6,
) -> PCASeriesPlot:
    """Render one PC1/PC2 panel per processing state, colored by ``color_by``.

    Parameters
    ----------
    panels:
        Ordered processing states (left -> right). All must share samples and order.
    color_by:
        Categorical sample-metadata column coloring the points.
    label_by:
        Optional metadata column whose value labels each point directly.
    category:
        Registry namespace for ``color_by`` values; defaults to ``color_by``.
    title:
        Optional short figure suptitle.
    legend_title:
        Title of the separate legend image; defaults to ``color_by``.
    standardize, random_state:
        Forwarded to :func:`compute_pca`.
    registry_path:
        The color registry (read-only).
    panel_size:
        Panel width/height in inches.
    """
    columns = [color_by] + ([label_by] if label_by is not None else [])
    _validate_panels(panels, columns)
    meta = panels[0].dataset.metadata
    groups = meta[color_by].astype(str).to_numpy()
    order = list(dict.fromkeys(groups.tolist()))
    color_map = registry_colors(
        category if category is not None else color_by,
        order,
        registry_path=registry_path,
    )
    point_labels = meta[label_by].astype(str).tolist() if label_by is not None else None

    results: dict[str, PCAResult] = {}
    with publication_style():
        n = len(panels)
        fig, axes_arr = plt.subplots(
            1, n, figsize=(panel_size * n + 0.6, panel_size + 0.9), squeeze=False
        )
        try:
            axes: list[Axes] = list(axes_arr[0])
            for ax, panel in zip(axes, panels, strict=True):
                result = compute_pca(
                    panel.dataset, standardize=standardize, random_state=random_state
                )
                results[panel.label] = result
                _draw_panel(ax, result, groups, order, color_map)
                pct = result.explained_variance_ratio * 100.0
                ax.set_xlabel(f"PC1 ({pct[0]:.1f}% variance)")
                ax.set_ylabel(f"PC2 ({pct[1]:.1f}% variance)")
                ax.set_title(panel.label, fontsize=12)
            if title is not None:
                fig.suptitle(title, fontsize=14, weight="bold")
            fig.tight_layout()
            if point_labels is not None:
                for ax, panel in zip(axes, panels, strict=True):
                    _place_labels(ax, fig, results[panel.label].scores, point_labels)
            legend_figure = _legend_figure(
                order, color_map, legend_title if legend_title else color_by
            )
        except BaseException:
            plt.close(fig)
            raise
    return PCASeriesPlot(
        figure=fig, legend_figure=legend_figure, results=results, color_map=color_map
    )


def save_pca_state_series(
    panels: Sequence[StatePanel],
    color_by: str,
    output_dir: str | Path,
    base_name: str,
    *,
    label_by: str | None = None,
    category: str | None = None,
    title: str | None = None,
    legend_title: str | None = None,
    standardize: bool = True,
    random_state: int = 0,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    dpi: int = 300,
) -> tuple[FigureArtifacts, PCASeriesPlot]:
    """Render (:func:`plot_pca_state_series`) and dual-export with a legend image."""
    with publication_style():
        plot = plot_pca_state_series(
            panels,
            color_by,
            label_by=label_by,
            category=category,
            title=title,
            legend_title=legend_title,
            standardize=standardize,
            random_state=random_state,
            registry_path=registry_path,
        )
        artifacts = save_figure(
            plot.figure, output_dir, base_name, legend_fig=plot.legend_figure, dpi=dpi
        )
    return artifacts, plot


def _draw_panel(
    ax: Axes,
    result: PCAResult,
    groups: np.ndarray,
    order: Sequence[str],
    color_map: dict[str, str],
) -> None:
    """Scatter PC1/PC2 per group and pad limits so direct labels fit inside."""
    scores = result.scores
    for group in order:
        mask = groups == group
        ax.scatter(
            scores[mask, 0],
            scores[mask, 1],
            color=color_map[group],
            edgecolor="k",
            linewidth=0.7,
            alpha=0.9,
            s=_MARKER_SIZE,
            zorder=3,
        )
    for dim, setter in ((0, ax.set_xlim), (1, ax.set_ylim)):
        lo = float(scores[:, dim].min())
        hi = float(scores[:, dim].max())
        span = hi - lo if hi > lo else 1.0
        pad = 0.22 * span
        setter(lo - pad, hi + pad)


# Candidate label anchors: (dx, dy) direction in units of the offset, ha, va.
_CANDIDATES: tuple[tuple[float, float, str, str], ...] = (
    (1.0, 0.6, "left", "bottom"),
    (1.0, -0.6, "left", "top"),
    (-1.0, 0.6, "right", "bottom"),
    (-1.0, -0.6, "right", "top"),
    (1.2, 0.0, "left", "center"),
    (-1.2, 0.0, "right", "center"),
    (0.0, 1.2, "center", "bottom"),
    (0.0, -1.2, "center", "top"),
)


def _overlap(a: Bbox, b: Bbox) -> float:
    """Overlap area (display px^2) of two bounding boxes."""
    w = min(a.x1, b.x1) - max(a.x0, b.x0)
    h = min(a.y1, b.y1) - max(a.y0, b.y0)
    return w * h if w > 0 and h > 0 else 0.0


def _place_labels(
    ax: Axes, fig: Figure, scores: np.ndarray, labels: Sequence[str]
) -> None:
    """Greedy collision-avoiding placement of one text label per point.

    For each point, tries eight anchor positions around the marker and keeps the one
    with the least overlap against already-placed labels, all markers, and the area
    outside the axes. Deterministic (fixed point order and candidate order).
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()  # type: ignore[attr-defined]
    axes_box = ax.get_window_extent(renderer)
    pts_to_px = fig.dpi / 72.0
    radius_px = np.sqrt(_MARKER_SIZE) / 2.0 * pts_to_px
    xy_disp = ax.transData.transform(scores[:, :2])
    marker_boxes = [
        Bbox.from_extents(x - radius_px, y - radius_px, x + radius_px, y + radius_px)
        for x, y in xy_disp
    ]
    placed: list[Bbox] = []
    for i, text in enumerate(labels):
        best_cost = np.inf
        best_artist = None
        best_box: Bbox | None = None
        for dx, dy, ha, va in _CANDIDATES:
            artist = ax.annotate(
                text,
                xy=(float(scores[i, 0]), float(scores[i, 1])),
                xytext=(dx * _LABEL_OFFSET_PT, dy * _LABEL_OFFSET_PT),
                textcoords="offset points",
                ha=ha,
                va=va,
                fontsize=_LABEL_FONTSIZE,
                zorder=4,
            )
            box = artist.get_window_extent(renderer)
            cost = sum(_overlap(box, other) for other in placed)
            cost += sum(_overlap(box, mb) for mb in marker_boxes)
            area = box.width * box.height
            cost += 10.0 * (area - _overlap(box, axes_box))  # outside the axes
            if cost < best_cost:
                if best_artist is not None:
                    best_artist.remove()
                best_cost = cost
                best_artist = artist
                best_box = box
            else:
                artist.remove()
            if cost == 0.0:
                break
        if best_box is None:  # pragma: no cover - candidates is non-empty
            raise RuntimeError("label placement produced no candidate")
        placed.append(best_box)


def _legend_figure(
    labels: Sequence[str], color_map: dict[str, str], legend_title: str
) -> Figure:
    """Standalone swatch legend, one marker per label (matches the scatter styling)."""
    height = max(1.4, 0.35 * len(labels) + 0.9)
    fig, ax = plt.subplots(figsize=(3.2, height))
    ax.axis("off")
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=10,
            markerfacecolor=color_map[label],
            markeredgecolor="k",
        )
        for label in labels
    ]
    ax.legend(
        handles,
        list(labels),
        title=legend_title,
        loc="center",
        frameon=False,
        fontsize=12,
        title_fontsize=13,
    )
    return fig
