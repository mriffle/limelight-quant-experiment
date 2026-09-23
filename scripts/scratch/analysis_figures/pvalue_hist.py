"""p-value histogram — calibration diagnostic (project copy of pvalue-hist@0.2).

Seeded from the plugin template ``lib/figures/pvalue_hist.py`` (v0.2). Draws the
distribution of **raw** per-feature p-values over ``[0, 1]`` as a density, with the
uniform null (density 1) as a reference line, and reports Storey's π0 per distribution.
Values outside ``[0, 1]`` raise (the "passed q by mistake" slip). Colors come from the
Okabe-Ito palette via :mod:`common.figures.colors`, **not persisted** by default
(``persist_colors=False``) because the overlay labels are figure-local.

Project deviations from the template (each deliberate):

* **Imports** point at the project packages.
* **Second channel for overlays** — ``linestyles`` gives each overlaid step histogram a
  line style as well as a color, so near-coincident distributions (and grayscale
  prints) stay distinguishable; the legend draws line swatches to match.
* **Uncapped π0 reported** — the template clamps π0 to 1; when the raw Storey
  estimate exceeds 1 (a deficit of small p / hump near 1) the legend says so
  (``π0 = 1.00, capped; raw 1.07``) instead of silently showing 1.00.
* **Legend keys the uniform reference** line as well as the distributions.
* :func:`plot_pvalue_small_multiples` — one bar histogram per panel (e.g. one per
  quantity), shared axes, with the uniform line, a dashed π0 line and a terse on-panel
  ``π0 ≈ x`` label, for a compact cross-quantity overview.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from common.figures.colors import DEFAULT_REGISTRY_PATH, assign_colors
from common.figures.figure_io import publication_style
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.typing import LineStyleType

__script_meta__: dict[str, object] = {
    "template": {"name": "pvalue-hist", "version": "0.2"},
    "kind": "module",
    "provides": [
        "PValueHistogramResult",
        "PValueHistogramPlot",
        "estimate_pi0",
        "estimate_pi0_raw",
        "plot_pvalue_histogram",
        "plot_pvalue_small_multiples",
    ],
    "uses": ["common.figures.colors", "common.figures.figure_io"],
    "seeded_from": "pvalue-hist@0.2",
    "description": (
        "Raw p-value calibration histogram (density over [0, 1], uniform reference, "
        "Storey pi0 per distribution; rejects values outside [0, 1]). Project "
        "additions: per-overlay line styles (second channel), uncapped-pi0 note when "
        "the Storey estimate exceeds 1, uniform line keyed in the legend, and a "
        "small-multiples variant with per-panel pi0 lines/labels. Figure-local colors "
        "(persist_colors=False). Separate legend image. Fail-loud."
    ),
}

DEFAULT_PVALUE_CATEGORY = "PValueDistribution"
_DEFAULT_PI0_LAMBDA = 0.5
_UNIFORM_COLOR = "0.25"


@dataclass(frozen=True)
class PValueHistogramResult:
    """Per-distribution summary underlying the figure.

    Attributes
    ----------
    counts:
        ``{label: number of finite p-values}``.
    pi0:
        ``{label: Storey π0 clamped to (0, 1]}``.
    pi0_raw:
        ``{label: unclamped Storey estimate}`` (may exceed 1).
    """

    counts: dict[str, int]
    pi0: dict[str, float]
    pi0_raw: dict[str, float]


@dataclass
class PValueHistogramPlot:
    """A rendered p-value histogram plus its companion legend figure."""

    figure: Figure
    legend_figure: Figure
    result: PValueHistogramResult
    color_map: dict[str, str]


def estimate_pi0_raw(pvalues: np.ndarray, lam: float = _DEFAULT_PI0_LAMBDA) -> float:
    """Unclamped Storey ``#{p > lam} / (m (1 - lam))`` over finite p; NaN if none."""
    if not 0.0 < lam < 1.0:
        raise ValueError(f"lam must be in (0, 1); got {lam}.")
    p = np.asarray(pvalues, dtype=float)
    finite = p[np.isfinite(p)]
    if finite.size == 0:
        return float("nan")
    return float(np.sum(finite > lam)) / (finite.size * (1.0 - lam))


def estimate_pi0(pvalues: np.ndarray, lam: float = _DEFAULT_PI0_LAMBDA) -> float:
    """Storey's fixed-``lam`` π0, clamped to ``(0, 1]`` (template behaviour)."""
    raw = estimate_pi0_raw(pvalues, lam)
    return raw if math.isnan(raw) else float(min(raw, 1.0))


def _validate_pvalues(label: str, p: np.ndarray) -> np.ndarray:
    """Return the finite p-values, raising if any finite value is outside ``[0, 1]``."""
    if p.ndim != 1:
        raise ValueError(f"{label!r}: p-values must be 1D; got {p.ndim}D.")
    finite = p[np.isfinite(p)]
    if finite.size == 0:
        raise ValueError(f"{label!r}: no finite p-values.")
    if finite.min() < 0.0 or finite.max() > 1.0:
        raise ValueError(
            f"{label!r}: p-values must lie in [0, 1]; got "
            f"[{finite.min():.4g}, {finite.max():.4g}]. Pass raw p-values, not BH "
            f"q-values or test statistics."
        )
    return np.asarray(finite, dtype=float)


def _summarize(
    data: Mapping[str, np.ndarray], lam: float
) -> tuple[dict[str, np.ndarray], PValueHistogramResult]:
    if len(data) == 0:
        raise ValueError("pvalues mapping is empty; pass >= 1 distribution.")
    finite = {
        str(k): _validate_pvalues(str(k), np.asarray(v, dtype=float))
        for k, v in data.items()
    }
    result = PValueHistogramResult(
        counts={k: int(v.size) for k, v in finite.items()},
        pi0={k: estimate_pi0(v, lam) for k, v in finite.items()},
        pi0_raw={k: estimate_pi0_raw(v, lam) for k, v in finite.items()},
    )
    return finite, result


def pi0_text(result: PValueHistogramResult, label: str) -> str:
    """``π0 = 0.86`` or ``π0 = 1.00, capped; raw 1.07`` for one distribution."""
    pi0 = result.pi0[label]
    raw = result.pi0_raw[label]
    if raw > 1.0:
        return f"π0 = {pi0:.2f}, capped; raw {raw:.2f}"
    return f"π0 = {pi0:.2f}"


def _bins(n_bins: int) -> list[float]:
    if n_bins < 1:
        raise ValueError(f"n_bins must be >= 1; got {n_bins}.")
    return [float(edge) for edge in np.linspace(0.0, 1.0, n_bins + 1)]


# --------------------------------------------------------------------------- #
# Overlay (one axes, several distributions)
# --------------------------------------------------------------------------- #
def plot_pvalue_histogram(
    pvalues: Mapping[str, np.ndarray],
    *,
    n_bins: int = 20,
    pi0_lambda: float = _DEFAULT_PI0_LAMBDA,
    linestyles: Mapping[str, LineStyleType] | None = None,
    category: str = DEFAULT_PVALUE_CATEGORY,
    title: str | None = None,
    xlabel: str = "raw p-value",
    legend_title: str = "Distribution",
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = False,
) -> PValueHistogramPlot:
    """Overlay raw-p distributions (step histograms) with the uniform reference.

    ``pvalues`` is an ordered ``{label: array}``; ``NaN`` is dropped, a finite value
    outside ``[0, 1]`` raises. ``linestyles`` (optional) maps each label to a
    matplotlib line style (a second channel beside color). π0 is reported per label in
    the separate legend.
    """
    bins = _bins(n_bins)
    finite, result = _summarize(pvalues, pi0_lambda)
    labels = list(finite)
    styles: dict[str, LineStyleType] = dict.fromkeys(labels, "-")
    if linestyles is not None:
        unknown = set(linestyles) - set(labels)
        if unknown:
            raise ValueError(f"linestyles has labels not in pvalues: {sorted(unknown)}")
        styles.update(linestyles)

    with publication_style():
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        try:
            color_map = assign_colors(
                category, labels, registry_path=registry_path, persist=persist_colors
            )
            for label in labels:
                # stairs with no baseline: no vertical drop lines at p = 0 / p = 1.
                density, edges = np.histogram(finite[label], bins=bins, density=True)
                ax.stairs(
                    density,
                    edges,
                    baseline=None,
                    color=color_map[label],
                    linewidth=2.0,
                    linestyle=styles[label],
                )
            ax.axhline(1.0, color=_UNIFORM_COLOR, linestyle=(0, (1, 2)), linewidth=1.2)
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(bottom=0.0)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("density (uniform = 1)")
            if title is not None:
                fig.suptitle(title, fontsize=12, weight="bold")
            handles: list[Artist] = [
                Line2D(
                    [0],
                    [0],
                    color=color_map[label],
                    linestyle=styles[label],
                    linewidth=2.0,
                )
                for label in labels
            ]
            texts = [
                f"{label} (n={result.counts[label]:,}; {pi0_text(result, label)})"
                for label in labels
            ]
            handles.append(
                Line2D([0], [0], color=_UNIFORM_COLOR, linestyle=(0, (1, 2)))
            )
            texts.append("uniform null (density 1)")
            legend_figure = _legend_figure(handles, texts, legend_title)
        except BaseException:
            plt.close(fig)
            raise

    return PValueHistogramPlot(
        figure=fig, legend_figure=legend_figure, result=result, color_map=color_map
    )


# --------------------------------------------------------------------------- #
# Small multiples (one distribution per panel)
# --------------------------------------------------------------------------- #
def plot_pvalue_small_multiples(
    pvalues: Mapping[str, np.ndarray],
    *,
    series_label: str,
    ncols: int = 2,
    n_bins: int = 20,
    pi0_lambda: float = _DEFAULT_PI0_LAMBDA,
    category: str = DEFAULT_PVALUE_CATEGORY,
    color_order: tuple[str, ...] | None = None,
    title: str | None = None,
    xlabel: str = "raw p-value",
    legend_title: str = "Key",
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = False,
) -> PValueHistogramPlot:
    """One bar histogram per panel (``{panel title: p array}``), shared axes.

    Every panel is the same series (``series_label``, one color). ``color_order``, when
    given, is the label list to assign colors over (so ``series_label`` gets the same
    palette slot it has in a companion overlay figure); it must contain
    ``series_label``. Each panel carries the uniform line, a dotted line at its Storey
    π0 and a terse ``π0 ≈ x`` label.
    """
    if ncols < 1:
        raise ValueError(f"ncols must be >= 1; got {ncols}.")
    bins = _bins(n_bins)
    finite, result = _summarize(pvalues, pi0_lambda)
    panels = list(finite)
    order = list(color_order) if color_order is not None else [series_label]
    if series_label not in order:
        raise ValueError(f"color_order must contain series_label {series_label!r}.")
    nrows = math.ceil(len(panels) / ncols)

    with publication_style():
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(4.0 * ncols, 3.1 * nrows + 0.5),
            sharex=True,
            sharey=True,
            squeeze=False,
        )
        try:
            color_map = assign_colors(
                category, order, registry_path=registry_path, persist=persist_colors
            )
            color = color_map[series_label]
            flat: list[Axes] = [a for row in axes for a in row]
            for ax, panel in zip(flat, panels, strict=False):
                _draw_panel(ax, finite[panel], bins, color, result, panel)
                ax.set_title(panel, fontsize=11)
            for ax in flat[len(panels) :]:
                ax.set_visible(False)
            for row in axes:
                row[0].set_ylabel("density (uniform = 1)")
            for ax in axes[-1]:
                ax.set_xlabel(xlabel)
            flat[0].set_xlim(0.0, 1.0)
            flat[0].set_ylim(bottom=0.0)
            if title is not None:
                fig.suptitle(title, fontsize=12, weight="bold")
            fig.tight_layout()
            handles: list[Artist] = [
                Patch(facecolor=color, alpha=0.6, edgecolor="none"),
                Line2D([0], [0], color=_UNIFORM_COLOR, linestyle=(0, (1, 2))),
                Line2D([0], [0], color="black", linestyle=(0, (5, 3)), linewidth=1.3),
            ]
            texts = [
                f"{series_label}: raw p",
                "uniform null (density 1)",
                f"Storey π0 (λ = {pi0_lambda:g})",
            ]
            legend_figure = _legend_figure(handles, texts, legend_title)
        except BaseException:
            plt.close(fig)
            raise

    return PValueHistogramPlot(
        figure=fig,
        legend_figure=legend_figure,
        result=result,
        color_map={series_label: color},
    )


def _draw_panel(
    ax: Axes,
    values: np.ndarray,
    bins: list[float],
    color: str,
    result: PValueHistogramResult,
    label: str,
) -> None:
    ax.hist(
        values,
        bins=bins,
        density=True,
        color=color,
        alpha=0.6,
        edgecolor="white",
        linewidth=0.5,
    )
    ax.axhline(1.0, color=_UNIFORM_COLOR, linestyle=(0, (1, 2)), linewidth=1.2)
    pi0 = result.pi0[label]
    ax.axhline(pi0, color="black", linestyle=(0, (5, 3)), linewidth=1.3)
    ax.text(
        0.98,
        0.96,
        f"π0 ≈ {pi0:.2f} · n = {result.counts[label]:,}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.5},
    )


def _legend_figure(
    handles: list[Artist], texts: list[str], legend_title: str
) -> Figure:
    """Standalone legend figure."""
    height = 0.3 * len(texts) + 0.45
    width = max(3.6, 0.105 * max(len(t) for t in texts) + 1.2)
    fig, ax = plt.subplots(figsize=(width, height))
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
