"""Supplementary DE figures for findings 0004 / 0005 (raloxifene-vs-control).

From-scratch Stage-4 family (no ``lib/figures`` template exists for these views).
Reads the differential-abundance tables written by
``scripts/scratch/de_raloxifene_vs_control.py`` and renders three figures for claims
that no existing figure shows:

  * :func:`plot_pvalue_by_tercile` -- LFQ protein raw-p histograms as small multiples
    by mean-log2-abundance tercile (equal-count), Storey pi0 per panel;
  * :func:`plot_residual_sd_by_design` -- per-protein residual SD distributions for
    the paired / batch / unadjusted designs (step densities on a log10 axis, medians
    marked);
  * :func:`plot_abundance_agreement` -- LFQ mean log2 intensity vs the spectral
    quantities' mean log2 abundance on the common protein set (points, LOWESS trend,
    Spearman rho).

Colors: the design colors are the same figure-local ``PValueDistribution`` mapping the
0004 p-value figures use (``assign_colors(..., persist=False)`` over the design labels
in paired/batch/unadjusted order), so a design keeps its color across finding 0004;
the registry is never written. The abundance-agreement points are a single neutral
(non-categorical) series. All inputs are validated fail-loud.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common.figures.colors import assign_colors
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator
from matplotlib.typing import LineStyleType
from numpy.typing import NDArray
from statsmodels.nonparametric.smoothers_lowess import lowess

from analysis_figures.quant_comparison import (
    MINUS,
    POINT_GRAY,
    RHO,
    TREND_COLOR,
    ComparisonData,
    spearman,
)

__script_meta__: dict[str, object] = {
    "task": "fig-de-supplementary",
    "kind": "module",
    "provides": [
        "DesignSpec",
        "DESIGNS",
        "PVALUE_CATEGORY",
        "read_table",
        "storey_pi0",
        "equal_count_terciles",
        "design_colors",
        "plot_pvalue_by_tercile",
        "plot_residual_sd_by_design",
        "plot_abundance_agreement",
    ],
    "uses": ["common.figures.colors", "analysis_figures.quant_comparison"],
    "seeded_from": None,
    "description": (
        "Supplementary DE figures (findings 0004/0005): p-value histograms by "
        "abundance tercile, residual SD by design, LFQ vs spectral mean-abundance "
        "agreement. Fail-loud input checks; registry read-only."
    ),
}

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.intp]

PVALUE_CATEGORY = "PValueDistribution"  # the 0004 p-value figures' color namespace
UNIFORM_COLOR = "#808080"
PI0_COLOR = "#000000"
TERCILE_NAMES = ("low", "mid", "high")


@dataclass(frozen=True)
class DesignSpec:
    """One fitted design: table suffix, display label (== 0004 label), line style."""

    key: str
    label: str
    linestyle: LineStyleType
    covariates: str
    residual_df: int


# Labels and order are identical to scripts/scratch/fig_de_raloxifene_vs_control.py
# DESIGNS, so assign_colors() over them reproduces the 0004 design colors.
DESIGNS: tuple[DesignSpec, ...] = (
    DesignSpec("paired", "paired: condition + pair", "-", "condition + pair", 3),
    DesignSpec("batch", "batch: condition + batch", "--", "condition + batch", 5),
    DesignSpec("unadjusted", "unadjusted: condition only", "-.", "condition", 6),
)


# --------------------------------------------------------------------------- #
# Loading + statistics
# --------------------------------------------------------------------------- #
def read_table(path: Path, columns: Sequence[str], id_column: str) -> pd.DataFrame:
    """Read a DE TSV; require ``columns`` finite and ``id_column`` unique."""
    if not path.is_file():
        raise FileNotFoundError(f"Required input missing: {path}")
    frame = pd.read_csv(path, sep="\t")
    missing = [c for c in (id_column, *columns) if c not in frame.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}.")
    if not frame[id_column].is_unique:
        raise ValueError(f"{path}: duplicate ids in {id_column!r}.")
    for col in columns:
        arr = pd.to_numeric(frame[col], errors="raise").to_numpy(dtype=np.float64)
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{path}: {col} has {int((~np.isfinite(arr)).sum())} NaN.")
    return frame


def storey_pi0(pvalues: FloatArray, lam: float = 0.5) -> float:
    """Storey fixed-lambda pi0 = #{p > lam} / (m (1 - lam)), clamped to <= 1.

    Same estimator as ``de_raloxifene_vs_control.storey_pi0`` (the summary.json value).
    """
    if not 0.0 < lam < 1.0:
        raise ValueError(f"lambda must be in (0, 1); got {lam}.")
    if pvalues.size == 0 or not np.all(np.isfinite(pvalues)):
        raise ValueError("storey_pi0 needs a non-empty array of finite p-values.")
    if np.any((pvalues < 0.0) | (pvalues > 1.0)):
        raise ValueError("p-values outside [0, 1].")
    return float(min(np.sum(pvalues > lam) / (pvalues.size * (1.0 - lam)), 1.0))


def equal_count_terciles(values: FloatArray) -> list[IntArray]:
    """Row indices of the low / mid / high equal-count terciles of ``values``.

    Stable sort then ``np.array_split`` (sizes differ by at most one; the larger
    group(s) come first). Deterministic for tied values.
    """
    if values.ndim != 1 or values.size < 3:
        raise ValueError("equal_count_terciles needs a 1-D array of >= 3 values.")
    order = np.argsort(values, kind="stable")
    return [np.asarray(c, dtype=np.intp) for c in np.array_split(order, 3)]


def design_colors(registry_path: Path) -> dict[str, str]:
    """0004 design colors: palette slots over the design labels (registry unchanged)."""
    colors = assign_colors(
        PVALUE_CATEGORY,
        [d.label for d in DESIGNS],
        registry_path=registry_path,
        persist=False,
    )
    if len(set(colors.values())) != len(DESIGNS):
        raise ValueError(f"design colors collide: {colors}.")
    return colors


# --------------------------------------------------------------------------- #
# Shared drawing helpers
# --------------------------------------------------------------------------- #
def _legend_figure(
    handles: Sequence[Artist], labels: Sequence[str], title: str | None = None
) -> Figure:
    width = max(3.2, 0.085 * max(len(t) for t in labels) + 1.0)
    fig = plt.figure(figsize=(width, 0.32 * len(handles) + 0.45))
    fig.legend(
        handles,
        labels,
        loc="center",
        frameon=False,
        handlelength=2.4,
        title=title,
        title_fontsize=11,
    )
    return fig


def _fmt(value: float, digits: int) -> str:
    return f"{value:.{digits}f}".replace("-", MINUS)


# --------------------------------------------------------------------------- #
# Figure 1 -- raw-p histograms by abundance tercile
# --------------------------------------------------------------------------- #
def plot_pvalue_by_tercile(
    pvalues: FloatArray,
    abundance: FloatArray,
    *,
    color: str,
    title: str,
    subtitle: str,
    abundance_name: str,
    n_bins: int = 20,
    lam: float = 0.5,
) -> tuple[Figure, Figure, dict[str, dict[str, float | int]]]:
    """Three panels (low/mid/high abundance tercile): density histogram + pi0."""
    if pvalues.shape != abundance.shape:
        raise ValueError("pvalues and abundance must have equal shapes.")
    if n_bins < 2:
        raise ValueError(f"n_bins must be >= 2; got {n_bins}.")
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    groups = equal_count_terciles(abundance)
    densities = [np.histogram(pvalues[g], bins=edges, density=True)[0] for g in groups]
    ymax = math.ceil(max(float(np.max(d)) for d in densities) * 1.12 / 0.25) * 0.25

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.3), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.14, top=0.72, wspace=0.08)
    stats: dict[str, dict[str, float | int]] = {}
    for ax, name, idx, dens in zip(axes, TERCILE_NAMES, groups, densities, strict=True):
        p = pvalues[idx]
        a = abundance[idx]
        pi0 = storey_pi0(p, lam)
        stats[name] = {
            "n": int(idx.size),
            "abundance_min": float(np.min(a)),
            "abundance_max": float(np.max(a)),
            "pi0": pi0,
            "n_p_gt_lambda": int(np.sum(p > lam)),
            "density_first_bin": float(dens[0]),
        }
        ax.stairs(dens, edges, fill=True, color=color, alpha=0.6)
        ax.stairs(dens, edges, color="white", linewidth=0.5)
        ax.axhline(1.0, color=UNIFORM_COLOR, linestyle=(0, (1, 2)), linewidth=1.3)
        ax.axhline(pi0, color=PI0_COLOR, linestyle=(0, (5, 3)), linewidth=1.3)
        ax.text(
            0.97,
            0.96,
            f"π0 = {pi0:.3f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=11.5,
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
        )
        ax.set_title(
            f"{name} tercile · n = {idx.size:,}\n"
            f"{abundance_name} {np.min(a):.1f}\u2013{np.max(a):.1f}",
            fontsize=11,
        )
        ax.set_xlabel("raw p-value")
    axes[0].set_xlim(0.0, 1.0)
    axes[0].set_ylim(0.0, ymax)
    axes[0].set_ylabel("density (uniform = 1)")
    fig.suptitle(title, fontsize=13.5, fontweight="bold", y=0.995, va="top")
    fig.text(0.5, 0.925, subtitle, ha="center", va="top", fontsize=10.5)

    handles: list[Artist] = [
        Patch(facecolor=color, alpha=0.6, edgecolor="none"),
        Line2D([], [], color=UNIFORM_COLOR, linestyle=(0, (1, 2)), linewidth=1.3),
        Line2D([], [], color=PI0_COLOR, linestyle=(0, (5, 3)), linewidth=1.3),
    ]
    labels = [
        f"raw p, {n_bins} bins (density)",
        "uniform null (density 1)",
        f"Storey π0 (λ = {lam:g})",
    ]
    return fig, _legend_figure(handles, labels), stats


# --------------------------------------------------------------------------- #
# Figure 2 -- residual SD by design
# --------------------------------------------------------------------------- #
def plot_residual_sd_by_design(
    residual_sd: Mapping[str, FloatArray],
    colors: Mapping[str, str],
    *,
    title: str,
    subtitle: str,
    bins_per_decade: int = 20,
) -> tuple[Figure, Figure, dict[str, dict[str, float | int]]]:
    """Overlaid step densities of log10 residual SD per design, medians marked.

    ``residual_sd`` is keyed by design key (``DESIGNS`` order). The histogram is taken
    on log10(SD) (density per log10 unit) and drawn on a log-scaled SD axis.
    """
    if set(residual_sd) != {d.key for d in DESIGNS}:
        raise ValueError(f"residual_sd keys {sorted(residual_sd)} != designs.")
    for key, arr in residual_sd.items():
        if arr.size == 0 or not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
            raise ValueError(f"{key}: residual SD must be finite and > 0.")
    all_log = np.log10(np.concatenate([residual_sd[d.key] for d in DESIGNS]))
    lo = math.floor(float(np.min(all_log)) * bins_per_decade) / bins_per_decade
    hi = math.ceil(float(np.max(all_log)) * bins_per_decade) / bins_per_decade
    log_edges = np.linspace(lo, hi, round((hi - lo) * bins_per_decade) + 1)

    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.13, top=0.84)
    stats: dict[str, dict[str, float | int]] = {}
    for spec in DESIGNS:
        sd = residual_sd[spec.key]
        dens, _ = np.histogram(np.log10(sd), bins=log_edges, density=True)
        color = colors[spec.label]
        ax.stairs(
            dens,
            10.0**log_edges,
            color=color,
            linewidth=2.0,
            linestyle=spec.linestyle,
            baseline=None,
        )
        med = float(np.median(sd))
        q1, q3 = (float(v) for v in np.quantile(sd, [0.25, 0.75]))
        stats[spec.key] = {
            "n": int(sd.size),
            "median": med,
            "q1": q1,
            "q3": q3,
            "min": float(np.min(sd)),
            "max": float(np.max(sd)),
            "residual_df": spec.residual_df,
        }
        ax.axvline(med, color=color, linestyle=(0, (1, 1.5)), linewidth=1.6)
    ax.set_xscale("log")
    ax.set_xlim(10.0**lo, 10.0**hi)
    ax.set_ylim(bottom=0.0)
    ticks = [
        t
        for t in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0)
        if lo <= math.log10(t) <= hi
    ]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("per-protein residual SD (log2 units; log scale)")
    ax.set_ylabel("density (per log10 unit)")
    ax.text(
        0.98,
        0.97,
        "median",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10.5,
        fontweight="bold",
    )
    for i, spec in enumerate(DESIGNS):
        ax.text(
            0.98,
            0.97 - 0.075 * (i + 1),
            f"{spec.key} {stats[spec.key]['median']:.3f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=10.5,
            color=colors[spec.label],
        )
    fig.suptitle(title, fontsize=13.5, fontweight="bold", y=0.995, va="top")
    fig.text(0.5, 0.93, subtitle, ha="center", va="top", fontsize=10.5)

    handles: list[Artist] = [
        Line2D([], [], color=colors[d.label], linestyle=d.linestyle, linewidth=2.0)
        for d in DESIGNS
    ]
    labels = [f"{d.label} (residual df {d.residual_df})" for d in DESIGNS]
    handles.append(
        Line2D([], [], color=TREND_COLOR, linestyle=(0, (1, 1.5)), linewidth=1.6)
    )
    labels.append("median (dotted, design color)")
    return fig, _legend_figure(handles, labels, "Design"), stats


# --------------------------------------------------------------------------- #
# Figure 3 -- LFQ vs spectral mean abundance
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class AgreementPanel:
    """One head-to-head panel: y quantity key (x is always LFQ) and its axis label."""

    key: str
    label: str
    ylabel: str


def _lowess(x: FloatArray, y: FloatArray, frac: float) -> tuple[FloatArray, FloatArray]:
    fitted = np.asarray(lowess(y, x, frac=frac, it=3, return_sorted=True), dtype=float)
    if fitted.ndim != 2 or fitted.shape[1] != 2 or not np.all(np.isfinite(fitted)):
        raise ValueError("LOWESS returned a malformed or non-finite fit.")
    return fitted[:, 0], fitted[:, 1]


def plot_abundance_agreement(
    data: ComparisonData,
    panels: Sequence[AgreementPanel],
    *,
    x_key: str,
    xlabel: str,
    title: str,
    subtitle: str,
    lowess_frac: float = 0.3,
) -> tuple[Figure, Figure, dict[str, dict[str, float | int]]]:
    """LFQ mean log2 abundance (x) vs each spectral quantity's (y); rho + LOWESS."""
    x = data.column(x_key, "mean_log2_abundance")
    fig, axes = plt.subplots(1, len(panels), figsize=(5.4 * len(panels), 5.0))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.13, top=0.8, wspace=0.26)
    flat: list[Axes] = list(np.atleast_1d(axes))
    stats: dict[str, dict[str, float | int]] = {}
    for ax, panel in zip(flat, panels, strict=True):
        y = data.column(panel.key, "mean_log2_abundance")
        rho = spearman(x, y)
        pearson = float(np.corrcoef(x, y)[0, 1])
        lx, ly = _lowess(x, y, lowess_frac)
        stats[f"LFQ-{panel.label}"] = {
            "n": int(x.size),
            "spearman": rho,
            "pearson": pearson,
            "y_min": float(np.min(y)),
            "y_max": float(np.max(y)),
        }
        ax.scatter(x, y, s=7, color=POINT_GRAY, alpha=0.35, linewidths=0, zorder=2)
        ax.plot(lx, ly, color=TREND_COLOR, linewidth=2.0, zorder=4)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(panel.ylabel)
        ax.set_title(f"LFQ vs {panel.label}", fontsize=12)
        ax.text(
            0.04,
            0.96,
            f"Spearman {RHO} = {_fmt(rho, 3)}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=11.5,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
        )
    fig.suptitle(title, fontsize=13.5, fontweight="bold", y=0.995, va="top")
    fig.text(0.5, 0.915, subtitle, ha="center", va="top", fontsize=10.5)

    handles: list[Artist] = [
        Line2D([], [], ls="none", marker="o", ms=5, color=POINT_GRAY, alpha=0.6),
        Line2D([], [], color=TREND_COLOR, linewidth=2.0),
    ]
    labels = [
        f"protein group (n = {data.n_features:,})",
        f"LOWESS trend (frac = {lowess_frac:g})",
    ]
    return fig, _legend_figure(handles, labels), stats
