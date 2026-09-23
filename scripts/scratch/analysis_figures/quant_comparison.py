"""Quant-comparison figure family: LFQ vs NSAF vs PSM on the DE common feature set.

From-scratch Stage-4 family (no ``lib/figures`` template exists). Reads the
differential-abundance outputs written by
``scripts/scratch/de_raloxifene_vs_control.py`` and renders three figures comparing
the three protein-level quantities on the 1,640-group common feature set under the
PRIMARY (paired) design:

  * :func:`plot_log2fc_scatter` -- per-protein log2FC, pairwise scatter (3 panels,
    equal axes, ``y = x``, Spearman rho per panel);
  * :func:`plot_precision_distributions` -- per-protein residual SD and 95% CI
    half-width distributions per quantity (step histograms, medians marked);
  * :func:`plot_sd_vs_abundance` -- residual SD vs mean log2 abundance, one panel per
    quantity on its own x scale, with an equal-count binned-median trend.

Inputs are validated fail-loud (:func:`load_comparison`): the common-set table's
``log2fc`` / ``residual_sd`` must equal the per-quantity ``<q>_paired.tsv`` tables
exactly (confirming the primary paired design; ``residual_sd`` must also *differ*
from the batch-design table so the check discriminates), and the mean log2 abundance
is joined from those per-quantity tables. Quantity colors are read (never written)
from the ``quantity`` category of ``state/color_registry.json``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common.figures.colors import load_palette, load_registry
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from numpy.typing import NDArray
from scipy import stats

__script_meta__: dict[str, object] = {
    "task": "fig-quant-comparison",
    "kind": "module",
    "provides": [
        "QUANTITIES",
        "QuantitySpec",
        "ComparisonData",
        "load_comparison",
        "quantity_colors",
        "spearman",
        "binned_median",
        "plot_log2fc_scatter",
        "plot_precision_distributions",
        "plot_sd_vs_abundance",
    ],
    "uses": ["common.figures.colors"],
    "seeded_from": None,
    "description": (
        "Quant-comparison figures (LFQ vs NSAF vs PSM) on the DE common feature set, "
        "primary paired design: pairwise log2FC scatter, residual-SD / CI-half-width "
        "distributions, residual SD vs mean log2 abundance. Fail-loud input checks; "
        "registry colors read-only."
    ),
}

FloatArray = NDArray[np.float64]

REGISTRY_CATEGORY = "quantity"
POINT_GRAY = "#404040"  # neutral (non-categorical) point color for the FC scatter
REFERENCE_GRAY = "#808080"  # y = x / zero reference lines
TREND_COLOR = "#000000"  # binned-median trend (palette black; not a category)
TABLE_NAME = "common_set_comparison.tsv"
_EXACT_ATOL = 1e-12
RHO = "\u03c1"  # Greek rho, the conventional Spearman symbol
MINUS = "\u2212"  # typographic minus for on-canvas negative numbers


@dataclass(frozen=True)
class QuantitySpec:
    """One compared quantity: its table prefix, display label and axis wording."""

    key: str  # column prefix in the common-set table / per-quantity table stem
    label: str  # display label == registry value in the ``quantity`` category
    join_column: str  # column in ``<key>_<design>.tsv`` holding the common-set id
    abundance_label: str  # x-axis label for the mean log2 abundance


QUANTITIES: tuple[QuantitySpec, ...] = (
    QuantitySpec("protein", "LFQ", "feature", "mean log2 LFQ intensity (a.u.)"),
    QuantitySpec("nsaf", "NSAF", "first_member_id", "mean log2 NSAF"),
    QuantitySpec("psm_log2", "PSM", "first_member_id", "mean log2 PSM count"),
)

_PER_QUANTITY_COLUMNS = ("log2fc", "ci_low", "ci_high", "residual_sd")


@dataclass(frozen=True)
class ComparisonData:
    """Validated wide table: one row per common-set protein group.

    ``table`` holds, per quantity key ``k``: ``k_log2fc``, ``k_residual_sd``,
    ``k_ci_half_width`` and ``k_mean_log2_abundance`` (all finite floats).
    """

    table: pd.DataFrame
    n_features: int
    design: str

    def column(self, key: str, what: str) -> FloatArray:
        """Return the ``<key>_<what>`` column as a float array."""
        arr: FloatArray = self.table[f"{key}_{what}"].to_numpy(dtype=np.float64)
        return arr


# --------------------------------------------------------------------------- #
# Loading + validation
# --------------------------------------------------------------------------- #
def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required input missing: {path}")
    return pd.read_csv(path, sep="\t")


def _require_columns(frame: pd.DataFrame, columns: Sequence[str], where: Path) -> None:
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise ValueError(f"{where}: missing columns {missing}.")


def _finite(values: pd.Series, name: str) -> FloatArray:
    arr: FloatArray = pd.to_numeric(values, errors="raise").to_numpy(dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise ValueError(
            f"{name}: {int((~np.isfinite(arr)).sum())} non-finite value(s)."
        )
    return arr


def _per_quantity_block(
    de_dir: Path, spec: QuantitySpec, design: str, ids: pd.Series
) -> pd.DataFrame:
    """Rows of ``<key>_<design>.tsv`` aligned to ``ids`` (fail loud on any gap)."""
    path = de_dir / f"{spec.key}_{design}.tsv"
    frame = _read_tsv(path)
    _require_columns(
        frame, [spec.join_column, "log2fc", "residual_sd", "mean_log2_abundance"], path
    )
    if not frame[spec.join_column].is_unique:
        raise ValueError(f"{path}: duplicate ids in {spec.join_column!r}.")
    indexed = frame.set_index(spec.join_column)
    absent = ~ids.isin(indexed.index)
    if absent.any():
        raise ValueError(
            f"{path}: {int(absent.sum())} common-set ids absent (e.g. "
            f"{ids[absent].iloc[0]!r})."
        )
    return indexed.loc[ids.to_numpy()].reset_index()


def load_comparison(
    de_dir: Path,
    *,
    design: str = "paired",
    discriminating_design: str = "batch",
    expected_n: int | None = None,
) -> ComparisonData:
    """Load + validate the common-set table and join mean log2 abundance.

    Confirms that the common-set table is the ``design`` fit: every quantity's
    ``log2fc`` and ``residual_sd`` must equal those of ``<key>_<design>.tsv``
    (tolerance 1e-12), and ``residual_sd`` must differ from
    ``<key>_<discriminating_design>.tsv`` (log2FC is design-invariant here, the
    residual SD is not). Raises on any missing column, non-finite value, duplicate or
    missing id, or mismatch.
    """
    path = de_dir / TABLE_NAME
    table = _read_tsv(path)
    columns = ["feature"] + [
        f"{q.key}_{c}" for q in QUANTITIES for c in _PER_QUANTITY_COLUMNS
    ]
    _require_columns(table, columns, path)
    if not table["feature"].is_unique:
        raise ValueError(f"{path}: duplicate feature ids.")
    if expected_n is not None and len(table) != expected_n:
        raise ValueError(f"{path}: expected {expected_n} rows, got {len(table)}.")
    ids = table["feature"].astype(str)

    out = pd.DataFrame({"feature": ids.to_numpy()})
    for spec in QUANTITIES:
        fc = _finite(table[f"{spec.key}_log2fc"], f"{spec.key}_log2fc")
        sd = _finite(table[f"{spec.key}_residual_sd"], f"{spec.key}_residual_sd")
        lo = _finite(table[f"{spec.key}_ci_low"], f"{spec.key}_ci_low")
        hi = _finite(table[f"{spec.key}_ci_high"], f"{spec.key}_ci_high")
        if np.any(sd < 0) or np.any(hi < lo):
            raise ValueError(f"{spec.key}: negative residual SD or inverted CI.")
        if not np.allclose((lo + hi) / 2.0, fc, rtol=0.0, atol=1e-9):
            raise ValueError(f"{spec.key}: CI not centred on log2FC.")

        block = _per_quantity_block(de_dir, spec, design, ids)
        if not (
            np.allclose(block["log2fc"], fc, rtol=0.0, atol=_EXACT_ATOL)
            and np.allclose(block["residual_sd"], sd, rtol=0.0, atol=_EXACT_ATOL)
        ):
            raise ValueError(
                f"{spec.key}: common-set log2fc/residual_sd do not match "
                f"{spec.key}_{design}.tsv -- table is not the {design!r} fit."
            )
        other = _per_quantity_block(de_dir, spec, discriminating_design, ids)
        if np.allclose(other["residual_sd"], sd, rtol=0.0, atol=_EXACT_ATOL):
            raise ValueError(
                f"{spec.key}: residual_sd identical under {design!r} and "
                f"{discriminating_design!r}; design check is not discriminating."
            )
        out[f"{spec.key}_log2fc"] = fc
        out[f"{spec.key}_residual_sd"] = sd
        out[f"{spec.key}_ci_half_width"] = (hi - lo) / 2.0
        out[f"{spec.key}_mean_log2_abundance"] = _finite(
            block["mean_log2_abundance"], f"{spec.key}_mean_log2_abundance"
        )
    return ComparisonData(table=out, n_features=len(out), design=design)


def quantity_colors(registry_path: Path) -> dict[str, str]:
    """Registry colors for the quantity labels (read-only; fail loud if absent)."""
    registry = load_registry(registry_path)
    palette = load_palette(registry_path)
    entry = registry.get(REGISTRY_CATEGORY)
    if not isinstance(entry, dict) or not isinstance(entry.get("values"), dict):
        raise ValueError(
            f"{registry_path}: no {REGISTRY_CATEGORY!r} category with a values map."
        )
    values: dict[str, object] = entry["values"]
    colors: dict[str, str] = {}
    for spec in QUANTITIES:
        color = values.get(spec.label)
        if not isinstance(color, str):
            raise ValueError(
                f"{registry_path}: {REGISTRY_CATEGORY}/{spec.label} has no color."
            )
        if color not in palette.colors:
            raise ValueError(
                f"{REGISTRY_CATEGORY}/{spec.label}={color} not in palette."
            )
        colors[spec.label] = color
    if len(set(colors.values())) != len(colors):
        raise ValueError(f"{REGISTRY_CATEGORY}: quantities share a color: {colors}.")
    return colors


# --------------------------------------------------------------------------- #
# Statistics shown on the canvas
# --------------------------------------------------------------------------- #
def spearman(x: FloatArray, y: FloatArray) -> float:
    """Spearman rank correlation (average ranks for ties)."""
    if x.shape != y.shape or x.size < 3:
        raise ValueError("spearman needs two equal-length arrays of >= 3 values.")
    rho = float(stats.spearmanr(x, y).statistic)
    if not math.isfinite(rho):
        raise ValueError("spearman: undefined (constant input).")
    return rho


def binned_median(
    x: FloatArray, y: FloatArray, n_bins: int
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Equal-count bins of ``x``; return (median x, median y, bin count) per bin."""
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("binned_median: x and y must be equal-length 1-D arrays.")
    if n_bins < 2 or x.size < 2 * n_bins:
        raise ValueError("binned_median: need >= 2 bins and >= 2 points per bin.")
    order = np.argsort(x, kind="stable")
    chunks = np.array_split(order, n_bins)
    mx = np.array([np.median(x[c]) for c in chunks], dtype=np.float64)
    my = np.array([np.median(y[c]) for c in chunks], dtype=np.float64)
    counts = np.array([c.size for c in chunks], dtype=np.float64)
    return mx, my, counts


# --------------------------------------------------------------------------- #
# Shared drawing helpers
# --------------------------------------------------------------------------- #
def _titles(fig: Figure, title: str, subtitle: str, top: float) -> None:
    fig.text(
        0.5,
        top + 0.085,
        title,
        ha="center",
        va="bottom",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(0.5, top + 0.035, subtitle, ha="center", va="bottom", fontsize=10.5)


def _legend_figure(handles: Sequence[Line2D], labels: Sequence[str]) -> Figure:
    fig = plt.figure(figsize=(3.2, 0.32 * len(handles) + 0.2))
    fig.legend(handles, labels, loc="center", frameon=False, handlelength=2.2)
    return fig


def _nice_limit(value: float, step: float) -> float:
    return math.ceil(value / step) * step


def _nice_step(span: float) -> float:
    """A tick step giving ~4-6 ticks over ``span``."""
    for step in (0.1, 0.2, 0.25, 0.5, 1.0, 2.0, 2.5, 5.0):
        if span / step <= 6:
            return step
    return 10.0


# --------------------------------------------------------------------------- #
# Figure 1 -- pairwise log2FC scatter
# --------------------------------------------------------------------------- #
PAIRS: tuple[tuple[str, str], ...] = (
    ("protein", "nsaf"),
    ("protein", "psm_log2"),
    ("nsaf", "psm_log2"),
)


def _spec(key: str) -> QuantitySpec:
    for spec in QUANTITIES:
        if spec.key == key:
            return spec
    raise KeyError(key)


def plot_log2fc_scatter(
    data: ComparisonData, *, title: str, subtitle: str
) -> tuple[Figure, Figure, dict[str, float]]:
    """Three equal-axis panels of per-protein log2FC, ``y = x`` and Spearman rho."""
    fcs = {q.key: data.column(q.key, "log2fc") for q in QUANTITIES}
    lim = _nice_limit(max(float(np.max(np.abs(v))) for v in fcs.values()) * 1.04, 0.25)
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 5.0))
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.12, top=0.84, wspace=0.28)
    rhos: dict[str, float] = {}
    for ax, (kx, ky) in zip(axes, PAIRS, strict=True):
        x, y = fcs[kx], fcs[ky]
        rho = spearman(x, y)
        rhos[f"{_spec(kx).label}-{_spec(ky).label}"] = rho
        ax.axhline(0.0, color=REFERENCE_GRAY, lw=0.6, alpha=0.5, zorder=1)
        ax.axvline(0.0, color=REFERENCE_GRAY, lw=0.6, alpha=0.5, zorder=1)
        ax.plot(
            [-lim, lim], [-lim, lim], ls="--", color=REFERENCE_GRAY, lw=1.1, zorder=2
        )
        ax.scatter(x, y, s=6, color=POINT_GRAY, alpha=0.35, linewidths=0, zorder=3)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal", adjustable="box")
        step = _nice_step(2 * lim)
        half = math.floor(lim / step + 1e-9) * step
        ticks = np.arange(-half, half + step / 2, step)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlabel(f"{_spec(kx).label} log2 fold change")
        ax.set_ylabel(f"{_spec(ky).label} log2 fold change")
        ax.text(
            0.04,
            0.96,
            f"Spearman {RHO} = {rho:.3f}".replace("-", MINUS),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=11,
        )
    _titles(fig, title, subtitle, 0.84)
    handles = [
        Line2D([], [], ls="none", marker="o", ms=5, color=POINT_GRAY, alpha=0.6),
        Line2D([], [], ls="--", color=REFERENCE_GRAY, lw=1.1),
    ]
    labels = [f"protein group (n = {data.n_features:,})", "y = x"]
    return fig, _legend_figure(handles, labels), rhos


# --------------------------------------------------------------------------- #
# Figure 2 -- precision distributions
# --------------------------------------------------------------------------- #
def _distribution_panel(
    ax: Axes,
    data: ComparisonData,
    what: str,
    colors: dict[str, str],
    xlabel: str,
) -> dict[str, float]:
    arrays = {q.key: data.column(q.key, what) for q in QUANTITIES}
    xmax = _nice_limit(max(float(np.max(a)) for a in arrays.values()), 0.1)
    bins = [float(b) for b in np.linspace(0.0, xmax, round(xmax / 0.025) + 1)]
    medians: dict[str, float] = {}
    for spec in QUANTITIES:
        values = arrays[spec.key]
        color = colors[spec.label]
        ax.hist(values, bins=bins, density=True, histtype="step", lw=1.6, color=color)
        med = float(np.median(values))
        medians[spec.label] = med
        ax.axvline(med, color=color, ls="--", lw=1.3)
    ax.set_xlim(0.0, xmax)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    ax.text(
        0.97,
        0.96,
        "median",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10.5,
        fontweight="bold",
    )
    for i, spec in enumerate(QUANTITIES):
        ax.text(
            0.97,
            0.96 - 0.075 * (i + 1),
            f"{spec.label} {medians[spec.label]:.3f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=10.5,
            color=colors[spec.label],
        )
    return medians


def plot_precision_distributions(
    data: ComparisonData, colors: dict[str, str], *, title: str, subtitle: str
) -> tuple[Figure, Figure, dict[str, dict[str, float]]]:
    """Residual SD and 95% CI half-width distributions per quantity, medians marked."""
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.14, top=0.82, wspace=0.2)
    medians = {
        "residual_sd": _distribution_panel(
            axes[0], data, "residual_sd", colors, "per-protein residual SD (log2 units)"
        ),
        "ci_half_width": _distribution_panel(
            axes[1],
            data,
            "ci_half_width",
            colors,
            "per-protein 95% CI half-width of log2FC (log2 units)",
        ),
    }
    _titles(fig, title, subtitle, 0.82)
    handles = [Line2D([], [], color=colors[q.label], lw=1.8) for q in QUANTITIES]
    handles.append(Line2D([], [], color=TREND_COLOR, ls="--", lw=1.3))
    labels = [q.label for q in QUANTITIES] + ["median (dashed, per quantity color)"]
    return fig, _legend_figure(handles, labels), medians


# --------------------------------------------------------------------------- #
# Figure 3 -- residual SD vs mean abundance
# --------------------------------------------------------------------------- #
def plot_sd_vs_abundance(
    data: ComparisonData,
    colors: dict[str, str],
    *,
    title: str,
    subtitle: str,
    n_bins: int = 10,
) -> tuple[Figure, Figure, dict[str, dict[str, object]]]:
    """Residual SD vs mean log2 abundance per quantity + binned-median trend."""
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.6), sharey=True)
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.14, top=0.8, wspace=0.12)
    ymax = _nice_limit(
        max(float(np.max(data.column(q.key, "residual_sd"))) for q in QUANTITIES), 0.1
    )
    trend: dict[str, dict[str, object]] = {}
    for ax, spec in zip(axes, QUANTITIES, strict=True):
        x = data.column(spec.key, "mean_log2_abundance")
        y = data.column(spec.key, "residual_sd")
        rho = spearman(x, y)
        mx, my, counts = binned_median(x, y, n_bins)
        trend[spec.label] = {
            "spearman_abundance_vs_sd": rho,
            "bin_median_abundance": [round(float(v), 4) for v in mx],
            "bin_median_residual_sd": [round(float(v), 4) for v in my],
            "bin_counts": [int(v) for v in counts],
        }
        ax.scatter(
            x, y, s=6, color=colors[spec.label], alpha=0.35, linewidths=0, zorder=2
        )
        ax.plot(
            mx,
            my,
            color=TREND_COLOR,
            lw=1.8,
            marker="o",
            ms=5,
            mfc=TREND_COLOR,
            mec="white",
            mew=0.8,
            zorder=4,
        )
        ax.set_ylim(0.0, ymax)
        ax.set_xlabel(spec.abundance_label)
        ax.set_title(spec.label, color=colors[spec.label], fontsize=13)
        ax.text(
            0.97,
            0.96,
            f"Spearman {RHO} = {rho:.2f}".replace("-", MINUS),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=10.5,
        )
    axes[0].set_ylabel("per-protein residual SD (log2 units)")
    _titles(fig, title, subtitle, 0.84)
    handles = [
        Line2D([], [], ls="none", marker="o", ms=5, color=colors[q.label])
        for q in QUANTITIES
    ]
    handles.append(
        Line2D([], [], color=TREND_COLOR, lw=1.8, marker="o", ms=5, mec="white")
    )
    labels = [f"{q.label} protein group" for q in QUANTITIES]
    labels.append(f"binned median ({n_bins} equal-count bins)")
    return fig, _legend_figure(handles, labels), trend
