"""Abundance-agreement figure family: FlashLFQ LFQ vs Limelight PSM / NSAF.

From-scratch Stage-4 family (no ``lib/figures`` template fits). Reads the tables
written by ``scripts/scratch/abundance_agreement.py`` to
``results/abundance-agreement/``
(indexed, with their columns, in that directory's ``summary.json``) and renders:

  finding 0010 (``figures/analysis/quant-comparison/abundance-agreement/``)
    * :func:`plot_density`             -- per-cell 2-D density, LFQ vs log2 PSM / NSAF;
    * :func:`plot_per_sample`          -- per-run Spearman with CIs, acquisition order;
    * :func:`plot_slopes`              -- slope estimators + LOWESS local slope;
    * :func:`plot_by_psm_count`        -- spread at fixed PSM count, by LFQ decile;
    * :func:`plot_residual_drivers`    -- partial-residual plots of the two drivers;
    * :func:`plot_detection_classes`   -- both / PSM-only / LFQ-only per run + ECDF;
    * :func:`plot_no_lfq_reasons`      -- the 897 PSM-but-never-LFQ proteins;
    * :func:`plot_within_protein`      -- per-protein run-to-run tracking vs nulls;
  finding 0011 (``figures/analysis/quant-comparison/contaminant-sharing/``)
    * :func:`plot_contaminant_sharing` -- real no-LFQ proteins that would regain
      unique peptides without the contaminant-list copies (table from
      :func:`contaminant_sharing_table`), with the contaminant copy's LFQ run count.

Every plotted number is read from the tables; each on-canvas statistic that
``summary.json`` also records is cross-checked against it (:func:`check`, fail loud).
Colors: the ``quantity`` and ``batch`` registry categories are read (never written);
figure-local labels (detection class, no-LFQ reason) get palette colors through
``assign_colors(..., persist=False)`` so the registry file is not modified.
"""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common.figures.colors import assign_colors, load_palette, load_registry
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
from numpy.typing import NDArray
from scipy import stats

__script_meta__: dict[str, object] = {
    "task": "fig-abundance-agreement",
    "kind": "module",
    "provides": [
        "AgreementData",
        "load_agreement",
        "check",
        "plot_density",
        "plot_per_sample",
        "plot_slopes",
        "plot_by_psm_count",
        "plot_residual_drivers",
        "plot_detection_classes",
        "plot_no_lfq_reasons",
        "plot_within_protein",
        "contaminant_sharing_table",
        "plot_contaminant_sharing",
    ],
    "uses": ["common.figures.colors"],
    "seeded_from": None,
    "description": (
        "Findings 0010/0011 figures: LFQ vs PSM/NSAF agreement (density, per-run "
        "correlation, slope estimators, count-noise breakdown, residual drivers, "
        "detection classes, no-LFQ reasons, within-protein tracking) and contaminant "
        "peptide sharing. Reads results/abundance-agreement tables; on-canvas "
        "numbers cross-checked against summary.json; registry read-only."
    ),
}

FloatArray = NDArray[np.float64]
Result = tuple[Figure, Figure, dict[str, Any]]

RHO = "\u03c1"
MINUS = "\u2212"
LE = "≤"
GE = "≥"
DOT = " · "
NEUTRAL = "#404040"  # non-categorical marks (points, reference text)
REF_GRAY = "#808080"  # reference lines
BAND_GRAY = "#d9d9d9"  # null / CI bands
LFQ_AXIS = "log2 LFQ intensity, median-normalized (a.u.)"
PROCESSING = (
    "LFQ: FlashLFQ, median-normalized log2" + DOT + "PSM, NSAF: Limelight, log2"
)
_ATOL = 6e-6  # summary.json rounds to 6 decimals

# Explicit Okabe-Ito members used for figure-local fit/trend lines (not categories);
# asserted to be palette members at load time.
FIT_LOWESS = "#D55E00"
FIT_LINE = "#56B4E9"
FIT_BLACK = "#000000"

DENSITY_CMAP = LinearSegmentedColormap.from_list(
    "density_greys", plt.get_cmap("Greys")(np.linspace(0.18, 0.92, 256))
)

PSM_BINS = ("1", "2-3", "4-7", "8-15", "16-31", ">=32")
# Slope estimators shown (the replicate-noise-only estimators -- Deming, attenuation-
# corrected OLS -- are deliberately omitted: their error model does not hold here).
ESTIMATORS = (
    ("ols_y_on_x", "OLS (y on x)"),
    ("sma", "SMA"),
    ("orthogonal", "orthogonal"),
    ("inverse_ols_x_on_y", "inverse OLS (x on y)"),
)
LFQ_SCALES = (
    ("protein", "protein LFQ", "o"),
    ("per_peptide", "LFQ per unique peptide", "D"),
)

DETECTION_CLASSES = (
    "both",
    "PSM only, no LFQ in any run",
    "PSM only, LFQ in other runs",
    "LFQ only (all MBR)",
)
REASONS = (
    ("shared_peptides_only", "shared peptides only"),
    ("unique_peptides_never_quantified", "unique peptides, never quantified"),
    ("flashlfq_nan_not_quantifiable", "FlashLFQ NaN token"),
)


# --------------------------------------------------------------------------- #
# Loading + checks
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class AgreementData:
    """All abundance-agreement tables, validated against ``summary.json``."""

    result_dir: Path
    summary: dict[str, Any]
    tables: dict[str, pd.DataFrame]
    samples: tuple[str, ...]  # acquisition order
    batch_of: dict[str, str]
    condition_of: dict[str, str]

    def t(self, name: str) -> pd.DataFrame:
        """Return the table ``<name>.tsv``."""
        return self.tables[name]


def load_agreement(result_dir: Path, acquisition_order: Sequence[str]) -> AgreementData:
    """Read ``summary.json`` and every table it indexes; fail loud on any drift.

    Each table must exist with exactly the columns ``summary.json`` lists, and the
    summary's sample list must equal ``acquisition_order`` (from samples.tsv).
    """
    summary_path = result_dir / "summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"Missing {summary_path}")
    summary: dict[str, Any] = json.loads(summary_path.read_text())
    tables: dict[str, pd.DataFrame] = {}
    for fname, columns in summary["tables"].items():
        path = result_dir / fname
        if not path.is_file():
            raise FileNotFoundError(f"summary.json lists {fname} but {path} is absent.")
        frame = pd.read_csv(path, sep="\t")
        if list(frame.columns) != list(columns):
            raise ValueError(
                f"{path}: columns {list(frame.columns)} != summary {columns}"
            )
        tables[fname.removesuffix(".tsv")] = frame
    samples = tuple(str(s) for s in summary["samples"])
    if samples != tuple(acquisition_order):
        raise ValueError(
            f"summary samples {samples} != acquisition order {tuple(acquisition_order)}"
        )
    det = tables["detection_per_sample"]
    if tuple(det.sample_id) != samples:
        raise ValueError("detection_per_sample.tsv is not in acquisition order.")
    batch_of = dict(zip(det.sample_id, det.batch, strict=True))
    condition_of = dict(zip(det.sample_id, det.condition, strict=True))
    palette = load_palette()
    for color in (FIT_LOWESS, FIT_LINE, FIT_BLACK):
        if color not in palette.colors:
            raise ValueError(f"{color} is not an Okabe-Ito palette member.")
    return AgreementData(result_dir, summary, tables, samples, batch_of, condition_of)


def check(label: str, computed: float, reported: object, atol: float = _ATOL) -> None:
    """Fail loud unless ``computed`` matches the summary.json value ``reported``."""
    if not isinstance(reported, int | float) or isinstance(reported, bool):
        raise ValueError(f"summary.json {label}: not a number ({reported!r}).")
    if not math.isclose(computed, float(reported), abs_tol=atol, rel_tol=0.0):
        raise ValueError(
            f"{label}: figure value {computed!r} != summary.json {reported!r}."
        )


def registry_colors(
    category: str, values: Sequence[str], registry_path: Path
) -> dict[str, str]:
    """Read (never write) registry colors for ``values`` of ``category``."""
    entry = load_registry(registry_path).get(category)
    if not isinstance(entry, dict) or not isinstance(entry.get("values"), dict):
        raise ValueError(f"{registry_path}: no {category!r} category.")
    out: dict[str, str] = {}
    for value in values:
        color = entry["values"].get(value)
        if not isinstance(color, str):
            raise ValueError(f"{registry_path}: {category}/{value} has no color.")
        out[value] = color
    return out


def local_colors(
    category: str, values: Sequence[str], registry_path: Path
) -> dict[str, str]:
    """Palette colors for figure-local labels, NOT persisted to the registry."""
    return assign_colors(
        category, list(values), registry_path=registry_path, persist=False
    )


# --------------------------------------------------------------------------- #
# Small drawing helpers
# --------------------------------------------------------------------------- #
def _titles(fig: Figure, title: str, subtitle: str, top: float) -> None:
    height = fig.get_figheight()
    fig.text(
        0.5,
        top + 0.98 / height,
        title,
        ha="center",
        va="bottom",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5, top + 0.58 / height, subtitle, ha="center", va="bottom", fontsize=10.5
    )


def _panel(ax: Axes, letter: str, title: str) -> None:
    ax.set_title(title, fontsize=11.5, fontweight="normal", loc="left", pad=8)
    ax.text(
        -0.02,
        1.02,
        letter,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="right",
        va="bottom",
    )


def _stat_box(ax: Axes, text: str, loc: str = "upper left") -> None:
    x, ha = (0.03, "left") if "left" in loc else (0.97, "right")
    y, va = (0.97, "top") if "upper" in loc else (0.03, "bottom")
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha=ha,
        va=va,
        fontsize=9.5,
        linespacing=1.35,
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "none", "alpha": 0.85},
    )


def _legend_figure(
    handles: Sequence[Any], labels: Sequence[str], ncol: int = 1, width: float = 4.0
) -> Figure:
    rows = math.ceil(len(handles) / ncol)
    fig = plt.figure(figsize=(width, 0.3 * rows + 0.25))
    fig.legend(
        handles, labels, loc="center", frameon=False, ncol=ncol, handlelength=2.2
    )
    return fig


def _plain_log(ax: Axes) -> None:
    """Plain-number tick labels (1, 10, 100) on a log x axis."""
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:,.0f}"))


def _batch_separator(ax: Axes, data: AgreementData) -> None:
    batches = [data.batch_of[s] for s in data.samples]
    for i in range(1, len(batches)):
        if batches[i] != batches[i - 1]:
            ax.axvline(i - 0.5, color=REF_GRAY, lw=0.9, ls=":")


def _fmt_ci(v: float, lo: float, hi: float, digits: int = 2) -> str:
    def f(z: float) -> str:
        return f"{z:.{digits}f}".replace("-", MINUS)

    return f"{f(v)} [{f(lo)}, {f(hi)}]"


def _row(frame: pd.DataFrame, **match: object) -> pd.Series:
    mask = np.ones(len(frame), dtype=bool)
    for key, value in match.items():
        mask &= (frame[key] == value).to_numpy()
    if int(mask.sum()) != 1:
        raise ValueError(f"expected exactly one row for {match}; got {int(mask.sum())}")
    row: pd.Series = frame.loc[mask].iloc[0]
    return row


def _f(value: object) -> float:
    return float(value)  # type: ignore[arg-type]


def _short(entry: str) -> str:
    return entry.removesuffix("_HUMAN")


def _sample_axis(ax: Axes, data: AgreementData) -> None:
    """Acquisition-order x axis with a dotted separator between batches."""
    ax.set_xticks(range(len(data.samples)), list(data.samples), rotation=45, ha="right")
    ax.set_xlim(-0.6, len(data.samples) - 0.4)
    _batch_separator(ax, data)
    ax.set_xlabel("run (acquisition order)")


# --------------------------------------------------------------------------- #
# 1. Density: LFQ vs log2 PSM / log2 NSAF
# --------------------------------------------------------------------------- #
def plot_density(data: AgreementData, qcol: dict[str, str]) -> Result:
    """Two hexbin panels of per-cell LFQ vs log2 PSM / NSAF, LOWESS + SMA lines."""
    del qcol  # quantities are identified by panel, not color
    cells = data.t("cells")
    common = cells[cells["in_common"]]
    agree = data.summary["1_agreement"]
    n_cells, n_prot = len(common), int(common.protein_id.nunique())
    if n_cells != agree["n_common_cells"]:
        raise ValueError(f"common cells {n_cells} != {agree['n_common_cells']}")
    x = common.lfq.to_numpy(dtype=np.float64)
    trend = data.t("trend_lowess")
    slopes = data.t("slopes")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4))
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.12, top=0.83, wspace=0.2)
    specs = (
        ("PSM", "log2_psm", "pooled_psm", "lfq_vs_psm", "log2 PSM count per run"),
        ("NSAF", "log2_nsaf", "pooled_nsaf", "lfq_vs_nsaf", "log2 NSAF per run"),
    )
    colls: list[PolyCollection] = []
    numbers: dict[str, Any] = {"n_cells": n_cells, "n_proteins": n_prot}
    xlo, xhi = math.floor(x.min()), math.ceil(x.max())
    for ax, (q, ycol, tname, comp, ylabel), letter in zip(
        axes, specs, "ab", strict=True
    ):
        y = common[ycol].to_numpy(dtype=np.float64)
        ylo, yhi = math.floor(y.min()) - 0.5, math.ceil(y.max()) + 0.5
        coll = ax.hexbin(
            x,
            y,
            gridsize=(60, 42),
            extent=(xlo, xhi, ylo, yhi),
            mincnt=1,
            cmap=DENSITY_CMAP,
            linewidths=0.0,
        )
        colls.append(coll)
        g = trend[(trend.trend == tname) & (trend.scope == "pooled_cells_all_common")]
        if len(g) != 200 or set(g.x_quantity) != {"lfq"}:
            raise ValueError(f"LOWESS grid {tname}: unexpected shape/axis.")
        ax.plot(g.x, g.fitted, color=FIT_LOWESS, lw=2.2, solid_capstyle="round")
        sma = _row(
            slopes, comparison=comp, set="pooled_cells_all_common", estimator="sma"
        )
        xs = np.array([x.min(), x.max()])
        ax.plot(
            xs,
            _f(sma.intercept) + _f(sma.slope) * xs,
            color=FIT_LINE,
            lw=2.0,
            ls=(0, (5, 2.5)),
        )
        # Recompute what the canvas states from the plotted cells.
        rho = float(stats.spearmanr(x, y).statistic)
        pooled = agree["correlations"][f"{comp}|all_common|spearman"]["pooled"]
        check(f"{comp} pooled spearman", rho, pooled[0])
        r = float(np.corrcoef(x, y)[0, 1])
        sma_calc = math.copysign(float(np.std(y, ddof=1) / np.std(x, ddof=1)), r)
        check(f"{comp} SMA slope", sma_calc, _f(sma.slope), atol=1e-5)
        check(
            f"{comp} SMA slope summary",
            _f(sma.slope),
            data.summary["2_trend"]["slopes"][f"{comp}|pooled_cells_all_common"]["sma"][
                0
            ],
        )
        _stat_box(
            ax,
            f"Spearman {RHO} = {_fmt_ci(rho, pooled[1], pooled[2], 3)}\n"
            f"SMA slope = {_f(sma.slope):.2f} (descriptive)",
        )
        numbers[q] = {
            "spearman": [round(rho, 4), pooled[1], pooled[2]],
            "sma_slope": [
                round(_f(sma.slope), 4),
                round(_f(sma.ci_low), 4),
                round(_f(sma.ci_high), 4),
            ],
            "sma_intercept": round(_f(sma.intercept), 4),
        }
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
        ax.set_xlabel(LFQ_AXIS)
        ax.set_ylabel(ylabel)
        _panel(ax, letter, f"LFQ vs {q}")
    vmax = max(float(np.max(c.get_array())) for c in colls)  # type: ignore[arg-type]
    norm = LogNorm(vmin=1.0, vmax=vmax)
    for c in colls:
        c.set_norm(norm)
    numbers["max_cells_per_hexagon"] = int(vmax)
    _titles(
        fig,
        "Per-run protein abundance: FlashLFQ LFQ vs Limelight spectral counts",
        f"{n_cells:,} run\u00d7protein cells quantified by both ({n_prot:,} proteins,"
        f" 8 runs){DOT}{PROCESSING}",
        0.83,
    )

    legend = plt.figure(figsize=(4.6, 1.55))
    cax = legend.add_axes((0.08, 0.66, 0.84, 0.13))
    cb = legend.colorbar(
        ScalarMappable(norm=norm, cmap=DENSITY_CMAP), cax=cax, orientation="horizontal"
    )
    cb.set_label("cells per hexagon (log scale)", fontsize=10)
    cb.set_ticks(
        [t for t in (1, 3, 10, 30, 100, 300) if t <= vmax],
        labels=[str(t) for t in (1, 3, 10, 30, 100, 300) if t <= vmax],
    )
    legend.legend(
        [
            Line2D([], [], color=FIT_LOWESS, lw=2.2),
            Line2D([], [], color=FIT_LINE, lw=2.0, ls=(0, (5, 2.5))),
        ],
        ["LOWESS trend (frac 0.3)", "SMA (symmetric descriptive fit)"],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=2,
        frameon=False,
    )
    return fig, legend, numbers


# --------------------------------------------------------------------------- #
# 2. Per-sample correlations
# --------------------------------------------------------------------------- #
def plot_per_sample(data: AgreementData, bcol: dict[str, str]) -> Result:
    """Per-run Spearman LFQ-PSM / LFQ-NSAF with bootstrap CIs, by acquisition order."""
    corr = data.t("correlations")
    agree = data.summary["1_agreement"]["correlations"]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.4), sharex=True)
    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.14, top=0.86, hspace=0.28)
    markers = {"control": "o", "raloxifene-d0": "^"}
    numbers: dict[str, Any] = {}
    for ax, (comp, q), letter in zip(
        axes, (("lfq_vs_psm", "PSM"), ("lfq_vs_nsaf", "NSAF")), "ab", strict=True
    ):
        sub = corr[
            (corr.comparison == comp)
            & (corr.subset == "all_common")
            & (corr.method == "spearman")
        ]
        per = (
            sub[sub.scope == "per_sample"]
            .set_index("sample_id")
            .loc[list(data.samples)]
        )
        pooled = _row(sub, scope="pooled_cells")
        block = agree[f"{comp}|all_common|spearman"]
        check(f"{comp} pooled", _f(pooled.r), block["pooled"][0])
        check(
            f"{comp} per-sample min", float(per.r.min()), block["per_sample_min_max"][0]
        )
        check(
            f"{comp} per-sample max", float(per.r.max()), block["per_sample_min_max"][1]
        )
        ax.axhspan(_f(pooled.ci_low), _f(pooled.ci_high), color=BAND_GRAY, lw=0)
        ax.axhline(_f(pooled.r), color=REF_GRAY, lw=1.2, ls="--")
        for i, sid in enumerate(data.samples):
            row = per.loc[sid]
            color = bcol[data.batch_of[sid]]
            ax.errorbar(
                i,
                row.r,
                yerr=[[row.r - row.ci_low], [row.ci_high - row.r]],
                fmt="none",
                ecolor=color,
                elinewidth=1.6,
                capsize=3,
            )
            ax.plot(
                i,
                row.r,
                marker=markers[data.condition_of[sid]],
                ms=8,
                color=color,
                mec="black",
                mew=0.6,
                ls="none",
            )
        mid = _f(pooled.r)
        ax.set_ylim(mid - 0.05, mid + 0.05)
        if letter == "a":
            _batch_separator(ax, data)
        ax.set_ylabel(f"Spearman {RHO}, LFQ vs {q}")
        _panel(
            ax,
            letter,
            f"LFQ vs {q}: pooled {RHO} = {_f(pooled.r):.3f}{DOT}per run "
            f"{per.r.min():.3f}\u2013{per.r.max():.3f}",
        )
        numbers[q] = {
            "pooled": [
                round(_f(pooled.r), 4),
                round(_f(pooled.ci_low), 4),
                round(_f(pooled.ci_high), 4),
            ],
            "per_sample": {s: round(float(per.loc[s].r), 4) for s in data.samples},
            "n_cells_min_max": [int(per.n_cells.min()), int(per.n_cells.max())],
        }
    _sample_axis(axes[1], data)
    _titles(
        fig,
        "Per-run agreement of LFQ with spectral counts",
        f"Spearman {RHO} over each run's proteins quantified by both "
        f"(n = {numbers['PSM']['n_cells_min_max'][0]:,}\u2013"
        f"{numbers['PSM']['n_cells_min_max'][1]:,}){DOT}95% bootstrap CI",
        0.86,
    )
    handles: list[Any] = [
        Patch(color=bcol[b]) for b in sorted(set(data.batch_of.values()))
    ]
    labels = [f"batch {b[1:]}" for b in sorted(set(data.batch_of.values()))]
    handles += [
        Line2D([], [], marker=m, color=NEUTRAL, ls="none", ms=8, mec="black", mew=0.6)
        for m in markers.values()
    ]
    labels += list(markers)
    handles += [Line2D([], [], color=REF_GRAY, ls="--", lw=1.2), Patch(color=BAND_GRAY)]
    labels += [f"pooled {RHO} (all runs)", "pooled 95% CI"]
    return fig, _legend_figure(handles, labels, ncol=2, width=5.4), numbers


# --------------------------------------------------------------------------- #
# 3. Slope estimators on two LFQ scales + LOWESS local slope
# --------------------------------------------------------------------------- #
def slope_estimates(x: FloatArray, y: FloatArray) -> dict[str, float]:
    """OLS (y on x), SMA, orthogonal (Deming, lambda = 1) and inverse OLS slopes."""
    sxx = float(np.var(x, ddof=1))
    syy = float(np.var(y, ddof=1))
    sxy = float(np.cov(x, y)[0, 1])
    r = sxy / math.sqrt(sxx * syy)
    return {
        "ols_y_on_x": sxy / sxx,
        "sma": math.copysign(math.sqrt(syy / sxx), r),
        "orthogonal": (syy - sxx + math.sqrt((syy - sxx) ** 2 + 4 * sxy**2))
        / (2 * sxy),
        "inverse_ols_x_on_y": syy / sxy,
        "r_pearson": r,
    }


def plot_slopes(data: AgreementData, qcol: dict[str, str]) -> Result:
    """Slope of log2 PSM / NSAF on LFQ by estimator and LFQ scale (protein means,
    point estimates), plus the LOWESS local slope along the LFQ range."""
    slopes = data.t("slopes")
    props = data.t("protein_properties")
    sub = props[props.complete8]
    n = len(sub)
    if n != data.summary["4_residuals"]["n_complete8"]:
        raise ValueError("complete8 protein count mismatch.")
    lfq = sub.mean_lfq.to_numpy(dtype=np.float64)
    lpep = sub.log2_n_unique_peptides.to_numpy(dtype=np.float64)
    if not np.all(np.isfinite(lpep)) or int(sub.n_unique_peptides.min()) < 1:
        raise ValueError("complete8 protein without a unique peptide.")
    xs = {"protein": lfq, "per_peptide": lfq - lpep}
    fig = plt.figure(figsize=(14.0, 5.4))
    outer = fig.add_gridspec(
        1,
        2,
        width_ratios=(2.1, 1.15),
        wspace=0.2,
        left=0.12,
        right=0.985,
        bottom=0.13,
        top=0.8,
    )
    inner = outer[0].subgridspec(1, 2, wspace=0.1)
    ax_p = fig.add_subplot(inner[0])
    ax_n = fig.add_subplot(inner[1], sharey=ax_p)
    ax_l = fig.add_subplot(outer[1])
    offsets = {"protein": -0.14, "per_peptide": 0.14}
    numbers: dict[str, Any] = {"n_proteins": n}
    for ax, (comp, q, ycol), letter in zip(
        (ax_p, ax_n),
        (
            ("lfq_vs_psm", "PSM", "mean_log2_psm"),
            ("lfq_vs_nsaf", "NSAF", "mean_log2_nsaf"),
        ),
        "ab",
        strict=True,
    ):
        y = sub[ycol].to_numpy(dtype=np.float64)
        ax.axvline(1.0, color=REF_GRAY, lw=1.1, ls="--")
        numbers[q] = {}
        for scale, _label, marker in LFQ_SCALES:
            est = slope_estimates(xs[scale], y)
            numbers[q][scale] = {k: round(v, 4) for k, v in est.items()}
            for i, (key, _lab) in enumerate(ESTIMATORS):
                if scale == "protein":  # the analysis table covers the protein scale
                    check(
                        f"{comp} {key}",
                        est[key],
                        _f(
                            _row(
                                slopes,
                                comparison=comp,
                                set="protein_mean_complete8",
                                estimator=key,
                            ).slope
                        ),
                        1e-5,
                    )
                yy = i + offsets[scale]
                filled = scale == "protein"
                ax.plot(
                    est[key],
                    yy,
                    marker=marker,
                    ms=8,
                    ls="none",
                    color=qcol[q],
                    mfc=qcol[q] if filled else "white",
                    mec=qcol[q] if not filled else "black",
                    mew=1.6 if not filled else 0.5,
                )
                ax.text(
                    est[key] + 0.05,
                    yy,
                    f"{est[key]:.2f}",
                    va="center",
                    fontsize=8.5,
                    color=NEUTRAL,
                )
        ax.set_xlim(0.5, 2.4)
        ax.set_xlabel(f"slope, log2 {q} per log2 LFQ unit")
        _panel(ax, letter, f"LFQ vs {q}")
    ax_p.set_yticks(range(len(ESTIMATORS)), [lab for _e, lab in ESTIMATORS])
    ax_p.set_ylim(len(ESTIMATORS) - 0.5, -0.5)
    plt.setp(ax_n.get_yticklabels(), visible=False)

    xm = np.sort(lfq)
    trend = data.t("trend_lowess")
    ax_l.axhline(1.0, color=REF_GRAY, lw=1.1, ls="--")
    for tname, q in (("mean_psm", "PSM"), ("mean_nsaf", "NSAF")):
        g = trend[(trend.trend == tname) & (trend.scope == "protein_mean_complete8")]
        gx = g.x.to_numpy(dtype=np.float64)
        pct = 100.0 * np.searchsorted(xm, gx, side="right") / xm.size
        ax_l.plot(pct, g.local_slope, color=qcol[q], lw=2.2)
        if q == "PSM":
            for key, qq in (
                ("10%", 0.1),
                ("25%", 0.25),
                ("50%", 0.5),
                ("75%", 0.75),
                ("90%", 0.9),
            ):
                val = float(np.interp(np.quantile(xm, qq), gx, g.local_slope))
                check(
                    f"local slope {key}",
                    val,
                    data.summary["2_trend"]["lowess_local_slope_mean_psm_by_lfq"][key],
                    1e-4,
                )  # grid TSV carries 6 significant digits
            numbers["local_slope_psm_by_pct"] = data.summary["2_trend"][
                "lowess_local_slope_mean_psm_by_lfq"
            ]
    ax_l.set_xlim(0, 100)
    ax_l.set_xlabel("protein-mean LFQ percentile (%)")
    ax_l.set_ylabel("local slope of LOWESS trend (protein LFQ scale)")
    _panel(ax_l, "c", "Non-linearity: LOWESS local slope")
    _titles(
        fig,
        "The slope of spectral count on LFQ is not identified",
        f"it depends on the error model (estimator) and the LFQ scale{DOT}protein "
        f"means over 8 runs, n = {n:,}{DOT}point estimates{DOT}log2\u2013log2",
        0.8,
    )
    handles: list[Any] = [
        Line2D(
            [], [], marker="o", ls="none", color=NEUTRAL, mec="black", mew=0.5, ms=8
        ),
        Line2D(
            [], [], marker="D", ls="none", color=NEUTRAL, mfc="white", mew=1.6, ms=8
        ),
        Line2D([], [], color=qcol["PSM"], lw=2.2),
        Line2D([], [], color=qcol["NSAF"], lw=2.2),
        Line2D([], [], color=REF_GRAY, lw=1.1, ls="--"),
    ]
    labels = [
        "protein LFQ (a, b)",
        f"LFQ per unique peptide: LFQ {MINUS} log2 n unique peptides (a, b)",
        "PSM",
        "NSAF",
        "slope 1",
    ]
    return fig, _legend_figure(handles, labels, ncol=1, width=6.4), numbers


# --------------------------------------------------------------------------- #
# 4. Agreement by PSM count / LFQ decile
# --------------------------------------------------------------------------- #
def plot_by_psm_count(data: AgreementData, qcol: dict[str, str]) -> Result:
    """LFQ spread at fixed PSM count; residual SD vs Poisson floor by LFQ decile."""
    cells = data.t("cells")
    common = cells[cells["in_common"]]
    bp = data.t("binned_by_psm").set_index("psm_bin").loc[list(PSM_BINS)]
    bl = data.t("binned_by_lfq_decile")
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 9.0))
    fig.subplots_adjust(
        left=0.08, right=0.98, bottom=0.08, top=0.88, hspace=0.42, wspace=0.24
    )
    pos = np.arange(len(PSM_BINS))
    xt = [f"{b.replace('>=', GE)}\n{int(bp.loc[b].n_cells):,}" for b in PSM_BINS]
    numbers: dict[str, Any] = {}

    ax = axes[0, 0]
    groups = [
        common.loc[common.psm_bin == b, "lfq"].to_numpy(dtype=np.float64)
        for b in PSM_BINS
    ]
    for b, gvals in zip(PSM_BINS, groups, strict=True):
        row = bp.loc[b]
        if gvals.size != int(row.n_cells):
            raise ValueError(f"bin {b}: {gvals.size} cells vs table {row.n_cells}")
        check(f"bin {b} median", float(np.median(gvals)), _f(row.median_lfq), 1e-4)
        check(f"bin {b} q10", float(np.quantile(gvals, 0.1)), _f(row.lfq_q10), 1e-4)
        check(f"bin {b} q90", float(np.quantile(gvals, 0.9)), _f(row.lfq_q90), 1e-4)
    parts = ax.violinplot(groups, positions=pos, widths=0.8, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor(qcol["LFQ"])
        body.set_edgecolor("none")
        body.set_alpha(0.3)
    ax.vlines(pos, bp.lfq_q10, bp.lfq_q90, color=qcol["LFQ"], lw=2.4)
    ax.plot(pos, bp.median_lfq, "o", color="white", mec=qcol["LFQ"], mew=1.8, ms=6)
    ax.set_xticks(pos, xt)
    ax.set_xlabel("PSMs in the run (bin; n cells)")
    ax.set_ylabel(LFQ_AXIS)
    _panel(ax, "a", "LFQ distribution at fixed PSM count")
    width = (bp.lfq_q90 - bp.lfq_q10).to_numpy(dtype=np.float64)
    numbers["lfq_q10_q90_width_by_bin"] = dict(
        zip(PSM_BINS, np.round(width, 2).tolist(), strict=True)
    )

    ax = axes[0, 1]
    ax.errorbar(
        pos,
        bp.lfq_resid_sd,
        yerr=[
            bp.lfq_resid_sd - bp.lfq_resid_sd_ci_low,
            bp.lfq_resid_sd_ci_high - bp.lfq_resid_sd,
        ],
        fmt="o-",
        color=qcol["LFQ"],
        capsize=3,
        lw=1.6,
        ms=6,
    )
    pooled_sd = data.summary["3_breakdown"]["pooled_resid_sd"]["lfq_given_log2_psm"]
    ax.axhline(pooled_sd, color=REF_GRAY, lw=1.1, ls="--")
    ax.text(
        3.5,
        pooled_sd + 0.02,
        f"pooled SD {pooled_sd:.2f}",
        ha="center",
        va="bottom",
        fontsize=9,
        color=NEUTRAL,
    )
    ax.set_xticks(pos, xt)
    ax.set_xlim(-0.5, len(PSM_BINS) - 0.5)
    ax.set_ylim(0, None)
    ax.set_xlabel("PSMs in the run (bin; n cells)")
    ax.set_ylabel("SD of LFQ residual from trend (log2)")
    _panel(ax, "b", "LFQ scatter around the LFQ-on-PSM trend")
    numbers["lfq_resid_sd_by_bin"] = dict(
        zip(
            PSM_BINS,
            np.round(bp.lfq_resid_sd.to_numpy(dtype=np.float64), 3).tolist(),
            strict=True,
        )
    )

    ax = axes[1, 0]
    ok = bp.spearman_lfq_psm_within_bin.notna().to_numpy()
    rs = bp.spearman_lfq_psm_within_bin.to_numpy(dtype=np.float64)
    lo = bp.spearman_ci_low.to_numpy(dtype=np.float64)
    hi = bp.spearman_ci_high.to_numpy(dtype=np.float64)
    ax.errorbar(
        pos[ok],
        rs[ok],
        yerr=[rs[ok] - lo[ok], hi[ok] - rs[ok]],
        fmt="o-",
        color=qcol["PSM"],
        capsize=3,
        lw=1.6,
        ms=6,
    )
    for i in np.flatnonzero(~ok):
        ax.text(
            pos[i],
            0.03,
            "n/a\n(PSM\nconstant)",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color=NEUTRAL,
        )
    ax.axhline(0, color=REF_GRAY, lw=0.8)
    ax.set_xticks(pos, xt)
    ax.set_xlim(-0.5, len(PSM_BINS) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel("PSMs in the run (bin; n cells)")
    ax.set_ylabel(f"Spearman {RHO} within bin")
    _panel(ax, "c", "LFQ vs PSM rank agreement inside each bin")
    numbers["within_bin_spearman"] = {
        b: (None if not k else round(float(v), 3))
        for b, k, v in zip(PSM_BINS, ok, rs, strict=True)
    }

    ax = axes[1, 1]
    dec = bl.lfq_decile.to_numpy()
    for q, col in (("PSM", "psm"), ("NSAF", "nsaf")):
        sd = bl[f"{col}_resid_sd"].to_numpy(dtype=np.float64)
        ax.errorbar(
            dec,
            sd,
            yerr=[
                sd - bl[f"{col}_resid_sd_ci_low"],
                bl[f"{col}_resid_sd_ci_high"] - sd,
            ],
            fmt="o-",
            color=qcol[q],
            capsize=3,
            lw=1.6,
            ms=6,
        )
    floor = bl.psm_poisson_floor_sd.to_numpy(dtype=np.float64)
    ax.plot(dec, floor, "D--", color=NEUTRAL, lw=1.3, ms=5)
    share = bl.psm_poisson_floor_var_share.to_numpy(dtype=np.float64)
    check(
        "floor share recomputed",
        float(
            np.max(
                np.abs(
                    floor**2 / bl.psm_resid_sd.to_numpy(dtype=np.float64) ** 2 - share
                )
            )
        ),
        0.0,
        1e-5,
    )
    for d_, f_, s_ in zip(dec, floor, share, strict=True):
        ax.text(
            d_,
            f_ - 0.06,
            f"{100 * s_:.0f}%",
            ha="center",
            va="top",
            fontsize=8,
            color=NEUTRAL,
        )
    overall = data.summary["3_breakdown"]["poisson_floor_share_of_pooled_psm_resid_var"]
    ax.set_xticks(dec)
    ax.set_ylim(0, None)
    ax.set_xlabel("LFQ decile of cells (1 = lowest; ~2,011 cells each)")
    ax.set_ylabel("SD of residual from LFQ trend (log2)")
    _panel(
        ax,
        "d",
        f"Spectral residual SD vs Poisson floor (overall {100 * overall:.0f}%"
        " of PSM var.)",
    )
    numbers["decile"] = {
        "psm_resid_sd": np.round(
            bl.psm_resid_sd.to_numpy(dtype=np.float64), 3
        ).tolist(),
        "nsaf_resid_sd": np.round(
            bl.nsaf_resid_sd.to_numpy(dtype=np.float64), 3
        ).tolist(),
        "poisson_floor_sd": np.round(floor, 3).tolist(),
        "floor_var_share": np.round(share, 3).tolist(),
        "overall_floor_share": overall,
    }
    _titles(
        fig,
        "Where LFQ and spectral counts disagree: count level and abundance",
        f"{len(common):,} cells quantified by both{DOT}95% bootstrap CIs"
        f" (proteins resampled){DOT}{PROCESSING}",
        0.88,
    )
    handles: list[Any] = [
        Patch(facecolor=qcol["LFQ"], alpha=0.3),
        Line2D(
            [],
            [],
            color=qcol["LFQ"],
            lw=2.4,
            marker="o",
            mfc="white",
            mec=qcol["LFQ"],
            mew=1.8,
        ),
        Line2D([], [], color=qcol["LFQ"], marker="o", lw=1.6),
        Line2D([], [], color=qcol["PSM"], marker="o", lw=1.6),
        Line2D([], [], color=qcol["NSAF"], marker="o", lw=1.6),
        Line2D([], [], color=NEUTRAL, marker="D", ls="--", lw=1.3, ms=5),
        Line2D([], [], color=REF_GRAY, ls="--", lw=1.1),
    ]
    labels = [
        "LFQ distribution (a)",
        "LFQ median, 10\u201390% range (a)",
        "LFQ residual SD (b)",
        f"PSM (c: {RHO}; d: residual SD)",
        "NSAF residual SD (d)",
        "Poisson floor SD for PSM (d; % = floor share of variance)",
        "pooled SD (b)",
    ]
    return fig, _legend_figure(handles, labels, ncol=1, width=5.2), numbers


# --------------------------------------------------------------------------- #
# 5. Residual drivers (partial-residual plots)
# --------------------------------------------------------------------------- #
PROPERTIES = (
    "log2_length_rel",
    "log2_n_unique_peptides",
    "mbr_intensity_frac_mean",
    "log2_shared_psm_inflation",
)


def plot_residual_drivers(data: AgreementData, qcol: dict[str, str]) -> Result:
    """Component-plus-residual plots: PSM residual vs shared-PSM inflation; NSAF
    residual vs relative length; fitted partial slope vs the theoretical slope."""
    props = data.t("protein_properties")
    sub = props[props.complete8].reset_index(drop=True)
    reg = data.t("residual_regression")
    n = len(sub)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4))
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.13, top=0.8, wspace=0.24)
    specs = (
        (
            "PSM",
            "resid_log2_psm",
            "log2_psm",
            "log2_shared_psm_inflation",
            1.0,
            f"shared-PSM inflation, {MINUS}log2(1 {MINUS} f)",
            "log2 PSM residual, adjusted (log2)",
        ),
        (
            "NSAF",
            "resid_log2_nsaf",
            "log2_nsaf",
            "log2_length_rel",
            -1.0,
            "log2 relative protein length",
            "log2 NSAF residual, adjusted (log2)",
        ),
    )
    numbers: dict[str, Any] = {"n_proteins": n}
    X = sub[list(PROPERTIES)].to_numpy(dtype=np.float64)
    design = np.column_stack([np.ones(n), X])
    for ax, (q, ycol, rname, term, theory, xlabel, ylabel), letter in zip(
        axes, specs, "ab", strict=True
    ):
        y = sub[ycol].to_numpy(dtype=np.float64)
        beta, *_ = np.linalg.lstsq(design, y, rcond=None)
        rows = reg[reg.residual == rname].set_index("term")
        for k, prop in enumerate(PROPERTIES):
            check(
                f"{rname} coef {prop}",
                float(beta[k + 1]),
                _f(rows.loc[prop].coef_per_unit),
                1e-5,
            )
        k = PROPERTIES.index(term)
        xk = X[:, k]
        others = [j for j in range(len(PROPERTIES)) if j != k]
        partial = y - (X[:, others] - X[:, others].mean(axis=0)) @ beta[1:][others]
        slope = float(beta[k + 1])
        a = float(partial.mean() - slope * xk.mean())
        rr = rows.loc[term]
        ax.scatter(xk, partial, s=9, color=qcol[q], alpha=0.45, lw=0, rasterized=True)
        xs = np.array([xk.min(), xk.max()])
        ax.plot(xs, a + theory * xs, color=REF_GRAY, lw=1.6, ls="--")
        ax.plot(xs, a + slope * xs, color=FIT_BLACK, lw=2.0)
        n_zero = int(np.sum(xk == 0))
        extra = f"\n{n_zero:,} proteins with f = 0 (x = 0)" if q == "PSM" else ""
        _stat_box(
            ax,
            f"slope {_fmt_ci(slope, _f(rr.ci_low_per_unit), _f(rr.ci_high_per_unit))}"
            f" vs theory {theory:+.0f}".replace("-", MINUS).replace("+1", "1")
            + f"\nΔR² (drop term) = {_f(rr.delta_r2_drop):.2f}{DOT}"
            f"n = {n:,}" + extra,
            loc="upper left" if q == "PSM" else "upper right",
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        name = "shared-peptide PSMs" if q == "PSM" else "protein length"
        _panel(ax, letter, f"{q} residual vs {name}")
        numbers[q] = {
            "term": term,
            "slope": round(slope, 4),
            "ci": [round(_f(rr.ci_low_per_unit), 4), round(_f(rr.ci_high_per_unit), 4)],
            "theory": theory,
            "delta_r2_drop": round(_f(rr.delta_r2_drop), 4),
            "model_r2": round(_f(rr.model_r2), 4),
            "n_x_zero": n_zero,
        }
    check(
        "psm model r2",
        numbers["PSM"]["model_r2"],
        data.summary["4_residuals"]["regression"]["log2_psm"]["r2"],
        1e-4,
    )
    _titles(
        fig,
        "What drives the spectral-count residual from the LFQ trend",
        f"protein means, quantified in all 8 runs (n = {n:,}){DOT}residual from "
        f"LOWESS on mean LFQ, other 3 terms partialled out{DOT}95% CI (HC3)",
        0.8,
    )
    handles: list[Any] = [
        Line2D([], [], marker="o", ls="none", color=qcol["PSM"], alpha=0.6),
        Line2D([], [], marker="o", ls="none", color=qcol["NSAF"], alpha=0.6),
        Line2D([], [], color=FIT_BLACK, lw=2.0),
        Line2D([], [], color=REF_GRAY, lw=1.6, ls="--"),
    ]
    labels = [
        "protein (PSM residual)",
        "protein (NSAF residual)",
        "fitted partial slope (4-term OLS)",
        "theoretical slope",
    ]
    return fig, _legend_figure(handles, labels, ncol=2, width=6.2), numbers


# --------------------------------------------------------------------------- #
# 6. Detection classes
# --------------------------------------------------------------------------- #
def _cell_classes(cells: pd.DataFrame) -> pd.Series:
    lfq_any = cells.groupby("protein_id")["has_lfq"].transform("any")
    cls = np.select(
        [
            cells.has_lfq & cells.has_psm,
            cells.has_psm & ~cells.has_lfq & ~lfq_any,
            cells.has_psm & ~cells.has_lfq & lfq_any,
            cells.has_lfq & ~cells.has_psm,
        ],
        list(DETECTION_CLASSES),
        default="",
    )
    if np.any(cls == ""):
        raise ValueError("cells.tsv row with neither LFQ nor PSM.")
    return pd.Series(cls, index=cells.index)


def plot_detection_classes(data: AgreementData, registry: Path) -> Result:
    """Per-run stacked detection classes + ECDF of PSM counts by class."""
    det = data.t("detection_per_sample").set_index("sample_id").loc[list(data.samples)]
    cells = data.t("cells").copy()
    cells["cls"] = _cell_classes(cells)
    colors = local_colors("detection_class", list(DETECTION_CLASSES), registry)
    col_of = {
        DETECTION_CLASSES[0]: "n_both",
        DETECTION_CLASSES[1]: "n_psm_only_never_lfq",
        DETECTION_CLASSES[2]: "n_psm_only_lfq_other_runs",
        DETECTION_CLASSES[3]: "n_lfq_only",
    }
    counts = cells.pivot_table(
        index="sample_id", columns="cls", values="psm", aggfunc="size", fill_value=0
    )
    for sid in data.samples:
        for cls, col in col_of.items():
            if int(counts.loc[sid, cls]) != int(det.loc[sid, col]):
                raise ValueError(f"{sid}/{cls}: cells {counts.loc[sid, cls]} != table")
    if not (det.n_lfq_only == det.n_lfq_only_mbr_only).all():
        raise ValueError("LFQ-only cells are not all MBR-only.")
    summ = data.summary["5_detection"]
    check("lfq-only total", float(det.n_lfq_only.sum()), summ["lfq_only_cells_total"])
    check(
        "psm-only min",
        float((det.n_psm_only).min()),
        summ["per_sample_psm_only_range"][0],
    )
    check(
        "psm-only max",
        float((det.n_psm_only).max()),
        summ["per_sample_psm_only_range"][1],
    )
    dist = data.t("detection_psm_distribution")
    tab = {
        "both": DETECTION_CLASSES[0],
        "psm_only_never_lfq": DETECTION_CLASSES[1],
        "psm_only_lfq_other_runs": DETECTION_CLASSES[2],
    }
    binned = cells[cells.has_psm].groupby(["sample_id", "cls", "psm_bin"]).size()
    for r in dist.itertuples():
        if int(binned.loc[(r.sample_id, tab[r[2]], r.psm_bin)]) != int(r.n_proteins):
            raise ValueError(f"PSM-bin distribution mismatch at {r}")

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(13.0, 5.6), gridspec_kw={"width_ratios": (1.25, 1.0)}
    )
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.2, top=0.8, wspace=0.22)
    bottom = np.zeros(len(data.samples))
    for cls, col in col_of.items():
        v = det[col].to_numpy(dtype=np.float64)
        ax_a.bar(
            range(len(v)),
            v,
            bottom=bottom,
            color=colors[cls],
            width=0.72,
            edgecolor="white",
            lw=0.6,
        )
        bottom += v
    _sample_axis(ax_a, data)
    ax_a.set_ylabel("proteins in the run")
    _panel(ax_a, "a", "Proteins per run by what quantified them")

    for cls in DETECTION_CLASSES[:3]:
        v = np.sort(cells.loc[cells.cls == cls, "psm"].to_numpy(dtype=np.float64))
        ax_b.step(
            np.concatenate([[1.0], v]),
            np.concatenate([[0.0], np.arange(1, v.size + 1) / v.size]),
            where="post",
            color=colors[cls],
            lw=2.0,
        )
    ax_b.set_xscale("log")
    ax_b.set_xlim(1, None)
    _plain_log(ax_b)
    ax_b.set_ylim(0, 1.01)
    ax_b.set_xlabel("PSMs for the protein in the run (log scale)")
    ax_b.set_ylabel("cumulative fraction of run\u00d7protein cells")
    auc_lo, auc_hi = (
        det.auc_psm_both_gt_psm_only.min(),
        det.auc_psm_both_gt_psm_only.max(),
    )
    _panel(ax_b, "b", "PSM counts by class, 8 runs pooled")
    med = {
        c: float(np.median(cells.loc[cells.cls == c, "psm"]))
        for c in DETECTION_CLASSES[:3]
    }
    _stat_box(
        ax_b,
        "median PSMs per run:\n"
        f"both {det.median_psm_both.min():.0f}\u2013{det.median_psm_both.max():.0f}"
        f"{DOT}PSM-only {det.median_psm_psm_only.min():.0f}\n"
        f"P(both > PSM-only) {auc_lo:.2f}\u2013{auc_hi:.2f}",
        loc="lower right",
    )
    _titles(
        fig,
        "Detection classes per run: LFQ vs spectral identification",
        f"{len(data.samples)} runs{DOT}{cells.protein_id.nunique():,} proteins with "
        f"LFQ or PSMs in any run{DOT}contaminants excluded",
        0.8,
    )
    numbers: dict[str, Any] = {
        "per_run": {
            sid: {c: int(det.loc[sid, col]) for c, col in col_of.items()}
            for sid in data.samples
        },
        "pooled_median_psm": med,
        "auc_range": [round(float(auc_lo), 3), round(float(auc_hi), 3)],
        "colors": colors,
    }
    handles = [Patch(color=colors[c]) for c in DETECTION_CLASSES]
    return (
        fig,
        _legend_figure(handles, list(DETECTION_CLASSES), ncol=2, width=6.6),
        numbers,
    )


# --------------------------------------------------------------------------- #
# 7. No-LFQ reasons
# --------------------------------------------------------------------------- #
def plot_no_lfq_reasons(data: AgreementData, registry: Path, top_n: int = 20) -> Result:
    """The PSM-but-never-LFQ proteins: counts by reason, PSM distribution, top N."""
    nl = data.t("no_lfq_proteins")
    summ = data.summary["5_detection"]
    labels = dict(REASONS)
    colors = local_colors("no_lfq_reason", [lab for _k, lab in REASONS], registry)
    if len(nl) != summ["n_proteins_psm_never_lfq"]:
        raise ValueError("no-LFQ protein count mismatch.")
    for key, _lab in REASONS:
        check(
            f"reason {key}", float((nl.reason == key).sum()), summ["reason_counts"][key]
        )
    cells = data.t("cells")
    lfq_any = cells.groupby("protein_id").has_lfq.any()
    tot = cells[cells.has_psm].groupby("protein_id").psm.sum()
    ref = tot[lfq_any.reindex(tot.index).to_numpy()].to_numpy(dtype=np.float64)
    check(
        "lfq-any median total psm",
        float(np.median(ref)),
        summ["lfq_any_total_psm_median"],
    )

    fig = plt.figure(figsize=(14.0, 6.6))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=(1.0, 1.05),
        height_ratios=(0.75, 1.1),
        left=0.14,
        right=0.985,
        bottom=0.1,
        top=0.84,
        wspace=0.32,
        hspace=0.55,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[:, 1])

    keys = [k for k, _l in REASONS]
    ncount = [int((nl.reason == k).sum()) for k in keys]
    psums = [int(nl.loc[nl.reason == k, "total_psm"].sum()) for k in keys]
    ax_a.barh(range(3), ncount, color=[colors[labels[k]] for k in keys], height=0.62)
    for i, (c, p) in enumerate(zip(ncount, psums, strict=True)):
        ax_a.text(
            c + 12, i, f"{c:,} proteins{DOT}{p:,} PSMs", va="center", fontsize=9.5
        )
    ax_a.set_yticks(range(3), [labels[k] for k in keys])
    ax_a.set_ylim(2.5, -0.5)
    ax_a.set_xlim(0, max(ncount) * 1.75)
    ax_a.set_xlabel("proteins with PSMs but no LFQ in any run")
    _panel(ax_a, "a", f"By reason (n = {len(nl):,})")

    rng = np.random.default_rng(0)  # jitter only; deterministic
    rows = [("LFQ in ≥ 1 run (reference)", ref, NEUTRAL)] + [
        (
            labels[k],
            nl.loc[nl.reason == k, "total_psm"].to_numpy(dtype=np.float64),
            colors[labels[k]],
        )
        for k in keys
    ]
    for i, (_lab, vals, color) in enumerate(rows):
        ax_b.scatter(
            vals,
            i + rng.uniform(-0.28, 0.28, vals.size),
            s=5,
            color=color,
            alpha=0.35,
            lw=0,
            rasterized=True,
        )
        m = float(np.median(vals))
        ax_b.plot([m, m], [i - 0.36, i + 0.36], color="black", lw=2.0)
    ax_b.set_xscale("log")
    _plain_log(ax_b)
    ax_b.set_yticks(
        range(len(rows)),
        [f"{lab}\nn = {v.size:,}{DOT}median {np.median(v):.0f}" for lab, v, _c in rows],
    )
    ax_b.set_ylim(len(rows) - 0.5, -0.5)
    ax_b.set_xlabel("total PSMs over 8 runs (log scale; bar = median)")
    _panel(ax_b, "b", "PSM totals per protein")

    top = nl.sort_values(["total_psm", "entry"], ascending=[False, True]).head(top_n)
    listed = [r["entry"] for r in summ["top_never_lfq"]][:top_n]
    if list(top.entry) != listed and sorted(top.entry) != sorted(listed):
        raise ValueError(f"top-{top_n} differs from summary: {list(top.entry)}")
    yy = np.arange(len(top))
    ax_c.barh(
        yy, top.total_psm, color=[colors[labels[r]] for r in top.reason], height=0.7
    )
    ax_c.set_yticks(yy, [_short(e) for e in top.entry])
    ax_c.set_ylim(len(top) - 0.5, -0.5)
    for y_, v in zip(yy, top.total_psm, strict=True):
        ax_c.text(v + 25, float(y_), f"{v:,}", va="center", fontsize=8.5)
    ax_c.set_xlim(0, float(top.total_psm.max()) * 1.15)
    ax_c.set_xlabel("total PSMs over 8 runs")
    _panel(ax_c, "c", f"Top {top_n} by PSMs")
    _titles(
        fig,
        "Proteins identified by PSMs but never quantified by FlashLFQ",
        f"{len(nl):,} of {data.summary['inputs']['n_analysis_proteins']:,} proteins"
        f"{DOT}contaminants excluded{DOT}8 runs",
        0.84,
    )
    numbers: dict[str, Any] = {
        "reason_counts": dict(zip(keys, ncount, strict=True)),
        "reason_psm_totals": dict(zip(keys, psums, strict=True)),
        "median_total_psm": {
            **{k: float(nl.loc[nl.reason == k, "total_psm"].median()) for k in keys},
            "lfq_any_reference": float(np.median(ref)),
        },
        "n_lfq_any_reference": int(ref.size),
        "top": [
            [_short(e), int(p), r]
            for e, p, r in zip(top.entry, top.total_psm, top.reason, strict=True)
        ],
        "colors": colors,
    }
    handles = [Patch(color=colors[labels[k]]) for k in keys]
    handles.append(Patch(color=NEUTRAL))
    leg_labels = [labels[k] for k in keys] + ["LFQ in ≥ 1 run (reference, b)"]
    return fig, _legend_figure(handles, leg_labels, ncol=2, width=7.2), numbers


# --------------------------------------------------------------------------- #
# 8. Within-protein tracking
# --------------------------------------------------------------------------- #
def _hist(ax: Axes, r: FloatArray, color: str) -> None:
    ax.hist(
        r,
        bins=np.linspace(-1, 1, 41).tolist(),
        color=color,
        alpha=0.85,
        ec="white",
        lw=0.4,
    )
    ax.axvline(float(np.mean(r)), color="black", lw=1.6, ls="--")
    ax.axvline(float(np.median(r)), color="black", lw=1.8)
    ax.axvline(0, color=REF_GRAY, lw=0.8)
    ax.set_xlim(-1, 1)


def plot_within_protein(data: AgreementData, qcol: dict[str, str]) -> Result:
    """Per-protein run-to-run Pearson r (LFQ vs log2 PSM): 8 runs, between pair means
    and within pairs (pair-respecting nulls only), and median r by count level.

    Tests that treat the 8 runs as exchangeable (per-protein 8! permutations, the
    shared 8-run relabeling, the independent per-protein shuffle) are NOT shown:
    runs come in 4 pairs, so exchangeability is violated and effective n = 4 pairs.
    """
    wp = data.t("within_protein_correlation")
    nulls = data.t("within_protein_null")
    by_bin = data.t("within_protein_by_psm_bin")
    six = data.summary["6_within_protein"]
    color = qcol["PSM"]
    fig = plt.figure(figsize=(14.5, 6.0))
    gs = fig.add_gridspec(
        2,
        3,
        width_ratios=(1.15, 1.0, 1.0),
        wspace=0.3,
        hspace=0.62,
        left=0.06,
        right=0.985,
        bottom=0.11,
        top=0.8,
    )
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b1 = fig.add_subplot(gs[0, 1])
    ax_b2 = fig.add_subplot(gs[1, 1], sharex=ax_b1)
    ax_c = fig.add_subplot(gs[:, 2])
    numbers: dict[str, Any] = {}
    specs = (
        (
            ax_a,
            "r_psm_raw_pearson",
            "pearson",
            "psm_raw_pearson",
            "a",
            "Across 8 runs",
            False,
        ),
        (
            ax_b1,
            "r_psm_raw_pair_means",
            "pearson_between_pair_means",
            "psm_raw_pearson_between_pair_means",
            "b",
            "Between the 4 pair means",
            True,
        ),
        (
            ax_b2,
            "r_psm_raw_within_pair",
            "pearson_within_pair",
            "psm_raw_pearson_within_pair",
            "",
            "Within pairs (4 differences)",
            True,
        ),
    )
    for ax, col, method, skey, letter, title, pair_null in specs:
        r = wp[col].dropna().to_numpy(dtype=np.float64)
        row = _row(nulls, comparison="lfq_vs_psm_raw", method=method)
        block = six[skey]
        if r.size != int(row.n_proteins) or r.size != block["n_proteins"]:
            raise ValueError(f"{col}: n {r.size} != {row.n_proteins}")
        check(f"{col} median", float(np.median(r)), block["median_r"])
        check(f"{col} mean", float(np.mean(r)), block["mean_r"])
        check(f"{col} frac>0", float(np.mean(r > 0)), block["frac_r_positive"])
        text = (
            f"median {np.median(r):.2f}{DOT}mean {np.mean(r):.2f}\n"
            f"{100 * np.mean(r > 0):.0f}% > 0{DOT}n = {r.size:,}"
        )
        entry: dict[str, Any] = {
            "n": int(r.size),
            "median": round(float(np.median(r)), 4),
            "mean": round(float(np.mean(r)), 4),
            "frac_pos": round(float(np.mean(r > 0)), 4),
        }
        if pair_null:
            lo, hi = _f(row.global_null_q025), _f(row.global_null_q975)
            ax.axvspan(lo, hi, fill=False, hatch="///", ec="#9a9a9a", lw=0)
            p2 = _f(row.global_perm_p_two_sided)
            n_perm = int(row.global_n_perm)
            what = "pair relabelings" if "means" in method else "swap patterns"
            text += f"\nmean r: p = {p2:.3f} ({n_perm} {what})"
            entry |= {
                "null_95_of_mean_r": [round(lo, 4), round(hi, 4)],
                "p_two_sided": round(p2, 4),
                "n_perm": n_perm,
                "min_attainable_p": round((1 if "means" in method else 2) / n_perm, 4),
            }
        _hist(ax, r, color)
        _stat_box(ax, text, loc="upper left")
        ax.set_ylabel("proteins")
        if letter:
            _panel(ax, letter, title)
        else:
            ax.set_title(title, fontsize=11.5, fontweight="normal", loc="left", pad=8)
        ax.set_ylim(0, ax.get_ylim()[1] * (1.75 if pair_null else 1.3))
        numbers[col] = entry
    ax_a.set_xlabel("per-protein Pearson r, LFQ vs log2 PSM")
    plt.setp(ax_b1.get_xticklabels(), visible=False)
    ax_b2.set_xlabel("per-protein Pearson r, LFQ vs log2 PSM")

    sub = by_bin[(by_bin.comparison == "lfq_vs_psm_raw") & (by_bin.method == "pearson")]
    sub = sub.set_index("mean_psm_bin").loc[list(PSM_BINS)]
    pos = np.arange(len(PSM_BINS))
    ax_c.plot(pos, sub.median_r, "o-", color=color, lw=1.8, ms=7, mec="black", mew=0.5)
    ax_c.axhline(0, color=REF_GRAY, lw=0.8)
    ax_c.set_xticks(
        pos,
        [
            f"{b.replace('>=', GE)}\n{int(n):,}"
            for b, n in zip(PSM_BINS, sub.n_proteins, strict=True)
        ],
    )
    ax_c.set_ylim(-0.3, 1.0)
    ax_c.set_xlabel("mean PSMs per run (bin; n proteins)")
    ax_c.set_ylabel("median per-protein Pearson r (8 runs)")
    _panel(ax_c, "c", "Tracking vs count level")
    numbers["median_r_by_mean_psm_bin"] = dict(
        zip(
            PSM_BINS,
            np.round(sub.median_r.to_numpy(dtype=np.float64), 3).tolist(),
            strict=True,
        )
    )
    _titles(
        fig,
        "Does a protein's spectral count follow its LFQ from run to run?",
        f"proteins quantified by both in all 8 runs{DOT}LFQ median-normalized log2 "
        f"vs raw log2 PSM{DOT}runs are 4 pairs: effective n = 4",
        0.8,
    )
    handles: list[Any] = [
        Patch(facecolor=color, alpha=0.85),
        Line2D([], [], color="black", lw=1.8),
        Line2D([], [], color="black", lw=1.6, ls="--"),
        Patch(fill=False, hatch="///", ec="#9a9a9a"),
        Line2D([], [], color=color, marker="o", lw=1.8, mec="black", mew=0.5),
    ]
    labels = [
        "proteins (histogram, a\u2013b)",
        "median r",
        "mean r",
        "95% null range of the mean r, pair-respecting (b)",
        "median r per bin (c)",
    ]
    return fig, _legend_figure(handles, labels, ncol=2, width=8.0), numbers


# --------------------------------------------------------------------------- #
# 9. Contaminant sharing (finding 0011)
# --------------------------------------------------------------------------- #
def _id_entry(protein_id: str) -> str:
    """Display entry for a protein-quants id (real ``sp|ACC|ENTRY`` or no-accession
    ``sp|ENTRY|`` contaminant form)."""
    parts = protein_id.split("|")
    if len(parts) != 3:
        raise ValueError(f"unexpected protein id {protein_id!r}")
    return parts[2] if parts[2] else parts[1]


def contaminant_sharing_table(
    protein_groups: Sequence[str],
    peptide_quantified: NDArray[np.bool_],
    contaminant_ids: frozenset[str],
    contaminant_lfq_runs: dict[str, int],
    no_lfq: pd.DataFrame,
) -> pd.DataFrame:
    """Real no-LFQ proteins that would gain >= 1 unique peptide if the contaminant-list
    copies were removed from the peptide protein groups.

    A peptide row becomes unique to real protein P when P is its only
    non-contaminant member and it has >= 1 contaminant member. ``exclusive`` = every
    peptide row of P becomes unique (P shares only with contaminant copies).
    ``carrier`` = P's contaminant partner with LFQ in the most runs (tie: most
    shared rows). One row per gaining protein, sorted by total PSMs.
    """
    if len(protein_groups) != peptide_quantified.size:
        raise ValueError("protein_groups / peptide_quantified length mismatch.")
    rows_of: dict[str, list[int]] = {}
    gain_of: dict[str, list[int]] = {}
    partners_of: dict[str, dict[str, int]] = {}
    for i, groups in enumerate(protein_groups):
        members = groups.split(";")
        real = [m for m in members if m not in contaminant_ids]
        conts = [m for m in members if m in contaminant_ids]
        for m in real:
            rows_of.setdefault(m, []).append(i)
            for c in conts:
                partners_of.setdefault(m, {})
                partners_of[m][c] = partners_of[m].get(c, 0) + 1
        if len(real) == 1 and conts:
            gain_of.setdefault(real[0], []).append(i)
    out = []
    for r in no_lfq.itertuples():
        gained = gain_of.get(str(r.protein_id), [])
        if not gained:
            continue
        partners = partners_of[str(r.protein_id)]
        ranked = sorted(partners, key=lambda c: (-partners[c], _id_entry(c)))
        carrier = max(ranked, key=lambda c: (contaminant_lfq_runs[c], partners[c]))
        n_rows = len(rows_of[str(r.protein_id)])
        out.append(
            {
                "protein_id": r.protein_id,
                "entry": r.entry,
                "reason": r.reason,
                "total_psm": int(r.total_psm),
                "n_runs_with_psm": int(r.n_runs_with_psm),
                "n_peptide_rows": n_rows,
                "n_gain_unique_rows": len(gained),
                "n_gain_unique_quantified": int(peptide_quantified[gained].sum()),
                "exclusive_contaminant_sharing": len(gained) == n_rows,
                "contaminant_partners": ";".join(
                    f"{_id_entry(c)}({partners[c]})" for c in ranked
                ),
                "partner_lfq_runs": ";".join(
                    f"{_id_entry(c)}:{contaminant_lfq_runs[c]}" for c in ranked
                ),
                "carrier_contaminant": _id_entry(carrier),
                "carrier_lfq_runs": contaminant_lfq_runs[carrier],
            }
        )
    table = pd.DataFrame(out)
    if table.empty:
        raise ValueError("no protein gains a unique peptide -- unexpected.")
    return table.sort_values(
        ["total_psm", "entry"], ascending=[False, True], ignore_index=True
    )


def plot_contaminant_sharing(
    data: AgreementData, qcol: dict[str, str], table: pd.DataFrame
) -> Result:
    """The real no-LFQ proteins that would regain unique peptides without the
    contaminant-list copies: PSMs (exclusive sharers solid) and the number of runs
    in which the contaminant copy carries the LFQ."""
    n_all = len(data.t("no_lfq_proteins"))
    n68 = data.summary["5_detection"]["n_shared_only_sharing_with_contaminant_entry"]
    n, n_excl = len(table), int(table.exclusive_contaminant_sharing.sum())
    psm_total = int(table.total_psm.sum())
    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(11.5, 7.0),
        sharey=True,
        gridspec_kw={"width_ratios": (1.35, 1.0)},
    )
    fig.subplots_adjust(left=0.1, right=0.86, bottom=0.09, top=0.84, wspace=0.08)
    yy = np.arange(n)
    excl = table.exclusive_contaminant_sharing.to_numpy(dtype=bool)
    psm = table.total_psm.to_numpy(dtype=np.float64)
    ax_a.barh(
        yy[excl], psm[excl], color=qcol["PSM"], height=0.68, ec=qcol["PSM"], lw=1.0
    )
    ax_a.barh(
        yy[~excl],
        psm[~excl],
        color="white",
        height=0.68,
        ec=qcol["PSM"],
        lw=1.0,
        hatch="////",
    )
    for y_, v, g, t in zip(
        yy, psm, table.n_gain_unique_rows, table.n_peptide_rows, strict=True
    ):
        ax_a.text(v + 25, y_, f"{int(v):,}{DOT}{g}/{t}", va="center", fontsize=8.5)
    ax_a.set_xlim(0, float(psm.max()) * 1.3)
    ax_a.set_yticks(yy, [_short(e) for e in table.entry])
    ax_a.set_ylim(n - 0.5, -0.5)
    ax_a.set_xlabel("PSMs of the real protein, 8 runs")
    _panel(ax_a, "a", "Real protein: PSMs, no LFQ")
    runs = table.carrier_lfq_runs.to_numpy(dtype=np.float64)
    ax_b.barh(yy, runs, color=qcol["LFQ"], height=0.68)
    for y_, v, c in zip(yy, runs, table.carrier_contaminant, strict=True):
        ax_b.text(v + 0.15, y_, f"{int(v)}{DOT}{_short(c)}", va="center", fontsize=8.5)
    ax_b.set_xlim(0, 8)
    ax_b.set_xticks(range(0, 9, 2))
    ax_b.set_xlabel("runs with LFQ on the contaminant copy (of 8)")
    _panel(ax_b, "b", "Contaminant-list copy: LFQ runs")
    _titles(
        fig,
        "Proteins that would regain unique peptides without contaminant copies",
        f"{n} of {n_all:,} PSM-but-no-LFQ proteins ({psm_total:,} PSMs){DOT}"
        f"{n_excl} share only with contaminant copies{DOT}"
        f"contaminant list = 33 refined entries",
        0.84,
    )
    numbers: dict[str, Any] = {
        "n_proteins": n,
        "n_exclusive": n_excl,
        "total_psm": psm_total,
        "total_psm_exclusive": int(psm[excl].sum()),
        "n_sharing_any_contaminant_entry_68": n68,
        "n_carrier_lfq_any_run": int(np.sum(runs > 0)),
        "n_carrier_lfq_all_8": int(np.sum(runs == 8)),
        "rows": [
            [_short(e), int(p), int(g), int(t), bool(x), _short(c), int(k)]
            for e, p, g, t, x, c, k in zip(
                table.entry,
                table.total_psm,
                table.n_gain_unique_rows,
                table.n_peptide_rows,
                table.exclusive_contaminant_sharing,
                table.carrier_contaminant,
                table.carrier_lfq_runs,
                strict=True,
            )
        ],
    }
    handles = [
        Patch(fc=qcol["PSM"], ec=qcol["PSM"]),
        Patch(fc="white", ec=qcol["PSM"], hatch="////"),
        Patch(color=qcol["LFQ"]),
    ]
    labels = [
        "PSMs; shares only with contaminant copies",
        "PSMs; also shares with real proteins",
        "runs with LFQ on the contaminant copy",
    ]
    return fig, _legend_figure(handles, labels, ncol=1, width=5.2), numbers
