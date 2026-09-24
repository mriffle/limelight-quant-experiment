"""LFQ protein abundance vs spectral counts (NSAF, PSM counts): agreement toolkit.

Pure, study-agnostic functions used by ``scripts/scratch/abundance_agreement.py``. Not
seeded from a ``lib/`` template: nothing in ``lib/`` covers method-agreement analysis
(correlation with cluster-bootstrap CIs, errors-in-variables slopes, the zero-truncated
Poisson floor of log counts, the NSAF length identity, within-feature permutation
nulls). Written to ``conventions/coding.md`` standards; a ``lib/`` "method-agreement"
template may be worth contributing.

Conventions used throughout
---------------------------
* Matrices are ``(n_samples, n_features)`` like the project ``Dataset`` contract;
  ``NaN`` = not observed (no LFQ value / no PSMs).
* All agreement is computed on the log2 scale.
* Confidence intervals are percentile bootstrap intervals that resample **proteins**
  (clusters). When a statistic pools several runs of one protein, whole proteins are
  resampled so within-protein dependence is respected.

The NSAF length identity
------------------------
NSAF for protein ``i`` in run ``r`` is ``(PSM_ir / L_i) / S_r`` with
``S_r = sum_j PSM_jr / L_j``, so ``log2(PSM_ir / NSAF_ir) = log2 L_i + log2 S_r``. This
is an additive two-way layout: :func:`estimate_relative_length` recovers ``log2 L_i``
up to one global constant. Limelight exports NSAF < 1e-3 to three significant figures
(relative error <= 0.5 %) but NSAF >= 1e-3 to three **decimals** (up to 50 % error at
0.001), so only the former cells are used for the point fit. Proteins with no such cell
get an interval estimate: the intersection over runs of the per-run rounding intervals.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import optimize, stats
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.multitest import multipletests

__script_meta__: dict[str, object] = {
    "task": "abundance-agreement",
    "kind": "module",
    "provides": [
        "align_by_first_member",
        "run_normalization_factors",
        "round_like_limelight_nsaf",
        "estimate_relative_length",
        "correlation",
        "bootstrap_ci",
        "cluster_index_sampler",
        "LowessTrend",
        "fit_lowess",
        "slope_estimates",
        "dynamic_range",
        "psm_bin_labels",
        "zt_poisson_mu_from_mean",
        "zt_poisson_log2_var",
        "rowwise_corr",
        "all_permutations",
        "exact_rowwise_perm_pvalues",
        "global_perm_mean_corr",
        "independent_perm_null_mean_corr",
        "all_sign_vectors",
        "within_pair_corr",
        "global_signflip_mean_corr",
        "peptide_protein_summary",
        "shared_psm_fraction",
        "median_polish",
        "bh",
        "prob_superiority",
    ],
    "uses": [],
    "seeded_from": None,
    "description": (
        "Agreement toolkit for LFQ vs NSAF / PSM counts on log2: first-member join, "
        "run-factor recovery, NSAF length-identity estimate (precise-cell two-way "
        "median fit + rounding-interval intersection), correlations with cluster "
        "bootstrap CIs, LOWESS trends, OLS/inverse-OLS/SMA/orthogonal/Deming slopes, "
        "dynamic range, zero-truncated Poisson log2 variance, within-protein "
        "correlations with exact per-protein and global run-label permutation nulls, "
        "peptide-level summaries (unique peptides, MBR fraction, shared-PSM fraction)."
    ),
}

FloatArray = npt.NDArray[np.float64]
BoolArray = npt.NDArray[np.bool_]
IntArray = npt.NDArray[np.int64]
CorrMethod = Literal["pearson", "spearman"]

NSAF_PRECISE_BELOW = 1e-3
NSAF_DECIMALS = 3
NSAF_SIG_FIGS = 3
PSM_BIN_EDGES: tuple[int, ...] = (1, 2, 4, 8, 16, 32)
DETECTION_TYPES: tuple[str, ...] = (
    "MSMS",
    "MBR",
    "NotDetected",
    "MSMSIdentifiedButNotQuantified",
    "MSMSAmbiguousPeakfinding",
)


# --------------------------------------------------------------------------- #
# Joins and scale recovery
# --------------------------------------------------------------------------- #
def align_by_first_member(
    protein_ids: Sequence[str], first_member_ids: Sequence[str]
) -> IntArray:
    """Return ``j`` with ``first_member_ids[j[i]] == protein_ids[i]`` for every ``i``.

    Fails loud unless the relation is a bijection (no duplicates on either side, every
    protein id present exactly once among the first members, same length).
    """
    prot = pd.Index(list(protein_ids))
    first = pd.Index(list(first_member_ids))
    if prot.has_duplicates:
        raise ValueError("protein_ids contain duplicates.")
    if first.has_duplicates:
        raise ValueError("first_member_ids contain duplicates.")
    if len(prot) != len(first):
        raise ValueError(
            f"Length mismatch: {len(prot)} protein ids vs {len(first)} first members."
        )
    indexer = np.asarray(first.get_indexer(prot), dtype=np.int64)
    if (indexer < 0).any():
        missing = [str(p) for p, k in zip(prot, indexer, strict=True) if k < 0][:5]
        raise ValueError(f"{int((indexer < 0).sum())} protein ids unmatched: {missing}")
    return indexer


def run_normalization_factors(
    raw: FloatArray, normalized: FloatArray, rtol: float = 1e-9
) -> FloatArray:
    """Per-run multiplicative factor ``normalized / raw`` (both linear, same columns).

    Fails loud unless the factor is constant within each run to ``rtol`` (i.e. the
    normalization really was a per-run scaling) and finite/positive.
    """
    if raw.shape != normalized.shape:
        raise ValueError(f"Shape mismatch {raw.shape} vs {normalized.shape}.")
    ratio = normalized / raw
    if not np.isfinite(ratio).all() or (ratio <= 0).any():
        raise ValueError("Non-finite or non-positive normalization ratio.")
    factors = np.median(ratio, axis=1)
    spread = np.max(np.abs(ratio / factors[:, None] - 1.0))
    if spread > rtol:
        raise ValueError(f"Normalization is not a per-run scaling (spread {spread}).")
    return np.asarray(factors, dtype=np.float64)


# --------------------------------------------------------------------------- #
# NSAF export rounding and the length identity
# --------------------------------------------------------------------------- #
def round_like_limelight_nsaf(
    nsaf: FloatArray,
    precise_below: float = NSAF_PRECISE_BELOW,
    decimals: int = NSAF_DECIMALS,
    sig_figs: int = NSAF_SIG_FIGS,
) -> FloatArray:
    """Mimic the Limelight display rounding of NSAF (for synthetic tests / docs).

    Values below ``precise_below`` keep ``sig_figs`` significant figures; values at or
    above it are rounded to ``decimals`` decimal places. NaN passes through.
    """
    out = np.full_like(nsaf, np.nan, dtype=np.float64)
    ok = np.isfinite(nsaf) & (nsaf > 0)
    small = ok & (nsaf < precise_below)
    big = ok & ~small
    exponent = np.floor(np.log10(nsaf[small]))
    scale = 10.0 ** (sig_figs - 1 - exponent)
    out[small] = np.round(nsaf[small] * scale) / scale
    out[big] = np.round(nsaf[big], decimals)
    return out


@dataclass(frozen=True)
class LengthEstimate:
    """Relative protein length from the NSAF identity (log2, median protein = 0).

    Attributes
    ----------
    log2_length:
        ``(n_features,)`` relative log2 length; NaN where no PSMs at all.
    method:
        ``"precise"`` (>= 1 cell with NSAF < 1e-3), ``"interval"`` (only rounded cells;
        midpoint of the intersected rounding intervals), ``"interval_empty"`` (the
        intervals do not intersect -> midpoint median, flagged) or ``"none"``.
    n_precise_cells:
        Number of precise cells used per protein.
    max_abs_dev_precise:
        Max |cell estimate - protein estimate| over precise cells (cross-run
        consistency; ~0.005 expected from 3-significant-figure rounding).
    interval_low, interval_high:
        Intersected log2-length bounds for interval proteins (NaN otherwise).
    run_offset:
        ``(n_samples,)`` estimated ``log2 S_r`` (sum-to-zero).
    """

    log2_length: FloatArray
    method: npt.NDArray[np.str_]
    n_precise_cells: IntArray
    max_abs_dev_precise: FloatArray
    interval_low: FloatArray
    interval_high: FloatArray
    run_offset: FloatArray


def estimate_relative_length(
    psm: FloatArray,
    nsaf: FloatArray,
    *,
    precise_below: float = NSAF_PRECISE_BELOW,
    decimals: int = NSAF_DECIMALS,
    max_iter: int = 100,
    tol: float = 1e-12,
) -> LengthEstimate:
    """Estimate relative log2 protein length from ``log2(PSM / NSAF)`` (see module doc).

    ``psm`` and ``nsaf`` are ``(n_samples, n_features)``; NaN = no PSMs. The two-way
    additive fit ``y_ir = a_i + c_r`` uses alternating medians over precise cells
    only. Interval proteins use the fitted run offsets ``c_r`` and, per run, the
    rounding interval ``[max(x - h, precise_below), x + h]`` with ``h = 0.5 *
    10**-decimals`` for the true NSAF.
    """
    if psm.shape != nsaf.shape:
        raise ValueError("psm and nsaf shapes differ.")
    if not np.array_equal(np.isnan(psm), np.isnan(nsaf)):
        raise ValueError("psm and nsaf missingness patterns differ.")
    if np.nanmin(psm) < 1 or np.nanmin(nsaf) <= 0:
        raise ValueError("Observed PSM counts must be >= 1 and NSAF > 0.")
    y = np.log2(psm / nsaf)
    precise = np.isfinite(y) & (nsaf < precise_below)
    yp = np.where(precise, y, np.nan)
    has_precise = precise.any(axis=0)
    n_s, n_f = y.shape
    c = np.zeros(n_s)
    a = np.full(n_f, np.nan)
    for iteration in range(max_iter):
        a_new = np.full(n_f, np.nan)
        a_new[has_precise] = np.nanmedian(yp[:, has_precise] - c[:, None], axis=0)
        c_new = np.nanmedian(yp[:, has_precise] - a_new[None, has_precise], axis=1)
        c_new = c_new - c_new.mean()
        delta = np.nanmax(np.abs(a_new - np.where(np.isnan(a), 0.0, a)))
        a, c = a_new, c_new
        if delta < tol and iteration > 0:
            break
    dev = np.where(precise, y - a[None, :] - c[:, None], np.nan)
    max_dev = np.full(n_f, np.nan)
    max_dev[has_precise] = np.nanmax(np.abs(dev[:, has_precise]), axis=0)

    h = 0.5 * 10.0 ** (-decimals)
    method = np.full(n_f, "none", dtype=object)
    method[has_precise] = "precise"
    lo_all = np.full(n_f, np.nan)
    hi_all = np.full(n_f, np.nan)
    rounded = np.isfinite(y) & ~precise
    need = rounded.any(axis=0) & ~has_precise
    for i in np.flatnonzero(need):
        rows = np.flatnonzero(rounded[:, i])
        x = nsaf[rows, i]
        base = np.log2(psm[rows, i]) - c[rows]
        lows = base - np.log2(x + h)
        highs = base - np.log2(np.maximum(x - h, precise_below))
        lo, hi = float(lows.max()), float(highs.min())
        if lo <= hi:
            a[i] = 0.5 * (lo + hi)
            method[i] = "interval"
        else:
            a[i] = float(np.median(0.5 * (lows + highs)))
            method[i] = "interval_empty"
        lo_all[i], hi_all[i] = lo, hi
    centre = float(np.nanmedian(a))
    return LengthEstimate(
        log2_length=np.asarray(a - centre, dtype=np.float64),
        method=np.asarray(method, dtype=str),
        n_precise_cells=np.asarray(precise.sum(axis=0), dtype=np.int64),
        max_abs_dev_precise=np.asarray(max_dev, dtype=np.float64),
        interval_low=np.asarray(lo_all - centre, dtype=np.float64),
        interval_high=np.asarray(hi_all - centre, dtype=np.float64),
        run_offset=np.asarray(c, dtype=np.float64),
    )


# --------------------------------------------------------------------------- #
# Correlation + bootstrap
# --------------------------------------------------------------------------- #
def correlation(x: FloatArray, y: FloatArray, method: CorrMethod) -> float:
    """Pearson or Spearman correlation of two equal-length finite vectors."""
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("x and y must be equal-length 1-D arrays.")
    if len(x) < 3:
        return math.nan
    if method == "spearman":
        x = np.asarray(stats.rankdata(x), dtype=np.float64)
        y = np.asarray(stats.rankdata(y), dtype=np.float64)
    xc = x - x.mean()
    yc = y - y.mean()
    denom = math.sqrt(float(xc @ xc) * float(yc @ yc))
    return math.nan if denom == 0 else float(xc @ yc) / denom


def cluster_index_sampler(
    groups: npt.NDArray[np.generic],
) -> Callable[[np.random.Generator], IntArray]:
    """Return a sampler drawing a cluster-bootstrap row index (resampling groups).

    When every group is a singleton this reduces to the ordinary bootstrap.
    """
    codes, uniques = pd.factorize(pd.Series(groups), sort=False)
    codes = np.asarray(codes, dtype=np.int64)
    n_groups = len(uniques)
    if n_groups == len(codes):

        def draw_simple(rng: np.random.Generator) -> IntArray:
            return np.asarray(rng.integers(0, n_groups, n_groups), dtype=np.int64)

        return draw_simple
    order = np.argsort(codes, kind="stable")
    counts = np.bincount(codes, minlength=n_groups)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])

    def draw_cluster(rng: np.random.Generator) -> IntArray:
        picked = rng.integers(0, n_groups, n_groups)
        sizes = counts[picked]
        offsets = np.repeat(starts[picked] - np.cumsum(sizes) + sizes, sizes)
        within = np.arange(int(sizes.sum())) + offsets
        return np.asarray(order[within], dtype=np.int64)

    return draw_cluster


def bootstrap_ci(
    stat: Callable[[IntArray], float],
    sampler: Callable[[np.random.Generator], IntArray],
    n_boot: int,
    rng: np.random.Generator,
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Percentile bootstrap CI of ``stat(index)`` over ``n_boot`` resamples."""
    if n_boot < 1:
        return (math.nan, math.nan)
    draws = np.array([stat(sampler(rng)) for _ in range(n_boot)], dtype=np.float64)
    draws = draws[np.isfinite(draws)]
    if len(draws) == 0:
        return (math.nan, math.nan)
    lo, hi = np.quantile(draws, [alpha / 2, 1 - alpha / 2])
    return (float(lo), float(hi))


def bh(pvalues: FloatArray) -> FloatArray:
    """Benjamini-Hochberg adjusted p-values (NaN kept as NaN)."""
    out = np.full_like(pvalues, np.nan, dtype=np.float64)
    ok = np.isfinite(pvalues)
    if ok.any():
        out[ok] = multipletests(pvalues[ok], method="fdr_bh")[1]
    return out


def prob_superiority(a: FloatArray, b: FloatArray) -> float:
    """P(A > B) + 0.5 P(A = B) (the Mann-Whitney AUC); 0.5 = no shift."""
    if len(a) == 0 or len(b) == 0:
        return math.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(u) / (len(a) * len(b))


# --------------------------------------------------------------------------- #
# Trend + slopes + dynamic range
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class LowessTrend:
    """A LOWESS fit, evaluable anywhere in its x-range by linear interpolation."""

    x_sorted: FloatArray
    y_fit: FloatArray
    frac: float

    def __call__(self, x: FloatArray) -> FloatArray:
        return np.asarray(np.interp(x, self.x_sorted, self.y_fit), dtype=np.float64)

    def grid(self, n: int = 200) -> pd.DataFrame:
        """Fit on an even grid with the local slope (finite differences)."""
        gx = np.linspace(self.x_sorted[0], self.x_sorted[-1], n)
        gy = self(gx)
        return pd.DataFrame({"x": gx, "fitted": gy, "local_slope": np.gradient(gy, gx)})


def fit_lowess(
    x: FloatArray, y: FloatArray, frac: float, it: int = 3, delta_frac: float = 0.005
) -> LowessTrend:
    """statsmodels LOWESS (``frac``, ``it`` robustifying iterations) of y on x."""
    if len(x) < 10:
        raise ValueError("Too few points for LOWESS.")
    delta = delta_frac * float(np.ptp(x))
    fit = lowess(y, x, frac=frac, it=it, delta=delta, return_sorted=True)
    xs = np.asarray(fit[:, 0], dtype=np.float64)
    ys = np.asarray(fit[:, 1], dtype=np.float64)
    ux, first = np.unique(xs, return_index=True)
    return LowessTrend(x_sorted=ux, y_fit=ys[first], frac=frac)


def slope_estimates(x: FloatArray, y: FloatArray, lam: float) -> dict[str, float]:
    """Line slopes of y on x under different error models.

    * ``ols_y_on_x``: all scatter attributed to y (attenuated by noise in x).
    * ``inverse_ols_x_on_y``: all scatter attributed to x (1 / slope of x on y).
    * ``sma``: standardized / reduced major axis, ``sign(r) * sd_y / sd_x`` -- the
      ratio of spreads, i.e. the relative dynamic range of the two measures.
    * ``orthogonal``: Deming with error-variance ratio 1 (same log2 units).
    * ``deming``: Deming with ``lam = var(err_y) / var(err_x)`` supplied.
    """
    xc = x - x.mean()
    yc = y - y.mean()
    n = len(x)
    sxx = float(xc @ xc) / (n - 1)
    syy = float(yc @ yc) / (n - 1)
    sxy = float(xc @ yc) / (n - 1)

    def deming(d: float) -> float:
        return (syy - d * sxx + math.sqrt((syy - d * sxx) ** 2 + 4 * d * sxy**2)) / (
            2 * sxy
        )

    return {
        "ols_y_on_x": sxy / sxx,
        "inverse_ols_x_on_y": syy / sxy,
        "sma": math.copysign(math.sqrt(syy / sxx), sxy),
        "orthogonal": deming(1.0),
        "deming": deming(lam),
    }


def dynamic_range(v: FloatArray) -> dict[str, float]:
    """Spread summaries of a log2 vector (log2 units; ``*_log10`` = orders of mag.)."""
    q = np.quantile(v, [0.01, 0.05, 0.95, 0.99])
    out = {
        "n": float(len(v)),
        "min": float(v.min()),
        "max": float(v.max()),
        "range": float(np.ptp(v)),
        "span_q01_q99": float(q[3] - q[0]),
        "span_q05_q95": float(q[2] - q[1]),
        "sd": float(v.std(ddof=1)),
    }
    log10_2 = math.log10(2.0)
    for key in ("range", "span_q01_q99", "span_q05_q95"):
        out[f"{key}_log10"] = out[key] * log10_2
    return out


def psm_bin_labels(
    counts: FloatArray, edges: Sequence[int] = PSM_BIN_EDGES
) -> npt.NDArray[np.object_]:
    """Label each count with its PSM bin (``"1"``, ``"2-3"``, ..., ``">=32"``)."""
    labels = []
    for lo, hi in itertools.pairwise(edges):
        labels.append(str(lo) if hi - lo == 1 else f"{lo}-{hi - 1}")
    labels.append(f">={edges[-1]}")
    idx = np.searchsorted(np.asarray(edges), counts, side="right") - 1
    out = np.full(len(counts), None, dtype=object)
    ok = np.isfinite(counts) & (idx >= 0)
    out[ok] = np.asarray(labels, dtype=object)[idx[ok]]
    return out


# --------------------------------------------------------------------------- #
# Zero-truncated Poisson floor of log2 counts
# --------------------------------------------------------------------------- #
def zt_poisson_mu_from_mean(mean: float) -> float:
    """Poisson rate whose zero-truncated mean equals ``mean`` (>= 1)."""
    if mean < 1:
        raise ValueError("A zero-truncated mean is >= 1.")
    if mean - 1 < 1e-9:
        return 0.0
    return float(optimize.brentq(lambda mu: mu / -math.expm1(-mu) - mean, 1e-12, mean))


def zt_poisson_log2_var(mu: float) -> float:
    """Variance of ``log2 K`` for ``K ~ Poisson(mu)`` conditioned on ``K >= 1``."""
    if mu <= 0:
        return 0.0
    kmax = int(mu + 12 * math.sqrt(mu) + 30)
    k = np.arange(1, kmax + 1)
    p = stats.poisson.pmf(k, mu)
    p = p / p.sum()
    lk = np.log2(k)
    m = float(p @ lk)
    return max(float(p @ (lk - m) ** 2), 0.0)


# --------------------------------------------------------------------------- #
# Within-protein (across-run) correlation and permutation nulls
# --------------------------------------------------------------------------- #
def _row_standardize(m: FloatArray) -> FloatArray:
    c = m - m.mean(axis=1, keepdims=True)
    sd = np.sqrt((c**2).mean(axis=1, keepdims=True))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.asarray(np.where(sd > 0, c / sd, np.nan), dtype=np.float64)


def _prep(m: FloatArray, method: CorrMethod) -> FloatArray:
    if method == "spearman":
        m = np.asarray(stats.rankdata(m, axis=1), dtype=np.float64)
    return _row_standardize(m)


def rowwise_corr(x: FloatArray, y: FloatArray, method: CorrMethod) -> FloatArray:
    """Per-row correlation of ``(n_features, n_runs)`` matrices (NaN if constant)."""
    if x.shape != y.shape:
        raise ValueError("Shape mismatch.")
    zx, zy = _prep(x, method), _prep(y, method)
    return np.asarray((zx * zy).mean(axis=1), dtype=np.float64)


def all_permutations(n: int) -> IntArray:
    """All ``n!`` permutations of ``range(n)`` (identity first)."""
    return np.asarray(list(itertools.permutations(range(n))), dtype=np.int64)


def exact_rowwise_perm_pvalues(
    x: FloatArray, y: FloatArray, method: CorrMethod, chunk: int = 16
) -> FloatArray:
    """Exact two-sided permutation p per row: share of the ``n!`` run relabelings of
    ``y`` with ``|r| >= |r_obs|`` (identity included; ties handled exactly)."""
    zx, zy = _prep(x, method), _prep(y, method)
    perms = all_permutations(x.shape[1])
    n_runs = x.shape[1]
    out = np.full(x.shape[0], np.nan)
    ok = np.flatnonzero(np.isfinite(zx).all(axis=1) & np.isfinite(zy).all(axis=1))
    for start in range(0, len(ok), chunk):
        rows = ok[start : start + chunk]
        zyp = zy[rows][:, perms]  # (chunk, n_perm, n_runs)
        r = np.einsum("pk,pqk->pq", zx[rows], zyp) / n_runs
        obs = np.abs(r[:, 0])[:, None]
        out[rows] = (np.abs(r) >= obs - 1e-12).mean(axis=1)
    return out


def global_perm_mean_corr(
    x: FloatArray, y: FloatArray, method: CorrMethod
) -> tuple[float, FloatArray]:
    """Mean per-row correlation under every **shared** run relabeling of ``y``.

    The same permutation is applied to all rows, so between-protein dependence
    (shared run effects) is preserved. Returns the observed mean (identity) and the
    exact null distribution over all ``n!`` permutations (identity included).
    Rows constant in either matrix are dropped.
    """
    zx, zy = _prep(x, method), _prep(y, method)
    keep = np.isfinite(zx).all(axis=1) & np.isfinite(zy).all(axis=1)
    zx, zy = zx[keep], zy[keep]
    n_runs = x.shape[1]
    cross = zx.T @ zy / (n_runs * len(zx))  # (n_runs, n_runs)
    perms = all_permutations(n_runs)
    null = cross[np.arange(n_runs)[None, :], perms].sum(axis=1)
    return float(null[0]), np.asarray(null, dtype=np.float64)


def independent_perm_null_mean_corr(
    x: FloatArray,
    y: FloatArray,
    method: CorrMethod,
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[FloatArray, FloatArray]:
    """Null of the mean and median per-row correlation when ``y``'s run labels are
    shuffled **independently per row** (destroys shared run structure too)."""
    zx, zy = _prep(x, method), _prep(y, method)
    keep = np.isfinite(zx).all(axis=1) & np.isfinite(zy).all(axis=1)
    zx, zy = zx[keep], zy[keep]
    means = np.empty(n_perm)
    medians = np.empty(n_perm)
    for b in range(n_perm):
        perm = np.argsort(rng.random(zy.shape), axis=1)
        r = (zx * np.take_along_axis(zy, perm, axis=1)).mean(axis=1)
        means[b] = r.mean()
        medians[b] = np.median(r)
    return means, medians


def all_sign_vectors(n: int) -> FloatArray:
    """All ``2**n`` vectors of +/-1 (all +1 first)."""
    return np.asarray(list(itertools.product((1.0, -1.0), repeat=n)), dtype=np.float64)


def within_pair_corr(dx: FloatArray, dy: FloatArray) -> FloatArray:
    """Per-row correlation of pair-centered values from within-pair differences.

    For two runs per pair, pair-centering leaves ``+/- d/2``; the Pearson correlation
    over the pair-centered runs equals the **uncentered** correlation of the
    differences ``dx``, ``dy`` (``(n_features, n_pairs)``). NaN if a row is all zero.
    """
    norm = np.sqrt((dx**2).sum(axis=1) * (dy**2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.asarray(
            np.where(norm > 0, (dx * dy).sum(axis=1) / norm, np.nan), dtype=np.float64
        )


def global_signflip_mean_corr(
    dx: FloatArray, dy: FloatArray
) -> tuple[float, FloatArray]:
    """Mean :func:`within_pair_corr` under every shared swap pattern of the two runs
    within each pair (a sign flip of ``dy``'s pair differences, identical for all
    rows). Exact null over ``2**n_pairs`` patterns, identity first."""
    norm = np.sqrt((dx**2).sum(axis=1) * (dy**2).sum(axis=1))
    keep = norm > 0
    prod = (dx * dy)[keep] / norm[keep, None]
    null = (prod @ all_sign_vectors(dx.shape[1]).T).mean(axis=0)
    return float(null[0]), np.asarray(null, dtype=np.float64)


# --------------------------------------------------------------------------- #
# Peptide-level summaries
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PeptideSummary:
    """Per-protein and per-(run, protein) peptide descriptors (protein order given).

    ``per_protein`` columns: n_unique_rows, n_unique_quantified, n_shared_rows,
    n_shared_quantified, det_<type> (counts over unique-peptide cells, 8 runs), and
    sharing_partners (up to 3 most frequent co-members on shared rows),
    shares_with_flagged.
    Matrices are ``(n_samples, n_proteins)``: unique peptides MSMS / MBR per run and
    the MBR share of the summed unique-peptide intensity (NaN when none quantified).
    """

    per_protein: pd.DataFrame
    n_unique_msms: FloatArray
    n_unique_mbr: FloatArray
    mbr_intensity_frac: FloatArray


def peptide_protein_summary(
    protein_ids: Sequence[str],
    peptide_protein_groups: Sequence[str],
    intensities: FloatArray,
    detection_type: npt.NDArray[np.str_],
    entry_of: dict[str, str] | None = None,
    flag_ids: frozenset[str] = frozenset(),
) -> PeptideSummary:
    """Summarize quantified peptides per protein.

    ``peptide_protein_groups`` holds each peptide's ``;``-joined protein-group ids;
    a peptide is *unique* when it lists exactly one id. ``intensities`` /
    ``detection_type`` are ``(n_samples, n_peptides)``, NaN / non-positive =
    unquantified. Protein ids absent from every peptide get zero counts.
    ``shares_with_flagged`` marks proteins with >= 1 shared peptide whose group also
    lists an id in ``flag_ids`` (e.g. contaminant-FASTA entries).
    """
    n_s, n_p = intensities.shape
    if detection_type.shape != intensities.shape or len(peptide_protein_groups) != n_p:
        raise ValueError("Peptide inputs are misaligned.")
    prot_index = pd.Index(list(protein_ids))
    members = [str(g).split(";") for g in peptide_protein_groups]
    is_unique = np.array([len(m) == 1 for m in members])
    quant = np.nan_to_num(intensities, nan=0.0) > 0
    quant_any = quant.any(axis=0)

    uniq_owner = np.array(
        [prot_index.get_loc(m[0]) if m[0] in prot_index else -1 for m in members]
    )
    uniq_owner = np.where(is_unique, uniq_owner, -1)
    n_prot = len(prot_index)
    sel = uniq_owner >= 0

    def per_owner(values: FloatArray) -> FloatArray:
        out = np.zeros((n_s, n_prot))
        for r in range(n_s):
            out[r] = np.bincount(
                uniq_owner[sel], weights=values[r, sel], minlength=n_prot
            )
        return out

    msms = (detection_type == "MSMS").astype(float)
    mbr = (detection_type == "MBR").astype(float)
    inten = np.nan_to_num(intensities, nan=0.0)
    n_msms = per_owner(msms)
    n_mbr = per_owner(mbr)
    tot_int = per_owner(inten)
    mbr_int = per_owner(inten * mbr)
    with np.errstate(invalid="ignore", divide="ignore"):
        mbr_frac = np.where(tot_int > 0, mbr_int / tot_int, np.nan)

    table = pd.DataFrame(index=prot_index)
    table["n_unique_rows"] = np.bincount(uniq_owner[sel], minlength=n_prot)
    table["n_unique_quantified"] = np.bincount(
        uniq_owner[sel], weights=quant_any[sel].astype(float), minlength=n_prot
    ).astype(int)
    for det in DETECTION_TYPES:
        counts = (detection_type == det).sum(axis=0).astype(float)
        table[f"det_{det}"] = np.bincount(
            uniq_owner[sel], weights=counts[sel], minlength=n_prot
        ).astype(int)

    shared_rows = np.zeros(n_prot, dtype=int)
    shared_quant = np.zeros(n_prot, dtype=int)
    flagged = np.zeros(n_prot, dtype=bool)
    partners: dict[int, dict[str, int]] = {}
    for k in np.flatnonzero(~is_unique):
        ids = [prot_index.get_loc(m) for m in members[k] if m in prot_index]
        has_flag = any(m in flag_ids for m in members[k])
        for i in ids:
            flagged[i] |= has_flag
            shared_rows[i] += 1
            shared_quant[i] += int(quant_any[k])
            bucket = partners.setdefault(i, {})
            for other in members[k]:
                if other != prot_index[i]:
                    label = entry_of.get(other, other) if entry_of else other
                    bucket[label] = bucket.get(label, 0) + 1
    table["n_shared_rows"] = shared_rows
    table["n_shared_quantified"] = shared_quant
    table["shares_with_flagged"] = flagged
    table["sharing_partners"] = [
        ";".join(
            f"{name}({cnt})"
            for name, cnt in sorted(
                partners.get(i, {}).items(), key=lambda kv: (-kv[1], kv[0])
            )[:3]
        )
        for i in range(n_prot)
    ]
    return PeptideSummary(
        per_protein=table,
        n_unique_msms=n_msms,
        n_unique_mbr=n_mbr,
        mbr_intensity_frac=np.asarray(mbr_frac, dtype=np.float64),
    )


@dataclass(frozen=True)
class SharedPsmSummary:
    """Per-(run, protein-row) PSMs on peptides unique to / shared beyond the row."""

    psm_total: FloatArray
    psm_shared: FloatArray
    n_members_unmatched: int


def shared_psm_fraction(
    row_members: Sequence[Sequence[str]],
    peptide_proteins: Sequence[Sequence[str]],
    peptide_psms: FloatArray,
) -> SharedPsmSummary:
    """Split each protein row's PSMs into peptides unique to the row vs shared.

    ``row_members[i]`` lists the accession keys of protein row ``i`` (a Limelight
    group can hold several indistinguishable members); ``peptide_proteins[k]`` the
    keys a peptide maps to; ``peptide_psms`` is ``(n_samples, n_peptides)`` with 0 for
    none. A peptide is *shared* when its keys map to more than one row.
    """
    key_to_row: dict[str, int] = {}
    for i, keys in enumerate(row_members):
        for key in keys:
            if key in key_to_row:
                raise ValueError(f"Accession key {key!r} belongs to two rows.")
            key_to_row[key] = i
    n_s = peptide_psms.shape[0]
    n_rows = len(row_members)
    total = np.zeros((n_s, n_rows))
    shared = np.zeros((n_s, n_rows))
    unmatched = 0
    for k, keys in enumerate(peptide_proteins):
        rows = set()
        for key in keys:
            if key in key_to_row:
                rows.add(key_to_row[key])
            else:
                unmatched += 1
        is_shared = len(rows) > 1 or any(key not in key_to_row for key in keys)
        for i in rows:
            total[:, i] += peptide_psms[:, k]
            if is_shared:
                shared[:, i] += peptide_psms[:, k]
    return SharedPsmSummary(
        psm_total=total, psm_shared=shared, n_members_unmatched=unmatched
    )


def median_polish(
    y: FloatArray, n_iter: int = 20
) -> tuple[float, FloatArray, FloatArray]:
    """Tukey median polish of ``(n_rows, n_cols)`` (NaN-aware): overall, row, col."""
    resid = y.copy()
    overall = 0.0
    row = np.zeros(y.shape[0])
    col = np.zeros(y.shape[1])
    for _ in range(n_iter):
        rm = np.nan_to_num(np.nanmedian(resid, axis=1))
        resid -= rm[:, None]
        row += rm
        m_col = float(np.median(col))
        col -= m_col
        overall += m_col
        cm = np.nan_to_num(np.nanmedian(resid, axis=0))
        resid -= cm[None, :]
        col += cm
        m_row = float(np.median(row))
        row -= m_row
        overall += m_row
    return overall, row, col
