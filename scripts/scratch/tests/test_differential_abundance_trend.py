"""Tests for the limma-trend option of the DE module and the trend sensitivity runner.

Planted truth: an 8-sample, 4-pair design whose per-feature residual variances are
drawn from a scaled-inverse-chi2 prior whose scale *falls with mean abundance* (a known
mean-variance trend, known d0). ``trend=True`` must recover the trend and d0; with no
trend planted it must reduce to (almost) the constant prior; ``trend=False`` must be the
unchanged constant-prior path. When ``Rscript`` + limma are available the fit is also
checked against R ``eBayes(trend=TRUE)`` directly.
"""

from __future__ import annotations

import shutil
import subprocess
import warnings
from pathlib import Path
from typing import cast

import de_raloxifene_vs_control as base_runner
import de_trend_sensitivity as trend_runner
import numpy as np
import pandas as pd
import pytest
from analysis import differential_abundance as da
from analysis.result_io import load_result, save_result
from loaders.data_loading import Dataset

N_FEATURES = 3000
D0_TRUE = 6.0
REFERENCE, TREATED = "control", "raloxifene-d0"


def _metadata() -> pd.DataFrame:
    pairs = ["P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4"]
    return pd.DataFrame(
        {
            "sample_id": [f"S{i}" for i in range(8)],
            "condition": [REFERENCE, TREATED] * 4,
            "candidate_pair": pairs,
            "batch": ["A"] * 6 + ["B"] * 2,
            "sample_role": ["experimental"] * 8,
        }
    ).set_index("sample_id", drop=False)


def true_prior_sd(mean: np.ndarray, *, trend: bool) -> np.ndarray:
    """Planted prior SD: 0.30 at mean 18 falling log-linearly to 0.075 at mean 30."""
    if not trend:
        return np.full_like(mean, 0.15)
    return np.asarray(0.30 * np.exp(np.log(0.25) * (mean - 18.0) / 12.0), dtype=float)


def planted(
    seed: int = 3, *, trend: bool = True, n: int = N_FEATURES, shift_idx: int = 0
) -> tuple[Dataset, np.ndarray]:
    """Paired dataset with a planted (or flat) mean-variance trend in the prior.

    Returns the dataset and the true per-feature prior SD. Feature means are uniform
    on [18, 30]; residual variance ~ s0(mean)^2 * d0 / chi2(d0); the pair effect is
    large (sd 1) so only the paired design sees the small within-pair noise.
    ``shift_idx`` features get a +1 log2 shift in the treated arm.
    """
    rng = np.random.default_rng(seed)
    md = _metadata()
    mean = rng.uniform(18.0, 30.0, size=n)
    s0 = true_prior_sd(mean, trend=trend)
    sigma2 = s0**2 * D0_TRUE / rng.chisquare(D0_TRUE, size=n)
    pair_idx = md["candidate_pair"].map({"P1": 0, "P2": 1, "P3": 2, "P4": 3})
    pair_effect = rng.normal(0.0, 1.0, size=(4, n))
    x = mean + pair_effect[pair_idx.to_numpy()] - pair_effect.mean(axis=0)
    x = x + rng.normal(size=(8, n)) * np.sqrt(sigma2)
    treated = (md["condition"] == TREATED).to_numpy()
    x[np.ix_(treated, list(range(shift_idx)))] += 1.0
    ds = Dataset(
        abundances=x,
        feature_names=np.array([f"F{i:04d}" for i in range(n)], dtype=object),
        feature_metadata=pd.DataFrame({"accession": [f"A{i}" for i in range(n)]}),
        metadata=md,
        scale="log2",
    )
    return ds, s0


def _fit(ds: Dataset, *, trend: bool) -> da.DifferentialAbundanceResult:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", da.ZeroResidualVarianceWarning)
        return da.differential_abundance(
            ds,
            "condition",
            covariates=["candidate_pair"],
            reference={"condition": REFERENCE},
            method="moderated",
            trend=trend,
        )


# --------------------------------------------------------------------------- #
# Planted trend recovered
# --------------------------------------------------------------------------- #
def test_trend_recovers_planted_mean_variance_trend_and_d0() -> None:
    ds, s0_true = planted(trend=True)
    res = _fit(ds, trend=True)
    assert isinstance(res, da.DifferentialAbundanceTrendResult)
    tab = res.contrast_table.set_index("feature").reindex(ds.feature_names)
    prior_sd = np.sqrt(tab["prior_variance"].to_numpy())
    rel_err = np.abs(prior_sd / s0_true - 1.0)
    assert np.median(rel_err) < 0.10
    assert np.max(rel_err) < 0.35
    # Four-fold range planted; the fitted prior spans most of it, decreasing.
    assert prior_sd.max() / prior_sd.min() > 3.0
    order = np.argsort(tab["mean_abundance"].to_numpy())
    assert prior_sd[order[:100]].mean() > prior_sd[order[-100:]].mean()
    assert res.prior_df is not None and 4.0 < res.prior_df < 9.0
    assert res.trend_spline_df == 4
    assert len(res.trend_knots) == 4
    assert res.prior_variance is None


def test_constant_prior_misses_the_trend_and_underestimates_d0() -> None:
    """Without the trend the spread of the trend is absorbed as prior variance."""
    ds, _ = planted(trend=True)
    nt = _fit(ds, trend=False)
    tr = _fit(ds, trend=True)
    assert nt.prior_df is not None and tr.prior_df is not None
    assert nt.prior_df < tr.prior_df


def test_flat_truth_gives_flat_trend_close_to_constant_prior() -> None:
    ds, _ = planted(trend=False, seed=11)
    tr = _fit(ds, trend=True)
    nt = _fit(ds, trend=False)
    tab = tr.contrast_table
    prior_sd = np.sqrt(tab["prior_variance"].to_numpy())
    assert prior_sd.max() / prior_sd.min() < 1.3
    assert tr.prior_df is not None and nt.prior_df is not None
    assert abs(tr.prior_df - nt.prior_df) / nt.prior_df < 0.1
    # p-values nearly identical under a flat truth.
    merged = tab.merge(nt.contrast_table, on="feature", suffixes=("_t", "_n"))
    assert np.corrcoef(np.log(merged["p_t"]), np.log(merged["p_n"]))[0, 1] > 0.99


def test_trend_detects_planted_shift_at_low_variance_abundance() -> None:
    """A shifted high-abundance (low-variance) feature ranks higher under trend."""
    ds, _ = planted(trend=True, seed=5, shift_idx=0)
    # Plant a 0.35 log2 shift in the 20 highest-abundance features.
    md = ds.metadata
    treated = (md["condition"] == TREATED).to_numpy()
    top = np.argsort(ds.abundances.mean(axis=0))[-20:]
    x = ds.abundances.copy()
    x[np.ix_(treated, top)] += 0.35
    ds2 = Dataset(x, ds.feature_names, ds.feature_metadata, md, "log2")
    tr = _fit(ds2, trend=True).contrast_table.set_index("feature")
    nt = _fit(ds2, trend=False).contrast_table.set_index("feature")
    names = ds.feature_names[top]
    assert tr.loc[names, "p"].median() < nt.loc[names, "p"].median()


# --------------------------------------------------------------------------- #
# trend=False unchanged
# --------------------------------------------------------------------------- #
def test_trend_false_is_the_unchanged_constant_prior_path() -> None:
    ds, _ = planted(trend=True, n=600)
    default = da.differential_abundance(
        ds,
        "condition",
        covariates=["candidate_pair"],
        reference={"condition": REFERENCE},
    )
    explicit = _fit(ds, trend=False)
    assert type(default) is da.DifferentialAbundanceResult
    assert type(explicit) is da.DifferentialAbundanceResult
    assert "prior_variance" not in default.table.columns
    pd.testing.assert_frame_equal(default.table, explicit.table)
    assert isinstance(default.prior_variance, float)
    # Constant prior, recomputed by hand from the module's no-trend estimator.
    s0_sq, d0, _ = da._fit_f_distribution_prior(
        np.asarray(
            (default.table.loc[default.table["is_contrast"], "sigma"] ** 2).to_numpy()
        ),
        3,
    )
    assert default.prior_variance == pytest.approx(s0_sq, rel=1e-12)
    assert default.prior_df == pytest.approx(d0, rel=1e-12)


def test_trend_requires_moderated() -> None:
    ds, _ = planted(n=200)
    with pytest.raises(ValueError, match="trend=True applies only"):
        da.differential_abundance(ds, "condition", method="ols", trend=True)


def test_trend_result_round_trips(tmp_path: Path) -> None:
    ds, _ = planted(n=400)
    res = _fit(ds, trend=True)
    save_result(res, tmp_path / "r")
    back = load_result(tmp_path / "r", da.DifferentialAbundanceTrendResult)
    assert isinstance(res, da.DifferentialAbundanceTrendResult)
    pd.testing.assert_frame_equal(back.table, res.table)
    assert back.trend is True
    assert back.trend_knots == pytest.approx(res.trend_knots)
    assert back.prior_df == res.prior_df


# --------------------------------------------------------------------------- #
# Internals: spline basis + infinite-d0 branch
# --------------------------------------------------------------------------- #
def test_natural_spline_basis_spans_lines_and_is_linear_outside() -> None:
    knots = (0.0, 1.0, 2.5, 4.0)
    x = np.linspace(0.0, 4.0, 50)
    b = da._natural_spline_basis(x, knots)
    assert b.shape == (50, 4)
    coef, *_ = np.linalg.lstsq(b, 2.0 + 3.0 * x, rcond=None)
    assert np.allclose(b @ coef, 2.0 + 3.0 * x)
    outside = np.array([4.5, 5.0, 5.5, -1.0, -0.5, 0.0])
    bo = da._natural_spline_basis(outside, knots)
    # Second differences vanish beyond each boundary (linear extrapolation).
    assert np.allclose(bo[:3, :][2] - 2 * bo[:3, :][1] + bo[:3, :][0], 0.0, atol=1e-9)
    assert np.allclose(bo[3:, :][2] - 2 * bo[3:, :][1] + bo[3:, :][0], 0.0, atol=1e-9)


def test_trend_prior_infinite_d0_when_variances_lie_on_the_trend() -> None:
    cov = np.linspace(18.0, 30.0, 400)
    sigma2 = (0.3 * np.exp(-0.1 * (cov - 18.0))) ** 2
    tp = da._fit_f_distribution_prior_trend(sigma2, 3, cov)
    assert np.isinf(tp.d0)
    ratio = tp.s0_sq / sigma2
    assert np.allclose(ratio, ratio[0], rtol=1e-3)


# --------------------------------------------------------------------------- #
# Runner pieces
# --------------------------------------------------------------------------- #
def test_relabel_diagnostic_ranks_planted_labelling_first_under_both_models() -> None:
    ds, _ = planted(n=800, seed=9, shift_idx=150)
    diag = trend_runner.pair_relabel_diagnostic(ds)
    assert diag["n_labellings"] == 8
    standing = cast(dict[str, dict[str, int]], diag["observed_standing"])
    for model in ("trend", "notrend"):
        assert standing[model]["observed_pi0_rank_lowest_first"] == 1
        assert standing[model]["n_labellings_with_hits_q<0.10_ge_observed"] == 1


def test_trend_table_joins_notrend_and_flags_crossings() -> None:
    ds, _ = planted(n=800, seed=9, shift_idx=40)
    spec = base_runner.QuantitySpec("x", ("accession",))
    tr = trend_runner.run_design(ds, ("candidate_pair",), trend=True)
    nt = trend_runner.run_design(ds, ("candidate_pair",), trend=False)
    nt_table = base_runner.contrast_table(nt, ds, spec)
    table = trend_runner.trend_table(tr, nt_table, ds, spec)
    assert len(table) == 800
    for col in ("prior_sd", "notrend_q", "hit_q05_trend", "hit_q10_notrend"):
        assert col in table.columns
    changes = trend_runner.hit_changes(table, spec)
    block = cast(dict[str, int], changes["q<0.05"])
    assert block["n_crossing"] == (
        block["n_gained_under_trend"] + block["n_lost_under_trend"]
    )


# --------------------------------------------------------------------------- #
# Direct R limma cross-check (skipped without Rscript + limma)
# --------------------------------------------------------------------------- #
def _has_limma() -> bool:
    if shutil.which("Rscript") is None:
        return False
    proc = subprocess.run(
        ["Rscript", "-e", "suppressPackageStartupMessages(library(limma))"],
        capture_output=True,
        check=False,
    )
    return proc.returncode == 0


@pytest.mark.skipif(not _has_limma(), reason="Rscript with limma not available")
def test_matches_r_limma_trend(tmp_path: Path) -> None:
    ds, _ = planted(n=500, seed=21, shift_idx=10)
    in_dir, out_dir = tmp_path / "in", tmp_path / "out"
    in_dir.mkdir()
    ds.metadata.loc[:, ["sample_id", "condition", "candidate_pair", "batch"]].to_csv(
        in_dir / "samples.tsv", sep="\t", index=False
    )
    frame = pd.DataFrame(ds.abundances.T, columns=list(ds.metadata["sample_id"]))
    frame.insert(0, "feature", ds.feature_names.astype(str))
    frame.to_csv(in_dir / "planted.tsv", sep="\t", index=False, float_format="%.17g")
    script = Path(trend_runner.__file__).resolve().parent / "limma_trend_check.R"
    subprocess.run(
        ["Rscript", str(script), str(in_dir), str(out_dir)],
        check=True,
        capture_output=True,
    )
    r = pd.read_csv(out_dir / "planted_paired.tsv", sep="\t")
    py = _fit(ds, trend=True)
    ct = py.contrast_table.set_index("feature").reindex(r["feature"])
    assert py.prior_df == pytest.approx(float(r["df_prior"].iloc[0]), rel=1e-9)
    assert np.allclose(ct["statistic"], r["t"], rtol=1e-9, atol=1e-10)
    assert np.allclose(ct["p"], r["p"], rtol=1e-8, atol=1e-300)
    assert np.allclose(ct["prior_variance"], r["s2_prior"], rtol=1e-9)
