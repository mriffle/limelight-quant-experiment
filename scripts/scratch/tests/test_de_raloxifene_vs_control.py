"""Tests for the raloxifene-vs-control DE runner and the DE-module behaviour it uses.

Planted truth: a synthetic 4-pair (8-sample) dataset with a strong pair effect and a
known shift in a few features -- the paired design must recover them, and the
unpaired designs must have less power. Plus the template behaviour the analysis relies
on (paired design == paired t-test, BH, the limma zero-variance floor, result_io
round-trip) and the runner's diagnostics.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import de_raloxifene_vs_control as runner
import numpy as np
import pandas as pd
import pytest
from analysis import differential_abundance as da
from analysis.result_io import load_result, save_result
from loaders.data_loading import Dataset
from scipy import stats

N_FEATURES = 600
UP = list(range(8))
DOWN = list(range(8, 15))
SHIFT = 1.0


def _metadata() -> pd.DataFrame:
    """Study layout: 4 pairs (one control + one raloxifene each); pair 4 = batch B."""
    pairs = ["P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4"]
    return pd.DataFrame(
        {
            "sample_id": [f"S{i}" for i in range(8)],
            "condition": [runner.REFERENCE, runner.TREATED] * 4,
            "candidate_pair": pairs,
            "batch": ["A"] * 6 + ["B"] * 2,
            "sample_role": ["experimental"] * 8,
        }
    ).set_index("sample_id", drop=False)


def _planted_dataset(seed: int = 7, n_broad: int = 0) -> Dataset:
    """Pair effect (sd 1) >> within-pair noise (scaled-inv-chi2, s0 0.15, d0 6).

    Features ``UP``/``DOWN`` are shifted by +/-``SHIFT`` in raloxifene; ``n_broad``
    further features get a small +0.2 shift (a broad, weak signal).
    """
    rng = np.random.default_rng(seed)
    md = _metadata()
    sigma2 = 0.15**2 * 6 / rng.chisquare(6, size=N_FEATURES)
    pair_idx = md["candidate_pair"].map({"P1": 0, "P2": 1, "P3": 2, "P4": 3})
    pair_effect = rng.normal(0.0, 1.0, size=(4, N_FEATURES))
    x = 22.0 + pair_effect[pair_idx.to_numpy()]
    x = x + rng.normal(size=(8, N_FEATURES)) * np.sqrt(sigma2)
    treated = (md["condition"] == runner.TREATED).to_numpy()
    x[np.ix_(treated, UP)] += SHIFT
    x[np.ix_(treated, DOWN)] -= SHIFT
    broad = list(range(len(UP) + len(DOWN), len(UP) + len(DOWN) + n_broad))
    x[np.ix_(treated, broad)] += 0.2
    names = np.array([f"f{j:03d}" for j in range(N_FEATURES)], dtype=object)
    fmeta = pd.DataFrame({"accession": names, "entry": names})
    return Dataset(
        abundances=x,
        feature_names=names,
        feature_metadata=fmeta,
        metadata=md,
        scale="log2",
    )


SPEC = runner.QuantitySpec("synthetic", ("accession", "entry"))


def _tables(ds: Dataset) -> dict[str, pd.DataFrame]:
    return {
        d: runner.contrast_table(runner.run_design(ds, cov), ds, SPEC)
        for d, cov in runner.DESIGNS.items()
    }


# --------------------------------------------------------------------------- #
# Planted truth
# --------------------------------------------------------------------------- #
def test_paired_design_recovers_planted_shifts_with_correct_sign() -> None:
    table = _tables(_planted_dataset())["paired"].set_index("feature")
    planted = [f"f{j:03d}" for j in UP + DOWN]
    hits = set(table.index[table["q"] < 0.05])
    assert len(hits & set(planted)) >= 14  # recall >= 14/15
    assert len(hits - set(planted)) <= 2  # few false positives among 585 nulls
    assert (table.loc[[f"f{j:03d}" for j in UP], "log2fc"] > 0.7).all()
    assert (table.loc[[f"f{j:03d}" for j in DOWN], "log2fc"] < -0.7).all()
    covered = (table["ci_low"] <= SHIFT) & (table["ci_high"] >= SHIFT)
    assert covered.loc[[f"f{j:03d}" for j in UP]].mean() >= 0.8


def test_unpaired_designs_have_less_power() -> None:
    tables = _tables(_planted_dataset())
    planted = [f"f{j:03d}" for j in UP + DOWN]

    def n_hits(t: pd.DataFrame) -> int:
        return int(t.set_index("feature").loc[planted, "q"].lt(0.05).sum())

    def median_p(t: pd.DataFrame) -> float:
        return float(t.set_index("feature").loc[planted, "p"].median())

    assert n_hits(tables["unadjusted"]) < n_hits(tables["paired"])
    assert n_hits(tables["batch"]) < n_hits(tables["paired"])
    assert median_p(tables["unadjusted"]) > 10 * median_p(tables["paired"])


def test_residual_df_and_identical_effects_across_balanced_designs() -> None:
    ds = _planted_dataset()
    results = {d: runner.run_design(ds, c) for d, c in runner.DESIGNS.items()}
    assert {d: r.residual_df for d, r in results.items()} == {
        "paired": 3,
        "batch": 5,
        "unadjusted": 6,
    }
    # Condition is balanced within every pair and batch, so the OLS effect is the
    # mean within-pair difference under all three designs (only the SE differs).
    effects = [
        r.contrast_table.sort_values("feature")["effect"].to_numpy()
        for r in results.values()
    ]
    np.testing.assert_allclose(effects[0], effects[1], atol=1e-10)
    np.testing.assert_allclose(effects[0], effects[2], atol=1e-10)


def test_relabel_diagnostic_ranks_planted_labelling_first() -> None:
    diag = runner.pair_relabel_diagnostic(_planted_dataset(n_broad=200))
    assert diag["n_labellings"] == 8
    labellings = diag["labellings"]
    assert isinstance(labellings, list)
    assert sum(bool(r["observed"]) for r in labellings) == 1
    assert diag["observed_pi0_rank_lowest_first"] == 1


# --------------------------------------------------------------------------- #
# Template behaviour relied on
# --------------------------------------------------------------------------- #
def test_pair_covariate_ols_equals_paired_t_test() -> None:
    ds = _planted_dataset()
    res = da.differential_abundance(
        ds,
        "condition",
        covariates=["candidate_pair"],
        reference={"condition": runner.REFERENCE},
        method="ols",
    )
    ct = res.contrast_table.set_index("feature").loc[ds.feature_names]
    treated = ds.metadata["condition"].to_numpy() == runner.TREATED
    expected = stats.ttest_rel(ds.abundances[treated], ds.abundances[~treated], axis=0)
    np.testing.assert_allclose(ct["statistic"], expected.statistic, rtol=1e-8)
    np.testing.assert_allclose(ct["p"], expected.pvalue, rtol=1e-8)


def test_bh_matches_statsmodels() -> None:
    from statsmodels.stats.multitest import multipletests

    p = np.random.default_rng(1).uniform(size=500) ** 2
    np.testing.assert_allclose(
        da.bh_adjust(p), multipletests(p, method="fdr_bh")[1], rtol=1e-12
    )


def test_prior_fit_recovers_d0_and_floor_inactive_without_zeros() -> None:
    rng = np.random.default_rng(3)
    s2 = 0.04 * 8 / rng.chisquare(8, size=20000) * rng.chisquare(3, 20000) / 3
    s0_sq, d0, n_floored = da._fit_f_distribution_prior(s2, 3)
    assert n_floored == 0
    assert 6.0 < d0 < 10.5
    assert 0.03 < s0_sq < 0.05


def test_floor_keeps_roundoff_zero_variances_from_collapsing_d0() -> None:
    rng = np.random.default_rng(4)
    s2 = 0.04 * 8 / rng.chisquare(8, size=2000) * rng.chisquare(3, 2000) / 3
    _, d0_clean, _ = da._fit_f_distribution_prior(s2, 3)
    with_zeros = np.concatenate([s2, np.full(20, 1e-31)])
    with pytest.warns(da.ZeroResidualVarianceWarning, match="20 feature"):
        _, d0_zero, n_floored = da._fit_f_distribution_prior(with_zeros, 3)
    assert n_floored == 20
    # Unfloored, log(1e-31) ~ -71 dominates the spread of log-variances: d0 -> ~0.3
    # (what the unmodified template gave on the real NSAF / PSM paired fits).
    # Floored (limma), 1% zeros still pull d0 down, but it stays a real prior.
    assert 1.5 < d0_zero < d0_clean
    floored, n = da._floor_variances(with_zeros)
    assert n == 20
    assert floored.min() == pytest.approx(1e-5 * np.median(with_zeros))


def test_floor_raises_when_most_variances_are_zero() -> None:
    s2 = np.concatenate([np.zeros(60), np.ones(40)])
    with pytest.raises(ValueError, match="More than half"):
        da._fit_f_distribution_prior(s2, 3)


def test_result_round_trips_through_result_io(tmp_path: Path) -> None:
    res = runner.run_design(_planted_dataset(), ("candidate_pair",))
    save_result(res, tmp_path / "r")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        back = load_result(tmp_path / "r", da.DifferentialAbundanceResult)
    pd.testing.assert_frame_equal(back.table, res.table)
    assert back.residual_df == 3
    assert back.prior_df == res.prior_df
    assert back.contrast_terms == res.contrast_terms


# --------------------------------------------------------------------------- #
# Runner helpers
# --------------------------------------------------------------------------- #
def test_calibration_read_labels() -> None:
    rng = np.random.default_rng(5)
    uniform = rng.uniform(size=5000)
    assert runner.calibration_read(uniform)["read"] == "approximately uniform"
    spiked = np.concatenate([uniform[:4000], rng.uniform(0, 0.01, size=1000)])
    assert "spike at 0" in str(runner.calibration_read(spiked)["read"])
    conservative = rng.beta(2.0, 1.0, size=5000)
    read = str(runner.calibration_read(conservative)["read"])
    assert "hump near 1" in read
    assert "deficit near 0" in read


def test_storey_pi0() -> None:
    rng = np.random.default_rng(6)
    assert runner.storey_pi0(rng.uniform(size=20000)) == pytest.approx(1.0, abs=0.03)
    mix = np.concatenate([rng.uniform(size=10000), np.full(10000, 1e-6)])
    assert runner.storey_pi0(mix) == pytest.approx(0.5, abs=0.03)


def test_drop_constant_features_and_design_validation() -> None:
    ds = _planted_dataset()
    ds.abundances[:, 20] = 5.0
    kept, dropped = runner.drop_constant_features(ds)
    assert dropped == ["f020"]
    assert kept.abundances.shape == (8, N_FEATURES - 1)
    assert len(kept.feature_metadata) == N_FEATURES - 1

    runner.validate_design(ds.metadata)
    broken = ds.metadata.copy()
    broken.loc["S1", "condition"] = runner.REFERENCE
    with pytest.raises(ValueError, match="Expected 4"):
        runner.validate_design(broken)
