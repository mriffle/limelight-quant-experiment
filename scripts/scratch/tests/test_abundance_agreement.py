"""Tests for the LFQ-vs-spectral-count agreement module and runner.

Focus: (1) the Limelight -> protein-quants join is a verified bijection and lands the
right PSM row on the right protein (synthetic + real-data spot check); (2) the NSAF
length identity recovers planted protein lengths through Limelight-style rounding;
plus invariants of the statistical helpers the analysis leans on (slopes, the
zero-truncated Poisson floor, exact permutation nulls, the within-pair identity,
the cluster bootstrap, peptide sharing).
"""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import abundance_agreement as runner
import numpy as np
import pandas as pd
import pytest
from analysis import abundance_agreement as aa

ROOT = Path(__file__).resolve().parents[3]


# --------------------------------------------------------------------------- #
# Join
# --------------------------------------------------------------------------- #
def test_align_by_first_member_is_a_bijection() -> None:
    prot = ["p_a", "p_b", "p_c", "p_d"]
    first = ["p_c", "p_a", "p_d", "p_b"]
    j = aa.align_by_first_member(prot, first)
    assert [first[k] for k in j] == prot


@pytest.mark.parametrize(
    ("prot", "first"),
    [
        (["a", "b"], ["a", "a"]),  # duplicate first member
        (["a", "a"], ["a", "b"]),  # duplicate protein id
        (["a", "b"], ["a", "c"]),  # unmatched
        (["a", "b"], ["a", "b", "c"]),  # length mismatch
    ],
)
def test_align_by_first_member_fails_loud(prot: list[str], first: list[str]) -> None:
    with pytest.raises(ValueError):
        aa.align_by_first_member(prot, first)


@pytest.mark.skipif(
    not (ROOT / "results/qc_states/manifest.json").is_file(),
    reason="real processing states not materialized",
)
def test_real_join_lands_psm_rows_on_the_right_protein() -> None:
    d, cont, info = runner.load_all(ROOT, ROOT / "results/qc_states")
    assert len(d.ids) + len(cont.ids) == 4344
    assert len(cont.ids) == 33
    # Independent of the loaders: parse the raw Limelight protein dump + label map.
    dump = pd.read_csv(
        ROOT / "data/protein-limelight-table-dump.txt",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    label_map = pd.read_csv(ROOT / "results/stage2/limelight_label_map.tsv", sep="\t")
    label_of = dict(zip(label_map.sample_id, label_map.label, strict=True))
    first_key = dump["Protein(s)"].str.split(",").str[0].str.strip()
    by_key = dict(zip(first_key, range(len(dump)), strict=True))
    raw = np.column_stack(
        [
            dump[f"PSMs ({label_of[s]})"]
            .str.replace(",", "", regex=False)
            .replace("", "0")
            .astype(float)
            .to_numpy()
            for s in d.samples
        ]
    )
    for col, pid in enumerate(d.ids):
        row = by_key[pid.split("_", 2)[2]]  # strip "psvid_<n>_"
        got = np.nan_to_num(d.psm[:, col], nan=0.0)
        np.testing.assert_array_equal(got, raw[row], err_msg=pid)
    # The recovered normalized LFQ equals the stored complete-case state.
    dev = info["lfq_matches_normalized_log_max_abs_dev"]
    assert isinstance(dev, float)
    assert dev < runner.LFQ_STATE_TOL


# --------------------------------------------------------------------------- #
# NSAF rounding + length identity
# --------------------------------------------------------------------------- #
def test_round_like_limelight_nsaf() -> None:
    x = np.array([2.3449e-5, 1.2e-7, 0.01849, 0.0014, np.nan])
    out = aa.round_like_limelight_nsaf(x)
    np.testing.assert_allclose(out[:4], [2.34e-5, 1.2e-7, 0.018, 0.001])
    assert np.isnan(out[4])


def _synthetic_nsaf(
    seed: int = 1, n_prot: int = 400, n_runs: int = 8
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    length = rng.integers(80, 2500, n_prot).astype(float)
    abundance = np.exp(rng.normal(0, 2.0, n_prot))
    abundance[:6] *= 400  # a few very abundant proteins -> NSAF >= 1e-3 only
    loading = np.exp(rng.normal(0, 0.2, n_runs))
    mu = loading[:, None] * abundance[None, :] * length[None, :] / 150.0
    psm = rng.poisson(mu).astype(float)
    psm[psm == 0] = np.nan
    saf = psm / length[None, :]
    nsaf = saf / np.nansum(saf, axis=1, keepdims=True)
    return length, psm, aa.round_like_limelight_nsaf(nsaf)


def test_length_identity_recovers_planted_lengths() -> None:
    length, psm, nsaf = _synthetic_nsaf()
    est = aa.estimate_relative_length(psm, nsaf)
    truth = np.log2(length)
    precise = est.method == "precise"
    interval = est.method == "interval"
    assert precise.sum() > 300
    assert interval.sum() >= 1, "planted abundant proteins must use the interval path"
    offset = np.median((est.log2_length - truth)[precise])
    err = est.log2_length - truth - offset
    # 3-significant-figure rounding -> <= ~0.01 log2 error on precise proteins.
    assert np.nanmax(np.abs(err[precise])) < 0.01
    # Interval proteins: the truth lies inside the intersected rounding interval.
    lo = est.interval_low[interval] - offset
    hi = est.interval_high[interval] - offset
    assert np.all(truth[interval] >= lo - 0.01)
    assert np.all(truth[interval] <= hi + 0.01)
    # Cross-run consistency on precise cells is at the rounding level.
    assert np.nanmax(est.max_abs_dev_precise) < 0.015


def test_length_identity_fails_loud_on_misaligned_missingness() -> None:
    _, psm, nsaf = _synthetic_nsaf(n_prot=50)
    nsaf = nsaf.copy()
    col = int(np.flatnonzero(~np.isnan(nsaf[0]))[0])
    nsaf[0, col] = np.nan
    with pytest.raises(ValueError, match="missingness"):
        aa.estimate_relative_length(psm, nsaf)


# --------------------------------------------------------------------------- #
# Statistical helpers
# --------------------------------------------------------------------------- #
def test_slopes_recover_planted_line_and_deming_limits() -> None:
    rng = np.random.default_rng(3)
    x = rng.normal(0, 2, 5000)
    y = 0.8 * x + 1.0
    est = aa.slope_estimates(x, y, lam=1.0)
    for value in est.values():
        assert value == pytest.approx(0.8, rel=1e-9)
    y_noisy = 0.8 * x + rng.normal(0, 0.5, len(x))
    est = aa.slope_estimates(x, y_noisy, lam=1e9)
    assert est["deming"] == pytest.approx(est["ols_y_on_x"], rel=1e-4)
    assert est["ols_y_on_x"] < est["sma"] < est["inverse_ols_x_on_y"]


def test_zero_truncated_poisson_floor() -> None:
    assert aa.zt_poisson_mu_from_mean(1.0) == 0.0
    assert aa.zt_poisson_log2_var(0.0) == 0.0
    mu = aa.zt_poisson_mu_from_mean(3.0)
    assert mu / -math.expm1(-mu) == pytest.approx(3.0)
    # Large-count delta method: var(log2 K) ~ 1 / (mu ln2^2).
    assert aa.zt_poisson_log2_var(400.0) == pytest.approx(
        1 / (400 * math.log(2) ** 2), rel=0.02
    )


def test_exact_perm_pvalue_and_global_null() -> None:
    x = np.arange(8, dtype=float)[None, :].repeat(3, axis=0)
    y = np.vstack([np.arange(8.0), np.arange(8.0)[::-1], [3, 1, 4, 1, 5, 9, 2, 6]])
    p = aa.exact_rowwise_perm_pvalues(x, y, "pearson")
    # |r| = 1 is attained only by the identity and the reversal: 2 / 8!.
    assert p[0] == pytest.approx(2 / math.factorial(8))
    assert p[1] == pytest.approx(2 / math.factorial(8))
    assert 0 < p[2] <= 1
    obs, null = aa.global_perm_mean_corr(x, y, "pearson")
    assert len(null) == math.factorial(8)
    assert null[0] == pytest.approx(obs)
    assert abs(float(null.mean())) < 1e-12  # centered rows -> permutation mean 0


def test_within_pair_identity_and_signflip_null() -> None:
    rng = np.random.default_rng(5)
    ctrl_x, trt_x = rng.normal(size=(20, 4)), rng.normal(size=(20, 4))
    ctrl_y, trt_y = rng.normal(size=(20, 4)), rng.normal(size=(20, 4))
    dx, dy = trt_x - ctrl_x, trt_y - ctrl_y
    r = aa.within_pair_corr(dx, dy)
    for i in range(20):
        xs = np.concatenate([ctrl_x[i], trt_x[i]])
        ys = np.concatenate([ctrl_y[i], trt_y[i]])
        pm_x = np.tile(0.5 * (ctrl_x[i] + trt_x[i]), 2)
        pm_y = np.tile(0.5 * (ctrl_y[i] + trt_y[i]), 2)
        expected = np.corrcoef(xs - pm_x, ys - pm_y)[0, 1]
        assert r[i] == pytest.approx(expected)
    obs, null = aa.global_signflip_mean_corr(dx, dy)
    assert len(null) == 16
    assert null[0] == pytest.approx(float(np.mean(r)))
    assert obs == pytest.approx(float(np.mean(r)))


def test_cluster_sampler_draws_whole_groups() -> None:
    groups = np.array(["a", "a", "b", "c", "c", "c"])
    draw = aa.cluster_index_sampler(groups)
    rng = np.random.default_rng(0)
    for _ in range(50):
        idx = draw(rng)
        picked = groups[idx]
        for g, n in pd.Series(picked).value_counts().items():
            assert n % int((groups == g).sum()) == 0


def test_peptide_summary_and_shared_psm() -> None:
    prot = ["P1", "P2", "P3"]
    groups = ["P1", "P1", "P1;P2", "P2", "P3;CONT"]
    inten = np.array([[10.0, 0, 5, 2, 7], [0, 4, 5, 0, 7]])
    det = np.array(
        [
            ["MSMS", "NotDetected", "MSMS", "MBR", "MSMS"],
            ["NotDetected", "MBR", "MSMS", "NotDetected", "MSMS"],
        ]
    )
    s = aa.peptide_protein_summary(
        prot, groups, inten, det, flag_ids=frozenset({"CONT"})
    )
    t = s.per_protein
    assert t.loc["P1", "n_unique_rows"] == 2
    assert t.loc["P1", "n_unique_quantified"] == 2
    assert t.loc["P3", "n_unique_rows"] == 0
    assert bool(t.loc["P3", "shares_with_flagged"])
    assert not bool(t.loc["P1", "shares_with_flagged"])
    np.testing.assert_array_equal(s.n_unique_msms[:, 0], [1, 0])
    np.testing.assert_array_equal(s.n_unique_mbr[:, 0], [0, 1])
    assert s.mbr_intensity_frac[0, 1] == pytest.approx(1.0)  # P2 run 0: only MBR
    assert np.isnan(s.mbr_intensity_frac[1, 1])

    shared = aa.shared_psm_fraction(
        [["A"], ["B", "B2"]],
        [["A"], ["A", "B"], ["B2"], ["B", "B2"]],
        np.array([[5.0, 3.0, 2.0, 7.0]]),
    )
    np.testing.assert_array_equal(shared.psm_total[0], [8.0, 12.0])
    np.testing.assert_array_equal(shared.psm_shared[0], [3.0, 3.0])
    assert shared.n_members_unmatched == 0


def test_psm_bin_labels() -> None:
    labels = aa.psm_bin_labels(np.array([1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 500.0]))
    assert list(labels) == [
        "1",
        "2-3",
        "2-3",
        "4-7",
        "4-7",
        "8-15",
        "8-15",
        "16-31",
        "16-31",
        ">=32",
        ">=32",
    ]


def test_all_permutations_identity_first() -> None:
    perms = aa.all_permutations(4)
    assert perms.shape == (24, 4)
    assert list(perms[0]) == [0, 1, 2, 3]
    assert {tuple(p) for p in perms} == set(itertools.permutations(range(4)))
