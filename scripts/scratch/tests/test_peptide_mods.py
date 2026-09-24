"""Tests for ``analysis.peptide_mods`` and the peptide-mod runner's helpers.

Focus (as requested): mass decomposition, the dump <-> quant join, discordance
counting; plus the statistics helpers the analysis relies on, with hand-verified
fixtures (limma 3.58.1 reference values for the moderated one-sample t) and
property/planted-truth checks.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import peptide_mod_analysis as runner
import pytest
from analysis import peptide_mods as pm
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy import stats
from statsmodels.stats.multitest import multipletests

# --------------------------------------------------------------------------- #
# Masses and decomposition
# --------------------------------------------------------------------------- #


def test_raloxifene_adduct_delta_arithmetic() -> None:
    """C28H27NO4S monoisotopic 473.16608; minus 2 H = 471.15043 (hand-computed)."""
    assert pytest.approx(473.166079, abs=1e-6) == pm.RALOXIFENE_MONOISOTOPIC
    assert pytest.approx(471.150429, abs=1e-6) == pm.RALOXIFENE_ADDUCT_DELTA
    # the quant file's appended mass agrees within display precision
    assert abs(471.1504 - pm.RALOXIFENE_ADDUCT_DELTA) < 1e-4


@pytest.mark.parametrize(
    ("mass", "expected"),
    [
        (0.0, (0, 0, 0)),
        (57.021464, (1, 0, 0)),
        (114.042928, (2, 0, 0)),
        (15.9949, (0, 1, 0)),
        (73.016364, (1, 1, 0)),
        (130.037828, (2, 1, 0)),
        (31.9898, (0, 2, 0)),
        (47.9847, (0, 3, 0)),
        (89.011264, (1, 2, 0)),
        (471.1504, (0, 0, 1)),
        (528.171864, (1, 0, 1)),
        (585.193328, (2, 0, 1)),
        (487.14529999999996, (0, 1, 1)),
        (544.166764, (1, 1, 1)),
        (942.3008, (0, 0, 2)),
    ],
)
def test_decompose_every_observed_mass(
    mass: float, expected: tuple[int, int, int]
) -> None:
    """Every total mass present in peptide-quants.tsv has exactly one composition."""
    sols = pm.decompose_mass(mass)
    assert len(sols) == 1
    comp, err = sols[0]
    assert (comp.n_cam, comp.n_ox, comp.n_ralox) == expected
    assert abs(err) < 1e-4


def test_decompose_unexplained_and_invalid() -> None:
    assert pm.decompose_mass(42.010565) == []  # acetyl: not a searched delta
    with pytest.raises(ValueError, match=">= 0"):
        pm.decompose_mass(-5.0)
    with pytest.raises(ValueError, match="finite"):
        pm.decompose_mass(float("nan"))
    with pytest.raises(ValueError, match="tolerance"):
        pm.decompose_mass(57.0, tolerance=0.0)


@settings(max_examples=150, deadline=None)
@given(n_cam=st.integers(0, 6), n_ox=st.integers(0, 4), n_ralox=st.integers(0, 3))
def test_decompose_roundtrip_unique(n_cam: int, n_ox: int, n_ralox: int) -> None:
    """Property: a composition's own mass decomposes back to it, uniquely."""
    comp = pm.ModComposition(n_cam, n_ox, n_ralox)
    sols = pm.decompose_mass(round(comp.mass, 6))
    assert [s[0] for s in sols] == [comp]


def test_class_precedence() -> None:
    assert pm.ModComposition(0, 0, 0).mod_class == pm.UNMODIFIED
    assert pm.ModComposition(2, 0, 0).mod_class == pm.CAM_ONLY
    assert pm.ModComposition(1, 1, 0).mod_class == pm.OXIDIZED
    assert pm.ModComposition(1, 1, 1).mod_class == pm.ADDUCT
    assert pm.ModComposition(2, 1, 0).label == "2CAM+1Ox"
    assert pm.ModComposition(0, 0, 0).label == "none"
    with pytest.raises(ValueError, match="Negative"):
        pm.ModComposition(-1, 0, 0)


def test_decompose_feature_masses_flags() -> None:
    out = pm.decompose_feature_masses(
        [114.042928, 42.010565, 528.171864], ["PEPCK", "PEPTIDE", "ACDCK"]
    )
    assert out["mod_class"].tolist() == [pm.CAM_ONLY, pm.UNEXPLAINED, pm.ADDUCT]
    assert bool(out["cam_exceeds_cys"].iloc[0])  # 2 CAM, 1 Cys
    assert not bool(out["cam_exceeds_cys"].iloc[2])
    assert pd.isna(out["n_cam"].iloc[1])
    assert out["composition"].iloc[1] == "unexplained"
    with pytest.raises(ValueError, match="length"):
        pm.decompose_feature_masses([0.0], ["A", "B"])


# --------------------------------------------------------------------------- #
# Localized sequences, dump reading, join
# --------------------------------------------------------------------------- #


def test_parse_localized_sequence() -> None:
    p = pm.parse_localized_sequence("LPEC[471.15]EADDGC[57.02]PKM[15.99]R")
    assert p.base_sequence == "LPECEADDGCPKMR"
    assert [s.label for s in p.sites] == ["C4:Ralox", "C10:CAM", "M13:Ox"]
    assert p.composition == pm.ModComposition(1, 1, 1)
    assert [s.label for s in p.sites_of("Ralox")] == ["C4:Ralox"]
    assert pm.parse_localized_sequence("PEPTIDEK").sites == ()


@pytest.mark.parametrize("bad", ["", "PEPc", "PEP[99.99]K", "PE[57.02", "[57.02]PEP"])
def test_parse_localized_sequence_fails_loud(bad: str) -> None:
    with pytest.raises(ValueError):
        pm.parse_localized_sequence(bad)


def test_parse_psm_count() -> None:
    assert pm.parse_psm_count("") == 0
    assert pm.parse_psm_count("7") == 7
    assert pm.parse_psm_count("1,234") == 1234
    for bad in ("1.5", "12,34", "x"):
        with pytest.raises(ValueError):
            pm.parse_psm_count(bad)


def _write_dump(path: Path, mbr_cell: str = "(MBR)") -> None:
    header = ["Peptide Sequence", "Unique", "Pre", "Post", "Protein(s)"]
    for label in ("L1", "L2"):
        header += [
            f"PSMs ({label})",
            f"Quant ({label})",
            f"MBR ({label})",
            f"Shared group [column content string: (SHARED GROUP)] ({label})",
        ]
    rows = [
        [
            "AM[15.99]MK",
            "*",
            "K",
            "A",
            "sp|P1|A_HUMAN",
            "2",
            "1e5",
            "",
            "",
            "",
            "",
            mbr_cell,
            "",
        ],
        [
            "AMM[15.99]K",
            "*",
            "K",
            "A",
            "sp|P1|A_HUMAN",
            "1,001",
            "1e5",
            "",
            "",
            "3",
            "1e5",
            "",
            "",
        ],
        [
            "GC[471.15]K",
            "",
            "K",
            "A",
            "sp|P2|B_HUMAN",
            "",
            "",
            "",
            "",
            "4",
            "1e5",
            "",
            "",
        ],
        ["WWK", "*", "K", "A", "sp|P3|C_HUMAN", "1", "1e5", "", "", "", "", "", ""],
    ]
    text = "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n"
    path.write_text(text, encoding="utf-8")


def test_read_dump_join_and_aggregate(tmp_path: Path) -> None:
    dump_file = tmp_path / "dump.txt"
    _write_dump(dump_file)
    dump = pm.read_peptide_dump(dump_file, {"L1": "S_b", "L2": "S_a"}, ["S_a", "S_b"])
    # sample order follows sample_order, not file order: L2 -> S_a is column 0
    assert dump.psms.tolist() == [[0, 2], [3, 1001], [4, 0], [0, 1]]
    assert dump.mbr.tolist() == [
        [True, False],
        [False, False],
        [False, False],
        [False, False],
    ]
    assert dump.rows["key"].tolist() == [
        "AMMK|0|1|0",
        "AMMK|0|1|0",
        "GCK|0|0|1",
        "WWK|0|0|0",
    ]

    feature_keys = ["AMMK|0|1|0", "GCK|0|0|1", "NOTINDUMP|0|0|0"]
    join = pm.join_dump_to_features(dump.rows["key"].tolist(), feature_keys)
    assert join.dump_to_feature.tolist() == [0, 0, 1, -1]
    assert join.n_isomers.tolist() == [2, 1, 0]
    assert join.unmatched_dump_keys == ("WWK|0|0|0",)
    assert join.unmatched_feature_keys == ("NOTINDUMP|0|0|0",)

    psms, mbr, text = pm.aggregate_dump_to_features(join, dump, 3)
    assert psms[0].tolist() == [3.0, 1003.0]  # summed over the two Ox isoforms
    assert psms[1].tolist() == [4.0, 0.0]
    assert np.isnan(psms[2]).all()  # no dump row
    assert mbr[0].tolist() == [True, False]  # OR over isoforms
    assert text["dump_ox_sites"].iloc[0] == "M2:Ox | M3:Ox"
    assert text["dump_ralox_residues"].iloc[1] == "C"
    assert bool(text["dump_unique"].iloc[0])
    assert not bool(text["dump_unique"].iloc[1])


def test_read_dump_fails_loud(tmp_path: Path) -> None:
    dump_file = tmp_path / "dump.txt"
    _write_dump(dump_file, mbr_cell="yes")
    with pytest.raises(ValueError, match="MBR"):
        pm.read_peptide_dump(dump_file, {"L1": "S_b", "L2": "S_a"}, ["S_a", "S_b"])
    _write_dump(dump_file)
    with pytest.raises(ValueError, match="without a sample mapping"):
        pm.read_peptide_dump(dump_file, {"L1": "S_b"}, ["S_a", "S_b"])
    with pytest.raises(ValueError, match="1:1"):
        pm.read_peptide_dump(dump_file, {"L1": "S_a", "L2": "S_a"}, ["S_a", "S_b"])


def test_join_duplicate_feature_key_fails() -> None:
    with pytest.raises(ValueError, match="Duplicate"):
        pm.join_dump_to_features(["A|0|0|0"], ["A|0|0|0", "A|0|0|0"])


def test_parse_protein_members() -> None:
    out = pm.parse_protein_members(
        "psvid_1_sp|P08684|CP3A4_HUMAN;psvid_2_sp|CATD_HUMAN|"
    )
    assert out == [
        ("psvid_1_sp|P08684|CP3A4_HUMAN", "P08684", "CP3A4_HUMAN"),
        ("psvid_2_sp|CATD_HUMAN|", "", "CATD_HUMAN"),
    ]
    with pytest.raises(ValueError):
        pm.parse_protein_members("P08684")


# --------------------------------------------------------------------------- #
# Pairs and discordance
# --------------------------------------------------------------------------- #


def _meta() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "condition": [
                "control",
                "control",
                "treated",
                "treated",
                "control",
                "treated",
            ],
            "pair": ["P1", "P2", "P1", "P2", "P3", "P3"],
        }
    )


def test_pair_layout_and_counts_hand_fixture() -> None:
    layout = pm.pair_layout(
        _meta(),
        pair_col="pair",
        condition_col="condition",
        control="control",
        treated="treated",
    )
    assert layout.pair_ids == ("P1", "P2", "P3")
    assert layout.control_idx.tolist() == [0, 1, 4]
    assert layout.treated_idx.tolist() == [2, 3, 5]
    # rows = samples (C1, C2, T1, T2, C3, T3); columns = 3 features
    det = np.array(
        [
            [False, True, True],  # C1
            [False, False, True],  # C2
            [True, True, False],  # T1
            [True, False, False],  # T2
            [False, True, False],  # C3
            [True, False, False],  # T3
        ]
    )
    counts = pm.discordance_counts(det, layout)
    assert counts["n_treated_only"].tolist() == [3, 0, 0]
    assert counts["n_control_only"].tolist() == [0, 1, 2]
    assert counts["n_both"].tolist() == [0, 1, 0]
    assert counts["n_neither"].tolist() == [0, 1, 1]


def test_pair_layout_rejects_unbalanced_pair() -> None:
    meta = _meta()
    meta.loc[5, "condition"] = "control"
    with pytest.raises(ValueError, match="exactly one"):
        pm.pair_layout(
            meta,
            pair_col="pair",
            condition_col="condition",
            control="control",
            treated="treated",
        )


def test_pair_states_rejects_non_bool() -> None:
    layout = pm.pair_layout(
        _meta(),
        pair_col="pair",
        condition_col="condition",
        control="control",
        treated="treated",
    )
    with pytest.raises(ValueError, match="boolean"):
        pm.pair_states(np.zeros((6, 2)), layout)


@settings(max_examples=60, deadline=None)
@given(st.integers(0, 2**31 - 1))
def test_discordance_invariants(seed: int) -> None:
    """Property: states partition the pairs; swapping conditions swaps directions."""
    rng = np.random.default_rng(seed)
    det = rng.random((6, 25)) < 0.5
    meta = _meta()
    layout = pm.pair_layout(
        meta,
        pair_col="pair",
        condition_col="condition",
        control="control",
        treated="treated",
    )
    counts = pm.discordance_counts(det, layout)
    assert (counts.sum(axis=1) == 3).all()
    swapped = meta.assign(
        condition=meta["condition"].map({"control": "treated", "treated": "control"})
    )
    layout_s = pm.pair_layout(
        swapped,
        pair_col="pair",
        condition_col="condition",
        control="control",
        treated="treated",
    )
    counts_s = pm.discordance_counts(det, layout_s)
    assert counts_s["n_treated_only"].tolist() == counts["n_control_only"].tolist()
    assert counts_s["n_both"].tolist() == counts["n_both"].tolist()


def test_detection_masks_definitions() -> None:
    det = np.array([["MSMS", "MBR", "NotDetected", "MSMSIdentifiedButNotQuantified"]])
    psms = np.array([[2.0], [0.0], [np.nan], [1.0]])  # (n_features, n_samples)
    masks = runner.detection_masks(det, psms)
    assert masks["quantified"].tolist() == [[True, True, False, False]]
    assert masks["msms"].tolist() == [[True, False, False, False]]
    assert masks["psm"].tolist() == [[True, False, False, True]]


def test_reference_index() -> None:
    feats = pd.DataFrame(
        {
            "base_sequence": ["GCK", "GCK", "GCK", "AMK", "AMK", "YCK", "YCK"],
            "mod_class": [
                pm.CAM_ONLY,
                pm.ADDUCT,
                pm.UNMODIFIED,
                pm.UNMODIFIED,
                pm.OXIDIZED,
                pm.ADDUCT,
                pm.CAM_ONLY,
            ],
            "n_cam": [1, 0, 0, 0, 0, 1, 1],
            "n_ox": [0, 0, 0, 0, 1, 0, 0],
            "n_ralox": [0, 1, 0, 0, 0, 1, 0],
            "dump_ralox_residues": ["", "C", "", "", "", "Y", ""],
        }
    )
    feats["key"] = [
        pm.composition_key(b, pm.ModComposition(int(c), int(o), int(r)))
        for b, c, o, r in zip(
            feats["base_sequence"],
            feats["n_cam"],
            feats["n_ox"],
            feats["n_ralox"],
            strict=True,
        )
    ]
    ref = runner._reference_index(feats)
    # CAM-only GCK -> unmodified GCK; Cys adduct -> the CAM form (Cys would be CAM'd);
    # Ox -> unmodified AMK; Tyr adduct with 1 CAM -> 1-CAM form; unmodified -> -1
    assert ref.tolist() == [2, 0, -1, -1, 3, 6, -1]


# --------------------------------------------------------------------------- #
# Statistics helpers
# --------------------------------------------------------------------------- #


def test_bh_matches_statsmodels_and_passes_nan() -> None:
    rng = np.random.default_rng(3)
    p = rng.random(50)
    ref = multipletests(p, method="fdr_bh")[1]
    np.testing.assert_allclose(pm.bh_adjust(p), ref, rtol=1e-12)
    q = pm.bh_adjust(np.array([0.01, np.nan, 0.04]))
    assert np.isnan(q[1])
    np.testing.assert_allclose(q[[0, 2]], [0.02, 0.04])
    with pytest.raises(ValueError):
        pm.bh_adjust(np.array([1.2]))


def test_binomial_proportion_hand_values() -> None:
    res = pm.binomial_proportion(9, 10)
    assert res["prop"] == pytest.approx(0.9)
    assert res["p"] == pytest.approx(22 / 1024)  # 2 * (1 + 10) / 2^10
    assert res["ci_low"] == pytest.approx(0.554983, abs=1e-5)  # Clopper-Pearson
    assert res["ci_high"] == pytest.approx(0.997471, abs=1e-5)
    assert np.isnan(pm.binomial_proportion(0, 0)["prop"])


def test_fisher_and_mantel_haenszel() -> None:
    tab = np.array([[12, 3], [40, 45]])
    f = pm.fisher_odds_ratio(tab)
    assert f["p"] == pytest.approx(stats.fisher_exact(tab).pvalue)
    assert f["ci_low"] < f["odds_ratio"] < f["ci_high"]
    mh = pm.mantel_haenszel(np.stack([tab, tab]))
    assert mh["odds_ratio"] == pytest.approx((12 * 45) / (3 * 40))
    assert np.isnan(pm.fisher_odds_ratio(np.array([[0, 0], [1, 2]]))["p"])


def test_pair_level_log_odds_ratios() -> None:
    tabs = np.array([[[4, 2], [10, 10]], [[3, 0], [5, 5]], [[0, 0], [1, 1]]])
    lor = pm.pair_level_log_odds_ratios(tabs)
    assert lor[0] == pytest.approx(np.log(2.0))
    assert lor[1] == pytest.approx(np.log((3.5 * 5.5) / (0.5 * 5.5)))  # Haldane
    assert np.isnan(lor[2])


def test_welch_matches_scipy() -> None:
    rng = np.random.default_rng(5)
    x, y = rng.normal(0.3, 1, 40), rng.normal(0, 2, 60)
    w = pm.welch_difference(x, y)
    assert w["p"] == pytest.approx(stats.ttest_ind(x, y, equal_var=False).pvalue)
    assert w["ci_low"] < w["diff"] < w["ci_high"]


def test_mann_whitney_shift_recovers_planted_shift() -> None:
    rng = np.random.default_rng(9)
    y = rng.normal(0, 1, 400)
    x = rng.normal(0, 1, 300) + 0.5
    res = pm.mann_whitney_shift(x, y)
    assert res["ci_low"] < 0.5 < res["ci_high"]
    assert res["rank_biserial"] > 0
    assert res["p"] < 1e-6


def test_paired_t_and_sign_test() -> None:
    res = pm.paired_t(np.array([1.0, 2.0, 3.0, 4.0]))
    assert res["mean"] == pytest.approx(2.5)
    assert res["p"] == pytest.approx(stats.ttest_1samp([1, 2, 3, 4], 0).pvalue)
    assert pm.sign_test(4, 0) == pytest.approx(0.125)
    with pytest.raises(ValueError):
        pm.paired_t(np.array([1.0, np.nan]))


def test_moderated_one_sample_matches_limma() -> None:
    """Reference: R limma 3.58.1 lmFit(d, intercept) + eBayes on the same matrix
    (written with 10 decimals): s2.prior 0.0276085593983933, df.prior
    3.29630091045506; first rows of topTable(confint=TRUE)."""
    rng = np.random.default_rng(11)
    s2 = 0.2**2 * 5 / rng.chisquare(5, size=40)
    d = rng.normal(size=(4, 40)) * np.sqrt(s2)
    d[:, :3] += 0.8
    d = np.round(d, 10)
    res = pm.moderated_one_sample(d)
    assert res.s0_sq == pytest.approx(0.0276085593983933, rel=1e-6)
    assert res.d0 == pytest.approx(3.29630091045506, rel=1e-6)
    np.testing.assert_allclose(
        res.t[:6],
        [
            10.0803398029743594,
            8.9516385580805142,
            8.4786308688438901,
            1.1305192636039123,
            1.3619185575345252,
            -0.0295519553554654,
        ],
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        res.p[:6],
        [
            4.07702270675825e-05,
            8.24836088966969e-05,
            1.13446427663349e-04,
            2.99474870830977e-01,
            2.19936960965250e-01,
            9.77338979892505e-01,
        ],
        rtol=1e-5,
    )
    np.testing.assert_allclose(res.q[:3], [0.00151261903551132] * 3, rtol=1e-5)
    np.testing.assert_allclose(
        res.ci_low[:3],
        [0.6289187719081164, 0.5684531550627743, 0.5430694876841933],
        rtol=1e-6,
    )


def test_moderated_one_sample_fails_on_nan() -> None:
    with pytest.raises(ValueError, match="non-finite"):
        pm.moderated_one_sample(np.array([[0.1, np.nan], [0.2, 0.3], [0.1, 0.2]]))


def test_pair_level_class_shift_planted() -> None:
    rng = np.random.default_rng(2)
    classes = np.array(["unmodified"] * 200 + ["CAM-only"] * 100)
    d = rng.normal(0, 0.2, size=(4, 300))
    d[:, 200:] -= 0.3
    out = pm.pair_level_class_shift(d, classes, "unmodified")
    r = out["CAM-only"]
    assert r["n_pairs_positive"] == 0
    ci_low, ci_high = r["ci_low"], r["ci_high"]
    assert isinstance(ci_low, float)
    assert isinstance(ci_high, float)
    assert ci_low < -0.3 < ci_high
    with pytest.raises(ValueError, match="baseline"):
        pm.pair_level_class_shift(d, classes, "oxidized")
