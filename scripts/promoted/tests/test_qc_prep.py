"""Tests for qc_prep.py's pure helper functions (subset_features, run_batch_correction).

The full pipeline (loaders -> handle_missing -> normalize -> ComBat -> save_dataset)
is exercised end-to-end by running qc_prep.py on the real data (Stage-3 gate); these
tests cover the two units of non-trivial logic this script adds on top of the
already-tested loaders/normalize/batch_correct modules.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from loaders.data_loading import Dataset
from loaders.limelight_loader import LimelightProteinCounts
from qc_prep import (
    build_nsaf_dataset,
    build_psm_dataset,
    process_psm_level,
    run_batch_correction,
    subset_features,
)


def _dataset(scale: str = "linear") -> Dataset:
    abundances = np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
    )
    return Dataset(
        abundances=abundances,
        feature_names=np.array(["F1", "F2", "F3"], dtype=str),
        feature_metadata=pd.DataFrame({"is_contaminant": [False, True, False]}),
        metadata=pd.DataFrame(
            {"batch": ["B1", "B1", "B2", "B2"], "condition": ["a", "b", "a", "b"]},
            index=pd.Index(["s0", "s1", "s2", "s3"], name="sample_id"),
        ),
        scale=scale,  # type: ignore[arg-type]
    )


def test_subset_features_keeps_only_masked_columns() -> None:
    ds = _dataset()
    mask = np.array([True, False, True])
    subset = subset_features(ds, mask)
    assert subset.abundances.shape == (4, 2)
    np.testing.assert_array_equal(subset.abundances, ds.abundances[:, [0, 2]])
    assert list(subset.feature_names) == ["F1", "F3"]
    assert list(subset.feature_metadata["is_contaminant"]) == [False, False]
    # Independent of the input (no aliasing): mutating the subset must not touch ds.
    subset.abundances[0, 0] = 999.0
    assert ds.abundances[0, 0] == 1.0


def test_subset_features_shape_mismatch_fails_loud() -> None:
    ds = _dataset()
    with pytest.raises(ValueError, match="keep_mask shape"):
        subset_features(ds, np.array([True, False]))  # wrong length


def test_run_batch_correction_success() -> None:
    # A perfectly patterned ramp (as in _dataset()) is numerically degenerate for
    # ComBat's per-batch variance estimate with only 2 samples/batch (produces NaN);
    # use seeded random data instead, which is well-conditioned.
    rng = np.random.default_rng(20260922)
    abundances = rng.normal(loc=10.0, scale=1.0, size=(6, 20))
    ds = Dataset(
        abundances=abundances,
        feature_names=np.array([f"F{i}" for i in range(20)], dtype=str),
        feature_metadata=pd.DataFrame({"is_contaminant": [False] * 20}),
        metadata=pd.DataFrame(
            {"batch": ["B1", "B1", "B1", "B2", "B2", "B2"]},
            index=pd.Index([f"s{i}" for i in range(6)], name="sample_id"),
        ),
        scale="log2",
    )
    outcome = run_batch_correction(ds)
    assert outcome.status == "ok"
    assert outcome.dataset is not None
    assert outcome.dataset.abundances.shape == ds.abundances.shape
    assert outcome.error is None


def _protein_feature_metadata() -> pd.DataFrame:
    """3 protein-quants rows: 2 normal, 1 pure contaminant (dump-refined already)."""
    return pd.DataFrame(
        {
            "protein_group": [
                "psvid_1_sp|ACC1|E1_HUMAN",
                "psvid_2_sp|ACC2|E2_HUMAN",
                "psvid_3_sp|E3_HUMAN|",
            ],
            "is_contaminant": [False, False, True],
            "contaminant_grouped_with_real": [False, False, False],
        }
    )


def _limelight(
    nsaf: np.ndarray | None = None,
    psms: np.ndarray | None = None,
    sample_ids: tuple[str, ...] = ("S1", "S2"),
) -> LimelightProteinCounts:
    """3 dump groups: 2 single-member (1:1 with protein-quants), 1 comma-joined
    (first member ACC2, matching protein-quants row 2; second member has no
    protein-quants counterpart, per the verified 1:1 relation)."""
    protein_group = np.array(
        ["sp|ACC1|E1_HUMAN", "sp|ACC2|E2_HUMAN,sp|ACCX|EX_HUMAN", "sp|E3_HUMAN|"],
        dtype=str,
    )
    nsaf = (
        nsaf if nsaf is not None else np.array([[0.01, 0.02, 0.0], [0.03, 0.0, 0.001]])
    )
    psms = psms if psms is not None else np.array([[5, 3, 0], [4, 0, 1]])
    return LimelightProteinCounts(
        protein_group=protein_group,
        member_accession_keys=[g.split(",") for g in protein_group],
        psms=psms.astype(np.int64),
        nsaf=nsaf.astype(float),
        sample_ids=np.array(sample_ids, dtype=str),
        label_map=pd.DataFrame({"label": ["l1", "l2"], "sample_id": list(sample_ids)}),
    )


def _metadata(sample_ids: tuple[str, ...] = ("S1", "S2")) -> pd.DataFrame:
    df = pd.DataFrame(
        {"sample_id": list(sample_ids), "batch": ["B1", "B1"], "condition": ["a", "b"]}
    )
    df.index = pd.Index(df["sample_id"], name="sample_id")
    return df


def test_build_nsaf_dataset_zero_to_nan_and_feature_metadata() -> None:
    limelight = _limelight()
    pfm = _protein_feature_metadata()
    ds = build_nsaf_dataset(limelight, pfm, _metadata())
    assert ds.abundances.shape == (2, 3)
    assert ds.scale == "linear"
    assert math.isnan(ds.abundances[0, 2])  # NSAF 0 -> NaN
    assert ds.abundances[0, 0] == pytest.approx(0.01)
    fm = ds.feature_metadata
    assert list(fm["first_member_id"]) == [
        "psvid_1_sp|ACC1|E1_HUMAN",
        "psvid_2_sp|ACC2|E2_HUMAN",
        "psvid_3_sp|E3_HUMAN|",
    ]
    assert list(fm["is_contaminant"]) == [False, False, True]


def test_build_psm_dataset_zero_to_nan_and_feature_metadata() -> None:
    limelight = _limelight()
    pfm = _protein_feature_metadata()
    ds = build_psm_dataset(limelight, pfm, _metadata())
    assert ds.abundances.shape == (2, 3)
    assert math.isnan(ds.abundances[1, 1])  # PSM 0 -> NaN
    assert ds.abundances[0, 0] == pytest.approx(5.0)
    assert list(ds.feature_metadata["is_contaminant"]) == [False, False, True]


def test_build_limelight_dataset_sample_order_mismatch_fails_loud() -> None:
    limelight = _limelight(sample_ids=("S1", "S2"))
    pfm = _protein_feature_metadata()
    reversed_metadata = _metadata(sample_ids=("S2", "S1"))
    with pytest.raises(ValueError, match="sample_id"):
        build_nsaf_dataset(limelight, pfm, reversed_metadata)


def test_process_psm_level_feature_set_mismatch_fails_loud(tmp_path: Path) -> None:
    """psm's own raw_linear_complete (PSM>=1 in every run) must match the nsaf
    level's; a deliberately wrong expected set must raise, not silently pass."""
    limelight = _limelight()
    pfm = _protein_feature_metadata()
    metadata = _metadata()
    wrong_expected = np.array(["not-a-real-group"], dtype=str)
    with pytest.raises(ValueError, match="does not match"):
        process_psm_level(limelight, pfm, metadata, tmp_path, wrong_expected)


def test_run_batch_correction_records_failure_not_swallowed() -> None:
    """A batch with < 2 samples makes combat_correct raise; the outcome must record
    the failure (status='failed', error set) rather than letting the exception
    propagate or silently working around it."""
    abundances = np.array([[1.0, 2.0], [4.0, 5.0], [7.0, 8.0]])
    ds = Dataset(
        abundances=abundances,
        feature_names=np.array(["F1", "F2"], dtype=str),
        feature_metadata=pd.DataFrame({"is_contaminant": [False, False]}),
        metadata=pd.DataFrame(
            {"batch": ["B1", "B2", "B2"]}, index=pd.Index(["s0", "s1", "s2"])
        ),
        scale="log2",
    )
    outcome = run_batch_correction(ds)
    assert outcome.status == "failed"
    assert outcome.dataset is None
    assert outcome.error is not None
    assert "1 sample" in outcome.error or ">= 2" in outcome.error
