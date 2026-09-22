"""Tests for ``scripts/scratch/common/design.py``."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from common import design


def test_run_half_label_even_batch_splits_evenly() -> None:
    # batch_size=4 -> half_size=2: positions 1,2 early; 3,4 late.
    assert [design.run_half_label(p, 4) for p in (1, 2, 3, 4)] == [
        "early",
        "early",
        "late",
        "late",
    ]


def test_run_half_label_odd_batch_middle_goes_early() -> None:
    # batch_size=5 -> half_size=ceil(5/2)=3: positions 1,2,3 early; 4,5 late.
    # The middle sample (position 3) is the load-bearing case: it must be
    # "early", not "late".
    assert [design.run_half_label(p, 5) for p in (1, 2, 3, 4, 5)] == [
        "early",
        "early",
        "early",
        "late",
        "late",
    ]


def test_run_half_label_single_sample_batch_is_early() -> None:
    assert design.run_half_label(1, 1) == "early"


def test_run_half_label_rejects_nonpositive_batch_size() -> None:
    with pytest.raises(ValueError, match="positive"):
        design.run_half_label(1, 0)


@pytest.mark.parametrize(
    ("position", "batch_size"),
    [(0, 5), (6, 5), (-1, 5)],
)
def test_run_half_label_rejects_out_of_range_position(
    position: int, batch_size: int
) -> None:
    with pytest.raises(ValueError, match="out of range"):
        design.run_half_label(position, batch_size)


@given(st.integers(min_value=1, max_value=50))
def test_run_half_label_exactly_ceil_half_are_early(batch_size: int) -> None:
    """Property: exactly ceil(batch_size / 2) positions are labeled "early",
    for every batch size, and "early" positions are always a prefix (1..k)."""
    labels = [design.run_half_label(p, batch_size) for p in range(1, batch_size + 1)]
    expected_early = -(-batch_size // 2)
    assert labels.count("early") == expected_early
    assert labels.count("late") == batch_size - expected_early
    # "early" positions form a contiguous prefix.
    assert labels == ["early"] * expected_early + ["late"] * (
        batch_size - expected_early
    )


def test_failure_marker_name_is_a_json_filename() -> None:
    assert design.FAILURE_MARKER_NAME == "FAILED.json"
