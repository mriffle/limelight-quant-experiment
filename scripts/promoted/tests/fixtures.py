"""Shared tiny hand-verified fixtures for the loader tests."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# 4 samples: 2 control, 2 raloxifene-d0; one batch of 4 (acquisition order 1..4).
SAMPLE_IDS = ["AZ001", "AZ002", "AZ003", "AZ004"]
SEARCH_IDS = ["101", "102", "103", "104"]
CONDITIONS = ["control", "raloxifene-d0", "control", "raloxifene-d0"]


def write_samples_tsv(path: Path) -> pd.DataFrame:
    """Write a tiny 4-sample ``samples.tsv`` fixture; return it as a DataFrame."""
    df = pd.DataFrame(
        {
            "sample_id": SAMPLE_IDS,
            "search_scan_file_id": SEARCH_IDS,
            "batch": ["B1", "B1", "B1", "B1"],
            "run_position_within_batch": ["1", "2", "3", "4"],
            "condition": CONDITIONS,
            "candidate_pair": ["P1", "P1", "P2", "P2"],
            "sample_role": ["experimental"] * 4,
        }
    )
    df.to_csv(path, sep="\t", index=False)
    return df


def intensity_col(search_id: str) -> str:
    return f"Intensity_search_scan_file_id_{search_id}"


def detection_col(search_id: str) -> str:
    return f"Detection Type_search_scan_file_id_{search_id}"
