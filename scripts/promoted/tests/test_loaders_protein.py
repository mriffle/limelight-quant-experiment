"""Tests for loaders/protein_loader.py: the two Stage-3 integrity obligations.

(A) test the loader — hand-verified fixtures, property/invariant checks, edge cases.
(B) verify-on-real-data lives in verify_loaders.py, not here.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from fixtures import (
    CONDITIONS,
    SAMPLE_IDS,
    SEARCH_IDS,
    intensity_col,
    write_samples_tsv,
)
from loaders.protein_loader import load_protein_dataset, refine_contaminants


def _write_protein_quants(
    path: Path, rows: list[dict[str, str]], search_ids: list[str] = SEARCH_IDS
) -> None:
    columns = [
        "Protein Groups",
        "Gene Name",
        "Organism",
        *[intensity_col(s) for s in search_ids],
    ]
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row[c] for c in columns) + "\n")


def _basic_rows() -> list[dict[str, str]]:
    return [
        {
            "Protein Groups": "psvid_1_sp|P00001|AAA_HUMAN",
            "Gene Name": "",
            "Organism": "",
            intensity_col("101"): "100.0",
            intensity_col("102"): "0",
            intensity_col("103"): "NaN",
            intensity_col("104"): "400.0",
        },
        {
            "Protein Groups": "psvid_2_sp|BBB_HUMAN|",  # contaminant, no accession
            "Gene Name": "",
            "Organism": "",
            intensity_col("101"): "10.0",
            intensity_col("102"): "20.0",
            intensity_col("103"): "30.0",
            intensity_col("104"): "40.0",
        },
        {
            "Protein Groups": "psvid_4_sp|P00004|DDD_HUMAN",  # all-missing feature
            "Gene Name": "",
            "Organism": "",
            intensity_col("101"): "0",
            intensity_col("102"): "0",
            intensity_col("103"): "NaN",
            intensity_col("104"): "0",
        },
    ]


@pytest.fixture
def basic_project(tmp_path: Path) -> tuple[Path, Path]:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    _write_protein_quants(quants_path, _basic_rows())
    return quants_path, samples_path


def test_shape_and_order(basic_project: tuple[Path, Path]) -> None:
    quants_path, samples_path = basic_project
    result = load_protein_dataset(quants_path, samples_path)
    ds = result.dataset
    assert ds.abundances.shape == (4, 3)
    assert ds.scale == "linear"
    assert list(ds.feature_names) == [
        "psvid_1_sp|P00001|AAA_HUMAN",
        "psvid_2_sp|BBB_HUMAN|",
        "psvid_4_sp|P00004|DDD_HUMAN",
    ]
    # Samples ordered by acquisition order (batch, run_position): AZ001..AZ004.
    assert list(ds.metadata["sample_id"]) == SAMPLE_IDS
    assert list(ds.metadata["condition"]) == CONDITIONS


def test_nan_and_zero_tokens_both_become_nan_but_stay_distinguishable(
    basic_project: tuple[Path, Path],
) -> None:
    quants_path, samples_path = basic_project
    result = load_protein_dataset(quants_path, samples_path)
    ds = result.dataset
    row0 = ds.abundances[:, 0]  # feature AAA_HUMAN, samples AZ001..AZ004
    assert row0[0] == pytest.approx(100.0)
    assert math.isnan(row0[1])  # raw "0"
    assert math.isnan(row0[2])  # raw "NaN"
    assert row0[3] == pytest.approx(400.0)

    assert result.was_zero_token[1, 0]
    assert not result.was_nan_token[1, 0]
    assert result.was_nan_token[2, 0]
    assert not result.was_zero_token[2, 0]
    # Positive cells untouched, bit-for-bit.
    assert not result.was_nan_token[0, 0]
    assert not result.was_zero_token[0, 0]
    assert not result.was_nan_token[3, 0]
    assert not result.was_zero_token[3, 0]


def test_all_missing_feature(basic_project: tuple[Path, Path]) -> None:
    quants_path, samples_path = basic_project
    result = load_protein_dataset(quants_path, samples_path)
    ds = result.dataset
    col = list(ds.feature_names).index("psvid_4_sp|P00004|DDD_HUMAN")
    assert np.isnan(ds.abundances[:, col]).all()


def test_contaminant_flag_and_accession_parsing(
    basic_project: tuple[Path, Path],
) -> None:
    quants_path, samples_path = basic_project
    result = load_protein_dataset(quants_path, samples_path)
    fm = result.dataset.feature_metadata.set_index("protein_group")
    normal = fm.loc["psvid_1_sp|P00001|AAA_HUMAN"]
    assert not normal["is_contaminant"]
    assert normal["accession"] == "P00001"
    assert normal["entry"] == "AAA_HUMAN"
    assert normal["psvid"] == 1

    contaminant = fm.loc["psvid_2_sp|BBB_HUMAN|"]
    assert contaminant["is_contaminant"]
    assert contaminant["accession"] == ""
    assert contaminant["entry"] == "BBB_HUMAN"
    # Without protein_limelight_file, no dump refinement is applied: the flag
    # column exists but is never set.
    assert not fm["contaminant_grouped_with_real"].any()


def _write_dump_groups(path: Path, groups: list[str]) -> None:
    """Write a minimal dump file with only the 'Protein(s)' column -- enough for
    refine_contaminants(), which only reads that one column."""
    with path.open("w", encoding="utf-8") as fh:
        fh.write("Protein(s)\n")
        for group in groups:
            fh.write(group + "\n")


def test_refine_contaminants_unit(tmp_path: Path) -> None:
    """Direct unit test of refine_contaminants() on a hand-built feature_metadata."""
    feature_metadata = pd.DataFrame(
        {
            "protein_group": [
                "psvid_10_sp|PURE_HUMAN|",
                "psvid_11_sp|GROUPED_HUMAN|",
                "psvid_12_sp|P99999|NORMAL_HUMAN",
            ],
            "is_contaminant": [True, True, False],
        }
    )
    dump_path = tmp_path / "dump.txt"
    _write_dump_groups(
        dump_path,
        [
            "sp|PURE_HUMAN|",  # single member -> pure contaminant, no refinement
            "sp|GROUPED_HUMAN|,sp|P12345|GROUPED_HUMAN",  # real sibling -> kept
            "sp|P99999|NORMAL_HUMAN",  # not a contaminant candidate; untouched
        ],
    )
    refined = refine_contaminants(feature_metadata, dump_path)
    by_id = refined.set_index("protein_group")
    assert by_id.loc["psvid_10_sp|PURE_HUMAN|", "is_contaminant"]
    assert not by_id.loc["psvid_10_sp|PURE_HUMAN|", "contaminant_grouped_with_real"]
    assert not by_id.loc["psvid_11_sp|GROUPED_HUMAN|", "is_contaminant"]
    assert by_id.loc["psvid_11_sp|GROUPED_HUMAN|", "contaminant_grouped_with_real"]
    assert not by_id.loc["psvid_12_sp|P99999|NORMAL_HUMAN", "is_contaminant"]
    assert not by_id.loc[
        "psvid_12_sp|P99999|NORMAL_HUMAN", "contaminant_grouped_with_real"
    ]
    # The input is not mutated (independent output): row 1 (GROUPED_HUMAN) flips to
    # False in the refined copy but the original DataFrame keeps its original True.
    assert feature_metadata.loc[1, "is_contaminant"]
    assert not refined.loc[1, "is_contaminant"]


def test_refine_contaminants_missing_dump_row_fails_loud(tmp_path: Path) -> None:
    feature_metadata = pd.DataFrame(
        {
            "protein_group": ["psvid_2_sp|BBB_HUMAN|"],
            "is_contaminant": [True],
        }
    )
    dump_path = tmp_path / "dump.txt"
    _write_dump_groups(dump_path, ["sp|SOMETHING_ELSE_HUMAN|"])  # no matching row
    with pytest.raises(ValueError, match="No Limelight dump row"):
        refine_contaminants(feature_metadata, dump_path)


def test_load_protein_dataset_applies_dump_refinement(tmp_path: Path) -> None:
    """End-to-end: load_protein_dataset(protein_limelight_file=...) applies the
    scientist-confirmed refinement (33-true-contaminant rule)."""
    quants_path = tmp_path / "protein-quants.tsv"
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    rows = [
        {
            "Protein Groups": "psvid_10_sp|PURE_HUMAN|",
            "Gene Name": "",
            "Organism": "",
            **{intensity_col(s): "100.0" for s in SEARCH_IDS},
        },
        {
            "Protein Groups": "psvid_11_sp|GROUPED_HUMAN|",
            "Gene Name": "",
            "Organism": "",
            **{intensity_col(s): "200.0" for s in SEARCH_IDS},
        },
    ]
    _write_protein_quants(quants_path, rows)
    dump_path = tmp_path / "dump.txt"
    _write_dump_groups(
        dump_path,
        ["sp|PURE_HUMAN|", "sp|GROUPED_HUMAN|,sp|P12345|GROUPED_HUMAN"],
    )
    result = load_protein_dataset(
        quants_path, samples_path, protein_limelight_file=dump_path
    )
    fm = result.dataset.feature_metadata.set_index("protein_group")
    assert fm.loc["psvid_10_sp|PURE_HUMAN|", "is_contaminant"]
    assert not fm.loc["psvid_11_sp|GROUPED_HUMAN|", "is_contaminant"]
    assert fm.loc["psvid_11_sp|GROUPED_HUMAN|", "contaminant_grouped_with_real"]


def test_positive_values_untouched_bit_for_bit(
    basic_project: tuple[Path, Path],
) -> None:
    quants_path, samples_path = basic_project
    result = load_protein_dataset(quants_path, samples_path)
    col = list(result.dataset.feature_names).index("psvid_2_sp|BBB_HUMAN|")
    np.testing.assert_array_equal(
        result.dataset.abundances[:, col], np.array([10.0, 20.0, 30.0, 40.0])
    )


def test_no_value_appears_that_was_not_in_source(
    basic_project: tuple[Path, Path],
) -> None:
    """Property: every finite abundance value equals some raw source cell exactly."""
    quants_path, samples_path = basic_project
    result = load_protein_dataset(quants_path, samples_path)
    raw = pd.read_csv(quants_path, sep="\t", dtype=str)
    raw_values = {
        float(v)
        for v in raw[[intensity_col(s) for s in SEARCH_IDS]].to_numpy().ravel()
        if v not in ("NaN",)
    }
    finite = result.dataset.abundances[np.isfinite(result.dataset.abundances)]
    assert set(finite.tolist()) <= raw_values


def test_counts_preserved(basic_project: tuple[Path, Path]) -> None:
    """Property: n_features/n_samples loaded == source counts."""
    quants_path, samples_path = basic_project
    result = load_protein_dataset(quants_path, samples_path)
    n_source_rows = sum(1 for _ in quants_path.open()) - 1
    assert result.dataset.abundances.shape[1] == n_source_rows
    assert result.dataset.abundances.shape[0] == len(SAMPLE_IDS)


def test_duplicate_feature_id_fails_loud(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    rows = _basic_rows()[:1] * 2  # duplicate the same id
    _write_protein_quants(quants_path, rows)
    with pytest.raises(ValueError, match="duplicate"):
        load_protein_dataset(quants_path, samples_path)


def test_orphan_sample_column_fails_loud(tmp_path: Path) -> None:
    """A data column whose search_scan_file_id is not in samples.tsv must raise."""
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    extra_ids = [*SEARCH_IDS, "999"]
    rows = _basic_rows()
    for row in rows:
        row[intensity_col("999")] = "1.0"
    _write_protein_quants(quants_path, rows, search_ids=extra_ids)
    with pytest.raises(ValueError, match="bijection"):
        load_protein_dataset(quants_path, samples_path)


def test_missing_sample_column_fails_loud(tmp_path: Path) -> None:
    """A sample in samples.tsv with no matching data column must raise."""
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    fewer_ids = SEARCH_IDS[:-1]
    rows = _basic_rows()
    for row in rows:
        del row[intensity_col(SEARCH_IDS[-1])]
    _write_protein_quants(quants_path, rows, search_ids=fewer_ids)
    with pytest.raises(ValueError, match="bijection"):
        load_protein_dataset(quants_path, samples_path)


def test_duplicate_header_column_fails_loud(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    header = [
        "Protein Groups",
        "Gene Name",
        "Organism",
        *[intensity_col(s) for s in SEARCH_IDS],
        intensity_col(SEARCH_IDS[0]),  # duplicate column name
    ]
    with quants_path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write(
            "\t".join(["psvid_1_sp|P00001|AAA_HUMAN", "", "", "1", "2", "3", "4", "5"])
            + "\n"
        )
    with pytest.raises(ValueError, match="duplicate column header"):
        load_protein_dataset(quants_path, samples_path)


def test_empty_file_fails_loud(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    header = [
        "Protein Groups",
        "Gene Name",
        "Organism",
        *[intensity_col(s) for s in SEARCH_IDS],
    ]
    quants_path.write_text("\t".join(header) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no data rows"):
        load_protein_dataset(quants_path, samples_path)


def test_negative_value_fails_loud(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    rows = _basic_rows()
    rows[0][intensity_col("101")] = "-5.0"
    _write_protein_quants(quants_path, rows)
    with pytest.raises(ValueError, match="negative"):
        load_protein_dataset(quants_path, samples_path)


def test_non_numeric_non_nan_value_fails_loud(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    rows = _basic_rows()
    rows[0][intensity_col("101")] = "Filtered"
    _write_protein_quants(quants_path, rows)
    with pytest.raises(ValueError, match="non-numeric"):
        load_protein_dataset(quants_path, samples_path)


def test_single_sample(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    df = pd.DataFrame(
        {
            "sample_id": ["AZ001"],
            "search_scan_file_id": ["101"],
            "batch": ["B1"],
            "run_position_within_batch": ["1"],
            "condition": ["control"],
        }
    )
    df.to_csv(samples_path, sep="\t", index=False)
    quants_path = tmp_path / "protein-quants.tsv"
    quants_path.write_text(
        "Protein Groups\tGene Name\tOrganism\t" + intensity_col("101") + "\n"
        "psvid_1_sp|P00001|AAA_HUMAN\t\t\t123.0\n",
        encoding="utf-8",
    )
    result = load_protein_dataset(quants_path, samples_path)
    assert result.dataset.abundances.shape == (1, 1)
    assert result.dataset.abundances[0, 0] == pytest.approx(123.0)


def test_missing_data_file_raises(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    with pytest.raises(FileNotFoundError):
        load_protein_dataset(tmp_path / "does-not-exist.tsv", samples_path)


def test_planted_truth_2x_fold_change_recovered(tmp_path: Path) -> None:
    """Planted-truth check: loader + median-normalize must recover log2FC ~= 1.

    4 samples (2 control @ baseline, 2 raloxifene-d0 @ baseline*2, for one feature);
    every other feature is held flat across samples so median normalization is a
    no-op (all per-sample medians equal), isolating the planted effect.
    """
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    control_ids = ["101", "103"]  # AZ001, AZ003 (see fixtures.CONDITIONS)
    ralox_ids = ["102", "104"]  # AZ002, AZ004
    rows = []
    # The planted feature: doubled in raloxifene-d0.
    rows.append(
        {
            "Protein Groups": "psvid_1_sp|P00001|PLANTED_HUMAN",
            "Gene Name": "",
            "Organism": "",
            **{intensity_col(s): "1000.0" for s in control_ids},
            **{intensity_col(s): "2000.0" for s in ralox_ids},
        }
    )
    # Ten flat filler features (same value in every sample) so the per-sample
    # median is identical across samples and median normalization is a no-op.
    for i in range(10):
        rows.append(
            {
                "Protein Groups": f"psvid_{i + 2}_sp|P{i:05d}|FLAT{i}_HUMAN",
                "Gene Name": "",
                "Organism": "",
                **{intensity_col(s): "500.0" for s in SEARCH_IDS},
            }
        )
    quants_path = tmp_path / "protein-quants.tsv"
    _write_protein_quants(quants_path, rows)

    from loaders.normalize import normalize

    result = load_protein_dataset(quants_path, samples_path)
    normalized = normalize(result.dataset, method="median")
    # log2FC computed directly from the normalized linear values (no pseudocount
    # bias): median normalization is a no-op here (every per-sample median is
    # identical, from the flat filler features), so this isolates the planted 2x.
    col = list(normalized.feature_names).index("psvid_1_sp|P00001|PLANTED_HUMAN")
    control_mask = normalized.metadata["condition"] == "control"
    ralox_mask = normalized.metadata["condition"] == "raloxifene-d0"
    control_mean = normalized.abundances[control_mask.to_numpy(), col].mean()
    ralox_mean = normalized.abundances[ralox_mask.to_numpy(), col].mean()
    log2fc = math.log2(ralox_mean / control_mean)
    assert log2fc == pytest.approx(1.0, abs=1e-6)
