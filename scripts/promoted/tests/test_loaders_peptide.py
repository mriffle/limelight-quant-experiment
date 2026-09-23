"""Tests for loaders/peptide_loader.py."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from fixtures import (
    CONDITIONS,
    SAMPLE_IDS,
    SEARCH_IDS,
    detection_col,
    intensity_col,
    write_samples_tsv,
)
from loaders.peptide_loader import load_peptide_dataset


def _write_peptide_quants(path: Path, rows: list[dict[str, str]]) -> None:
    columns = [
        "Sequence",
        "Base Sequence",
        "Protein Groups",
        "Gene Names",
        "Organism",
        *[intensity_col(s) for s in SEARCH_IDS],
        *[detection_col(s) for s in SEARCH_IDS],
    ]
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row[c] for c in columns) + "\n")


def _row(
    sequence: str,
    base: str,
    protein_groups: str,
    intensities: dict[str, str],
    types: dict[str, str],
) -> dict[str, str]:
    row = {
        "Sequence": sequence,
        "Base Sequence": base,
        "Protein Groups": protein_groups,
        "Gene Names": "",
        "Organism": "",
    }
    for sid in SEARCH_IDS:
        row[intensity_col(sid)] = intensities[sid]
        row[detection_col(sid)] = types[sid]
    return row


def _basic_rows() -> list[dict[str, str]]:
    return [
        _row(
            "AAAK[+0.0]",
            "AAAK",
            "psvid_1_sp|P00001|AAA_HUMAN",
            {"101": "100.0", "102": "0", "103": "50.0", "104": "0"},
            {
                "101": "MSMS",
                "102": "NotDetected",
                "103": "MBR",
                "104": "MSMSIdentifiedButNotQuantified",
            },
        ),
        _row(
            "CCCK[+57.02]",
            "CCCK",
            "psvid_2_sp|BBB_HUMAN|",  # contaminant single member
            {"101": "10.0", "102": "20.0", "103": "30.0", "104": "40.0"},
            {"101": "MSMS", "102": "MSMS", "103": "MSMS", "104": "MSMS"},
        ),
        _row(
            "DDDK[+0.0]",
            "DDDK",
            "psvid_3_sp|P00003|CCC_HUMAN;psvid_4_sp|EEE_HUMAN|",  # shared;1 contam.
            {"101": "0", "102": "0", "103": "0", "104": "0"},
            {
                "101": "NotDetected",
                "102": "NotDetected",
                "103": "NotDetected",
                "104": "MSMSAmbiguousPeakfinding",
            },
        ),
    ]


@pytest.fixture
def basic_project(tmp_path: Path) -> tuple[Path, Path]:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "peptide-quants.tsv"
    _write_peptide_quants(quants_path, _basic_rows())
    return quants_path, samples_path


def test_shape_and_order(basic_project: tuple[Path, Path]) -> None:
    quants_path, samples_path = basic_project
    result = load_peptide_dataset(quants_path, samples_path)
    ds = result.dataset
    assert ds.abundances.shape == (4, 3)
    assert list(ds.feature_names) == ["AAAK[+0.0]", "CCCK[+57.02]", "DDDK[+0.0]"]
    assert list(ds.metadata["sample_id"]) == SAMPLE_IDS
    assert list(ds.metadata["condition"]) == CONDITIONS


def test_zero_becomes_nan_positive_untouched(basic_project: tuple[Path, Path]) -> None:
    quants_path, samples_path = basic_project
    result = load_peptide_dataset(quants_path, samples_path)
    col = list(result.dataset.feature_names).index("AAAK[+0.0]")
    row = result.dataset.abundances[:, col]
    assert row[0] == pytest.approx(100.0)
    assert math.isnan(row[1])
    assert row[2] == pytest.approx(50.0)
    assert math.isnan(row[3])


def test_all_zero_row_is_all_nan(basic_project: tuple[Path, Path]) -> None:
    quants_path, samples_path = basic_project
    result = load_peptide_dataset(quants_path, samples_path)
    col = list(result.dataset.feature_names).index("DDDK[+0.0]")
    assert np.isnan(result.dataset.abundances[:, col]).all()


def test_detection_type_matrix_returned(basic_project: tuple[Path, Path]) -> None:
    quants_path, samples_path = basic_project
    result = load_peptide_dataset(quants_path, samples_path)
    col = list(result.dataset.feature_names).index("AAAK[+0.0]")
    assert list(result.detection_type[:, col]) == [
        "MSMS",
        "NotDetected",
        "MBR",
        "MSMSIdentifiedButNotQuantified",
    ]


def test_mod_mass_and_contaminant_and_base_sequence(
    basic_project: tuple[Path, Path],
) -> None:
    quants_path, samples_path = basic_project
    result = load_peptide_dataset(quants_path, samples_path)
    fm = result.dataset.feature_metadata.set_index("sequence")
    assert fm.loc["AAAK[+0.0]", "total_mod_mass"] == pytest.approx(0.0)
    assert fm.loc["CCCK[+57.02]", "total_mod_mass"] == pytest.approx(57.02)
    assert fm.loc["AAAK[+0.0]", "base_sequence"] == "AAAK"
    assert not fm.loc["AAAK[+0.0]", "is_contaminant"]
    assert fm.loc["CCCK[+57.02]", "is_contaminant"]  # single contaminant member
    assert fm.loc["DDDK[+0.0]", "is_contaminant"]  # one of two ';'-joined members


def test_contaminant_ids_uses_refined_set_not_coarse_regex(
    basic_project: tuple[Path, Path],
) -> None:
    """With ``contaminant_ids`` given, only members literally IN that set count --
    a no-accession-form id that the dump refined to "kept" (not in the pure set,
    e.g. because it turned out grouped_with_real) must NOT flag a peptide as
    contaminant, even though it still matches the coarser regex."""
    quants_path, samples_path = basic_project
    # "psvid_4_sp|EEE_HUMAN|" is simulated as dump-refined to "kept" (NOT pure);
    # only "psvid_2_sp|BBB_HUMAN|" is a true (pure) contaminant here.
    contaminant_ids = frozenset({"psvid_2_sp|BBB_HUMAN|"})
    result = load_peptide_dataset(
        quants_path, samples_path, contaminant_ids=contaminant_ids
    )
    fm = result.dataset.feature_metadata.set_index("sequence")
    assert fm.loc["CCCK[+57.02]", "is_contaminant"]  # its sole member is pure
    assert not fm.loc["DDDK[+0.0]", "is_contaminant"]  # only member is NOT pure
    assert not fm.loc["AAAK[+0.0]", "is_contaminant"]  # no contaminant member at all


def test_detection_type_invariant_violation_fails_loud(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "peptide-quants.tsv"
    rows = _basic_rows()
    # Violate: positive intensity but a zero-type Detection Type.
    rows[0]["Detection Type_search_scan_file_id_101"] = "NotDetected"
    _write_peptide_quants(quants_path, rows)
    with pytest.raises(ValueError, match="invariant"):
        load_peptide_dataset(quants_path, samples_path)


def test_unknown_detection_type_fails_loud(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "peptide-quants.tsv"
    rows = _basic_rows()
    rows[0]["Detection Type_search_scan_file_id_101"] = "SomethingElse"
    _write_peptide_quants(quants_path, rows)
    with pytest.raises(ValueError, match="undocumented Detection Type"):
        load_peptide_dataset(quants_path, samples_path)


def test_sequence_base_mismatch_fails_loud(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "peptide-quants.tsv"
    rows = _basic_rows()
    rows[0]["Base Sequence"] = "WRONG"
    _write_peptide_quants(quants_path, rows)
    with pytest.raises(ValueError, match="Base Sequence"):
        load_peptide_dataset(quants_path, samples_path)


def test_literal_nan_token_fails_loud(tmp_path: Path) -> None:
    """peptide-quants documents only '0' as missing; a literal 'NaN' is unexpected."""
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "peptide-quants.tsv"
    rows = _basic_rows()
    rows[0]["Intensity_search_scan_file_id_101"] = "NaN"
    rows[0]["Detection Type_search_scan_file_id_101"] = "NotDetected"
    _write_peptide_quants(quants_path, rows)
    with pytest.raises(ValueError, match="literal 'NaN'"):
        load_peptide_dataset(quants_path, samples_path)


def test_duplicate_sequence_fails_loud(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "peptide-quants.tsv"
    rows = _basic_rows()[:1] * 2
    _write_peptide_quants(quants_path, rows)
    with pytest.raises(ValueError, match="duplicate"):
        load_peptide_dataset(quants_path, samples_path)


def test_empty_file_fails_loud(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "peptide-quants.tsv"
    columns = [
        "Sequence",
        "Base Sequence",
        "Protein Groups",
        "Gene Names",
        "Organism",
        *[intensity_col(s) for s in SEARCH_IDS],
        *[detection_col(s) for s in SEARCH_IDS],
    ]
    quants_path.write_text("\t".join(columns) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no data rows"):
        load_peptide_dataset(quants_path, samples_path)


def test_counts_preserved(basic_project: tuple[Path, Path]) -> None:
    quants_path, samples_path = basic_project
    result = load_peptide_dataset(quants_path, samples_path)
    n_source_rows = sum(1 for _ in quants_path.open()) - 1
    assert result.dataset.abundances.shape[1] == n_source_rows
    assert result.dataset.abundances.shape[0] == len(SAMPLE_IDS)
