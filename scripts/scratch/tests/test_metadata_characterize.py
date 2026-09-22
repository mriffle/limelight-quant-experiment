"""Tests for ``scripts/scratch/metadata_characterize.py``.

Covers: filename-regex parsing (good + malformed + invalid-date), mapping-file
parsing, the bijection/fail-loud validity checks (synthetic fixtures, never the
real data) including the quant-header id-list checks (duplicates/length/leading
zeros) and interior-blank-line detection, unit tests for each hypothesis H1-H5
(including the inconclusive branch), the exact permutation p-value on
hand-checkable cases (3v3 and 2v2), a planted-truth check (fully confounded run
order recovered; interleaved design not flagged), a brute-force oracle for the
stratified permutation test, bias-corrected Cramer's V on a known table, the
role-keyed data-version stamp, the output-dir-under-data guard, and an
end-to-end ``main()`` check (outputs written, byte-identical across two runs,
failure marker written and no other outputs on a failing run).
"""

from __future__ import annotations

import itertools
import json
import sys
from collections.abc import Callable
from datetime import date
from pathlib import Path

import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import metadata_characterize as mc

# ---------------------------------------------------------------------------
# parse_filename: good + malformed filenames
# ---------------------------------------------------------------------------


def test_parse_filename_good() -> None:
    parsed = mc.parse_filename("UWPRExp480_2021_0506_AZ_034_AZ905_AZ_complex.mzML")
    assert parsed.year == 2021
    assert parsed.month == 5
    assert parsed.day == 6
    assert parsed.acq_date == date(2021, 5, 6)
    assert parsed.seq_number == 34
    assert parsed.sample_id == "AZ905"
    assert parsed.sample_number == 905


@pytest.mark.parametrize(
    "bad_name",
    [
        "UWPRExp480_2021_0506_AZ_034_AZ905_AZ_complex.raw",  # wrong extension
        "UWPRExp480_2021_506_AZ_034_AZ905_AZ_complex.mzML",  # MMDD not 4 digits
        "UWPRExp480_2021_0506_AZ_34_AZ905_AZ_complex.mzML",  # seq not 3 digits
        "UWPRExp480_2021_0506_AZ_034_905_AZ_complex.mzML",  # sample_id missing 'AZ'
        "UWPRExp999_2021_0506_AZ_034_AZ905_AZ_complex.mzML",  # wrong prefix
        "UWPRExp480_2021_0506_AZ_034_AZ905_AZ_complex.mzml",  # wrong case extension
        "",
    ],
)
def test_parse_filename_malformed_raises_pattern_error(bad_name: str) -> None:
    with pytest.raises(mc.FilenamePatternError):
        mc.parse_filename(bad_name)


def test_parse_filename_invalid_calendar_date_raises_date_error() -> None:
    # Month 13 does not exist.
    with pytest.raises(mc.FilenameDateError):
        mc.parse_filename("UWPRExp480_2021_1301_AZ_034_AZ905_AZ_complex.mzML")
    # Feb 30 does not exist.
    with pytest.raises(mc.FilenameDateError):
        mc.parse_filename("UWPRExp480_2021_0230_AZ_034_AZ905_AZ_complex.mzML")


# ---------------------------------------------------------------------------
# parse_mapping_file
# ---------------------------------------------------------------------------


def test_parse_mapping_file_good(tmp_path: Path) -> None:
    path = tmp_path / "mapping.txt"
    path.write_text("18458 (a.mzML), 18461 (b.mzML)", encoding="utf-8")
    entries = mc.parse_mapping_file(path)
    assert entries == [
        mc.MappingEntry(18458, "a.mzML"),
        mc.MappingEntry(18461, "b.mzML"),
    ]


def test_parse_mapping_file_empty_raises(tmp_path: Path) -> None:
    path = tmp_path / "mapping.txt"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="empty"):
        mc.parse_mapping_file(path)


def test_parse_mapping_file_malformed_entry_raises(tmp_path: Path) -> None:
    path = tmp_path / "mapping.txt"
    path.write_text("18458 a.mzML, 18461 (b.mzML)", encoding="utf-8")
    with pytest.raises(ValueError, match="does not match"):
        mc.parse_mapping_file(path)


# ---------------------------------------------------------------------------
# Synthetic fixtures for the fail-loud validity checks (never the real data)
# ---------------------------------------------------------------------------


_GOOD_FILES = [
    "UWPRExp480_2021_0506_AZ_034_AZ905_AZ_complex.mzML",
    "UWPRExp480_2021_0506_AZ_045_AZ906_AZ_complex.mzML",
]
_GOOD_CONDITIONS = ["control", "raloxifene-d0"]


def _write_fixture_files(
    tmp_path: Path,
    *,
    metadata_rows: list[tuple[str, str]],
    mapping_text: str,
    protein_quants_header: str = (
        "Protein Groups\tIntensity_search_scan_file_id_1\t"
        "Intensity_search_scan_file_id_2\n"
    ),
    peptide_quants_header: str = (
        "Sequence\tIntensity_search_scan_file_id_1\tIntensity_search_scan_file_id_2\t"
        "Detection Type_search_scan_file_id_1\tDetection Type_search_scan_file_id_2\n"
    ),
    protein_limelight_header: str = (
        "Protein(s)\tPSMs (1)\tPSMs (2)\tNSAF (1)\tNSAF (2)\n"
    ),
    peptide_limelight_header: str = (
        "Peptide Sequence\tPSMs (1)\tPSMs (2)\tQuant (1)\tQuant (2)\n"
    ),
) -> dict[str, Path]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    metadata_file = tmp_path / "metadata.tsv"
    lines = ["Replicate\tcondition"]
    lines.extend(f"{f}\t{c}" for f, c in metadata_rows)
    metadata_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    mapping_file = tmp_path / "mapping.txt"
    mapping_file.write_text(mapping_text, encoding="utf-8")

    protein_quants_file = tmp_path / "protein-quants.tsv"
    protein_quants_file.write_text(protein_quants_header, encoding="utf-8")

    peptide_quants_file = tmp_path / "peptide-quants.tsv"
    peptide_quants_file.write_text(peptide_quants_header, encoding="utf-8")

    protein_limelight_file = tmp_path / "protein-limelight.txt"
    protein_limelight_file.write_text(protein_limelight_header, encoding="utf-8")

    peptide_limelight_file = tmp_path / "peptide-limelight.txt"
    peptide_limelight_file.write_text(peptide_limelight_header, encoding="utf-8")

    return {
        "metadata_file": metadata_file,
        "mapping_file": mapping_file,
        "protein_quants_file": protein_quants_file,
        "peptide_quants_file": peptide_quants_file,
        "protein_limelight_file": protein_limelight_file,
        "peptide_limelight_file": peptide_limelight_file,
    }


def _two_sample_mapping() -> str:
    return f"1 ({_GOOD_FILES[0]}), 2 ({_GOOD_FILES[1]})"


_ARGV_FLAG_BY_ROLE = {
    "metadata_file": "--metadata-file",
    "mapping_file": "--mapping-file",
    "protein_quants_file": "--protein-quants-file",
    "peptide_quants_file": "--peptide-quants-file",
    "protein_limelight_file": "--protein-limelight-file",
    "peptide_limelight_file": "--peptide-limelight-file",
}


def _argv_from_paths(paths: dict[str, Path], output_dir: Path) -> list[str]:
    argv: list[str] = []
    for role, flag in _ARGV_FLAG_BY_ROLE.items():
        argv += [flag, str(paths[role])]
    argv += ["--output-dir", str(output_dir), "--log-level", "WARNING"]
    return argv


def test_validity_checks_pass_on_a_clean_2_sample_fixture(tmp_path: Path) -> None:
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=_two_sample_mapping(),
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    mc.raise_on_any_failure(checks)  # must not raise
    assert all(c.passed for c in checks)


def test_validity_checks_fail_loud_on_bijection_mismatch(tmp_path: Path) -> None:
    """Mapping file references a file not present in metadata.tsv -> raises."""
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=(
            f"1 ({_GOOD_FILES[0]}), "
            "2 (UWPRExp480_2021_0506_AZ_099_AZ999_AZ_complex.mzML)"
        ),
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    with pytest.raises(ValueError, match="validity check"):
        mc.raise_on_any_failure(checks)
    failed_names = {c.name for c in checks if not c.passed}
    assert "mapping_bijection_with_metadata_replicates" in failed_names


def test_validity_checks_fail_loud_on_duplicate_replicate(tmp_path: Path) -> None:
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=[
            (_GOOD_FILES[0], "control"),
            (_GOOD_FILES[0], "raloxifene-d0"),
        ],
        mapping_text=_two_sample_mapping(),
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    assert not next(c for c in checks if c.name == "replicate_unique").passed
    with pytest.raises(ValueError):
        mc.raise_on_any_failure(checks)


def test_validity_checks_fail_loud_on_bad_condition(tmp_path: Path) -> None:
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=[(_GOOD_FILES[0], "TREATED"), (_GOOD_FILES[1], "control")],
        mapping_text=_two_sample_mapping(),
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    assert not next(c for c in checks if c.name == "condition_in_allowed_set").passed
    with pytest.raises(ValueError):
        mc.raise_on_any_failure(checks)


def test_validity_checks_fail_loud_on_invalid_date(tmp_path: Path) -> None:
    """A replicate whose filename structurally matches but encodes month 13.

    The regex-match check must still PASS (the pattern's digit groups are all the
    right shape); only the date-validity check must fail, and it must report the
    actual bad row rather than being silently skipped.
    """
    bad_file = "UWPRExp480_2021_1301_AZ_034_AZ905_AZ_complex.mzML"
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=[(bad_file, "control"), (_GOOD_FILES[1], "raloxifene-d0")],
        mapping_text=f"1 ({bad_file}), 2 ({_GOOD_FILES[1]})",
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    by_name = {c.name: c for c in checks}
    assert by_name["replicate_matches_filename_regex"].passed
    assert not by_name["replicate_dates_valid_calendar_dates"].passed
    assert (
        "invalid calendar date"
        in by_name["replicate_dates_valid_calendar_dates"].detail
    )
    with pytest.raises(ValueError):
        mc.raise_on_any_failure(checks)


def test_validity_checks_fail_loud_on_mismatched_quant_ids(tmp_path: Path) -> None:
    """protein-quants header references an id not in the mapping file."""
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=_two_sample_mapping(),
        protein_quants_header=(
            "Protein Groups\tIntensity_search_scan_file_id_1\t"
            "Intensity_search_scan_file_id_999\n"
        ),
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    assert not next(
        c for c in checks if c.name == "protein_quants_intensity_ids_match_mapping"
    ).passed
    with pytest.raises(ValueError):
        mc.raise_on_any_failure(checks)


def test_validity_checks_fail_loud_on_duplicate_quant_id_column(tmp_path: Path) -> None:
    """Same id in two columns (a duplicated per-run column) must be caught as a
    LIST-level defect (length + duplicate-id), not silently absorbed by a set."""
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=_two_sample_mapping(),
        protein_quants_header=(
            "Protein Groups\tIntensity_search_scan_file_id_1\t"
            "Intensity_search_scan_file_id_1\n"
        ),
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    check = next(
        c for c in checks if c.name == "protein_quants_intensity_ids_match_mapping"
    )
    assert not check.passed
    assert "duplicate id" in check.detail
    with pytest.raises(ValueError):
        mc.raise_on_any_failure(checks)


def test_validity_checks_fail_loud_on_wrong_quant_column_count(tmp_path: Path) -> None:
    """3 Intensity columns for a 2-sample fixture (one id duplicated, extra
    column) must be caught by the length check even though the SET of ids
    still matches the mapping ids exactly."""
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=_two_sample_mapping(),
        protein_quants_header=(
            "Protein Groups\tIntensity_search_scan_file_id_1\t"
            "Intensity_search_scan_file_id_2\tIntensity_search_scan_file_id_2\n"
        ),
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    check = next(
        c for c in checks if c.name == "protein_quants_intensity_ids_match_mapping"
    )
    assert not check.passed
    assert "found 3 matching column(s), expected 2" in check.detail
    with pytest.raises(ValueError):
        mc.raise_on_any_failure(checks)


def test_validity_checks_fail_loud_on_non_canonical_quant_id_text(
    tmp_path: Path,
) -> None:
    """A leading-zero id (``_01``) must be rejected even though ``int('01')``
    happens to equal the right mapping id."""
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=_two_sample_mapping(),
        protein_quants_header=(
            "Protein Groups\tIntensity_search_scan_file_id_01\t"
            "Intensity_search_scan_file_id_2\n"
        ),
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    check = next(
        c for c in checks if c.name == "protein_quants_intensity_ids_match_mapping"
    )
    assert not check.passed
    assert "non-canonical" in check.detail
    assert "01" in check.detail
    with pytest.raises(ValueError):
        mc.raise_on_any_failure(checks)


def test_validity_checks_fail_loud_on_limelight_label_mismatch(tmp_path: Path) -> None:
    """peptide-limelight dump has only 1 label instead of matching protein's 2."""
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=_two_sample_mapping(),
        peptide_limelight_header="Peptide Sequence\tPSMs (1)\tQuant (1)\n",
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    names = {c.name for c in checks if not c.passed}
    assert "peptide_limelight_n_distinct_labels_per_family_consistent" in names
    assert "limelight_label_sets_match_across_files" in names
    with pytest.raises(ValueError):
        mc.raise_on_any_failure(checks)


def test_validity_checks_fail_loud_on_blank_cell(tmp_path: Path) -> None:
    """A blank ``condition`` cell must fail ``metadata_no_blank_cells`` and the
    overall run must raise (fail loud, not warn)."""
    metadata_file = tmp_path / "metadata.tsv"
    metadata_file.write_text(
        f"Replicate\tcondition\n{_GOOD_FILES[0]}\t\n{_GOOD_FILES[1]}\tcontrol\n",
        encoding="utf-8",
    )
    header, rows = mc.read_tsv_rows(metadata_file)
    mapping_file = tmp_path / "mapping.txt"
    mapping_file.write_text(_two_sample_mapping(), encoding="utf-8")
    raw = mc.RawInputs(
        metadata_header=header,
        metadata_rows=rows,
        mapping_entries=mc.parse_mapping_file(mapping_file),
        protein_quants_header=[
            "Protein Groups",
            "Intensity_search_scan_file_id_1",
            "Intensity_search_scan_file_id_2",
        ],
        peptide_quants_header=[
            "Sequence",
            "Intensity_search_scan_file_id_1",
            "Intensity_search_scan_file_id_2",
            "Detection Type_search_scan_file_id_1",
            "Detection Type_search_scan_file_id_2",
        ],
        protein_limelight_header=["Protein(s)", "PSMs (1)", "PSMs (2)"],
        peptide_limelight_header=["Peptide Sequence", "PSMs (1)", "PSMs (2)"],
    )
    checks = mc.run_validity_checks(raw)
    assert not next(c for c in checks if c.name == "metadata_no_blank_cells").passed
    with pytest.raises(ValueError, match="validity check"):
        mc.raise_on_any_failure(checks)


def test_read_tsv_rows_empty_file_raises(tmp_path: Path) -> None:
    path = tmp_path / "empty.tsv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="empty"):
        mc.read_tsv_rows(path)


def test_read_tsv_rows_trailing_blank_line_tolerated(tmp_path: Path) -> None:
    path = tmp_path / "metadata.tsv"
    path.write_text("Replicate\tcondition\nfileA\tcontrol\n\n", encoding="utf-8")
    header, rows = mc.read_tsv_rows(path)
    assert header == ["Replicate", "condition"]
    assert rows == [["fileA", "control"]]


def test_read_tsv_rows_interior_blank_line_raises(tmp_path: Path) -> None:
    path = tmp_path / "metadata.tsv"
    path.write_text(
        "Replicate\tcondition\nfileA\tcontrol\n\nfileB\tcontrol\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="blank line"):
        mc.read_tsv_rows(path)


def test_read_header_empty_file_raises(tmp_path: Path) -> None:
    path = tmp_path / "empty.tsv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="empty"):
        mc.read_header(path)


# ---------------------------------------------------------------------------
# extract_family_labels: ambiguous / truncated Limelight-dump-style headers
# ---------------------------------------------------------------------------


def test_extract_family_labels_handles_nested_parens() -> None:
    header = [
        "Protein(s)",
        "PSMs (1_0)",
        "PSMs (2_0318_)",
        "Shared group [column content string: (SHARED GROUP)] (1_0)",
        "Shared group [column content string: (SHARED GROUP)] (2_0318_)",
    ]
    families = mc.extract_family_labels(header)
    assert families["PSMs"] == ["1_0", "2_0318_"]
    assert families["Shared group [column content string: (SHARED GROUP)]"] == [
        "1_0",
        "2_0318_",
    ]
    # "Protein(s)" (an id column) has no trailing " (label)" and is skipped.
    assert "Protein(s)" not in families


def test_extract_ids_raw_preserves_duplicates_and_raw_text() -> None:
    header = [
        "Intensity_search_scan_file_id_1",
        "Intensity_search_scan_file_id_01",
        "Intensity_search_scan_file_id_2",
    ]
    raw = mc.extract_ids_raw(header, mc.INTENSITY_COL_RE)
    assert raw == [("1", 1), ("01", 1), ("2", 2)]


# ---------------------------------------------------------------------------
# Exact permutation test: hand-checkable cases (3v3 and 2v2)
# ---------------------------------------------------------------------------


def test_exact_permutation_3v3_fully_separated_one_sided_is_1_over_20() -> None:
    # Controls at the 3 smallest values, other group at the 3 largest -> full
    # separation. C(6,3) = 20 equally likely assignments; only 1 achieves this.
    seq_numbers = [10, 20, 30, 40, 50, 60]
    is_control = [True, True, True, False, False, False]
    result = mc.exact_permutation_test(seq_numbers, is_control)
    assert result.n_permutations == 20
    assert result.u_observed == 9
    assert result.u_max == 9
    assert result.p_one_sided_observed_direction == pytest.approx(1 / 20)
    assert result.p_two_sided == pytest.approx(2 / 20)
    assert result.rank_biserial == pytest.approx(1.0)
    assert result.direction == "control_earlier"


def test_exact_permutation_requires_both_groups_present() -> None:
    with pytest.raises(ValueError, match="at least one member"):
        mc.exact_permutation_test([1, 2, 3], [True, True, True])


def test_exact_permutation_rejects_ties() -> None:
    with pytest.raises(ValueError, match="distinct"):
        mc.exact_permutation_test([1, 1, 2], [True, False, True])


def test_exact_permutation_2v2_exact_distribution_and_p_values() -> None:
    # Interleaved order (not separated): control, other, control, other.
    seq_numbers = [10, 20, 30, 40]
    is_control = [True, False, True, False]
    result = mc.exact_permutation_test(seq_numbers, is_control)

    # Hand-enumerated U over all C(4,2)=6 combinations (control-position choice):
    # positions (0,1)->4, (0,2)->3, (0,3)->2, (1,2)->2, (1,3)->1, (2,3)->0.
    combos = list(itertools.combinations(range(4), 2))
    values = seq_numbers
    us = []
    for combo in combos:
        control_vals = [values[i] for i in combo]
        other_vals = [values[i] for i in range(4) if i not in combo]
        us.append(sum(1 for cv in control_vals for ov in other_vals if cv < ov))
    assert sorted(us) == [0, 1, 2, 2, 3, 4]

    assert result.n_permutations == 6
    assert result.u_observed == 3
    assert result.u_max == 4
    assert result.p_one_sided_observed_direction == pytest.approx(2 / 6)
    assert result.p_two_sided == pytest.approx(4 / 6)


def test_exact_stratified_permutation_matches_hand_calculation() -> None:
    # Batch 1: 3 control at smallest seq, 3 other at largest -> fully separated.
    # Batch 2: 1 control at smallest seq, 1 other at largest -> fully separated.
    batches = [
        ([10, 20, 30, 40, 50, 60], [True, True, True, False, False, False]),
        ([100, 200], [True, False]),
    ]
    result = mc.exact_stratified_permutation_test(batches)
    assert result.n_permutations == 20 * 2
    assert result.u_observed == 9 + 1
    assert result.u_max == 9 + 1
    assert result.p_one_sided_observed_direction == pytest.approx(1 / 40)
    assert result.p_two_sided == pytest.approx(2 / 40)


def test_exact_stratified_permutation_refuses_oversized_enumeration() -> None:
    batches = [([10, 20, 30, 40, 50, 60], [True, True, True, False, False, False])]
    with pytest.raises(ValueError, match="max_total_permutations"):
        mc.exact_stratified_permutation_test(batches, max_total_permutations=5)


def test_exact_stratified_permutation_requires_at_least_one_batch() -> None:
    with pytest.raises(ValueError, match="at least one batch"):
        mc.exact_stratified_permutation_test([])


def test_exact_stratified_permutation_single_batch_matches_single_batch_test() -> None:
    """A stratified test over exactly one batch must reproduce
    ``exact_permutation_test`` on that batch field-for-field."""
    seq_numbers = [3, 1, 9, 4, 2]
    is_control = [True, False, True, False, True]
    single = mc.exact_permutation_test(seq_numbers, is_control)
    stratified = mc.exact_stratified_permutation_test([(seq_numbers, is_control)])
    assert stratified.n_total == single.n_total
    assert stratified.n_control == single.n_control
    assert stratified.n_other == single.n_other
    assert stratified.u_observed == single.u_observed
    assert stratified.u_max == single.u_max
    assert stratified.n_permutations == single.n_permutations
    assert stratified.rank_biserial == pytest.approx(single.rank_biserial)
    assert stratified.direction == single.direction
    assert stratified.p_one_sided_observed_direction == pytest.approx(
        single.p_one_sided_observed_direction
    )
    assert stratified.p_two_sided == pytest.approx(single.p_two_sided)


def _brute_force_u(values: list[int], control_positions: tuple[int, ...]) -> int:
    """Independent (non-module) reference U statistic for the oracle test."""
    control_set = set(control_positions)
    control_vals = [values[i] for i in control_positions]
    other_vals = [values[i] for i in range(len(values)) if i not in control_set]
    return sum(1 for cv in control_vals for ov in other_vals if cv < ov)


def _brute_force_stratified(
    batches: list[tuple[list[int], list[bool]]],
) -> tuple[int, int, int, float, float]:
    """Independent (non-module) reference implementation of the stratified exact
    test, built from scratch (does not call any of the module's helpers) as an
    oracle to check ``exact_stratified_permutation_test`` against."""
    per_batch_all_u: list[list[int]] = []
    obs_total = 0
    max_total = 0
    n_total = 0
    for values, is_control in batches:
        n = len(values)
        n_control = sum(is_control)
        n_other = n - n_control
        obs_positions = tuple(i for i, c in enumerate(is_control) if c)
        obs_total += _brute_force_u(values, obs_positions)
        max_total += n_control * n_other
        n_total += n
        all_u = [
            _brute_force_u(values, combo)
            for combo in itertools.combinations(range(n), n_control)
        ]
        per_batch_all_u.append(all_u)

    joint_sums = [sum(combo) for combo in itertools.product(*per_batch_all_u)]
    n_perm = len(joint_sums)
    center2 = max_total
    obs_dist2 = abs(2 * obs_total - center2)
    p_two = sum(1 for u in joint_sums if abs(2 * u - center2) >= obs_dist2) / n_perm
    if 2 * obs_total >= center2:
        p_one = sum(1 for u in joint_sums if u >= obs_total) / n_perm
    else:
        p_one = sum(1 for u in joint_sums if u <= obs_total) / n_perm
    return obs_total, max_total, n_perm, p_one, p_two


@st.composite
def _small_multi_batch_case(draw: st.DrawFn) -> list[tuple[list[int], list[bool]]]:
    n_batches = draw(st.integers(min_value=1, max_value=2))
    batches: list[tuple[list[int], list[bool]]] = []
    offset = 0
    for _ in range(n_batches):
        n = draw(st.integers(min_value=2, max_value=4))
        n_control = draw(st.integers(min_value=1, max_value=n - 1))
        values = list(range(offset, offset + n))
        offset += n + 1
        positions = list(range(n))
        control_positions = draw(
            st.permutations(positions).map(_take_sorted_prefix(n_control))
        )
        is_control = [i in control_positions for i in positions]
        batches.append((values, is_control))
    return batches


def _take_sorted_prefix(k: int) -> Callable[[list[int]], list[int]]:
    """Bind ``k`` explicitly (avoids a late-binding closure over the loop var)."""

    def _take(values: list[int]) -> list[int]:
        return sorted(values[:k])

    return _take


@given(_small_multi_batch_case())
def test_exact_stratified_permutation_matches_brute_force_oracle(
    batches: list[tuple[list[int], list[bool]]],
) -> None:
    result = mc.exact_stratified_permutation_test(batches)
    obs, u_max, n_perm, p_one, p_two = _brute_force_stratified(batches)
    assert result.u_observed == obs
    assert result.u_max == u_max
    assert result.n_permutations == n_perm
    assert result.p_one_sided_observed_direction == pytest.approx(p_one)
    assert result.p_two_sided == pytest.approx(p_two)


# ---------------------------------------------------------------------------
# Cramer's V: known tables
# ---------------------------------------------------------------------------


def test_cramers_v_perfect_association_2x2_is_one() -> None:
    crosstab = pd.DataFrame([[4, 0], [0, 4]], index=["a", "b"], columns=["x", "y"])
    assert mc.cramers_v_uncorrected(crosstab) == pytest.approx(1.0)
    assert mc.cramers_v_bias_corrected(crosstab) == pytest.approx(1.0)


def test_cramers_v_independent_table_is_zero() -> None:
    # Proportional rows -> independence -> V = 0 exactly (both variants).
    crosstab = pd.DataFrame([[3, 1], [3, 1]], index=["a", "b"], columns=["x", "y"])
    assert mc.cramers_v_uncorrected(crosstab) == pytest.approx(0.0, abs=1e-9)
    assert mc.cramers_v_bias_corrected(crosstab) == pytest.approx(0.0, abs=1e-9)


def test_cramers_v_bias_corrected_below_raw_under_noisy_association() -> None:
    # A noisy (imperfect) 3x3 association: bias correction should pull V down
    # relative to the uncorrected estimate at this small n (textbook regime the
    # Bergsma correction targets).
    crosstab = pd.DataFrame(
        [[5, 2, 1], [1, 5, 2], [2, 1, 5]],
        index=["a", "b", "c"],
        columns=["x", "y", "z"],
    )
    raw = mc.cramers_v_uncorrected(crosstab)
    corrected = mc.cramers_v_bias_corrected(crosstab)
    assert corrected <= raw + 1e-9


def test_cramers_v_degenerate_table_is_zero() -> None:
    crosstab = pd.DataFrame([[4], [4]], index=["a", "b"], columns=["x"])
    assert mc.cramers_v_uncorrected(crosstab) == 0.0
    assert mc.cramers_v_bias_corrected(crosstab) == 0.0


# ---------------------------------------------------------------------------
# build_samples_table: derived fields on a small synthetic fixture
# ---------------------------------------------------------------------------


def test_build_samples_table_derives_expected_fields(tmp_path: Path) -> None:
    files = [
        "UWPRExp480_2021_0101_AZ_001_AZ101_AZ_complex.mzML",
        "UWPRExp480_2021_0101_AZ_002_AZ102_AZ_complex.mzML",
    ]
    conditions = ["control", "raloxifene-d0"]
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(files, conditions, strict=True)),
        mapping_text=f"1 ({files[0]}), 2 ({files[1]})",
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    mc.raise_on_any_failure(checks)
    samples = mc.build_samples_table(raw)

    assert list(samples["sample_id"]) == ["AZ101", "AZ102"]
    assert list(samples["candidate_pair"]) == ["P101_102", "P101_102"]
    assert list(samples["run_position_within_batch"]) == [1, 2]
    assert list(samples["run_half"]) == ["early", "late"]
    assert list(samples["sample_role"]) == ["experimental", "experimental"]
    assert list(samples["search_scan_file_id"]) == [1, 2]
    assert list(samples["intensity_column"]) == [
        "Intensity_search_scan_file_id_1",
        "Intensity_search_scan_file_id_2",
    ]


def test_build_run_layout_and_table1_shapes(tmp_path: Path) -> None:
    files = [
        "UWPRExp480_2021_0101_AZ_001_AZ101_AZ_complex.mzML",
        "UWPRExp480_2021_0101_AZ_002_AZ102_AZ_complex.mzML",
    ]
    conditions = ["control", "raloxifene-d0"]
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(files, conditions, strict=True)),
        mapping_text=f"1 ({files[0]}), 2 ({files[1]})",
    )
    raw = mc.load_raw_inputs(**paths)
    mc.raise_on_any_failure(mc.run_validity_checks(raw))
    samples = mc.build_samples_table(raw)

    run_layout = mc.build_run_layout(samples)
    assert run_layout.shape == (1, 3)  # 1 batch, batch-col + 2 run positions
    assert run_layout.loc[0, "1"] == "control/AZ101"
    assert run_layout.loc[0, "2"] == "raloxifene-d0/AZ102"

    table1 = mc.build_table1(samples)
    assert set(table1["condition"]) == {"control", "raloxifene-d0"}
    assert (table1["n"] == 1).all()

    join_key = mc.build_join_key(samples)
    assert (join_key["limelight_dump_label"] == "UNRESOLVED").all()
    assert set(join_key["sample_id"]) == {"AZ101", "AZ102"}


# ---------------------------------------------------------------------------
# Hypotheses H1-H5: direct unit tests against small hand-built samples tables
# ---------------------------------------------------------------------------


def _samples_df(rows: list[dict[str, object]]) -> pd.DataFrame:
    """Build a minimal samples-table-shaped DataFrame for hypothesis unit tests."""
    return pd.DataFrame(rows)


def test_hypothesis_h1_balanced_passes() -> None:
    # H1 checks the study's specific design target (4 vs 4), not "balanced" in
    # general, per the Stage-1 spec ("H1 balanced arms: 4 control, 4
    # raloxifene-d0") — matching the fixed 8-sample cohort this project has.
    df = _samples_df(
        [{"condition": c} for c in (["control"] * 4 + ["raloxifene-d0"] * 4)]
    )
    row = mc.hypothesis_h1(df)
    assert row.passed is True
    assert row.statistic_value == "4,4"


def test_hypothesis_h1_unbalanced_fails() -> None:
    df = _samples_df(
        [{"condition": c} for c in ["control", "control", "control", "raloxifene-d0"]]
    )
    row = mc.hypothesis_h1(df)
    assert row.passed is False


def test_hypothesis_h2_batch_aliased_with_condition_fails() -> None:
    # Batch B2 has only "control" -> condition IS aliased with batch there.
    df = _samples_df(
        [
            {"batch": "B1", "condition": "control"},
            {"batch": "B1", "condition": "raloxifene-d0"},
            {"batch": "B2", "condition": "control"},
            {"batch": "B2", "condition": "control"},
        ]
    )
    row = mc.hypothesis_h2(df)
    assert row.passed is False


def test_hypothesis_h2_every_batch_has_both_conditions_passes() -> None:
    df = _samples_df(
        [
            {"batch": "B1", "condition": "control"},
            {"batch": "B1", "condition": "raloxifene-d0"},
            {"batch": "B2", "condition": "control"},
            {"batch": "B2", "condition": "raloxifene-d0"},
        ]
    )
    row = mc.hypothesis_h2(df)
    assert row.passed is True


def test_hypothesis_h3_parity_and_pairing_passes() -> None:
    df = _samples_df(
        [
            {"candidate_pair": "P1_2", "condition": "control", "sample_number": 1},
            {
                "candidate_pair": "P1_2",
                "condition": "raloxifene-d0",
                "sample_number": 2,
            },
            {"candidate_pair": "P3_4", "condition": "control", "sample_number": 3},
            {
                "candidate_pair": "P3_4",
                "condition": "raloxifene-d0",
                "sample_number": 4,
            },
        ]
    )
    row = mc.hypothesis_h3(df)
    assert row.passed is True
    assert row.statistic_value == "True,True"


def test_hypothesis_h3_parity_violation_fails() -> None:
    # sample_number=1 (odd) is raloxifene-d0, not control -> parity rule broken.
    df = _samples_df(
        [
            {
                "candidate_pair": "P1_2",
                "condition": "raloxifene-d0",
                "sample_number": 1,
            },
            {"candidate_pair": "P1_2", "condition": "control", "sample_number": 2},
        ]
    )
    row = mc.hypothesis_h3(df)
    assert row.passed is False
    assert row.statistic_value == "True,False"


def test_hypothesis_h3_pair_composition_violation_fails() -> None:
    # Pair P1_2 has two controls, not one of each.
    df = _samples_df(
        [
            {"candidate_pair": "P1_2", "condition": "control", "sample_number": 1},
            {"candidate_pair": "P1_2", "condition": "control", "sample_number": 2},
        ]
    )
    row = mc.hypothesis_h3(df)
    assert row.passed is False
    assert row.statistic_value == "False,False"


def test_hypothesis_h4_inconclusive_when_batch_has_lt2_of_one_condition() -> None:
    df = _samples_df(
        [
            {"batch": "B1", "seq_number": 1, "condition": "control"},
            {"batch": "B1", "seq_number": 2, "condition": "raloxifene-d0"},
        ]
    )
    rows = mc.hypothesis_h4(df, max_total_permutations=1000)
    batch_row = next(r for r in rows if r.scope == "batch=B1")
    assert batch_row.passed is None
    assert batch_row.p_one_sided_observed_direction is None
    assert batch_row.p_two_sided is None
    assert "No information on its own" in batch_row.note


def test_hypothesis_h5_order_matches() -> None:
    df = _samples_df(
        [
            {"batch": "B1", "seq_number": 1, "sample_number": 101},
            {"batch": "B1", "seq_number": 2, "sample_number": 102},
        ]
    )
    rows = mc.hypothesis_h5(df)
    assert len(rows) == 1
    assert rows[0].passed is True


def test_hypothesis_h5_order_mismatches() -> None:
    df = _samples_df(
        [
            {"batch": "B1", "seq_number": 1, "sample_number": 102},
            {"batch": "B1", "seq_number": 2, "sample_number": 101},
        ]
    )
    rows = mc.hypothesis_h5(df)
    assert len(rows) == 1
    assert rows[0].passed is False


# ---------------------------------------------------------------------------
# Planted-truth check: fully confounded run order recovered; interleaved not
# ---------------------------------------------------------------------------


def test_h4_planted_truth_fully_confounded_run_order_is_flagged(
    tmp_path: Path,
) -> None:
    """4 samples, 1 batch: the two controls run first, the two raloxifene-d0
    run after -> a fully planted confound the pipeline must recover, with the
    exact p-value 1/C(4,2) = 1/6."""
    files = [
        "UWPRExp480_2021_0101_AZ_001_AZ201_AZ_complex.mzML",
        "UWPRExp480_2021_0101_AZ_002_AZ203_AZ_complex.mzML",
        "UWPRExp480_2021_0101_AZ_003_AZ202_AZ_complex.mzML",
        "UWPRExp480_2021_0101_AZ_004_AZ204_AZ_complex.mzML",
    ]
    conditions = ["control", "control", "raloxifene-d0", "raloxifene-d0"]
    mapping_text = ", ".join(f"{i + 1} ({f})" for i, f in enumerate(files))
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(files, conditions, strict=True)),
        mapping_text=mapping_text,
        protein_quants_header=(
            "Protein Groups\t"
            + "\t".join(f"Intensity_search_scan_file_id_{i + 1}" for i in range(4))
            + "\n"
        ),
        peptide_quants_header=(
            "Sequence\t"
            + "\t".join(f"Intensity_search_scan_file_id_{i + 1}" for i in range(4))
            + "\t"
            + "\t".join(f"Detection Type_search_scan_file_id_{i + 1}" for i in range(4))
            + "\n"
        ),
        protein_limelight_header=(
            "Protein(s)\t"
            + "\t".join(f"PSMs ({i + 1})" for i in range(4))
            + "\t"
            + "\t".join(f"NSAF ({i + 1})" for i in range(4))
            + "\n"
        ),
        peptide_limelight_header=(
            "Peptide Sequence\t"
            + "\t".join(f"PSMs ({i + 1})" for i in range(4))
            + "\t"
            + "\t".join(f"Quant ({i + 1})" for i in range(4))
            + "\n"
        ),
    )
    raw = mc.load_raw_inputs(**paths)
    mc.raise_on_any_failure(mc.run_validity_checks(raw))
    samples = mc.build_samples_table(raw)

    h4_rows = mc.hypothesis_h4(samples, max_total_permutations=1000)
    batch_row = next(r for r in h4_rows if r.scope.startswith("batch="))
    assert batch_row.passed is True
    assert batch_row.p_one_sided_observed_direction == pytest.approx(1 / 6)
    assert batch_row.direction == "control_earlier"


def test_h4_interleaved_design_not_flagged(tmp_path: Path) -> None:
    """Same 4 samples, but run order alternates control/raloxifene-d0 -> the
    pipeline must NOT flag this as (fully) confounded."""
    files = [
        "UWPRExp480_2021_0101_AZ_001_AZ201_AZ_complex.mzML",
        "UWPRExp480_2021_0101_AZ_002_AZ202_AZ_complex.mzML",
        "UWPRExp480_2021_0101_AZ_003_AZ203_AZ_complex.mzML",
        "UWPRExp480_2021_0101_AZ_004_AZ204_AZ_complex.mzML",
    ]
    conditions = ["control", "raloxifene-d0", "control", "raloxifene-d0"]
    mapping_text = ", ".join(f"{i + 1} ({f})" for i, f in enumerate(files))
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(files, conditions, strict=True)),
        mapping_text=mapping_text,
        protein_quants_header=(
            "Protein Groups\t"
            + "\t".join(f"Intensity_search_scan_file_id_{i + 1}" for i in range(4))
            + "\n"
        ),
        peptide_quants_header=(
            "Sequence\t"
            + "\t".join(f"Intensity_search_scan_file_id_{i + 1}" for i in range(4))
            + "\t"
            + "\t".join(f"Detection Type_search_scan_file_id_{i + 1}" for i in range(4))
            + "\n"
        ),
        protein_limelight_header=(
            "Protein(s)\t"
            + "\t".join(f"PSMs ({i + 1})" for i in range(4))
            + "\t"
            + "\t".join(f"NSAF ({i + 1})" for i in range(4))
            + "\n"
        ),
        peptide_limelight_header=(
            "Peptide Sequence\t"
            + "\t".join(f"PSMs ({i + 1})" for i in range(4))
            + "\t"
            + "\t".join(f"Quant ({i + 1})" for i in range(4))
            + "\n"
        ),
    )
    raw = mc.load_raw_inputs(**paths)
    mc.raise_on_any_failure(mc.run_validity_checks(raw))
    samples = mc.build_samples_table(raw)

    h4_rows = mc.hypothesis_h4(samples, max_total_permutations=1000)
    batch_row = next(r for r in h4_rows if r.scope.startswith("batch="))
    assert batch_row.passed is False


# ---------------------------------------------------------------------------
# Data-version stamp (role-keyed)
# ---------------------------------------------------------------------------


def test_compute_file_hashes_and_data_version_role_order_independent(
    tmp_path: Path,
) -> None:
    f1 = tmp_path / "a.txt"
    f1.write_text("hello", encoding="utf-8")
    f2 = tmp_path / "b.txt"
    f2.write_text("world", encoding="utf-8")

    hashes_1 = mc.compute_file_hashes([("role_a", f1), ("role_b", f2)])
    hashes_2 = mc.compute_file_hashes([("role_b", f2), ("role_a", f1)])

    assert hashes_1["role_a"]["path"] == str(f1)
    assert hashes_1["role_a"]["sha256"] == mc.sha256_of_file(f1)
    # Order of the (role, path) pairs must not affect the combined stamp.
    assert mc.compute_data_version(hashes_1) == mc.compute_data_version(hashes_2)
    assert mc.compute_data_version(hashes_1).startswith("sha256:")

    f1.write_text("HELLO", encoding="utf-8")
    hashes_3 = mc.compute_file_hashes([("role_a", f1), ("role_b", f2)])
    assert mc.compute_data_version(hashes_3) != mc.compute_data_version(hashes_1)


def test_compute_file_hashes_rejects_duplicate_roles(tmp_path: Path) -> None:
    f1 = tmp_path / "a.txt"
    f1.write_text("hello", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate"):
        mc.compute_file_hashes([("role_a", f1), ("role_a", f1)])


# ---------------------------------------------------------------------------
# --output-dir-under-data guard
# ---------------------------------------------------------------------------


def test_assert_output_dir_not_under_data_refuses_nested_or_equal_dir(
    tmp_path: Path,
) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    metadata_file = data_dir / "metadata.tsv"
    metadata_file.write_text("x", encoding="utf-8")

    with pytest.raises(ValueError, match="read-only"):
        mc._assert_output_dir_not_under_data(
            data_dir / "results", [("metadata_file", metadata_file)]
        )
    with pytest.raises(ValueError, match="read-only"):
        mc._assert_output_dir_not_under_data(
            data_dir, [("metadata_file", metadata_file)]
        )


def test_assert_output_dir_not_under_data_allows_sibling_dir(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    metadata_file = data_dir / "metadata.tsv"
    metadata_file.write_text("x", encoding="utf-8")
    mc._assert_output_dir_not_under_data(
        tmp_path / "results", [("metadata_file", metadata_file)]
    )  # must not raise


# ---------------------------------------------------------------------------
# Edge cases: empty / all-missing / single-sample / duplicate-id inputs
# ---------------------------------------------------------------------------


def test_metadata_with_only_header_fails_gracefully(tmp_path: Path) -> None:
    path = tmp_path / "metadata.tsv"
    path.write_text("Replicate\tcondition\n", encoding="utf-8")
    header, rows = mc.read_tsv_rows(path)
    assert header == ["Replicate", "condition"]
    assert rows == []


def test_single_sample_batch_cannot_run_h4_single_batch_test() -> None:
    with pytest.raises(ValueError, match="at least one member of each group"):
        mc.exact_permutation_test([1], [True])


def test_duplicate_ids_in_mapping_detected_as_check_failure(tmp_path: Path) -> None:
    paths = _write_fixture_files(
        tmp_path,
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=f"1 ({_GOOD_FILES[0]}), 1 ({_GOOD_FILES[1]})",
    )
    raw = mc.load_raw_inputs(**paths)
    checks = mc.run_validity_checks(raw)
    assert not next(
        c for c in checks if c.name == "mapping_exactly_n_unique_ids_and_files"
    ).passed
    with pytest.raises(ValueError):
        mc.raise_on_any_failure(checks)


# ---------------------------------------------------------------------------
# End-to-end main(): outputs written, byte-identical, failure marker
# ---------------------------------------------------------------------------


_OUTPUT_FILENAMES = [
    "samples.tsv",
    "validity_checks.tsv",
    "hypotheses.tsv",
    "associations.tsv",
    "crosstab_condition_batch.tsv",
    "crosstab_condition_run_half.tsv",
    "crosstab_condition_pair.tsv",
    "crosstab_batch_run_position.tsv",
    "run_layout.tsv",
    "table1.tsv",
    "table1.md",
    "join_key.tsv",
    "data_version.json",
    "summary.json",
]


def test_main_end_to_end_writes_expected_outputs_and_is_byte_identical(
    tmp_path: Path,
) -> None:
    paths = _write_fixture_files(
        tmp_path / "in",
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=_two_sample_mapping(),
    )
    out1 = tmp_path / "out1"
    out2 = tmp_path / "out2"

    mc.main(_argv_from_paths(paths, out1))
    mc.main(_argv_from_paths(paths, out2))

    for filename in _OUTPUT_FILENAMES:
        assert (out1 / filename).is_file(), f"{filename} missing from run 1"
        assert (out2 / filename).is_file(), f"{filename} missing from run 2"
        bytes1 = (out1 / filename).read_bytes()
        bytes2 = (out2 / filename).read_bytes()
        assert bytes1 == bytes2, f"{filename} differs between two identical runs"

    assert not (out1 / mc.FAILURE_MARKER_NAME).exists()
    assert not (out2 / mc.FAILURE_MARKER_NAME).exists()

    summary = json.loads((out1 / "summary.json").read_text(encoding="utf-8"))
    assert "output_dir" not in summary["params"]
    assert summary["validity_checks"]["all_passed"] is True


def test_main_writes_failure_marker_and_raises_and_no_other_outputs(
    tmp_path: Path,
) -> None:
    paths = _write_fixture_files(
        tmp_path / "in",
        metadata_rows=[(_GOOD_FILES[0], "BOGUS"), (_GOOD_FILES[1], "control")],
        mapping_text=_two_sample_mapping(),
    )
    out = tmp_path / "out"

    with pytest.raises(ValueError, match="validity check"):
        mc.main(_argv_from_paths(paths, out))

    assert (out / mc.FAILURE_MARKER_NAME).is_file()
    marker = json.loads((out / mc.FAILURE_MARKER_NAME).read_text(encoding="utf-8"))
    assert marker["failed"] is True
    assert "condition_in_allowed_set" in marker["error"]

    for filename in _OUTPUT_FILENAMES:
        assert not (out / filename).exists(), f"{filename} unexpectedly written"


def test_main_clears_stale_failure_marker_on_a_subsequent_success(
    tmp_path: Path,
) -> None:
    bad_paths = _write_fixture_files(
        tmp_path / "bad",
        metadata_rows=[(_GOOD_FILES[0], "BOGUS"), (_GOOD_FILES[1], "control")],
        mapping_text=_two_sample_mapping(),
    )
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        mc.main(_argv_from_paths(bad_paths, out))
    assert (out / mc.FAILURE_MARKER_NAME).is_file()

    good_paths = _write_fixture_files(
        tmp_path / "good",
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=_two_sample_mapping(),
    )
    mc.main(_argv_from_paths(good_paths, out))
    assert not (out / mc.FAILURE_MARKER_NAME).exists()
    assert (out / "samples.tsv").is_file()


def test_main_refuses_output_dir_under_data(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    paths = _write_fixture_files(
        data_dir,
        metadata_rows=list(zip(_GOOD_FILES, _GOOD_CONDITIONS, strict=True)),
        mapping_text=_two_sample_mapping(),
    )
    with pytest.raises(ValueError, match="read-only"):
        mc.main(_argv_from_paths(paths, data_dir / "results"))


# ---------------------------------------------------------------------------
# Property tests (hypothesis): exact permutation test invariants
# ---------------------------------------------------------------------------


@st.composite
def _batch_case(draw: st.DrawFn) -> tuple[list[int], list[bool]]:
    n = draw(st.integers(min_value=2, max_value=7))
    n_control = draw(st.integers(min_value=1, max_value=n - 1))
    values = draw(
        st.lists(
            st.integers(min_value=-1000, max_value=1000),
            min_size=n,
            max_size=n,
            unique=True,
        )
    )
    positions = list(range(n))
    # Choose which positions are control via a random subset of the right size.
    control_positions = draw(
        st.permutations(positions).map(lambda p: sorted(p[:n_control]))
    )
    is_control = [i in control_positions for i in positions]
    return values, is_control


@given(_batch_case())
def test_exact_permutation_p_values_are_valid_probabilities(
    case: tuple[list[int], list[bool]],
) -> None:
    values, is_control = case
    result = mc.exact_permutation_test(values, is_control)
    assert 0 < result.p_one_sided_observed_direction <= 1
    assert 0 < result.p_two_sided <= 1
    assert -1.0 - 1e-9 <= result.rank_biserial <= 1.0 + 1e-9
    assert 0 <= result.u_observed <= result.u_max


@given(_batch_case())
def test_exact_permutation_symmetric_under_label_swap(
    case: tuple[list[int], list[bool]],
) -> None:
    values, is_control = case
    result = mc.exact_permutation_test(values, is_control)
    swapped = [not c for c in is_control]
    result_swapped = mc.exact_permutation_test(values, swapped)

    assert result.p_one_sided_observed_direction == pytest.approx(
        result_swapped.p_one_sided_observed_direction
    )
    assert result.p_two_sided == pytest.approx(result_swapped.p_two_sided)
    assert result.rank_biserial == pytest.approx(-result_swapped.rank_biserial)
    assert result.u_observed + result_swapped.u_observed == result.u_max
