"""Stage-1 metadata characterization for the raloxifene/HLM LFQ experiment.

Parses the study's raw metadata (``metadata.tsv``, the scan-file<->search-id mapping
file, and the *header lines only* of the quant and Limelight table-dump files),
validates their structural integrity (fail loud — raise, do not warn), derives a tidy
per-sample table, tests a set of design-invariant hypotheses proposed at Stage 0/1
(batch/run-order/pairing structure that was *not* recorded by the scientist), computes
confounding statistics (Cramer's V, raw and bias-corrected) and cohort-characterization
tables, and writes everything as machine-readable TSV/JSON to ``results/metadata/`` for
downstream figure generation and for findings to cite.

Deliberately does **not** read the quantitative VALUES of any data file — only header
lines of the quant/Limelight-dump files (Stage 2 has not happened yet) — and does not
plot anything (a separate figure-generator consumes these tables).

Design-invariant hypotheses (H1-H5) are NOT validity checks: a failing hypothesis is a
*finding*, recorded with its evidence, not a reason to raise. Validity checks (file
structure, bijections, id-set equality) ARE fail-loud: any failure raises, and no output
is written for a run that fails them.

No stochastic step: every hypothesis test below is *exact* (full enumeration of the
combinatorial permutation space), so there is no RNG and no seed to record. Enumeration
size is bounded (``--max-stratified-permutations``) and the script raises rather than
silently truncating or taking an unbounded time if the bound would be exceeded.

Run:
    ./.venv/bin/python scripts/scratch/metadata_characterize.py \\
        --metadata-file data/metadata.tsv \\
        --mapping-file data/scan-file-search-id-mapping.txt \\
        --protein-quants-file data/protein-quants.tsv \\
        --peptide-quants-file data/peptide-quants.tsv \\
        --protein-limelight-file data/protein-limelight-table-dump.txt \\
        --peptide-limelight-file data/peptide-limelight-table-dump.txt \\
        --output-dir results/metadata
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import logging
import math
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd

__script_meta__: dict[str, object] = {
    "task": "metadata-characterize",
    "kind": "analysis",
    "provides": [],
    "uses": [],
    "seeded_from": {"template": "batch-correct-combat", "version": "0.2"},
    "description": (
        "Stage-1 metadata validity checks, design-invariant hypothesis tests "
        "(batch/run-order/pairing), confounding statistics (Cramer's V raw + "
        "bias-corrected), and cohort-characterization tables for the raloxifene/HLM "
        "LFQ metadata. Reads header lines only; makes no figures."
    ),
}

LOGGER = logging.getLogger("metadata_characterize")

# ---------------------------------------------------------------------------
# Constants: filename pattern, allowed values, column-label patterns
# ---------------------------------------------------------------------------

# UWPRExp480_<YYYY>_<MMDD>_AZ_<NNN>_AZ<nnn>_AZ_complex.mzML
FILENAME_RE = re.compile(
    r"^UWPRExp480_(?P<year>\d{4})_(?P<month>\d{2})(?P<day>\d{2})_AZ_"
    r"(?P<seq>\d{3})_(?P<sample_id>AZ\d{3})_AZ_complex\.mzML$"
)

ALLOWED_CONDITIONS: frozenset[str] = frozenset({"control", "raloxifene-d0"})

# A per-run label at the very end of a Limelight-dump column header, e.g.
# "PSMs (1_0506)" -> family "PSMs", label "1_0506". Anchored to the *end* of the
# string so a family name that itself contains parentheses (peptide dump's
# "Shared group [column content string: (SHARED GROUP)] (1_0)") is not confused.
LABEL_RE = re.compile(r"\s\(([^()]+)\)$")

MAPPING_ENTRY_RE = re.compile(r"^(\d+)\s*\(([^()]+)\)$")

INTENSITY_COL_RE = re.compile(r"^Intensity_search_scan_file_id_(\d+)$")
DETECTION_COL_RE = re.compile(r"^Detection Type_search_scan_file_id_(\d+)$")


# ---------------------------------------------------------------------------
# Parsing primitives
# ---------------------------------------------------------------------------


class FilenamePatternError(ValueError):
    """Raised when a Replicate filename does not match the expected pattern."""


class FilenameDateError(ValueError):
    """Raised when a Replicate filename encodes an invalid calendar date."""


@dataclass(frozen=True)
class ParsedFilename:
    """Fields derived from one ``Replicate`` mzML filename."""

    replicate: str
    year: int
    month: int
    day: int
    acq_date: date
    seq_number: int
    sample_id: str
    sample_number: int


def parse_filename(name: str) -> ParsedFilename:
    """Parse a Replicate filename against the scientist-confirmed pattern.

    Raises :class:`FilenamePatternError` if ``name`` does not match the pattern
    ``UWPRExp480_<YYYY>_<MMDD>_AZ_<NNN>_AZ<nnn>_AZ_complex.mzML``, or
    :class:`FilenameDateError` if the encoded year/month/day is not a valid calendar
    date. Both are `ValueError` subclasses so callers that only care about "fail loud"
    can catch `ValueError`.
    """
    match = FILENAME_RE.match(name)
    if match is None:
        raise FilenamePatternError(
            f"Replicate filename {name!r} does not match the expected pattern "
            f"'UWPRExp480_<YYYY>_<MMDD>_AZ_<NNN>_AZ<nnn>_AZ_complex.mzML'."
        )
    year = int(match.group("year"))
    month = int(match.group("month"))
    day = int(match.group("day"))
    try:
        acq_date = date(year, month, day)
    except ValueError as exc:
        raise FilenameDateError(
            f"Replicate filename {name!r} encodes an invalid calendar date "
            f"{year:04d}-{month:02d}-{day:02d}: {exc}"
        ) from exc
    seq_number = int(match.group("seq"))
    sample_id = match.group("sample_id")
    sample_number = int(sample_id[2:])
    return ParsedFilename(
        replicate=name,
        year=year,
        month=month,
        day=day,
        acq_date=acq_date,
        seq_number=seq_number,
        sample_id=sample_id,
        sample_number=sample_number,
    )


@dataclass(frozen=True)
class MappingEntry:
    """One ``<search_scan_file_id> (<mzML file name>)`` entry."""

    search_scan_file_id: int
    file_name: str


def parse_mapping_file(path: Path) -> list[MappingEntry]:
    """Parse the single-line ``id (file), id (file), ...`` mapping file.

    Raises ``ValueError`` if the file is empty or any comma-separated entry does not
    match ``<digits> (<file name>)``. Does not itself enforce uniqueness of ids or
    files (or a count of 8) — that is a validity *check*, evaluated by the caller so a
    structural violation is reported alongside the others rather than raising deep in
    the parser.
    """
    text = path.read_text(encoding="utf-8")
    stripped = text.strip()
    if not stripped:
        raise ValueError(f"Mapping file {path} is empty.")
    raw_entries = [entry.strip() for entry in stripped.split(",")]
    entries: list[MappingEntry] = []
    for raw in raw_entries:
        match = MAPPING_ENTRY_RE.match(raw)
        if match is None:
            raise ValueError(
                f"Mapping file {path} entry {raw!r} does not match the expected "
                f"'<search_scan_file_id> (<mzML file name>)' pattern."
            )
        entries.append(MappingEntry(int(match.group(1)), match.group(2)))
    return entries


def read_header(path: Path, *, sep: str = "\t") -> list[str]:
    """Return just the first (header) line of ``path``, split on ``sep``.

    Never reads beyond the first line — the contract for the quant/Limelight-dump
    files, whose quantitative VALUES this script must not read.
    """
    with path.open(encoding="utf-8") as handle:
        first_line = handle.readline()
    if not first_line:
        raise ValueError(f"{path} is empty (no header line).")
    return first_line.rstrip("\r\n").split(sep)


def read_tsv_rows(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read a small TSV file fully as raw strings: ``(header, data_rows)``.

    Used only for ``metadata.tsv`` (8 rows) — small enough, and control over exact
    cell contents (to check for blanks) matters more here than streaming.

    A wholly blank *line* (``csv.reader`` yields ``[]``) is only tolerated when it
    trails the last real row (e.g. a trailing newline at EOF) — it is silently
    dropped in that position only. A blank line *between* real rows raises: silently
    skipping an interior blank line would drop a data row without any record of it,
    exactly the "silent row drop" the correctness charter forbids.
    """
    with path.open(newline="", encoding="utf-8") as handle:
        raw_rows = list(csv.reader(handle, delimiter="\t"))
    nonblank_indices = [i for i, row in enumerate(raw_rows) if row]
    if not nonblank_indices:
        raise ValueError(f"{path} is empty.")
    last_nonblank = nonblank_indices[-1]
    interior_blanks = [
        i + 1 for i, row in enumerate(raw_rows[:last_nonblank]) if not row
    ]
    if interior_blanks:
        raise ValueError(
            f"{path} has blank line(s) at line number(s) {interior_blanks} "
            f"in between data rows; only a *trailing* blank line is tolerated "
            f"(e.g. a trailing newline at EOF). An interior blank line would "
            f"otherwise silently drop a data row."
        )
    rows = raw_rows[: last_nonblank + 1]
    header, *data_rows = rows
    return header, data_rows


def extract_ids_raw(
    header: Sequence[str], pattern: re.Pattern[str]
) -> list[tuple[str, int]]:
    """Return ``(raw_captured_text, int_value)`` for each matching header column.

    Returns a **list**, in header order, with duplicates preserved (never a
    ``set``) — a header with the same id in two columns, or with one column
    fewer/more than expected, must stay visible to the caller rather than being
    silently collapsed away by set deduplication.
    """
    results: list[tuple[str, int]] = []
    for col in header:
        match = pattern.match(col)
        if match is not None:
            text = match.group(1)
            results.append((text, int(text)))
    return results


def _check_quant_header_ids(
    name: str,
    header: Sequence[str],
    pattern: re.Pattern[str],
    mapping_ids: Sequence[int],
    expected_n: int,
    header_label: str,
) -> CheckResult:
    """Validate one quant-file header's per-run id columns against the mapping ids.

    Checks, as **lists** (not sets, so a duplicated or missing column is visible
    rather than silently absorbed by set equality):
      * exactly ``expected_n`` matching columns;
      * no duplicate id across those columns;
      * every captured id's raw text is already canonical (``text == str(int(text))``
        — rejects a leading-zero variant like ``_018454`` that would silently parse
        to the right integer but is not the id as actually written);
      * the resulting ids, as a sorted list, equal the mapping ids' sorted list.
    """
    raw = extract_ids_raw(header, pattern)
    texts = [text for text, _ in raw]
    ids = [value for _, value in raw]
    problems: list[str] = []

    if len(ids) != expected_n:
        problems.append(f"found {len(ids)} matching column(s), expected {expected_n}")

    dup_ids = _duplicates(str(i) for i in ids)
    if dup_ids:
        problems.append(f"duplicate id(s) across columns: {dup_ids}")

    non_canonical = [t for t in texts if t != str(int(t))]
    if non_canonical:
        problems.append(f"non-canonical id text (e.g. leading zeros): {non_canonical}")

    if sorted(ids) != sorted(mapping_ids):
        problems.append(
            f"{header_label} ids {sorted(ids)} != mapping ids {sorted(mapping_ids)}"
        )

    passed = not problems
    detail = "OK" if passed else "; ".join(problems)
    return _check(name, passed, detail)


def extract_family_labels(header: Sequence[str]) -> dict[str, list[str]]:
    """Group Limelight-dump per-run columns by measure family -> list of run labels.

    A column with no trailing ``" (label)"`` (e.g. an id column like ``"Protein(s)"``)
    is skipped. Family names are used only to group; they are not otherwise validated.
    """
    families: dict[str, list[str]] = {}
    for col in header:
        match = LABEL_RE.search(col)
        if match is None:
            continue
        label = match.group(1)
        family = col[: match.start()].strip()
        families.setdefault(family, []).append(label)
    return families


def sha256_of_file(path: Path, *, chunk_size: int = 1 << 20) -> str:
    """Stream-hash ``path`` (raw bytes; never parsed/interpreted) with sha256."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compute_file_hashes(
    data_files: Sequence[tuple[str, Path]],
) -> dict[str, dict[str, str]]:
    """Hash each ``(role, path)`` pair; returns ``{role: {"path": ..., "sha256": ..}}``.

    Keyed by ROLE (``metadata_file``, ``mapping_file``, ...), not by filename, so
    two data files that happen to share a basename cannot collide/overwrite each
    other's hash. Raises if a role appears more than once (the combined stamp below
    assumes exactly one hash per role).
    """
    roles = [role for role, _ in data_files]
    dup_roles = _duplicates(roles)
    if dup_roles:
        raise ValueError(f"Duplicate data-file role(s): {dup_roles}")
    return {
        role: {"path": str(path), "sha256": sha256_of_file(path)}
        for role, path in data_files
    }


def compute_data_version(file_hashes: dict[str, dict[str, str]]) -> str:
    """Combine per-role hashes into one ``sha256:<hex>`` data-version stamp.

    The combined stamp is computed over sorted ``"<role>:<sha256>"`` lines (role,
    not path/filename, so the stamp is stable regardless of where the input files
    happen to live on disk).
    """
    lines = sorted(f"{role}:{info['sha256']}\n" for role, info in file_hashes.items())
    combined = hashlib.sha256("".join(lines).encode("utf-8")).hexdigest()
    return f"sha256:{combined}"


def _duplicates(items: Iterable[str]) -> list[str]:
    """Return the values that appear more than once in ``items`` (sorted, unique)."""
    seen: set[str] = set()
    dupes: set[str] = set()
    for item in items:
        if item in seen:
            dupes.add(item)
        seen.add(item)
    return sorted(dupes)


# ---------------------------------------------------------------------------
# Raw inputs + validity checks
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RawInputs:
    """Everything read from disk, unvalidated (header lines / mapping / metadata)."""

    metadata_header: list[str]
    metadata_rows: list[list[str]]
    mapping_entries: list[MappingEntry]
    protein_quants_header: list[str]
    peptide_quants_header: list[str]
    protein_limelight_header: list[str]
    peptide_limelight_header: list[str]


def load_raw_inputs(
    *,
    metadata_file: Path,
    mapping_file: Path,
    protein_quants_file: Path,
    peptide_quants_file: Path,
    protein_limelight_file: Path,
    peptide_limelight_file: Path,
) -> RawInputs:
    """Read all raw inputs: full ``metadata.tsv``/mapping file, header lines only
    for the quant and Limelight-dump files."""
    metadata_header, metadata_rows = read_tsv_rows(metadata_file)
    mapping_entries = parse_mapping_file(mapping_file)
    return RawInputs(
        metadata_header=metadata_header,
        metadata_rows=metadata_rows,
        mapping_entries=mapping_entries,
        protein_quants_header=read_header(protein_quants_file),
        peptide_quants_header=read_header(peptide_quants_file),
        protein_limelight_header=read_header(protein_limelight_file),
        peptide_limelight_header=read_header(peptide_limelight_file),
    )


@dataclass(frozen=True)
class CheckResult:
    """One validity-check outcome: name, pass/fail, human-readable evidence."""

    name: str
    passed: bool
    detail: str


def _check(name: str, passed: bool, detail: str) -> CheckResult:
    return CheckResult(name=name, passed=passed, detail=detail)


def run_validity_checks(raw: RawInputs) -> list[CheckResult]:
    """Evaluate every Stage-1 structural validity check; never raises itself.

    Every check is evaluated (even after an earlier one fails) so a failing run
    reports the full picture in one shot; the caller decides whether to raise.
    """
    results: list[CheckResult] = []

    # -- metadata.tsv structural checks --------------------------------------------
    expected_header = ["Replicate", "condition"]
    results.append(
        _check(
            "metadata_header_exact",
            raw.metadata_header == expected_header,
            f"expected {expected_header}, got {raw.metadata_header}",
        )
    )

    blanks = [
        (r_idx + 2, c_idx)
        for r_idx, row in enumerate(raw.metadata_rows)
        for c_idx, cell in enumerate(row)
        if cell.strip() == ""
    ]
    results.append(
        _check(
            "metadata_no_blank_cells",
            len(blanks) == 0,
            "no blank cells"
            if not blanks
            else f"blank cell(s) at (row, col): {blanks}",
        )
    )

    n_cols_ok = all(len(row) == len(expected_header) for row in raw.metadata_rows)
    results.append(
        _check(
            "metadata_row_width_consistent",
            n_cols_ok,
            f"all {len(raw.metadata_rows)} rows have {len(expected_header)} cells"
            if n_cols_ok
            else "one or more rows do not have exactly 2 cells",
        )
    )

    replicates = [row[0] for row in raw.metadata_rows if len(row) > 0]
    dup_replicates = _duplicates(replicates)
    results.append(
        _check(
            "replicate_unique",
            len(dup_replicates) == 0,
            "all Replicate values unique"
            if not dup_replicates
            else f"duplicate Replicate value(s): {dup_replicates}",
        )
    )

    parse_errors: list[str] = []
    date_errors: list[str] = []
    parsed_by_replicate: dict[str, ParsedFilename] = {}
    for replicate in replicates:
        try:
            parsed_by_replicate[replicate] = parse_filename(replicate)
        except FilenamePatternError as exc:
            parse_errors.append(str(exc))
        except FilenameDateError as exc:
            # Pattern matched (right digit groups) but the encoded year/month/day is
            # not a real calendar date (e.g. month 13, Feb 30) — a distinct failure
            # mode from a structural pattern mismatch, reported by its own check
            # below rather than left to propagate and abort the whole function.
            date_errors.append(str(exc))
    results.append(
        _check(
            "replicate_matches_filename_regex",
            len(parse_errors) == 0,
            "all Replicate values match the filename pattern"
            if not parse_errors
            else "; ".join(parse_errors),
        )
    )

    results.append(
        _check(
            "replicate_dates_valid_calendar_dates",
            len(date_errors) == 0,
            "all encoded dates are valid calendar dates"
            if not date_errors
            else "; ".join(date_errors),
        )
    )

    conditions = [row[1] for row in raw.metadata_rows if len(row) > 1]
    bad_conditions = sorted({c for c in conditions if c not in ALLOWED_CONDITIONS})
    results.append(
        _check(
            "condition_in_allowed_set",
            len(bad_conditions) == 0,
            f"all conditions in {sorted(ALLOWED_CONDITIONS)}"
            if not bad_conditions
            else f"unexpected condition value(s): {bad_conditions}",
        )
    )

    if not parse_errors and not date_errors:
        sample_ids = [p.sample_id for p in parsed_by_replicate.values()]
        dup_sample_ids = _duplicates(sample_ids)
        results.append(
            _check(
                "sample_id_unique",
                len(dup_sample_ids) == 0,
                "all sample_id values unique"
                if not dup_sample_ids
                else f"duplicate sample_id value(s): {dup_sample_ids}",
            )
        )

        seq_by_date: dict[date, list[int]] = {}
        for parsed in parsed_by_replicate.values():
            seq_by_date.setdefault(parsed.acq_date, []).append(parsed.seq_number)
        seq_dupes = {
            d: _duplicates(str(s) for s in seqs)
            for d, seqs in seq_by_date.items()
            if len(seqs) != len(set(seqs))
        }
        results.append(
            _check(
                "seq_number_unique_within_date",
                len(seq_dupes) == 0,
                "seq_number unique within each acquisition date"
                if not seq_dupes
                else f"duplicate seq_number(s) within date: {seq_dupes}",
            )
        )
    else:
        skip_detail = "skipped: filename regex/date check(s) failed first"
        results.append(_check("sample_id_unique", False, skip_detail))
        results.append(_check("seq_number_unique_within_date", False, skip_detail))

    # -- mapping-file checks ----------------------------------------------------------
    # Expected count is derived from the metadata itself (== number of Replicate rows,
    # 8 for this study), not hardcoded, so the same check is meaningfully testable on
    # smaller synthetic fixtures while still enforcing "exactly 8" on the real data.
    expected_n = len(replicates)
    mapping_ids = [e.search_scan_file_id for e in raw.mapping_entries]
    mapping_files = [e.file_name for e in raw.mapping_entries]
    dup_ids = _duplicates(str(i) for i in mapping_ids)
    dup_files = _duplicates(mapping_files)
    results.append(
        _check(
            "mapping_exactly_n_unique_ids_and_files",
            len(mapping_ids) == expected_n
            and len(mapping_files) == expected_n
            and not dup_ids
            and not dup_files,
            f"expected {expected_n} (from metadata row count); got "
            f"{len(mapping_ids)} ids ({len(set(mapping_ids))} unique), "
            f"{len(mapping_files)} files ({len(set(mapping_files))} unique)",
        )
    )

    replicate_set = set(replicates)
    mapping_file_set = set(mapping_files)
    missing_in_mapping = sorted(replicate_set - mapping_file_set)
    extra_in_mapping = sorted(mapping_file_set - replicate_set)
    results.append(
        _check(
            "mapping_bijection_with_metadata_replicates",
            not missing_in_mapping and not extra_in_mapping,
            "mapping file names == metadata Replicate values"
            if not (missing_in_mapping or extra_in_mapping)
            else (
                f"in metadata but not mapping: {missing_in_mapping}; "
                f"in mapping but not metadata: {extra_in_mapping}"
            ),
        )
    )

    results.append(
        _check_quant_header_ids(
            "protein_quants_intensity_ids_match_mapping",
            raw.protein_quants_header,
            INTENSITY_COL_RE,
            mapping_ids,
            expected_n,
            "protein-quants Intensity",
        )
    )

    results.append(
        _check_quant_header_ids(
            "peptide_quants_intensity_ids_match_mapping",
            raw.peptide_quants_header,
            INTENSITY_COL_RE,
            mapping_ids,
            expected_n,
            "peptide-quants Intensity",
        )
    )

    results.append(
        _check_quant_header_ids(
            "peptide_quants_detection_type_ids_match_mapping",
            raw.peptide_quants_header,
            DETECTION_COL_RE,
            mapping_ids,
            expected_n,
            "peptide-quants 'Detection Type'",
        )
    )

    # -- Limelight dump checks --------------------------------------------------------
    protein_families = extract_family_labels(raw.protein_limelight_header)
    peptide_families = extract_family_labels(raw.peptide_limelight_header)

    for dump_name, families in (
        ("protein_limelight", protein_families),
        ("peptide_limelight", peptide_families),
    ):
        label_sets = {family: set(labels) for family, labels in families.items()}
        counts_ok = all(
            len(labels) == expected_n and len(set(labels)) == expected_n
            for labels in families.values()
        )
        distinct_sets = {frozenset(s) for s in label_sets.values()}
        consistent = len(distinct_sets) <= 1
        results.append(
            _check(
                f"{dump_name}_n_distinct_labels_per_family_consistent",
                counts_ok and consistent and len(families) > 0,
                f"expected {expected_n} distinct labels per family (from metadata "
                f"row count); families={list(families)}; per-family label counts="
                f"{ {f: len(v) for f, v in families.items()} }; "
                f"consistent across families={consistent}",
            )
        )

    protein_label_set = next(
        iter({frozenset(s) for s in protein_families.values()}), frozenset()
    )
    peptide_label_set = next(
        iter({frozenset(s) for s in peptide_families.values()}), frozenset()
    )
    results.append(
        _check(
            "limelight_label_sets_match_across_files",
            protein_label_set == peptide_label_set
            and len(protein_label_set) == expected_n,
            f"protein labels={sorted(protein_label_set)}; "
            f"peptide labels={sorted(peptide_label_set)}",
        )
    )

    return results


def raise_on_any_failure(checks: list[CheckResult]) -> None:
    """Fail loud: raise ``ValueError`` listing every failed check, or return quietly."""
    failed = [c for c in checks if not c.passed]
    if failed:
        detail = "\n".join(f"  - {c.name}: {c.detail}" for c in failed)
        raise ValueError(
            f"{len(failed)} of {len(checks)} Stage-1 metadata validity check(s) "
            f"FAILED:\n{detail}"
        )


# ---------------------------------------------------------------------------
# Phase B: derived tidy samples table (only reached once all checks pass)
# ---------------------------------------------------------------------------


def build_samples_table(raw: RawInputs) -> pd.DataFrame:
    """Build the tidy per-sample table. Assumes ``run_validity_checks`` all passed."""
    mapping_by_file = {e.file_name: e.search_scan_file_id for e in raw.mapping_entries}

    records: list[dict[str, Any]] = []
    for replicate, condition in ((row[0], row[1]) for row in raw.metadata_rows):
        parsed = parse_filename(replicate)
        pair_low = (
            parsed.sample_number
            if parsed.sample_number % 2 == 1
            else parsed.sample_number - 1
        )
        pair_high = pair_low + 1
        search_scan_file_id = mapping_by_file[replicate]
        records.append(
            {
                "file": replicate,
                "condition": condition,
                "acq_date": parsed.acq_date.isoformat(),
                "batch": f"B{parsed.acq_date.isoformat()}",
                "seq_number": parsed.seq_number,
                "sample_id": parsed.sample_id,
                "sample_number": parsed.sample_number,
                "candidate_pair": f"P{pair_low}_{pair_high}",
                "sample_role": "experimental",
                "search_scan_file_id": search_scan_file_id,
                "intensity_column": (
                    f"Intensity_search_scan_file_id_{search_scan_file_id}"
                ),
                "detection_type_column": (
                    f"Detection Type_search_scan_file_id_{search_scan_file_id}"
                ),
            }
        )
    samples = pd.DataFrame.from_records(records)
    samples = samples.sort_values(["batch", "seq_number"], kind="stable").reset_index(
        drop=True
    )

    samples["run_position_within_batch"] = (
        samples.groupby("batch")["seq_number"].rank(method="first").astype(int)
    )
    batch_size = samples.groupby("batch")["batch"].transform("size")
    # ceil(size / 2): for an ODD-sized batch, the extra (middle) sample goes to
    # "early", not "late" — a deliberate, documented tie-break (also called out in
    # the condition x run_half association's output note), not an accident of
    # integer division.
    half_size = -(-batch_size // 2)
    is_early = samples["run_position_within_batch"] <= half_size
    samples["run_half"] = is_early.map({True: "early", False: "late"})

    column_order = [
        "file",
        "condition",
        "acq_date",
        "batch",
        "seq_number",
        "run_position_within_batch",
        "run_half",
        "sample_id",
        "sample_number",
        "candidate_pair",
        "sample_role",
        "search_scan_file_id",
        "intensity_column",
        "detection_type_column",
    ]
    return samples[column_order]


# ---------------------------------------------------------------------------
# Exact permutation test machinery (H4)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PermutationResult:
    """Exact Mann-Whitney-style permutation test result (integer U, no ties).

    ``values`` (``seq_numbers`` in the caller) are treated as an order proxy where
    smaller = earlier; ``direction`` and ``rank_biserial`` describe whether the
    "control" group tends to sit earlier or later in that order.

    ``p_one_sided_observed_direction`` is the one-sided p-value *in whichever
    direction was actually observed* — that direction was **not pre-specified**
    before looking at the data, so this number is reported for reference only and
    is not a pre-specification-honest p-value (picking the direction post hoc after
    seeing which way the data leaned inflates it, the classic one-sided-after-the-
    fact problem). ``p_two_sided`` — which does not depend on which direction was
    observed — is the primary, pre-specification-honest value; treat it as such.
    """

    n_total: int
    n_control: int
    n_other: int
    u_observed: int
    u_max: int
    rank_biserial: float
    direction: str
    n_permutations: int
    p_one_sided_observed_direction: float
    p_two_sided: float


def _u_statistic(values: Sequence[int], control_positions: Sequence[int]) -> int:
    """Count pairs (control value, other value) with control value < other value."""
    control_set = set(control_positions)
    control_vals = [values[i] for i in control_positions]
    other_vals = [values[i] for i in range(len(values)) if i not in control_set]
    return sum(1 for cv in control_vals for ov in other_vals if cv < ov)


def _batch_u_distribution(
    seq_numbers: Sequence[int], is_control: Sequence[bool]
) -> tuple[int, list[int], int]:
    """Return ``(u_observed, all_u_values_over_combinations, u_max)`` for one batch."""
    n = len(seq_numbers)
    if len(seq_numbers) != len(is_control):
        raise ValueError("seq_numbers and is_control must be the same length.")
    if len(set(seq_numbers)) != n:
        raise ValueError("Batch seq_numbers must be distinct (no ties supported).")
    n_control = int(sum(is_control))
    n_other = n - n_control
    u_max = n_control * n_other
    observed_positions = [i for i, c in enumerate(is_control) if c]
    u_obs = _u_statistic(seq_numbers, observed_positions)
    if n_control == 0 or n_other == 0:
        return u_obs, [u_obs], u_max
    combos = itertools.combinations(range(n), n_control)
    us = [_u_statistic(seq_numbers, combo) for combo in combos]
    return u_obs, us, u_max


def _permutation_result_from_u(
    *,
    n_total: int,
    n_control: int,
    n_other: int,
    u_obs: int,
    u_max: int,
    u_values: Sequence[int],
) -> PermutationResult:
    """Shared one-/two-sided exact p-value math given an enumerated U distribution."""
    n_perm = len(u_values)
    center2 = u_max  # 2 * center, since center = u_max / 2
    obs_dist2 = abs(2 * u_obs - center2)
    n_two = sum(1 for u in u_values if abs(2 * u - center2) >= obs_dist2)
    p_two = n_two / n_perm

    if 2 * u_obs >= center2:
        n_one = sum(1 for u in u_values if u >= u_obs)
    else:
        n_one = sum(1 for u in u_values if u <= u_obs)
    p_one = n_one / n_perm

    rank_biserial = (2.0 * u_obs / u_max) - 1.0 if u_max > 0 else 0.0
    if u_max == 0:
        direction = "undetermined (u_max=0: one group is empty)"
    elif rank_biserial > 0:
        direction = "control_earlier"
    elif rank_biserial < 0:
        direction = "control_later"
    else:
        direction = "no_difference"

    return PermutationResult(
        n_total=n_total,
        n_control=n_control,
        n_other=n_other,
        u_observed=u_obs,
        u_max=u_max,
        rank_biserial=rank_biserial,
        direction=direction,
        n_permutations=n_perm,
        p_one_sided_observed_direction=p_one,
        p_two_sided=p_two,
    )


def exact_permutation_test(
    seq_numbers: Sequence[int], is_control: Sequence[bool]
) -> PermutationResult:
    """Exact one-batch permutation test of seq_number order vs. condition label.

    Enumerates every way of choosing which positions are "control" among the fixed,
    distinct ``seq_numbers`` (``C(n, n_control)`` equally likely assignments) and
    computes the Mann-Whitney U (count of control-before-other pairs). Requires at
    least one member of each group.
    """
    n_control = int(sum(is_control))
    n_other = len(is_control) - n_control
    if n_control == 0 or n_other == 0:
        raise ValueError(
            "exact_permutation_test requires at least one member of each group; got "
            f"n_control={n_control}, n_other={n_other}."
        )
    u_obs, u_values, u_max = _batch_u_distribution(seq_numbers, is_control)
    return _permutation_result_from_u(
        n_total=len(seq_numbers),
        n_control=n_control,
        n_other=n_other,
        u_obs=u_obs,
        u_max=u_max,
        u_values=u_values,
    )


def exact_stratified_permutation_test(
    batches: Sequence[tuple[Sequence[int], Sequence[bool]]],
    *,
    max_total_permutations: int = 2_000_000,
) -> PermutationResult:
    """Exact stratified permutation test: permute condition independently within
    each batch, combined statistic = sum of per-batch U. Raises rather than silently
    truncating if the joint enumeration would exceed ``max_total_permutations``."""
    if not batches:
        raise ValueError(
            "exact_stratified_permutation_test requires at least one batch."
        )
    per_batch = [_batch_u_distribution(seq, ctrl) for seq, ctrl in batches]

    total_perms = 1
    for _, us, _ in per_batch:
        total_perms *= len(us)
    if total_perms > max_total_permutations:
        raise ValueError(
            f"Stratified enumeration would require {total_perms} joint permutations "
            f"(> max_total_permutations={max_total_permutations}); refusing rather "
            f"than silently taking a long time."
        )

    u_obs_total = sum(u for u, _, _ in per_batch)
    u_max_total = sum(m for _, _, m in per_batch)
    n_total = sum(len(seq) for seq, _ in batches)
    n_control_total = sum(int(sum(ctrl)) for _, ctrl in batches)
    n_other_total = n_total - n_control_total

    combined_u_values = [
        sum(combo) for combo in itertools.product(*(us for _, us, _ in per_batch))
    ]

    return _permutation_result_from_u(
        n_total=n_total,
        n_control=n_control_total,
        n_other=n_other_total,
        u_obs=u_obs_total,
        u_max=u_max_total,
        u_values=combined_u_values,
    )


# ---------------------------------------------------------------------------
# Association statistics (Cramer's V raw + bias-corrected)
# ---------------------------------------------------------------------------


def cramers_v_uncorrected(crosstab: pd.DataFrame) -> float:
    """Textbook (uncorrected) Cramer's V; ``0.0`` for a degenerate table.

    ``expected`` is always strictly positive here: every row/column of a
    ``pd.crosstab`` corresponds to a category with at least one observation, so
    both its row-sum and column-sum are > 0, and ``expected = row_sum * col_sum
    / n`` can never be zero. No zero-division guard is needed (matches
    :func:`cramers_v_bias_corrected`, which relies on the same fact).
    """
    observed = crosstab.to_numpy(dtype=float)
    n = float(observed.sum())
    n_rows, n_cols = observed.shape
    if n <= 0 or min(n_rows, n_cols) < 2:
        return 0.0
    row_sums = observed.sum(axis=1, keepdims=True)
    col_sums = observed.sum(axis=0, keepdims=True)
    expected = row_sums @ col_sums / n
    chi2 = float(((observed - expected) ** 2 / expected).sum())
    phi2 = chi2 / n
    denom = min(n_rows, n_cols) - 1
    if denom <= 0:
        return 0.0
    return math.sqrt(phi2 / denom)


def cramers_v_bias_corrected(crosstab: pd.DataFrame) -> float:
    """Bias-corrected Cramer's V (Bergsma 2013).

    Adapted from ``lib/common/batch_correct.py::_cramers_v`` (seeded_from lineage,
    see module ``__script_meta__``). Reimplemented here rather than imported so this
    Stage-1 script does not depend on that module's ``pycombat`` runtime dependency,
    which this project does not otherwise need. Returns ``0.0`` for a degenerate table
    (one row/col, or n <= 1).
    """
    observed = crosstab.to_numpy(dtype=float)
    n = float(observed.sum())
    n_rows, n_cols = observed.shape
    if n <= 1 or min(n_rows, n_cols) < 2:
        return 0.0
    row_sums = observed.sum(axis=1, keepdims=True)
    col_sums = observed.sum(axis=0, keepdims=True)
    expected = row_sums @ col_sums / n
    chi2 = float(((observed - expected) ** 2 / expected).sum())
    phi2 = chi2 / n
    phi2_corr = max(0.0, phi2 - (n_rows - 1) * (n_cols - 1) / (n - 1))
    rows_corr = n_rows - (n_rows - 1) ** 2 / (n - 1)
    cols_corr = n_cols - (n_cols - 1) ** 2 / (n - 1)
    denom = min(rows_corr - 1, cols_corr - 1)
    if denom <= 0:
        return 0.0
    return math.sqrt(phi2_corr / denom)


@dataclass(frozen=True)
class AssociationResult:
    """Cramer's V (raw + bias-corrected) for one pair of categorical columns."""

    variable_a: str
    variable_b: str
    n: int
    cramers_v_raw: float
    cramers_v_bias_corrected: float
    crosstab: pd.DataFrame


def compute_association(
    samples: pd.DataFrame, col_a: str, col_b: str
) -> AssociationResult:
    """Build the ``col_a`` x ``col_b`` crosstab and both Cramer's V variants."""
    crosstab = pd.crosstab(samples[col_a], samples[col_b])
    return AssociationResult(
        variable_a=col_a,
        variable_b=col_b,
        n=int(crosstab.to_numpy().sum()),
        cramers_v_raw=cramers_v_uncorrected(crosstab),
        cramers_v_bias_corrected=cramers_v_bias_corrected(crosstab),
        crosstab=crosstab,
    )


# ---------------------------------------------------------------------------
# Hypotheses H1-H5
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HypothesisRow:
    """One row of the ``hypotheses.tsv`` output.

    ``p_one_sided_observed_direction`` / ``direction``: the direction (e.g.
    "control earlier") was NOT pre-specified before looking at the data, so the
    one-sided p-value is reported for reference only, not as a pre-specification-
    honest result. ``p_two_sided`` is the primary, pre-specification-honest value.
    """

    hypothesis_id: str
    scope: str
    description: str
    passed: bool | None  # None = inconclusive / no information
    statistic_name: str
    statistic_value: str
    direction: str | None
    p_one_sided_observed_direction: float | None
    p_two_sided: float | None
    note: str


def hypothesis_h1(samples: pd.DataFrame) -> HypothesisRow:
    counts = samples["condition"].value_counts()
    n_control = int(counts.get("control", 0))
    n_ralox = int(counts.get("raloxifene-d0", 0))
    passed = n_control == 4 and n_ralox == 4
    return HypothesisRow(
        hypothesis_id="H1",
        scope="overall",
        description="Balanced arms: 4 control, 4 raloxifene-d0.",
        passed=passed,
        statistic_name="n_control,n_raloxifene-d0",
        statistic_value=f"{n_control},{n_ralox}",
        direction=None,
        p_one_sided_observed_direction=None,
        p_two_sided=None,
        note="See crosstab_condition_batch.tsv for the per-batch breakdown.",
    )


def hypothesis_h2(samples: pd.DataFrame) -> HypothesisRow:
    crosstab = pd.crosstab(samples["batch"], samples["condition"])
    both_present = bool((crosstab > 0).all(axis=1).all())
    return HypothesisRow(
        hypothesis_id="H2",
        scope="overall",
        description=(
            "Every batch contains both conditions (condition not aliased with batch)."
        ),
        passed=both_present,
        statistic_name="per_batch_condition_counts",
        statistic_value=str(crosstab.to_dict()),
        direction=None,
        p_one_sided_observed_direction=None,
        p_two_sided=None,
        note="See crosstab_condition_batch.tsv for full counts.",
    )


def hypothesis_h3(samples: pd.DataFrame) -> HypothesisRow:
    pair_crosstab = pd.crosstab(samples["candidate_pair"], samples["condition"])
    one_each = (
        bool((pair_crosstab.to_numpy() == 1).all()) and pair_crosstab.shape[1] == 2
    )
    is_odd = samples["sample_number"] % 2 == 1
    is_control = samples["condition"] == "control"
    parity_holds = bool((is_odd == is_control).all())
    passed = one_each and parity_holds
    return HypothesisRow(
        hypothesis_id="H3",
        scope="overall",
        description=(
            "Each candidate pair (consecutive sample-ID pair) contains exactly one "
            "control and one raloxifene-d0; parity rule (control <-> odd sample "
            "number) holds for all 8."
        ),
        passed=passed,
        statistic_name="one_control_one_ralox_per_pair,parity_rule_holds",
        statistic_value=f"{one_each},{parity_holds}",
        direction=None,
        p_one_sided_observed_direction=None,
        p_two_sided=None,
        note="See crosstab_condition_pair.tsv for the pair-level counts.",
    )


def hypothesis_h4(
    samples: pd.DataFrame, *, max_total_permutations: int
) -> list[HypothesisRow]:
    rows: list[HypothesisRow] = []
    batches_for_stratified: list[tuple[list[int], list[bool]]] = []
    for batch, group in samples.groupby("batch", sort=True):
        seq_numbers = group["seq_number"].tolist()
        is_control = (group["condition"] == "control").tolist()
        n_control = sum(is_control)
        n_other = len(is_control) - n_control
        batches_for_stratified.append((seq_numbers, is_control))
        if n_control >= 2 and n_other >= 2:
            result = exact_permutation_test(seq_numbers, is_control)
            fully_separated = result.u_observed in (0, result.u_max)
            rows.append(
                HypothesisRow(
                    hypothesis_id="H4",
                    scope=f"batch={batch}",
                    description=(
                        "Within-batch run order (seq_number) vs. condition: exact "
                        "permutation test of separation, this batch only."
                    ),
                    passed=fully_separated,
                    statistic_name="U_observed/U_max,rank_biserial",
                    statistic_value=f"{result.u_observed}/{result.u_max},{result.rank_biserial:.4f}",
                    direction=result.direction,
                    p_one_sided_observed_direction=result.p_one_sided_observed_direction,
                    p_two_sided=result.p_two_sided,
                    note=(
                        f"n={result.n_total} ({result.n_control} control, "
                        f"{result.n_other} raloxifene-d0), "
                        f"{result.n_permutations} exact permutations enumerated. "
                        f"direction={result.direction!r} was NOT pre-specified "
                        "before inspecting the data — p_one_sided_observed_direction "
                        "is for reference only; p_two_sided is the primary, "
                        "pre-specification-honest value."
                    ),
                )
            )
        else:
            rows.append(
                HypothesisRow(
                    hypothesis_id="H4",
                    scope=f"batch={batch}",
                    description=(
                        "Within-batch run order (seq_number) vs. condition: this "
                        "batch has fewer than 2 samples of one condition."
                    ),
                    passed=None,
                    statistic_name="n_control,n_other",
                    statistic_value=f"{n_control},{n_other}",
                    direction=None,
                    p_one_sided_observed_direction=None,
                    p_two_sided=None,
                    note=(
                        "No information on its own: with e.g. 1 vs 1 there are only "
                        "2 possible label assignments, so any observed order is "
                        "trivially 'fully separated' and uninformative in isolation."
                    ),
                )
            )

    stratified = exact_stratified_permutation_test(
        batches_for_stratified, max_total_permutations=max_total_permutations
    )
    rows.append(
        HypothesisRow(
            hypothesis_id="H4",
            scope="stratified_across_batches",
            description=(
                "Stratified exact permutation test: condition permuted independently "
                "within each batch; combined statistic = sum of per-batch U."
            ),
            passed=stratified.u_observed in (0, stratified.u_max),
            statistic_name="U_observed/U_max,rank_biserial",
            statistic_value=(
                f"{stratified.u_observed}/{stratified.u_max},{stratified.rank_biserial:.4f}"
            ),
            direction=stratified.direction,
            p_one_sided_observed_direction=stratified.p_one_sided_observed_direction,
            p_two_sided=stratified.p_two_sided,
            note=(
                f"n={stratified.n_total} across {len(batches_for_stratified)} "
                f"batch(es), {stratified.n_permutations} joint exact permutations "
                "enumerated. This is the primary evidence for run-order confounding "
                "(the single 1-vs-1 2022 batch carries no information alone but "
                "contributes to the joint enumeration). "
                f"direction={stratified.direction!r} was NOT pre-specified before "
                "inspecting the data — p_one_sided_observed_direction is for "
                "reference only; p_two_sided is the primary, "
                "pre-specification-honest value."
            ),
        )
    )
    return rows


def hypothesis_h5(samples: pd.DataFrame) -> list[HypothesisRow]:
    rows: list[HypothesisRow] = []
    for batch, group in samples.groupby("batch", sort=True):
        ordered_by_seq = group.sort_values("seq_number")["sample_number"].tolist()
        ordered_by_sample = sorted(group["sample_number"].tolist())
        matches = ordered_by_seq == ordered_by_sample
        rows.append(
            HypothesisRow(
                hypothesis_id="H5",
                scope=f"batch={batch}",
                description=(
                    "Sample-number order equals run (seq_number) order within batch "
                    "(were samples run in ID order?)."
                ),
                passed=matches,
                statistic_name="run_order_sample_numbers,id_order_sample_numbers",
                statistic_value=f"{ordered_by_seq},{ordered_by_sample}",
                direction=None,
                p_one_sided_observed_direction=None,
                p_two_sided=None,
                note=(
                    "Exact match required (n is too small here for a meaningful "
                    "correlation coefficient)."
                ),
            )
        )
    return rows


# ---------------------------------------------------------------------------
# Cohort characterization: crosstabs, run layout, Table 1, join key
# ---------------------------------------------------------------------------


def build_run_layout(samples: pd.DataFrame) -> pd.DataFrame:
    """Rows = batch, columns = run position (1..max), cell = 'condition/sample_id'."""
    max_position = int(samples["run_position_within_batch"].max())
    batches = sorted(samples["batch"].unique())
    layout = pd.DataFrame(
        "",
        index=pd.Index(batches, name="batch"),
        columns=[str(p) for p in range(1, max_position + 1)],
    )
    for row in samples.itertuples():
        layout.loc[row.batch, str(row.run_position_within_batch)] = (
            f"{row.condition}/{row.sample_id}"
        )
    return layout.reset_index()


def build_table1(samples: pd.DataFrame) -> pd.DataFrame:
    """Publication-ready Table 1, one row per condition."""
    rows = []
    for condition, group in samples.groupby("condition", sort=True):
        batch_counts = group["batch"].value_counts().sort_index()
        batch_counts_str = "; ".join(f"{b}={c}" for b, c in batch_counts.items())
        rows.append(
            {
                "condition": condition,
                "n": len(group),
                "batch_counts": batch_counts_str,
                "seq_number_min": int(group["seq_number"].min()),
                "seq_number_max": int(group["seq_number"].max()),
                "median_run_position": float(
                    group["run_position_within_batch"].median()
                ),
                "sample_ids": ",".join(sorted(group["sample_id"])),
                "note": (
                    "No covariates (e.g. donor, sex, age) are recorded in the metadata."
                ),
            }
        )
    return pd.DataFrame(rows)


def render_table1_markdown(table1: pd.DataFrame) -> str:
    """Render ``table1`` as a Markdown table (no external dependency needed)."""
    columns = list(table1.columns)
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    lines = [header, separator]
    for _, row in table1.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in columns) + " |")
    return "\n".join(lines) + "\n"


def build_join_key(samples: pd.DataFrame) -> pd.DataFrame:
    """File <-> search_scan_file_id <-> quant column names <-> sample_id/condition."""
    join_key = samples[
        [
            "file",
            "search_scan_file_id",
            "intensity_column",
            "detection_type_column",
            "sample_id",
            "condition",
        ]
    ].copy()
    join_key["limelight_dump_label"] = "UNRESOLVED"
    return join_key


# ---------------------------------------------------------------------------
# CLI / orchestration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Args:
    metadata_file: Path
    mapping_file: Path
    protein_quants_file: Path
    peptide_quants_file: Path
    protein_limelight_file: Path
    peptide_limelight_file: Path
    output_dir: Path
    max_stratified_permutations: int
    log_level: str


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-file", type=Path, default=Path("data/metadata.tsv"))
    parser.add_argument(
        "--mapping-file",
        type=Path,
        default=Path("data/scan-file-search-id-mapping.txt"),
    )
    parser.add_argument(
        "--protein-quants-file", type=Path, default=Path("data/protein-quants.tsv")
    )
    parser.add_argument(
        "--peptide-quants-file", type=Path, default=Path("data/peptide-quants.tsv")
    )
    parser.add_argument(
        "--protein-limelight-file",
        type=Path,
        default=Path("data/protein-limelight-table-dump.txt"),
    )
    parser.add_argument(
        "--peptide-limelight-file",
        type=Path,
        default=Path("data/peptide-limelight-table-dump.txt"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/metadata"))
    parser.add_argument(
        "--max-stratified-permutations",
        type=int,
        default=2_000_000,
        help="Refuse the stratified H4 test rather than enumerate more than this many "
        "joint permutations.",
    )
    parser.add_argument("--log-level", default="INFO")
    ns = parser.parse_args(argv)
    return Args(
        metadata_file=ns.metadata_file,
        mapping_file=ns.mapping_file,
        protein_quants_file=ns.protein_quants_file,
        peptide_quants_file=ns.peptide_quants_file,
        protein_limelight_file=ns.protein_limelight_file,
        peptide_limelight_file=ns.peptide_limelight_file,
        output_dir=ns.output_dir,
        max_stratified_permutations=ns.max_stratified_permutations,
        log_level=ns.log_level,
    )


def _association_row(result: AssociationResult) -> dict[str, object]:
    note = (
        "n is small (<=8); the bias-corrected V is unstable in this regime "
        "(can be driven to 0 even under apparent association) and the raw V is "
        "upward-biased. Treat the exact-permutation result (H4) as the primary "
        "evidence for run-order confounding; use both V's only as descriptive "
        "context."
    )
    if result.variable_b == "run_half":
        note += (
            " run_half splits each batch at ceil(batch_size / 2): for an "
            "odd-sized batch the extra (middle) sample is assigned to the "
            "'early' half, not 'late'."
        )
    return {
        "variable_a": result.variable_a,
        "variable_b": result.variable_b,
        "n": result.n,
        "cramers_v_raw": round(result.cramers_v_raw, 6),
        "cramers_v_bias_corrected": round(result.cramers_v_bias_corrected, 6),
        "note": note,
    }


def _hypothesis_row_to_dict(row: HypothesisRow) -> dict[str, object]:
    return {
        "hypothesis_id": row.hypothesis_id,
        "scope": row.scope,
        "description": row.description,
        "passed": "" if row.passed is None else row.passed,
        "statistic_name": row.statistic_name,
        "statistic_value": str(row.statistic_value),
        "direction": "" if row.direction is None else row.direction,
        "p_one_sided_observed_direction": (
            ""
            if row.p_one_sided_observed_direction is None
            else row.p_one_sided_observed_direction
        ),
        "p_two_sided": "" if row.p_two_sided is None else row.p_two_sided,
        "note": row.note,
    }


def _assert_output_dir_not_under_data(
    output_dir: Path, data_files: Sequence[tuple[str, Path]]
) -> None:
    """Refuse an ``output_dir`` that resolves inside a read-only data directory.

    ``data/`` is read-only by project convention (enforced elsewhere by a hook);
    this is a fail-loud, script-level guard against pointing ``--output-dir`` at
    it (or a subdirectory of it) so a misconfigured run fails clearly here rather
    than mysteriously at the write step.
    """
    resolved_output = output_dir.resolve()
    for role, path in data_files:
        data_dir = path.resolve().parent
        if resolved_output == data_dir or data_dir in resolved_output.parents:
            raise ValueError(
                f"--output-dir {output_dir} resolves under the same directory as "
                f"input data file {role!r} ({data_dir}), which is read-only by "
                "project convention. Choose an output directory outside the data "
                "directory (e.g. results/metadata)."
            )


FAILURE_MARKER_NAME = "FAILED.json"


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    logging.basicConfig(
        level=args.log_level, format="%(asctime)s %(levelname)s %(message)s"
    )
    LOGGER.info("metadata_characterize starting")
    LOGGER.info("params: %s", args)

    data_files: dict[str, Path] = {
        "metadata_file": args.metadata_file,
        "mapping_file": args.mapping_file,
        "protein_quants_file": args.protein_quants_file,
        "peptide_quants_file": args.peptide_quants_file,
        "protein_limelight_file": args.protein_limelight_file,
        "peptide_limelight_file": args.peptide_limelight_file,
    }
    for name, path in data_files.items():
        if not path.is_file():
            raise FileNotFoundError(f"{name} not found: {path}")

    data_file_pairs = list(data_files.items())
    _assert_output_dir_not_under_data(args.output_dir, data_file_pairs)

    # Create the output directory and clear any stale failure marker from a
    # previous FAILED attempt up front, so that (a) the marker is written into an
    # existing directory if this attempt also fails, and (b) a successful run
    # never leaves a stale marker behind claiming the (fresh) outputs are bad.
    args.output_dir.mkdir(parents=True, exist_ok=True)
    marker_path = args.output_dir / FAILURE_MARKER_NAME
    marker_path.unlink(missing_ok=True)

    try:
        LOGGER.info(
            "hashing raw data files for the data-version stamp (bytes only, no parsing)"
        )
        file_hashes = compute_file_hashes(data_file_pairs)
        data_version = compute_data_version(file_hashes)
        LOGGER.info("data_version = %s", data_version)

        LOGGER.info(
            "reading raw inputs (header lines only for quant/Limelight-dump files)"
        )
        raw = load_raw_inputs(
            metadata_file=args.metadata_file,
            mapping_file=args.mapping_file,
            protein_quants_file=args.protein_quants_file,
            peptide_quants_file=args.peptide_quants_file,
            protein_limelight_file=args.protein_limelight_file,
            peptide_limelight_file=args.peptide_limelight_file,
        )
        LOGGER.info(
            "metadata.tsv: %d data row(s); mapping file: %d entries",
            len(raw.metadata_rows),
            len(raw.mapping_entries),
        )

        checks = run_validity_checks(raw)
        n_passed = sum(c.passed for c in checks)
        LOGGER.info("validity checks: %d/%d passed", n_passed, len(checks))
        for c in checks:
            (LOGGER.info if c.passed else LOGGER.error)(
                "  [%s] %s: %s", "PASS" if c.passed else "FAIL", c.name, c.detail
            )
        raise_on_any_failure(checks)

        samples = build_samples_table(raw)
        LOGGER.info("built samples table: shape=%s", samples.shape)

        h1 = hypothesis_h1(samples)
        h2 = hypothesis_h2(samples)
        h3 = hypothesis_h3(samples)
        h4_rows = hypothesis_h4(
            samples, max_total_permutations=args.max_stratified_permutations
        )
        h5_rows = hypothesis_h5(samples)
        hypotheses = [h1, h2, h3, *h4_rows, *h5_rows]
        for h in hypotheses:
            LOGGER.info(
                "[%s/%s] passed=%s %s",
                h.hypothesis_id,
                h.scope,
                h.passed,
                h.statistic_value,
            )

        assoc_condition_batch = compute_association(samples, "condition", "batch")
        assoc_condition_run_half = compute_association(samples, "condition", "run_half")
        assoc_condition_pair = compute_association(
            samples, "condition", "candidate_pair"
        )
        associations = [
            assoc_condition_batch,
            assoc_condition_run_half,
            assoc_condition_pair,
        ]
        for a in associations:
            LOGGER.info(
                "association %s x %s: V_raw=%.4f V_corrected=%.4f",
                a.variable_a,
                a.variable_b,
                a.cramers_v_raw,
                a.cramers_v_bias_corrected,
            )

        batch_run_position = pd.crosstab(
            samples["batch"], samples["run_position_within_batch"]
        )
        run_layout = build_run_layout(samples)
        table1 = build_table1(samples)
        table1_md = render_table1_markdown(table1)
        join_key = build_join_key(samples)

        samples.to_csv(args.output_dir / "samples.tsv", sep="\t", index=False)

        validity_df = pd.DataFrame(
            [{"name": c.name, "passed": c.passed, "detail": c.detail} for c in checks]
        )
        validity_df.to_csv(
            args.output_dir / "validity_checks.tsv", sep="\t", index=False
        )

        hypotheses_df = pd.DataFrame([_hypothesis_row_to_dict(h) for h in hypotheses])
        hypotheses_df.to_csv(args.output_dir / "hypotheses.tsv", sep="\t", index=False)

        associations_df = pd.DataFrame([_association_row(a) for a in associations])
        associations_df.to_csv(
            args.output_dir / "associations.tsv", sep="\t", index=False
        )

        assoc_condition_batch.crosstab.to_csv(
            args.output_dir / "crosstab_condition_batch.tsv", sep="\t"
        )
        assoc_condition_run_half.crosstab.to_csv(
            args.output_dir / "crosstab_condition_run_half.tsv", sep="\t"
        )
        assoc_condition_pair.crosstab.to_csv(
            args.output_dir / "crosstab_condition_pair.tsv", sep="\t"
        )
        batch_run_position.to_csv(
            args.output_dir / "crosstab_batch_run_position.tsv", sep="\t"
        )

        run_layout.to_csv(args.output_dir / "run_layout.tsv", sep="\t", index=False)
        table1.to_csv(args.output_dir / "table1.tsv", sep="\t", index=False)
        (args.output_dir / "table1.md").write_text(table1_md, encoding="utf-8")
        join_key.to_csv(args.output_dir / "join_key.tsv", sep="\t", index=False)

        data_version_payload = {
            "file_hashes": file_hashes,
            "data_version": data_version,
        }
        (args.output_dir / "data_version.json").write_text(
            json.dumps(data_version_payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )

        summary = {
            "params": {
                "metadata_file": str(args.metadata_file),
                "mapping_file": str(args.mapping_file),
                "protein_quants_file": str(args.protein_quants_file),
                "peptide_quants_file": str(args.peptide_quants_file),
                "protein_limelight_file": str(args.protein_limelight_file),
                "peptide_limelight_file": str(args.peptide_limelight_file),
                # output_dir is deliberately NOT recorded here: it is the one
                # parameter that varies purely with *where* you chose to write
                # results, never with *what* was computed, and including it would
                # make summary.json differ byte-for-byte between two runs of the
                # identical inputs written to two different --output-dir values.
                "max_stratified_permutations": args.max_stratified_permutations,
            },
            "data_version": data_version,
            "file_hashes": file_hashes,
            "n_samples": len(samples),
            "validity_checks": {
                "n_passed": int(n_passed),
                "n_total": len(checks),
                "all_passed": True,
            },
            "hypotheses": [_hypothesis_row_to_dict(h) for h in hypotheses],
            "associations": [_association_row(a) for a in associations],
            "table1": table1.to_dict(orient="records"),
            "join_key_n_rows": len(join_key),
        }
        (args.output_dir / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True, default=str),
            encoding="utf-8",
        )
    except Exception as exc:
        marker_path.write_text(
            json.dumps(
                {
                    "failed": True,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        LOGGER.error("metadata_characterize FAILED: %s", exc)
        raise

    LOGGER.info("wrote outputs to %s", args.output_dir)
    LOGGER.info("metadata_characterize done")


if __name__ == "__main__":
    main()
