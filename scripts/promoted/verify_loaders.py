"""Stage-3 integrity gate: verify the protein/peptide/Limelight loaders on real data.

Obligation (B) of the two-part integrity gate (conventions/correctness.md) — obligation
(A), testing the loaders themselves, lives in ``tests/test_loaders_*.py``. This script
loads the real project data through the three loaders and checks the result against the
raw files by an INDEPENDENT parsing path (Python's ``csv`` module, never pandas/the
loader's own code) so a loader bug cannot be common-mode with its own verification.

Checks performed (see ``run_checks`` for the full list): raw row/column counts vs. an
independent header parse; random seeded spot-cell reconciliation (>= 200 cells per
matrix) against the raw text; orientation; dtypes; value ranges; identifier round-trip;
missing-token counts (0 / literal "NaN"); the peptide Detection-Type invariant,
independently re-tallied; contaminant counts (dump-refined 33/449, scientist decision)
and complete non-contaminant feature counts; exact sample<->metadata pairing
(condition/batch/sample_role breakdowns), INDEPENDENTLY RE-DERIVED row-by-row from
``data/metadata.tsv`` + ``data/scan-file-search-id-mapping.txt`` (never from
``samples.tsv``); the Limelight label map (reproduced by the loader itself, which
cross-checks ``results/stage2/limelight_label_map.tsv``); a FULL (not sampled) PSMs/
NSAF reconciliation of every one of the 4,344x8 cells against an independent csv-module
parse of the dump; a two-way derivation of per-sample total/median intensity (Dataset
vs. the independent raw parse); and that the combined data_version over the same six
roles as ``results/metadata/data_version.json`` reproduces that file's stamp exactly.

Fails loud: any check failure raises after all checks have run (so a single run reports
every failure, not just the first), writes a FAILED marker into ``--output-dir``, and
exits non-zero.

Run:
    ./.venv/bin/python scripts/promoted/verify_loaders.py \\
        --metadata-file data/metadata.tsv \\
        --mapping-file data/scan-file-search-id-mapping.txt \\
        --protein-quants-file data/protein-quants.tsv \\
        --peptide-quants-file data/peptide-quants.tsv \\
        --protein-limelight-file data/protein-limelight-table-dump.txt \\
        --peptide-limelight-file data/peptide-limelight-table-dump.txt \\
        --samples-file results/metadata/samples.tsv \\
        --cross-check-file results/stage2/limelight_label_map.tsv \\
        --output-dir results/stage3
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import random
import re
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from common.design import FAILURE_MARKER_NAME  # noqa: E402
from common.hashing import compute_data_version, compute_file_hashes  # noqa: E402
from loaders.limelight_loader import load_limelight_protein_counts  # noqa: E402
from loaders.peptide_loader import (  # noqa: E402
    POSITIVE_DETECTION_TYPES,
    ZERO_DETECTION_TYPES,
    load_peptide_dataset,
)
from loaders.protein_loader import load_protein_dataset  # noqa: E402
from loaders.samples import read_samples  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "verify-loaders",
    "kind": "analysis",
    "provides": [],
    "uses": [
        "common.hashing",
        "common.design",
        "loaders.protein_loader",
        "loaders.peptide_loader",
        "loaders.limelight_loader",
        "loaders.samples",
    ],
    "seeded_from": None,
    "description": (
        "Stage-3 integrity-gate verification: loads the real protein/peptide/"
        "Limelight data through the project loaders and cross-checks counts, spot "
        "cells, dtypes, value ranges, identifier round-trip, missing-token counts, "
        "the Detection-Type invariant, contaminant/complete-feature counts, exact "
        "sample<->metadata pairing, the Limelight label map, and a two-way total/"
        "median derivation -- all against an INDEPENDENT raw csv-module parse, "
        "never the loaders' own pandas code. Fail-loud (FAILED marker + nonzero "
        "exit) on any mismatch."
    ),
}

LOGGER = logging.getLogger("verify_loaders")

EXPECTED_N_PROTEINS = 4344
EXPECTED_N_PEPTIDES = 33358
EXPECTED_N_SAMPLES = 8
EXPECTED_PROTEIN_ZEROS = 13084
EXPECTED_PROTEIN_NAN_TOKENS = 14
EXPECTED_PEPTIDE_ZEROS = 89181
# Scientist decision (Stage-3 code review): a no-accession protein-quants row is a
# TRUE contaminant only if its Limelight dump group has no real-accession sibling
# (33 of the 52 no-accession candidates); the other 19 are kept
# (contaminant_grouped_with_real=True). Peptides are excluded only via the 33.
EXPECTED_PROTEIN_CONTAMINANTS = 33
EXPECTED_PEPTIDE_CONTAMINANTS = 449
EXPECTED_PROTEIN_COMPLETE_NONCONTAMINANT = 1801
EXPECTED_PEPTIDE_COMPLETE_NONCONTAMINANT = 10229
EXPECTED_DETECTION_TYPE_COUNTS = {
    "MSMS": 152578,
    "MBR": 25105,
    "NotDetected": 75363,
    "MSMSIdentifiedButNotQuantified": 10135,
    "MSMSAmbiguousPeakfinding": 3683,
}
EXPECTED_N_COMMA_JOINED_GROUPS = 27
# The combined data_version stamp results/metadata/data_version.json records over
# these same six roles (state/METADATA.md); this script must reproduce it exactly.
EXPECTED_DATA_VERSION = (
    "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
)

MAPPING_ENTRY_RE = re.compile(r"^(\d+)\s*\(([^()]+)\)$")


@dataclass(frozen=True)
class CheckResult:
    name: str
    passed: bool
    detail: str


def _check(name: str, passed: bool, detail: str) -> CheckResult:
    return CheckResult(name=name, passed=passed, detail=detail)


# --------------------------------------------------------------------------- #
# Independent raw parsing (csv module, never pandas) -- the whole point of (B).
# --------------------------------------------------------------------------- #


def read_raw_table(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read a TSV wholly via ``csv.reader``: ``(header, data_rows)``. No pandas."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path} is empty.")
    header, *data_rows = rows
    return header, data_rows


def _col_index(header: Sequence[str], name: str) -> int:
    return header.index(name)


def _raw_cell_to_expected(raw: str) -> float:
    """The loader's documented missing policy, applied to one raw cell."""
    if raw == "NaN":
        return float("nan")
    value = float(raw)
    return float("nan") if value == 0.0 else value


# --------------------------------------------------------------------------- #
# Independent pairing re-derivation (item D): data/metadata.tsv +
# data/scan-file-search-id-mapping.txt via csv/re -- NEVER results/metadata/
# samples.tsv, so a bug shared between metadata_characterize.py and this script
# cannot hide behind agreement with its own derived artifact.
# --------------------------------------------------------------------------- #

REPLICATE_RE = re.compile(
    r"^UWPRExp480_(?P<year>\d{4})_(?P<month>\d{2})(?P<day>\d{2})_AZ_"
    r"(?P<seq>\d{3})_(?P<sample_id>AZ\d{3})_AZ_complex\.mzML$"
)
# Substrings that would mark a run as a QC/pool/reference/blank sample rather than
# an experimental one (state/METADATA.md "Experimental vs. control samples" -- none
# are expected in this study, checked directly rather than assumed).
QC_MARKER_SUBSTRINGS = ("pool", "qc", "reference", "blank", "standard")


def _parse_mapping_file(path: Path) -> dict[str, str]:
    """Independent parse of the ``id (file), id (file), ...`` mapping file."""
    text = path.read_text(encoding="utf-8").strip()
    out: dict[str, str] = {}
    for raw_entry in text.split(","):
        raw_entry = raw_entry.strip()
        match = MAPPING_ENTRY_RE.match(raw_entry)
        if match is None:
            raise ValueError(
                f"{path}: entry {raw_entry!r} does not match '<id> (<file>)'."
            )
        out[match.group(2)] = match.group(1)
    return out


def derive_pairing_independently(
    metadata_file: Path, mapping_file: Path
) -> pd.DataFrame:
    """Re-derive sample_id/search_scan_file_id/condition/batch from the two raw
    source files only (csv module + regex; never pandas, never samples.tsv)."""
    with metadata_file.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"{metadata_file} is empty.")
    header, *data_rows = rows
    if header != ["Replicate", "condition"]:
        raise ValueError(f"{metadata_file}: unexpected header {header}.")
    mapping = _parse_mapping_file(mapping_file)

    records: list[dict[str, str]] = []
    for row in data_rows:
        if not row:
            continue
        replicate, condition = row
        match = REPLICATE_RE.match(replicate)
        if match is None:
            raise ValueError(f"{metadata_file}: Replicate {replicate!r} unparseable.")
        if replicate not in mapping:
            raise ValueError(f"{replicate!r} is missing from {mapping_file}.")
        batch = f"B{match.group('year')}-{match.group('month')}-{match.group('day')}"
        records.append(
            {
                "file": replicate,
                "condition": condition,
                "batch": batch,
                "sample_id": match.group("sample_id"),
                "search_scan_file_id": mapping[replicate],
            }
        )
    if len(records) != len(data_rows):
        raise ValueError(
            f"{metadata_file}: {len(data_rows)} data row(s) but only "
            f"{len(records)} parsed (blank line(s)?)."
        )
    return pd.DataFrame.from_records(records)


# --------------------------------------------------------------------------- #
# Individual check groups
# --------------------------------------------------------------------------- #


def check_raw_counts(
    protein_header: list[str],
    protein_rows: list[list[str]],
    peptide_header: list[str],
    peptide_rows: list[list[str]],
) -> list[CheckResult]:
    checks = []
    checks.append(
        _check(
            "protein_raw_row_count",
            len(protein_rows) == EXPECTED_N_PROTEINS,
            f"expected {EXPECTED_N_PROTEINS}, got {len(protein_rows)}",
        )
    )
    checks.append(
        _check(
            "peptide_raw_row_count",
            len(peptide_rows) == EXPECTED_N_PEPTIDES,
            f"expected {EXPECTED_N_PEPTIDES}, got {len(peptide_rows)}",
        )
    )
    n_protein_intensity = sum(
        1 for c in protein_header if c.startswith("Intensity_search_scan_file_id_")
    )
    n_peptide_intensity = sum(
        1 for c in peptide_header if c.startswith("Intensity_search_scan_file_id_")
    )
    n_peptide_detection = sum(
        1 for c in peptide_header if c.startswith("Detection Type_search_scan_file_id_")
    )
    checks.append(
        _check(
            "protein_header_n_intensity_columns",
            n_protein_intensity == EXPECTED_N_SAMPLES,
            f"expected {EXPECTED_N_SAMPLES}, got {n_protein_intensity}",
        )
    )
    checks.append(
        _check(
            "peptide_header_n_intensity_and_detection_columns",
            n_peptide_intensity == EXPECTED_N_SAMPLES
            and n_peptide_detection == EXPECTED_N_SAMPLES,
            f"expected {EXPECTED_N_SAMPLES}/{EXPECTED_N_SAMPLES}, got "
            f"{n_peptide_intensity}/{n_peptide_detection}",
        )
    )
    return checks


def check_orientation_and_dtype(
    name: str,
    ds_abundances: np.ndarray,
    ds_feature_names: np.ndarray,
    n_expected_features: int,
) -> list[CheckResult]:
    checks = []
    checks.append(
        _check(
            f"{name}_orientation_rows_are_samples",
            ds_abundances.shape == (EXPECTED_N_SAMPLES, n_expected_features),
            f"abundances.shape={ds_abundances.shape}, expected "
            f"({EXPECTED_N_SAMPLES}, {n_expected_features})",
        )
    )
    checks.append(
        _check(
            f"{name}_abundances_dtype_float",
            np.issubdtype(ds_abundances.dtype, np.floating),
            f"dtype={ds_abundances.dtype}",
        )
    )
    checks.append(
        _check(
            f"{name}_feature_names_dtype_str",
            ds_feature_names.dtype.kind in ("U", "O"),
            f"dtype={ds_feature_names.dtype}",
        )
    )
    finite = ds_abundances[np.isfinite(ds_abundances)]
    all_positive = bool((finite > 0).all()) if finite.size else True
    checks.append(
        _check(
            f"{name}_finite_values_all_positive",
            all_positive,
            f"n_finite={finite.size}, n_nonpositive={int((finite <= 0).sum())}",
        )
    )
    checks.append(
        _check(
            f"{name}_finite_values_all_finite",
            bool(np.isfinite(finite).all()),
            "no inf among finite-flagged values"
            if finite.size
            else "n/a (no finite values)",
        )
    )
    return checks


def check_identifier_roundtrip(
    name: str,
    ds_feature_names: np.ndarray,
    raw_header: list[str],
    raw_rows: list[list[str]],
    id_col_name: str,
) -> CheckResult:
    idx = _col_index(raw_header, id_col_name)
    raw_ids = [row[idx] for row in raw_rows]
    matches = list(ds_feature_names) == raw_ids
    if matches:
        detail = "feature_names == raw column order/text exactly"
    else:
        first_diff = next(
            (
                i
                for i, (a, b) in enumerate(zip(ds_feature_names, raw_ids, strict=True))
                if a != b
            ),
            None,
        )
        detail = f"mismatch at first differing index {first_diff}"
    return _check(f"{name}_identifier_roundtrip", matches, detail)


def check_missing_counts_protein(
    was_nan_token: np.ndarray, was_zero_token: np.ndarray
) -> list[CheckResult]:
    n_nan = int(was_nan_token.sum())
    n_zero = int(was_zero_token.sum())
    return [
        _check(
            "protein_nan_token_count",
            n_nan == EXPECTED_PROTEIN_NAN_TOKENS,
            f"expected {EXPECTED_PROTEIN_NAN_TOKENS}, got {n_nan}",
        ),
        _check(
            "protein_zero_token_count",
            n_zero == EXPECTED_PROTEIN_ZEROS,
            f"expected {EXPECTED_PROTEIN_ZEROS}, got {n_zero}",
        ),
    ]


def check_missing_counts_peptide(abundances: np.ndarray) -> CheckResult:
    n_missing = int(np.isnan(abundances).sum())
    return _check(
        "peptide_zero_count",
        n_missing == EXPECTED_PEPTIDE_ZEROS,
        f"expected {EXPECTED_PEPTIDE_ZEROS}, got {n_missing}",
    )


def check_detection_type_tally(
    peptide_header: list[str], peptide_rows: list[list[str]]
) -> CheckResult:
    """Independently re-tally Detection Type counts from the raw csv parse."""
    detection_cols = [
        i
        for i, c in enumerate(peptide_header)
        if c.startswith("Detection Type_search_scan_file_id_")
    ]
    counts: dict[str, int] = dict.fromkeys(EXPECTED_DETECTION_TYPE_COUNTS, 0)
    unknown = 0
    for row in peptide_rows:
        for i in detection_cols:
            value = row[i]
            if value in counts:
                counts[value] += 1
            else:
                unknown += 1
    # Cross-check the hardcoded expectation against the loader's own type taxonomy
    # (loaders.peptide_loader.POSITIVE_DETECTION_TYPES / ZERO_DETECTION_TYPES), so a
    # future change to that taxonomy cannot silently drift from this script's
    # constant without a check noticing.
    taxonomy_agrees = set(EXPECTED_DETECTION_TYPE_COUNTS) == (
        POSITIVE_DETECTION_TYPES | ZERO_DETECTION_TYPES
    )
    passed = (
        counts == EXPECTED_DETECTION_TYPE_COUNTS and unknown == 0 and taxonomy_agrees
    )
    detail = (
        f"counts={counts}, unknown={unknown}, "
        f"expected={EXPECTED_DETECTION_TYPE_COUNTS}, "
        f"taxonomy_agrees_with_loader={taxonomy_agrees}"
    )
    return _check("peptide_detection_type_tally", passed, detail)


def check_contaminant_and_complete_counts(
    name: str,
    feature_metadata: pd.DataFrame,
    abundances: np.ndarray,
    expected_contaminants: int,
    expected_complete: int,
) -> list[CheckResult]:
    is_contaminant = feature_metadata["is_contaminant"].to_numpy(dtype=bool)
    n_contaminant = int(is_contaminant.sum())
    non_contaminant_abundances = abundances[:, ~is_contaminant]
    complete = np.isfinite(non_contaminant_abundances).all(axis=0)
    n_complete = int(complete.sum())
    return [
        _check(
            f"{name}_contaminant_count",
            n_contaminant == expected_contaminants,
            f"expected {expected_contaminants}, got {n_contaminant}",
        ),
        _check(
            f"{name}_complete_noncontaminant_count",
            n_complete == expected_complete,
            f"expected {expected_complete}, got {n_complete}",
        ),
    ]


def check_pairing(
    independent: pd.DataFrame,
    samples: pd.DataFrame,
    protein_header: list[str],
    peptide_header: list[str],
) -> list[CheckResult]:
    """Pairing checks driven by ``independent`` (re-derived from the two raw source
    files, never ``samples.tsv``); ``samples`` (the loaders' actual input) is only
    used for the row-by-row cross-check, so a bug shared between
    ``metadata_characterize.py`` and this script cannot hide as agreement."""
    checks = []
    checks.append(
        _check(
            "pairing_n_samples",
            len(independent) == EXPECTED_N_SAMPLES,
            f"expected {EXPECTED_N_SAMPLES}, got {len(independent)}",
        )
    )
    condition_counts = independent["condition"].value_counts().to_dict()
    checks.append(
        _check(
            "pairing_condition_4_4",
            condition_counts.get("control") == 4
            and condition_counts.get("raloxifene-d0") == 4,
            f"condition_counts={condition_counts}",
        )
    )
    batch_counts = independent["batch"].value_counts().to_dict()
    checks.append(
        _check(
            "pairing_batch_6_2",
            batch_counts.get("B2021-05-06") == 6
            and batch_counts.get("B2022-03-18") == 2,
            f"batch_counts={batch_counts}",
        )
    )
    qc_hits = [
        f
        for f in independent["file"]
        if any(marker in f.lower() for marker in QC_MARKER_SUBSTRINGS)
    ]
    checks.append(
        _check(
            "pairing_experimental_8_control_0",
            len(qc_hits) == 0,
            "no QC/pool/reference/blank marker in any of the 8 Replicate filenames"
            if not qc_hits
            else f"QC-marker filename(s): {qc_hits}",
        )
    )

    merged = independent.merge(
        samples[["sample_id", "search_scan_file_id", "condition", "batch"]],
        on="sample_id",
        how="outer",
        suffixes=("_independent", "_samples_tsv"),
        indicator=True,
    )
    only_one_side = merged[merged["_merge"] != "both"]
    disagreements = merged[
        (merged["_merge"] == "both")
        & (
            (
                merged["search_scan_file_id_independent"]
                != merged["search_scan_file_id_samples_tsv"]
            )
            | (merged["condition_independent"] != merged["condition_samples_tsv"])
            | (merged["batch_independent"] != merged["batch_samples_tsv"])
        )
    ]
    checks.append(
        _check(
            "pairing_matches_samples_tsv_row_by_row",
            only_one_side.empty and disagreements.empty,
            "independent re-derivation == samples.tsv for every sample_id "
            "(search_scan_file_id, condition, batch)"
            if only_one_side.empty and disagreements.empty
            else f"only-one-side: {only_one_side['sample_id'].tolist()}; "
            f"disagreements: {disagreements['sample_id'].tolist()}",
        )
    )

    protein_ids = {
        m.group(1)
        for c in protein_header
        if (m := re.match(r"^Intensity_search_scan_file_id_(\d+)$", c))
    }
    peptide_ids = {
        m.group(1)
        for c in peptide_header
        if (m := re.match(r"^Intensity_search_scan_file_id_(\d+)$", c))
    }
    missing_in_protein = sorted(set(independent["search_scan_file_id"]) - protein_ids)
    missing_in_peptide = sorted(set(independent["search_scan_file_id"]) - peptide_ids)
    checks.append(
        _check(
            "pairing_search_id_has_data_columns",
            not missing_in_protein and not missing_in_peptide,
            "every independently-derived search_scan_file_id has an Intensity "
            "column in both quant files"
            if not missing_in_protein and not missing_in_peptide
            else f"missing in protein header: {missing_in_protein}; missing in "
            f"peptide header: {missing_in_peptide}",
        )
    )
    return checks


def check_spot_reconciliation(
    name: str,
    ds_abundances: np.ndarray,
    ds_sample_ids: list[str],
    raw_header: list[str],
    raw_rows: list[list[str]],
    id_col_name: str,
    intensity_col_by_sample: dict[str, str],
    *,
    n_checks: int,
    seed: int,
) -> CheckResult:
    """Re-derive >= n_checks random cells from the raw file (csv module) and compare
    to the loaded Dataset (NaN==NaN; positive values to within a 1e-12 relative
    tolerance).

    Not bit-exact: Python's ``float()`` (used for ``expected`` here) and pandas'
    ``pd.to_numeric`` (used by the loader's ``parse_abundance_matrix``) can each
    round the very same decimal string to an adjacent float64 (a documented, benign
    1-ULP difference between the two parsers -- verified directly, e.g.
    ``float("501066.73557059787")`` vs ``pd.to_numeric(...)`` differ by 1 ULP; both
    are within 1e-15 relative of the true value). ``rtol=1e-12`` is far tighter than
    that noise floor yet many orders of magnitude looser than it, so it still catches
    any real corruption (wrong row/column, truncation, unit error) while not flagging
    parser-precision noise as a data-integrity failure.
    """
    rng = random.Random(seed)
    n_features = len(raw_rows)
    mismatches: list[str] = []
    n_checked = 0
    for _ in range(n_checks):
        feature_idx = rng.randrange(n_features)
        sample_id = rng.choice(ds_sample_ids)
        col_name = intensity_col_by_sample[sample_id]
        col_idx = _col_index(raw_header, col_name)
        raw_value = raw_rows[feature_idx][col_idx]
        expected = _raw_cell_to_expected(raw_value)
        sample_row = ds_sample_ids.index(sample_id)
        actual = ds_abundances[sample_row, feature_idx]
        n_checked += 1
        both_nan = np.isnan(expected) and np.isnan(actual)
        close_enough = bool(np.isclose(expected, actual, rtol=1e-12, atol=0.0))
        if not (both_nan or close_enough):
            mismatches.append(
                f"feature={feature_idx} sample={sample_id} raw={raw_value!r} "
                f"expected={expected} actual={actual}"
            )
    return _check(
        f"{name}_spot_reconciliation",
        len(mismatches) == 0 and n_checked >= 200,
        f"checked={n_checked}, mismatches={len(mismatches)}"
        + (f"; e.g. {mismatches[:3]}" if mismatches else ""),
    )


def check_two_way_totals(
    name: str,
    ds_abundances: np.ndarray,
    ds_sample_ids: list[str],
    raw_header: list[str],
    raw_rows: list[list[str]],
    intensity_col_by_sample: dict[str, str],
) -> CheckResult:
    """Per-sample total/median from the Dataset vs. an independent raw-parse sum."""
    mismatches: list[str] = []
    for row_idx, sample_id in enumerate(ds_sample_ids):
        col_name = intensity_col_by_sample[sample_id]
        col_idx = _col_index(raw_header, col_name)
        raw_values = [_raw_cell_to_expected(row[col_idx]) for row in raw_rows]
        raw_arr = np.array(raw_values, dtype=float)
        raw_total = float(np.nansum(raw_arr))
        raw_median = float(np.nanmedian(raw_arr))
        ds_total = float(np.nansum(ds_abundances[row_idx, :]))
        ds_median = float(np.nanmedian(ds_abundances[row_idx, :]))
        if not (
            np.isclose(raw_total, ds_total, rtol=1e-9)
            and np.isclose(raw_median, ds_median, rtol=1e-9)
        ):
            mismatches.append(
                f"sample={sample_id} raw_total={raw_total} ds_total={ds_total} "
                f"raw_median={raw_median} ds_median={ds_median}"
            )
    return _check(
        f"{name}_two_way_total_median",
        len(mismatches) == 0,
        "raw-parse and Dataset totals/medians agree for every sample"
        if not mismatches
        else f"mismatches: {mismatches}",
    )


PSMS_COL_RE = re.compile(r"^PSMs \((.+)\)$")
NSAF_COL_RE = re.compile(r"^NSAF \((.+)\)$")


def check_psms_nsaf_full(
    protein_limelight_file: Path,
    cross_check_file: Path,
    sample_ids: list[str],
    samples_tsv_sample_ids: list[str],
    psms: np.ndarray,
    nsaf: np.ndarray,
) -> list[CheckResult]:
    """Item B (blocking): full, non-sampled reconciliation of every PSMs/NSAF cell.

    Parses the dump with the ``csv`` module (never pandas), maps label -> sample via
    ``cross_check_file`` (also read with ``csv``), and compares every one of the
    4,344 x 8 PSMs cells and NSAF cells against ``limelight.psms`` / ``limelight.nsaf``
    -- no sampling. Also asserts the dump's structural facts (4,344 group rows, 27
    comma-joined groups) and that the NSAF column sums land in the documented
    [0.98, 1.0] range (state/DATA_DESCRIPTION.md "Known data-quality issues").
    """
    header, rows = read_raw_table(protein_limelight_file)
    checks: list[CheckResult] = []

    protein_col_idx = _col_index(header, "Protein(s)")
    n_groups = len(rows)
    checks.append(
        _check(
            "psms_nsaf_n_group_rows",
            n_groups == EXPECTED_N_PROTEINS,
            f"expected {EXPECTED_N_PROTEINS}, got {n_groups}",
        )
    )
    n_comma_joined = sum(1 for row in rows if "," in row[protein_col_idx])
    checks.append(
        _check(
            "psms_nsaf_n_comma_joined_groups",
            n_comma_joined == EXPECTED_N_COMMA_JOINED_GROUPS,
            f"expected {EXPECTED_N_COMMA_JOINED_GROUPS}, got {n_comma_joined}",
        )
    )
    checks.append(
        _check(
            "psms_nsaf_sample_ids_match_samples_tsv_order",
            sample_ids == samples_tsv_sample_ids,
            f"limelight.sample_ids={sample_ids} vs samples.tsv order "
            f"{samples_tsv_sample_ids}",
        )
    )

    with cross_check_file.open(newline="", encoding="utf-8") as handle:
        cc_rows = list(csv.reader(handle, delimiter="\t"))
    cc_header, *cc_data = cc_rows
    label_col = cc_header.index("label")
    sample_col = cc_header.index("sample_id")
    label_to_sample = {row[label_col]: row[sample_col] for row in cc_data if row}

    psms_cols: dict[str, int] = {}
    nsaf_cols: dict[str, int] = {}
    for i, col in enumerate(header):
        m_psms = PSMS_COL_RE.match(col)
        if m_psms is not None:
            psms_cols[m_psms.group(1)] = i
        m_nsaf = NSAF_COL_RE.match(col)
        if m_nsaf is not None:
            nsaf_cols[m_nsaf.group(1)] = i

    mismatches_psms: list[str] = []
    mismatches_nsaf: list[str] = []
    nsaf_col_sums: dict[str, float] = {}
    n_psms_checked = 0
    n_nsaf_checked = 0
    for label, sample_id in label_to_sample.items():
        if (
            sample_id not in sample_ids
            or label not in psms_cols
            or label not in nsaf_cols
        ):
            continue
        row_idx = sample_ids.index(sample_id)
        psms_idx = psms_cols[label]
        nsaf_idx = nsaf_cols[label]
        col_sum = 0.0
        for group_idx, row in enumerate(rows):
            raw_psms = row[psms_idx].replace(",", "").strip()
            expected_psms = 0 if raw_psms == "" else int(raw_psms)
            actual_psms = int(psms[row_idx, group_idx])
            n_psms_checked += 1
            if expected_psms != actual_psms:
                mismatches_psms.append(
                    f"label={label} group={group_idx} expected={expected_psms} "
                    f"actual={actual_psms}"
                )
            raw_nsaf = row[nsaf_idx].replace(",", "").strip()
            expected_nsaf = 0.0 if raw_nsaf == "" else float(raw_nsaf)
            actual_nsaf = float(nsaf[row_idx, group_idx])
            n_nsaf_checked += 1
            if not bool(np.isclose(expected_nsaf, actual_nsaf, rtol=1e-9, atol=1e-12)):
                mismatches_nsaf.append(
                    f"label={label} group={group_idx} expected={expected_nsaf} "
                    f"actual={actual_nsaf}"
                )
            col_sum += expected_nsaf
        nsaf_col_sums[label] = col_sum

    checks.append(
        _check(
            "psms_full_reconciliation",
            len(mismatches_psms) == 0 and n_psms_checked == n_groups * len(sample_ids),
            f"checked={n_psms_checked} (expected {n_groups * len(sample_ids)}), "
            f"mismatches={len(mismatches_psms)}"
            + (f"; e.g. {mismatches_psms[:3]}" if mismatches_psms else ""),
        )
    )
    checks.append(
        _check(
            "nsaf_full_reconciliation",
            len(mismatches_nsaf) == 0 and n_nsaf_checked == n_groups * len(sample_ids),
            f"checked={n_nsaf_checked} (expected {n_groups * len(sample_ids)}), "
            f"mismatches={len(mismatches_nsaf)}"
            + (f"; e.g. {mismatches_nsaf[:3]}" if mismatches_nsaf else ""),
        )
    )
    sums_in_range = all(0.98 <= s <= 1.0 for s in nsaf_col_sums.values())
    checks.append(
        _check(
            "nsaf_column_sums_in_range",
            sums_in_range and len(nsaf_col_sums) == len(sample_ids),
            f"sums={ {k: round(v, 4) for k, v in nsaf_col_sums.items()} }",
        )
    )
    return checks


# --------------------------------------------------------------------------- #
# CLI / orchestration
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Args:
    metadata_file: Path
    mapping_file: Path
    protein_quants_file: Path
    peptide_quants_file: Path
    protein_limelight_file: Path
    peptide_limelight_file: Path
    samples_file: Path
    cross_check_file: Path
    output_dir: Path
    n_spot_checks: int
    seed: int
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
    parser.add_argument(
        "--samples-file", type=Path, default=Path("results/metadata/samples.tsv")
    )
    parser.add_argument(
        "--cross-check-file",
        type=Path,
        default=Path("results/stage2/limelight_label_map.tsv"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/stage3"))
    parser.add_argument("--n-spot-checks", type=int, default=250)
    parser.add_argument("--seed", type=int, default=20260922)
    parser.add_argument("--log-level", type=str.upper, default="INFO")
    ns = parser.parse_args(argv)
    return Args(
        metadata_file=ns.metadata_file,
        mapping_file=ns.mapping_file,
        protein_quants_file=ns.protein_quants_file,
        peptide_quants_file=ns.peptide_quants_file,
        protein_limelight_file=ns.protein_limelight_file,
        peptide_limelight_file=ns.peptide_limelight_file,
        samples_file=ns.samples_file,
        cross_check_file=ns.cross_check_file,
        output_dir=ns.output_dir,
        n_spot_checks=ns.n_spot_checks,
        seed=ns.seed,
        log_level=ns.log_level,
    )


def run_checks(args: Args) -> list[CheckResult]:
    """Run every Stage-3 verification check; never raises itself (returns results)."""
    checks: list[CheckResult] = []

    LOGGER.info("reading raw files via the independent csv-module path")
    protein_header, protein_rows = read_raw_table(args.protein_quants_file)
    peptide_header, peptide_rows = read_raw_table(args.peptide_quants_file)
    checks += check_raw_counts(
        protein_header, protein_rows, peptide_header, peptide_rows
    )

    LOGGER.info("re-deriving pairing independently (metadata.tsv + mapping file)")
    independent_pairing = derive_pairing_independently(
        args.metadata_file, args.mapping_file
    )
    samples = read_samples(args.samples_file)
    checks += check_pairing(
        independent_pairing, samples, protein_header, peptide_header
    )

    LOGGER.info("loading via the project loaders")
    protein_result = load_protein_dataset(
        args.protein_quants_file,
        args.samples_file,
        protein_limelight_file=args.protein_limelight_file,
    )
    pds = protein_result.dataset
    contaminant_ids = frozenset(
        pds.feature_metadata.loc[
            pds.feature_metadata["is_contaminant"], "protein_group"
        ]
    )
    peptide_result = load_peptide_dataset(
        args.peptide_quants_file, args.samples_file, contaminant_ids=contaminant_ids
    )
    eds = peptide_result.dataset

    checks += check_orientation_and_dtype(
        "protein", pds.abundances, pds.feature_names, EXPECTED_N_PROTEINS
    )
    checks += check_orientation_and_dtype(
        "peptide", eds.abundances, eds.feature_names, EXPECTED_N_PEPTIDES
    )
    checks.append(
        check_identifier_roundtrip(
            "protein", pds.feature_names, protein_header, protein_rows, "Protein Groups"
        )
    )
    checks.append(
        check_identifier_roundtrip(
            "peptide", eds.feature_names, peptide_header, peptide_rows, "Sequence"
        )
    )
    checks += check_missing_counts_protein(
        protein_result.was_nan_token, protein_result.was_zero_token
    )
    checks.append(check_missing_counts_peptide(eds.abundances))
    checks.append(check_detection_type_tally(peptide_header, peptide_rows))
    checks += check_contaminant_and_complete_counts(
        "protein",
        pds.feature_metadata,
        pds.abundances,
        EXPECTED_PROTEIN_CONTAMINANTS,
        EXPECTED_PROTEIN_COMPLETE_NONCONTAMINANT,
    )
    checks += check_contaminant_and_complete_counts(
        "peptide",
        eds.feature_metadata,
        eds.abundances,
        EXPECTED_PEPTIDE_CONTAMINANTS,
        EXPECTED_PEPTIDE_COMPLETE_NONCONTAMINANT,
    )

    protein_sample_ids = list(pds.metadata["sample_id"])
    intensity_col_by_sample = (
        dict(zip(samples["sample_id"], samples["intensity_column"], strict=True))
        if "intensity_column" in samples.columns
        else {
            row.sample_id: f"Intensity_search_scan_file_id_{row.search_scan_file_id}"
            for row in samples.itertuples()
        }
    )
    checks.append(
        check_spot_reconciliation(
            "protein",
            pds.abundances,
            protein_sample_ids,
            protein_header,
            protein_rows,
            "Protein Groups",
            intensity_col_by_sample,
            n_checks=args.n_spot_checks,
            seed=args.seed,
        )
    )
    checks.append(
        check_spot_reconciliation(
            "peptide",
            eds.abundances,
            list(eds.metadata["sample_id"]),
            peptide_header,
            peptide_rows,
            "Sequence",
            intensity_col_by_sample,
            n_checks=args.n_spot_checks,
            seed=args.seed + 1,
        )
    )
    checks.append(
        check_two_way_totals(
            "protein",
            pds.abundances,
            protein_sample_ids,
            protein_header,
            protein_rows,
            intensity_col_by_sample,
        )
    )
    checks.append(
        check_two_way_totals(
            "peptide",
            eds.abundances,
            list(eds.metadata["sample_id"]),
            peptide_header,
            peptide_rows,
            intensity_col_by_sample,
        )
    )

    LOGGER.info(
        "resolving Limelight labels (loader cross-checks the Stage-2 map itself)"
    )
    try:
        limelight = load_limelight_protein_counts(
            args.protein_limelight_file,
            args.protein_quants_file,
            args.samples_file,
            cross_check_file=args.cross_check_file,
        )
        checks.append(
            _check(
                "limelight_label_map_reproduced",
                True,
                f"resolved {len(limelight.label_map)} label(s), all matched "
                f"{args.cross_check_file}",
            )
        )
        LOGGER.info("full (non-sampled) PSMs/NSAF reconciliation (item B)")
        checks += check_psms_nsaf_full(
            args.protein_limelight_file,
            args.cross_check_file,
            list(limelight.sample_ids),
            list(samples["sample_id"]),
            limelight.psms,
            limelight.nsaf,
        )
    except Exception as exc:
        # Deliberately broad: caught so this one check reports FAIL rather than
        # aborting the whole verification run before every other check has had a
        # chance to report (fail-loud-but-complete, see module docstring).
        checks.append(_check("limelight_label_map_reproduced", False, str(exc)))

    return checks


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    logging.basicConfig(
        level=args.log_level, format="%(asctime)s %(levelname)s %(message)s"
    )
    LOGGER.info("verify_loaders starting")
    LOGGER.info("params: %s", args)

    # Same six roles, same names, as results/metadata/data_version.json
    # (state/METADATA.md) -- so this script's data_version is directly comparable.
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

    args.output_dir.mkdir(parents=True, exist_ok=True)
    marker_path = args.output_dir / FAILURE_MARKER_NAME
    marker_path.unlink(missing_ok=True)

    try:
        file_hashes = compute_file_hashes(list(data_files.items()))
        data_version = compute_data_version(file_hashes)
        LOGGER.info("data_version = %s", data_version)

        checks = run_checks(args)
        checks.append(
            _check(
                "data_version_matches_metadata_stamp",
                data_version == EXPECTED_DATA_VERSION,
                f"computed={data_version}, expected (results/metadata/"
                f"data_version.json)={EXPECTED_DATA_VERSION}",
            )
        )
        n_passed = sum(c.passed for c in checks)
        for c in checks:
            (LOGGER.info if c.passed else LOGGER.error)(
                "  [%s] %s: %s", "PASS" if c.passed else "FAIL", c.name, c.detail
            )
        LOGGER.info("checks: %d/%d passed", n_passed, len(checks))

        payload: dict[str, Any] = {
            "params": {
                "metadata_file": str(args.metadata_file),
                "mapping_file": str(args.mapping_file),
                "protein_quants_file": str(args.protein_quants_file),
                "peptide_quants_file": str(args.peptide_quants_file),
                "protein_limelight_file": str(args.protein_limelight_file),
                "peptide_limelight_file": str(args.peptide_limelight_file),
                "samples_file": str(args.samples_file),
                "cross_check_file": str(args.cross_check_file),
                "n_spot_checks": args.n_spot_checks,
                "seed": args.seed,
            },
            "expected_data_version": EXPECTED_DATA_VERSION,
            "data_version": data_version,
            "file_hashes": file_hashes,
            "n_checks_passed": int(n_passed),
            "n_checks_total": len(checks),
            "all_passed": n_passed == len(checks),
            "checks": [
                {"name": c.name, "passed": c.passed, "detail": c.detail} for c in checks
            ],
        }
        (args.output_dir / "verification.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8"
        )

        md_lines = [
            "# Stage-3 loader verification",
            "",
            f"data_version: `{data_version}`",
            "",
            f"**{n_passed}/{len(checks)} checks passed.**",
            "",
            "| Check | Result | Detail |",
            "|---|---|---|",
        ]
        for c in checks:
            md_lines.append(
                f"| {c.name} | {'PASS' if c.passed else 'FAIL'} | {c.detail} |"
            )
        (args.output_dir / "verification.md").write_text(
            "\n".join(md_lines) + "\n", encoding="utf-8"
        )

        failed = [c for c in checks if not c.passed]
        if failed:
            detail = "\n".join(f"  - {c.name}: {c.detail}" for c in failed)
            raise ValueError(
                f"{len(failed)} of {len(checks)} Stage-3 verification check(s) "
                f"FAILED:\n{detail}"
            )
    except Exception as exc:
        marker_path.write_text(
            json.dumps(
                {"failed": True, "error_type": type(exc).__name__, "error": str(exc)},
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        LOGGER.error("verify_loaders FAILED: %s", exc)
        raise

    LOGGER.info("wrote outputs to %s", args.output_dir)
    LOGGER.info("verify_loaders done")


if __name__ == "__main__":
    main()
