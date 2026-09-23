"""Verified loader for the FlashLFQ peptide-quants matrix (Stage 3).

PROJECT LOADER seeded from the plugin template ``lib/common/data_loading.py``
(``wide-data-loader`` v0.5) as a guide, not a direct call — same rationale as
``protein_loader.py`` (see its module docstring): ``data/peptide-quants.tsv`` needs
custom feature-id parsing (base sequence + trailing total-modification-mass bracket),
contaminant propagation from the semicolon-joined ``Protein Groups`` members, and a
per-cell ``Detection Type`` matrix the template has no concept of. Reuses
``protein_loader.parse_abundance_matrix`` (the token-aware text->float parser) and
``protein_loader.CONTAMINANT_RE`` (the scientist-confirmed contaminant regex) rather
than duplicating them — one implementation, per the script-registry convention.

Feature id: ``Sequence`` = base sequence + one trailing bracketed total modification
mass (e.g. ``LCR[+114.042928]``; ``state/DATA_DESCRIPTION.md`` "Feature
identifiers"). Cross-checked against the file's own ``Base Sequence`` column
(fail loud on disagreement).

Missing values: unlike protein-quants, peptide-quants documents only ``"0"`` as a
missing token (no literal ``"NaN"``; ``state/DATA_DESCRIPTION.md`` "Missing-value
semantics"). The loader still scans for a literal ``"NaN"`` token (reusing the same
parser as the protein loader) but treats its presence as a violated assumption and
fails loud, rather than silently reusing the protein policy on undocumented ground.

Detection Type invariant (Stage 2, re-asserted here as a loader-level integrity
check, not just documentation): ``intensity > 0 <=> Detection Type in {MSMS, MBR}``
and ``intensity == 0 <=> Detection Type in {NotDetected,
MSMSIdentifiedButNotQuantified, MSMSAmbiguousPeakfinding}``. Any other Detection
Type value, or any violation of the bijection, fails loud.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from loaders.data_loading import Dataset
from loaders.protein_loader import CONTAMINANT_RE, parse_abundance_matrix
from loaders.samples import read_samples
from loaders.text_io import assert_no_duplicate_headers, read_raw_header

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": [
        "PeptideLoadResult",
        "load_peptide_dataset",
        "PEPTIDE_ID_RE",
        "POSITIVE_DETECTION_TYPES",
        "ZERO_DETECTION_TYPES",
    ],
    "uses": [
        "loaders.data_loading",
        "loaders.protein_loader",
        "loaders.samples",
        "loaders.text_io",
    ],
    "seeded_from": {"template": "wide-data-loader", "version": "0.5"},
    "description": (
        "Verified loader for data/peptide-quants.tsv: feature id = Sequence (base "
        "sequence + trailing total-mod-mass bracket, cross-checked against Base "
        "Sequence); sample columns paired 1:1 to results/metadata/samples.tsv; '0' "
        "-> NaN in the returned Dataset (linear scale); feature_metadata carries "
        "base_sequence/protein_groups/is_contaminant/total_mod_mass. is_contaminant "
        "uses the scientist-confirmed dump-refined 33-id contaminant set when "
        "contaminant_ids is given (any ;-member in that set), else the coarser "
        "no-accession-form-only regex match. Also returns the per-cell Detection "
        "Type matrix and asserts the Stage-2 intensity<->Detection Type invariant, "
        "fail-loud. Reuses protein_loader's token-aware parser and contaminant "
        "regex (no duplication)."
    ),
}

INTENSITY_COL_RE = re.compile(r"^Intensity_search_scan_file_id_(\d+)$")
DETECTION_COL_RE = re.compile(r"^Detection Type_search_scan_file_id_(\d+)$")

# base sequence (upper-case amino acids) + exactly one trailing "[+<mass>]".
PEPTIDE_ID_RE = re.compile(r"^(?P<base>[A-Z]+)\[\+(?P<mass>[0-9.]+)\]$")

POSITIVE_DETECTION_TYPES: frozenset[str] = frozenset({"MSMS", "MBR"})
ZERO_DETECTION_TYPES: frozenset[str] = frozenset(
    {"NotDetected", "MSMSIdentifiedButNotQuantified", "MSMSAmbiguousPeakfinding"}
)
VALID_DETECTION_TYPES: frozenset[str] = POSITIVE_DETECTION_TYPES | ZERO_DETECTION_TYPES


@dataclass(frozen=True)
class PeptideLoadResult:
    """The peptide :class:`Dataset` plus the per-cell Detection Type matrix.

    Attributes
    ----------
    dataset:
        Linear-scale Dataset; ``abundances`` has ``NaN`` where the raw cell was ``0``
        (not quantified).
    detection_type:
        ``(n_samples, n_features)`` str array, same orientation as
        ``dataset.abundances``, carrying the raw ``Detection Type`` value for every
        cell (including the quantified ones).
    """

    dataset: Dataset
    detection_type: np.ndarray


def _header_ids(header: list[str], pattern: re.Pattern[str]) -> dict[str, str]:
    """Map ``search_scan_file_id`` text -> matching column name, from a header list."""
    out: dict[str, str] = {}
    for col in header:
        match = pattern.match(col)
        if match is not None:
            out[match.group(1)] = col
    return out


def _check_bijection(
    data_file: Path,
    samples_file: Path,
    header_ids: set[str],
    sample_ids: set[str],
    column_family: str,
) -> None:
    """Fail loud on any orphan/duplicate/missing sample<->column pairing."""
    orphan_in_data = sorted(header_ids - sample_ids)
    orphan_in_samples = sorted(sample_ids - header_ids)
    if orphan_in_data or orphan_in_samples:
        raise ValueError(
            f"{column_family} sample<->column pairing between {data_file} and "
            f"{samples_file} is not an exact bijection. search_scan_file_id present "
            f"in the data file but not in samples.tsv: {orphan_in_data}; present in "
            f"samples.tsv but not in the data file: {orphan_in_samples}."
        )


def _parse_feature_metadata(
    sequences: pd.Series,
    base_sequence_col: pd.Series,
    protein_groups_col: pd.Series,
    contaminant_ids: frozenset[str] | None,
) -> pd.DataFrame:
    """Parse each ``Sequence`` into base/mass; validate against Base Sequence; flag
    contaminants from the semicolon-joined ``Protein Groups`` members.

    ``contaminant_ids``, when given, is the scientist-confirmed set of TRUE
    (dump-refined) contaminant ``Protein Groups`` ids (see
    ``protein_loader.refine_contaminants``): a peptide is a contaminant iff any
    ``;``-member is IN this set (excludes only the 33 pure contaminants, keeping
    peptides whose sole "contaminant-looking" member is one of the 19 kept
    grouped-with-real rows). ``None`` falls back to the coarser
    no-accession-form-only regex match (any member matches ``CONTAMINANT_RE``).
    """
    base_parsed: list[str] = []
    mass: list[float] = []
    for seq, base_expected in zip(sequences, base_sequence_col, strict=True):
        match = PEPTIDE_ID_RE.match(seq)
        if match is None:
            raise ValueError(
                f"Sequence {seq!r} does not match the expected "
                f"'<BASE>[+<mass>]' pattern."
            )
        base = match.group("base")
        if base != base_expected:
            raise ValueError(
                f"Sequence {seq!r} parses to base {base!r} but the file's own "
                f"'Base Sequence' column says {base_expected!r}."
            )
        base_parsed.append(base)
        mass.append(float(match.group("mass")))

    if contaminant_ids is None:
        is_contaminant = [
            any(CONTAMINANT_RE.match(member) for member in groups.split(";"))
            for groups in protein_groups_col
        ]
    else:
        is_contaminant = [
            any(member in contaminant_ids for member in groups.split(";"))
            for groups in protein_groups_col
        ]
    return pd.DataFrame(
        {
            "sequence": sequences.to_numpy(dtype=str),
            "base_sequence": base_parsed,
            "protein_groups": protein_groups_col.to_numpy(dtype=str),
            "is_contaminant": is_contaminant,
            "total_mod_mass": pd.array(mass, dtype="float64"),
        }
    )


def _assert_detection_invariant(
    intensity: np.ndarray, detection_type: np.ndarray, data_file: Path
) -> None:
    """Fail loud unless ``intensity > 0 <=> detection_type in POSITIVE_DETECTION_TYPES``
    and ``intensity == 0 <=> detection_type in ZERO_DETECTION_TYPES`` hold everywhere,
    and every Detection Type value is one of the five documented types."""
    unknown = set(np.unique(detection_type)) - VALID_DETECTION_TYPES
    if unknown:
        raise ValueError(
            f"{data_file} has undocumented Detection Type value(s): {sorted(unknown)}."
        )
    is_positive_type = np.isin(detection_type, list(POSITIVE_DETECTION_TYPES))
    is_zero_type = np.isin(detection_type, list(ZERO_DETECTION_TYPES))
    is_positive_intensity = intensity > 0

    violation_a = is_positive_type & ~is_positive_intensity
    violation_b = is_positive_intensity & ~is_positive_type
    violation_c = is_zero_type & (intensity != 0)
    violation_d = (intensity == 0) & ~is_zero_type
    n_bad = int((violation_a | violation_b | violation_c | violation_d).sum())
    if n_bad:
        raise ValueError(
            f"{data_file}: {n_bad} cell(s) violate the intensity<->Detection Type "
            f"invariant (intensity > 0 <=> type in {sorted(POSITIVE_DETECTION_TYPES)}; "
            f"intensity == 0 <=> type in {sorted(ZERO_DETECTION_TYPES)})."
        )


def load_peptide_dataset(
    peptide_quants_file: str | Path,
    samples_file: str | Path = "results/metadata/samples.tsv",
    *,
    id_column: str = "Sequence",
    contaminant_ids: frozenset[str] | None = None,
) -> PeptideLoadResult:
    """Load ``peptide_quants_file`` into the standard :class:`Dataset` contract.

    See the module docstring for the missing-value, feature-metadata, and
    Detection-Type-invariant policy. ``contaminant_ids``, when given, applies the
    scientist-confirmed dump-refined contaminant set (see
    ``protein_loader.refine_contaminants`` -- pass the ``protein_group`` ids where
    its ``is_contaminant`` is ``True``); otherwise the coarser no-accession-form-only
    regex rule is used. Fails loud on: a non-unique feature id, an
    orphaned/duplicated/missing sample column (either the Intensity or the Detection
    Type family), a non-numeric abundance cell, an undocumented literal ``"NaN"``
    token, a negative abundance, a Sequence that disagrees with Base Sequence, an
    undocumented Detection Type value, or an intensity<->Detection-Type violation.
    """
    data_file = Path(peptide_quants_file)
    samp_file = Path(samples_file)
    if not data_file.is_file():
        raise FileNotFoundError(f"peptide_quants_file not found: {data_file}")
    samples = read_samples(samp_file)
    sample_ids = samples["search_scan_file_id"].tolist()

    raw_header = read_raw_header(data_file)
    assert_no_duplicate_headers(data_file, raw_header)
    frame = pd.read_csv(data_file, sep="\t", dtype=str, keep_default_na=False)
    for required in (id_column, "Base Sequence", "Protein Groups"):
        if required not in frame.columns:
            raise ValueError(f"{data_file} has no column {required!r}.")
    if len(frame) == 0:
        raise ValueError(f"{data_file} has no data rows (header only).")

    intensity_ids = _header_ids(list(frame.columns), INTENSITY_COL_RE)
    detection_ids = _header_ids(list(frame.columns), DETECTION_COL_RE)
    _check_bijection(
        data_file, samp_file, set(intensity_ids), set(sample_ids), "Intensity"
    )
    _check_bijection(
        data_file, samp_file, set(detection_ids), set(sample_ids), "Detection Type"
    )

    sequences = frame[id_column]
    blank = sequences.str.strip() == ""
    if bool(blank.any()):
        rows = (sequences.index[blank][:5] + 2).tolist()
        raise ValueError(
            f"{data_file} has blank {id_column!r} value(s) at row(s) {rows}."
        )
    dup = sequences[sequences.duplicated(keep=False)].unique().tolist()
    if dup:
        raise ValueError(
            f"{data_file} has duplicate {id_column!r} value(s): {dup[:5]}."
        )

    n_features = len(frame)
    n_samples = len(samples)
    intensity_columns = [intensity_ids[sid] for sid in sample_ids]
    detection_columns = [detection_ids[sid] for sid in sample_ids]

    raw_intensity = frame[intensity_columns].to_numpy(dtype=str)  # (n_features, n_s)
    detection_type = frame[detection_columns].to_numpy(dtype=str)  # (n_features, n_s)

    parsed, was_nan_token = parse_abundance_matrix(raw_intensity, data_file)
    if bool(was_nan_token.any()):
        raise ValueError(
            f"{data_file} has {int(was_nan_token.sum())} literal 'NaN' abundance "
            f"token(s); peptide-quants is documented to use only '0' for "
            f"not-quantified (state/DATA_DESCRIPTION.md). This is an unexpected, "
            f"undocumented token — resolve before loading."
        )
    finite_negative = np.isfinite(parsed) & (parsed < 0)
    if bool(finite_negative.any()):
        raise ValueError(
            f"{data_file} has {int(finite_negative.sum())} negative abundance "
            f"value(s); intensities must be >= 0."
        )

    _assert_detection_invariant(parsed, detection_type, data_file)

    abundances = np.where(parsed == 0.0, np.nan, parsed)

    # Transpose to the Dataset convention: (n_samples, n_features).
    abundances_t = abundances.T
    detection_type_t = detection_type.T
    if abundances_t.shape != (n_samples, n_features):
        raise ValueError(
            f"Assembled abundances shape {abundances_t.shape} != expected "
            f"({n_samples}, {n_features})."
        )

    feature_metadata = _parse_feature_metadata(
        sequences, frame["Base Sequence"], frame["Protein Groups"], contaminant_ids
    )
    metadata = samples.copy()
    metadata.index = pd.Index(samples["sample_id"].to_numpy(), name="sample_id")

    dataset = Dataset(
        abundances=abundances_t,
        feature_names=sequences.to_numpy(dtype=str),
        feature_metadata=feature_metadata,
        metadata=metadata,
        scale="linear",
    )
    return PeptideLoadResult(dataset=dataset, detection_type=detection_type_t)
