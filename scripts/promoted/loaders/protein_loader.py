"""Verified loader for the FlashLFQ protein-quants matrix (Stage 3).

PROJECT LOADER seeded from the plugin template ``lib/common/data_loading.py``
(``wide-data-loader`` v0.5) as a guide, not a direct call: ``data/protein-quants.tsv``
does not fit the template's generic ``id_columns`` shape because (a) the missing-value
convention is two distinct raw tokens (``"0"`` and literal text ``"NaN"``) that must
both become ``NaN`` in the returned :class:`~loaders.data_loading.Dataset` while
staying *distinguishable* for provenance, and (b) the feature id needs regex parsing
into structured metadata (accession / entry / contaminant flag / psvid), which the
template does not do. Per the template's own guidance, this is a custom loader that
returns the same ``Dataset`` structure (``state/DATA_DESCRIPTION.md``,
``state/METADATA.md`` §"Join key to the data matrices").

Sample pairing: reads ``results/metadata/samples.tsv`` (the promoted
``metadata_characterize.py`` output — already an exact, verified bijection between
``Replicate`` mzML file names and ``search_scan_file_id``; see
``state/METADATA.md``). Samples are ordered by acquisition order (batch, then
``run_position_within_batch``), re-derived here (not trusted blindly from the file's
row order) so the returned ``Dataset`` is correctly ordered even if ``samples.tsv``
were regenerated in a different order.

Feature id: ``Protein Groups`` (text, read as ``str`` — never type-inferred).
Feature metadata is parsed from the id via the two scientist-confirmed regexes
(``state/DATA_DESCRIPTION.md`` "Contaminants and decoys"):
  * normal:      ``psvid_<n>_sp|<accession>|<entry>``
  * no-accession form: ``psvid_<n>_sp|<entry>|``      (regex
                  ``^psvid_\\d+_sp\\|[^|]+\\|$``)

CONTAMINANT RULE (scientist-confirmed, revised): matching the no-accession form is
only a *candidate*. 19 of the 52 no-accession rows are the Limelight dump's 1:1
counterpart of an indistinguishable-group row whose OTHER member is a real UniProt
accession (e.g. dump row ``sp|CATD_HUMAN|,sp|P07339|CATD_HUMAN``) -- they stand in
for that real protein and must be KEPT, not excluded. A protein is a TRUE
contaminant iff it matches the no-accession form AND its dump group has no other
member with a real accession (33 of the 52). ``load_protein_dataset``'s optional
``protein_limelight_file`` argument applies this refinement (see
:func:`refine_contaminants`); without it, ``is_contaminant`` is the coarser
no-accession-form-only flag (52) and ``contaminant_grouped_with_real`` is always
``False``.

Missing values: raw cells are read as text so the two tokens are told apart before
any numeric coercion. Both ``"0"`` and literal ``"NaN"`` become ``NaN`` in
``Dataset.abundances`` (the project's missing-value policy,
``state/DATA_DESCRIPTION.md`` "Missing-value semantics") — but the loader also
returns two boolean ``(n_samples, n_features)`` masks (``was_nan_token`` /
``was_zero_token``) so the two raw tokens stay distinguishable for provenance/QC,
per this task's requirement.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from loaders.data_loading import Dataset
from loaders.samples import read_samples
from loaders.text_io import assert_no_duplicate_headers, read_raw_header

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": [
        "ProteinLoadResult",
        "load_protein_dataset",
        "PROTEIN_ID_RE",
        "CONTAMINANT_RE",
        "refine_contaminants",
        "parse_abundance_matrix",
    ],
    "uses": ["loaders.data_loading", "loaders.samples", "loaders.text_io"],
    "seeded_from": {"template": "wide-data-loader", "version": "0.5"},
    "description": (
        "Verified loader for data/protein-quants.tsv: feature id = Protein Groups "
        "(text); sample columns paired 1:1 (fail loud on orphan/duplicate) to "
        "results/metadata/samples.tsv via search_scan_file_id; both '0' and literal "
        "'NaN' raw tokens -> NaN in the returned Dataset (linear scale), with a "
        "per-cell boolean token mask kept distinguishable; feature_metadata parsed "
        "into accession/entry/is_contaminant/contaminant_grouped_with_real/psvid. "
        "is_contaminant is the scientist-confirmed dump-refined rule (33 true "
        "contaminants) when protein_limelight_file is given to "
        "load_protein_dataset/refine_contaminants, else the coarser no-accession-"
        "form-only flag (52). Study-specific (not the generic wide-data-loader, "
        "per its own guidance for files that don't fit that shape)."
    ),
}

# psvid_<n>_sp|<mid>|<entry> ; no-accession rows have mid=entry-name and an empty
# trailing entry (the id ends right after the second '|').
PROTEIN_ID_RE = re.compile(r"^psvid_(?P<psvid>\d+)_sp\|(?P<mid>[^|]*)\|(?P<entry>.*)$")
# No-accession candidate form: id ends at the second '|' (nothing after it). A
# candidate is a TRUE contaminant only after refine_contaminants() checks its dump
# group (see the module docstring's CONTAMINANT RULE).
CONTAMINANT_RE = re.compile(r"^psvid_\d+_sp\|[^|]+\|$")
# Strips the leading "psvid_<n>_" from a protein-quants id -> the Limelight dump's
# bare accession-key form (also used by loaders.limelight_loader).
PSVID_PREFIX_RE = re.compile(r"^psvid_\d+_")
# A dump group member with no real accession (bare "sp|<entry>|", nothing after the
# second '|') -- the dump-row analogue of CONTAMINANT_RE, without the psvid prefix.
BARE_NO_ACCESSION_RE = re.compile(r"^sp\|[^|]+\|$")

INTENSITY_COL_RE = re.compile(r"^Intensity_search_scan_file_id_(\d+)$")


@dataclass(frozen=True)
class ProteinLoadResult:
    """The protein :class:`Dataset` plus the raw-token provenance masks.

    Attributes
    ----------
    dataset:
        Linear-scale Dataset; ``abundances`` has ``NaN`` for both raw-missing tokens.
    was_nan_token, was_zero_token:
        ``(n_samples, n_features)`` boolean masks, mutually exclusive, aligned to
        ``dataset.abundances``: ``True`` where the *raw* source cell was literal text
        ``"NaN"`` / numeric ``"0"`` respectively. Both positions are ``NaN`` in
        ``dataset.abundances``; the masks are what keeps the two tokens
        distinguishable (per-cell token record).
    """

    dataset: Dataset
    was_nan_token: np.ndarray
    was_zero_token: np.ndarray


def _read_raw_protein_table(
    data_file: Path, id_column: str
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Read protein-quants.tsv as text; return ``(frame, {sample_id: column_name})``.

    Every cell is read as ``str`` (``dtype=str``, ``keep_default_na=False``) so the
    raw ``"0"`` / literal ``"NaN"`` tokens are told apart from real numeric text
    before any coercion — pandas' own NA-sniffing would otherwise silently convert
    the literal ``"NaN"`` cells to an actual NaN float during the read and the two
    tokens could no longer be distinguished.
    """
    raw_header = read_raw_header(data_file)
    assert_no_duplicate_headers(data_file, raw_header)
    frame = pd.read_csv(data_file, sep="\t", dtype=str, keep_default_na=False)
    if id_column not in frame.columns:
        raise ValueError(f"{data_file} has no id column {id_column!r}.")
    header_ids: dict[str, str] = {}
    for col in frame.columns:
        match = INTENSITY_COL_RE.match(col)
        if match is not None:
            header_ids[match.group(1)] = col
    return frame, header_ids


def _check_bijection(
    data_file: Path, samples_file: Path, header_ids: set[str], sample_ids: set[str]
) -> None:
    """Fail loud on any orphan/duplicate/missing sample<->column pairing."""
    orphan_in_data = sorted(header_ids - sample_ids)
    orphan_in_samples = sorted(sample_ids - header_ids)
    if orphan_in_data or orphan_in_samples:
        raise ValueError(
            f"Sample<->column pairing between {data_file} and {samples_file} is not "
            f"an exact bijection. search_scan_file_id present in the data file but "
            f"not in samples.tsv: {orphan_in_data}; present in samples.tsv but not "
            f"in the data file: {orphan_in_samples}."
        )


def _parse_feature_metadata(protein_ids: pd.Series) -> pd.DataFrame:
    """Parse each ``Protein Groups`` id into accession/entry/is_contaminant/psvid.

    ``is_contaminant`` here is the coarse no-accession-form-only candidate flag;
    call :func:`refine_contaminants` (or pass ``protein_limelight_file`` to
    :func:`load_protein_dataset`) to apply the scientist-confirmed dump-based
    refinement. ``contaminant_grouped_with_real`` defaults to ``False`` and is only
    ever set by that refinement.
    """
    psvid: list[int] = []
    accession: list[str] = []
    entry: list[str] = []
    is_contaminant: list[bool] = []
    for pid in protein_ids:
        match = PROTEIN_ID_RE.match(pid)
        if match is None:
            raise ValueError(
                f"Protein Groups id {pid!r} does not match the expected "
                f"'psvid_<n>_sp|<accession-or-entry>|<entry-or-empty>' pattern."
            )
        contaminant = bool(CONTAMINANT_RE.match(pid))
        psvid.append(int(match.group("psvid")))
        is_contaminant.append(contaminant)
        if contaminant:
            accession.append("")
            entry.append(match.group("mid"))
        else:
            accession.append(match.group("mid"))
            entry.append(match.group("entry"))
    return pd.DataFrame(
        {
            "protein_group": protein_ids.to_numpy(dtype=str),
            "psvid": pd.array(psvid, dtype="int64"),
            "accession": accession,
            "entry": entry,
            "is_contaminant": is_contaminant,
            "contaminant_grouped_with_real": [False] * len(protein_ids),
        }
    )


def _read_dump_group_members(protein_limelight_file: Path) -> dict[str, list[str]]:
    """Map each Limelight dump row's first group member (bare accession key) to the
    row's full comma-split member list.

    The dump<->protein-quants relation is 1:1: each dump row's FIRST member is
    exactly the one protein-quants row it corresponds to (verified: all 4,344 dump
    first-members equal the 4,344 protein-quants accession keys, a bijection; the 28
    other, non-first group members appear in no protein-quants row at all). This is
    a correction of an earlier, wrong "one group -> many quant rows" assumption --
    see loaders/limelight_loader.py's module docstring for the same correction.
    """
    frame = pd.read_csv(
        protein_limelight_file, sep="\t", dtype=str, usecols=["Protein(s)"]
    )
    members = frame["Protein(s)"].str.split(",")
    first = members.apply(lambda m: m[0])
    dup = first[first.duplicated(keep=False)].unique().tolist()
    if dup:
        raise ValueError(
            f"{protein_limelight_file}: {len(dup)} first-group-member value(s) are "
            f"not unique (breaks the assumed 1:1 dump<->protein-quants relation): "
            f"{dup[:5]}."
        )
    return dict(zip(first, members, strict=True))


def refine_contaminants(
    feature_metadata: pd.DataFrame, protein_limelight_file: str | Path
) -> pd.DataFrame:
    """Apply the scientist-confirmed dump-based contaminant refinement.

    A no-accession-form candidate (``feature_metadata["is_contaminant"]`` already
    ``True``, from :func:`_parse_feature_metadata`'s regex-only pass) is a TRUE
    contaminant only if its Limelight dump group (looked up via the verified 1:1
    dump<->protein-quants relation, see :func:`_read_dump_group_members`) has no
    OTHER member with a real UniProt accession. Candidates whose group DOES have a
    real-accession sibling are flipped to ``is_contaminant=False`` and flagged
    ``contaminant_grouped_with_real=True`` (kept, not excluded -- they stand in for
    that real protein; state/DATA_DESCRIPTION.md "Contaminants and decoys", scientist
    decision revising the original all-no-accession-form-is-contaminant rule).

    Returns a new, independent DataFrame (the input is not mutated). Fails loud if a
    candidate's accession key has no corresponding dump row (would break the
    assumed 1:1 relation) or if ``protein_limelight_file`` doesn't parse as expected.
    """
    groups = _read_dump_group_members(Path(protein_limelight_file))
    refined = feature_metadata.copy()
    is_contaminant = refined["is_contaminant"].to_numpy(dtype=bool).copy()
    grouped_with_real = np.zeros(len(refined), dtype=bool)
    keys = refined["protein_group"].apply(lambda x: PSVID_PREFIX_RE.sub("", x))
    for i, (key, candidate) in enumerate(zip(keys, is_contaminant, strict=True)):
        if not candidate:
            continue
        members = groups.get(key)
        if members is None:
            raise ValueError(
                f"No Limelight dump row's first group member equals {key!r} "
                f"(protein-quants row {refined['protein_group'].iloc[i]!r}); cannot "
                f"apply the dump-based contaminant refinement."
            )
        has_real_sibling = any(
            not BARE_NO_ACCESSION_RE.match(member) for member in members[1:]
        )
        if has_real_sibling:
            is_contaminant[i] = False
            grouped_with_real[i] = True
    refined["is_contaminant"] = is_contaminant
    refined["contaminant_grouped_with_real"] = grouped_with_real
    return refined


def parse_abundance_matrix(
    raw: np.ndarray, data_file: Path, nan_token: str = "NaN"
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized text->float parse of a raw string matrix.

    Returns ``(parsed, was_nan_token)``.

    ``raw`` is any-shape array of ``str`` cells. ``nan_token`` cells parse to
    ``float("nan")`` and are flagged in the returned mask; every other cell must be
    numeric-coercible (fail loud, with an example, otherwise) — no other token is a
    silently-accepted missing value.
    """
    was_nan_token = raw == nan_token
    to_parse = np.where(was_nan_token, "nan", raw)
    parsed = pd.to_numeric(pd.Series(to_parse.ravel()), errors="coerce").to_numpy(
        dtype=float
    )
    bad = np.isnan(parsed) & ~was_nan_token.ravel()
    if bool(bad.any()):
        offenders = raw.ravel()[bad]
        raise ValueError(
            f"{data_file} has {int(bad.sum())} non-numeric, non-{nan_token!r} "
            f"abundance value(s), e.g. {offenders[0]!r}."
        )
    return parsed.reshape(raw.shape), was_nan_token


def load_protein_dataset(
    protein_quants_file: str | Path,
    samples_file: str | Path = "results/metadata/samples.tsv",
    *,
    id_column: str = "Protein Groups",
    protein_limelight_file: str | Path | None = None,
) -> ProteinLoadResult:
    """Load ``protein_quants_file`` into the standard :class:`Dataset` contract.

    See the module docstring for the missing-value and feature-metadata policy.
    ``protein_limelight_file``, when given, applies the scientist-confirmed
    dump-based contaminant refinement (:func:`refine_contaminants`) so
    ``is_contaminant`` reflects the 33 TRUE contaminants (not the coarser 52
    no-accession-form candidates) and ``contaminant_grouped_with_real`` is set for
    the 19 kept rows. Fails loud on: a non-unique feature id, an
    orphaned/duplicated/missing sample column, a non-numeric non-"NaN" abundance
    cell, a negative abundance, a Protein Groups id that does not match the
    expected pattern, or (with ``protein_limelight_file``) a contaminant candidate
    whose dump group cannot be found.
    """
    data_file = Path(protein_quants_file)
    samp_file = Path(samples_file)
    if not data_file.is_file():
        raise FileNotFoundError(f"protein_quants_file not found: {data_file}")
    samples = read_samples(samp_file)
    sample_ids = samples["search_scan_file_id"].tolist()

    frame, header_ids = _read_raw_protein_table(data_file, id_column)
    if len(frame) == 0:
        raise ValueError(f"{data_file} has no data rows (header only).")
    _check_bijection(data_file, samp_file, set(header_ids), set(sample_ids))

    protein_ids = frame[id_column]
    blank = protein_ids.str.strip() == ""
    if bool(blank.any()):
        rows = (protein_ids.index[blank][:5] + 2).tolist()
        raise ValueError(
            f"{data_file} has blank {id_column!r} value(s) at row(s) {rows}."
        )
    dup = protein_ids[protein_ids.duplicated(keep=False)].unique().tolist()
    if dup:
        raise ValueError(
            f"{data_file} has duplicate {id_column!r} value(s): {dup[:5]}."
        )

    n_features = len(frame)
    n_samples = len(samples)
    ordered_columns = [header_ids[sid] for sid in sample_ids]
    raw = frame[ordered_columns].to_numpy(dtype=str)  # (n_features, n_samples)

    parsed, was_nan_token = parse_abundance_matrix(raw, data_file)
    finite_negative = np.isfinite(parsed) & (parsed < 0)
    if bool(finite_negative.any()):
        n_bad = int(finite_negative.sum())
        raise ValueError(
            f"{data_file} has {n_bad} negative abundance value(s); intensities must "
            f"be >= 0."
        )
    was_zero_token = (~was_nan_token) & (parsed == 0.0)
    abundances = np.where(was_nan_token | was_zero_token, np.nan, parsed)

    # Transpose to the Dataset convention: (n_samples, n_features).
    abundances_t = abundances.T
    was_nan_token_t = was_nan_token.T
    was_zero_token_t = was_zero_token.T
    if abundances_t.shape != (n_samples, n_features):
        raise ValueError(
            f"Assembled abundances shape {abundances_t.shape} != expected "
            f"({n_samples}, {n_features})."
        )

    feature_metadata = _parse_feature_metadata(protein_ids)
    if protein_limelight_file is not None:
        feature_metadata = refine_contaminants(feature_metadata, protein_limelight_file)
    metadata = samples.copy()
    metadata.index = pd.Index(samples["sample_id"].to_numpy(), name="sample_id")

    dataset = Dataset(
        abundances=abundances_t,
        feature_names=protein_ids.to_numpy(dtype=str),
        feature_metadata=feature_metadata,
        metadata=metadata,
        scale="linear",
    )
    return ProteinLoadResult(
        dataset=dataset,
        was_nan_token=was_nan_token_t,
        was_zero_token=was_zero_token_t,
    )
