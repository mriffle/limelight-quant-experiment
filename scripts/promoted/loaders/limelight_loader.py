"""Verified loader for the Limelight protein table dump (PSMs / NSAF), Stage 3.

PROJECT LOADER for ``data/protein-limelight-table-dump.txt``. Not seeded from a
``lib/`` template — nothing in ``lib/`` addresses resolving ambiguous per-run column
labels by VALUE CORRESPONDENCE against an already-loaded matrix, which is the whole
point of this loader (``state/DATA_DESCRIPTION.md`` "Sample identifiers ->
metadata"). Reuses ``protein_loader.load_protein_dataset`` (for the reference
FlashLFQ intensities the labels are matched against) rather than re-parsing
``protein-quants.tsv`` a second time.

Why the return type is NOT the project's ``Dataset`` contract: PSMs/NSAF are keyed
by Limelight **protein group rows**. The dump<->protein-quants relation is actually
1:1 (verified): each of the 4,344 dump rows' FIRST comma-separated member is exactly
the one ``protein-quants`` row it corresponds to, and this is a bijection over all
4,344 rows of both files. (An earlier draft of this loader/``state/
DATA_DESCRIPTION.md`` wrongly described this as "one group -> many quant rows" /
"the quant file splits indistinguishable groups into separate rows" -- that was
incorrect and has been corrected: the 28 OTHER, non-first members across the 27
comma-joined groups appear in no ``protein-quants`` row at all.) Despite the 1:1
relation, PSMs/NSAF are still returned as a dedicated tidy container
(:class:`LimelightProteinCounts`) rather than the ``Dataset`` contract: they are
count/fraction measures joined to ``protein_loader``'s feature space by accession
key (not the primary feature space here), and ``member_accession_keys`` carries the
extra, quant-less group members for provenance/QC visibility -- information the
``Dataset`` contract has no slot for.

Label resolution (``state/DATA_DESCRIPTION.md`` "Limelight dumps: resolved"): the
dump's per-run column suffix (e.g. ``"1_0506_A"``) carries no direct sample id. Each
label is assigned to the one ``protein-quants`` sample whose sig3-rounded
(3-significant-figure, Limelight's display precision) intensities match the dump's
own ``Quant (FlashLFQ)`` column for that label, over every single-member protein row
in common. A label is accepted only if its best match_fraction is >= 0.999 AND its
runner-up is <= 0.9 (an unambiguous winner); otherwise this loader fails loud rather
than guessing. The result is additionally cross-checked, byte-for-byte on the
label -> sample_id assignment, against the Stage-2 artifact
``results/stage2/limelight_label_map.tsv`` (skippable via ``cross_check_file=None``,
e.g. for synthetic test fixtures that have no Stage-2 counterpart).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from loaders.protein_loader import PSVID_PREFIX_RE, load_protein_dataset
from loaders.samples import read_samples
from loaders.text_io import assert_no_duplicate_headers, read_raw_header

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": [
        "LimelightProteinCounts",
        "load_limelight_protein_counts",
        "sig3",
    ],
    "uses": ["loaders.protein_loader", "loaders.samples", "loaders.text_io"],
    # PSVID_PREFIX_RE is imported from loaders.protein_loader (single implementation).
    "seeded_from": None,
    "description": (
        "Verified loader for data/protein-limelight-table-dump.txt: PSMs (thousands- "
        "separator int; blank -> 0) and NSAF (3 significant figures; blank -> 0) per "
        "sample x protein-group, plus the label map. Per-run column labels are "
        "resolved to samples by VALUE CORRESPONDENCE against protein-quants "
        "intensities (sig3-rounded exact match, >= 0.999 best / <= 0.9 runner-up), "
        "cross-checked against results/stage2/limelight_label_map.tsv. 27 comma-"
        "joined indistinguishable-group rows keep the group string and an explicit "
        "member accession-key list; the dump<->protein-quants relation is 1:1 via "
        "each row's first member (verified bijection; the 27 groups' 28 other "
        "members appear in no protein-quants row). Returns a dedicated tidy "
        "container, not the Dataset contract -- see the module docstring."
    ),
}

LABEL_RE = re.compile(r"\(([^()]+)\)$")

DEFAULT_CROSS_CHECK_FILE = "results/stage2/limelight_label_map.tsv"

# The dump's own literal text for FlashLFQ's "not quantifiable across runs" case
# (state/DATA_DESCRIPTION.md "Missing-value semantics"), matching the 14 literal
# "NaN" cells in protein-quants -- a Quant (FlashLFQ) token, never blank.
FLASHLFQ_NOT_QUANTIFIABLE_TOKEN = "FlashLFQ: NaN — not quantifiable across runs"


@dataclass(frozen=True)
class LimelightProteinCounts:
    """PSMs/NSAF per sample x Limelight protein-group row, plus the label map.

    Attributes
    ----------
    protein_group:
        ``(n_groups,)`` str array — the dump's raw ``Protein(s)`` value (the join
        key), verbatim (a single accession-key string, or a comma-joined group).
    member_accession_keys:
        ``(n_groups,)`` list of lists — each group's members, split on ``","``, in
        the ``sp|<accession>|<entry>`` form matching ``protein-quants`` ids with
        their ``psvid_<n>_`` prefix stripped (for joining to
        ``protein_loader``-loaded protein rows).
    psms:
        ``(n_samples, n_groups)`` int array. Blank raw cells are ``0`` PSMs.
    nsaf:
        ``(n_samples, n_groups)`` float array (3-significant-figure display
        precision, as shipped). Blank raw cells are ``0.0``.
    sample_ids:
        ``(n_samples,)`` str array, in acquisition order (matches ``psms``/``nsaf``
        row order).
    label_map:
        One row per resolved dump label: ``label``, ``sample_id``, ``match_fraction``,
        ``second_best``, ``condition``, ``batch``, ``search_scan_file_id``.
    """

    protein_group: np.ndarray
    member_accession_keys: list[list[str]]
    psms: np.ndarray
    nsaf: np.ndarray
    sample_ids: np.ndarray
    label_map: pd.DataFrame


def sig3(x: np.ndarray) -> np.ndarray:
    """Round to 3 significant figures (Limelight's display precision).

    ``0`` stays exactly ``0``; ``NaN`` propagates as ``NaN`` (never silently becomes
    ``0``, unlike a naive "``x>0`` mask, ``else`` leave the rest" implementation).
    """
    x = np.asarray(x, dtype=float)
    out = np.full(x.shape, np.nan, dtype=float)
    out[x == 0] = 0.0
    positive = np.isfinite(x) & (x > 0)
    if bool(positive.any()):
        magnitude = np.floor(np.log10(x[positive]))
        out[positive] = np.round(x[positive] / 10**magnitude, 2) * 10**magnitude
    return out


def _extract_labels(header: list[str], family_prefix: str) -> dict[str, str]:
    """Map label -> column name for columns named ``"<family_prefix> (<label>)"``."""
    out: dict[str, str] = {}
    prefix = f"{family_prefix} ("
    for col in header:
        if not col.startswith(prefix):
            continue
        match = LABEL_RE.search(col)
        if match is None:
            raise ValueError(f"Column {col!r} does not end in '(<label>)'.")
        out[match.group(1)] = col
    return out


def _parse_numeric_column(
    raw: pd.Series,
    data_file: Path,
    column: str,
    *,
    nan_tokens: frozenset[str] = frozenset(),
) -> np.ndarray:
    """Parse a thousands-separated numeric column; blank -> 0.0; fail loud otherwise.

    ``nan_tokens`` declares additional literal strings that parse to ``NaN`` rather
    than raising -- e.g. the dump's ``Quant (FlashLFQ)`` column spells FlashLFQ's own
    "not quantifiable across runs" case as literal text (matching the 14
    protein-quants literal-``"NaN"`` cells; ``state/DATA_DESCRIPTION.md`` "Missing-
    value semantics"), not a blank.
    """
    cleaned = raw.str.replace(",", "", regex=False)
    cleaned = cleaned.where(cleaned.str.strip() != "", "0")
    is_nan_token = cleaned.isin(nan_tokens)
    cleaned = cleaned.where(~is_nan_token, "nan")
    parsed = pd.to_numeric(cleaned, errors="coerce")
    bad = parsed.isna() & ~is_nan_token
    if bool(bad.any()):
        offenders = raw[bad].unique().tolist()
        raise ValueError(
            f"{data_file}, column {column!r}: {int(bad.sum())} value(s) are not "
            f"numeric (after stripping ',' and treating blank as 0), e.g. "
            f"{offenders[:3]}."
        )
    return np.asarray(parsed.to_numpy(dtype=float))


def _parse_psms_column(raw: pd.Series, data_file: Path, column: str) -> np.ndarray:
    """Parse a PSMs column to non-negative integers (blank -> 0 PSMs)."""
    values = _parse_numeric_column(raw, data_file, column)
    non_integer = values % 1 != 0
    if bool(non_integer.any()):
        raise ValueError(
            f"{data_file}, column {column!r}: {int(non_integer.sum())} non-integer "
            f"PSMs value(s)."
        )
    negative = values < 0
    if bool(negative.any()):
        raise ValueError(
            f"{data_file}, column {column!r}: {int(negative.sum())} negative PSMs "
            f"value(s)."
        )
    return values.astype(np.int64)


def _accession_key(protein_group_id: str) -> str:
    """Strip the ``psvid_<n>_`` prefix: protein-quants id -> dump accession key."""
    return PSVID_PREFIX_RE.sub("", protein_group_id)


def _resolve_labels(
    dump: pd.DataFrame,
    dump_path: Path,
    labels: list[str],
    ff_quant_ids: dict[str, str],
    pq_by_key: pd.DataFrame,
    *,
    match_threshold: float,
    second_best_threshold: float,
) -> pd.DataFrame:
    """Resolve each dump label to the one sample whose sig3(intensity) it matches.

    ``pq_by_key`` is (accession_key x sample_id), the reference intensities.
    Returns a frame with columns ``label``, ``sample_id``, ``match_fraction``,
    ``second_best``, indexed by label. Fails loud if any label's best match is below
    ``match_threshold`` or its runner-up exceeds ``second_best_threshold``, or if the
    resulting label -> sample_id assignment is not a bijection.
    """
    common = dump.index.intersection(pq_by_key.index)
    if len(common) == 0:
        raise ValueError(
            "No dump 'Protein(s)' rows (single-member) match any protein-quants "
            "accession key; cannot resolve Limelight labels by value correspondence."
        )
    sample_ids = list(pq_by_key.columns)
    score = pd.DataFrame(index=labels, columns=sample_ids, dtype=float)
    for label in labels:
        d = _parse_numeric_column(
            dump.loc[common, ff_quant_ids[label]],
            dump_path,
            ff_quant_ids[label],
            nan_tokens=frozenset({FLASHLFQ_NOT_QUANTIFIABLE_TOKEN}),
        )
        ok = ~np.isnan(d)
        for sample_id in sample_ids:
            q = sig3(pq_by_key.loc[common, sample_id].to_numpy(dtype=float))
            score.loc[label, sample_id] = float(
                np.isclose(d[ok], q[ok], rtol=1e-6).mean()
            )

    best_sample = score.idxmax(axis=1)
    best_score = score.max(axis=1)
    second_best = score.apply(lambda r: sorted(r)[-2], axis=1)

    ambiguous = (best_score < match_threshold) | (second_best > second_best_threshold)
    if bool(ambiguous.any()):
        bad_labels = score.index[ambiguous].tolist()
        raise ValueError(
            f"Limelight label(s) {bad_labels} could not be resolved unambiguously: "
            f"best match_fraction must be >= {match_threshold} and runner-up <= "
            f"{second_best_threshold}. Scores:\n{score.loc[bad_labels].round(4)}"
        )

    resolved = pd.DataFrame(
        {
            "label": labels,
            "sample_id": best_sample.to_numpy(),
            "match_fraction": best_score.to_numpy(),
            "second_best": second_best.to_numpy(),
        }
    )
    dup_samples = (
        resolved.loc[resolved["sample_id"].duplicated(keep=False), "sample_id"]
        .unique()
        .tolist()
    )
    if dup_samples:
        raise ValueError(
            f"Resolved label -> sample_id assignment is not a bijection; sample(s) "
            f"claimed by more than one label: {dup_samples}."
        )
    return resolved.set_index("label")


def _cross_check(resolved: pd.DataFrame, cross_check_file: str | Path) -> None:
    """Fail loud if ``resolved`` disagrees with the Stage-2 label map on disk."""
    path = Path(cross_check_file)
    if not path.is_file():
        raise FileNotFoundError(f"cross_check_file not found: {path}")
    reference = pd.read_csv(path, sep="\t", dtype=str).set_index("label")
    ours = resolved["sample_id"].astype(str)
    theirs = reference["sample_id"].astype(str)
    if set(ours.index) != set(theirs.index):
        raise ValueError(
            f"Label set resolved here {sorted(ours.index)} does not match "
            f"{path} {sorted(theirs.index)}."
        )
    disagreements = {
        label: (ours[label], theirs[label])
        for label in ours.index
        if ours[label] != theirs[label]
    }
    if disagreements:
        raise ValueError(
            f"Resolved Limelight label map disagrees with {path} for label(s): "
            f"{disagreements}."
        )


def load_limelight_protein_counts(
    protein_limelight_file: str | Path,
    protein_quants_file: str | Path,
    samples_file: str | Path = "results/metadata/samples.tsv",
    *,
    cross_check_file: str | Path | None = DEFAULT_CROSS_CHECK_FILE,
    match_threshold: float = 0.999,
    second_best_threshold: float = 0.9,
) -> LimelightProteinCounts:
    """Load PSMs/NSAF per sample x Limelight protein-group, resolving labels by value.

    See the module docstring for the resolution rule and why the return type is a
    dedicated container rather than the project's ``Dataset`` contract. Fails loud
    on: an unresolvable/ambiguous label, a non-bijective label map, disagreement
    with ``cross_check_file`` (unless ``None``), a non-numeric PSMs/NSAF/quant cell,
    or a structurally inconsistent dump header (missing PSMs/NSAF/Quant(FlashLFQ)
    family for a label).
    """
    dump_path = Path(protein_limelight_file)
    if not dump_path.is_file():
        raise FileNotFoundError(f"protein_limelight_file not found: {dump_path}")

    samples = read_samples(samples_file)
    sample_ids = samples["sample_id"].tolist()

    protein_result = load_protein_dataset(protein_quants_file, samples_file)
    ds = protein_result.dataset
    # Reconstruct the pre-missing-conversion values (real zeros distinguishable from
    # NaN tokens) from the token masks, since sig3-matching needs the raw zero.
    pre_conversion = np.where(
        protein_result.was_nan_token,
        np.nan,
        np.where(protein_result.was_zero_token, 0.0, ds.abundances),
    )
    accession_keys = [_accession_key(pid) for pid in ds.feature_names]
    if len(set(accession_keys)) != len(accession_keys):
        raise ValueError(
            "protein-quants accession keys (psvid prefix stripped) collide."
        )
    # ds.metadata is indexed by sample_id in ds.abundances row order; pre_conversion
    # is (n_samples, n_features) in that same order.
    pq_by_key = pd.DataFrame(
        pre_conversion.T, index=accession_keys, columns=list(ds.metadata["sample_id"])
    )
    # Re-order columns to this loader's own sample_ids order (identical set; asserted).
    if set(pq_by_key.columns) != set(sample_ids):
        raise ValueError("protein_loader sample_ids disagree with samples.tsv.")
    pq_by_key = pq_by_key[sample_ids]

    raw_header = read_raw_header(dump_path)
    assert_no_duplicate_headers(dump_path, raw_header)
    dump = pd.read_csv(dump_path, sep="\t", dtype=str, keep_default_na=False)
    protein_group_col = "Protein(s)"
    if protein_group_col not in dump.columns:
        raise ValueError(f"{dump_path} has no column {protein_group_col!r}.")
    dup_groups = (
        dump[protein_group_col][dump[protein_group_col].duplicated(keep=False)]
        .unique()
        .tolist()
    )
    if dup_groups:
        raise ValueError(
            f"{dump_path} has duplicate {protein_group_col!r} row(s): {dup_groups[:5]}."
        )
    dump = dump.set_index(protein_group_col, drop=False)

    header = list(dump.columns)
    psms_ids = _extract_labels(header, "PSMs")
    nsaf_ids = _extract_labels(header, "NSAF")
    ff_quant_ids = _extract_labels(header, "Quant (FlashLFQ)")
    labels = sorted(psms_ids)
    if not labels:
        raise ValueError(f"{dump_path} has no 'PSMs (<label>)' columns.")
    for name, ids in (
        ("PSMs", psms_ids),
        ("NSAF", nsaf_ids),
        ("Quant (FlashLFQ)", ff_quant_ids),
    ):
        if set(ids) != set(labels):
            raise ValueError(
                f"{dump_path}: {name!r} column labels {sorted(ids)} != PSMs labels "
                f"{labels}."
            )
    if len(labels) != len(sample_ids):
        raise ValueError(
            f"{dump_path} has {len(labels)} distinct label(s) but samples.tsv has "
            f"{len(sample_ids)} sample(s)."
        )

    resolved = _resolve_labels(
        dump,
        dump_path,
        labels,
        ff_quant_ids,
        pq_by_key,
        match_threshold=match_threshold,
        second_best_threshold=second_best_threshold,
    )
    if cross_check_file is not None:
        _cross_check(resolved, cross_check_file)

    label_map = resolved.reset_index().merge(
        samples[["sample_id", "condition", "batch", "search_scan_file_id"]],
        on="sample_id",
        how="left",
    )

    n_groups = len(dump)
    psms = np.empty((len(sample_ids), n_groups), dtype=np.int64)
    nsaf = np.empty((len(sample_ids), n_groups), dtype=float)
    for row_idx, sample_id in enumerate(sample_ids):
        label = resolved.index[resolved["sample_id"] == sample_id][0]
        psms[row_idx, :] = _parse_psms_column(
            dump[psms_ids[label]], dump_path, psms_ids[label]
        )
        nsaf[row_idx, :] = _parse_numeric_column(
            dump[nsaf_ids[label]], dump_path, nsaf_ids[label]
        )

    protein_group = dump[protein_group_col].to_numpy(dtype=str)
    member_accession_keys = [row.split(",") for row in protein_group]

    return LimelightProteinCounts(
        protein_group=protein_group,
        member_accession_keys=member_accession_keys,
        psms=psms,
        nsaf=nsaf,
        sample_ids=np.array(sample_ids, dtype=str),
        label_map=label_map,
    )
