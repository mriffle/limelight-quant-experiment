"""Shared ``results/metadata/samples.tsv`` reader for every project loader.

One implementation, imported by ``protein_loader.py``, ``peptide_loader.py``, and
``limelight_loader.py`` (script-registry convention: shared code lives in one module
and is imported, never copy-pasted). Not seeded from a ``lib/`` template — the
``samples.tsv`` shape (columns, acquisition-order sort keys) is entirely
project-specific (produced by the promoted ``metadata_characterize.py``).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": ["read_samples", "REQUIRED_SAMPLE_COLUMNS"],
    "uses": [],
    "seeded_from": None,
    "description": (
        "Reads + validates results/metadata/samples.tsv (required columns present, "
        "sample_id and search_scan_file_id unique) and returns it sorted into "
        "acquisition order (batch, then run_position_within_batch) -- re-derived "
        "here rather than trusted from the file's row order. Shared by every "
        "project loader that pairs a quant/dump file to samples.tsv."
    ),
}

REQUIRED_SAMPLE_COLUMNS = (
    "sample_id",
    "search_scan_file_id",
    "batch",
    "run_position_within_batch",
    "condition",
)


def read_samples(samples_file: str | Path) -> pd.DataFrame:
    """Read + validate ``samples_file``; return it ordered by acquisition order.

    Raises ``FileNotFoundError`` if the file is absent, ``ValueError`` if a required
    column is missing, the table is empty, or ``sample_id`` / ``search_scan_file_id``
    is not unique.
    """
    path = Path(samples_file)
    if not path.is_file():
        raise FileNotFoundError(f"samples_file not found: {path}")
    samples = pd.read_csv(path, sep="\t", dtype=str)
    missing_cols = [c for c in REQUIRED_SAMPLE_COLUMNS if c not in samples.columns]
    if missing_cols:
        raise ValueError(f"{path} is missing required column(s): {missing_cols}")
    if samples.empty:
        raise ValueError(f"{path} has no data rows.")
    for col in ("sample_id", "search_scan_file_id"):
        dup = samples.loc[samples[col].duplicated(keep=False), col].unique().tolist()
        if dup:
            raise ValueError(f"{path} has duplicate {col!r} value(s): {dup}")
    ordered = samples.assign(
        _run_position=pd.to_numeric(
            samples["run_position_within_batch"], errors="raise"
        )
    ).sort_values(["batch", "_run_position"], kind="stable")
    return ordered.drop(columns="_run_position").reset_index(drop=True)
