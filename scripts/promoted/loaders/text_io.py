"""Shared raw-header reading for the project's wide data-file loaders.

Factored out (script-registry: shared code lives in one module, imported, never
copy-pasted) because both ``protein_loader.py`` and ``peptide_loader.py`` need the
same fail-loud guard: pandas silently renames a duplicated column header
(``"Intensity_..._101"`` -> ``"Intensity_..._101.1"``) rather than raising, which
would otherwise make ``INTENSITY_COL_RE`` simply fail to match the *second* column
and its data vanish with no error at all -- exactly the "silent row/column drop"
the correctness charter forbids. Reading the raw header line first, before pandas
gets to it, catches this.
"""

from __future__ import annotations

from pathlib import Path

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": ["read_raw_header", "assert_no_duplicate_headers"],
    "uses": [],
    "seeded_from": None,
    "description": (
        "Reads just the raw first (header) line of a wide data file and asserts no "
        "column name repeats, before pandas can silently rename a duplicate and "
        "drop its data. Shared by protein_loader.py and peptide_loader.py."
    ),
}


def read_raw_header(path: Path, *, sep: str = "\t") -> list[str]:
    """Return the raw first-line column names of ``path``, before pandas parses it."""
    with path.open(encoding="utf-8") as handle:
        first_line = handle.readline()
    if not first_line:
        raise ValueError(f"{path} is empty (no header line).")
    return first_line.rstrip("\r\n").split(sep)


def assert_no_duplicate_headers(path: Path, header: list[str]) -> None:
    """Raise if ``header`` (the raw column names) contains a repeated name."""
    seen: set[str] = set()
    dupes: set[str] = set()
    for name in header:
        if name in seen:
            dupes.add(name)
        seen.add(name)
    if dupes:
        raise ValueError(
            f"{path} has duplicate column header(s): {sorted(dupes)}. pandas would "
            f"silently rename these on read (e.g. 'col' -> 'col.1'), dropping the "
            f"duplicate's data from any regex-based column match."
        )
