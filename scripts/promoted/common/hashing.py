"""File-hashing and data-version helpers shared by every Stage-1 script.

One implementation, used by both ``metadata_characterize.py`` (the loader/analysis
script, which hashes the raw ``data/`` inputs) and ``metadata_figures.py`` (which
hashes its own inputs — the precomputed ``results/metadata`` tables, the color
registry, and its own code modules — for per-figure provenance). Held to the
correctness charter (conventions/correctness.md): **assume nothing, verify
everything, fail loud.**

Design note: hashes are keyed by a caller-supplied ROLE (e.g. ``metadata_file``,
or a project-relative path string), never by bare filename — two inputs that
happen to share a basename must never collide/overwrite each other's hash, and
the combined ``data_version`` stamp is computed over sorted ``"<role>:<sha256>"``
lines so it is stable regardless of where the files happen to live on disk.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable, Sequence
from pathlib import Path

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": [
        "sha256_of_file",
        "compute_file_hashes",
        "compute_data_version",
        "duplicates",
    ],
    "uses": [],
    "seeded_from": None,
    "description": (
        "sha256 file hashing + role-keyed data-version stamp, shared by "
        "metadata_characterize.py and metadata_figures.py so there is exactly "
        "one hashing implementation in the project."
    ),
}


def sha256_of_file(path: Path, *, chunk_size: int = 1 << 20) -> str:
    """Stream-hash ``path`` (raw bytes; never parsed/interpreted) with sha256."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def duplicates(items: Iterable[str]) -> list[str]:
    """Return the values that appear more than once in ``items`` (sorted, unique)."""
    seen: set[str] = set()
    dupes: set[str] = set()
    for item in items:
        if item in seen:
            dupes.add(item)
        seen.add(item)
    return sorted(dupes)


def compute_file_hashes(
    data_files: Sequence[tuple[str, Path]],
) -> dict[str, dict[str, str]]:
    """Hash each ``(role, path)`` pair; returns ``{role: {"path": ..., "sha256": ..}}``.

    Keyed by ROLE, not by filename, so two files that happen to share a basename
    cannot collide/overwrite each other's hash. Raises if a role appears more than
    once (the combined stamp below assumes exactly one hash per role).
    """
    roles = [role for role, _ in data_files]
    dup_roles = duplicates(roles)
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
