"""Tests for ``scripts/scratch/common/hashing.py``."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from common import hashing


def test_sha256_of_file_matches_known_digest(tmp_path: Path) -> None:
    path = tmp_path / "a.txt"
    path.write_text("hello", encoding="utf-8")
    # sha256("hello") — a well-known test vector.
    assert (
        hashing.sha256_of_file(path)
        == "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
    )


def test_compute_file_hashes_keyed_by_role_with_path(tmp_path: Path) -> None:
    f1 = tmp_path / "a.txt"
    f1.write_text("hello", encoding="utf-8")
    f2 = tmp_path / "b.txt"
    f2.write_text("world", encoding="utf-8")
    out = hashing.compute_file_hashes([("role_a", f1), ("role_b", f2)])
    assert out["role_a"] == {
        "path": str(f1),
        "sha256": hashing.sha256_of_file(f1),
    }
    assert out["role_b"]["sha256"] == hashing.sha256_of_file(f2)


def test_compute_file_hashes_rejects_duplicate_roles(tmp_path: Path) -> None:
    f1 = tmp_path / "a.txt"
    f1.write_text("hello", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate"):
        hashing.compute_file_hashes([("role_a", f1), ("role_a", f1)])


def test_compute_data_version_role_order_independent_and_content_sensitive(
    tmp_path: Path,
) -> None:
    f1 = tmp_path / "a.txt"
    f1.write_text("hello", encoding="utf-8")
    f2 = tmp_path / "b.txt"
    f2.write_text("world", encoding="utf-8")

    hashes_1 = hashing.compute_file_hashes([("role_a", f1), ("role_b", f2)])
    hashes_2 = hashing.compute_file_hashes([("role_b", f2), ("role_a", f1)])
    assert hashing.compute_data_version(hashes_1) == hashing.compute_data_version(
        hashes_2
    )
    assert hashing.compute_data_version(hashes_1).startswith("sha256:")

    f1.write_text("HELLO", encoding="utf-8")
    hashes_3 = hashing.compute_file_hashes([("role_a", f1), ("role_b", f2)])
    assert hashing.compute_data_version(hashes_3) != hashing.compute_data_version(
        hashes_1
    )


def test_compute_data_version_keyed_by_role_not_path(tmp_path: Path) -> None:
    """Two different paths under the same role produce the same stamp component
    key (role), so the stamp depends on role + content, never on where a file
    happens to live on disk."""
    dir_a = tmp_path / "a"
    dir_a.mkdir()
    dir_b = tmp_path / "b"
    dir_b.mkdir()
    f_a = dir_a / "x.txt"
    f_a.write_text("same content", encoding="utf-8")
    f_b = dir_b / "y.txt"
    f_b.write_text("same content", encoding="utf-8")

    hashes_a = hashing.compute_file_hashes([("role", f_a)])
    hashes_b = hashing.compute_file_hashes([("role", f_b)])
    assert hashing.compute_data_version(hashes_a) == hashing.compute_data_version(
        hashes_b
    )
