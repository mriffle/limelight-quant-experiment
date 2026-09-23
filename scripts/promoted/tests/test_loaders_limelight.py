"""Tests for loaders/limelight_loader.py."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from fixtures import SAMPLE_IDS, SEARCH_IDS, intensity_col, write_samples_tsv
from loaders.limelight_loader import load_limelight_protein_counts, sig3

# search_id -> (accession -> raw protein-quants intensity), one vector per sample.
PROTEINS = ["ACC1", "ACC2", "ACC3", "ACC4"]
VALUES: dict[str, dict[str, float]] = {
    "101": {"ACC1": 123000.0, "ACC2": 56700.0, "ACC3": 34500.0, "ACC4": 78100.0},
    "102": {"ACC1": 45600.0, "ACC2": 89100.0, "ACC3": 91200.0, "ACC4": 23400.0},
    "103": {"ACC1": 78900.0, "ACC2": 12300.0, "ACC3": 45600.0, "ACC4": 91000.0},
    "104": {"ACC1": 234000.0, "ACC2": 67800.0, "ACC3": 11100.0, "ACC4": 34400.0},
}
# Deliberately scrambled label -> search_id map the loader must recover.
LABEL_TO_SEARCH_ID = {"L_a": "103", "L_b": "101", "L_c": "104", "L_d": "102"}


def _protein_id(acc: str) -> str:
    return f"psvid_{PROTEINS.index(acc) + 1}_sp|{acc}|E{acc}_HUMAN"


def _write_protein_quants(path: Path) -> None:
    columns = [
        "Protein Groups",
        "Gene Name",
        "Organism",
        *[intensity_col(s) for s in SEARCH_IDS],
    ]
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for acc in PROTEINS:
            row = [_protein_id(acc), "", "", *[str(VALUES[s][acc]) for s in SEARCH_IDS]]
            fh.write("\t".join(row) + "\n")


def _write_dump(
    path: Path,
    labels: dict[str, str],
    *,
    blank_psms_cell: tuple[str, str] | None = None,
    not_quantifiable_cell: tuple[str, str] | None = None,
    include_group_row: bool = True,
) -> None:
    """Write a synthetic protein-limelight-table-dump.txt.

    ``labels``: label -> search_id it truly corresponds to (for the FlashLFQ quant
    values, which the loader must recover by value correspondence).
    ``blank_psms_cell``: optional (label, accession) whose PSMs cell is left blank.
    ``not_quantifiable_cell``: optional (label, accession) whose Quant (FlashLFQ)
    cell is set to the dump's literal "not quantifiable across runs" text instead
    of a number (mirrors the real dump for the 14 literal-"NaN" protein-quants
    cells).
    """
    label_names = list(labels)
    columns = ["Protein(s)", "Protein Description(s)"]
    columns += [f"PSMs ({label})" for label in label_names]
    columns += [f"NSAF ({label})" for label in label_names]
    columns += [f"Quant (FlashLFQ) ({label})" for label in label_names]
    columns += [f"Quant (Limelight) ({label})" for label in label_names]
    columns += ["Protein Group Number"]

    rows: list[list[str]] = []
    for i, acc in enumerate(PROTEINS):
        row: dict[str, str] = {
            "Protein(s)": f"sp|{acc}|E{acc}_HUMAN",
            "Protein Description(s)": "desc",
            "Protein Group Number": str(i + 1),
        }
        for label, search_id in labels.items():
            psms_val = (
                "1,234" if (label, acc) == (label_names[0], PROTEINS[0]) else "56"
            )
            if blank_psms_cell == (label, acc):
                psms_val = ""
            row[f"PSMs ({label})"] = psms_val
            row[f"NSAF ({label})"] = "0.018"
            ff_quant = f"{sig3(np.array([VALUES[search_id][acc]]))[0]:.6g}"
            if not_quantifiable_cell == (label, acc):
                ff_quant = "FlashLFQ: NaN — not quantifiable across runs"
            row[f"Quant (FlashLFQ) ({label})"] = ff_quant
            row[f"Quant (Limelight) ({label})"] = "overlapping signal"
        rows.append([row[c] for c in columns])

    if include_group_row:
        group_row: dict[str, str] = {
            "Protein(s)": "sp|GX1|G1_HUMAN,sp|GX2|G2_HUMAN",
            "Protein Description(s)": "group desc",
            "Protein Group Number": str(len(PROTEINS) + 1),
        }
        for label in label_names:
            group_row[f"PSMs ({label})"] = "10"
            group_row[f"NSAF ({label})"] = "0.001"
            group_row[f"Quant (FlashLFQ) ({label})"] = "1.00e+6"
            group_row[f"Quant (Limelight) ({label})"] = "overlapping signal"
        rows.append([group_row[c] for c in columns])

    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for data_row in rows:
            fh.write("\t".join(data_row) + "\n")


@pytest.fixture
def basic_project(tmp_path: Path) -> tuple[Path, Path, Path]:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    _write_protein_quants(quants_path)
    dump_path = tmp_path / "dump.txt"
    _write_dump(dump_path, LABEL_TO_SEARCH_ID)
    return dump_path, quants_path, samples_path


def test_sig3_rounding() -> None:
    x = np.array([123456.0, 0.0, np.nan, 9.994e-5])
    out = sig3(x)
    assert out[0] == pytest.approx(123000.0)
    assert out[1] == 0.0
    assert math.isnan(out[2])
    # 9.994e-5 -> 3 significant figures -> 9.99e-5 (rounds down at the 4th digit).
    assert out[3] == pytest.approx(9.99e-5, rel=1e-9)


def test_labels_resolved_by_value_correspondence(
    basic_project: tuple[Path, Path, Path],
) -> None:
    dump_path, quants_path, samples_path = basic_project
    result = load_limelight_protein_counts(
        dump_path, quants_path, samples_path, cross_check_file=None
    )
    resolved = dict(
        zip(result.label_map["label"], result.label_map["sample_id"], strict=True)
    )

    search_to_sample = dict(zip(SEARCH_IDS, SAMPLE_IDS, strict=True))
    expected = {
        label: search_to_sample[sid] for label, sid in LABEL_TO_SEARCH_ID.items()
    }
    assert resolved == expected
    assert (result.label_map["match_fraction"] >= 0.999).all()


def test_group_row_kept_with_member_list(
    basic_project: tuple[Path, Path, Path],
) -> None:
    dump_path, quants_path, samples_path = basic_project
    result = load_limelight_protein_counts(
        dump_path, quants_path, samples_path, cross_check_file=None
    )
    idx = list(result.protein_group).index("sp|GX1|G1_HUMAN,sp|GX2|G2_HUMAN")
    assert result.member_accession_keys[idx] == ["sp|GX1|G1_HUMAN", "sp|GX2|G2_HUMAN"]


def test_flashlfq_not_quantifiable_token_does_not_break_resolution(
    tmp_path: Path,
) -> None:
    """One Quant (FlashLFQ) cell is the dump's literal 'not quantifiable' text
    (mirrors FlashLFQ's own not-quantifiable-across-runs case, not a blank); it
    must be skipped from label-matching scoring, not raise or corrupt the result."""
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    _write_protein_quants(quants_path)
    dump_path = tmp_path / "dump.txt"
    label_a = next(iter(LABEL_TO_SEARCH_ID))
    _write_dump(
        dump_path, LABEL_TO_SEARCH_ID, not_quantifiable_cell=(label_a, PROTEINS[0])
    )
    result = load_limelight_protein_counts(
        dump_path, quants_path, samples_path, cross_check_file=None
    )
    resolved = dict(
        zip(result.label_map["label"], result.label_map["sample_id"], strict=True)
    )
    search_to_sample = dict(zip(SEARCH_IDS, SAMPLE_IDS, strict=True))
    expected = {
        label: search_to_sample[sid] for label, sid in LABEL_TO_SEARCH_ID.items()
    }
    assert resolved == expected


def test_psms_thousands_separator_and_blank(
    basic_project: tuple[Path, Path, Path],
) -> None:
    dump_path, quants_path, samples_path = basic_project
    result = load_limelight_protein_counts(
        dump_path, quants_path, samples_path, cross_check_file=None
    )
    label_a = next(iter(LABEL_TO_SEARCH_ID))  # "L_a", the label that got "1,234" PSMs
    sample_row = list(result.sample_ids).index(
        dict(
            zip(result.label_map["label"], result.label_map["sample_id"], strict=True)
        )[label_a]
    )
    protein_col = list(result.protein_group).index(
        f"sp|{PROTEINS[0]}|E{PROTEINS[0]}_HUMAN"
    )
    assert result.psms[sample_row, protein_col] == 1234


def test_blank_psms_cell_is_zero(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    _write_protein_quants(quants_path)
    dump_path = tmp_path / "dump.txt"
    label_a = "L_a"
    _write_dump(dump_path, LABEL_TO_SEARCH_ID, blank_psms_cell=(label_a, PROTEINS[1]))
    result = load_limelight_protein_counts(
        dump_path, quants_path, samples_path, cross_check_file=None
    )
    sample_id = dict(
        zip(result.label_map["label"], result.label_map["sample_id"], strict=True)
    )[label_a]
    sample_row = list(result.sample_ids).index(sample_id)
    protein_col = list(result.protein_group).index(
        f"sp|{PROTEINS[1]}|E{PROTEINS[1]}_HUMAN"
    )
    assert result.psms[sample_row, protein_col] == 0


def test_ambiguous_label_fails_loud(tmp_path: Path) -> None:
    """Two samples with identical protein-quants vectors: no label can be resolved."""
    samples_path = tmp_path / "samples_2.tsv"
    df = pd.DataFrame(
        {
            "sample_id": ["AZ001", "AZ002"],
            "search_scan_file_id": ["101", "102"],
            "batch": ["B1", "B1"],
            "run_position_within_batch": ["1", "2"],
            "condition": ["control", "raloxifene-d0"],
        }
    )
    df.to_csv(samples_path, sep="\t", index=False)

    quants_path = tmp_path / "protein-quants_2.tsv"
    columns = [
        "Protein Groups",
        "Gene Name",
        "Organism",
        intensity_col("101"),
        intensity_col("102"),
    ]
    with quants_path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for acc in PROTEINS:
            fh.write(
                "\t".join([_protein_id(acc), "", "", "1000.0", "1000.0"]) + "\n"
            )  # identical across both samples -> ambiguous

    dump_path = tmp_path / "dump_2.txt"
    label_names = ["L_x", "L_y"]  # both must match both samples equally -> ambiguous
    columns2 = ["Protein(s)", "Protein Description(s)"]
    columns2 += [f"PSMs ({label})" for label in label_names]
    columns2 += [f"NSAF ({label})" for label in label_names]
    columns2 += [f"Quant (FlashLFQ) ({label})" for label in label_names]
    columns2 += [f"Quant (Limelight) ({label})" for label in label_names]
    with dump_path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns2) + "\n")
        for acc in PROTEINS:
            row = {
                "Protein(s)": f"sp|{acc}|E{acc}_HUMAN",
                "Protein Description(s)": "desc",
            }
            for label in label_names:
                row[f"PSMs ({label})"] = "10"
                row[f"NSAF ({label})"] = "0.01"
                row[f"Quant (FlashLFQ) ({label})"] = (
                    f"{sig3(np.array([1000.0]))[0]:.6g}"
                )
                row[f"Quant (Limelight) ({label})"] = "overlapping signal"
            fh.write("\t".join(row[c] for c in columns2) + "\n")

    with pytest.raises(ValueError, match="unambiguously"):
        load_limelight_protein_counts(
            dump_path, quants_path, samples_path, cross_check_file=None
        )


def test_cross_check_mismatch_fails_loud(
    basic_project: tuple[Path, Path, Path],
) -> None:
    dump_path, quants_path, samples_path = basic_project

    wrong_map = pd.DataFrame(
        {
            "label": list(LABEL_TO_SEARCH_ID),
            "sample_id": [SAMPLE_IDS[0]] * len(LABEL_TO_SEARCH_ID),
        }
    )
    cross_check_path = dump_path.parent / "wrong_label_map.tsv"
    wrong_map.to_csv(cross_check_path, sep="\t", index=False)
    with pytest.raises(ValueError, match="disagrees"):
        load_limelight_protein_counts(
            dump_path, quants_path, samples_path, cross_check_file=cross_check_path
        )


def test_cross_check_agreement_passes(basic_project: tuple[Path, Path, Path]) -> None:
    dump_path, quants_path, samples_path = basic_project
    result = load_limelight_protein_counts(
        dump_path, quants_path, samples_path, cross_check_file=None
    )
    cross_check_path = dump_path.parent / "correct_label_map.tsv"
    result.label_map[["label", "sample_id"]].to_csv(
        cross_check_path, sep="\t", index=False
    )
    # Should not raise.
    load_limelight_protein_counts(
        dump_path, quants_path, samples_path, cross_check_file=cross_check_path
    )


def test_missing_dump_file_raises(tmp_path: Path) -> None:
    samples_path = tmp_path / "samples.tsv"
    write_samples_tsv(samples_path)
    quants_path = tmp_path / "protein-quants.tsv"
    _write_protein_quants(quants_path)
    with pytest.raises(FileNotFoundError):
        load_limelight_protein_counts(
            tmp_path / "no-such-dump.txt",
            quants_path,
            samples_path,
            cross_check_file=None,
        )
