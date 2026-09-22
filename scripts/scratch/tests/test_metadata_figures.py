"""Smoke + guard tests for ``scripts/scratch/metadata_figures.py``.

Uses synthetic fixtures (never the real data): a small two-batch design written as the
same precomputed tables ``metadata_characterize.py`` emits, plus a temporary color
registry. Checks that all four figures + legends + provenance are written, that a
precomputed table disagreeing with ``samples.tsv`` fails loud, that a level missing from
the registry fails loud (no ad-hoc colors), and the separation-boundary logic.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import metadata_figures as mf

# batch, seq, position, half, sample_id, condition, pair
_ROWS = [
    ("BA", 10, 1, "early", "S1", "control", "P1"),
    ("BA", 12, 2, "early", "S3", "control", "P2"),
    ("BA", 15, 3, "late", "S2", "raloxifene-d0", "P1"),
    ("BA", 18, 4, "late", "S4", "raloxifene-d0", "P2"),
    ("BB", 3, 1, "early", "S5", "control", "P3"),
    ("BB", 7, 2, "late", "S6", "raloxifene-d0", "P3"),
]

_PALETTE = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"]


def _write_inputs(root: Path) -> tuple[Path, Path]:
    meta = root / "results" / "metadata"
    meta.mkdir(parents=True)
    samples = pd.DataFrame(
        _ROWS,
        columns=[
            "batch",
            "seq_number",
            "run_position_within_batch",
            "run_half",
            "sample_id",
            "condition",
            "candidate_pair",
        ],
    )
    samples.to_csv(meta / "samples.tsv", sep="\t", index=False)
    for column, name in (
        ("batch", "crosstab_condition_batch.tsv"),
        ("run_half", "crosstab_condition_run_half.tsv"),
    ):
        ct = pd.crosstab(samples["condition"], samples[column])
        ct.to_csv(meta / name, sep="\t")
    layout = samples.assign(
        cell=samples["condition"] + "/" + samples["sample_id"]
    ).pivot_table(
        index="batch",
        columns="run_position_within_batch",
        values="cell",
        aggfunc="first",
    )
    layout.to_csv(meta / "run_layout.tsv", sep="\t")
    pd.DataFrame(
        [
            {
                "hypothesis_id": "H4",
                "scope": "batch=BA",
                "p_one_sided": "0.1667",
                "p_two_sided": "0.3333",
            },
            {
                "hypothesis_id": "H4",
                "scope": mf.STRATIFIED_SCOPE,
                "p_one_sided": "0.0833",
                "p_two_sided": "0.1667",
            },
        ]
    ).to_csv(meta / "hypotheses.tsv", sep="\t", index=False)
    (meta / "data_version.json").write_text(json.dumps({"data_version": "sha256:test"}))
    registry = root / "state" / "color_registry.json"
    registry.parent.mkdir(parents=True)
    registry.write_text(
        json.dumps(
            {
                "_palette": {"name": "Okabe-Ito", "colors": _PALETTE},
                "condition": {
                    "scope": "project",
                    "values": {"control": "#0072B2", "raloxifene-d0": "#D55E00"},
                },
                "batch": {
                    "scope": "project",
                    "values": {"BA": "#009E73", "BB": "#E69F00"},
                },
                "candidate_pair": {
                    "scope": "project",
                    "values": {"P1": "#E69F00", "P2": "#56B4E9", "P3": "#F0E442"},
                },
            }
        )
    )
    return meta, registry


def _dirs(root: Path) -> mf.OutputDirs:
    base = root / "figures" / "metadata"
    return mf.OutputDirs(
        base / "distributions", base / "crosstabs", base / "run-layout"
    )


def test_render_all_writes_every_artifact(tmp_path: Path) -> None:
    meta, registry = _write_inputs(tmp_path)
    before = registry.read_text()
    prov_file = tmp_path / "prov.json"
    records = mf.render_all(
        meta, registry, _dirs(tmp_path), prov_file, dpi=50, project_root=tmp_path
    )
    assert set(records) == {
        mf.STEM_COUNTS,
        mf.STEM_BATCH,
        mf.STEM_RUN_HALF,
        mf.STEM_RUN_LAYOUT,
    }
    for rec in records.values():
        for key in ("svg", "png", "legend_svg", "legend_png"):
            path = tmp_path / str(rec[key])
            assert path.is_file() and path.stat().st_size > 0, key
        assert rec["data_version"] == "sha256:test"
    assert json.loads(prov_file.read_text()) == records
    # The registry is read-only for this script.
    assert registry.read_text() == before
    # The run-layout subtitle shows the two-sided stratified p, never a one-sided one.
    layout_svg = (tmp_path / str(records[mf.STEM_RUN_LAYOUT]["svg"])).read_text()
    assert "Exact permutation p = 0.167 (two-sided)" in layout_svg
    assert "one-sided" not in layout_svg


def test_crosstab_disagreement_fails_loud(tmp_path: Path) -> None:
    meta, registry = _write_inputs(tmp_path)
    path = meta / "crosstab_condition_batch.tsv"
    ct = pd.read_csv(path, sep="\t", index_col="condition")
    ct.loc["control", "BA"] += 1
    ct.to_csv(path, sep="\t")
    with pytest.raises(mf.MetadataConsistencyError, match="crosstab_condition_batch"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_run_layout_disagreement_fails_loud(tmp_path: Path) -> None:
    meta, registry = _write_inputs(tmp_path)
    path = meta / "run_layout.tsv"
    text = path.read_text().replace("control/S1", "control/S9")
    path.write_text(text)
    with pytest.raises(mf.MetadataConsistencyError, match="run_layout"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_level_missing_from_registry_fails_loud(tmp_path: Path) -> None:
    meta, registry = _write_inputs(tmp_path)
    reg = json.loads(registry.read_text())
    del reg["candidate_pair"]["values"]["P3"]
    registry.write_text(json.dumps(reg))
    with pytest.raises(KeyError, match="P3"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


@pytest.mark.parametrize(
    ("layout", "expected"),
    [
        ({1: "control", 2: "control", 3: "raloxifene-d0"}, 2.5),
        ({1: "control", 2: "raloxifene-d0", 3: "control"}, None),
        ({1: "raloxifene-d0", 2: "control"}, None),
        ({1: "control", 2: "control"}, None),
    ],
)
def test_separation_boundary(layout: dict[int, str], expected: float | None) -> None:
    assert mf.separation_boundary(layout) == expected
