"""Tests for ``scripts/scratch/metadata_figures.py``.

Two fixture families:

* ``_build_characterize_output`` — runs the REAL ``metadata_characterize.main()``
  on synthetic (never real-data) inputs and returns the ``results/metadata``
  directory it wrote. Used for every scenario the real pipeline can actually
  produce (a clean render; an odd-sized batch; a failed upstream run's
  ``FAILED.json``; a degenerate single-sample cohort's blank stratified
  p-value; provenance's ``data_version`` matching ``data_version.json``).
* ``_write_raw_tables`` — writes the precomputed tables directly (bypassing
  ``metadata_characterize.py`` entirely), for scenarios the real pipeline's own
  validity checks would never let through (a hand-tampered/corrupted table: a
  crosstab or run_layout disagreeing with samples.tsv, a duplicate sample_id,
  an unknown condition/run_half level, a missing/duplicate stratified H4 row,
  an empty samples.tsv) — i.e. simulating a bug or manual edit downstream of a
  trustworthy characterize run, which this script's own defenses must catch.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import metadata_characterize as mc
import metadata_figures as mf
from common import design
from common.figures.colors import CategoricalPaletteExceededError

_PALETTE = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#000000",
]


# --------------------------------------------------------------------------- #
# Fixture family 1: run the real metadata_characterize.main()
# --------------------------------------------------------------------------- #


def _filename(year: int, month: int, day: int, seq: int, sample_num: int) -> str:
    return (
        f"UWPRExp480_{year:04d}_{month:02d}{day:02d}_AZ_{seq:03d}_"
        f"AZ{sample_num:03d}_AZ_complex.mzML"
    )


def _build_characterize_output(
    tmp_path: Path, rows: list[tuple[int, int, int, int, int, str]]
) -> Path:
    """Run ``metadata_characterize.main()`` on a synthetic cohort; return its
    output dir (``<tmp_path>/results/metadata``).

    ``rows`` is ``(year, month, day, seq, sample_number, condition)`` tuples.
    """
    in_dir = tmp_path / "in"
    in_dir.mkdir(parents=True, exist_ok=True)
    files = [_filename(y, m, d, seq, sn) for (y, m, d, seq, sn, _c) in rows]
    conditions = [c for (*_rest, c) in rows]
    lines = ["Replicate\tcondition"] + [
        f"{f}\t{c}" for f, c in zip(files, conditions, strict=True)
    ]
    (in_dir / "metadata.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")

    ids = list(range(1, len(files) + 1))
    mapping_text = ", ".join(f"{i} ({f})" for i, f in zip(ids, files, strict=True))
    (in_dir / "mapping.txt").write_text(mapping_text, encoding="utf-8")

    intensity_cols = "\t".join(f"Intensity_search_scan_file_id_{i}" for i in ids)
    detection_cols = "\t".join(f"Detection Type_search_scan_file_id_{i}" for i in ids)
    (in_dir / "protein-quants.tsv").write_text(
        f"Protein Groups\t{intensity_cols}\n", encoding="utf-8"
    )
    (in_dir / "peptide-quants.tsv").write_text(
        f"Sequence\t{intensity_cols}\t{detection_cols}\n", encoding="utf-8"
    )

    psms_cols = "\t".join(f"PSMs ({i})" for i in ids)
    nsaf_cols = "\t".join(f"NSAF ({i})" for i in ids)
    (in_dir / "protein-limelight.txt").write_text(
        f"Protein(s)\t{psms_cols}\t{nsaf_cols}\n", encoding="utf-8"
    )
    quant_cols = "\t".join(f"Quant ({i})" for i in ids)
    (in_dir / "peptide-limelight.txt").write_text(
        f"Peptide Sequence\t{psms_cols}\t{quant_cols}\n", encoding="utf-8"
    )

    out_dir = tmp_path / "results" / "metadata"
    argv = [
        "--metadata-file",
        str(in_dir / "metadata.tsv"),
        "--mapping-file",
        str(in_dir / "mapping.txt"),
        "--protein-quants-file",
        str(in_dir / "protein-quants.tsv"),
        "--peptide-quants-file",
        str(in_dir / "peptide-quants.tsv"),
        "--protein-limelight-file",
        str(in_dir / "protein-limelight.txt"),
        "--peptide-limelight-file",
        str(in_dir / "peptide-limelight.txt"),
        "--output-dir",
        str(out_dir),
        "--log-level",
        "WARNING",
    ]
    mc.main(argv)
    return out_dir


def _write_registry_for(metadata_dir: Path, registry_path: Path) -> None:
    """Seed a registry with exactly the levels present in ``metadata_dir``'s
    ``samples.tsv`` (plus both canonical conditions, so a single-condition
    cohort still has both colors available)."""
    samples = pd.read_csv(metadata_dir / "samples.tsv", sep="\t", dtype=str)
    batches = sorted(samples["batch"].unique())
    pairs = sorted(samples["candidate_pair"].unique())
    conditions = sorted(set(samples["condition"].unique()) | set(mf.CONDITION_ORDER))
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    registry_path.write_text(
        json.dumps(
            {
                "_palette": {"name": "Okabe-Ito", "colors": _PALETTE},
                "condition": {
                    "scope": "project",
                    "values": {
                        c: _PALETTE[i % len(_PALETTE)] for i, c in enumerate(conditions)
                    },
                },
                "batch": {
                    "scope": "project",
                    "values": {
                        b: _PALETTE[i % len(_PALETTE)] for i, b in enumerate(batches)
                    },
                },
                "candidate_pair": {
                    "scope": "project",
                    "values": {
                        p: _PALETTE[i % len(_PALETTE)] for i, p in enumerate(pairs)
                    },
                },
            }
        ),
        encoding="utf-8",
    )


def _dirs(root: Path) -> mf.OutputDirs:
    base = root / "figures" / "metadata"
    return mf.OutputDirs(
        base / "distributions", base / "crosstabs", base / "run-layout"
    )


def _default_cohort() -> list[tuple[int, int, int, int, int, str]]:
    """4 samples, 1 batch: a fully confounded run order (like the real design)."""
    return [
        (2021, 1, 1, 1, 201, "control"),
        (2021, 1, 1, 2, 203, "control"),
        (2021, 1, 1, 3, 202, "raloxifene-d0"),
        (2021, 1, 1, 4, 204, "raloxifene-d0"),
    ]


# --------------------------------------------------------------------------- #
# render_all — the clean/successful path
# --------------------------------------------------------------------------- #


def test_render_all_writes_every_artifact(tmp_path: Path) -> None:
    metadata_dir = _build_characterize_output(tmp_path, _default_cohort())
    registry = tmp_path / "state" / "color_registry.json"
    _write_registry_for(metadata_dir, registry)
    before = registry.read_text(encoding="utf-8")
    prov_file = tmp_path / "prov.json"

    records = mf.render_all(
        metadata_dir,
        registry,
        _dirs(tmp_path),
        prov_file,
        dpi=50,
        project_root=tmp_path,
    )
    assert set(records) == {
        mf.STEM_COUNTS,
        mf.STEM_BATCH,
        mf.STEM_RUN_HALF,
        mf.STEM_RUN_LAYOUT,
    }
    data_version = json.loads(
        (metadata_dir / "data_version.json").read_text(encoding="utf-8")
    )["data_version"]
    for rec in records.values():
        for key in ("svg", "png", "legend_svg", "legend_png"):
            path = tmp_path / str(rec[key])
            assert path.is_file() and path.stat().st_size > 0, key
        # provenance data_version must match data_version.json exactly.
        assert rec["data_version"] == data_version
        # inputs always carries the registry + both figure-code modules' hashes,
        # in addition to whichever data tables that figure used.
        inputs = rec["inputs"]
        assert isinstance(inputs, dict)
        assert any(str(k).endswith("color_registry.json") for k in inputs)
        assert any(str(k).endswith("colors.py") for k in inputs)
        assert any(str(k).endswith("figure_io.py") for k in inputs)
    assert json.loads(prov_file.read_text(encoding="utf-8")) == records
    # The registry is read-only for this script.
    assert registry.read_text(encoding="utf-8") == before
    # The run-layout subtitle names the test as stratified-by-batch and shows
    # the two-sided p, never a one-sided one.
    layout_svg = (tmp_path / str(records[mf.STEM_RUN_LAYOUT]["svg"])).read_text(
        encoding="utf-8"
    )
    # _default_cohort() is 1 batch, 2v2, fully separated -> u_max=4, u_obs=4,
    # p_two_sided = 2/6 (the 2v2 exact-permutation case): 2 of the 6 possible
    # label assignments are at least this separated (this one, and its mirror).
    assert "Stratified-by-batch exact permutation p = 0.333 (two-sided)" in layout_svg
    assert "one-sided" not in layout_svg
    run_half_svg = (tmp_path / str(records[mf.STEM_RUN_HALF]["svg"])).read_text(
        encoding="utf-8"
    )
    assert "presumed acquisition order" in run_half_svg


def test_render_all_git_commit_is_populated_in_a_repo(tmp_path: Path) -> None:
    """The project is a git repo; script.git_commit must be a real commit hash,
    not null, when rendering from within it."""
    metadata_dir = _build_characterize_output(tmp_path, _default_cohort())
    registry = tmp_path / "state" / "color_registry.json"
    _write_registry_for(metadata_dir, registry)
    project_root = Path(__file__).resolve().parents[3]
    records = mf.render_all(
        metadata_dir,
        registry,
        _dirs(tmp_path),
        tmp_path / "prov.json",
        dpi=50,
        project_root=project_root,
    )
    script = next(iter(records.values()))["script"]
    assert isinstance(script, dict)
    commit = script["git_commit"]
    assert isinstance(commit, str) and len(commit) == 40


# --------------------------------------------------------------------------- #
# B1 — odd-sized batch: run_half must use common.design (middle -> early)
# --------------------------------------------------------------------------- #


def test_render_all_succeeds_for_odd_sized_batch(tmp_path: Path) -> None:
    rows = [
        (2021, 1, 1, 1, 201, "control"),
        (2021, 1, 1, 2, 202, "control"),
        (2021, 1, 1, 3, 203, "raloxifene-d0"),
    ]
    metadata_dir = _build_characterize_output(tmp_path, rows)
    samples = pd.read_csv(metadata_dir / "samples.tsv", sep="\t")
    # Middle position (2 of 3) must be "early" per common.design.run_half_label.
    assert design.run_half_label(2, 3) == "early"
    by_pos = dict(
        zip(samples["run_position_within_batch"], samples["run_half"], strict=True)
    )
    assert by_pos == {1: "early", 2: "early", 3: "late"}

    registry = tmp_path / "state" / "color_registry.json"
    _write_registry_for(metadata_dir, registry)
    # Must not raise: metadata_figures' own re-derivation of run_half (via the
    # same common.design function) must agree with what characterize wrote.
    mf.render_all(
        metadata_dir, registry, _dirs(tmp_path), tmp_path / "prov.json", dpi=50
    )


# --------------------------------------------------------------------------- #
# B2 — refuse to render if FAILED.json is present
# --------------------------------------------------------------------------- #


def test_render_all_refuses_when_failed_marker_present(tmp_path: Path) -> None:
    metadata_dir = _build_characterize_output(tmp_path, _default_cohort())
    registry = tmp_path / "state" / "color_registry.json"
    _write_registry_for(metadata_dir, registry)
    (metadata_dir / design.FAILURE_MARKER_NAME).write_text(
        json.dumps({"failed": True, "error": "boom"}), encoding="utf-8"
    )
    with pytest.raises(mf.UpstreamRunFailedError, match=r"FAILED\.json"):
        mf.render_all(
            metadata_dir, registry, _dirs(tmp_path), tmp_path / "prov.json", dpi=50
        )


def test_render_all_refuses_from_a_genuinely_failed_characterize_run(
    tmp_path: Path,
) -> None:
    """End-to-end: metadata_characterize.main() itself writes FAILED.json on a
    bad cohort, and metadata_figures.py refuses to render from it."""
    bad_rows = [(2021, 1, 1, 1, 201, "BOGUS"), (2021, 1, 1, 2, 202, "control")]
    metadata_dir = tmp_path / "results" / "metadata"
    with pytest.raises(ValueError, match="validity check"):
        _build_characterize_output(tmp_path, bad_rows)
    assert (metadata_dir / design.FAILURE_MARKER_NAME).is_file()
    registry = tmp_path / "state" / "color_registry.json"
    registry.parent.mkdir(parents=True, exist_ok=True)
    registry.write_text(
        json.dumps({"_palette": {"name": "Okabe-Ito", "colors": _PALETTE}}),
        encoding="utf-8",
    )
    with pytest.raises(mf.UpstreamRunFailedError):
        mf.render_all(
            metadata_dir, registry, _dirs(tmp_path), tmp_path / "prov.json", dpi=50
        )


# --------------------------------------------------------------------------- #
# B5 — blank p_two_sided / single-sample cohort (A1's u_max==0 branch)
# --------------------------------------------------------------------------- #


def test_render_all_rejects_blank_stratified_p_two_sided(tmp_path: Path) -> None:
    """Direct unit test of the guard: a hand-blanked p_two_sided must raise a
    clear error, not crash on ``float("")``."""
    metadata_dir = _build_characterize_output(tmp_path, _default_cohort())
    hyp_path = metadata_dir / "hypotheses.tsv"
    hyp = pd.read_csv(hyp_path, sep="\t", dtype=str, keep_default_na=False)
    mask = (hyp["hypothesis_id"] == "H4") & (hyp["scope"] == mf.STRATIFIED_SCOPE)
    assert mask.any()
    hyp.loc[mask, "p_two_sided"] = ""
    hyp.to_csv(hyp_path, sep="\t", index=False)

    registry = tmp_path / "state" / "color_registry.json"
    _write_registry_for(metadata_dir, registry)
    with pytest.raises(mf.MetadataConsistencyError, match="blank"):
        mf.render_all(
            metadata_dir, registry, _dirs(tmp_path), tmp_path / "prov.json", dpi=50
        )


def test_render_all_rejects_single_sample_cohort_blank_p(tmp_path: Path) -> None:
    """A genuine 1-sample cohort: u_max=0 (A1) -> hypotheses.tsv's stratified
    p_two_sided is blank -> metadata_figures.py must fail loud, not render a
    figure with a fabricated p-value."""
    metadata_dir = _build_characterize_output(
        tmp_path, [(2021, 1, 1, 1, 201, "control")]
    )
    hyp = pd.read_csv(
        metadata_dir / "hypotheses.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    strat = hyp[hyp["scope"] == mf.STRATIFIED_SCOPE]
    assert len(strat) == 1
    assert strat["p_two_sided"].iloc[0] == ""

    registry = tmp_path / "state" / "color_registry.json"
    _write_registry_for(metadata_dir, registry)
    with pytest.raises(mf.MetadataConsistencyError, match="blank"):
        mf.render_all(
            metadata_dir, registry, _dirs(tmp_path), tmp_path / "prov.json", dpi=50
        )


def test_render_all_rejects_missing_stratified_row(tmp_path: Path) -> None:
    metadata_dir = _build_characterize_output(tmp_path, _default_cohort())
    hyp_path = metadata_dir / "hypotheses.tsv"
    hyp = pd.read_csv(hyp_path, sep="\t", dtype=str, keep_default_na=False)
    hyp = hyp[hyp["scope"] != mf.STRATIFIED_SCOPE]
    hyp.to_csv(hyp_path, sep="\t", index=False)

    registry = tmp_path / "state" / "color_registry.json"
    _write_registry_for(metadata_dir, registry)
    with pytest.raises(mf.MetadataConsistencyError, match="expected exactly one"):
        mf.render_all(
            metadata_dir, registry, _dirs(tmp_path), tmp_path / "prov.json", dpi=50
        )


def test_render_all_rejects_duplicate_stratified_row(tmp_path: Path) -> None:
    metadata_dir = _build_characterize_output(tmp_path, _default_cohort())
    hyp_path = metadata_dir / "hypotheses.tsv"
    hyp = pd.read_csv(hyp_path, sep="\t", dtype=str, keep_default_na=False)
    strat_rows = hyp[hyp["scope"] == mf.STRATIFIED_SCOPE]
    hyp = pd.concat([hyp, strat_rows], ignore_index=True)
    hyp.to_csv(hyp_path, sep="\t", index=False)

    registry = tmp_path / "state" / "color_registry.json"
    _write_registry_for(metadata_dir, registry)
    with pytest.raises(mf.MetadataConsistencyError, match="expected exactly one"):
        mf.render_all(
            metadata_dir, registry, _dirs(tmp_path), tmp_path / "prov.json", dpi=50
        )


# --------------------------------------------------------------------------- #
# >8 categories -> CategoricalPaletteExceededError
# --------------------------------------------------------------------------- #


def test_more_than_eight_candidate_pair_levels_raises(tmp_path: Path) -> None:
    # 201 (odd) and 202 (even) are a consecutive pair -> both map to the SAME
    # candidate_pair "P201_202"; one of each condition keeps the stratified
    # p_two_sided non-blank (u_max=1 > 0), so load_tables succeeds and we reach
    # the color-assignment step this test targets.
    metadata_dir = _build_characterize_output(
        tmp_path,
        [
            (2021, 1, 1, 1, 201, "control"),
            (2021, 1, 1, 2, 202, "raloxifene-d0"),
        ],
    )
    samples = pd.read_csv(metadata_dir / "samples.tsv", sep="\t", dtype=str)
    (actual_pair,) = samples["candidate_pair"].unique()

    # Pre-seed the registry's candidate_pair category with 9 DISTINCT labels
    # (including the one this cohort actually uses) -- assign_colors' >8 guard
    # fires on the count of already-registered distinct labels alone.
    nine_pairs = {
        f"P{900 + 2 * i}_{901 + 2 * i}": _PALETTE[i % len(_PALETTE)] for i in range(9)
    }
    nine_pairs[actual_pair] = _PALETTE[0]
    registry = tmp_path / "state" / "color_registry.json"
    registry.parent.mkdir(parents=True, exist_ok=True)
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
                    "values": dict.fromkeys(samples["batch"].unique(), "#009E73"),
                },
                "candidate_pair": {"scope": "project", "values": nine_pairs},
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(CategoricalPaletteExceededError):
        mf.render_all(
            metadata_dir, registry, _dirs(tmp_path), tmp_path / "prov.json", dpi=50
        )


# --------------------------------------------------------------------------- #
# Corrupted/tampered precomputed tables (never producible by the real
# pipeline) -- direct table-writing fixtures.
# --------------------------------------------------------------------------- #

# batch, seq, position, half, sample_id, condition, pair
_ROWS = [
    ("BA", 10, 1, "early", "S1", "control", "P1"),
    ("BA", 12, 2, "early", "S3", "control", "P2"),
    ("BA", 15, 3, "late", "S2", "raloxifene-d0", "P1"),
    ("BA", 18, 4, "late", "S4", "raloxifene-d0", "P2"),
    ("BB", 3, 1, "early", "S5", "control", "P3"),
    ("BB", 7, 2, "late", "S6", "raloxifene-d0", "P3"),
]


def _write_raw_tables(
    root: Path, rows: list[tuple[str, int, int, str, str, str, str]] | None = None
) -> tuple[Path, Path]:
    """Write the precomputed tables DIRECTLY (bypassing metadata_characterize.py),
    for scenarios its own validity checks would never let through."""
    meta = root / "results" / "metadata"
    meta.mkdir(parents=True)
    samples = pd.DataFrame(
        _ROWS if rows is None else rows,
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
                "p_one_sided_observed_direction": "0.1667",
                "p_two_sided": "0.3333",
            },
            {
                "hypothesis_id": "H4",
                "scope": mf.STRATIFIED_SCOPE,
                "p_one_sided_observed_direction": "0.0833",
                "p_two_sided": "0.1667",
            },
        ]
    ).to_csv(meta / "hypotheses.tsv", sep="\t", index=False)
    (meta / "data_version.json").write_text(
        json.dumps({"data_version": "sha256:test"}), encoding="utf-8"
    )
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
        ),
        encoding="utf-8",
    )
    return meta, registry


def test_crosstab_disagreement_fails_loud(tmp_path: Path) -> None:
    meta, registry = _write_raw_tables(tmp_path)
    path = meta / "crosstab_condition_batch.tsv"
    ct = pd.read_csv(path, sep="\t", index_col="condition")
    ct.loc["control", "BA"] += 1
    ct.to_csv(path, sep="\t")
    with pytest.raises(mf.MetadataConsistencyError, match="crosstab_condition_batch"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_run_layout_disagreement_fails_loud(tmp_path: Path) -> None:
    meta, registry = _write_raw_tables(tmp_path)
    path = meta / "run_layout.tsv"
    text = path.read_text(encoding="utf-8").replace("control/S1", "control/S9")
    path.write_text(text, encoding="utf-8")
    with pytest.raises(mf.MetadataConsistencyError, match="run_layout"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_level_missing_from_registry_fails_loud(tmp_path: Path) -> None:
    meta, registry = _write_raw_tables(tmp_path)
    reg = json.loads(registry.read_text(encoding="utf-8"))
    del reg["candidate_pair"]["values"]["P3"]
    registry.write_text(json.dumps(reg), encoding="utf-8")
    with pytest.raises(KeyError, match="P3"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_duplicate_sample_id_fails_loud(tmp_path: Path) -> None:
    meta, registry = _write_raw_tables(tmp_path)
    path = meta / "samples.tsv"
    text = path.read_text(encoding="utf-8").replace("S6", "S5")
    path.write_text(text, encoding="utf-8")
    with pytest.raises(mf.MetadataConsistencyError, match="duplicate sample_id"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_unknown_condition_level_fails_loud(tmp_path: Path) -> None:
    meta, registry = _write_raw_tables(tmp_path)
    path = meta / "samples.tsv"
    text = path.read_text(encoding="utf-8").replace("control", "TREATED")
    path.write_text(text, encoding="utf-8")
    with pytest.raises(mf.MetadataConsistencyError, match="condition levels"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_unknown_run_half_level_fails_loud(tmp_path: Path) -> None:
    meta, registry = _write_raw_tables(tmp_path)
    path = meta / "samples.tsv"
    text = path.read_text(encoding="utf-8").replace("early", "middle")
    path.write_text(text, encoding="utf-8")
    with pytest.raises(mf.MetadataConsistencyError, match="run_half levels"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_empty_samples_tsv_fails_loud(tmp_path: Path) -> None:
    meta, registry = _write_raw_tables(tmp_path)
    header = (meta / "samples.tsv").read_text(encoding="utf-8").splitlines()[0]
    (meta / "samples.tsv").write_text(header + "\n", encoding="utf-8")
    with pytest.raises(mf.MetadataConsistencyError, match="no rows"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_run_half_violating_shared_rule_fails_loud(tmp_path: Path) -> None:
    # Odd batch of 3: the shared rule puts the middle run in "early"; label it
    # "late" (crosstabs/layout are derived from the same rows, so only the rule
    # check can fire).
    rows = [
        ("BA", 10, 1, "early", "S1", "control", "P1"),
        ("BA", 12, 2, "late", "S3", "control", "P2"),
        ("BA", 15, 3, "late", "S2", "raloxifene-d0", "P1"),
        ("BB", 3, 1, "early", "S5", "control", "P3"),
        ("BB", 7, 2, "late", "S6", "raloxifene-d0", "P3"),
    ]
    meta, registry = _write_raw_tables(tmp_path, rows)
    with pytest.raises(mf.MetadataConsistencyError, match="run_half"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_run_position_not_seq_rank_fails_loud(tmp_path: Path) -> None:
    rows = list(_ROWS)
    rows[0] = ("BA", 10, 2, "early", "S1", "control", "P1")
    rows[1] = ("BA", 12, 1, "early", "S3", "control", "P2")
    meta, registry = _write_raw_tables(tmp_path, rows)
    with pytest.raises(mf.MetadataConsistencyError, match="run positions"):
        mf.render_all(meta, registry, _dirs(tmp_path), tmp_path / "p.json")


def test_svgs_keep_text_as_text(tmp_path: Path) -> None:
    # Regression: save_figure must run inside publication_style() so
    # svg.fonttype="none" is in effect (text stays <text>, not glyph paths).
    metadata_dir = _build_characterize_output(tmp_path, _default_cohort())
    registry = tmp_path / "state" / "color_registry.json"
    _write_registry_for(metadata_dir, registry)
    records = mf.render_all(
        metadata_dir,
        registry,
        _dirs(tmp_path),
        tmp_path / "prov.json",
        dpi=50,
        project_root=tmp_path,
    )
    for record in records.values():
        svg = (tmp_path / str(record["svg"])).read_text(encoding="utf-8")
        assert "<text" in svg
