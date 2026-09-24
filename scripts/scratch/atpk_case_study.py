"""Stage 4: ATPK_HUMAN (P56134) case study -- why LFQ and spectral counts disagree.

ATPK is the one protein that becomes a q < 0.05 hit under limma-trend (paired LFQ,
log2FC -0.52), while its spectral-count fold changes point the other way (NSAF / log2
PSM about +0.8). This script lays out every measurement behind those numbers:

  * per-sample and per-pair (raloxifene - control, within candidate_pair) values for
    protein LFQ (raw log2, median-normalized log2, and the per-sample normalization
    offset), log2 NSAF, and PSM counts (count and log2),
  * the LFQ peptide evidence: every FlashLFQ peptide assigned to the ATPK group, its
    per-run intensity and FlashLFQ Detection Type (MSMS / MBR / ambiguous / not
    detected; via the verified ``loaders.peptide_loader``), whether it is shared, its
    normalized log2 per-pair differences and its peptide-level DE result,
  * the per-peptide PSM counts and MBR flags from the Limelight peptide dump (runs
    mapped by the Stage-2 label map ``results/stage2/limelight_label_map.tsv``),
  * a reconciliation of the FlashLFQ protein intensity with the sum of its reported
    peptide intensities, and the Poisson count-noise scale of the PSM differences.

Outputs (default ``results/de/raloxifene-vs-control/trend``): ``atpk_case.json``,
``atpk_case_samples.tsv`` (long: one row per quantity x feature x sample) and
``atpk_case_pairs.tsv`` (one row per quantity x feature x pair; the per-pair plot
input). Descriptive only: no new test is run (the DE numbers are read from the base
and trend tables).

Run:
    ./.venv/bin/python scripts/scratch/atpk_case_study.py
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import de_raloxifene_vs_control as base  # noqa: E402
from common.hashing import sha256_of_file  # noqa: E402
from loaders.dataset_io import load_dataset  # noqa: E402
from loaders.peptide_loader import load_peptide_dataset  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "de-raloxifene-vs-control-trend",
    "kind": "analysis",
    "provides": [],
    "uses": [
        "de_raloxifene_vs_control",
        "loaders.dataset_io",
        "loaders.peptide_loader",
        "common.hashing",
    ],
    "seeded_from": None,
    "description": (
        "Descriptive ATPK_HUMAN (P56134) case study: per-sample / per-pair LFQ, NSAF "
        "and PSM values, peptide-level LFQ evidence with FlashLFQ detection types, "
        "per-peptide PSMs from the Limelight dump, protein-vs-peptide reconciliation."
    ),
}

LOG = logging.getLogger("atpk_case_study")

ACCESSION = "P56134"
PROTEIN_GROUP = "psvid_86283_sp|P56134|ATPK_HUMAN"
META_COLUMNS = ("condition", "candidate_pair", "batch", "seq_number")
DUMP_PREFIXES = ("PSMs", "MBR", "Quant")


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def _one_feature(ds_path: Path, key: str, value: str) -> tuple[np.ndarray, pd.Index]:
    """The single feature of a qc state whose ``key`` equals ``value`` (fail loud)."""
    ds = load_dataset(ds_path)
    ids = ds.feature_names if key == "feature" else ds.feature_metadata[key].to_numpy()
    idx = np.flatnonzero(np.asarray(ids).astype(str) == value)
    if idx.size != 1:
        raise ValueError(f"{ds_path}: expected one {key} == {value!r}; got {idx.size}.")
    return np.asarray(ds.abundances[:, idx[0]], dtype=float), ds.metadata.index


def per_sample_offsets(qc_root: Path, level: str) -> pd.Series:
    """Per-sample median-normalization offset: normalized_log - log2(raw_linear).

    Checked constant across features (within 1e-3 log2) -- the normalization is one
    multiplicative factor per sample.
    """
    raw = load_dataset(qc_root / level / "raw_linear")
    nl = load_dataset(qc_root / level / "normalized_log")
    raw_df = pd.DataFrame(
        raw.abundances.T, index=raw.feature_names, columns=raw.metadata.index
    )
    nl_df = pd.DataFrame(
        nl.abundances.T, index=nl.feature_names, columns=nl.metadata.index
    )
    diff = nl_df - np.log2(raw_df.loc[nl_df.index, nl_df.columns])
    spread = float((diff.max() - diff.min()).max())
    if spread > 1e-3:
        raise ValueError(f"{level}: normalization offset not constant ({spread:.2e}).")
    return diff.median()


def read_peptide_dump(
    dump_file: Path, label_map_file: Path, sample_ids: pd.Index
) -> pd.DataFrame:
    """ATPK rows of the Limelight peptide dump, long format (one row per run).

    Columns: dump_sequence, base_sequence, sample_id, psms (int), mbr (bool), quant
    (the display string: a 3-significant-figure number or e.g. 'overlapping signal').
    """
    frame = pd.read_csv(
        dump_file, sep="\t", dtype=str, keep_default_na=False, quoting=3
    )
    label_map = pd.read_csv(label_map_file, sep="\t", dtype=str)
    if set(label_map["sample_id"]) != set(sample_ids.astype(str)):
        raise ValueError("Label map samples differ from the analysis samples.")
    rows = frame[frame["Protein(s)"].str.contains(ACCESSION, regex=False)]
    if rows.empty:
        raise ValueError(f"No {ACCESSION} rows in {dump_file}.")
    out: list[dict[str, object]] = []
    for _, r in rows.iterrows():
        seq = str(r["Peptide Sequence"])
        base_seq = re.sub(r"\[[^\]]*\]", "", seq)
        for _, lm in label_map.iterrows():
            label = str(lm["label"])
            cols = {p: f"{p} ({label})" for p in DUMP_PREFIXES}
            missing = [c for c in cols.values() if c not in frame.columns]
            if missing:
                raise ValueError(f"Dump lacks columns {missing}.")
            psm_text = str(r[cols["PSMs"]]).replace(",", "").strip()
            if not psm_text.isdigit():
                raise ValueError(f"Non-integer PSMs {psm_text!r} for {seq}/{label}.")
            mbr_text = str(r[cols["MBR"]]).strip()
            if mbr_text not in ("", "(MBR)"):
                raise ValueError(f"Unexpected MBR token {mbr_text!r}.")
            out.append(
                {
                    "dump_sequence": seq,
                    "base_sequence": base_seq,
                    "shared_in_dump": str(r["Unique"]).strip() != "*",
                    "sample_id": str(lm["sample_id"]),
                    "psms": int(psm_text),
                    "mbr": mbr_text == "(MBR)",
                    "quant": str(r[cols["Quant"]]).strip(),
                }
            )
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- #
# Pair arithmetic
# --------------------------------------------------------------------------- #
def pair_rows(
    quantity: str,
    feature: str,
    values: pd.Series,
    meta: pd.DataFrame,
    *,
    ratio: bool = False,
) -> list[dict[str, object]]:
    """Within-pair raloxifene - control rows for one quantity x feature.

    ``ratio=True`` (counts): the difference is ``log2(ralox / control)`` when both are
    positive, else ``NaN``. ``NaN`` values (not quantified) give a ``NaN`` difference.
    """
    rows: list[dict[str, object]] = []
    for pair, sub in meta.groupby("candidate_pair", sort=True):
        ctrl = sub.index[sub["condition"] == base.REFERENCE][0]
        ralo = sub.index[sub["condition"] == base.TREATED][0]
        c, t = float(values[ctrl]), float(values[ralo])
        if ratio:
            diff = float(np.log2(t / c)) if c > 0 and t > 0 else float("nan")
        else:
            diff = t - c
        rows.append(
            {
                "quantity": quantity,
                "feature": feature,
                "pair": str(pair),
                "batch": str(sub["batch"].iloc[0]),
                "control_sample": str(ctrl),
                "raloxifene_sample": str(ralo),
                "control_seq_number": int(meta.loc[ctrl, "seq_number"]),
                "raloxifene_seq_number": int(meta.loc[ralo, "seq_number"]),
                "control_value": c,
                "raloxifene_value": t,
                "difference": diff,
                "difference_kind": "log2 ratio" if ratio else "difference (log2)",
            }
        )
    return rows


def _summ(rows: list[dict[str, object]]) -> dict[str, object]:
    d = np.array([float(str(r["difference"])) for r in rows])
    fin = d[np.isfinite(d)]
    return {
        "per_pair": {str(r["pair"]): _r(float(str(r["difference"]))) for r in rows},
        "mean": _r(float(fin.mean())) if fin.size else None,
        "sd": _r(float(fin.std(ddof=1))) if fin.size > 1 else None,
        "n_pairs": int(fin.size),
        "n_pairs_negative": int((fin < 0).sum()),
        "n_pairs_positive": int((fin > 0).sum()),
    }


def _r(x: float, nd: int = 4) -> float | None:
    return round(x, nd) if np.isfinite(x) else None


def de_row(path: Path, feature: str, key: str = "feature") -> dict[str, object]:
    """ATPK's row of a DE table (trend tables carry the no-trend q as well).

    ``key`` is the id column to match (``first_member_id`` for NSAF / PSM tables).
    """
    t = pd.read_csv(path, sep="\t", dtype={key: str})
    row = t[t[key] == feature]
    if len(row) != 1:
        raise ValueError(f"{path}: {feature!r} not found exactly once.")
    r = row.iloc[0]
    out: dict[str, object] = {
        "log2fc": _r(float(r["log2fc"])),
        "ci95": [_r(float(r["ci_low"])), _r(float(r["ci_high"]))],
        "p": float(r["p"]),
        "q": _r(float(r["q"])),
        "residual_sd": _r(float(r["residual_sd"])),
        "rank_by_p": int(t["p"].rank(method="min")[row.index[0]]),
        "n_features": len(t),
    }
    if "notrend_q" in t.columns:
        out.update(
            {
                "p_notrend": float(r["notrend_p"]),
                "q_notrend": _r(float(r["notrend_q"])),
                "prior_sd_trend": _r(float(r["prior_sd"])),
            }
        )
    return out


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(argv: Sequence[str] | None = None) -> int:
    root = _SCRATCH.parent.parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--qc-root", type=Path, default=root / "results" / "qc_states")
    parser.add_argument("--data-dir", type=Path, default=root / "data")
    parser.add_argument(
        "--samples-file", type=Path, default=root / "results/metadata/samples.tsv"
    )
    parser.add_argument(
        "--label-map",
        type=Path,
        default=root / "results/stage2/limelight_label_map.tsv",
    )
    parser.add_argument(
        "--de-dir", type=Path, default=root / "results/de/raloxifene-vs-control"
    )
    parser.add_argument(
        "--out-dir", type=Path, default=root / "results/de/raloxifene-vs-control/trend"
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    qc: Path = args.qc_root
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((qc / "manifest.json").read_text(encoding="utf-8"))

    # ---- protein-level values ------------------------------------------------ #
    lfq_raw, samples = _one_feature(qc / "protein/raw_linear", "feature", PROTEIN_GROUP)
    lfq_norm, s2 = _one_feature(qc / "protein/normalized_log", "feature", PROTEIN_GROUP)
    nsaf_log2, s3 = _one_feature(qc / "nsaf/raw_log", "first_member_id", PROTEIN_GROUP)
    psm, s4 = _one_feature(
        qc / "psm/raw_linear_complete", "first_member_id", PROTEIN_GROUP
    )
    for other in (s2, s3, s4):
        if not other.equals(samples):
            raise ValueError("Sample order differs between qc states.")
    meta = (
        load_dataset(qc / "protein/normalized_log")
        .metadata.loc[:, list(META_COLUMNS)]
        .copy()
    )
    base.validate_design(load_dataset(qc / "protein/normalized_log").metadata)
    protein_offsets = per_sample_offsets(qc, "protein")
    peptide_offsets = per_sample_offsets(qc, "peptide")
    if not np.allclose(
        lfq_norm, np.log2(lfq_raw) + protein_offsets[samples].to_numpy(), atol=1e-4
    ):
        raise ValueError("ATPK normalized log2 != raw log2 + sample offset.")

    protein_values: dict[str, pd.Series] = {
        "lfq_protein_raw_log2": pd.Series(np.log2(lfq_raw), index=samples),
        "lfq_protein_norm_log2": pd.Series(lfq_norm, index=samples),
        "lfq_protein_norm_offset": protein_offsets[samples],
        "nsaf_log2": pd.Series(nsaf_log2, index=samples),
        "psm_count": pd.Series(psm, index=samples),
        "psm_log2": pd.Series(np.log2(psm), index=samples),
    }

    # ---- peptide-level LFQ evidence ------------------------------------------ #
    pep = load_peptide_dataset(args.data_dir / "peptide-quants.tsv", args.samples_file)
    pds = pep.dataset
    if not pds.metadata.index.astype(str).equals(samples.astype(str)):
        pep_order = [list(pds.metadata.index.astype(str)).index(s) for s in samples]
    else:
        pep_order = list(range(len(samples)))
    in_atpk = pds.feature_metadata["protein_groups"].str.contains(ACCESSION).to_numpy()
    complete_peptides = set(
        load_dataset(qc / "peptide/normalized_log").feature_names.astype(str)
    )
    pep_trend = pd.read_csv(
        args.de_dir / "trend/peptide_paired.tsv", sep="\t", dtype={"feature": str}
    )
    peptides: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    pair_table: list[dict[str, object]] = []
    peptide_intensity: dict[str, pd.Series] = {}
    for j in np.flatnonzero(in_atpk):
        fid = str(pds.feature_names[j])
        fmeta = pds.feature_metadata.iloc[j]
        inten = pd.Series(pds.abundances[pep_order, j], index=samples)
        det = pd.Series(pep.detection_type[pep_order, j], index=samples)
        peptide_intensity[fid] = inten.fillna(0.0)
        norm = np.log2(inten) + peptide_offsets[samples]
        rows = pair_rows("lfq_peptide_norm_log2", fid, norm, meta)
        pair_table.extend(rows)
        for sid in samples:
            sample_rows.append(
                {
                    "quantity": "lfq_peptide_norm_log2",
                    "feature": fid,
                    "sample_id": str(sid),
                    "value": float(norm[sid]),
                    "raw_intensity": float(inten[sid]),
                    "detection_type": str(det[sid]),
                }
            )
        de = pep_trend[pep_trend["feature"] == fid]
        peptides.append(
            {
                "peptide": fid,
                "base_sequence": str(fmeta["base_sequence"]),
                "total_mod_mass": float(fmeta["total_mod_mass"]),
                "protein_groups": str(fmeta["protein_groups"]),
                "shared": ";" in str(fmeta["protein_groups"]),
                "in_complete_de_set": fid in complete_peptides,
                "detection_type_counts": det.value_counts().to_dict(),
                "detection_type_by_sample": det.astype(str).to_dict(),
                "n_mbr": int((det == "MBR").sum()),
                "within_pair_norm_log2_difference": _summ(rows),
                "peptide_de_paired": None
                if de.empty
                else {
                    "log2fc": _r(float(de["log2fc"].iloc[0])),
                    "ci95": [
                        _r(float(de["ci_low"].iloc[0])),
                        _r(float(de["ci_high"].iloc[0])),
                    ],
                    "q_trend": _r(float(de["q"].iloc[0])),
                    "q_notrend": _r(float(de["notrend_q"].iloc[0])),
                },
            }
        )

    # ---- per-peptide PSMs from the Limelight dump ---------------------------- #
    dump = read_peptide_dump(
        args.data_dir / "peptide-limelight-table-dump.txt", args.label_map, samples
    )
    dump_summary: list[dict[str, object]] = []
    for seq, sub in dump.groupby("dump_sequence", sort=True):
        counts = sub.set_index("sample_id")["psms"].reindex(samples.astype(str))
        counts.index = samples
        rows = pair_rows("psm_count_peptide_dump", str(seq), counts, meta, ratio=True)
        pair_table.extend(rows)
        for sid in samples:
            r = sub[sub["sample_id"] == str(sid)].iloc[0]
            sample_rows.append(
                {
                    "quantity": "psm_count_peptide_dump",
                    "feature": str(seq),
                    "sample_id": str(sid),
                    "value": float(r["psms"]),
                    "mbr_flag": bool(r["mbr"]),
                    "dump_quant": str(r["quant"]),
                }
            )
        ctrl = meta.index[meta["condition"] == base.REFERENCE]
        ralo = meta.index[meta["condition"] == base.TREATED]
        dump_summary.append(
            {
                "dump_sequence": str(seq),
                "psms_by_sample": {str(k): int(v) for k, v in counts.items()},
                "psms_total_control": int(counts[ctrl].sum()),
                "psms_total_raloxifene": int(counts[ralo].sum()),
                "n_runs_mbr": int(sub["mbr"].sum()),
                "dump_quant_by_sample": dict(
                    zip(sub["sample_id"], sub["quant"], strict=True)
                ),
                "within_pair_log2_psm_ratio": _summ(rows),
            }
        )
    total_dump = dump.groupby("sample_id")["psms"].sum().reindex(samples.astype(str))
    if not np.array_equal(total_dump.to_numpy(), psm.astype(int)):
        LOG.warning(
            "Sum of peptide-dump PSMs %s != protein PSMs %s",
            total_dump.to_list(),
            psm.tolist(),
        )

    # ---- protein-level per-sample + per-pair rows ---------------------------- #
    for q, vals in protein_values.items():
        is_count = q == "psm_count"
        pair_table.extend(pair_rows(q, PROTEIN_GROUP, vals, meta, ratio=is_count))
        for sid in samples:
            sample_rows.append(
                {
                    "quantity": q,
                    "feature": PROTEIN_GROUP,
                    "sample_id": str(sid),
                    "value": float(vals[sid]),
                }
            )
    pairs_df = pd.DataFrame(pair_table)
    samples_df = pd.DataFrame(sample_rows).merge(
        meta.reset_index()
        .rename(columns={"index": "sample_id"})
        .astype({"sample_id": str}),
        on="sample_id",
        how="left",
        validate="m:1",
    )

    # ---- reconciliation + count noise ---------------------------------------- #
    pep_sum = sum(peptide_intensity.values())
    recon_ratio = pd.Series(lfq_raw, index=samples) / pep_sum
    psm_s = protein_values["psm_count"]
    poisson_sd = []
    for _, sub in meta.groupby("candidate_pair"):
        c = sub.index[sub["condition"] == base.REFERENCE][0]
        t = sub.index[sub["condition"] == base.TREATED][0]
        poisson_sd.append(float(np.sqrt(1 / psm_s[c] + 1 / psm_s[t]) / np.log(2)))

    protein_summary = {
        q: _summ([r for r in pair_table if r["quantity"] == q]) for q in protein_values
    }
    de_results = {
        "lfq_protein": {
            d: de_row(args.de_dir / f"trend/protein_{d}.tsv", PROTEIN_GROUP)
            for d in base.DESIGNS
        },
        "nsaf_paired": de_row(
            args.de_dir / "trend/nsaf_paired.tsv", PROTEIN_GROUP, "first_member_id"
        ),
        "psm_log2_paired": de_row(
            args.de_dir / "trend/psm_log2_paired.tsv", PROTEIN_GROUP, "first_member_id"
        ),
    }
    ox = [d for d in dump_summary if "[" in str(d["dump_sequence"])]
    lfq_pair = protein_summary["lfq_protein_norm_log2"]
    raw_pair = protein_summary["lfq_protein_raw_log2"]
    psm_pair = protein_summary["psm_count"]
    offs_pair = protein_summary["lfq_protein_norm_offset"]
    total_psm_c = sum(int(str(d["psms_total_control"])) for d in dump_summary)
    total_psm_r = sum(int(str(d["psms_total_raloxifene"])) for d in dump_summary)
    top = max(
        dump_summary,
        key=lambda d: (
            int(str(d["psms_total_raloxifene"])) - int(str(d["psms_total_control"]))
        ),
    )
    top_gain = int(str(top["psms_total_raloxifene"])) - int(
        str(top["psms_total_control"])
    )
    n_mbr = sum(int(str(pp["n_mbr"])) for pp in peptides) + sum(
        int(str(d["n_runs_mbr"])) for d in dump_summary
    )
    n_shared = sum(bool(pp["shared"]) for pp in peptides)
    pep_lines = [
        f"{pp['peptide']}: mean within-pair normalized log2 difference "
        f"{pp['within_pair_norm_log2_difference']['mean']} "  # type: ignore[index]
        f"(per pair {pp['within_pair_norm_log2_difference']['per_pair']}; "  # type: ignore[index]
        f"detection types {pp['detection_type_counts']})"
        for pp in peptides
    ]
    ox_lfq = [pp for pp in peptides if float(str(pp["total_mod_mass"])) != 0.0]
    ox_direction = (
        "falls"
        if ox_lfq
        and float(str(ox_lfq[0]["within_pair_norm_log2_difference"]["mean"]))  # type: ignore[index]
        < 0
        else "does not fall"
    )
    interpretation = [
        f"LFQ protein: normalized log2 within-pair differences "
        f"{lfq_pair['per_pair']} (mean {lfq_pair['mean']}, SD {lfq_pair['sd']}); "
        f"{lfq_pair['n_pairs_negative']} of {lfq_pair['n_pairs']} pairs negative and "
        f"nearly identical, so the residual SD is tiny; the trend prior (lower at high "
        f"abundance) then moderates it less towards a larger variance and it passes "
        f"q < 0.05.",
        f"Normalization does not create the LFQ decrease: raw log2 differences "
        f"{raw_pair['per_pair']} (mean {raw_pair['mean']}); the median-normalization "
        f"offsets add {offs_pair['per_pair']} per pair (raloxifene runs scaled up), "
        f"shrinking the effect and making it more uniform across pairs.",
        "LFQ peptide evidence (FlashLFQ peptide table): " + "; ".join(pep_lines) + ".",
        f"Spectral counts are small: protein PSMs per run "
        f"{ {str(k): int(v) for k, v in psm_s.items()} } (mean "
        f"{psm_s.mean():.2f}); within-pair log2 PSM ratios {psm_pair['per_pair']}. "
        f"Poisson noise alone gives a per-pair log2-ratio SD of about "
        f"{np.median(poisson_sd):.2f}.",
        "Per-peptide PSMs (Limelight dump), control vs raloxifene totals: "
        + "; ".join(
            f"{d['dump_sequence']} {d['psms_total_control']} vs "
            f"{d['psms_total_raloxifene']}"
            for d in dump_summary
        )
        + f" (all peptides {total_psm_c} vs {total_psm_r}). {top['dump_sequence']} "
        f"gains {top_gain:+d} PSMs (all other ATPK peptides together "
        f"{total_psm_r - total_psm_c - top_gain:+d}); "
        f"the modified LFQ peptide "
        f"({ox_lfq[0]['peptide'] if ox_lfq else 'none'}) MS1 intensity "
        f"{ox_direction} in the same pairs.",
        f"MBR / shared peptides: {n_mbr} MBR values among the ATPK peptides "
        f"(FlashLFQ detection types + dump MBR flags) and {n_shared} shared peptides.",
        f"The FlashLFQ protein intensity is not the sum of the reported peptide "
        f"intensities (ratio protein / peptide sum per sample "
        f"{ {str(k): round(float(v), 3) for k, v in recon_ratio.items()} }): the "
        f"protein roll-up uses signal not visible in the peptide table (e.g. peaks "
        f"reported as MSMSAmbiguousPeakfinding, intensity 0 in the peptide table).",
        "Run order is aliased with condition (control first in every pair, finding "
        "0001), so a run-order effect (e.g. Met oxidation or chromatographic drift "
        "along the queue) cannot be separated from treatment for any of these "
        "directions.",
    ]
    result = {
        "question": "Why does ATPK (P56134) go down in LFQ (log2FC -0.52, trend "
        "q 0.04) but up in NSAF / PSM counts?",
        "data_version": manifest["data_version"],
        "script": "scripts/scratch/atpk_case_study.py",
        "script_sha256": sha256_of_file(Path(__file__)),
        "protein_group": PROTEIN_GROUP,
        "gene": "ATP5MF (ATP synthase F(0) complex subunit f, mitochondrial)",
        "difference_convention": "raloxifene-d0 minus control within candidate_pair "
        "(log2 units; counts as log2 ratio)",
        "samples": {
            str(s): {c: str(meta.loc[s, c]) for c in META_COLUMNS} for s in samples
        },
        "protein_per_sample": {
            q: {str(k): _r(float(v), 5) for k, v in vals.items()}
            for q, vals in protein_values.items()
        },
        "protein_within_pair": protein_summary,
        "mean_psms_per_run": round(float(psm_s.mean()), 3),
        "poisson_sd_log2_psm_ratio_per_pair": [round(v, 3) for v in poisson_sd],
        "de_results": de_results,
        "lfq_peptides": peptides,
        "peptide_dump_psms": dump_summary,
        "oxidized_peptide_dump_rows": [str(d["dump_sequence"]) for d in ox],
        "protein_vs_peptide_sum": {
            "peptide_intensity_sum_by_sample": {
                str(k): float(v) for k, v in pep_sum.items()
            },
            "ratio_protein_over_peptide_sum": {
                str(k): round(float(v), 4) for k, v in recon_ratio.items()
            },
            "note": "Across all proteins only ~11% of protein x run cells equal the "
            "sum of their unique reported peptide intensities (checked 2026-09-24); "
            "the FlashLFQ protein roll-up is not reconstructable from peptide-quants.",
        },
        "interpretation": interpretation,
        "caveats": [
            "Descriptive; no test beyond the existing DE tables. One protein chosen "
            "because it is the top hit: selection makes its effect size optimistic.",
            "Dump numbers are display-rounded (3 s.f.); PSM counts are exact integers.",
            "Four pairs; the 2022 pair (P941_942) is a different acquisition batch.",
        ],
        "outputs": {
            "samples_tsv": "atpk_case_samples.tsv",
            "pairs_tsv": "atpk_case_pairs.tsv",
        },
    }
    (out_dir / "atpk_case.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    samples_df.to_csv(out_dir / "atpk_case_samples.tsv", sep="\t", index=False)
    pairs_df.to_csv(out_dir / "atpk_case_pairs.tsv", sep="\t", index=False)
    LOG.info("wrote %s", out_dir / "atpk_case.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
