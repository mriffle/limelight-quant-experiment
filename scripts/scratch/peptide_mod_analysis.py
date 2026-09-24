"""Stage 4: peptide-level modification analysis, raloxifene-d0 vs control.

Six analyses on the FlashLFQ peptide LFQ data (``state/DATA_DESCRIPTION.md``), all on
the 8 experimental runs (no QC/pool controls exist), with the 4 matched
control/raloxifene pairs (finding 0003) as the unit of replication:

1. **Modification decomposition** of every peptide's single appended total mod mass
   into carbamidomethyl (57.021464), oxidation (15.994915) and the raloxifene adduct
   (C28H27NO4S - 2H, derived in ``analysis.peptide_mods``); classes unmodified /
   CAM-only / oxidized / raloxifene-adduct / unexplained; the Limelight peptide dump's
   localized sequences joined on (base, composition) give sites.
2. **Raloxifene-adduct peptides**: per-sample detection type, intensity, dump PSMs /
   MBR flag; per-pair treated-vs-control detection; exact binomial (McNemar-style) on
   discordant pair x peptide cells; adducted proteins (CYPs?).
3. **Presence/absence for all peptides** (complete-case DE cannot see
   treatment-specific peptides): per-peptide pair discordance under three detection
   definitions (quantified = MSMS|MBR; msms = MSMS only; psm = dump PSMs > 0); per
   class, the treated-only share of discordant cells vs the unmodified baseline (Fisher
   exact, primary; pair-stratified Mantel-Haenszel and a peptide-level collapse as
   sensitivities; BH over classes), baseline per pair / per batch (run-order check).
4. **Intensity by class**: paired log2FC (existing ``peptide_paired.tsv``) per class,
   class vs unmodified (Welch mean difference + CI; Mann-Whitney with Hodges-Lehmann
   shift + CI; BH over classes); top peptides by q and by |log2FC|.
5. **Protein-adjusted change**: per pair (peptide log2FC - its protein's log2FC), a
   moderated one-sample t over pairs (BH over peptides), and modified - reference-form
   log2FC (occupancy-like); class tests (BH over classes).
6. **Oxidation index** per protein and globally (sum oxidized / sum all forms of the
   base sequences that have an oxidized form), paired log2 ratio.

Inputs are read-only: ``results/qc_states/peptide/{raw_linear,raw_linear_complete,
normalized_log}``, ``results/qc_states/protein/normalized_log``,
``results/de/raloxifene-vs-control/peptide_{paired,batch,unadjusted}.tsv``,
``data/peptide-quants.tsv`` (Detection Type, via ``loaders.peptide_loader``),
``data/peptide-limelight-table-dump.txt`` and
``results/stage2/limelight_label_map.tsv``.
Outputs (TSV, plot-ready) + ``summary.json`` under ``--out-dir``
(default ``results/peptide-mods``). Deterministic; no seed consumed. No figures.

Run:
    ./.venv/bin/python scripts/scratch/peptide_mod_analysis.py
"""

from __future__ import annotations

import argparse
import json
import logging
import platform
import re
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis import peptide_mods as pm  # noqa: E402
from common.hashing import sha256_of_file  # noqa: E402
from loaders.data_loading import Dataset  # noqa: E402
from loaders.dataset_io import load_dataset  # noqa: E402
from loaders.peptide_loader import load_peptide_dataset  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "peptide-mod-analysis",
    "kind": "analysis",
    "provides": [],
    "uses": [
        "analysis.peptide_mods",
        "loaders.dataset_io",
        "loaders.peptide_loader",
        "common.hashing",
    ],
    "seeded_from": None,
    "description": (
        "Peptide-level modification analysis raloxifene-d0 vs control: mass "
        "decomposition (CAM/Ox/raloxifene adduct), dump localization join, adduct "
        "peptide detection by pair (exact binomial), presence/absence discordance by "
        "class vs unmodified baseline (Fisher/MH, BH), paired log2FC by class "
        "(Welch/Mann-Whitney, BH), protein-adjusted and modified-minus-reference "
        "change (moderated one-sample t, BH), oxidation index (paired)."
    ),
}

LOG = logging.getLogger("peptide_mod_analysis")

CONTROL = "control"
TREATED = "raloxifene-d0"
PAIR_COL = "candidate_pair"
CONDITION_COL = "condition"
POSITIVE_TYPES = ("MSMS", "MBR")
MS2_ID_TYPES = ("MSMS", "MSMSIdentifiedButNotQuantified", "MSMSAmbiguousPeakfinding")
DETECTION_DEFINITIONS: tuple[str, ...] = ("quantified", "msms", "psm")
DEFINITION_TEXT: dict[str, str] = {
    "quantified": "FlashLFQ Detection Type MSMS or MBR (intensity > 0)",
    "msms": "FlashLFQ Detection Type MSMS only (MBR transfers excluded)",
    "psm": "Limelight dump PSMs > 0 in the run (MS2 identification, quant or not)",
}
NON_BASELINE_CLASSES: tuple[str, ...] = (
    pm.CAM_ONLY,
    pm.OXIDIZED,
    pm.ADDUCT,
    pm.UNEXPLAINED,
)
CYP_ENTRY_RE = re.compile(r"^CP\d")
CYP3A4_ACCESSION = "P08684"


@dataclass(frozen=True)
class Inputs:
    """Everything loaded and cross-checked, samples in one shared order."""

    raw: Dataset
    raw_complete: Dataset
    norm_log: Dataset
    protein_log: Dataset
    detection_type: np.ndarray
    dump: pm.PeptideDump
    de: dict[str, pd.DataFrame]
    data_version: str
    input_hashes: dict[str, str]


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def _sample_ids(ds: Dataset) -> list[str]:
    return [str(s) for s in ds.metadata.index]


def load_inputs(root: Path) -> Inputs:
    """Load and cross-check every input; fail loud on any disagreement."""
    qc = root / "results" / "qc_states"
    manifest = json.loads((qc / "manifest.json").read_text(encoding="utf-8"))
    raw = load_dataset(qc / "peptide" / "raw_linear")
    raw_complete = load_dataset(qc / "peptide" / "raw_linear_complete")
    norm_log = load_dataset(qc / "peptide" / "normalized_log")
    protein_log = load_dataset(qc / "protein" / "normalized_log")
    for name, ds, scale in (
        ("peptide raw_linear", raw, "linear"),
        ("peptide raw_linear_complete", raw_complete, "linear"),
        ("peptide normalized_log", norm_log, "log2"),
        ("protein normalized_log", protein_log, "log2"),
    ):
        if ds.scale != scale:
            raise ValueError(f"{name} has scale {ds.scale!r}; expected {scale!r}.")
        if _sample_ids(ds) != _sample_ids(raw):
            raise ValueError(f"{name} sample order differs from peptide raw_linear.")

    peptide_file = root / "data" / "peptide-quants.tsv"
    loaded = load_peptide_dataset(
        peptide_file, root / "results" / "metadata" / "samples.tsv"
    )
    if _sample_ids(loaded.dataset) != _sample_ids(raw):
        raise ValueError("Loader sample order differs from the saved raw_linear state.")
    if list(loaded.dataset.feature_names) != list(raw.feature_names):
        raise ValueError(
            "Loader feature order differs from the saved raw_linear state."
        )
    if not np.array_equal(
        loaded.dataset.abundances, raw.abundances, equal_nan=True
    ):  # same file, same parser: must be bit-identical
        raise ValueError("Loader abundances differ from the saved raw_linear state.")
    LOG.info(
        "Detection Type matrix %s verified against raw_linear state",
        loaded.detection_type.shape,
    )

    label_map = pd.read_csv(
        root / "results" / "stage2" / "limelight_label_map.tsv", sep="\t", dtype=str
    )
    label_to_sample = dict(zip(label_map["label"], label_map["sample_id"], strict=True))
    dump_file = root / "data" / "peptide-limelight-table-dump.txt"
    dump = pm.read_peptide_dump(dump_file, label_to_sample, _sample_ids(raw))

    de_dir = root / "results" / "de" / "raloxifene-vs-control"
    de = {
        design: pd.read_csv(de_dir / f"peptide_{design}.tsv", sep="\t")
        for design in ("paired", "batch", "unadjusted")
    }
    for design, table in de.items():
        if table["feature"].duplicated().any():
            raise ValueError(f"peptide_{design}.tsv has duplicate features.")
        if set(table["feature"]) != set(norm_log.feature_names):
            raise ValueError(
                f"peptide_{design}.tsv features differ from the normalized_log state."
            )

    hashes = {
        "peptide_quants": sha256_of_file(peptide_file),
        "peptide_dump": sha256_of_file(dump_file),
        **{f"peptide_{d}.tsv": sha256_of_file(de_dir / f"peptide_{d}.tsv") for d in de},
    }
    return Inputs(
        raw=raw,
        raw_complete=raw_complete,
        norm_log=norm_log,
        protein_log=protein_log,
        detection_type=loaded.detection_type,
        dump=dump,
        de=de,
        data_version=str(manifest["data_version"]),
        input_hashes=hashes,
    )


# --------------------------------------------------------------------------- #
# 1. Decomposition + join
# --------------------------------------------------------------------------- #
def build_feature_table(raw: Dataset, tolerance: float) -> pd.DataFrame:
    """One row per quant peptide: ids, protein members, decomposition, class."""
    fm = raw.feature_metadata.reset_index(drop=True)
    dec = pm.decompose_feature_masses(
        fm["total_mod_mass"].to_numpy(dtype=float).tolist(),
        fm["base_sequence"].astype(str).tolist(),
        tolerance=tolerance,
    )
    members = [pm.parse_protein_members(str(g)) for g in fm["protein_groups"]]
    table = pd.DataFrame(
        {
            "feature": np.asarray(raw.feature_names, dtype=str),
            "base_sequence": fm["base_sequence"].astype(str).to_numpy(),
            "protein_groups": fm["protein_groups"].astype(str).to_numpy(),
            "accessions": [";".join(m[1] for m in ms) for ms in members],
            "entries": [";".join(m[2] for m in ms) for ms in members],
            "n_protein_members": [len(ms) for ms in members],
            "is_contaminant": fm["is_contaminant"].to_numpy(dtype=bool),
        }
    )
    table = pd.concat([table, dec], axis=1)
    keys = []
    for base, row in zip(table["base_sequence"], dec.itertuples(), strict=True):
        if row.mod_class == pm.UNEXPLAINED:
            keys.append(f"{base}|unexplained|{row.total_mod_mass!r}")
        else:
            keys.append(
                pm.composition_key(
                    base,
                    pm.ModComposition(int(row.n_cam), int(row.n_ox), int(row.n_ralox)),
                )
            )
    table["key"] = keys
    return table


def mass_decomposition_table(features: pd.DataFrame) -> pd.DataFrame:
    """Per distinct total mass: composition, class, error, feature count."""
    grouped = (
        features.groupby("total_mod_mass", sort=True)
        .agg(
            composition=("composition", "first"),
            mod_class=("mod_class", "first"),
            mass_error=("mass_error", "first"),
            n_solutions=("n_solutions", "first"),
            n_features=("feature", "size"),
            n_non_contaminant=("is_contaminant", lambda s: int((~s).sum())),
            n_cam_exceeds_cys=("cam_exceeds_cys", "sum"),
        )
        .reset_index()
    )
    return grouped


# --------------------------------------------------------------------------- #
# Detection
# --------------------------------------------------------------------------- #
def detection_masks(
    detection_type: np.ndarray, dump_psms: np.ndarray
) -> dict[str, np.ndarray]:
    """``(n_samples, n_features)`` bool masks for each detection definition."""
    psm_t = np.nan_to_num(dump_psms.T, nan=0.0)
    return {
        "quantified": np.isin(detection_type, POSITIVE_TYPES),
        "msms": detection_type == "MSMS",
        "psm": psm_t > 0,
    }


def dump_crosschecks(
    detection_type: np.ndarray, dump_psms: np.ndarray, dump_mbr: np.ndarray
) -> dict[str, object]:
    """Agreement of dump PSMs / MBR flag with the FlashLFQ Detection Type per cell."""
    has_dump = np.isfinite(dump_psms).T  # (n_samples, n_features)
    det = detection_type[has_dump]
    psm_pos = (np.nan_to_num(dump_psms.T) > 0)[has_dump]
    mbr = dump_mbr.T[has_dump]
    out: dict[str, object] = {"n_cells_with_dump_row": int(has_dump.sum())}
    tab_psm = pd.crosstab(
        pd.Series(det, name="detection_type"), pd.Series(psm_pos, name="dump_psms_gt0")
    )
    tab_mbr = pd.crosstab(
        pd.Series(det, name="detection_type"), pd.Series(mbr, name="dump_mbr_flag")
    )
    out["detection_type_x_dump_psms_gt0"] = {
        str(k): {str(c): int(v) for c, v in row.items()}
        for k, row in tab_psm.iterrows()
    }
    out["detection_type_x_dump_mbr_flag"] = {
        str(k): {str(c): int(v) for c, v in row.items()}
        for k, row in tab_mbr.iterrows()
    }
    is_mbr_type = det == "MBR"
    out["mbr_flag_agreement"] = float(np.mean(is_mbr_type == mbr))
    return out


# --------------------------------------------------------------------------- #
# 2. Adduct peptides
# --------------------------------------------------------------------------- #
def _reference_index(features: pd.DataFrame) -> np.ndarray:
    """Index of each modified form's reference form (``-1`` if none / unmodified).

    Reference = same base sequence, no Ox, no adduct, and CAM count = own CAM count +
    adduct sites on Cys (the Cys would be CAM'd when not adducted) for oxidized /
    adduct forms; the fully unmodified form for CAM-only forms. When the dump
    isoforms disagree on the adducted residue the reference is undefined (``-1``).
    """
    lookup = {k: i for i, k in enumerate(features["key"])}
    ref = np.full(len(features), -1, dtype=np.int64)
    for i, row in enumerate(features.itertuples()):
        cls = row.mod_class
        if cls in (pm.UNMODIFIED, pm.UNEXPLAINED):
            continue
        n_cam, n_ralox = int(row.n_cam), int(row.n_ralox)
        if cls == pm.CAM_ONLY:
            ref_cam = 0
        elif n_ralox == 0:
            ref_cam = n_cam
        else:
            residues = str(row.dump_ralox_residues)
            if residues == "C":
                ref_cam = n_cam + n_ralox
            elif residues and "C" not in residues:
                ref_cam = n_cam
            else:
                continue
        key = pm.composition_key(
            str(row.base_sequence), pm.ModComposition(ref_cam, 0, 0)
        )
        ref[i] = lookup.get(key, -1)
    return ref


def adduct_tables(
    features: pd.DataFrame,
    inputs: Inputs,
    masks: dict[str, np.ndarray],
    dump_psms: np.ndarray,
    dump_mbr: np.ndarray,
    layout: pm.PairLayout,
    ref_idx: np.ndarray,
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Wide + long adduct tables, per-pair summary, protein list, and tests."""
    md = inputs.raw.metadata
    samples = _sample_ids(inputs.raw)
    idx = np.flatnonzero(features["n_ralox"].fillna(0).to_numpy(dtype=int) > 0)
    inten = inputs.raw.abundances
    det = inputs.detection_type

    wide = features.iloc[idx][
        [
            "feature",
            "base_sequence",
            "accessions",
            "entries",
            "protein_groups",
            "is_contaminant",
            "composition",
            "total_mod_mass",
            "n_cam",
            "n_ox",
            "n_ralox",
            "dump_localized_sequences",
            "dump_ralox_sites",
            "dump_ralox_residues",
            "dump_ox_sites",
            "dump_unique",
        ]
    ].reset_index(drop=True)
    for defn in DETECTION_DEFINITIONS:
        for col in ("n_treated_only", "n_control_only", "n_both", "n_neither"):
            wide[f"{defn}_{col}"] = features.iloc[idx][f"{defn}_{col}"].to_numpy()
    wide["n_samples_quantified"] = masks["quantified"][:, idx].sum(axis=0)
    wide["n_control_quantified"] = masks["quantified"][layout.control_idx][:, idx].sum(
        axis=0
    )
    wide["n_treated_quantified"] = masks["quantified"][layout.treated_idx][:, idx].sum(
        axis=0
    )
    wide["n_control_msms"] = masks["msms"][layout.control_idx][:, idx].sum(axis=0)
    wide["n_treated_msms"] = masks["msms"][layout.treated_idx][:, idx].sum(axis=0)
    wide["n_control_psm"] = masks["psm"][layout.control_idx][:, idx].sum(axis=0)
    wide["n_treated_psm"] = masks["psm"][layout.treated_idx][:, idx].sum(axis=0)
    wide["control_psms_total"] = np.nansum(
        dump_psms[idx][:, layout.control_idx], axis=1
    )
    wide["treated_psms_total"] = np.nansum(
        dump_psms[idx][:, layout.treated_idx], axis=1
    )
    wide["reference_feature"] = [
        str(features["feature"].iloc[r]) if r >= 0 else "" for r in ref_idx[idx]
    ]
    for j, s in enumerate(samples):
        wide[f"{s}_detection_type"] = det[j, idx]
        wide[f"{s}_intensity"] = inten[j, idx]
        wide[f"{s}_dump_psms"] = dump_psms[idx, j]
        wide[f"{s}_dump_mbr"] = dump_mbr[idx, j]
    # occupancy-like fraction in treated runs: adduct / (adduct + reference form)
    occ = np.full((len(idx), len(layout.treated_idx)), np.nan)
    for a, (f, r) in enumerate(zip(idx, ref_idx[idx], strict=True)):
        if r < 0:
            continue
        ia = inten[layout.treated_idx, f]
        ir = inten[layout.treated_idx, r]
        both = np.isfinite(ia) & np.isfinite(ir)
        occ[a, both] = ia[both] / (ia[both] + ir[both])
    for k, pair in enumerate(layout.pair_ids):
        wide[f"{pair}_treated_adduct_fraction"] = occ[:, k]
    has_occ = np.isfinite(occ).any(axis=1)
    med_occ = np.full(len(idx), np.nan)
    if has_occ.any():
        med_occ[has_occ] = np.nanmedian(occ[has_occ], axis=1)
    wide["median_treated_adduct_fraction"] = med_occ

    # long table (plot-ready)
    long_rows = []
    for a, f in enumerate(idx):
        for j, s in enumerate(samples):
            long_rows.append(
                {
                    "feature": features["feature"].iloc[f],
                    "entries": features["entries"].iloc[f],
                    "dump_ralox_sites": features["dump_ralox_sites"].iloc[f],
                    "sample_id": s,
                    "condition": md[CONDITION_COL].iloc[j],
                    "candidate_pair": md[PAIR_COL].iloc[j],
                    "batch": md["batch"].iloc[j],
                    "run_position_within_batch": int(
                        md["run_position_within_batch"].iloc[j]
                    ),
                    "detection_type": det[j, f],
                    "intensity": inten[j, f],
                    "log2_intensity": float(np.log2(inten[j, f]))
                    if np.isfinite(inten[j, f])
                    else np.nan,
                    "dump_psms": dump_psms[f, j],
                    "dump_mbr": bool(dump_mbr[f, j]),
                    "adduct_index": a,
                }
            )
    long = pd.DataFrame(long_rows)

    # per pair summary + class-level exact binomial (McNemar-style)
    pair_rows = []
    tests: dict[str, object] = {}
    for defn in DETECTION_DEFINITIONS:
        t_only, c_only, both, neither = pm.pair_states(masks[defn][:, idx], layout)
        for k, pair in enumerate(layout.pair_ids):
            c_sample = samples[int(layout.control_idx[k])]
            pair_rows.append(
                {
                    "definition": defn,
                    "candidate_pair": pair,
                    "batch": md["batch"].iloc[int(layout.control_idx[k])],
                    "control_sample": c_sample,
                    "treated_sample": samples[int(layout.treated_idx[k])],
                    "n_adduct_features": len(idx),
                    "n_control_detected": int((c_only[k] | both[k]).sum()),
                    "n_treated_detected": int((t_only[k] | both[k]).sum()),
                    "n_treated_only": int(t_only[k].sum()),
                    "n_control_only": int(c_only[k].sum()),
                    "n_both": int(both[k].sum()),
                    "n_neither": int(neither[k].sum()),
                }
            )
        k_t, k_c = int(t_only.sum()), int(c_only.sum())
        binom = pm.binomial_proportion(k_t, k_t + k_c, alpha=alpha)
        n_ctrl_det = [
            int((c_only[k] | both[k]).sum()) for k in range(len(layout.pair_ids))
        ]
        n_trt_det = [
            int((t_only[k] | both[k]).sum()) for k in range(len(layout.pair_ids))
        ]
        n_pos = sum(t > c for t, c in zip(n_trt_det, n_ctrl_det, strict=True))
        n_neg = sum(t < c for t, c in zip(n_trt_det, n_ctrl_det, strict=True))
        residues = features["dump_ralox_residues"].to_numpy(dtype=str)[idx]
        is_cys = residues == "C"
        by_residue: dict[str, object] = {}
        for label, sel in (("Cys", is_cys), ("non-Cys (Tyr/Trp)", ~is_cys)):
            kt, kc = int(t_only[:, sel].sum()), int(c_only[:, sel].sum())
            by_residue[label] = {
                "n_features": int(sel.sum()),
                "treated_only_cells": kt,
                "control_only_cells": kc,
                "both_cells": int(both[:, sel].sum()),
                "n_features_detected_any_control": int(
                    (c_only[:, sel] | both[:, sel]).any(axis=0).sum()
                ),
                **pm.binomial_proportion(kt, kt + kc, alpha=alpha),
            }
        fisher_res = pm.fisher_odds_ratio(
            np.array(
                [
                    [int(t_only[:, is_cys].sum()), int(c_only[:, is_cys].sum())],
                    [int(t_only[:, ~is_cys].sum()), int(c_only[:, ~is_cys].sum())],
                ]
            ),
            alpha=alpha,
        )
        by_residue["Cys_vs_nonCys_fisher"] = fisher_res
        ctrl_any = (c_only | both).any(axis=0)
        by_residue["Cys_vs_nonCys_features_with_any_control_detection_fisher"] = (
            pm.fisher_odds_ratio(
                np.array(
                    [
                        [
                            int((ctrl_any & is_cys).sum()),
                            int((~ctrl_any & is_cys).sum()),
                        ],
                        [
                            int((ctrl_any & ~is_cys).sum()),
                            int((~ctrl_any & ~is_cys).sum()),
                        ],
                    ]
                ),
                alpha=alpha,
            )
        )
        tests[defn] = {
            "definition": DEFINITION_TEXT[defn],
            "by_residue": by_residue,
            "n_adduct_features": len(idx),
            "n_cells": len(idx) * len(layout.pair_ids),
            "treated_only_cells": k_t,
            "control_only_cells": k_c,
            "both_cells": int(both.sum()),
            "neither_cells": int(neither.sum()),
            "prop_treated_only_among_discordant": binom["prop"],
            "ci_low": binom["ci_low"],
            "ci_high": binom["ci_high"],
            "p_exact_binomial_vs_0.5": binom["p"],
            "test": "exact two-sided binomial on discordant cells (McNemar-style); "
            "Clopper-Pearson 95% CI",
            "n_features_detected_any_control": int(
                masks[defn][layout.control_idx][:, idx].any(axis=0).sum()
            ),
            "n_features_detected_any_treated": int(
                masks[defn][layout.treated_idx][:, idx].any(axis=0).sum()
            ),
            "pair_level_detected_counts_treated": n_trt_det,
            "pair_level_detected_counts_control": n_ctrl_det,
            "pair_level_sign_test_p": pm.sign_test(n_pos, n_neg),
            "pair_level_note": "4 pairs: min attainable two-sided sign-test p = 0.125",
        }
    pair_summary = pd.DataFrame(pair_rows)

    # control detections: how were they made?
    ctrl_rows = []
    for f in idx:
        for k in range(len(layout.pair_ids)):
            cj, tj = int(layout.control_idx[k]), int(layout.treated_idx[k])
            if det[cj, f] in POSITIVE_TYPES or dump_psms[f, cj] > 0:
                ctrl_rows.append(
                    {
                        "feature": features["feature"].iloc[f],
                        "candidate_pair": layout.pair_ids[k],
                        "control_detection_type": det[cj, f],
                        "control_dump_psms": dump_psms[f, cj],
                        "treated_detection_type": det[tj, f],
                        "treated_dump_psms": dump_psms[f, tj],
                        "log2_control_over_treated": float(
                            np.log2(inten[cj, f] / inten[tj, f])
                        )
                        if np.isfinite(inten[cj, f]) and np.isfinite(inten[tj, f])
                        else np.nan,
                    }
                )
    ctrl_det = pd.DataFrame(ctrl_rows)
    if len(ctrl_det):
        tests["control_detections"] = {
            "n_control_cells_detected_or_psm": len(ctrl_det),
            "by_control_detection_type": {
                str(k): int(v)
                for k, v in ctrl_det["control_detection_type"].value_counts().items()
            },
            "n_with_control_psm": int((ctrl_det["control_dump_psms"] > 0).sum()),
            "median_log2_control_over_treated": float(
                np.nanmedian(ctrl_det["log2_control_over_treated"])
            )
            if ctrl_det["log2_control_over_treated"].notna().any()
            else None,
        }

    # proteins
    prot_rows = []
    for f in idx:
        for member, acc, entry in pm.parse_protein_members(
            str(features["protein_groups"].iloc[f])
        ):
            prot_rows.append(
                {
                    "protein_member": member,
                    "accession": acc,
                    "entry": entry,
                    "feature": features["feature"].iloc[f],
                    "ralox_sites": features["dump_ralox_sites"].iloc[f],
                    "ralox_residues": features["dump_ralox_residues"].iloc[f],
                    "shared_peptide": int(features["n_protein_members"].iloc[f]) > 1,
                    "n_treated_quantified": int(
                        masks["quantified"][layout.treated_idx, f].sum()
                    ),
                    "n_control_quantified": int(
                        masks["quantified"][layout.control_idx, f].sum()
                    ),
                    "n_treated_msms": int(masks["msms"][layout.treated_idx, f].sum()),
                    "n_control_msms": int(masks["msms"][layout.control_idx, f].sum()),
                }
            )
    prot_long = pd.DataFrame(prot_rows)
    proteins = (
        prot_long.groupby(["accession", "entry"], sort=True)
        .agg(
            n_adduct_features=("feature", "nunique"),
            features=("feature", lambda s: ";".join(sorted(set(s)))),
            ralox_sites=("ralox_sites", lambda s: " | ".join(sorted(set(s)))),
            ralox_residues=(
                "ralox_residues",
                lambda s: "".join(sorted(set("".join(s)))),
            ),
            any_shared_peptide=("shared_peptide", "any"),
            max_treated_quantified=("n_treated_quantified", "max"),
            max_control_quantified=("n_control_quantified", "max"),
            max_treated_msms=("n_treated_msms", "max"),
            max_control_msms=("n_control_msms", "max"),
        )
        .reset_index()
    )
    proteins["is_cyp"] = proteins["entry"].str.match(CYP_ENTRY_RE.pattern)
    proteins = proteins.sort_values(
        ["is_cyp", "n_adduct_features", "entry"], ascending=[False, False, True]
    )
    return wide, long, pair_summary, proteins, {"tests": tests, "ctrl": ctrl_det}


# --------------------------------------------------------------------------- #
# 3. Presence / absence by class
# --------------------------------------------------------------------------- #
def discordance_class_tests(
    features: pd.DataFrame,
    masks: dict[str, np.ndarray],
    layout: pm.PairLayout,
    batch_of_pair: dict[str, str],
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Per definition x class: discordant-cell asymmetry vs 0.5 and vs unmodified."""
    keep = ~features["is_contaminant"].to_numpy(dtype=bool)
    classes = features["mod_class"].to_numpy(dtype=str)
    rows = []
    pair_rows = []
    batch_tests: dict[str, object] = {}
    for defn in DETECTION_DEFINITIONS:
        t_only, c_only, _, _ = pm.pair_states(masks[defn], layout)
        per_class: dict[str, dict[str, object]] = {}
        for cls in pm.MOD_CLASSES:
            sel = keep & (classes == cls)
            tt = t_only[:, sel]
            cc = c_only[:, sel]
            k_t, k_c = int(tt.sum()), int(cc.sum())
            binom = pm.binomial_proportion(k_t, k_t + k_c, alpha=alpha)
            net = tt.sum(axis=0).astype(int) - cc.sum(axis=0).astype(int)
            per_class[cls] = {
                "tt": tt,
                "cc": cc,
                "n_pep_lean_t": int((net > 0).sum()),
                "n_pep_lean_c": int((net < 0).sum()),
            }
            detected_any = masks[defn][:, sel].any(axis=0)
            rows.append(
                {
                    "definition": defn,
                    "mod_class": cls,
                    "n_peptides": int(sel.sum()),
                    "n_peptides_detected_any": int(detected_any.sum()),
                    "treated_only_cells": k_t,
                    "control_only_cells": k_c,
                    "prop_treated_only": binom["prop"],
                    "prop_ci_low": binom["ci_low"],
                    "prop_ci_high": binom["ci_high"],
                    "p_binomial_vs_0.5": binom["p"],
                    "n_peptides_lean_treated": per_class[cls]["n_pep_lean_t"],
                    "n_peptides_lean_control": per_class[cls]["n_pep_lean_c"],
                    "n_peptides_treated_only_4of4": int(
                        (tt.sum(axis=0) == len(layout.pair_ids)).sum()
                    ),
                    "n_peptides_control_only_4of4": int(
                        (cc.sum(axis=0) == len(layout.pair_ids)).sum()
                    ),
                    "n_peptides_treated_only_ge3_none_opposite": int(
                        ((tt.sum(axis=0) >= 3) & (cc.sum(axis=0) == 0)).sum()
                    ),
                    "n_peptides_control_only_ge3_none_opposite": int(
                        ((cc.sum(axis=0) >= 3) & (tt.sum(axis=0) == 0)).sum()
                    ),
                }
            )
            for k, pair in enumerate(layout.pair_ids):
                b = pm.binomial_proportion(
                    int(tt[k].sum()), int(tt[k].sum() + cc[k].sum()), alpha=alpha
                )
                pair_rows.append(
                    {
                        "definition": defn,
                        "mod_class": cls,
                        "candidate_pair": pair,
                        "batch": batch_of_pair[pair],
                        "treated_only_cells": int(tt[k].sum()),
                        "control_only_cells": int(cc[k].sum()),
                        "prop_treated_only": b["prop"],
                        "prop_ci_low": b["ci_low"],
                        "prop_ci_high": b["ci_high"],
                        "n_detected_control": int(
                            masks[defn][layout.control_idx[k], sel].sum()
                        ),
                        "n_detected_treated": int(
                            masks[defn][layout.treated_idx[k], sel].sum()
                        ),
                    }
                )
        base = per_class[pm.UNMODIFIED]
        bt, bc = base["tt"], base["cc"]
        assert isinstance(bt, np.ndarray)
        assert isinstance(bc, np.ndarray)
        for row in rows:
            if row["definition"] != defn or row["mod_class"] == pm.UNMODIFIED:
                continue
            cls = str(row["mod_class"])
            ct, cc_ = per_class[cls]["tt"], per_class[cls]["cc"]
            assert isinstance(ct, np.ndarray)
            assert isinstance(cc_, np.ndarray)
            fisher = pm.fisher_odds_ratio(
                np.array(
                    [[int(ct.sum()), int(cc_.sum())], [int(bt.sum()), int(bc.sum())]]
                ),
                alpha=alpha,
            )
            strata = np.array(
                [
                    [
                        [int(ct[k].sum()), int(cc_[k].sum())],
                        [int(bt[k].sum()), int(bc[k].sum())],
                    ]
                    for k in range(len(layout.pair_ids))
                ]
            )
            mh = pm.mantel_haenszel(strata, alpha=alpha)
            log_or = pm.pair_level_log_odds_ratios(strata)
            ok = np.isfinite(log_or)
            pl = (
                pm.paired_t(log_or[ok], alpha=alpha)
                if ok.sum() >= 2
                else {"mean": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
            )
            pep = pm.fisher_odds_ratio(
                np.array(
                    [
                        [
                            per_class[cls]["n_pep_lean_t"],
                            per_class[cls]["n_pep_lean_c"],
                        ],
                        [base["n_pep_lean_t"], base["n_pep_lean_c"]],
                    ]
                ),
                alpha=alpha,
            )
            row.update(
                fisher_or_vs_unmodified=fisher["odds_ratio"],
                fisher_ci_low=fisher["ci_low"],
                fisher_ci_high=fisher["ci_high"],
                fisher_p=fisher["p"],
                mh_or_vs_unmodified=mh["odds_ratio"],
                mh_ci_low=mh["ci_low"],
                mh_ci_high=mh["ci_high"],
                mh_p=mh["p"],
                peptide_level_or=pep["odds_ratio"],
                peptide_level_ci_low=pep["ci_low"],
                peptide_level_ci_high=pep["ci_high"],
                peptide_level_p=pep["p"],
                pair_level_or=float(np.exp(pl["mean"])),
                pair_level_ci_low=float(np.exp(pl["ci_low"])),
                pair_level_ci_high=float(np.exp(pl["ci_high"])),
                pair_level_p=pl["p"],
                pair_level_n_pairs_or_gt1=int((log_or[ok] > 0).sum()),
                pair_level_per_pair_or=";".join(f"{v:.3g}" for v in np.exp(log_or)),
            )
        # baseline vs batch: 2022 pair vs pooled 2021 pairs (unmodified)
        batches = sorted(set(batch_of_pair.values()))
        if len(batches) == 2:
            in_b = np.array(
                [batch_of_pair[p] == batches[1] for p in layout.pair_ids], dtype=bool
            )
            tab = np.array(
                [
                    [int(bt[in_b].sum()), int(bc[in_b].sum())],
                    [int(bt[~in_b].sum()), int(bc[~in_b].sum())],
                ]
            )
            f = pm.fisher_odds_ratio(tab, alpha=alpha)
            batch_tests[defn] = {
                "comparison": f"unmodified treated-only share, {batches[1]} pair vs "
                f"{batches[0]} pairs (Fisher exact OR)",
                "table": tab.tolist(),
                **f,
            }
    table = pd.DataFrame(rows)
    for defn in DETECTION_DEFINITIONS:
        m = table["definition"] == defn
        table.loc[m, "q_binomial_vs_0.5_BH"] = pm.bh_adjust(
            table.loc[m, "p_binomial_vs_0.5"].to_numpy(dtype=float)
        )
        nb = m & (table["mod_class"] != pm.UNMODIFIED)
        for col in ("fisher", "mh", "peptide_level", "pair_level"):
            table.loc[nb, f"{col}_q_BH"] = pm.bh_adjust(
                table.loc[nb, f"{col}_p"].to_numpy(dtype=float)
            )
    return table, pd.DataFrame(pair_rows), batch_tests


def discordant_peptide_list(
    features: pd.DataFrame,
    masks: dict[str, np.ndarray],
    layout: pm.PairLayout,
    inputs: Inputs,
    min_same: int,
) -> pd.DataFrame:
    """Peptides discordant in >= ``min_same`` pairs one way and 0 the other way."""
    samples = _sample_ids(inputs.raw)
    keep = ~features["is_contaminant"].to_numpy(dtype=bool)
    frames = []
    for defn in DETECTION_DEFINITIONS:
        n_t = features[f"{defn}_n_treated_only"].to_numpy()
        n_c = features[f"{defn}_n_control_only"].to_numpy()
        up = keep & (n_t >= min_same) & (n_c == 0)
        down = keep & (n_c >= min_same) & (n_t == 0)
        sel = np.flatnonzero(up | down)
        sub = features.iloc[sel][
            [
                "feature",
                "base_sequence",
                "accessions",
                "entries",
                "mod_class",
                "composition",
                "dump_sites",
                f"{defn}_n_treated_only",
                f"{defn}_n_control_only",
                f"{defn}_n_both",
                f"{defn}_n_neither",
            ]
        ].rename(
            columns={
                f"{defn}_n_treated_only": "n_treated_only",
                f"{defn}_n_control_only": "n_control_only",
                f"{defn}_n_both": "n_both",
                f"{defn}_n_neither": "n_neither",
            }
        )
        sub.insert(0, "definition", defn)
        sub["direction"] = np.where(up[sel], "treated_only", "control_only")
        sub["all_pairs"] = np.maximum(
            sub["n_treated_only"], sub["n_control_only"]
        ) == len(layout.pair_ids)
        inten = inputs.raw.abundances[:, sel]
        for j, s in enumerate(samples):
            sub[f"{s}_intensity"] = inten[j]
            sub[f"{s}_detection_type"] = inputs.detection_type[j, sel]
        present = np.full(len(sel), np.nan)
        anyq = np.isfinite(inten).any(axis=0)
        if anyq.any():
            present[anyq] = np.nanmedian(inten[:, anyq], axis=0)
        sub["median_intensity_where_quantified"] = present
        frames.append(sub)
    out = pd.concat(frames, ignore_index=True)
    return out.sort_values(
        ["definition", "direction", "all_pairs", "median_intensity_where_quantified"],
        ascending=[True, False, False, False],
    )


# --------------------------------------------------------------------------- #
# 4. Intensity by class
# --------------------------------------------------------------------------- #
def _class_vs_unmodified(
    values: pd.Series, classes: pd.Series, alpha: float, effect: str
) -> pd.DataFrame:
    """Per class summary + Welch and Mann-Whitney vs unmodified; BH over classes."""
    base = values[classes == pm.UNMODIFIED].to_numpy(dtype=float)
    rows = []
    for cls in pm.MOD_CLASSES:
        x = values[classes == cls].to_numpy(dtype=float)
        row: dict[str, object] = {
            "mod_class": cls,
            "effect": effect,
            "n": len(x),
            "mean": float(x.mean()) if len(x) else np.nan,
            "median": float(np.median(x)) if len(x) else np.nan,
            "sd": float(x.std(ddof=1)) if len(x) > 1 else np.nan,
            "q25": float(np.quantile(x, 0.25)) if len(x) else np.nan,
            "q75": float(np.quantile(x, 0.75)) if len(x) else np.nan,
            "frac_positive": float((x > 0).mean()) if len(x) else np.nan,
        }
        if cls != pm.UNMODIFIED and len(x) >= 2:
            w = pm.welch_difference(x, base, alpha=alpha)
            mw = pm.mann_whitney_shift(x, base, alpha=alpha)
            row.update(
                welch_mean_diff_vs_unmodified=w["diff"],
                welch_ci_low=w["ci_low"],
                welch_ci_high=w["ci_high"],
                welch_p=w["p"],
                hl_shift_vs_unmodified=mw["hl_shift"],
                hl_ci_low=mw["ci_low"],
                hl_ci_high=mw["ci_high"],
                rank_biserial=mw["rank_biserial"],
                mannwhitney_p=mw["p"],
            )
        rows.append(row)
    out = pd.DataFrame(rows)
    for col in ("welch_p", "mannwhitney_p"):
        if col not in out:
            out[col] = np.nan
        out[col.replace("_p", "_q_BH")] = pm.bh_adjust(out[col].to_numpy(dtype=float))
    return out


def _add_pair_level(
    tests: pd.DataFrame,
    diffs: np.ndarray,
    classes: np.ndarray,
    pair_ids: Sequence[str],
    alpha: float,
) -> pd.DataFrame:
    """Append the pair-level class shift (t across pairs, df 3) + BH over classes."""
    shift = pm.pair_level_class_shift(diffs, classes, pm.UNMODIFIED, alpha=alpha)
    out = tests.copy()
    for col in (
        "pair_shift_mean",
        "pair_shift_ci_low",
        "pair_shift_ci_high",
        "pair_shift_p",
        "pair_shift_n_pairs_positive",
    ):
        out[col] = np.nan
    out["pair_shift_per_pair"] = ""
    for i, cls in enumerate(out["mod_class"]):
        if cls not in shift:
            continue
        r = shift[str(cls)]
        out.loc[i, "pair_shift_mean"] = r["mean"]
        out.loc[i, "pair_shift_ci_low"] = r["ci_low"]
        out.loc[i, "pair_shift_ci_high"] = r["ci_high"]
        out.loc[i, "pair_shift_p"] = r["p"]
        out.loc[i, "pair_shift_n_pairs_positive"] = r["n_pairs_positive"]
        per_pair = r["per_pair_shift"]
        assert isinstance(per_pair, list)
        out.loc[i, "pair_shift_per_pair"] = ";".join(
            f"{p}:{v:.4f}" for p, v in zip(pair_ids, per_pair, strict=True)
        )
    out["pair_shift_q_BH"] = pm.bh_adjust(out["pair_shift_p"].to_numpy(dtype=float))
    return out


def intensity_by_class(
    features: pd.DataFrame,
    inputs: Inputs,
    layout: pm.PairLayout,
    alpha: float,
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Paired log2FC per class, class tests, top peptides by q and by |log2FC|."""
    de = inputs.de
    paired = de["paired"].merge(
        features[
            [
                "feature",
                "accessions",
                "entries",
                "mod_class",
                "composition",
                "n_dump_isoforms",
                "dump_sites",
            ]
        ],
        on="feature",
        how="left",
        validate="one_to_one",
    )
    if paired["mod_class"].isna().any():
        raise ValueError("Paired DE features missing from the feature table.")
    for design in ("batch", "unadjusted"):
        paired = paired.merge(
            de[design][["feature", "q"]].rename(columns={"q": f"q_{design}"}),
            on="feature",
            how="left",
            validate="one_to_one",
        )
    tests = _class_vs_unmodified(
        paired["log2fc"], paired["mod_class"], alpha, "paired log2FC (ralox - control)"
    )
    pep_names = [str(n) for n in inputs.norm_log.feature_names]
    cls_by_feat = features.set_index("feature")["mod_class"]
    tests = _add_pair_level(
        tests,
        _pair_diffs(inputs.norm_log, layout),
        cls_by_feat.loc[pep_names].to_numpy(dtype=str),
        layout.pair_ids,
        alpha,
    )
    cols = [
        "feature",
        "base_sequence",
        "accessions",
        "entries",
        "mod_class",
        "composition",
        "log2fc",
        "ci_low",
        "ci_high",
        "se",
        "p",
        "q",
        "q_batch",
        "q_unadjusted",
        "mean_log2_abundance",
        "n_dump_isoforms",
        "dump_sites",
    ]
    by_q = paired.sort_values(["q", "p"]).head(top_n)[cols]
    by_fc = (
        paired.assign(abs_log2fc=paired["log2fc"].abs())
        .sort_values("abs_log2fc", ascending=False)
        .head(top_n)[cols]
    )
    return paired[cols], tests, by_q, by_fc


# --------------------------------------------------------------------------- #
# 5. Protein-adjusted and modified-minus-reference
# --------------------------------------------------------------------------- #
def _pair_diffs(ds: Dataset, layout: pm.PairLayout) -> np.ndarray:
    """``(n_pairs, n_features)`` treated - control, per pair (log2 input)."""
    x = ds.abundances
    return np.asarray(x[layout.treated_idx] - x[layout.control_idx], dtype=float)


def protein_adjusted(
    features: pd.DataFrame,
    inputs: Inputs,
    layout: pm.PairLayout,
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Peptide - protein paired log2FC per unique peptide; moderated one-sample t."""
    pep_names = [str(n) for n in inputs.norm_log.feature_names]
    prot_names = [str(n) for n in inputs.protein_log.feature_names]
    pep_d = _pair_diffs(inputs.norm_log, layout)
    prot_d = _pair_diffs(inputs.protein_log, layout)

    # reproduce the existing paired DE contrast (implementation check)
    check = pm.moderated_one_sample(pep_d, alpha=alpha)
    paired = inputs.de["paired"].set_index("feature").loc[pep_names]
    repro = {
        "max_abs_diff_log2fc": float(
            np.max(np.abs(check.mean - paired["log2fc"].to_numpy()))
        ),
        "max_abs_diff_p": float(np.max(np.abs(check.p - paired["p"].to_numpy()))),
        "max_abs_diff_q": float(np.max(np.abs(check.q - paired["q"].to_numpy()))),
        "max_abs_diff_ci_low": float(
            np.max(np.abs(check.ci_low - paired["ci_low"].to_numpy()))
        ),
    }
    if repro["max_abs_diff_p"] > 1e-8 or repro["max_abs_diff_log2fc"] > 1e-10:
        raise ValueError(
            f"Moderated one-sample t does not reproduce peptide_paired.tsv: {repro}."
        )

    feat = features.set_index("feature").loc[pep_names]
    prot_index = {n: i for i, n in enumerate(prot_names)}
    groups = feat["protein_groups"].astype(str).to_numpy()
    unique = feat["n_protein_members"].to_numpy(dtype=int) == 1
    col = np.array(
        [
            prot_index.get(g, -1) if u else -1
            for g, u in zip(groups, unique, strict=True)
        ]
    )
    usable = col >= 0
    adj_raw = pep_d[:, usable] - prot_d[:, col[usable]]
    # peptide and protein levels were median-normalized separately, so each pair
    # carries a level offset; centre each pair on its median peptide-minus-protein.
    offsets = np.median(adj_raw, axis=1)
    adj = adj_raw - offsets[:, None]
    res = pm.moderated_one_sample(adj, alpha=alpha)

    # peptide support per protein (raw, any quantified value, unique peptides)
    all_feat = features[
        (~features["is_contaminant"]) & (features["n_protein_members"] == 1)
    ]
    quant_any = np.isfinite(inputs.raw.abundances).any(axis=0)
    all_feat = all_feat[quant_any[all_feat.index.to_numpy()]]
    n_base = all_feat.groupby("protein_groups")["base_sequence"].nunique()
    n_complete = pd.Series(groups[usable]).value_counts()

    out = pd.DataFrame(
        {
            "feature": np.array(pep_names)[usable],
            "protein_group": groups[usable],
            "accessions": feat["accessions"].to_numpy()[usable],
            "entries": feat["entries"].to_numpy()[usable],
            "mod_class": feat["mod_class"].to_numpy()[usable],
            "composition": feat["composition"].to_numpy()[usable],
            "peptide_log2fc": pep_d[:, usable].mean(axis=0),
            "protein_log2fc": prot_d[:, col[usable]].mean(axis=0),
            "adjusted_log2fc": res.mean,
            "ci_low": res.ci_low,
            "ci_high": res.ci_high,
            "se": res.se,
            "t": res.t,
            "p": res.p,
            "q": res.q,
        }
    )
    for k, pair in enumerate(layout.pair_ids):
        out[f"adjusted_{pair}"] = adj[k]
    out["protein_n_base_sequences_quantified"] = (
        out["protein_group"].map(n_base).fillna(0).astype(int)
    )
    out["protein_n_peptides_complete"] = (
        out["protein_group"].map(n_complete).astype(int)
    )
    out["single_peptide_protein"] = out["protein_n_base_sequences_quantified"] <= 1
    out = out.sort_values(["q", "p"]).reset_index(drop=True)
    tests = _class_vs_unmodified(
        out["adjusted_log2fc"],
        out["mod_class"],
        alpha,
        "peptide - protein paired log2FC (pair-median centred)",
    )
    tests = _add_pair_level(
        tests,
        adj,
        feat["mod_class"].to_numpy(dtype=str)[usable],
        layout.pair_ids,
        alpha,
    )
    info: dict[str, object] = {
        "reproduces_peptide_paired_tsv": repro,
        "n_complete_peptides": len(pep_names),
        "n_shared_peptides_excluded": int((~unique).sum()),
        "n_unique_peptides_protein_not_complete": int((unique & (col < 0)).sum()),
        "n_tested": int(usable.sum()),
        "per_pair_centering_offset": dict(
            zip(layout.pair_ids, offsets.tolist(), strict=True)
        ),
        "n_single_peptide_protein": int(out["single_peptide_protein"].sum()),
        "prior_s0_sq": res.s0_sq,
        "prior_d0": res.d0,
        "n_q_lt_0.05": int((out["q"] < 0.05).sum()),
        "n_q_lt_0.10": int((out["q"] < 0.10).sum()),
        "n_q_lt_0.05_multi_peptide": int(
            ((out["q"] < 0.05) & ~out["single_peptide_protein"]).sum()
        ),
        "min_q": float(out["q"].min()),
    }
    return out, tests, info


def modified_vs_reference(
    features: pd.DataFrame,
    inputs: Inputs,
    layout: pm.PairLayout,
    ref_idx: np.ndarray,
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Modified - reference-form paired log2FC (occupancy-like), complete forms only."""
    pep_names = [str(n) for n in inputs.norm_log.feature_names]
    pos = {n: i for i, n in enumerate(pep_names)}
    names = features["feature"].astype(str).to_numpy()
    pairs_mr = [
        (i, int(r))
        for i, r in enumerate(ref_idx)
        if r >= 0
        and names[i] in pos
        and names[r] in pos
        and not features["is_contaminant"].iloc[i]
    ]
    n_modified_complete = int(
        sum(
            1
            for i, n in enumerate(names)
            if n in pos
            and features["mod_class"].iloc[i] not in (pm.UNMODIFIED, pm.UNEXPLAINED)
        )
    )
    pep_d = _pair_diffs(inputs.norm_log, layout)
    mi = np.array([pos[names[i]] for i, _ in pairs_mr], dtype=np.int64)
    ri = np.array([pos[names[r]] for _, r in pairs_mr], dtype=np.int64)
    diff = pep_d[:, mi] - pep_d[:, ri]
    res = pm.moderated_one_sample(diff, alpha=alpha)
    fi = np.array([i for i, _ in pairs_mr], dtype=np.int64)
    out = pd.DataFrame(
        {
            "feature": names[fi],
            "reference_feature": [names[r] for _, r in pairs_mr],
            "accessions": features["accessions"].to_numpy()[fi],
            "entries": features["entries"].to_numpy()[fi],
            "mod_class": features["mod_class"].to_numpy()[fi],
            "composition": features["composition"].to_numpy()[fi],
            "dump_sites": features["dump_sites"].to_numpy()[fi],
            "modified_log2fc": pep_d[:, mi].mean(axis=0),
            "reference_log2fc": pep_d[:, ri].mean(axis=0),
            "modified_minus_reference": res.mean,
            "ci_low": res.ci_low,
            "ci_high": res.ci_high,
            "se": res.se,
            "t": res.t,
            "p": res.p,
            "q": res.q,
        }
    )
    for k, pair in enumerate(layout.pair_ids):
        out[f"diff_{pair}"] = diff[k]
    out = out.sort_values(["q", "p"]).reset_index(drop=True)
    rows = []
    for cls in (pm.CAM_ONLY, pm.OXIDIZED, pm.ADDUCT):
        x = out.loc[out["mod_class"] == cls, "modified_minus_reference"].to_numpy(
            dtype=float
        )
        sel = (out["mod_class"] == cls).to_numpy()
        nan = float("nan")
        if len(x) >= 2:
            t = pm.paired_t(x, alpha=alpha)
            w = float(stats.wilcoxon(x).pvalue)
            per_pair = np.array(
                [out.loc[sel, f"diff_{p}"].mean() for p in layout.pair_ids]
            )
            pl = pm.paired_t(per_pair, alpha=alpha)
        else:
            t = {"mean": nan, "ci_low": nan, "ci_high": nan, "p": nan, "df": nan}
            w = nan
            per_pair = np.full(len(layout.pair_ids), nan)
            pl = {"mean": nan, "ci_low": nan, "ci_high": nan, "p": nan, "df": nan}
        rows.append(
            {
                "mod_class": cls,
                "effect": "mean over forms of (modified - reference) paired log2FC",
                "n_forms": len(x),
                "mean": t["mean"],
                "ci_low": t["ci_low"],
                "ci_high": t["ci_high"],
                "median": float(np.median(x)) if len(x) else np.nan,
                "t_test_p": t["p"],
                "wilcoxon_p": w,
                "pair_level_mean": pl["mean"],
                "pair_level_ci_low": pl["ci_low"],
                "pair_level_ci_high": pl["ci_high"],
                "pair_level_p": pl["p"],
                "pair_level_per_pair": ";".join(
                    f"{p}:{v:.4f}"
                    for p, v in zip(layout.pair_ids, per_pair, strict=True)
                ),
            }
        )
    tests = pd.DataFrame(rows)
    tests["t_test_q_BH"] = pm.bh_adjust(tests["t_test_p"].to_numpy(dtype=float))
    tests["wilcoxon_q_BH"] = pm.bh_adjust(tests["wilcoxon_p"].to_numpy(dtype=float))
    tests["pair_level_q_BH"] = pm.bh_adjust(tests["pair_level_p"].to_numpy(dtype=float))
    info: dict[str, object] = {
        "n_modified_complete": n_modified_complete,
        "n_with_complete_reference": len(pairs_mr),
        "prior_s0_sq": res.s0_sq,
        "prior_d0": res.d0,
        "n_q_lt_0.05": int((out["q"] < 0.05).sum()),
        "n_q_lt_0.10": int((out["q"] < 0.10).sum()),
        "min_q": float(out["q"].min()),
    }
    return out, tests, info


# --------------------------------------------------------------------------- #
# 6. Oxidation index
# --------------------------------------------------------------------------- #
def _ox_sums(
    features: pd.DataFrame, x: np.ndarray, eligible_feature: np.ndarray
) -> pd.DataFrame:
    """Per base sequence (with >=1 oxidized and >=1 non-oxidized eligible form):
    per-sample sums of oxidized and of all forms. ``x`` is (n_samples, n_features)
    linear with NaN -> treated as 0 by the caller."""
    ncls = features["mod_class"].to_numpy(dtype=str)
    is_ox = (ncls == pm.OXIDIZED) & eligible_feature
    is_non = np.isin(ncls, [pm.UNMODIFIED, pm.CAM_ONLY]) & eligible_feature
    base = features["base_sequence"].to_numpy(dtype=str)
    ok_bases = set(base[is_ox]) & set(base[is_non])
    use = (is_ox | is_non) & np.isin(base, list(ok_bases))
    frame = pd.DataFrame(np.nan_to_num(x[:, use]).T, columns=range(x.shape[0]))
    frame["base_sequence"] = base[use]
    frame["protein_groups"] = features["protein_groups"].to_numpy(dtype=str)[use]
    frame["n_protein_members"] = features["n_protein_members"].to_numpy()[use]
    frame["is_ox"] = is_ox[use]
    return frame


def oxidation_index(
    features: pd.DataFrame,
    inputs: Inputs,
    layout: pm.PairLayout,
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Per-protein and global oxidation index, raloxifene vs control (paired)."""
    samples = _sample_ids(inputs.raw)
    md = inputs.raw.metadata
    raw = inputs.raw.abundances
    non_contam = ~features["is_contaminant"].to_numpy(dtype=bool)
    complete = np.isfinite(raw).all(axis=0) & non_contam
    n_s = len(samples)

    def per_sample_index(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        ox = frame.loc[frame["is_ox"], list(range(n_s))].sum(axis=0).to_numpy()
        tot = frame[list(range(n_s))].sum(axis=0).to_numpy()
        return np.asarray(ox, dtype=float), np.asarray(tot, dtype=float)

    results: dict[str, object] = {}
    sample_rows = []
    for variant, eligible in (
        ("complete_case", complete),
        ("missing_as_zero", non_contam),
    ):
        frame = _ox_sums(features, raw, eligible)
        ox, tot = per_sample_index(frame)
        idx = ox / tot
        log_ratio = np.log2(idx[layout.treated_idx] / idx[layout.control_idx])
        t = pm.paired_t(log_ratio, alpha=alpha)
        n_pos = int((log_ratio > 0).sum())
        results[variant] = {
            "n_base_sequences": int(frame["base_sequence"].nunique()),
            "n_forms": len(frame),
            "per_pair_log2_ratio_treated_over_control": dict(
                zip(layout.pair_ids, log_ratio.tolist(), strict=True)
            ),
            "mean_log2_ratio": t["mean"],
            "ci_low": t["ci_low"],
            "ci_high": t["ci_high"],
            "paired_t_p": t["p"],
            "df": t["df"],
            "sign_test_p": pm.sign_test(n_pos, len(log_ratio) - n_pos),
            "test": "paired t on per-pair log2(global index treated / control), "
            "df 3; exact sign test (min p 0.125)",
        }
        for j, s in enumerate(samples):
            sample_rows.append(
                {
                    "variant": variant,
                    "sample_id": s,
                    "condition": md[CONDITION_COL].iloc[j],
                    "candidate_pair": md[PAIR_COL].iloc[j],
                    "batch": md["batch"].iloc[j],
                    "run_position_within_batch": int(
                        md["run_position_within_batch"].iloc[j]
                    ),
                    "oxidized_sum": ox[j],
                    "total_sum": tot[j],
                    "oxidation_index": idx[j],
                }
            )
    # per protein (complete case, unique peptides)
    frame = _ox_sums(features, raw, complete)
    frame = frame[frame["n_protein_members"] == 1]
    cols = list(range(n_s))
    ox_p = frame[frame["is_ox"]].groupby("protein_groups")[cols].sum()
    tot_p = frame.groupby("protein_groups")[cols].sum().loc[ox_p.index]
    nbase = frame.groupby("protein_groups")["base_sequence"].nunique().loc[ox_p.index]
    idx_p = (ox_p / tot_p).to_numpy()  # (n_prot, n_samples)
    lr = np.log2(idx_p[:, layout.treated_idx] / idx_p[:, layout.control_idx]).T
    prot = pd.DataFrame({"protein_group": ox_p.index.to_numpy(dtype=str)})
    members = [pm.parse_protein_members(g)[0] for g in prot["protein_group"]]
    prot["accession"] = [m[1] for m in members]
    prot["entry"] = [m[2] for m in members]
    prot["n_base_sequences"] = nbase.to_numpy()
    for j, s in enumerate(samples):
        prot[f"{s}_oxidation_index"] = idx_p[:, j]
    for k, pair in enumerate(layout.pair_ids):
        prot[f"log2_ratio_{pair}"] = lr[k]
    if lr.shape[1] >= 3:
        res = pm.moderated_one_sample(lr, alpha=alpha)
        prot["mean_log2_ratio"] = res.mean
        prot["ci_low"] = res.ci_low
        prot["ci_high"] = res.ci_high
        prot["t"] = res.t
        prot["p"] = res.p
        prot["q"] = res.q
        # pair-level summary: median over proteins per pair
        med = np.median(lr, axis=1)
        tm = pm.paired_t(med, alpha=alpha)
        results["per_protein"] = {
            "n_proteins": int(lr.shape[1]),
            "prior_s0_sq": res.s0_sq,
            "prior_d0": res.d0,
            "n_q_lt_0.05": int((res.q < 0.05).sum()),
            "n_q_lt_0.10": int((res.q < 0.10).sum()),
            "min_q": float(np.nanmin(res.q)),
            "frac_proteins_mean_log2_ratio_positive": float((res.mean > 0).mean()),
            "per_pair_median_log2_ratio": dict(
                zip(layout.pair_ids, med.tolist(), strict=True)
            ),
            "median_over_proteins_paired_t": tm,
        }
        prot = prot.sort_values(["q", "p"])
    return prot, pd.DataFrame(sample_rows), results


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def _json_default(obj: object) -> object:
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Not JSON serializable: {type(obj)}")


def _clean(obj: object) -> object:
    """Replace float NaN/inf with None recursively (strict JSON)."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, list | tuple):
        return [_clean(v) for v in obj]
    if isinstance(obj, float) and not np.isfinite(obj):
        return None if np.isnan(obj) else ("inf" if obj > 0 else "-inf")
    return obj


def _records(frame: pd.DataFrame) -> list[dict[str, object]]:
    return [
        {str(k): _clean(v) for k, v in row.items()}
        for row in json.loads(frame.to_json(orient="records"))
    ]


def run(
    root: Path, out_dir: Path, tolerance: float, alpha: float, top_n: int, min_same: int
) -> dict[str, object]:
    """Run all six analyses; write tables + summary.json; return the summary."""
    inputs = load_inputs(root)
    out_dir.mkdir(parents=True, exist_ok=True)
    md = inputs.raw.metadata
    layout = pm.pair_layout(
        md,
        pair_col=PAIR_COL,
        condition_col=CONDITION_COL,
        control=CONTROL,
        treated=TREATED,
    )
    batch_of_pair = {
        str(p): str(md["batch"].iloc[int(i)])
        for p, i in zip(layout.pair_ids, layout.control_idx, strict=True)
    }
    LOG.info("pairs %s; data_version %s", layout.pair_ids, inputs.data_version)

    # 1. decomposition + join
    features = build_feature_table(inputs.raw, tolerance)
    join = pm.join_dump_to_features(
        inputs.dump.rows["key"].tolist(), features["key"].tolist()
    )
    dump_psms, dump_mbr, dump_text = pm.aggregate_dump_to_features(
        join, inputs.dump, len(features)
    )
    features = pd.concat([features, dump_text], axis=1)
    features["n_dump_isoforms"] = join.n_isomers
    xcheck = dump_crosschecks(inputs.detection_type, dump_psms, dump_mbr)
    unmatched = pd.DataFrame(
        {
            "side": ["dump"] * len(join.unmatched_dump_keys)
            + ["quant"] * len(join.unmatched_feature_keys),
            "key": list(join.unmatched_dump_keys) + list(join.unmatched_feature_keys),
        }
    )
    unmatched.to_csv(out_dir / "dump_join_unmatched.tsv", sep="\t", index=False)
    LOG.info(
        "dump join: %d dump rows, %d unmatched dump, %d unmatched quant",
        len(inputs.dump.rows),
        len(join.unmatched_dump_keys),
        len(join.unmatched_feature_keys),
    )

    masks = detection_masks(inputs.detection_type, dump_psms)
    for defn in DETECTION_DEFINITIONS:
        counts = pm.discordance_counts(masks[defn], layout)
        for col in counts.columns:
            features[f"{defn}_{col}"] = counts[col].to_numpy()
        features[f"{defn}_n_samples"] = masks[defn].sum(axis=0)
    ref_idx = _reference_index(features)
    features["reference_feature"] = [
        str(features["feature"].iloc[r]) if r >= 0 else "" for r in ref_idx
    ]
    features["complete_case"] = features["feature"].isin(
        [str(n) for n in inputs.norm_log.feature_names]
    )
    mass_tab = mass_decomposition_table(features)
    mass_tab.to_csv(out_dir / "mass_decomposition.tsv", sep="\t", index=False)
    features.drop(columns=["key"]).to_csv(
        out_dir / "peptide_classes.tsv", sep="\t", index=False
    )

    # 2. adducts
    a_wide, a_long, a_pairs, a_prot, a_extra = adduct_tables(
        features, inputs, masks, dump_psms, dump_mbr, layout, ref_idx, alpha
    )
    a_wide.to_csv(out_dir / "adduct_peptides.tsv", sep="\t", index=False)
    a_long.to_csv(out_dir / "adduct_peptides_long.tsv", sep="\t", index=False)
    a_pairs.to_csv(out_dir / "adduct_pair_summary.tsv", sep="\t", index=False)
    a_prot.to_csv(out_dir / "adduct_proteins.tsv", sep="\t", index=False)
    ctrl_det = a_extra["ctrl"]
    assert isinstance(ctrl_det, pd.DataFrame)
    ctrl_det.to_csv(out_dir / "adduct_control_detections.tsv", sep="\t", index=False)

    # 3. presence/absence
    d_tests, d_pairs, d_batch = discordance_class_tests(
        features, masks, layout, batch_of_pair, alpha
    )
    d_tests.to_csv(out_dir / "discordance_class_tests.tsv", sep="\t", index=False)
    d_pairs.to_csv(out_dir / "discordance_by_pair.tsv", sep="\t", index=False)
    d_list = discordant_peptide_list(features, masks, layout, inputs, min_same)
    d_list.to_csv(out_dir / "discordant_peptides.tsv", sep="\t", index=False)

    # 4. intensity by class
    fc, fc_tests, top_q, top_fc = intensity_by_class(
        features, inputs, layout, alpha, top_n
    )
    fc.to_csv(out_dir / "class_log2fc.tsv", sep="\t", index=False)
    fc_tests.to_csv(out_dir / "class_log2fc_tests.tsv", sep="\t", index=False)
    top_q.to_csv(out_dir / "top_peptides_by_q.tsv", sep="\t", index=False)
    top_fc.to_csv(out_dir / "top_peptides_by_abs_log2fc.tsv", sep="\t", index=False)

    # 5. protein-adjusted + modified-vs-reference
    padj, padj_tests, padj_info = protein_adjusted(features, inputs, layout, alpha)
    padj.to_csv(out_dir / "protein_adjusted.tsv", sep="\t", index=False)
    padj_tests.to_csv(
        out_dir / "protein_adjusted_class_tests.tsv", sep="\t", index=False
    )
    mvr, mvr_tests, mvr_info = modified_vs_reference(
        features, inputs, layout, ref_idx, alpha
    )
    mvr.to_csv(out_dir / "modified_vs_reference.tsv", sep="\t", index=False)
    mvr_tests.to_csv(
        out_dir / "modified_vs_reference_class_tests.tsv", sep="\t", index=False
    )

    # 6. oxidation index
    ox_prot, ox_samples, ox_res = oxidation_index(features, inputs, layout, alpha)
    ox_prot.to_csv(out_dir / "oxidation_index_protein.tsv", sep="\t", index=False)
    ox_samples.to_csv(out_dir / "oxidation_index_sample.tsv", sep="\t", index=False)

    # CYP3A4 context
    cyp = features[features["accessions"].str.contains(CYP3A4_ACCESSION)]
    cyp3a4 = {
        "n_features": len(cyp),
        "by_class": {
            str(k): int(v) for k, v in cyp["mod_class"].value_counts().items()
        },
        "n_adduct_features": int((cyp["mod_class"] == pm.ADDUCT).sum()),
        "n_cys_containing_features": int((cyp["n_cys"] > 0).sum()),
    }

    classes_nc = features.loc[~features["is_contaminant"], "mod_class"].value_counts()
    summary: dict[str, object] = {
        "task": "peptide-mod-analysis",
        "data_version": inputs.data_version,
        "script": "scripts/scratch/peptide_mod_analysis.py",
        "script_sha256": sha256_of_file(Path(__file__)),
        "module_sha256": sha256_of_file(Path(pm.__file__)),
        "input_sha256": inputs.input_hashes,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "params": {
            "mass_tolerance": tolerance,
            "alpha": alpha,
            "top_n": top_n,
            "min_same_direction_pairs": min_same,
            "sample_set": "experimental (all 8 runs; no QC/pool controls exist)",
            "pairs": list(layout.pair_ids),
            "contrast": f"{TREATED} - {CONTROL} (positive = higher in raloxifene)",
            "detection_definitions": DEFINITION_TEXT,
            "seed": None,
        },
        "mass_arithmetic": {
            "raloxifene_formula": "C28H27NO4S",
            "raloxifene_monoisotopic": pm.RALOXIFENE_MONOISOTOPIC,
            "adduct_delta_M_minus_2H": pm.RALOXIFENE_ADDUCT_DELTA,
            "quant_file_adduct_mass": 471.1504,
            "error_quant_minus_derived": 471.1504 - pm.RALOXIFENE_ADDUCT_DELTA,
            "deltas": pm.MOD_DELTAS,
        },
        "decomposition": {
            "n_features": len(features),
            "n_distinct_masses": len(mass_tab),
            "n_unexplained_features": int(
                (features["mod_class"] == pm.UNEXPLAINED).sum()
            ),
            "n_ambiguous_masses": int((mass_tab["n_solutions"] > 1).sum()),
            "max_abs_mass_error": float(np.nanmax(np.abs(mass_tab["mass_error"]))),
            "n_cam_exceeds_cys": int(features["cam_exceeds_cys"].sum()),
            "class_counts_all": {
                str(k): int(v) for k, v in features["mod_class"].value_counts().items()
            },
            "class_counts_non_contaminant": {
                str(k): int(v) for k, v in classes_nc.items()
            },
        },
        "dump_join": {
            "n_dump_rows": len(inputs.dump.rows),
            "n_quant_features": len(features),
            "n_unmatched_dump_rows": len(join.unmatched_dump_keys),
            "n_unmatched_quant_features": len(join.unmatched_feature_keys),
            "isoforms_per_feature": {
                str(k): int(v)
                for k, v in pd.Series(join.n_isomers)
                .value_counts()
                .sort_index()
                .items()
            },
            "crosschecks": xcheck,
        },
        "adducts": {
            "n_adduct_features": len(a_wide),
            "n_adduct_features_non_contaminant": int((~a_wide["is_contaminant"]).sum()),
            "residue_counts": {
                str(k): int(v)
                for k, v in a_wide["dump_ralox_residues"].value_counts().items()
            },
            "n_proteins": int(a_prot["entry"].nunique()),
            "cyp_proteins": a_prot.loc[a_prot["is_cyp"], "entry"].tolist(),
            "cyp3a4_P08684_adducted": bool(
                (a_prot["accession"] == CYP3A4_ACCESSION).any()
            ),
            "cyp3a4_context": cyp3a4,
            "tests": a_extra["tests"],
        },
        "presence_absence": {
            "class_tests": _records(d_tests),
            "baseline_by_batch": d_batch,
            "n_discordant_list_rows": len(d_list),
            "correction": "BH over classes within each detection definition "
            "(vs-0.5 family: 5 classes; vs-unmodified families: 4 classes)",
        },
        "intensity_by_class": {
            "n_complete_peptides": len(fc),
            "class_tests": _records(fc_tests),
            "n_q_lt_0.10": int((fc["q"] < 0.10).sum()),
            "correction": "BH over the 4 non-baseline classes, per test",
        },
        "protein_adjusted": {
            **padj_info,
            "class_tests": _records(padj_tests),
        },
        "modified_vs_reference": {
            **mvr_info,
            "class_tests": _records(mvr_tests),
        },
        "oxidation": ox_res,
        "caveats": [
            "4 matched pairs, residual df 3: all results exploratory.",
            "Run order aliased with condition (control first in every pair; "
            "finding 0001): treated-only detections and any shift can reflect "
            "run position/column state as well as treatment.",
            "MBR can transfer identifications across runs; msms/psm definitions are "
            "the sensitivity analyses.",
            "Discordant cells are not independent across pairs of one peptide or "
            "peptides of one run; the pair-stratified MH and peptide-level collapse "
            "are the dependence sensitivities.",
            "Peptide positions within proteins are not derivable without a FASTA; "
            "sites are residue + position within the peptide.",
        ],
    }
    (out_dir / "summary.json").write_text(
        json.dumps(_clean(summary), indent=2, default=_json_default), encoding="utf-8"
    )
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""
    root_default = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--root", type=Path, default=root_default)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument(
        "--mass-tolerance", type=float, default=pm.DEFAULT_MASS_TOLERANCE
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=30)
    parser.add_argument("--min-same-direction", type=int, default=3)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    out_dir = args.out_dir or args.root / "results" / "peptide-mods"
    LOG.info("params: %s", vars(args))
    run(
        args.root,
        out_dir,
        args.mass_tolerance,
        args.alpha,
        args.top_n,
        args.min_same_direction,
    )
    LOG.info("wrote %s", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
