"""Stage-3 QC "prep-once": materialize the processing-state matrices to disk.

For each level in {protein, peptide}, loads the real data via the verified loaders
and writes a sequence of processing states (``lib/common/dataset_io.py``,
``dataset-io`` v0.1, seeded as ``loaders/dataset_io.py``) under
``results/qc_states/<level>/<state>/`` so downstream QC figure jobs each
``load_dataset`` the one state they need instead of re-loading / re-normalizing /
re-running ComBat (the workflow's "prep-once" pattern):

  * ``raw_linear``            -- full matrix, all features incl. contaminants
                                  (flagged in feature_metadata), NaN = missing.
  * ``raw_linear_complete``   -- non-contaminant features quantified in all 8 runs
                                  (the scientist-confirmed analysis set,
                                  ``state/DATA_DESCRIPTION.md`` "Preprocessing
                                  decisions" -- ``handle_missing`` with
                                  ``max_missing_fraction=0``).
  * ``raw_log``               -- log2(raw_linear_complete).
  * ``normalized_linear``     -- median-normalized raw_linear_complete (``lib/
                                  common/normalize.py`` ``normalize``, method
                                  ``"median"``; linear in, linear out --
                                  scale-rescaled to the mean of per-sample medians).
  * ``normalized_log``        -- log2(normalized_linear).
  * ``batch_corrected_log``   -- ComBat on normalized_log, BATCH LABEL ONLY
                                  (``lib/common/batch_correct.py``, never the
                                  condition). Recorded, not silently worked around,
                                  if it fails or warns (B2022 has only 1 control + 1
                                  raloxifene sample -- a fragile 2-sample batch).
  * ``batch_corrected_linear``-- to_linear(batch_corrected_log), when it succeeded.

Also, an ``nsaf`` level (scientist request, protein-group level, no normalization/
ComBat -- NSAF is already normalized): ``raw_linear`` (all 4,344 Limelight dump
groups x 8; NSAF ``0`` -> ``NaN``, meaning no PSMs in that run), ``raw_linear_complete``
(NSAF > 0 in all 8 runs, contaminants excluded), and ``raw_log`` (log2 of complete,
with a tiny pseudocount -- see :func:`process_nsaf_level`). The PSMs matrix is saved
alongside as ``results/qc_states/nsaf/psms.tsv`` for reference.

Contaminant rule (scientist decision, Stage-3 code review): a no-accession
protein-quants row is a TRUE contaminant only if its Limelight dump group has no
real-accession sibling (``protein_loader.refine_contaminants``); protein/peptide/nsaf
``is_contaminant`` all use this rule.

Also runs ``assess_batch_confounding`` once (batch vs. condition, batch vs.
candidate_pair) and saves it to ``results/stage3/batch_confounding.json``.

No stochastic step (deterministic loaders/transforms; pycombat's empirical-Bayes fit
is deterministic given its input), so there is no seed to record.

Fails loud: any exception writes a FAILED marker into ``--qc-states-dir`` (cleared
at the start of every run) and re-raises, so a crash never leaves a stale-looking
but actually-broken ``qc_states/`` tree with no indication anything went wrong.

Run:
    ./.venv/bin/python scripts/promoted/qc_prep.py \\
        --metadata-file data/metadata.tsv \\
        --mapping-file data/scan-file-search-id-mapping.txt \\
        --protein-quants-file data/protein-quants.tsv \\
        --peptide-quants-file data/peptide-quants.tsv \\
        --protein-limelight-file data/protein-limelight-table-dump.txt \\
        --peptide-limelight-file data/peptide-limelight-table-dump.txt \\
        --samples-file results/metadata/samples.tsv \\
        --cross-check-file results/stage2/limelight_label_map.tsv \\
        --qc-states-dir results/qc_states \\
        --stage3-dir results/stage3
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from common.design import FAILURE_MARKER_NAME  # noqa: E402
from common.hashing import compute_data_version, compute_file_hashes  # noqa: E402
from loaders.batch_correct import (  # noqa: E402
    BatchConfoundingWarning,
    BatchPassthroughWarning,
    ConfoundingReport,
    assess_batch_confounding,
    combat_correct,
)
from loaders.data_loading import Dataset  # noqa: E402
from loaders.dataset_io import load_dataset, save_dataset  # noqa: E402
from loaders.limelight_loader import (  # noqa: E402
    LimelightProteinCounts,
    load_limelight_protein_counts,
)
from loaders.missing_values import handle_missing  # noqa: E402
from loaders.normalize import (  # noqa: E402
    NormalizationMethod,
    log2_transform,
    normalize,
    to_linear,
)
from loaders.peptide_loader import load_peptide_dataset  # noqa: E402
from loaders.protein_loader import PSVID_PREFIX_RE, load_protein_dataset  # noqa: E402
from loaders.samples import read_samples  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "qc-prep",
    "kind": "analysis",
    "provides": [],
    "uses": [
        "common.hashing",
        "common.design",
        "loaders.protein_loader",
        "loaders.peptide_loader",
        "loaders.limelight_loader",
        "loaders.samples",
        "loaders.missing_values",
        "loaders.normalize",
        "loaders.batch_correct",
        "loaders.dataset_io",
        "loaders.data_loading",
    ],
    "seeded_from": None,
    "description": (
        "Stage-3 QC prep-once: materializes raw/complete/normalized/log/batch-"
        "corrected Dataset states for protein and peptide levels, plus a "
        "normalization/ComBat-free nsaf level (protein-group NSAF from "
        "load_limelight_protein_counts), under results/qc_states/<level>/<state>/ "
        "via dataset_io.save_dataset; runs assess_batch_confounding once (batch x "
        "condition, batch x candidate_pair); records normalization params + "
        "feature counts + data_version in a manifest; fails loud (FAILED marker) "
        "on any crash. is_contaminant uses the scientist-confirmed dump-refined "
        "rule (33 protein / 449 peptide) at every level. Batch correction is "
        "batch-label-only; ComBat failure/warnings on the 2-sample B2022 batch are "
        "recorded, not worked around; the data_version stamp is asserted to equal "
        "results/metadata/data_version.json's."
    ),
}

LOGGER = logging.getLogger("qc_prep")

NORMALIZATION_METHOD: NormalizationMethod = "median"
BATCH_COLUMN = "batch"
CONFOUNDING_COVARIATES = ("condition", "candidate_pair")
# NSAF values are ~1e-5..3e-2 (three orders of magnitude smaller than intensity);
# the intensity default pseudocount=1.0 would swamp every real value and destroy
# the log-fold-change scale. This pseudocount only guards the log2(0) edge case
# (unused in practice: raw_linear_complete already excludes NSAF==0 rows) and sits
# 4+ orders of magnitude below the smallest real NSAF value.
NSAF_LOG_PSEUDOCOUNT = 1e-9
# Same six roles, same names, as results/metadata/data_version.json
# (state/METADATA.md) -- this script's data_version must reproduce it exactly.
EXPECTED_DATA_VERSION = (
    "sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74"
)


def subset_features(dataset: Dataset, keep_mask: np.ndarray) -> Dataset:
    """Return a new Dataset keeping only the features where ``keep_mask`` is True.

    Not part of ``loaders`` because nothing else in the project needs a bare
    feature-boolean-mask subset (the missing_values module's filter is the general
    case; this is only for the contaminant exclusion step here).
    """
    keep_mask = np.asarray(keep_mask, dtype=bool)
    if keep_mask.shape != (dataset.abundances.shape[1],):
        raise ValueError(
            f"keep_mask shape {keep_mask.shape} does not match n_features "
            f"({dataset.abundances.shape[1]},)."
        )
    return replace(
        dataset,
        abundances=dataset.abundances[:, keep_mask].copy(),
        feature_names=dataset.feature_names[keep_mask].copy(),
        feature_metadata=dataset.feature_metadata.loc[keep_mask].reset_index(drop=True),
        metadata=dataset.metadata.copy(),
    )


@dataclass
class BatchCorrectionOutcome:
    status: str  # "ok" | "failed"
    warnings_raised: list[str]
    error: str | None
    dataset: Dataset | None


def run_batch_correction(dataset: Dataset) -> BatchCorrectionOutcome:
    """ComBat batch-only correct ``dataset``; record (not swallow) failures/warnings."""
    caught: list[str] = []
    try:
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always", BatchConfoundingWarning)
            warnings.simplefilter("always", BatchPassthroughWarning)
            corrected = combat_correct(dataset, BATCH_COLUMN)
            caught = [f"{w.category.__name__}: {w.message}" for w in recorded]
        return BatchCorrectionOutcome(
            status="ok", warnings_raised=caught, error=None, dataset=corrected
        )
    except Exception as exc:
        return BatchCorrectionOutcome(
            status="failed",
            warnings_raised=caught,
            error=f"{type(exc).__name__}: {exc}",
            dataset=None,
        )


def _dataset_summary(dataset: Dataset) -> dict[str, object]:
    return {
        "n_samples": int(dataset.abundances.shape[0]),
        "n_features": int(dataset.abundances.shape[1]),
        "scale": dataset.scale,
    }


def _confounding_report_to_dict(report: ConfoundingReport) -> dict[str, object]:
    return {
        "batch_column": report.batch_column,
        "covariate_column": report.covariate_column,
        "crosstab": report.crosstab.to_dict(),
        "cramers_v": report.cramers_v,
        "perfectly_confounded": report.perfectly_confounded,
        "message": report.message,
    }


def process_level(
    level: str,
    dataset: Dataset,
    qc_states_dir: Path,
) -> dict[str, object]:
    """Materialize every processing state for one level; return its manifest entry."""
    level_dir = qc_states_dir / level
    states: dict[str, object] = {}

    LOGGER.info("[%s] raw_linear: %s", level, dataset.abundances.shape)
    save_dataset(dataset, level_dir / "raw_linear")
    states["raw_linear"] = _dataset_summary(dataset)

    is_contaminant = dataset.feature_metadata["is_contaminant"].to_numpy(dtype=bool)
    non_contaminant = subset_features(dataset, ~is_contaminant)
    LOGGER.info(
        "[%s] excluded %d contaminant feature(s); %d remain before completeness filter",
        level,
        int(is_contaminant.sum()),
        non_contaminant.abundances.shape[1],
    )

    missing_result = handle_missing(
        non_contaminant, max_missing_fraction=0.0, impute=None
    )
    complete = missing_result.dataset
    LOGGER.info(
        "[%s] raw_linear_complete: %d -> %d feature(s) (dropped %d incomplete)",
        level,
        missing_result.n_features_in,
        missing_result.n_features_out,
        len(missing_result.dropped_features),
    )
    save_dataset(complete, level_dir / "raw_linear_complete")
    states["raw_linear_complete"] = _dataset_summary(complete)

    raw_log = log2_transform(complete)
    save_dataset(raw_log, level_dir / "raw_log")
    states["raw_log"] = _dataset_summary(raw_log)

    normalized_linear = normalize(complete, method=NORMALIZATION_METHOD)
    save_dataset(normalized_linear, level_dir / "normalized_linear")
    states["normalized_linear"] = _dataset_summary(normalized_linear)

    normalized_log = log2_transform(normalized_linear)
    save_dataset(normalized_log, level_dir / "normalized_log")
    states["normalized_log"] = _dataset_summary(normalized_log)

    outcome = run_batch_correction(normalized_log)
    if outcome.status == "ok" and outcome.dataset is not None:
        save_dataset(outcome.dataset, level_dir / "batch_corrected_log")
        states["batch_corrected_log"] = {
            **_dataset_summary(outcome.dataset),
            "warnings": outcome.warnings_raised,
        }
        linear_back = to_linear(outcome.dataset)
        save_dataset(linear_back, level_dir / "batch_corrected_linear")
        states["batch_corrected_linear"] = {
            **_dataset_summary(linear_back),
            "warnings": outcome.warnings_raised,
        }
        for w in outcome.warnings_raised:
            LOGGER.warning("[%s] ComBat: %s", level, w)
    else:
        LOGGER.error("[%s] ComBat batch correction FAILED: %s", level, outcome.error)
        states["batch_corrected_log"] = {"status": "failed", "error": outcome.error}
        states["batch_corrected_linear"] = {
            "status": "failed",
            "error": "upstream batch_corrected_log failed",
        }

    return {
        "states": states,
        "normalization": {
            "method": NORMALIZATION_METHOD,
            "note": (
                "Median normalization (pronoms MedianNormalizer): divide each "
                "sample by its median, rescaled to the mean of medians. Computed "
                "over non-contaminant features quantified in all 8 runs."
            ),
        },
        "missing_value_handling": {
            "max_missing_fraction": missing_result.max_missing_fraction,
            "impute": missing_result.impute,
            "n_features_in": missing_result.n_features_in,
            "n_features_out": missing_result.n_features_out,
            "n_dropped": len(missing_result.dropped_features),
        },
        "n_contaminant_features_excluded": int(is_contaminant.sum()),
    }


def _limelight_feature_metadata(
    limelight: LimelightProteinCounts, protein_feature_metadata: pd.DataFrame
) -> pd.DataFrame:
    """Shared feature_metadata for the nsaf/psm levels (item F/G): protein_group,
    first_member_id, is_contaminant, contaminant_grouped_with_real -- joined from
    the ALREADY dump-refined protein ``feature_metadata`` via each Limelight group's
    first member (never re-derived here -- one contaminant classification, reused
    everywhere). Fails loud if a group's first member has no matching
    protein-quants row (would break the verified 1:1 dump<->quant relation).
    """
    accession_keys = protein_feature_metadata["protein_group"].apply(
        lambda pid: PSVID_PREFIX_RE.sub("", pid)
    )
    protein_group_by_key = dict(
        zip(accession_keys, protein_feature_metadata["protein_group"], strict=True)
    )
    by_protein_group = protein_feature_metadata.set_index("protein_group")

    first_member_keys = [group.split(",")[0] for group in limelight.protein_group]
    missing = [k for k in first_member_keys if k not in protein_group_by_key]
    if missing:
        raise ValueError(
            f"{len(missing)} Limelight group(s) have a first member with no "
            f"matching protein-quants row (breaks the assumed 1:1 dump<->quant "
            f"relation): {missing[:5]}."
        )
    first_member_ids = [protein_group_by_key[k] for k in first_member_keys]
    is_contaminant = by_protein_group.loc[first_member_ids, "is_contaminant"].to_numpy(
        dtype=bool
    )
    grouped_with_real = by_protein_group.loc[
        first_member_ids, "contaminant_grouped_with_real"
    ].to_numpy(dtype=bool)
    return pd.DataFrame(
        {
            "protein_group": np.asarray(limelight.protein_group, dtype=str),
            "first_member_id": first_member_ids,
            "is_contaminant": is_contaminant,
            "contaminant_grouped_with_real": grouped_with_real,
        }
    )


def _build_limelight_dataset(
    limelight: LimelightProteinCounts,
    protein_feature_metadata: pd.DataFrame,
    metadata: pd.DataFrame,
    values: np.ndarray,
) -> Dataset:
    """Build a protein-group x sample ``Dataset`` from a limelight-derived matrix
    (``limelight.nsaf`` or ``limelight.psms``); ``0`` -> ``NaN``, ``scale="linear"``.
    Fails loud if the sample order disagrees with ``metadata``.
    """
    if list(limelight.sample_ids) != list(metadata["sample_id"]):
        raise ValueError(
            "limelight.sample_ids does not match the given metadata's sample_id "
            "order; the two must be row-aligned."
        )
    feature_metadata = _limelight_feature_metadata(limelight, protein_feature_metadata)
    abundances = np.where(values == 0.0, np.nan, values).astype(float)
    return Dataset(
        abundances=abundances,
        feature_names=np.asarray(limelight.protein_group, dtype=str),
        feature_metadata=feature_metadata,
        metadata=metadata.copy(),
        scale="linear",
    )


def build_nsaf_dataset(
    limelight: LimelightProteinCounts,
    protein_feature_metadata: pd.DataFrame,
    metadata: pd.DataFrame,
) -> Dataset:
    """Build the protein-group x sample NSAF ``Dataset`` (item F, scientist request).

    Feature id = the Limelight dump's raw group string (``limelight.protein_group``;
    unique; 1:1 with protein-quants via each group's first member -- verified
    bijection, see ``loaders/limelight_loader.py``). NSAF ``0`` -> ``NaN`` (no PSMs
    quantified in that run). See :func:`_limelight_feature_metadata` for the
    ``feature_metadata`` policy.
    """
    return _build_limelight_dataset(
        limelight, protein_feature_metadata, metadata, limelight.nsaf
    )


def build_psm_dataset(
    limelight: LimelightProteinCounts,
    protein_feature_metadata: pd.DataFrame,
    metadata: pd.DataFrame,
) -> Dataset:
    """Build the protein-group x sample PSM-count ``Dataset`` (item G, scientist
    request). PSM ``0`` -> ``NaN`` (not identified in that run). Same feature id /
    ``feature_metadata`` policy as :func:`build_nsaf_dataset`.
    """
    return _build_limelight_dataset(
        limelight,
        protein_feature_metadata,
        metadata,
        limelight.psms.astype(float),
    )


def process_nsaf_level(
    limelight: LimelightProteinCounts,
    protein_feature_metadata: pd.DataFrame,
    metadata: pd.DataFrame,
    qc_states_dir: Path,
) -> tuple[dict[str, object], np.ndarray]:
    """Materialize the nsaf level's states + the PSMs reference TSV (item F).

    NO normalization and NO batch correction (scientist: NSAF is already
    normalized) -- only ``raw_linear`` / ``raw_linear_complete`` / ``raw_log``.
    Returns ``(manifest_entry, raw_linear_complete.feature_names)`` -- the latter
    is what the psm level's ``raw_linear_complete`` must match exactly (item G).
    """
    dataset = build_nsaf_dataset(limelight, protein_feature_metadata, metadata)
    level_dir = qc_states_dir / "nsaf"
    states: dict[str, object] = {}

    LOGGER.info("[nsaf] raw_linear: %s", dataset.abundances.shape)
    save_dataset(dataset, level_dir / "raw_linear")
    states["raw_linear"] = _dataset_summary(dataset)

    is_contaminant = dataset.feature_metadata["is_contaminant"].to_numpy(dtype=bool)
    non_contaminant = subset_features(dataset, ~is_contaminant)
    missing_result = handle_missing(
        non_contaminant, max_missing_fraction=0.0, impute=None
    )
    complete = missing_result.dataset
    LOGGER.info(
        "[nsaf] raw_linear_complete: %d -> %d feature(s) (dropped %d incomplete)",
        missing_result.n_features_in,
        missing_result.n_features_out,
        len(missing_result.dropped_features),
    )
    save_dataset(complete, level_dir / "raw_linear_complete")
    states["raw_linear_complete"] = _dataset_summary(complete)

    raw_log = log2_transform(complete, pseudocount=NSAF_LOG_PSEUDOCOUNT)
    save_dataset(raw_log, level_dir / "raw_log")
    states["raw_log"] = _dataset_summary(raw_log)

    psms_path = level_dir / "psms.tsv"
    psms_df = pd.DataFrame(
        limelight.psms.T,
        index=pd.Index(limelight.protein_group, name="protein_group"),
        columns=list(limelight.sample_ids),
    )
    psms_df.to_csv(psms_path, sep="\t")
    LOGGER.info("[nsaf] wrote PSMs reference matrix to %s", psms_path)

    manifest_entry: dict[str, object] = {
        "states": states,
        "normalization": None,
        "batch_correction": None,
        "note": (
            "NSAF is already normalized (scientist decision); no normalize/ComBat "
            "step for this level."
        ),
        "missing_value_handling": {
            "max_missing_fraction": missing_result.max_missing_fraction,
            "impute": missing_result.impute,
            "n_features_in": missing_result.n_features_in,
            "n_features_out": missing_result.n_features_out,
            "n_dropped": len(missing_result.dropped_features),
        },
        "n_contaminant_features_excluded": int(is_contaminant.sum()),
        "psms_reference_file": str(psms_path),
    }
    return manifest_entry, complete.feature_names


def process_psm_level(
    limelight: LimelightProteinCounts,
    protein_feature_metadata: pd.DataFrame,
    metadata: pd.DataFrame,
    qc_states_dir: Path,
    nsaf_complete_feature_names: np.ndarray,
) -> dict[str, object]:
    """Materialize the psm level (item G, scientist request): PSM counts used
    AS-IS -- NO normalization, NO batch correction, NO log transform.

    ``raw_linear_complete``'s feature set (>= 1 PSM in all 8 runs, contaminants
    excluded) is asserted -- not assumed -- identical to the nsaf level's, since
    both were filtered independently (PSM>=1 vs. NSAF>0) and should coincide.
    """
    dataset = build_psm_dataset(limelight, protein_feature_metadata, metadata)
    level_dir = qc_states_dir / "psm"
    states: dict[str, object] = {}

    LOGGER.info("[psm] raw_linear: %s", dataset.abundances.shape)
    save_dataset(dataset, level_dir / "raw_linear")
    states["raw_linear"] = _dataset_summary(dataset)

    is_contaminant = dataset.feature_metadata["is_contaminant"].to_numpy(dtype=bool)
    non_contaminant = subset_features(dataset, ~is_contaminant)
    missing_result = handle_missing(
        non_contaminant, max_missing_fraction=0.0, impute=None
    )
    complete = missing_result.dataset
    LOGGER.info(
        "[psm] raw_linear_complete: %d -> %d feature(s) (dropped %d incomplete)",
        missing_result.n_features_in,
        missing_result.n_features_out,
        len(missing_result.dropped_features),
    )
    complete_names = set(complete.feature_names.tolist())
    nsaf_names = set(np.asarray(nsaf_complete_feature_names, dtype=str).tolist())
    if complete_names != nsaf_names:
        only_psm = sorted(complete_names - nsaf_names)[:5]
        only_nsaf = sorted(nsaf_names - complete_names)[:5]
        raise ValueError(
            f"psm raw_linear_complete feature set ({len(complete_names)}) does not "
            f"match the nsaf level's ({len(nsaf_names)}); only in psm: {only_psm}; "
            f"only in nsaf: {only_nsaf}. Expected PSM >= 1 in all 8 runs to "
            f"coincide exactly with NSAF > 0 in all 8 runs."
        )
    save_dataset(complete, level_dir / "raw_linear_complete")
    states["raw_linear_complete"] = _dataset_summary(complete)

    return {
        "states": states,
        "normalization": None,
        "batch_correction": None,
        "log_transform": None,
        "note": (
            "PSM counts used as-is (scientist's explicit choice): no normalize/"
            "ComBat/log2 step for this level. raw_linear_complete's feature set is "
            "asserted identical to the nsaf level's (PSM>=1 in all 8 runs "
            "coincides exactly with NSAF>0 in all 8 runs)."
        ),
        "missing_value_handling": {
            "max_missing_fraction": missing_result.max_missing_fraction,
            "impute": missing_result.impute,
            "n_features_in": missing_result.n_features_in,
            "n_features_out": missing_result.n_features_out,
            "n_dropped": len(missing_result.dropped_features),
        },
        "n_contaminant_features_excluded": int(is_contaminant.sum()),
    }


@dataclass(frozen=True)
class Args:
    metadata_file: Path
    mapping_file: Path
    protein_quants_file: Path
    peptide_quants_file: Path
    protein_limelight_file: Path
    peptide_limelight_file: Path
    samples_file: Path
    cross_check_file: Path
    qc_states_dir: Path
    stage3_dir: Path
    log_level: str


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-file", type=Path, default=Path("data/metadata.tsv"))
    parser.add_argument(
        "--mapping-file",
        type=Path,
        default=Path("data/scan-file-search-id-mapping.txt"),
    )
    parser.add_argument(
        "--protein-quants-file", type=Path, default=Path("data/protein-quants.tsv")
    )
    parser.add_argument(
        "--peptide-quants-file", type=Path, default=Path("data/peptide-quants.tsv")
    )
    parser.add_argument(
        "--protein-limelight-file",
        type=Path,
        default=Path("data/protein-limelight-table-dump.txt"),
    )
    parser.add_argument(
        "--peptide-limelight-file",
        type=Path,
        default=Path("data/peptide-limelight-table-dump.txt"),
    )
    parser.add_argument(
        "--samples-file", type=Path, default=Path("results/metadata/samples.tsv")
    )
    parser.add_argument(
        "--cross-check-file",
        type=Path,
        default=Path("results/stage2/limelight_label_map.tsv"),
    )
    parser.add_argument("--qc-states-dir", type=Path, default=Path("results/qc_states"))
    parser.add_argument("--stage3-dir", type=Path, default=Path("results/stage3"))
    parser.add_argument("--log-level", type=str.upper, default="INFO")
    ns = parser.parse_args(argv)
    return Args(
        metadata_file=ns.metadata_file,
        mapping_file=ns.mapping_file,
        protein_quants_file=ns.protein_quants_file,
        peptide_quants_file=ns.peptide_quants_file,
        protein_limelight_file=ns.protein_limelight_file,
        peptide_limelight_file=ns.peptide_limelight_file,
        samples_file=ns.samples_file,
        cross_check_file=ns.cross_check_file,
        qc_states_dir=ns.qc_states_dir,
        stage3_dir=ns.stage3_dir,
        log_level=ns.log_level,
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    logging.basicConfig(
        level=args.log_level, format="%(asctime)s %(levelname)s %(message)s"
    )
    LOGGER.info("qc_prep starting")
    LOGGER.info("params: %s", args)

    # Same six roles, same names, as results/metadata/data_version.json.
    data_files: dict[str, Path] = {
        "metadata_file": args.metadata_file,
        "mapping_file": args.mapping_file,
        "protein_quants_file": args.protein_quants_file,
        "peptide_quants_file": args.peptide_quants_file,
        "protein_limelight_file": args.protein_limelight_file,
        "peptide_limelight_file": args.peptide_limelight_file,
    }
    for name, path in data_files.items():
        if not path.is_file():
            raise FileNotFoundError(f"{name} not found: {path}")
    if not args.samples_file.is_file():
        raise FileNotFoundError(f"samples_file not found: {args.samples_file}")

    args.qc_states_dir.mkdir(parents=True, exist_ok=True)
    marker_path = args.qc_states_dir / FAILURE_MARKER_NAME
    marker_path.unlink(missing_ok=True)

    try:
        file_hashes = compute_file_hashes(list(data_files.items()))
        data_version = compute_data_version(file_hashes)
        LOGGER.info("data_version = %s", data_version)
        if data_version != EXPECTED_DATA_VERSION:
            raise ValueError(
                f"data_version {data_version} != results/metadata/data_version.json's "
                f"{EXPECTED_DATA_VERSION}; the raw data/ files have changed since "
                f"Stage 1/2 -- re-run the Stage 1/2 pipeline before Stage-3 QC prep."
            )

        protein_result = load_protein_dataset(
            args.protein_quants_file,
            args.samples_file,
            protein_limelight_file=args.protein_limelight_file,
        )
        contaminant_ids = frozenset(
            protein_result.dataset.feature_metadata.loc[
                protein_result.dataset.feature_metadata["is_contaminant"],
                "protein_group",
            ]
        )
        peptide_result = load_peptide_dataset(
            args.peptide_quants_file,
            args.samples_file,
            contaminant_ids=contaminant_ids,
        )
        limelight = load_limelight_protein_counts(
            args.protein_limelight_file,
            args.protein_quants_file,
            args.samples_file,
            cross_check_file=args.cross_check_file,
        )
        samples = read_samples(args.samples_file)
        nsaf_metadata = samples.copy()
        nsaf_metadata.index = pd.Index(
            samples["sample_id"].to_numpy(), name="sample_id"
        )

        levels: dict[str, object] = {}
        levels["protein"] = process_level(
            "protein", protein_result.dataset, args.qc_states_dir
        )
        levels["peptide"] = process_level(
            "peptide", peptide_result.dataset, args.qc_states_dir
        )
        nsaf_manifest, nsaf_complete_feature_names = process_nsaf_level(
            limelight,
            protein_result.dataset.feature_metadata,
            nsaf_metadata,
            args.qc_states_dir,
        )
        levels["nsaf"] = nsaf_manifest
        levels["psm"] = process_psm_level(
            limelight,
            protein_result.dataset.feature_metadata,
            nsaf_metadata,
            args.qc_states_dir,
            nsaf_complete_feature_names,
        )

        LOGGER.info(
            "assessing batch confounding (once): batch vs %s", CONFOUNDING_COVARIATES
        )
        # Any Dataset with the sample metadata will do (identical 8-sample metadata
        # at every level); the protein complete dataset is used as the reference.
        protein_complete_path = args.qc_states_dir / "protein" / "raw_linear_complete"
        reference = load_dataset(protein_complete_path)
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always", BatchConfoundingWarning)
            reports = assess_batch_confounding(
                reference, BATCH_COLUMN, list(CONFOUNDING_COVARIATES)
            )
            confounding_warnings = [str(w.message) for w in recorded]

        args.stage3_dir.mkdir(parents=True, exist_ok=True)
        confounding_payload: dict[str, Any] = {
            "data_version": data_version,
            "batch_column": BATCH_COLUMN,
            "covariates": list(CONFOUNDING_COVARIATES),
            "reports": [_confounding_report_to_dict(r) for r in reports],
            "warnings_raised": confounding_warnings,
        }
        (args.stage3_dir / "batch_confounding.json").write_text(
            json.dumps(confounding_payload, indent=2, sort_keys=True, default=str),
            encoding="utf-8",
        )
        for report in reports:
            LOGGER.info(
                "confounding batch x %s: V=%.4f perfectly_confounded=%s",
                report.covariate_column,
                report.cramers_v,
                report.perfectly_confounded,
            )

        manifest = {
            "data_version": data_version,
            "expected_data_version": EXPECTED_DATA_VERSION,
            "file_hashes": file_hashes,
            "params": {
                "metadata_file": str(args.metadata_file),
                "mapping_file": str(args.mapping_file),
                "protein_quants_file": str(args.protein_quants_file),
                "peptide_quants_file": str(args.peptide_quants_file),
                "protein_limelight_file": str(args.protein_limelight_file),
                "peptide_limelight_file": str(args.peptide_limelight_file),
                "samples_file": str(args.samples_file),
                "cross_check_file": str(args.cross_check_file),
                "normalization_method": NORMALIZATION_METHOD,
                "batch_column": BATCH_COLUMN,
                "nsaf_log_pseudocount": NSAF_LOG_PSEUDOCOUNT,
            },
            "levels": levels,
        }
        (args.qc_states_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True, default=str),
            encoding="utf-8",
        )
    except Exception as exc:
        marker_path.write_text(
            json.dumps(
                {"failed": True, "error_type": type(exc).__name__, "error": str(exc)},
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        LOGGER.error("qc_prep FAILED: %s", exc)
        raise

    LOGGER.info(
        "wrote qc_states to %s, batch_confounding.json to %s",
        args.qc_states_dir,
        args.stage3_dir,
    )
    LOGGER.info("qc_prep done")


if __name__ == "__main__":
    main()
