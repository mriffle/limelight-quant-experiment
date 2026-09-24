"""Stage 4: LFQ protein abundance vs spectral counting (NSAF, PSM counts) -- agreement.

Question (scientist): how well do FlashLFQ MS1 protein intensities agree with the two
spectral-count abundance measures exported by Limelight, per run and averaged over
runs, where does agreement break down, what explains the disagreement, how do the
measures differ in what they detect, and do they track each other across runs within
a protein?

Quantities (all log2):
  * ``lfq``   -- FlashLFQ protein intensity, median-normalized with the project's
                 per-run factors (recovered from ``protein/normalized_linear`` /
                 ``protein/raw_linear`` and applied to every quantified cell, so the
                 1,801 complete proteins equal ``protein/normalized_log`` exactly);
  * ``psm``   -- Limelight PSM count, raw (unnormalized, as the scientist uses it);
  * ``nsaf``  -- Limelight NSAF as exported (display-rounded; see caveats).
Contaminants (33, ``is_contaminant``) are excluded from every analysis and summarized
separately.

Sections (see ``summary.json`` for the key numbers and the table index):
  1. pairwise agreement (pooled cells, per sample, protein means; Pearson + Spearman,
     protein-cluster bootstrap CIs; all-common-cells and complete-in-all-8 subsets);
  2. LOWESS calibration trends, OLS / inverse-OLS / SMA / orthogonal / Deming slopes
     (bootstrap CIs), dynamic range;
  3. agreement by PSM-count bin and LFQ decile with the zero-truncated Poisson floor,
     residual variance partition, replicate count noise vs Poisson;
  4. residuals of log2 PSM / NSAF from the LFQ trend vs protein properties (relative
     length from the NSAF identity, unique quantified peptides, MBR intensity share,
     shared-peptide PSM share): Spearman + multiple regression, BH;
  5. detection: PSMs without LFQ and LFQ without PSMs, per sample and per protein,
     with the mechanism (shared-only peptides, unquantified peptides, MBR-only);
  6. within-protein across-run correlation with exact per-protein and global
     run-label permutation nulls (BH over proteins).

Deterministic given ``--seed`` (bootstrap / independent-shuffle null only).

Run:
    ./.venv/bin/python scripts/scratch/abundance_agreement.py
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
import time
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis import abundance_agreement as aa  # noqa: E402
from common.hashing import sha256_of_file  # noqa: E402
from loaders.dataset_io import load_dataset  # noqa: E402
from loaders.peptide_loader import load_peptide_dataset  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "abundance-agreement",
    "kind": "analysis",
    "provides": [],
    "uses": [
        "analysis.abundance_agreement",
        "loaders.dataset_io",
        "loaders.peptide_loader",
        "common.hashing",
    ],
    "seeded_from": None,
    "description": (
        "LFQ protein abundance vs NSAF and PSM counts: pooled/per-sample/protein-mean "
        "correlations with bootstrap CIs, LOWESS calibration + errors-in-variables "
        "slopes + dynamic range, agreement by PSM bin / LFQ decile vs the Poisson "
        "floor, residual-vs-property analysis (NSAF-identity length, peptides, MBR, "
        "shared PSMs), detection comparison with mechanism, within-protein across-run "
        "tracking with permutation nulls. Tables + summary.json; no figures."
    ),
}

LOG = logging.getLogger("abundance_agreement")

DATA_FILES: dict[str, str] = {
    "protein_quants_file": "data/protein-quants.tsv",
    "peptide_quants_file": "data/peptide-quants.tsv",
    "peptide_limelight_file": "data/peptide-limelight-table-dump.txt",
}
METHODS: tuple[aa.CorrMethod, ...] = ("spearman", "pearson")
N_TOP_NO_LFQ = 25
N_PERM_PER_BIN = 500  # independent-shuffle null draws per mean-PSM bin (section 6)
# log2 tolerance vs the stored normalized_log state (observed ~5e-6: float noise).
LFQ_STATE_TOL = 1e-4


@dataclass(frozen=True)
class Params:
    """Every knob that changes the numbers."""

    seed: int
    n_boot: int
    n_perm_independent: int
    lowess_frac: float
    lowess_it: int
    psm_bin_edges: tuple[int, ...]
    n_lfq_bins: int
    slope_min_mean_psm: float


@dataclass(frozen=True)
class Data:
    """Aligned matrices over the non-contaminant analysis proteins."""

    ids: np.ndarray  # protein-quants ids
    entry: np.ndarray
    accession: np.ndarray
    samples: list[str]
    batch: np.ndarray
    condition: np.ndarray
    pair: np.ndarray  # candidate_pair per sample
    run_factor_log2: np.ndarray  # log2 median-normalization factor per run
    lfq: np.ndarray  # log2 normalized, NaN = not quantified
    psm: np.ndarray  # counts, NaN = 0 PSMs
    nsaf: np.ndarray  # as exported, NaN = 0 PSMs
    flashlfq_nan: np.ndarray  # bool per protein: literal NaN token in protein-quants


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def verify_data_files(root: Path, manifest: dict[str, object]) -> None:
    """Fail loud if a raw file read here differs from the QC-state manifest."""
    hashes = manifest["file_hashes"]
    if not isinstance(hashes, dict):
        raise TypeError("manifest file_hashes is not a mapping.")
    for key, rel in DATA_FILES.items():
        entry = hashes[key]
        if entry["path"] != rel:
            raise ValueError(f"{key}: manifest path {entry['path']} != {rel}")
        got = sha256_of_file(root / rel)
        if got != entry["sha256"]:
            raise ValueError(f"{rel}: sha256 {got} != manifest {entry['sha256']}")


def flashlfq_nan_ids(protein_quants_file: Path) -> set[str]:
    """Protein ids carrying the literal ``NaN`` token (FlashLFQ 'not quantifiable')."""
    frame = pd.read_csv(protein_quants_file, sep="\t", dtype=str, keep_default_na=False)
    cols = [c for c in frame.columns if c.startswith("Intensity_")]
    if len(cols) != 8:
        raise ValueError(f"Expected 8 intensity columns, got {len(cols)}.")
    has_nan = (frame[cols] == "NaN").any(axis=1)
    return set(frame.loc[has_nan, "Protein Groups"].astype(str))


def load_all(root: Path, qc_root: Path) -> tuple[Data, Data, dict[str, object]]:
    """Load + align LFQ / PSM / NSAF; return (analysis set, contaminants, info)."""
    raw = load_dataset(qc_root / "protein/raw_linear")
    norm_lin = load_dataset(qc_root / "protein/normalized_linear")
    norm_log = load_dataset(qc_root / "protein/normalized_log")
    nsaf_ds = load_dataset(qc_root / "nsaf/raw_linear")
    psm_ds = load_dataset(qc_root / "psm/raw_linear")
    for ds in (norm_lin, norm_log, nsaf_ds, psm_ds):
        if not ds.metadata.index.equals(raw.metadata.index):
            raise ValueError("Sample order differs between processing states.")
    if raw.scale != "linear" or norm_log.scale != "log2":
        raise ValueError("Unexpected scale tags.")

    ids = np.asarray(raw.feature_names, dtype=str)
    col = pd.Index(ids).get_indexer(norm_lin.feature_names)
    if (col < 0).any():
        raise ValueError("normalized_linear features missing from raw_linear.")
    factors = aa.run_normalization_factors(raw.abundances[:, col], norm_lin.abundances)
    lfq = np.log2(raw.abundances * factors[:, None])
    col_log = pd.Index(ids).get_indexer(norm_log.feature_names)
    max_dev = float(np.max(np.abs(lfq[:, col_log] - norm_log.abundances)))
    if max_dev > LFQ_STATE_TOL:
        raise ValueError(f"Recovered LFQ != normalized_log (max dev {max_dev}).")

    if not np.array_equal(nsaf_ds.feature_names, psm_ds.feature_names):
        raise ValueError("NSAF and PSM feature orders differ.")
    first = nsaf_ds.feature_metadata["first_member_id"].astype(str).tolist()
    if first != psm_ds.feature_metadata["first_member_id"].astype(str).tolist():
        raise ValueError("NSAF and PSM first_member_id differ.")
    j = aa.align_by_first_member(ids.tolist(), first)
    psm = psm_ds.abundances[:, j]
    nsaf = nsaf_ds.abundances[:, j]
    if not np.array_equal(np.isnan(psm), np.isnan(nsaf)):
        raise ValueError("PSM and NSAF missingness differ.")
    obs = psm[~np.isnan(psm)]
    if (obs < 1).any() or (obs != np.round(obs)).any():
        raise ValueError("PSM counts must be integers >= 1 where present.")
    cont = raw.feature_metadata["is_contaminant"].to_numpy(dtype=bool)
    cont_spec = nsaf_ds.feature_metadata["is_contaminant"].to_numpy(dtype=bool)[j]
    if not np.array_equal(cont, cont_spec):
        raise ValueError("Contaminant flags disagree between LFQ and spectral rows.")
    nan_ids = flashlfq_nan_ids(root / DATA_FILES["protein_quants_file"])
    if not nan_ids <= set(ids):
        raise ValueError("NaN-token protein ids not in the protein matrix.")
    meta = raw.metadata
    fm = raw.feature_metadata

    def subset(mask: np.ndarray) -> Data:
        return Data(
            ids=ids[mask],
            entry=fm["entry"].astype(str).to_numpy()[mask],
            accession=fm["accession"].astype(str).to_numpy()[mask],
            samples=[str(s) for s in meta.index],
            batch=meta["batch"].astype(str).to_numpy(),
            condition=meta["condition"].astype(str).to_numpy(),
            pair=meta["candidate_pair"].astype(str).to_numpy(),
            run_factor_log2=np.log2(factors),
            lfq=lfq[:, mask],
            psm=psm[:, mask],
            nsaf=nsaf[:, mask],
            flashlfq_nan=np.isin(ids[mask], list(nan_ids)),
        )

    info: dict[str, object] = {
        "n_protein_rows": len(ids),
        "n_contaminants_excluded": int(cont.sum()),
        "join": "Limelight rows -> protein-quants rows via feature_metadata "
        "first_member_id; verified bijection (4,344 <-> 4,344), identical NSAF/PSM "
        "row order and missingness, identical contaminant flags",
        "lfq_normalization_factors": dict(
            zip(meta.index.astype(str), map(float, factors), strict=True)
        ),
        "lfq_matches_normalized_log_max_abs_dev": max_dev,
        "n_flashlfq_nan_token_proteins": len(nan_ids),
    }
    return subset(~cont), subset(cont), info


def load_peptides(
    root: Path, data: Data, cont: Data
) -> tuple[aa.PeptideSummary, object]:
    """Peptide-level summaries (unique quantified peptides, MBR share) per protein."""
    contaminant_ids = frozenset(cont.ids.tolist())
    res = load_peptide_dataset(
        root / DATA_FILES["peptide_quants_file"],
        root / "results/metadata/samples.tsv",
        contaminant_ids=contaminant_ids,
    )
    if [str(s) for s in res.dataset.metadata.index] != data.samples:
        raise ValueError("Peptide sample order differs from protein sample order.")
    groups = res.dataset.feature_metadata["protein_groups"].astype(str).tolist()
    entry_of = dict(zip(data.ids.tolist(), data.entry.tolist(), strict=True))
    entry_of |= {
        i: f"{e}[contaminant]"
        for i, e in zip(cont.ids.tolist(), cont.entry.tolist(), strict=True)
    }
    summary = aa.peptide_protein_summary(
        data.ids.tolist(),
        groups,
        np.asarray(res.dataset.abundances, dtype=float),
        np.asarray(res.detection_type, dtype=str),
        entry_of=entry_of,
        flag_ids=contaminant_ids,
    )
    return summary, res


def _parse_psm_text(col: pd.Series) -> np.ndarray:
    txt = col.astype(str).str.replace(",", "", regex=False).str.strip()
    txt = txt.mask(txt == "", "0")
    if not txt.str.fullmatch(r"\d+").all():
        raise ValueError("Non-integer PSM text in peptide dump.")
    return np.asarray(txt.astype(int).to_numpy(), dtype=float)


def load_shared_psm(root: Path, data: Data, qc_root: Path) -> aa.SharedPsmSummary:
    """PSMs on shared vs group-unique peptides per protein row (peptide dump)."""
    dump = pd.read_csv(
        root / DATA_FILES["peptide_limelight_file"],
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    label_map = pd.read_csv(root / "results/stage2/limelight_label_map.tsv", sep="\t")
    label_of = dict(
        zip(
            label_map["sample_id"].astype(str),
            label_map["label"].astype(str),
            strict=True,
        )
    )
    psms = np.vstack(
        [_parse_psm_text(dump[f"PSMs ({label_of[s]})"]) for s in data.samples]
    )
    nsaf_ds = load_dataset(qc_root / "nsaf/raw_linear")
    group_of = dict(
        zip(
            nsaf_ds.feature_metadata["first_member_id"].astype(str),
            nsaf_ds.feature_metadata["protein_group"].astype(str),
            strict=True,
        )
    )
    # Every Limelight group (incl. contaminants) is a row so shared-ness is judged
    # against the full protein space; analysis rows are the first len(data.ids).
    all_first = list(group_of)
    order = data.ids.tolist() + [
        f for f in all_first if f not in set(data.ids.tolist())
    ]
    members = [[k.strip() for k in group_of[f].split(",")] for f in order]
    peptide_proteins = [
        [k.strip() for k in p.split(",")] for p in dump["Protein(s)"].astype(str)
    ]
    summary = aa.shared_psm_fraction(members, peptide_proteins, psms)
    if summary.n_members_unmatched:
        raise ValueError(
            f"{summary.n_members_unmatched} peptide-dump protein keys match no group."
        )
    n = len(data.ids)
    return aa.SharedPsmSummary(
        psm_total=summary.psm_total[:, :n],
        psm_shared=summary.psm_shared[:, :n],
        n_members_unmatched=summary.n_members_unmatched,
    )


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _num(x: float) -> float | None:
    return None if x is None or not math.isfinite(x) else round(float(x), 6)


def corr_row(
    x: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    method: aa.CorrMethod,
    params: Params,
    rng: np.random.Generator,
) -> dict[str, object]:
    """Correlation + protein-cluster bootstrap CI as a table row."""
    r = aa.correlation(x, y, method)
    sampler = aa.cluster_index_sampler(groups)
    lo, hi = aa.bootstrap_ci(
        lambda idx: aa.correlation(x[idx], y[idx], method), sampler, params.n_boot, rng
    )
    return {
        "method": method,
        "n_cells": len(x),
        "n_proteins": len(np.unique(groups)),
        "r": r,
        "ci_low": lo,
        "ci_high": hi,
    }


def _masked_mean(m: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Column mean of ``m`` over ``mask``; NaN (silently) for all-masked columns."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.asarray(np.nanmean(np.where(mask, m, np.nan), axis=0), dtype=float)


def within_protein_center(m: np.ndarray) -> np.ndarray:
    """Subtract each protein's (column's) mean over its observed runs."""
    return m - np.nanmean(m, axis=0, keepdims=True)


def run_offsets(m: np.ndarray) -> np.ndarray:
    """Per-run median of (value - protein mean) over complete columns."""
    return np.asarray(np.median(m - m.mean(axis=0, keepdims=True), axis=1), dtype=float)


# --------------------------------------------------------------------------- #
# 1. Pairwise agreement
# --------------------------------------------------------------------------- #
def section_correlations(
    d: Data,
    log2_len: np.ndarray,
    params: Params,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, dict[str, object]]:
    lpsm = np.log2(d.psm)
    lnsaf = np.log2(d.nsaf)
    common = ~np.isnan(d.lfq) & ~np.isnan(d.psm)
    complete8 = common.all(axis=0)
    psm_rn = lpsm - run_offsets(lpsm[:, complete8])[:, None]
    pairs: dict[str, tuple[np.ndarray, np.ndarray]] = {
        "lfq_vs_psm": (d.lfq, lpsm),
        "lfq_vs_nsaf": (d.lfq, lnsaf),
        "lfq_vs_psm_runnorm": (d.lfq, psm_rn),
        "psm_vs_nsaf": (lpsm, lnsaf),
        "lfq_per_length_vs_nsaf": (d.lfq - log2_len[None, :], lnsaf),
    }
    pooled_only = {"lfq_vs_psm_runnorm"}
    rows: list[dict[str, object]] = []
    prot_idx = np.broadcast_to(np.arange(d.lfq.shape[1]), d.lfq.shape)
    subsets = {"all_common": common, "complete8": common & complete8[None, :]}
    for sub_name, mask in subsets.items():
        for name, (x, y) in pairs.items():
            for method in METHODS:
                base = {"comparison": name, "subset": sub_name}
                rows.append(
                    base
                    | {"scope": "pooled_cells", "sample_id": "all"}
                    | corr_row(x[mask], y[mask], prot_idx[mask], method, params, rng)
                )
                if name in pooled_only:
                    continue
                for s, sid in enumerate(d.samples):
                    m = mask[s]
                    rows.append(
                        base
                        | {"scope": "per_sample", "sample_id": sid}
                        | corr_row(
                            x[s, m], y[s, m], prot_idx[s, m], method, params, rng
                        )
                    )
                keep = mask.any(axis=0)
                xm = np.nanmean(np.where(mask, x, np.nan)[:, keep], axis=0)
                ym = np.nanmean(np.where(mask, y, np.nan)[:, keep], axis=0)
                rows.append(
                    base
                    | {"scope": "protein_mean", "sample_id": "mean_over_runs"}
                    | corr_row(xm, ym, np.flatnonzero(keep), method, params, rng)
                )
    table = pd.DataFrame(rows)
    batch_of = dict(zip(d.samples, d.batch, strict=True))
    table["batch"] = table["sample_id"].map(batch_of).fillna("")
    summ: dict[str, object] = {}
    for (comp, sub, method), g in table.groupby(["comparison", "subset", "method"]):
        pooled = g[g["scope"] == "pooled_cells"].iloc[0]
        entry: dict[str, object] = {
            "pooled": [_num(pooled.r), _num(pooled.ci_low), _num(pooled.ci_high)],
            "pooled_n_cells": int(pooled.n_cells),
            "pooled_n_proteins": int(pooled.n_proteins),
        }
        ps = g[g["scope"] == "per_sample"]
        if len(ps):
            entry["per_sample_min_max"] = [_num(ps.r.min()), _num(ps.r.max())]
            entry["per_sample_mean_by_batch"] = {
                b: _num(v) for b, v in ps.groupby("batch").r.mean().items()
            }
            entry["per_sample_n_cells_min_max"] = [
                int(ps.n_cells.min()),
                int(ps.n_cells.max()),
            ]
            pm = g[g["scope"] == "protein_mean"].iloc[0]
            entry["protein_mean"] = [_num(pm.r), _num(pm.ci_low), _num(pm.ci_high)]
            entry["protein_mean_n"] = int(pm.n_proteins)
        summ[f"{comp}|{sub}|{method}"] = entry
    return table, {
        "n_common_cells": int(common.sum()),
        "n_complete8_proteins": int(complete8.sum()),
        "psm_run_offsets_log2": dict(
            zip(d.samples, map(_num, run_offsets(lpsm[:, complete8])), strict=True)
        ),
        "correlations": summ,
    }


# --------------------------------------------------------------------------- #
# 2. Trends, slopes, dynamic range
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Trends:
    pooled_psm: aa.LowessTrend  # log2 PSM ~ lfq (all common cells)
    pooled_nsaf: aa.LowessTrend
    pooled_lfq_on_psm: aa.LowessTrend  # lfq ~ log2 PSM
    mean_psm: aa.LowessTrend  # protein means, complete8
    mean_nsaf: aa.LowessTrend


def section_trends(
    d: Data, params: Params, rng: np.random.Generator
) -> tuple[Trends, dict[str, pd.DataFrame], dict[str, object]]:
    lpsm, lnsaf = np.log2(d.psm), np.log2(d.nsaf)
    common = ~np.isnan(d.lfq) & ~np.isnan(d.psm)
    c8 = common.all(axis=0)
    x, yp, yn = d.lfq[common], lpsm[common], lnsaf[common]
    f, it = params.lowess_frac, params.lowess_it
    xm, ypm, ynm = d.lfq[:, c8].mean(0), lpsm[:, c8].mean(0), lnsaf[:, c8].mean(0)
    trends = Trends(
        pooled_psm=aa.fit_lowess(x, yp, f, it),
        pooled_nsaf=aa.fit_lowess(x, yn, f, it),
        pooled_lfq_on_psm=aa.fit_lowess(yp, x, f, it),
        mean_psm=aa.fit_lowess(xm, ypm, f, it),
        mean_nsaf=aa.fit_lowess(xm, ynm, f, it),
    )
    grids = []
    for name, tr, scope, xlab, ylab in (
        ("pooled_psm", trends.pooled_psm, "pooled_cells_all_common", "lfq", "log2_psm"),
        (
            "pooled_nsaf",
            trends.pooled_nsaf,
            "pooled_cells_all_common",
            "lfq",
            "log2_nsaf",
        ),
        (
            "pooled_lfq_on_psm",
            trends.pooled_lfq_on_psm,
            "pooled_cells_all_common",
            "log2_psm",
            "lfq",
        ),
        ("mean_psm", trends.mean_psm, "protein_mean_complete8", "lfq", "log2_psm"),
        ("mean_nsaf", trends.mean_nsaf, "protein_mean_complete8", "lfq", "log2_nsaf"),
    ):
        g = tr.grid()
        g.insert(0, "trend", name)
        g.insert(1, "scope", scope)
        g.insert(2, "x_quantity", xlab)
        g.insert(3, "y_quantity", ylab)
        grids.append(g)
    grid_table = pd.concat(grids, ignore_index=True)

    # Error-variance ratio from replicates (within-protein across-run variance).
    var_w = {
        "lfq": float(np.mean(d.lfq[:, c8].var(axis=0, ddof=1))),
        "psm": float(np.mean(lpsm[:, c8].var(axis=0, ddof=1))),
        "nsaf": float(np.mean(lnsaf[:, c8].var(axis=0, ddof=1))),
    }
    n_runs = d.lfq.shape[0]
    reliab = {
        k: 1
        - (v / n_runs) / float(np.var({"lfq": xm, "psm": ypm, "nsaf": ynm}[k], ddof=1))
        for k, v in var_w.items()
    }
    mean_psm_count = d.psm[:, c8].mean(axis=0)
    slope_rows = []
    prot_idx = np.broadcast_to(np.arange(d.lfq.shape[1]), d.lfq.shape)
    for comp, ymean, ypool in (("psm", ypm, yp), ("nsaf", ynm, yn)):
        lam = var_w[comp] / var_w["lfq"]
        hi_psm = mean_psm_count >= params.slope_min_mean_psm
        # (x, y, cluster ids, replicate error variance of one x value)
        sets = {
            "protein_mean_complete8": (
                xm,
                ymean,
                np.arange(len(xm)),
                var_w["lfq"] / n_runs,
            ),
            f"protein_mean_complete8_meanpsm_ge{params.slope_min_mean_psm:g}": (
                xm[hi_psm],
                ymean[hi_psm],
                np.flatnonzero(hi_psm),
                var_w["lfq"] / n_runs,
            ),
            "pooled_cells_all_common": (x, ypool, prot_idx[common], var_w["lfq"]),
        }
        for set_name, (sx, sy, grp, err_x) in sets.items():
            rel_x = 1 - err_x / float(np.var(sx, ddof=1))

            def with_attenuation(
                e: dict[str, float], r: float = rel_x
            ) -> dict[str, float]:
                # OLS corrected for replicate (measurement) error in x only; protein-
                # specific equation error is left in y, unlike Deming/orthogonal.
                return e | {"ols_attenuation_corrected": e["ols_y_on_x"] / r}

            est = with_attenuation(aa.slope_estimates(sx, sy, lam))
            sampler = aa.cluster_index_sampler(grp)
            draws = []
            for _ in range(params.n_boot):
                idx = sampler(rng)
                draws.append(
                    with_attenuation(aa.slope_estimates(sx[idx], sy[idx], lam))
                )
            for estimator, value in est.items():
                vals = np.array([dr[estimator] for dr in draws])
                lo, hi = np.quantile(vals, [0.025, 0.975])
                slope_rows.append(
                    {
                        "comparison": f"lfq_vs_{comp}",
                        "set": set_name,
                        "estimator": estimator,
                        "slope": value,
                        "ci_low": lo,
                        "ci_high": hi,
                        "intercept": float(sy.mean() - value * sx.mean()),
                        "n_points": len(sx),
                        "n_proteins": len(np.unique(grp)),
                        "deming_lambda": lam,
                        "reliability_x": rel_x,
                        "r_pearson": aa.correlation(sx, sy, "pearson"),
                    }
                )
    slopes = pd.DataFrame(slope_rows)

    dr_rows = []
    for name, v in (("lfq", xm), ("log2_psm", ypm), ("log2_nsaf", ynm)):
        dr_rows.append(
            {
                "scope": "protein_mean_complete8",
                "sample_id": "mean_over_runs",
                "quantity": name,
            }
            | aa.dynamic_range(v)
        )
    for s, sid in enumerate(d.samples):
        m = common[s]
        for name, v in (
            ("lfq", d.lfq[s, m]),
            ("log2_psm", lpsm[s, m]),
            ("log2_nsaf", lnsaf[s, m]),
        ):
            dr_rows.append(
                {"scope": "per_sample_common_cells", "sample_id": sid, "quantity": name}
                | aa.dynamic_range(v)
            )
    dyn = pd.DataFrame(dr_rows)

    def pick(comp: str, set_name: str, est: str) -> list[float | None]:
        r = slopes[
            (slopes.comparison == comp)
            & (slopes.set == set_name)
            & (slopes.estimator == est)
        ].iloc[0]
        return [_num(r.slope), _num(r.ci_low), _num(r.ci_high)]

    summ: dict[str, object] = {
        "lowess": {"frac": params.lowess_frac, "it": params.lowess_it},
        "within_protein_var_log2": {k: _num(v) for k, v in var_w.items()},
        "reliability_of_protein_means": {k: _num(v) for k, v in reliab.items()},
        "slopes": {
            f"{comp}|{st}": {e: pick(comp, st, e) for e in slopes.estimator.unique()}
            for comp in slopes.comparison.unique()
            for st in slopes.set.unique()
        },
        "dynamic_range_protein_mean_complete8": {
            r.quantity: {
                "range_log2": _num(r.range),
                "span_q01_q99_log2": _num(r.span_q01_q99),
                "sd_log2": _num(r.sd),
                "range_log10": _num(r.range_log10),
            }
            for r in dyn[dyn.scope == "protein_mean_complete8"].itertuples()
        },
        "lowess_local_slope_mean_psm_by_lfq": {
            f"{q:.0%}": _num(
                float(
                    np.interp(
                        np.quantile(xm, q),
                        trends.mean_psm.grid().x,
                        trends.mean_psm.grid().local_slope,
                    )
                )
            )
            for q in (0.1, 0.25, 0.5, 0.75, 0.9)
        },
    }
    return (
        trends,
        {"trend_lowess": grid_table, "slopes": slopes, "dynamic_range": dyn},
        summ,
    )


# --------------------------------------------------------------------------- #
# 3. Where agreement breaks down
# --------------------------------------------------------------------------- #
def _zt_floor_var(mean_counts: np.ndarray) -> np.ndarray:
    """Zero-truncated Poisson variance of log2 K for each (zero-truncated) mean."""
    cache: dict[float, float] = {}
    out = np.empty(len(mean_counts))
    for i, m in enumerate(np.round(mean_counts, 4)):
        key = float(m)
        if key not in cache:
            cache[key] = aa.zt_poisson_log2_var(
                aa.zt_poisson_mu_from_mean(max(key, 1.0))
            )
        out[i] = cache[key]
    return out


def _sd_ci(
    values: np.ndarray, groups: np.ndarray, params: Params, rng: np.random.Generator
) -> tuple[float, float, float]:
    sd = float(np.std(values, ddof=1)) if len(values) > 1 else math.nan
    lo, hi = aa.bootstrap_ci(
        lambda idx: float(np.std(values[idx], ddof=1)),
        aa.cluster_index_sampler(groups),
        params.n_boot,
        rng,
    )
    return sd, lo, hi


def section_breakdown(
    d: Data, trends: Trends, params: Params, rng: np.random.Generator
) -> tuple[dict[str, pd.DataFrame], dict[str, object], dict[str, np.ndarray]]:
    lpsm, lnsaf = np.log2(d.psm), np.log2(d.nsaf)
    common = ~np.isnan(d.lfq) & ~np.isnan(d.psm)
    prot = np.broadcast_to(np.arange(d.lfq.shape[1]), d.lfq.shape)
    x, yp, yn, g = d.lfq[common], lpsm[common], lnsaf[common], prot[common]
    counts = d.psm[common]
    resid_psm = yp - trends.pooled_psm(x)
    resid_nsaf = yn - trends.pooled_nsaf(x)
    resid_lfq = x - trends.pooled_lfq_on_psm(yp)
    # Protein-level expected count: zero-truncated mean over the protein's common cells.
    prot_mean_count = _masked_mean(d.psm, common)
    floor_var = _zt_floor_var(prot_mean_count[g])

    bins = aa.psm_bin_labels(counts, params.psm_bin_edges)
    bin_order = list(
        dict.fromkeys(
            aa.psm_bin_labels(
                np.array(params.psm_bin_edges, float), params.psm_bin_edges
            )
        )
    )
    rows = []
    for b in bin_order:
        m = bins == b
        sd, lo, hi = _sd_ci(resid_lfq[m], g[m], params, rng)
        rho = (
            corr_row(x[m], yp[m], g[m], "spearman", params, rng)
            if np.ptp(yp[m]) > 0
            else None
        )
        rows.append(
            {
                "psm_bin": b,
                "n_cells": int(m.sum()),
                "n_proteins": len(np.unique(g[m])),
                "median_psm": float(np.median(counts[m])),
                "median_lfq": float(np.median(x[m])),
                "lfq_q10": float(np.quantile(x[m], 0.1)),
                "lfq_q90": float(np.quantile(x[m], 0.9)),
                "lfq_resid_sd": sd,
                "lfq_resid_sd_ci_low": lo,
                "lfq_resid_sd_ci_high": hi,
                "lfq_resid_robust_sd": float(
                    1.4826 * stats.median_abs_deviation(resid_lfq[m])
                ),
                "spearman_lfq_psm_within_bin": None if rho is None else rho["r"],
                "spearman_ci_low": None if rho is None else rho["ci_low"],
                "spearman_ci_high": None if rho is None else rho["ci_high"],
            }
        )
    by_psm = pd.DataFrame(rows)

    edges = np.quantile(x, np.linspace(0, 1, params.n_lfq_bins + 1))
    dec = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, params.n_lfq_bins - 1)
    rows = []
    for k in range(params.n_lfq_bins):
        m = dec == k
        sd_p, lo_p, hi_p = _sd_ci(resid_psm[m], g[m], params, rng)
        sd_n, lo_n, hi_n = _sd_ci(resid_nsaf[m], g[m], params, rng)
        rho_p = corr_row(x[m], yp[m], g[m], "spearman", params, rng)
        rho_n = corr_row(x[m], yn[m], g[m], "spearman", params, rng)
        floor_sd = math.sqrt(float(np.mean(floor_var[m])))
        rows.append(
            {
                "lfq_decile": k + 1,
                "lfq_low": float(edges[k]),
                "lfq_high": float(edges[k + 1]),
                "n_cells": int(m.sum()),
                "n_proteins": len(np.unique(g[m])),
                "median_psm": float(np.median(counts[m])),
                "frac_psm_eq1": float(np.mean(counts[m] == 1)),
                "frac_psm_le3": float(np.mean(counts[m] <= 3)),
                "psm_resid_sd": sd_p,
                "psm_resid_sd_ci_low": lo_p,
                "psm_resid_sd_ci_high": hi_p,
                "psm_poisson_floor_sd": floor_sd,
                "psm_poisson_floor_var_share": floor_sd**2 / sd_p**2,
                "nsaf_resid_sd": sd_n,
                "nsaf_resid_sd_ci_low": lo_n,
                "nsaf_resid_sd_ci_high": hi_n,
                "spearman_lfq_psm": rho_p["r"],
                "spearman_lfq_psm_ci_low": rho_p["ci_low"],
                "spearman_lfq_psm_ci_high": rho_p["ci_high"],
                "spearman_lfq_nsaf": rho_n["r"],
                "spearman_lfq_nsaf_ci_low": rho_n["ci_low"],
                "spearman_lfq_nsaf_ci_high": rho_n["ci_high"],
            }
        )
    by_lfq = pd.DataFrame(rows)

    # Variance partition of the pooled-trend residual on complete8 cells.
    c8 = common.all(axis=0)
    part = {}
    for name, tr, ymat in (
        ("psm", trends.pooled_psm, lpsm),
        ("nsaf", trends.pooled_nsaf, lnsaf),
    ):
        r = ymat[:, c8] - tr(d.lfq[:, c8])
        total = float(r.var(ddof=1))
        between = float(r.mean(axis=0).var(ddof=1))
        within = float(r.var(axis=0, ddof=1).mean())
        part[name] = {
            "total_var": _num(total),
            "between_protein_var": _num(between),
            "within_protein_var": _num(within),
            "between_share": _num(between / (between + within)),
        }

    # Replicate count noise vs zero-truncated Poisson (complete8 proteins).
    lp8 = lpsm[:, c8] - run_offsets(lpsm[:, c8])[:, None]
    lf8 = d.lfq[:, c8] - run_offsets(d.lfq[:, c8])[:, None]
    ln8 = lnsaf[:, c8] - run_offsets(lnsaf[:, c8])[:, None]
    mean_c = d.psm[:, c8].mean(axis=0)
    zt_var = _zt_floor_var(mean_c)
    mbins = aa.psm_bin_labels(mean_c, params.psm_bin_edges)
    rows = []
    for b in bin_order:
        m = mbins == b
        if not m.any():
            continue
        vp, vl, vn = (
            lp8[:, m].var(axis=0, ddof=1),
            lf8[:, m].var(axis=0, ddof=1),
            ln8[:, m].var(axis=0, ddof=1),
        )
        rows.append(
            {
                "mean_psm_bin": b,
                "n_proteins": int(m.sum()),
                "psm_sd_median": float(np.median(np.sqrt(vp))),
                "psm_rms_sd": math.sqrt(float(vp.mean())),
                "poisson_expected_rms_sd": math.sqrt(float(zt_var[m].mean())),
                "psm_overdispersion_var_ratio": float(vp.mean() / zt_var[m].mean())
                if zt_var[m].mean() > 0
                else math.nan,
                "nsaf_sd_median": float(np.median(np.sqrt(vn))),
                "nsaf_rms_sd": math.sqrt(float(vn.mean())),
                "lfq_sd_median": float(np.median(np.sqrt(vl))),
                "lfq_rms_sd": math.sqrt(float(vl.mean())),
            }
        )
    count_noise = pd.DataFrame(rows)
    floor_share_all = float(np.mean(floor_var) / np.var(resid_psm, ddof=1))
    summ: dict[str, object] = {
        "psm_bins": bin_order,
        "pooled_resid_sd": {
            "log2_psm_given_lfq": _num(float(np.std(resid_psm, ddof=1))),
            "log2_nsaf_given_lfq": _num(float(np.std(resid_nsaf, ddof=1))),
            "lfq_given_log2_psm": _num(float(np.std(resid_lfq, ddof=1))),
        },
        "poisson_floor_share_of_pooled_psm_resid_var": _num(floor_share_all),
        "residual_variance_partition_complete8": part,
        "frac_common_cells_psm_eq1": _num(float(np.mean(counts == 1))),
        "frac_common_cells_psm_le3": _num(float(np.mean(counts <= 3))),
    }
    cell_cols = {
        "resid_log2_psm": np.where(
            common, lpsm - trends.pooled_psm(np.nan_to_num(d.lfq)), np.nan
        ),
        "resid_log2_nsaf": np.where(
            common, lnsaf - trends.pooled_nsaf(np.nan_to_num(d.lfq)), np.nan
        ),
        "resid_lfq_given_psm": np.where(
            common, d.lfq - trends.pooled_lfq_on_psm(np.nan_to_num(lpsm)), np.nan
        ),
        "lfq_decile": np.where(
            common,
            np.clip(
                np.searchsorted(edges, np.nan_to_num(d.lfq), side="right"),
                1,
                params.n_lfq_bins,
            ),
            np.nan,
        ),
    }
    return (
        {
            "binned_by_psm": by_psm,
            "binned_by_lfq_decile": by_lfq,
            "count_noise_replicates": count_noise,
        },
        summ,
        cell_cols,
    )


# --------------------------------------------------------------------------- #
# 4. Residual structure vs protein properties
# --------------------------------------------------------------------------- #
PROPERTIES: tuple[str, ...] = (
    "log2_length_rel",
    "log2_n_unique_peptides",
    "mbr_intensity_frac_mean",
    "log2_shared_psm_inflation",
)


def section_residuals(
    d: Data,
    trends: Trends,
    length: aa.LengthEstimate,
    pep: aa.PeptideSummary,
    shared: aa.SharedPsmSummary,
    params: Params,
    rng: np.random.Generator,
) -> tuple[dict[str, pd.DataFrame], dict[str, object]]:
    lpsm, lnsaf = np.log2(d.psm), np.log2(d.nsaf)
    common = ~np.isnan(d.lfq) & ~np.isnan(d.psm)
    c8 = common.all(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mbr_mean = _masked_mean(pep.mbr_intensity_frac, common)
    tot = shared.psm_total.sum(axis=0)
    shared_frac = np.where(
        tot > 0, shared.psm_shared.sum(axis=0) / np.where(tot > 0, tot, 1), np.nan
    )
    n_uq = pep.per_protein["n_unique_quantified"].to_numpy(dtype=float)
    props = pd.DataFrame(
        {
            "protein_id": d.ids,
            "entry": d.entry,
            "complete8": c8,
            "n_common_runs": common.sum(axis=0),
            "mean_lfq": _masked_mean(d.lfq, common),
            "mean_log2_psm": _masked_mean(lpsm, common),
            "mean_log2_nsaf": _masked_mean(lnsaf, common),
            "mean_psm": _masked_mean(d.psm, common),
            "log2_length_rel": length.log2_length,
            "length_method": length.method,
            "length_n_precise_cells": length.n_precise_cells,
            "length_max_abs_dev_precise": length.max_abs_dev_precise,
            "length_interval_low": length.interval_low,
            "length_interval_high": length.interval_high,
            "n_unique_peptides": n_uq,
            "log2_n_unique_peptides": np.log2(np.where(n_uq > 0, n_uq, np.nan)),
            "n_shared_peptides": pep.per_protein["n_shared_rows"].to_numpy(),
            "mbr_intensity_frac_mean": mbr_mean,
            "shared_psm_frac": shared_frac,
            # PSM_total = PSM_unique / (1 - f): the log2 factor by which shared-peptide
            # PSMs (excluded from FlashLFQ protein quant) inflate the count.
            # Undefined (NaN) when every PSM is on a shared peptide (f = 1).
            "log2_shared_psm_inflation": -np.log2(
                np.where(shared_frac < 1, 1 - shared_frac, np.nan)
            ),
        }
    )
    props["resid_log2_psm"] = np.where(
        c8, props.mean_log2_psm - trends.mean_psm(props.mean_lfq.to_numpy()), np.nan
    )
    props["resid_log2_nsaf"] = np.where(
        c8, props.mean_log2_nsaf - trends.mean_nsaf(props.mean_lfq.to_numpy()), np.nan
    )
    sub = props[c8].reset_index(drop=True)
    for c in PROPERTIES:
        if not np.isfinite(sub[c].to_numpy(dtype=float)).all():
            raise ValueError(f"Property {c} is non-finite on the complete8 set.")

    rows = []
    for comp in ("psm", "nsaf"):
        yv = sub[f"resid_log2_{comp}"].to_numpy()
        for prop in PROPERTIES:
            xv = sub[prop].to_numpy()
            res = corr_row(xv, yv, np.arange(len(xv)), "spearman", params, rng)
            p = float(stats.spearmanr(xv, yv).pvalue)
            rows.append(
                {
                    "residual": f"log2_{comp}",
                    "property": prop,
                    "scope": "protein_mean_complete8",
                    "p_value": p,
                }
                | res
            )
        # Within-protein: does a run's MBR share explain that run's residual?
        r_cell = np.log2(d.psm if comp == "psm" else d.nsaf)[:, c8] - (
            trends.pooled_psm if comp == "psm" else trends.pooled_nsaf
        )(d.lfq[:, c8])
        mbr_cell = pep.mbr_intensity_frac[:, c8]
        rc, mc = within_protein_center(r_cell), within_protein_center(mbr_cell)
        ok = np.isfinite(rc) & np.isfinite(mc)
        grp = np.broadcast_to(np.arange(rc.shape[1]), rc.shape)[ok]
        res = corr_row(mc[ok], rc[ok], grp, "spearman", params, rng)
        p = float(stats.spearmanr(mc[ok], rc[ok]).pvalue)
        rows.append(
            {
                "residual": f"log2_{comp}",
                "property": "mbr_intensity_frac_run",
                "scope": "cell_within_protein_centered",
                "p_value": p,
            }
            | res
        )
    assoc = pd.DataFrame(rows)
    assoc["q_bh"] = aa.bh(assoc.p_value.to_numpy())
    assoc["correction"] = f"BH over the {len(assoc)} residual x property Spearman tests"

    # Multiple regression on standardized properties.
    z = (sub[list(PROPERTIES)] - sub[list(PROPERTIES)].mean()) / sub[
        list(PROPERTIES)
    ].std(ddof=1)
    xmat = sm.add_constant(z.to_numpy())
    reg_rows = []
    r2: dict[str, object] = {}
    for comp in ("psm", "nsaf"):
        yv = sub[f"resid_log2_{comp}"].to_numpy()
        fit = sm.OLS(yv, xmat).fit(cov_type="HC3")
        boots = []
        sampler = aa.cluster_index_sampler(np.arange(len(yv)))
        for _ in range(params.n_boot):
            idx = sampler(rng)
            boots.append(np.linalg.lstsq(xmat[idx], yv[idx], rcond=None)[0])
        bmat = np.array(boots)
        full_r2 = float(fit.rsquared)
        drop = {}
        single = {}
        for k, prop in enumerate(PROPERTIES):
            keep = [0] + [kk + 1 for kk in range(len(PROPERTIES)) if kk != k]
            drop[prop] = full_r2 - float(sm.OLS(yv, xmat[:, keep]).fit().rsquared)
            single[prop] = float(sm.OLS(yv, xmat[:, [0, k + 1]]).fit().rsquared)
            lo, hi = np.quantile(bmat[:, k + 1], [0.025, 0.975])
            reg_rows.append(
                {
                    "residual": f"log2_{comp}",
                    "term": prop,
                    "coef_per_sd": float(fit.params[k + 1]),
                    "coef_per_unit": float(fit.params[k + 1])
                    / float(sub[prop].std(ddof=1)),
                    "ci_low_per_unit": lo / float(sub[prop].std(ddof=1)),
                    "ci_high_per_unit": hi / float(sub[prop].std(ddof=1)),
                    "ci_low": lo,
                    "ci_high": hi,
                    "p_value_hc3": float(fit.pvalues[k + 1]),
                    "delta_r2_drop": drop[prop],
                    "r2_alone": single[prop],
                    "n": len(yv),
                    "model_r2": full_r2,
                    "residual_sd": float(np.std(yv, ddof=1)),
                }
            )
        r2[f"log2_{comp}"] = {
            "r2": _num(full_r2),
            "delta_r2_drop": {k: _num(v) for k, v in drop.items()},
            "r2_alone": {k: _num(v) for k, v in single.items()},
        }
    reg = pd.DataFrame(reg_rows)
    reg["q_bh"] = aa.bh(reg.p_value_hc3.to_numpy())
    reg["correction"] = f"BH over the {len(reg)} regression coefficients (HC3 SEs)"

    prop_corr = sub[list(PROPERTIES)].corr(method="spearman")
    methods = pd.Series(length.method).value_counts().to_dict()
    precise = length.method == "precise"
    interval = np.isin(length.method, ["interval", "interval_empty"])
    summ: dict[str, object] = {
        "n_complete8": len(sub),
        "length_estimate": {
            "method_counts": {str(k): int(v) for k, v in methods.items()},
            "run_offsets_log2_S": dict(
                zip(d.samples, map(_num, length.run_offset), strict=True)
            ),
            "precise_max_abs_dev_quantiles_50_99_100": [
                _num(float(v))
                for v in np.nanquantile(
                    length.max_abs_dev_precise[precise], [0.5, 0.99, 1.0]
                )
            ],
            "interval_width_log2_quantiles_50_90_100": [
                _num(float(v))
                for v in np.quantile(
                    (length.interval_high - length.interval_low)[interval],
                    [0.5, 0.9, 1.0],
                )
            ]
            if interval.any()
            else None,
            "relative_length_spread_log2_q05_q95": [
                _num(float(v)) for v in np.nanquantile(length.log2_length, [0.05, 0.95])
            ],
        },
        "property_spearman_matrix_complete8": {
            a: {b: _num(float(prop_corr.loc[a, b])) for b in PROPERTIES}
            for a in PROPERTIES
        },
        "spearman_resid_vs_property": {
            f"{r.residual}|{r.property}|{r.scope}": [
                _num(r.r),
                _num(r.ci_low),
                _num(r.ci_high),
                _num(r.q_bh),
            ]
            for r in assoc.itertuples()
        },
        "regression": r2,
    }
    return {
        "protein_properties": props,
        "residual_property_assoc": assoc,
        "residual_regression": reg,
    }, summ


# --------------------------------------------------------------------------- #
# 5. Detection
# --------------------------------------------------------------------------- #
def section_detection(
    d: Data,
    pep: aa.PeptideSummary,
    shared: aa.SharedPsmSummary,
    params: Params,
    rng: np.random.Generator,
) -> tuple[dict[str, pd.DataFrame], dict[str, object]]:
    has_l, has_p = ~np.isnan(d.lfq), ~np.isnan(d.psm)
    lfq_any = has_l.any(axis=0)
    rows, dist_rows = [], []
    pvals = []
    for s, sid in enumerate(d.samples):
        both = has_l[s] & has_p[s]
        p_only = ~has_l[s] & has_p[s]
        l_only = has_l[s] & ~has_p[s]
        p_both, p_po = d.psm[s, both], d.psm[s, p_only]
        auc = aa.prob_superiority(p_both, p_po)
        lo, hi = aa.bootstrap_ci(
            _auc_stat(p_both, p_po),
            _two_sample_sampler(len(p_both), len(p_po)),
            params.n_boot,
            rng,
        )
        pval = float(stats.mannwhitneyu(p_both, p_po, alternative="two-sided").pvalue)
        pvals.append(pval)
        mbr_only = l_only & (pep.n_unique_msms[s] == 0) & (pep.n_unique_mbr[s] > 0)
        rows.append(
            {
                "sample_id": sid,
                "batch": d.batch[s],
                "condition": d.condition[s],
                "n_proteins": d.lfq.shape[1],
                "n_both": int(both.sum()),
                "n_psm_only": int(p_only.sum()),
                "n_psm_only_never_lfq": int((p_only & ~lfq_any).sum()),
                "n_psm_only_lfq_other_runs": int((p_only & lfq_any).sum()),
                "n_lfq_only": int(l_only.sum()),
                "n_lfq_only_mbr_only": int(mbr_only.sum()),
                "n_neither": int((~has_l[s] & ~has_p[s]).sum()),
                "median_psm_both": float(np.median(p_both)),
                "median_psm_psm_only": float(np.median(p_po)),
                "frac_psm_only_eq1": float(np.mean(p_po == 1)),
                "frac_both_eq1": float(np.mean(p_both == 1)),
                "auc_psm_both_gt_psm_only": auc,
                "auc_ci_low": lo,
                "auc_ci_high": hi,
                "mannwhitney_p": pval,
                "median_lfq_both": float(np.median(d.lfq[s, both])),
                "median_lfq_lfq_only": float(np.median(d.lfq[s, l_only]))
                if l_only.any()
                else math.nan,
            }
        )
        for cls, m in (
            ("both", both),
            ("psm_only_lfq_other_runs", p_only & lfq_any),
            ("psm_only_never_lfq", p_only & ~lfq_any),
        ):
            labels = aa.psm_bin_labels(d.psm[s, m], params.psm_bin_edges)
            for b, n in pd.Series(labels).value_counts().items():
                dist_rows.append(
                    {"sample_id": sid, "class": cls, "psm_bin": b, "n_proteins": int(n)}
                )
    per_sample = pd.DataFrame(rows)
    per_sample["mannwhitney_q_bh"] = aa.bh(per_sample.mannwhitney_p.to_numpy())
    per_sample["correction"] = "BH over the 8 per-sample Mann-Whitney tests"
    dist = pd.DataFrame(dist_rows)

    # Per-protein: never-LFQ vs LFQ-in-any-run, total PSMs.
    tot_psm = np.nansum(d.psm, axis=0)
    has_any_psm = has_p.any(axis=0)
    never = ~lfq_any & has_any_psm
    a_tot, b_tot = tot_psm[lfq_any & has_any_psm], tot_psm[never]
    auc_all = aa.prob_superiority(a_tot, b_tot)
    lo_all, hi_all = aa.bootstrap_ci(
        _auc_stat(a_tot, b_tot),
        _two_sample_sampler(len(a_tot), len(b_tot)),
        params.n_boot,
        rng,
    )
    pp = pep.per_protein
    cat = np.full(len(d.ids), "", dtype=object)
    n_uq_rows = pp["n_unique_rows"].to_numpy()
    n_uq_quant = pp["n_unique_quantified"].to_numpy()
    cat[never & d.flashlfq_nan] = "flashlfq_nan_not_quantifiable"
    rest = never & ~d.flashlfq_nan
    cat[rest & (n_uq_rows == 0)] = "shared_peptides_only"
    cat[rest & (n_uq_rows > 0) & (n_uq_quant == 0)] = "unique_peptides_never_quantified"
    cat[rest & (n_uq_quant > 0)] = "unique_peptide_quantified_but_no_protein_value"
    tot_sh = shared.psm_total.sum(axis=0)
    no_lfq = pd.DataFrame(
        {
            "protein_id": d.ids,
            "entry": d.entry,
            "accession": d.accession,
            "total_psm": tot_psm,
            "max_psm_run": np.nanmax(np.where(has_p, d.psm, 0), axis=0),
            "n_runs_with_psm": has_p.sum(axis=0),
            "mean_nsaf": _masked_mean(d.nsaf, has_p),
            "reason": cat,
            "n_unique_peptide_rows": n_uq_rows,
            "n_unique_quantified": n_uq_quant,
            "n_shared_peptide_rows": pp["n_shared_rows"].to_numpy(),
            "n_shared_quantified": pp["n_shared_quantified"].to_numpy(),
            **{c: pp[c].to_numpy() for c in pp.columns if c.startswith("det_")},
            "shared_psm_frac": np.where(
                tot_sh > 0,
                shared.psm_shared.sum(axis=0) / np.where(tot_sh > 0, tot_sh, 1),
                np.nan,
            ),
            "sharing_partners": pp["sharing_partners"].to_numpy(),
            "shares_with_contaminant_entry": pp["shares_with_flagged"].to_numpy(),
            "flashlfq_nan_token": d.flashlfq_nan,
        }
    )[never].sort_values("total_psm", ascending=False, kind="stable")

    # LFQ present, 0 PSMs in that run.
    s_idx, p_idx = np.nonzero(has_l & ~has_p)
    lfq_only = pd.DataFrame(
        {
            "protein_id": d.ids[p_idx],
            "entry": d.entry[p_idx],
            "sample_id": np.asarray(d.samples)[s_idx],
            "batch": d.batch[s_idx],
            "lfq": d.lfq[s_idx, p_idx],
            "n_unique_peptides_msms_run": pep.n_unique_msms[s_idx, p_idx],
            "n_unique_peptides_mbr_run": pep.n_unique_mbr[s_idx, p_idx],
            "mbr_intensity_frac_run": pep.mbr_intensity_frac[s_idx, p_idx],
            "n_runs_with_psm_any": has_p.sum(axis=0)[p_idx],
        }
    )
    top = no_lfq.head(N_TOP_NO_LFQ)
    summ: dict[str, object] = {
        "n_proteins_psm_never_lfq": int(never.sum()),
        "reason_counts": {
            str(k): int(v) for k, v in no_lfq.reason.value_counts().items()
        },
        "unique_never_quantified_shared_psm_frac_median_q10": [
            _num(float(v))
            for v in np.nanquantile(
                no_lfq.shared_psm_frac[
                    no_lfq.reason == "unique_peptides_never_quantified"
                ].to_numpy(dtype=float),
                [0.5, 0.1],
            )
        ],
        "unique_never_quantified_n_unique_rows_median": _num(
            float(
                no_lfq.n_unique_peptide_rows[
                    no_lfq.reason == "unique_peptides_never_quantified"
                ].median()
            )
        ),
        "n_shared_only_sharing_with_contaminant_entry": int(
            (
                no_lfq.shares_with_contaminant_entry
                & (no_lfq.reason == "shared_peptides_only")
            ).sum()
        ),
        "psm_of_shared_only_sharing_with_contaminant_entry": int(
            no_lfq.total_psm[
                no_lfq.shares_with_contaminant_entry
                & (no_lfq.reason == "shared_peptides_only")
            ].sum()
        ),
        "never_lfq_total_psm_median_max": [
            _num(float(np.median(b_tot))),
            _num(float(b_tot.max())),
        ],
        "lfq_any_total_psm_median": _num(float(np.median(a_tot))),
        "auc_total_psm_lfq_any_gt_never": [_num(auc_all), _num(lo_all), _num(hi_all)],
        "auc_total_psm_mannwhitney_p": _num(
            float(stats.mannwhitneyu(a_tot, b_tot, alternative="two-sided").pvalue)
        ),
        "per_sample_psm_only_range": [
            int(per_sample.n_psm_only.min()),
            int(per_sample.n_psm_only.max()),
        ],
        "per_sample_lfq_only_range": [
            int(per_sample.n_lfq_only.min()),
            int(per_sample.n_lfq_only.max()),
        ],
        "lfq_only_cells_total": len(lfq_only),
        "lfq_only_cells_mbr_only": int(
            (
                (lfq_only.n_unique_peptides_msms_run == 0)
                & (lfq_only.n_unique_peptides_mbr_run > 0)
            ).sum()
        ),
        "top_never_lfq": [
            {
                "entry": r.entry,
                "total_psm": int(r.total_psm),
                "reason": r.reason,
                "partners": r.sharing_partners,
            }
            for r in top.itertuples()
        ],
    }
    return {
        "detection_per_sample": per_sample,
        "detection_psm_distribution": dist,
        "no_lfq_proteins": no_lfq,
        "lfq_without_psm_cells": lfq_only,
    }, summ


def _auc_stat(a: np.ndarray, b: np.ndarray) -> Callable[[np.ndarray], float]:
    """AUC statistic on a stacked two-sample bootstrap index (see below)."""
    n_a = len(a)

    def stat(idx: np.ndarray) -> float:
        return aa.prob_superiority(a[idx[idx < n_a]], b[idx[idx >= n_a] - n_a])

    return stat


def _two_sample_sampler(
    n_a: int, n_b: int
) -> Callable[[np.random.Generator], np.ndarray]:
    """Stratified bootstrap: indices < n_a index sample A, >= n_a index sample B."""

    def draw(rng: np.random.Generator) -> np.ndarray:
        return np.concatenate(
            [rng.integers(0, n_a, n_a), n_a + rng.integers(0, n_b, n_b)]
        )

    return draw


# --------------------------------------------------------------------------- #
# 6. Within-protein tracking across runs
# --------------------------------------------------------------------------- #
def section_within(
    d: Data, params: Params, rng: np.random.Generator
) -> tuple[dict[str, pd.DataFrame], dict[str, object]]:
    common = ~np.isnan(d.lfq) & ~np.isnan(d.psm)
    c8 = common.all(axis=0)
    lfq = d.lfq[:, c8].T  # (n_prot, n_runs)
    lpsm = np.log2(d.psm[:, c8]).T
    variants = {
        "psm_raw": lpsm,
        "psm_runnorm": lpsm - run_offsets(lpsm.T)[None, :],
        "nsaf": np.log2(d.nsaf[:, c8]).T,
    }
    mean_psm = d.psm[:, c8].mean(axis=0)
    bins = aa.psm_bin_labels(mean_psm, params.psm_bin_edges)
    per = pd.DataFrame(
        {
            "protein_id": d.ids[c8],
            "entry": d.entry[c8],
            "mean_psm": mean_psm,
            "mean_psm_bin": bins,
            "lfq_sd_runs": lfq.std(axis=1, ddof=1),
            "log2_psm_sd_runs": lpsm.std(axis=1, ddof=1),
        }
    )
    null_rows, bin_rows = [], []
    summ: dict[str, object] = {}
    for vname, y in variants.items():
        for method in METHODS:
            key = f"{vname}_{method}"
            r = aa.rowwise_corr(lfq, y, method)
            p = aa.exact_rowwise_perm_pvalues(lfq, y, method)
            per[f"r_{key}"] = r
            per[f"p_{key}"] = p
            per[f"q_{key}"] = aa.bh(p)
            obs_mean, gnull = aa.global_perm_mean_corr(lfq, y, method)
            imeans, imeds = aa.independent_perm_null_mean_corr(
                lfq, y, method, params.n_perm_independent, rng
            )
            ok = np.isfinite(r)
            g_p = float(np.mean(gnull >= obs_mean - 1e-12))
            g_p2 = float(np.mean(np.abs(gnull) >= abs(obs_mean) - 1e-12))
            i_p = (1 + int(np.sum(np.abs(imeans) >= abs(obs_mean)))) / (1 + len(imeans))
            q = per[f"q_{key}"].to_numpy()
            row = {
                "comparison": f"lfq_vs_{vname}",
                "method": method,
                "n_proteins": int(ok.sum()),
                "n_constant_dropped": int((~ok).sum()),
                "mean_r": float(np.nanmean(r)),
                "median_r": float(np.nanmedian(r)),
                "frac_r_positive": float(np.mean(r[ok] > 0)),
                "global_null_mean": float(gnull.mean()),
                "global_null_q025": float(np.quantile(gnull, 0.025)),
                "global_null_q975": float(np.quantile(gnull, 0.975)),
                "global_perm_p_one_sided": g_p,
                "global_perm_p_two_sided": g_p2,
                "global_n_perm": len(gnull),
                "indep_null_mean_q025": float(np.quantile(imeans, 0.025)),
                "indep_null_mean_q975": float(np.quantile(imeans, 0.975)),
                "indep_null_median_q025": float(np.quantile(imeds, 0.025)),
                "indep_null_median_q975": float(np.quantile(imeds, 0.975)),
                "indep_perm_p_two_sided": i_p,
                "indep_n_perm": len(imeans),
                "n_protein_q_lt_0.05": int(np.sum(q < 0.05)),
                "n_protein_q_lt_0.10": int(np.sum(q < 0.10)),
                "n_protein_p_lt_0.05": int(np.sum(p < 0.05)),
            }
            null_rows.append(row)
            summ[key] = {
                k: (_num(v) if isinstance(v, float) else v)
                for k, v in row.items()
                if k not in ("comparison", "method")
            }
            for b in dict.fromkeys(
                aa.psm_bin_labels(
                    np.array(params.psm_bin_edges, float), params.psm_bin_edges
                )
            ):
                m = (bins == b) & ok
                if m.sum() < 3:
                    continue
                zx = lfq[m]
                zy = y[m]
                _, bm = aa.independent_perm_null_mean_corr(
                    zx, zy, method, N_PERM_PER_BIN, rng
                )
                bin_rows.append(
                    {
                        "comparison": f"lfq_vs_{vname}",
                        "method": method,
                        "mean_psm_bin": b,
                        "n_proteins": int(m.sum()),
                        "median_r": float(np.median(r[m])),
                        "mean_r": float(np.mean(r[m])),
                        "indep_null_median_q025": float(np.quantile(bm, 0.025)),
                        "indep_null_median_q975": float(np.quantile(bm, 0.975)),
                        "median_lfq_sd_runs": float(np.median(per.lfq_sd_runs[m])),
                        "median_log2_psm_sd_runs": float(
                            np.median(per.log2_psm_sd_runs[m])
                        ),
                    }
                )
    # Decomposition by candidate pair: between-pair means vs within-pair differences.
    pairs = sorted(set(d.pair.tolist()))
    ctrl, trt = [], []
    for pr_ in pairs:
        idx = np.flatnonzero(d.pair == pr_)
        if len(idx) != 2 or set(d.condition[idx]) != {"control", "raloxifene-d0"}:
            raise ValueError(f"Pair {pr_} is not one control + one raloxifene run.")
        ctrl.append(int(idx[d.condition[idx] == "control"][0]))
        trt.append(int(idx[d.condition[idx] == "raloxifene-d0"][0]))
    ci, ti = np.array(ctrl), np.array(trt)

    def ss_between_share(m: np.ndarray) -> np.ndarray:
        pm = 0.5 * (m[:, ci] + m[:, ti])
        tot = ((m - m.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
        betw = 2 * ((pm - pm.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(tot > 0, betw / tot, np.nan)

    share = {
        "lfq": ss_between_share(lfq),
        "log2_psm": ss_between_share(lpsm),
        "log2_nsaf": ss_between_share(variants["nsaf"]),
    }
    for vname, y in variants.items():
        dx, dy = lfq[:, ti] - lfq[:, ci], y[:, ti] - y[:, ci]
        r_wp = aa.within_pair_corr(dx, dy)
        obs, gnull = aa.global_signflip_mean_corr(dx, dy)
        per[f"r_{vname}_within_pair"] = r_wp
        pmx, pmy = 0.5 * (lfq[:, ci] + lfq[:, ti]), 0.5 * (y[:, ci] + y[:, ti])
        r_bp = aa.rowwise_corr(pmx, pmy, "pearson")
        obs_bp, gnull_bp = aa.global_perm_mean_corr(pmx, pmy, "pearson")
        per[f"r_{vname}_pair_means"] = r_bp
        for label, r_, o_, g_ in (
            ("pearson_within_pair", r_wp, obs, gnull),
            ("pearson_between_pair_means", r_bp, obs_bp, gnull_bp),
        ):
            ok = np.isfinite(r_)
            row = {
                "comparison": f"lfq_vs_{vname}",
                "method": label,
                "n_proteins": int(ok.sum()),
                "n_constant_dropped": int((~ok).sum()),
                "mean_r": float(np.nanmean(r_)),
                "median_r": float(np.nanmedian(r_)),
                "frac_r_positive": float(np.mean(r_[ok] > 0)),
                "global_null_mean": float(g_.mean()),
                "global_null_q025": float(np.quantile(g_, 0.025)),
                "global_null_q975": float(np.quantile(g_, 0.975)),
                "global_perm_p_one_sided": float(np.mean(g_ >= o_ - 1e-12)),
                "global_perm_p_two_sided": float(
                    np.mean(np.abs(g_) >= abs(o_) - 1e-12)
                ),
                "global_n_perm": len(g_),
            }
            null_rows.append(row)
            summ[f"{vname}_{label}"] = {
                k: (_num(v) if isinstance(v, float) else v)
                for k, v in row.items()
                if k not in ("comparison", "method")
            }
    summ["between_pair_ss_share_median"] = {
        k: _num(float(np.nanmedian(v))) for k, v in share.items()
    }
    summ["between_pair_ss_share_expected_if_exchangeable"] = _num(3 / 7)
    per["lfq_between_pair_ss_share"] = share["lfq"]
    per["log2_psm_between_pair_ss_share"] = share["log2_psm"]
    nulls = pd.DataFrame(null_rows)
    nulls["correction_per_protein"] = (
        "exact permutation p (all 8! run relabelings) per protein; "
        "BH over proteins within comparison x method"
    )
    return {
        "within_protein_correlation": per,
        "within_protein_null": nulls,
        "within_protein_by_psm_bin": pd.DataFrame(bin_rows),
    }, summ


# --------------------------------------------------------------------------- #
# FlashLFQ roll-up diagnostic, contaminants, cells table
# --------------------------------------------------------------------------- #
def rollup_check(
    d: Data, pep_res: object, raw_lfq_log2: np.ndarray
) -> dict[str, object]:
    """Is protein LFQ ~ a median polish of its unique peptides' log2 intensities?"""
    ds = pep_res.dataset  # type: ignore[attr-defined]  # PeptideLoadResult (untyped attr access)
    groups = ds.feature_metadata["protein_groups"].astype(str).to_numpy()
    uniq = np.array([";" not in g for g in groups])
    logi = np.where(
        ds.abundances > 0,
        np.log2(np.where(ds.abundances > 0, ds.abundances, 1.0)),
        np.nan,
    )
    owner = pd.Series(np.flatnonzero(uniq)).groupby(groups[uniq]).apply(list)
    complete = ~np.isnan(raw_lfq_log2).any(axis=0)
    prof_dev, offset, n_pep = [], [], []
    for i in np.flatnonzero(complete):
        rows = owner.get(d.ids[i])
        if rows is None:
            continue
        y = logi[:, rows].T
        y = y[~np.isnan(y).all(axis=1)]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            overall, _, col = aa.median_polish(y)
        li = raw_lfq_log2[:, i]
        prof_dev.append(float(np.std((li - li.mean()) - (col - col.mean()))))
        offset.append(float(np.mean(li - overall - col)))
        n_pep.append(len(y))
    off, npep = np.array(offset), np.array(n_pep)
    return {
        "n_proteins": len(prof_dev),
        "profile_dev_sd_quantiles_50_90_99": [
            _num(float(v)) for v in np.quantile(prof_dev, [0.5, 0.9, 0.99])
        ],
        "offset_vs_log2_n_peptides_spearman": _num(
            float(stats.spearmanr(off, np.log2(npep)).statistic)
        ),
        "offset_minus_log2_n_peptides_median_iqr": [
            _num(float(v)) for v in np.quantile(off - np.log2(npep), [0.25, 0.5, 0.75])
        ],
        "note": "Tukey median polish on unique-peptide log2 intensities (MSMS+MBR) per "
        "LFQ-complete protein; profile_dev = SD over runs of the difference between "
        "the protein's log2 LFQ profile and the polish column effects.",
    }


def contaminant_table(c: Data) -> pd.DataFrame:
    has_l, has_p = ~np.isnan(c.lfq), ~np.isnan(c.psm)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return pd.DataFrame(
            {
                "protein_id": c.ids,
                "entry": c.entry,
                "n_runs_lfq": has_l.sum(axis=0),
                "n_runs_psm": has_p.sum(axis=0),
                "total_psm": np.nansum(c.psm, axis=0),
                "mean_lfq": np.nanmean(c.lfq, axis=0),
                "mean_log2_psm": np.nanmean(np.log2(c.psm), axis=0),
            }
        ).sort_values("total_psm", ascending=False)


def cells_table(
    d: Data,
    length: aa.LengthEstimate,
    pep: aa.PeptideSummary,
    shared: aa.SharedPsmSummary,
    trends: Trends,
    extra: dict[str, np.ndarray],
    params: Params,
) -> pd.DataFrame:
    has_l, has_p = ~np.isnan(d.lfq), ~np.isnan(d.psm)
    keep = has_l | has_p
    s_idx, p_idx = np.nonzero(keep)
    common = has_l & has_p
    c8 = common.all(axis=0)
    lpsm = np.log2(d.psm)
    offs = run_offsets(lpsm[:, c8])
    with np.errstate(invalid="ignore", divide="ignore"):
        shared_frac = np.where(
            shared.psm_total > 0, shared.psm_shared / shared.psm_total, np.nan
        )
    table = pd.DataFrame(
        {
            "protein_id": d.ids[p_idx],
            "entry": d.entry[p_idx],
            "sample_id": np.asarray(d.samples)[s_idx],
            "batch": d.batch[s_idx],
            "condition": d.condition[s_idx],
            "has_lfq": has_l[s_idx, p_idx],
            "has_psm": has_p[s_idx, p_idx],
            "in_common": common[s_idx, p_idx],
            "complete8": c8[p_idx],
            "lfq": d.lfq[s_idx, p_idx],
            "psm": d.psm[s_idx, p_idx],
            "log2_psm": lpsm[s_idx, p_idx],
            "log2_psm_runnorm": (lpsm - offs[:, None])[s_idx, p_idx],
            "nsaf": d.nsaf[s_idx, p_idx],
            "log2_nsaf": np.log2(d.nsaf)[s_idx, p_idx],
            "nsaf_precise": (d.nsaf < aa.NSAF_PRECISE_BELOW)[s_idx, p_idx],
            "log2_length_rel": length.log2_length[p_idx],
            "psm_bin": aa.psm_bin_labels(d.psm[s_idx, p_idx], params.psm_bin_edges),
            **{k: v[s_idx, p_idx] for k, v in extra.items()},
            "n_unique_peptides_msms_run": pep.n_unique_msms[s_idx, p_idx],
            "n_unique_peptides_mbr_run": pep.n_unique_mbr[s_idx, p_idx],
            "mbr_intensity_frac_run": pep.mbr_intensity_frac[s_idx, p_idx],
            "shared_psm_frac_run": shared_frac[s_idx, p_idx],
        }
    )
    table["fit_log2_psm_given_lfq"] = np.where(
        table.in_common, trends.pooled_psm(table.lfq.fillna(0).to_numpy()), np.nan
    )
    table["fit_log2_nsaf_given_lfq"] = np.where(
        table.in_common, trends.pooled_nsaf(table.lfq.fillna(0).to_numpy()), np.nan
    )
    return table


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def _write(tables: dict[str, pd.DataFrame], out: Path) -> dict[str, list[str]]:
    index = {}
    for name, t in tables.items():
        t.to_csv(out / f"{name}.tsv", sep="\t", index=False, float_format="%.6g")
        index[f"{name}.tsv"] = [str(c) for c in t.columns]
        LOG.info("wrote %s (%d rows)", name, len(t))
    return index


def main(argv: Sequence[str] | None = None) -> int:
    root = _SCRATCH.parent.parent
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", type=Path, default=root)
    ap.add_argument("--qc-root", type=Path, default=root / "results" / "qc_states")
    ap.add_argument(
        "--out-dir", type=Path, default=root / "results" / "abundance-agreement"
    )
    ap.add_argument("--seed", type=int, default=20260924)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--n-perm-independent", type=int, default=2000)
    ap.add_argument("--lowess-frac", type=float, default=0.3)
    ap.add_argument("--lowess-it", type=int, default=3)
    ap.add_argument("--n-lfq-bins", type=int, default=10)
    ap.add_argument("--slope-min-mean-psm", type=float, default=8.0)
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    params = Params(
        seed=args.seed,
        n_boot=args.n_boot,
        n_perm_independent=args.n_perm_independent,
        lowess_frac=args.lowess_frac,
        lowess_it=args.lowess_it,
        psm_bin_edges=aa.PSM_BIN_EDGES,
        n_lfq_bins=args.n_lfq_bins,
        slope_min_mean_psm=args.slope_min_mean_psm,
    )
    LOG.info("params %s", params)
    rng = np.random.default_rng(params.seed)
    out: Path = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.qc_root / "manifest.json").read_text(encoding="utf-8"))
    verify_data_files(args.root, manifest)
    t0 = time.time()

    d, cont, info = load_all(args.root, args.qc_root)
    LOG.info("analysis proteins %d, contaminants %d", len(d.ids), len(cont.ids))
    pep, pep_res = load_peptides(args.root, d, cont)
    shared = load_shared_psm(args.root, d, args.qc_root)
    length = aa.estimate_relative_length(
        np.where(np.isnan(d.psm), np.nan, d.psm), d.nsaf
    )
    raw_log2 = d.lfq - d.run_factor_log2[:, None]

    tables: dict[str, pd.DataFrame] = {}
    summary: dict[str, object] = {}
    t, s = section_correlations(d, length.log2_length, params, rng)
    tables["correlations"] = t
    summary["1_agreement"] = s
    LOG.info("section 1 done (%.0fs)", time.time() - t0)
    trends, t2, s2 = section_trends(d, params, rng)
    tables |= t2
    summary["2_trend"] = s2
    LOG.info("section 2 done (%.0fs)", time.time() - t0)
    t3, s3, cell_extra = section_breakdown(d, trends, params, rng)
    tables |= t3
    summary["3_breakdown"] = s3
    LOG.info("section 3 done (%.0fs)", time.time() - t0)
    t4, s4 = section_residuals(d, trends, length, pep, shared, params, rng)
    tables |= t4
    summary["4_residuals"] = s4
    LOG.info("section 4 done (%.0fs)", time.time() - t0)
    t5, s5 = section_detection(d, pep, shared, params, rng)
    tables |= t5
    summary["5_detection"] = s5
    summary["5_detection"]["flashlfq_rollup_check"] = rollup_check(d, pep_res, raw_log2)  # type: ignore[index]  # s5 is a dict
    LOG.info("section 5 done (%.0fs)", time.time() - t0)
    t6, s6 = section_within(d, params, rng)
    tables |= t6
    summary["6_within_protein"] = s6
    LOG.info("section 6 done (%.0fs)", time.time() - t0)
    tables["cells"] = cells_table(d, length, pep, shared, trends, cell_extra, params)
    tables["contaminants"] = contaminant_table(cont)
    index = _write(tables, out)

    doc = {
        "question": "Agreement of FlashLFQ protein intensity (median-normalized log2) "
        "with Limelight NSAF and PSM counts (log2), per run and averaged over runs.",
        "data_version": manifest["data_version"],
        "script": "scripts/scratch/abundance_agreement.py",
        "module": "scripts/scratch/analysis/abundance_agreement.py",
        "script_sha256": sha256_of_file(Path(__file__)),
        "module_sha256": sha256_of_file(Path(aa.__file__)),
        "seeded_from": None,
        "params": {
            k: (list(v) if isinstance(v, tuple) else v)
            for k, v in params.__dict__.items()
        },
        "sample_set": "experimental (8 of 8; no QC/pool controls exist)",
        "samples": d.samples,
        "inputs": info
        | {
            "n_analysis_proteins": len(d.ids),
            "shared_psm_members_unmatched": shared.n_members_unmatched,
        },
        "ci_method": "percentile bootstrap resampling proteins (clusters of cells for "
        "pooled statistics), n_boot as in params",
        "corrections": {
            "residual_property_assoc": "BH over all residual x property Spearman tests",
            "residual_regression": "BH over regression coefficients (HC3 p)",
            "detection_per_sample": "BH over the 8 per-sample Mann-Whitney tests",
            "within_protein": "exact per-protein permutation p (8! relabelings), "
            "BH over proteins per comparison x method",
            "correlations": "no tests (descriptive correlations with CIs)",
        },
        **summary,
        "tables": index,
        "runtime_s": round(time.time() - t0, 1),
    }
    (out / "summary.json").write_text(
        json.dumps(doc, indent=2, default=str) + "\n", encoding="utf-8"
    )
    LOG.info("wrote %s", out / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
