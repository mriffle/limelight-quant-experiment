"""Stage-3 QC runner: spectral-count (NSAF, PSM count) versions of the QC families.

Loads the prep-once spectral states written by ``qc_prep.py`` --

  * ``results/qc_states/nsaf/{raw_linear, raw_linear_complete, raw_log}`` -- NSAF
    as-is (no normalization, no batch correction; scientist's choice). ``raw_log`` is
    log2(NSAF + 1e-9) over the 2,168 non-contaminant groups with NSAF > 0 in all 8.
  * ``results/qc_states/psm/{raw_linear, raw_linear_complete}`` -- PSM counts as-is
    (no normalization, no batch correction, no log; scientist's explicit choice).

-- and renders single-state figures through the existing, reviewed family modules in
``qc_figures`` (id_depth, missingness, dynamic_range, abundance_boxplot, cv,
correlation, pca) into the structured family directories under ``figures/qc/``.
Each figure is ``<stem>.{svg,png}`` + ``<stem>.legend.{svg,png}``. Where a module
hard-codes an "intensity" axis label, the label is replaced with the measure's own
(NSAF / PSM count) before saving; no computation is altered.

Scale warnings the modules raise (Pearson / PCA / box plot on linear PSM counts) are
expected here -- the scientist asked for raw counts -- and are recorded, not silenced,
in the provenance JSON ``results/stage3/figure_provenance/spectral.json`` (runner and
module sha256, git HEAD, data_version, input-file sha256s, params, key numbers).

Load -> plot only; deterministic (no stochastic step).

Run (from the project root):
    ./.venv/bin/python scripts/promoted/qc_fig_spectral.py
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
_PROJECT_ROOT = _SCRATCH.parent.parent
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402
from common.figures.figure_io import (  # noqa: E402
    FigureArtifacts,
    publication_style,
    save_figure,
)
from common.hashing import sha256_of_file  # noqa: E402
from loaders.data_loading import Dataset  # noqa: E402
from loaders.dataset_io import load_dataset  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from qc_figures.abundance_boxplot import (  # noqa: E402
    median_range,
    plot_abundance_boxplots,
)
from qc_figures.correlation import plot_sample_correlation  # noqa: E402
from qc_figures.cv import plot_cv_distribution  # noqa: E402
from qc_figures.dynamic_range import plot_dynamic_range  # noqa: E402
from qc_figures.id_depth import plot_id_depth  # noqa: E402
from qc_figures.missingness import (  # noqa: E402
    plot_missingness,
    summarize_completeness,
)
from qc_figures.pca import StatePanel, plot_pca_state_series  # noqa: E402

__script_meta__: dict[str, object] = {
    "task": "qc-fig-spectral",
    "kind": "figure",
    "provides": [],
    "uses": [
        "loaders.dataset_io",
        "common.figures.figure_io",
        "common.hashing",
        "qc_figures.id_depth",
        "qc_figures.missingness",
        "qc_figures.dynamic_range",
        "qc_figures.abundance_boxplot",
        "qc_figures.cv",
        "qc_figures.correlation",
        "qc_figures.pca",
    ],
    "seeded_from": None,
    "description": (
        "Stage-3 QC NSAF + PSM-count versions of the id-depth, missingness, "
        "dynamic-range, abundance-boxplot, cv, sample-correlation and pca families "
        "(single state per figure, as-is data, no normalization/batch correction), "
        "with a provenance JSON."
    ),
}

LOG = logging.getLogger("qc_fig_spectral")

PROCESSING = "as-is, no normalization or batch correction"
CONTAMINANT_GROUP = "contaminant"
LABEL_TOP_N = 5
EXPECTED_N_SAMPLES = 8
CV_STATE_LABEL = "Raw"  # registry NormalizationState value (read-only registry)
CORR_ANNOTATIONS: dict[str, str] = {
    "condition": "condition",
    "batch": "batch",
    "candidate_pair": "candidate pair",
}
PCA_COLORINGS: tuple[tuple[str, str], ...] = (
    ("condition", "Condition"),
    ("batch", "Batch"),
)
MODULES: tuple[str, ...] = (
    "id_depth",
    "missingness",
    "dynamic_range",
    "abundance_boxplot",
    "cv",
    "correlation",
    "pca",
)


@dataclass(frozen=True)
class Measure:
    """One spectral measure: its state directory names and display vocabulary."""

    key: str  # state root under results/qc_states/ and stem token
    name: str  # display name in titles / labels
    log_state: str | None  # log2 complete state, if the measure has one
    log_label: str  # axis label for log2 values of this measure


MEASURES: dict[str, Measure] = {
    "nsaf": Measure(key="nsaf", name="NSAF", log_state="raw_log", log_label="NSAF"),
    "psm": Measure(key="psm", name="PSM count", log_state=None, log_label="PSM count"),
}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def _rel(path: Path | None) -> str | None:
    if path is None:
        return None
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(_PROJECT_ROOT.resolve()))
    except ValueError:
        return str(resolved)


def _git_head(root: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def _state_hashes(state_dir: Path) -> dict[str, str]:
    return {
        str(_rel(p)): sha256_of_file(p)
        for p in sorted(state_dir.iterdir())
        if p.is_file()
    }


def _artifacts(artifacts: FigureArtifacts) -> dict[str, str | None]:
    return {
        "svg": _rel(artifacts.svg),
        "png": _rel(artifacts.png),
        "legend_svg": _rel(artifacts.legend_svg),
        "legend_png": _rel(artifacts.legend_png),
    }


def _relabel(fig: Figure, old: str, new: str) -> None:
    """Replace every axis label equal to ``old`` (keeps font properties); fail loud."""
    n = 0
    for ax in fig.axes:
        for label in (ax.xaxis.label, ax.yaxis.label):
            if label.get_text() == old:
                label.set_text(new)
                n += 1
    if n == 0:
        raise ValueError(f"no axis label {old!r} to replace with {new!r}.")


def _check_samples(ds: Dataset, where: str) -> None:
    """Fail loud unless 8 experimental samples in acquisition order."""
    meta = ds.metadata
    roles = sorted({str(r) for r in meta["sample_role"]})
    if roles != ["experimental"] or len(meta) != EXPECTED_N_SAMPLES:
        raise ValueError(f"{where}: expected 8 experimental samples; got {roles}.")
    keys = list(
        zip(
            meta["batch"].astype(str),
            meta["run_position_within_batch"].astype(int),
            strict=True,
        )
    )
    if keys != sorted(keys):
        raise ValueError(f"{where}: samples not in acquisition order: {keys}.")


def _entry(protein_group: str) -> str:
    """Short display name: first member's last non-empty '|' token."""
    first = protein_group.split(",")[0]
    tokens = [t for t in first.split("|") if t]
    return tokens[-1] if tokens else protein_group


@dataclass
class Rendered:
    """One saved figure plus its provenance/key-number record."""

    stem: str
    family: str
    record: dict[str, Any]


class Runner:
    """Loads states once and renders each family figure into ``figures_root``."""

    def __init__(self, states_root: Path, figures_root: Path, registry: Path, dpi: int):
        self.states_root = states_root
        self.figures_root = figures_root
        self.registry = registry
        self.dpi = dpi
        self._cache: dict[tuple[str, str], Dataset] = {}
        self.inputs: dict[str, dict[str, str]] = {}

    def load(self, measure: str, state: str) -> Dataset:
        key = (measure, state)
        if key not in self._cache:
            state_dir = self.states_root / measure / state
            ds = load_dataset(state_dir)
            _check_samples(ds, f"{measure}/{state}")
            self._cache[key] = ds
            self.inputs[f"{measure}/{state}"] = _state_hashes(state_dir)
            LOG.info(
                "loaded %s/%s %s (%s)", measure, state, ds.abundances.shape, ds.scale
            )
        return self._cache[key]

    def save(
        self,
        family: str,
        stem: str,
        build: Callable[[], tuple[Figure, Figure | None, dict[str, Any]]],
    ) -> Rendered:
        """Build inside publication_style, record warnings, save dual + legend."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with publication_style():
                fig, legend_fig, record = build()
                artifacts = save_figure(
                    fig,
                    self.figures_root / family,
                    stem,
                    legend_fig=legend_fig,
                    dpi=self.dpi,
                )
        msgs = [f"{w.category.__name__}: {w.message}" for w in caught]
        for m in msgs:
            LOG.warning("%s: %s", stem, m)
        record["artifacts"] = _artifacts(artifacts)
        record["warnings"] = msgs
        LOG.info("wrote %s", record["artifacts"]["png"])
        return Rendered(stem=stem, family=family, record=record)

    # ------------------------------------------------------------------ #
    # Families
    # ------------------------------------------------------------------ #

    def id_depth(self) -> Rendered:
        nsaf = self.load("nsaf", "raw_linear")
        psm = self.load("psm", "raw_linear")
        det_n = np.isfinite(nsaf.abundances) & (nsaf.abundances > 0)
        det_p = np.isfinite(psm.abundances) & (psm.abundances > 0)
        if not np.array_equal(
            np.asarray(nsaf.feature_names), np.asarray(psm.feature_names)
        ) or not np.array_equal(det_n, det_p):
            raise ValueError("NSAF>0 and PSM>=1 detection differ; need two panels.")
        # Short panel label: the module draws batch names on the same top band, so a
        # long left-aligned panel title would collide with the first batch's label.
        label = "Protein groups"

        def build() -> tuple[Figure, Figure | None, dict[str, Any]]:
            plot = plot_id_depth(
                {label: psm},
                color_by="condition",
                reference_mask=np.ones(EXPECTED_N_SAMPLES, dtype=bool),
                ylabel={label: "Protein groups identified (PSM ≥ 1)"},
                title=(
                    "Protein groups identified per run (PSM ≥ 1, identical to NSAF > 0)"
                    f"\nNSAF / PSM count, {PROCESSING} · linear"
                ),
                legend_title="Condition",
                annotate_counts=True,
                divider_by="batch",
                registry_path=self.registry,
                persist_colors=False,
            )
            keep = ~psm.feature_metadata["is_contaminant"].to_numpy(dtype=bool)
            noncontam = det_p[:, keep].sum(axis=1)
            ids = [str(s) for s in plot.result.sample_ids]
            rec: dict[str, Any] = {
                "inputs": ["psm/raw_linear", "nsaf/raw_linear"],
                "detection_rule": "finite and > 0 (PSM >= 1 <=> NSAF > 0; verified "
                "identical cell-by-cell over all 4,344 groups)",
                "color_map": plot.color_map,
                "n_features": int(psm.abundances.shape[1]),
                "counts": dict(
                    zip(ids, (int(v) for v in plot.result.counts[label]), strict=True)
                ),
                "counts_non_contaminant": dict(
                    zip(ids, (int(v) for v in noncontam), strict=True)
                ),
                "median": plot.result.reference_median[label],
            }
            return plot.figure, plot.legend_figure, rec

        return self.save("id-depth", "id-depth-nsaf-psm-raw-linear", build)

    def missingness(self, m: Measure) -> Rendered:
        ds = self.load(m.key, "raw_linear")
        n_samples, n_features = ds.abundances.shape

        def build() -> tuple[Figure, Figure | None, dict[str, Any]]:
            plot = plot_missingness(
                ds,
                color_by="batch",
                title=(
                    f"{m.name} missingness · {PROCESSING} · linear "
                    f"({n_samples} runs, {n_features:,} protein groups)"
                ),
                legend_title="batch",
                registry_path=self.registry,
            )
            _relabel(
                plot.figure,
                "mean log2 intensity of detected values",
                f"mean log2 {m.log_label} of detected values",
            )
            res = plot.result
            rec: dict[str, Any] = {
                "inputs": [f"{m.key}/raw_linear"],
                "color_map": plot.color_map,
                "n_features": int(n_features),
                "n_contaminant_features_included": int(
                    ds.feature_metadata["is_contaminant"].sum()
                ),
                "mnar_pearson_r": res.mnar_correlation,
                "n_mnar_features": res.n_mnar_features,
                "features_retained": summarize_completeness(
                    ds, fractions=(1.0, 0.5), group_by="batch"
                ),
            }
            return plot.figure, plot.legend_figure, rec

        return self.save("missingness", f"missingness-{m.key}-raw-linear", build)

    def dynamic_range(self, m: Measure) -> Rendered:
        ds = self.load(m.key, "raw_linear")
        fm = ds.feature_metadata
        names = np.asarray(ds.feature_names).astype(str)
        is_contam = fm["is_contaminant"].to_numpy(dtype=bool)
        entries = np.array([_entry(str(g)) for g in fm["protein_group"]])
        a = np.asarray(ds.abundances, dtype=float)
        detected_any = (np.isfinite(a) & (a > 0)).any(axis=0)
        mark = is_contam & detected_any
        highlight = {
            str(f): str(e) for f, e in zip(names[mark], entries[mark], strict=True)
        }

        def build() -> tuple[Figure, Figure | None, dict[str, Any]]:
            plot = plot_dynamic_range(
                ds,
                highlight_features=highlight,
                highlight_groups=dict.fromkeys(highlight, CONTAMINANT_GROUP),
                label_top_n=LABEL_TOP_N,
                title=(
                    f"{m.name} dynamic range · {PROCESSING}\n"
                    "log2 axis · all 8 runs · contaminants marked"
                ),
                registry_path=self.registry,
                persist_colors=False,
            )
            _relabel(
                plot.figure,
                "log2 intensity (median of detected, a.u.)",
                f"log2 {m.log_label} (median of detected runs)",
            )
            res = plot.result
            ranked = res.feature_names_ranked.astype(str)
            rank_of = {f: i + 1 for i, f in enumerate(ranked)}
            entry_of = dict(zip(names, entries, strict=True))
            contam_ranks = sorted(rank_of[f] for f in highlight)
            k = res.n_features_detected

            def row(i: int) -> dict[str, Any]:
                return {
                    "rank": i + 1,
                    "entry": str(entry_of[ranked[i]]),
                    "contaminant": bool(ranked[i] in highlight),
                    "log2_median": round(float(res.log2_median[i]), 3),
                    "n_detected": int(res.n_detected[i]),
                }

            rec: dict[str, Any] = {
                "inputs": [f"{m.key}/raw_linear"],
                "color_map": plot.color_map,
                "n_features_total": res.n_features_total,
                "n_features_detected": k,
                "dynamic_range_orders": round(res.dynamic_range_orders, 3),
                "log2_median_max": round(float(res.log2_median[0]), 3),
                "log2_median_min": round(float(res.log2_median[-1]), 3),
                "top10": [row(i) for i in range(min(10, k))],
                "contaminants": {
                    "n_flagged_is_contaminant": int(is_contam.sum()),
                    "n_marked": len(highlight),
                    "n_in_top10": sum(r <= 10 for r in contam_ranks),
                    "n_in_top100": sum(r <= 100 for r in contam_ranks),
                    "median_rank": float(np.median(contam_ranks))
                    if contam_ranks
                    else None,
                    "labelled": [row(r - 1) for r in contam_ranks[:LABEL_TOP_N]],
                },
            }
            return plot.figure, plot.legend_figure, rec

        return self.save("dynamic-range", f"dynamic-range-{m.key}-raw-linear", build)

    def boxplot(self, m: Measure) -> Rendered:
        if m.log_state is not None:
            state, scale_tok, scale_txt = m.log_state, "raw-log2", "log2"
            ylabel = f"log2 {m.log_label}"
        else:
            state, scale_tok = "raw_linear_complete", "raw-linear"
            scale_txt = "raw counts, linear"
            ylabel = f"{m.log_label} (raw counts, linear)"
        ds = self.load(m.key, state)
        n_features = ds.abundances.shape[1]
        panel = scale_txt

        def build() -> tuple[Figure, Figure | None, dict[str, Any]]:
            plot = plot_abundance_boxplots(
                {panel: ds},
                categorical_annotations=("condition", "batch"),
                feature_type="protein",
                title=None,
                sample_label_column="sample_id",
                annotate_median_range=True,
                registry_path=self.registry,
                persist_colors=False,
            )
            old = (
                "log2 protein abundance (a.u.)"
                if m.log_state
                else "protein abundance (linear)"
            )
            _relabel(plot.figure, old, ylabel)
            # The module's GridSpec is fixed at top=0.93, which leaves no room for a
            # suptitle above a single panel: add height and lower the grid top.
            fig = plot.figure
            width, height = fig.get_size_inches()
            fig.set_size_inches(width, height + 0.7)
            grid = fig.axes[0].get_subplotspec()
            if grid is None:
                raise ValueError("box-plot panel has no subplotspec.")
            gridspec = grid.get_gridspec()
            if not isinstance(gridspec, GridSpec):
                raise TypeError(f"unexpected box-plot gridspec {type(gridspec)}.")
            gridspec.update(top=0.86)
            fig.suptitle(
                f"{m.name} per sample · {PROCESSING}\n{n_features:,} complete "
                "non-contaminant protein groups",
                fontsize=13,
                weight="bold",
                y=0.98,
            )
            meds = plot.result.medians[panel]
            ids = [str(s) for s in plot.result.sample_ids]
            rec: dict[str, Any] = {
                "inputs": [f"{m.key}/{state}"],
                "panel_label": panel,
                "color_maps": plot.color_maps,
                "n_features": int(n_features),
                "per_sample_median": dict(
                    zip(ids, (round(float(v), 6) for v in meds), strict=True)
                ),
                "median_range": round(float(median_range(meds)), 6),
            }
            return plot.figure, plot.legend_figure, rec

        return self.save(
            "abundance-boxplot", f"abundance-boxplot-{m.key}-{scale_tok}", build
        )

    def cv(self, m: Measure) -> Rendered:
        ds = self.load(m.key, "raw_linear_complete")
        n_features = ds.abundances.shape[1]

        def build() -> tuple[Figure, Figure | None, dict[str, Any]]:
            plot = plot_cv_distribution(
                {CV_STATE_LABEL: ds},
                feature_type="protein",
                xlabel=f"{m.name} CV (std / mean), linear scale",
                title=(
                    f"{m.name} CV across 8 samples · {PROCESSING} · linear\n"
                    f"n = {n_features:,} complete protein groups"
                ),
                legend_title=f"{m.name}, as-is (linear)",
                registry_path=self.registry,
                persist_colors=False,
            )
            cvs = plot.result.cvs[CV_STATE_LABEL]
            rec: dict[str, Any] = {
                "inputs": [f"{m.key}/raw_linear_complete"],
                "color_map": plot.color_map,
                "n_features": int(n_features),
                "median_cv": float(plot.result.medians[CV_STATE_LABEL]),
                "cv_quartiles": [
                    round(float(q), 4) for q in np.nanpercentile(cvs, [25, 75])
                ],
                "n_finite_cv": int(np.isfinite(cvs).sum()),
                "cv_definition": "std(ddof=1)/mean per feature over all 8 samples",
            }
            return plot.figure, plot.legend_figure, rec

        return self.save("cv", f"cv-experimental-{m.key}-raw-linear", build)

    def correlation(self, m: Measure) -> Rendered:
        if m.log_state is not None:
            state, scale_tok, scale_txt = m.log_state, "raw-log2", "log2"
        else:
            state, scale_tok = "raw_linear_complete", "raw-linear"
            scale_txt = "raw linear counts (non-log)"
        ds = self.load(m.key, state)
        n_samples, n_features = ds.abundances.shape

        def build() -> tuple[Figure, Figure | None, dict[str, Any]]:
            plot = plot_sample_correlation(
                ds,
                annotations=CORR_ANNOTATIONS,
                method="pearson",
                title=(
                    f"Sample correlation, {m.name} ({n_features:,} protein groups, "
                    f"n = {n_samples})\n{PROCESSING} · {scale_txt}\n"
                    "Pearson · average linkage"
                ),
                registry_path=self.registry,
            )
            res = plot.result
            iu = np.triu_indices(n_samples, k=1)
            pairs = [
                (
                    str(res.sample_ids[i]),
                    str(res.sample_ids[j]),
                    float(res.matrix[i, j]),
                )
                for i, j in zip(iu[0], iu[1], strict=True)
            ]
            lo = min(pairs, key=lambda p: p[2])
            hi = max(pairs, key=lambda p: p[2])
            ids = [str(s) for s in res.sample_ids]
            rec: dict[str, Any] = {
                "inputs": [f"{m.key}/{state}"],
                "color_maps": plot.color_maps,
                "n_features": int(n_features),
                "clustered_order": [ids[i] for i in res.order],
                "off_diagonal_min": {"pair": list(lo[:2]), "r": round(lo[2], 4)},
                "off_diagonal_max": {"pair": list(hi[:2]), "r": round(hi[2], 4)},
                "off_diagonal_median": round(float(np.median(res.off_diagonal())), 4),
                "matrix": {
                    a: {b: round(float(res.matrix[i, j]), 4) for j, b in enumerate(ids)}
                    for i, a in enumerate(ids)
                },
            }
            return plot.figure, plot.legend_figure, rec

        return self.save(
            "sample-correlation", f"sample-correlation-{m.key}-{scale_tok}", build
        )

    def pca(self, m: Measure, coloring: str, legend_title: str) -> Rendered:
        if m.log_state is not None:
            state, scale_tok, scale_txt = m.log_state, "raw-log2", "log2"
        else:
            state, scale_tok, scale_txt = "raw_linear_complete", "raw-linear", "linear"
        full = self.load(m.key, state)
        n_samples, n_total = full.abundances.shape
        # Standardized PCA divides by each feature's SD: features identical in all 8
        # runs carry no information and would divide by zero, so they are dropped
        # (recorded). max == min, not std == 0: the SD of identical log2 values can
        # come out ~1e-15 from rounding while being 0 after re-computation.
        a_full = np.asarray(full.abundances, dtype=float)
        varies = a_full.max(axis=0) > a_full.min(axis=0)
        dropped = [str(f) for f in np.asarray(full.feature_names)[~varies]]
        ds = replace(
            full,
            abundances=np.asarray(full.abundances)[:, varies],
            feature_names=np.asarray(full.feature_names)[varies],
            feature_metadata=full.feature_metadata.loc[varies].reset_index(drop=True),
        )
        n_features = int(varies.sum())

        def build() -> tuple[Figure, Figure | None, dict[str, Any]]:
            plot = plot_pca_state_series(
                [StatePanel(label=f"{PROCESSING} · {scale_txt}", dataset=ds)],
                coloring,
                label_by="sample_id",
                title=(
                    f"PCA of {m.name} by {coloring}\n"
                    f"{n_features:,} of {n_total:,} complete protein groups "
                    f"({len(dropped)} constant dropped) · n = {n_samples}"
                ),
                legend_title=legend_title,
                registry_path=self.registry,
                panel_size=5.4,
            )
            res = next(iter(plot.results.values()))
            meta = ds.metadata
            ids = meta["sample_id"].astype(str).tolist()
            pc1 = res.scores[:, 0]
            rec: dict[str, Any] = {
                "inputs": [f"{m.key}/{state}"],
                "color_map": plot.color_map,
                "n_features": res.n_features,
                "n_complete_features": int(n_total),
                "zero_variance_dropped": dropped,
                "standardize": res.standardized,
                "pc1_pct": round(float(res.explained_variance_ratio[0]) * 100, 2),
                "pc2_pct": round(float(res.explained_variance_ratio[1]) * 100, 2),
                "scores": {
                    sid: [round(float(v), 3) for v in res.scores[i, :2]]
                    for i, sid in enumerate(ids)
                },
                "pc1_mean_by": {
                    col: {
                        str(g): round(float(pc1[meta[col].astype(str) == g].mean()), 3)
                        for g in dict.fromkeys(meta[col].astype(str))
                    }
                    for col in ("condition", "batch")
                },
            }
            return plot.figure, plot.legend_figure, rec

        return self.save("pca", f"pca-by-{coloring}-{m.key}-{scale_tok}", build)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--qc-states-dir", type=Path, default=_PROJECT_ROOT / "results" / "qc_states"
    )
    p.add_argument(
        "--figures-root",
        type=Path,
        default=_PROJECT_ROOT / "figures" / "qc",
        help="Parent of the family directories (figures/qc/<family>/).",
    )
    p.add_argument(
        "--registry", type=Path, default=_PROJECT_ROOT / "state" / "color_registry.json"
    )
    p.add_argument(
        "--provenance",
        type=Path,
        default=_PROJECT_ROOT
        / "results"
        / "stage3"
        / "figure_provenance"
        / "spectral.json",
    )
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args(argv)
    manifest_path: Path = args.qc_states_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    data_version = manifest.get("data_version")
    if not data_version:
        raise ValueError(f"{manifest_path} has no data_version.")

    runner = Runner(args.qc_states_dir, args.figures_root, args.registry, args.dpi)
    rendered: list[Rendered] = [runner.id_depth()]
    for m in MEASURES.values():
        rendered.append(runner.missingness(m))
        rendered.append(runner.dynamic_range(m))
        rendered.append(runner.boxplot(m))
        rendered.append(runner.cv(m))
        rendered.append(runner.correlation(m))
        for coloring, legend_title in PCA_COLORINGS:
            rendered.append(runner.pca(m, coloring, legend_title))

    script = Path(__file__).resolve()
    provenance: dict[str, Any] = {
        "family": "spectral (nsaf + psm) across qc families",
        "scripts": {
            "runner": {"path": _rel(script), "sha256": sha256_of_file(script)},
            "modules": {
                name: {
                    "path": _rel(_SCRATCH / "qc_figures" / f"{name}.py"),
                    "sha256": sha256_of_file(_SCRATCH / "qc_figures" / f"{name}.py"),
                }
                for name in MODULES
            },
        },
        "git_head": _git_head(_PROJECT_ROOT),
        "data_version": data_version,
        "qc_states_manifest": {
            "path": _rel(manifest_path),
            "sha256": sha256_of_file(manifest_path),
        },
        "color_registry_sha256": sha256_of_file(args.registry),
        "params": {
            "processing": PROCESSING,
            "nsaf_log_pseudocount": manifest.get("params", {}).get(
                "nsaf_log_pseudocount"
            ),
            "detection_rule": "finite and > 0",
            "contaminant_highlight": "feature_metadata.is_contaminant (pure "
            "contaminant groups), detected in >= 1 run",
            "label_top_n": LABEL_TOP_N,
            "missingness_color_by": "batch",
            "cv_registry_label": f"NormalizationState/{CV_STATE_LABEL}",
            "correlation": {"method": "pearson", "linkage": "average"},
            "pca": {"standardize": True, "svd_solver": "full", "n_components": 2},
            "dpi": args.dpi,
        },
        "inputs": runner.inputs,
        "figures": {r.stem: {"family": r.family, **r.record} for r in rendered},
    }
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    LOG.info("wrote provenance %s", _rel(args.provenance))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
