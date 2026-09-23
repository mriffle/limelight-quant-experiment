"""Stage-4 figure family ``differential-abundance`` for raloxifene-d0 vs control.

Loads the **saved** ``DifferentialAbundanceResult`` objects written by
``scripts/scratch/de_raloxifene_vs_control.py`` (never re-fits the models) and renders,
into ``figures/analysis/differential-abundance/raloxifene-vs-control/`` (finding 0004):

1. ``0004-volcano-<quantity>-paired`` — volcano of the primary (paired) design, one per
   quantity (protein, peptide, nsaf, psm_log2); the 8 smallest-q features are labelled
   (entry names; peptide = sequence · protein entry) even though none is a hit.
2. ``0004-pvalue-hist-<quantity>-designs`` — raw-p histograms of the three designs
   (paired / batch / unadjusted) overlaid, with Storey π0 per design in the legend.
3. ``0004-pvalue-hist-quantities-paired`` — the four quantities' paired-design raw-p
   histograms as small multiples with π0 on each panel.

Every figure gets a separate ``.legend.{svg,png}``. Before plotting, each result object
is cross-checked against its exported TSV (same features; effect/p/q equal) and against
``summary.json`` (feature count, hit count, π0) — fail loud on any mismatch. Provenance
(script + module sha256, git HEAD + dirty flags, data_version, input sha256s, params,
colors, labelled features, artifact paths) is written to
``results/de/raloxifene-vs-control/figure_provenance.json``.

Deterministic (no stochastic step; textalloc placement is deterministic).

Run (from the project root):
    ./.venv/bin/python scripts/scratch/fig_de_raloxifene_vs_control.py
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.typing import LineStyleType

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"
for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis.differential_abundance import (  # noqa: E402
    DifferentialAbundanceResult,
)
from analysis.result_io import load_result  # noqa: E402
from analysis_figures import pvalue_hist as ph  # noqa: E402
from analysis_figures import volcano as vo  # noqa: E402
from common.figures.figure_io import (  # noqa: E402
    FigureArtifacts,
    publication_style,
    save_figure,
)
from common.hashing import sha256_of_file  # noqa: E402

__script_meta__: dict[str, object] = {
    "kind": "script",
    "uses": [
        "analysis_figures.volcano",
        "analysis_figures.pvalue_hist",
        "analysis.result_io",
        "analysis.differential_abundance",
        "common.figures.figure_io",
        "common.hashing",
    ],
    "description": (
        "Stage-4 differential-abundance figures for raloxifene-d0 vs control "
        "(finding 0004): paired-design volcano per quantity, per-quantity raw-p "
        "histograms overlaying the three designs, and a paired-design raw-p small-"
        "multiples overview. Loads cached result objects; no re-fitting."
    ),
}

FINDING_ID = "0004"
FDR = 0.05
ANNOTATE_TOP = 8
N_BINS = 20
PI0_LAMBDA = 0.5
PRIMARY_DESIGN = "paired"
CONTRAST_TERM = "condition[raloxifene-d0 vs control]"
EFFECT_LABEL = "log2 fold change (raloxifene-d0 / control)"
DIRECTION_LABELS = ("higher in raloxifene-d0", "higher in control")
# Float tolerance for result-object vs TSV cross-checks (TSV is a %.17g round trip).
_TOL = 1e-9


@dataclass(frozen=True)
class Quantity:
    """Display metadata for one quantity."""

    key: str
    noun: str
    unit: str
    processing: str
    panel: str
    # textalloc longest leader (axes fraction); peptide's near-coincident top cluster
    # needs longer leaders to keep them off neighbouring labels.
    label_max_distance: float = 0.28


QUANTITIES: tuple[Quantity, ...] = (
    Quantity(
        "protein",
        "Protein",
        "proteins",
        "LFQ median-normalized log2",
        "Protein · LFQ median-normalized log2",
    ),
    Quantity(
        "peptide",
        "Peptide",
        "peptides",
        "LFQ median-normalized log2",
        "Peptide · LFQ median-normalized log2",
        0.35,
    ),
    Quantity(
        "nsaf",
        "NSAF",
        "protein groups",
        "NSAF as-is log2",
        "NSAF · as-is log2",
    ),
    Quantity(
        "psm_log2",
        "PSM count",
        "protein groups",
        "PSM counts log2, unnormalized",
        "PSM count · log2, unnormalized",
    ),
)

# (design key, display label, line style) in draw / legend order; paired first so it
# takes the first palette slot in every figure.
DESIGNS: tuple[tuple[str, str, LineStyleType], ...] = (
    ("paired", "paired: condition + pair", "-"),
    ("batch", "batch: condition + batch", "--"),
    ("unadjusted", "unadjusted: condition only", "-."),
)
DESIGN_LABEL = {k: label for k, label, _ in DESIGNS}

# Peptide sequences longer than this are abbreviated on the canvas (first 7 + "…" +
# last 5 residues) so labels fit crossing-free; the full label is kept in provenance.
MAX_LABEL_SEQ_LEN = 15
_MOD_RE = re.compile(r"\[([+-]?[0-9.]+)\]$")


# --------------------------------------------------------------------------- #
# Loading + cross-checks
# --------------------------------------------------------------------------- #
def _entry(protein_id: str) -> str:
    """Entry name of one protein id (the last non-empty ``|`` field).

    Handles both id forms present in the data: ``..._sp|P56134|ATPK_HUMAN`` and the
    no-accession form ``..._sp|CATD_HUMAN|`` (-> ``CATD_HUMAN``).
    """
    fields = [f for f in protein_id.split("|") if f]
    if len(fields) < 2:
        raise ValueError(f"unrecognized protein id {protein_id!r}.")
    return fields[-1]


def _group_label(group: str) -> str:
    """First member's entry name, ``+k`` for the k further members (``;``/``,``)."""
    members = [m for m in re.split(r"[;,]", group) if m]
    if not members:
        raise ValueError(f"empty protein group {group!r}.")
    first = _entry(members[0])
    return first + (f" +{len(members) - 1}" if len(members) > 1 else "")


def _abbrev_seq(seq: str) -> str:
    """``IDGNLVVRPYTPISSDDDKGFVDLVIK`` -> ``IDGNLVV…DLVIK`` beyond MAX_LABEL_SEQ_LEN."""
    return seq if len(seq) <= MAX_LABEL_SEQ_LEN else f"{seq[:7]}…{seq[-5:]}"


def _peptide_labels(tsv: pd.DataFrame, *, abbreviate: bool) -> dict[str, str]:
    """``SEQUENCE · ENTRY`` (``+k`` more groups); mod bracket iff sequence shared."""
    shared = set(tsv.loc[tsv["base_sequence"].duplicated(keep=False), "base_sequence"])
    out: dict[str, str] = {}
    for feature, seq, groups in zip(
        tsv["feature"], tsv["base_sequence"], tsv["protein_groups"], strict=True
    ):
        prot = _group_label(str(groups))
        seq_text = _abbrev_seq(str(seq)) if abbreviate else str(seq)
        if str(seq) in shared:
            m = _MOD_RE.search(str(feature))
            if m is None:
                raise ValueError(f"peptide {feature!r}: no trailing mod bracket.")
            seq_text += f"[{float(m.group(1)):+.2f}]"
        out[str(feature)] = f"{seq_text} · {prot}"
    return out


def label_map(
    quantity: str, tsv: pd.DataFrame, *, abbreviate: bool = True
) -> dict[str, str]:
    """Feature id -> display label for one quantity (peptides abbreviated if long)."""
    if quantity == "protein":
        return {
            str(f): str(e) for f, e in zip(tsv["feature"], tsv["entry"], strict=True)
        }
    if quantity == "peptide":
        return _peptide_labels(tsv, abbreviate=abbreviate)
    return {str(f): _group_label(str(f)) for f in tsv["feature"]}


def _check_against_tsv(
    name: str, result: DifferentialAbundanceResult, tsv: pd.DataFrame
) -> None:
    """Result-object contrast rows must equal the exported TSV (fail loud)."""
    if tuple(result.contrast_terms) != (CONTRAST_TERM,):
        raise ValueError(f"{name}: contrast terms {result.contrast_terms}.")
    rows = result.table[result.table["term"] == CONTRAST_TERM].set_index("feature")
    t = tsv.set_index("feature")
    if set(rows.index) != set(t.index) or len(rows) != len(t):
        raise ValueError(f"{name}: result features differ from the TSV features.")
    t = t.loc[rows.index]
    for col_r, col_t in (("effect", "log2fc"), ("p", "p"), ("q", "q")):
        diff = np.abs(rows[col_r].to_numpy(float) - t[col_t].to_numpy(float))
        if not np.all(diff <= _TOL):
            raise ValueError(f"{name}: {col_r} differs from TSV by {diff.max():.3g}.")


def _check_against_summary(
    name: str,
    design_summary: Mapping[str, Any],
    n_features: int,
    n_hits: int,
    pi0: float,
) -> None:
    if int(design_summary["n_features_tested"]) != n_features:
        raise ValueError(f"{name}: n_features {n_features} != summary.")
    if int(design_summary["hits"][f"q<{FDR:.2f}"]["total"]) != n_hits:
        raise ValueError(f"{name}: hit count {n_hits} != summary.")
    if abs(float(design_summary[f"storey_pi0_lambda{PI0_LAMBDA}"]) - pi0) > 5e-5:
        raise ValueError(f"{name}: pi0 {pi0} != summary.")


def _contrast_p(result: DifferentialAbundanceResult) -> np.ndarray:
    """Raw p-values of the contrast term (not q)."""
    rows = result.table[result.table["term"] == CONTRAST_TERM]
    return np.asarray(rows["p"].to_numpy(dtype=float), dtype=float)


def _n_hits(result: DifferentialAbundanceResult) -> int:
    """Contrast-term features with BH q <= FDR, counted from the result object."""
    rows = result.table[result.table["term"] == CONTRAST_TERM]
    return int((rows["q"].to_numpy(dtype=float) <= FDR).sum())


def _files_sha(directory: Path) -> dict[str, str]:
    return {
        p.name: sha256_of_file(p) for p in sorted(directory.iterdir()) if p.is_file()
    }


def _git(args: list[str], cwd: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", *args], cwd=cwd, capture_output=True, text=True, check=True
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.stdout.strip()


def _artifacts(a: FigureArtifacts) -> dict[str, str | None]:
    return {
        "svg": str(a.svg),
        "png": str(a.png),
        "legend_svg": str(a.legend_svg) if a.legend_svg else None,
        "legend_png": str(a.legend_png) if a.legend_png else None,
    }


# --------------------------------------------------------------------------- #
# Run
# --------------------------------------------------------------------------- #
def run(
    de_dir: Path,
    output_dir: Path,
    registry_path: Path,
    qc_manifest_path: Path,
    provenance_path: Path,
) -> dict[str, Any]:
    """Render the whole family; return (and write) the provenance record."""
    summary_path = de_dir / "summary.json"
    summary = json.loads(summary_path.read_text())
    qc_manifest = json.loads(qc_manifest_path.read_text())
    data_version = summary["data_version"]
    if qc_manifest.get("data_version") != data_version:
        raise ValueError(
            f"data_version mismatch: summary {data_version} vs qc manifest "
            f"{qc_manifest.get('data_version')}."
        )
    if summary["primary_design"] != PRIMARY_DESIGN:
        raise ValueError(f"primary design is {summary['primary_design']!r}.")

    project_root = Path.cwd()
    script = Path(__file__).resolve()
    modules = {
        "volcano": Path(vo.__file__).resolve(),
        "pvalue_hist": Path(ph.__file__).resolve(),
    }
    tracked = [script, *modules.values()]
    record: dict[str, Any] = {
        "script": {
            "path": str(script.relative_to(project_root)),
            "sha256": sha256_of_file(script),
        },
        "modules": {
            k: {
                "path": str(p.relative_to(project_root)),
                "sha256": sha256_of_file(p),
                "seeded_from": (vo if k == "volcano" else ph).__script_meta__[
                    "seeded_from"
                ],
            }
            for k, p in modules.items()
        },
        "git_head": _git(["rev-parse", "HEAD"], project_root),
        "git_status_of_scripts": _git(
            ["status", "--porcelain", "--", *[str(p) for p in tracked]], project_root
        ),
        "data_version": data_version,
        "summary_json": {
            "path": str(summary_path),
            "sha256": sha256_of_file(summary_path),
        },
        "qc_states_manifest": {
            "path": str(qc_manifest_path),
            "sha256": sha256_of_file(qc_manifest_path),
        },
        "color_registry": str(registry_path),
        "params": {
            "fdr": FDR,
            "annotate_top": ANNOTATE_TOP,
            "annotate_scope": "all (smallest q among all tested; ties by raw p)",
            "label_max_distance": {q.key: q.label_max_distance for q in QUANTITIES},
            "n_bins": N_BINS,
            "pi0_lambda": PI0_LAMBDA,
            "contrast_term": CONTRAST_TERM,
            "primary_design": PRIMARY_DESIGN,
            "pvalue_colors_persisted": False,
            "volcano_category": vo.DEFAULT_SIGNIFICANCE_CATEGORY,
        },
        "inputs": {},
        "figures": {},
    }

    results: dict[tuple[str, str], DifferentialAbundanceResult] = {}
    tsvs: dict[tuple[str, str], pd.DataFrame] = {}
    for qty in QUANTITIES:
        for design, _, _ in DESIGNS:
            name = f"{qty.key}_{design}"
            rdir = de_dir / "results" / name
            tsv_path = de_dir / f"{name}.tsv"
            res = load_result(rdir, DifferentialAbundanceResult)
            tsv = pd.read_csv(tsv_path, sep="\t")
            _check_against_tsv(name, res, tsv)
            results[(qty.key, design)] = res
            tsvs[(qty.key, design)] = tsv
            record["inputs"][name] = {
                "result_dir": str(rdir),
                "result_files_sha256": _files_sha(rdir),
                "tsv": {"path": str(tsv_path), "sha256": sha256_of_file(tsv_path)},
            }

    with publication_style():
        # 1. volcanoes (primary design)
        for qty in QUANTITIES:
            res = results[(qty.key, PRIMARY_DESIGN)]
            lmap = label_map(qty.key, tsvs[(qty.key, PRIMARY_DESIGN)])
            lmap_full = label_map(
                qty.key, tsvs[(qty.key, PRIMARY_DESIGN)], abbreviate=False
            )
            n = int((res.table["term"] == CONTRAST_TERM).sum())
            dsum = summary["quantities"][qty.key]["designs"][PRIMARY_DESIGN]
            n_hits = int(dsum["hits"][f"q<{FDR:.2f}"]["total"])
            title = (
                f"{qty.noun}: raloxifene-d0 vs control · "
                f"{DESIGN_LABEL[PRIMARY_DESIGN]}\n"
                f"{qty.processing} · n = {n:,} {qty.unit} · "
                f"{n_hits} hits at q < {FDR:g}"
            )
            plot, features = vo.volcano_from_result(
                res,
                term=CONTRAST_TERM,
                fdr=FDR,
                effect_label=EFFECT_LABEL,
                annotate_top=ANNOTATE_TOP,
                annotate_scope="all",
                label_map=lmap,
                label_max_distance=qty.label_max_distance,
                direction_labels=DIRECTION_LABELS,
                title=title,
                legend_title=f"Significance (BH, n = {n:,})",
                registry_path=registry_path,
                persist_colors=True,
            )
            if plot.counts.up + plot.counts.down != n_hits:
                raise ValueError(f"{qty.key}: volcano hit count != summary.")
            if len(set(plot.labelled_text)) != len(plot.labelled_text):
                raise ValueError(f"{qty.key}: duplicate labels {plot.labelled_text}.")
            stem = f"{FINDING_ID}-volcano-{qty.key.replace('_', '-')}-{PRIMARY_DESIGN}"
            arts = save_figure(
                plot.figure, output_dir, stem, legend_fig=plot.legend_figure
            )
            rows = res.table[res.table["term"] == CONTRAST_TERM].reset_index(drop=True)
            record["figures"][stem] = {
                "kind": "volcano",
                "quantity": qty.key,
                "design": PRIMARY_DESIGN,
                "n_features": n,
                "counts": {
                    "up": plot.counts.up,
                    "down": plot.counts.down,
                    "ns": plot.counts.ns,
                },
                "color_map": plot.color_map,
                "labelled": [
                    {
                        "feature": features[i],
                        "label": text,
                        "label_full": lmap_full[features[i]],
                        "log2fc": float(rows.loc[i, "effect"]),
                        "p": float(rows.loc[i, "p"]),
                        "q": float(rows.loc[i, "q"]),
                    }
                    for i, text in zip(
                        plot.labelled_index, plot.labelled_text, strict=True
                    )
                ],
                "n_labelled_hits": plot.n_labelled_hits,
                "q_threshold_neg_log10": float(-np.log10(FDR)),
                "label_box_top_neg_log10_q": plot.label_box_top,
                "label_layout_conflicts": plot.label_conflicts,
                "min_q": float(rows["q"].min()),
                "label_groups_merged_as_coincident": [
                    [plot.labelled_text[k] for k in g]
                    for g in plot.label_groups
                    if len(g) > 1
                ],
                "max_neg_log10_q": float(-np.log10(rows["q"].min())),
                "effect_range": [
                    float(rows["effect"].min()),
                    float(rows["effect"].max()),
                ],
                "artifacts": _artifacts(arts),
            }

        # 2. per-quantity design overlays
        for qty in QUANTITIES:
            data = {
                label: _contrast_p(results[(qty.key, d)]) for d, label, _ in DESIGNS
            }
            n = int(next(iter(data.values())).size)
            hplot = ph.plot_pvalue_histogram(
                data,
                n_bins=N_BINS,
                pi0_lambda=PI0_LAMBDA,
                linestyles={label: ls for _, label, ls in DESIGNS},
                title=(
                    f"{qty.noun}: raw p-value by design · raloxifene-d0 vs control\n"
                    f"{qty.processing} · n = {n:,} {qty.unit}"
                ),
                legend_title="Design (moderated t)",
                registry_path=registry_path,
                persist_colors=False,
            )
            for d, label, _ in DESIGNS:
                dsum = summary["quantities"][qty.key]["designs"][d]
                _check_against_summary(
                    f"{qty.key}_{d}",
                    dsum,
                    hplot.result.counts[label],
                    _n_hits(results[(qty.key, d)]),
                    hplot.result.pi0[label],
                )
            stem = f"{FINDING_ID}-pvalue-hist-{qty.key.replace('_', '-')}-designs"
            arts = save_figure(
                hplot.figure, output_dir, stem, legend_fig=hplot.legend_figure
            )
            record["figures"][stem] = {
                "kind": "pvalue-hist overlay",
                "quantity": qty.key,
                "designs": [d for d, _, _ in DESIGNS],
                "counts": hplot.result.counts,
                "pi0": hplot.result.pi0,
                "pi0_raw_uncapped": hplot.result.pi0_raw,
                "color_map": hplot.color_map,
                "linestyles": {label: ls for _, label, ls in DESIGNS},
                "artifacts": _artifacts(arts),
            }

        # 3. paired-design small multiples across quantities
        panels = {
            qty.panel: _contrast_p(results[(qty.key, PRIMARY_DESIGN)])
            for qty in QUANTITIES
        }
        splot = ph.plot_pvalue_small_multiples(
            panels,
            series_label=DESIGN_LABEL[PRIMARY_DESIGN],
            color_order=tuple(label for _, label, _ in DESIGNS),
            n_bins=N_BINS,
            pi0_lambda=PI0_LAMBDA,
            title=(
                "Raw p-values by quantity · raloxifene-d0 vs control · "
                f"{DESIGN_LABEL[PRIMARY_DESIGN]}"
            ),
            legend_title="Key",
            registry_path=registry_path,
            persist_colors=False,
        )
        for qty, panel in zip(QUANTITIES, panels, strict=True):
            dsum = summary["quantities"][qty.key]["designs"][PRIMARY_DESIGN]
            _check_against_summary(
                f"{qty.key}_{PRIMARY_DESIGN}",
                dsum,
                splot.result.counts[panel],
                _n_hits(results[(qty.key, PRIMARY_DESIGN)]),
                splot.result.pi0[panel],
            )
        stem = f"{FINDING_ID}-pvalue-hist-quantities-{PRIMARY_DESIGN}"
        arts = save_figure(
            splot.figure, output_dir, stem, legend_fig=splot.legend_figure
        )
        record["figures"][stem] = {
            "kind": "pvalue-hist small multiples",
            "design": PRIMARY_DESIGN,
            "panels": list(panels),
            "counts": splot.result.counts,
            "pi0": splot.result.pi0,
            "color_map": splot.color_map,
            "artifacts": _artifacts(arts),
        }

    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    parser.add_argument(
        "--de-dir", type=Path, default=Path("results/de/raloxifene-vs-control")
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figures/analysis/differential-abundance/raloxifene-vs-control"),
    )
    parser.add_argument(
        "--registry-path", type=Path, default=Path("state/color_registry.json")
    )
    parser.add_argument(
        "--qc-manifest-path", type=Path, default=Path("results/qc_states/manifest.json")
    )
    parser.add_argument(
        "--provenance-path",
        type=Path,
        default=Path("results/de/raloxifene-vs-control/figure_provenance.json"),
    )
    args = parser.parse_args(argv)
    record = run(
        args.de_dir,
        args.output_dir,
        args.registry_path,
        args.qc_manifest_path,
        args.provenance_path,
    )
    for stem, fig in record["figures"].items():
        extra = fig.get("counts")
        print(f"{stem}: {extra}")
        if "labelled" in fig:
            print("   labels:", [lab["label"] for lab in fig["labelled"]])
        if "pi0" in fig:
            print("   pi0:", fig["pi0"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
