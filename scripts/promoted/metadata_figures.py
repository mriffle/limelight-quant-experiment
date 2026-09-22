"""Stage-1 cohort / metadata figures for the raloxifene/HLM LFQ experiment.

Renders the Stage-1 descriptive figure family from the **precomputed** tables that
``scripts/scratch/metadata_characterize.py`` wrote to ``results/metadata/``; it reads no
raw data and recomputes no statistic:

1. ``figures/metadata/distributions/cohort-counts`` - sample counts per level of
   condition, batch (proxy = acquisition date) and candidate pair (small multiples).
2. ``figures/metadata/crosstabs/0002-condition-by-batch`` - condition x batch counts
   (embedded in caveat finding 0002).
3. ``figures/metadata/crosstabs/0001-condition-by-run-half`` - condition x within-batch
   run half (early / late).
4. ``figures/metadata/run-layout/0001-run-layout`` - one row per batch, x = within-batch
   run position (from the file sequence number, *presumed* acquisition order), one
   marker per run colored by condition and labeled with its sample id (key evidence
   figure of caveat finding 0001, run order aliased with condition).

Every figure is dual-exported (SVG + 300-DPI PNG) with a separate legend image
``<stem>.legend.{svg,png}`` through the project's ``common.figures.figure_io``; every
categorical color is read from ``state/color_registry.json`` through
``common.figures.colors`` (read-only here: a level missing from the registry raises
rather than being assigned a new color). The precomputed crosstabs and run layout are
cross-checked against ``samples.tsv`` before anything is drawn (fail loud on any
disagreement), using the SAME run-half rule (``common.design.run_half_label``) that
``metadata_characterize.py`` used to write it, so the two can never silently disagree.

This script refuses to render at all if ``results/metadata/FAILED.json`` is present —
that marker means the upstream ``metadata_characterize.py`` run failed, so anything
already in ``results/metadata/`` may be stale/inconsistent leftovers from an earlier,
different run.

Per-figure provenance (script path + sha256, git commit if any, data_version, input
table hashes — including the color registry and the figure-code modules that shaped the
render, not just the data tables — parameters, artifact paths) is written to a JSON file
(default ``results/metadata/figure_provenance.json``).

No stochastic step (no RNG, no seed). No processing state applies (metadata only).

Run (from the project root):
    ./.venv/bin/python scripts/scratch/metadata_figures.py
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common.design import FAILURE_MARKER_NAME, run_half_label
from common.figures import colors as colors_module
from common.figures import figure_io as figure_io_module
from common.figures.colors import DEFAULT_REGISTRY_PATH, assign_colors, load_registry
from common.figures.figure_io import FigureArtifacts, publication_style, save_figure
from common.hashing import sha256_of_file

__script_meta__: dict[str, object] = {
    "task": "metadata-figures",
    "kind": "analysis",
    "provides": [],
    "uses": [
        "common.design",
        "common.hashing",
        "common.figures.colors",
        "common.figures.figure_io",
    ],
    # Written from scratch for this project; it only *uses* the seeded figure-io /
    # okabe-ito-colors modules (imports above), it was not itself seeded from either
    # template.
    "seeded_from": None,
    "description": (
        "Stage-1 cohort figures (level counts, condition x batch, condition x run "
        "half, within-batch run layout) rendered from the precomputed "
        "results/metadata tables; registry colors, dual export + separate legend, "
        "per-figure provenance JSON. Refuses to render if the upstream run failed "
        "(FAILED.json present)."
    ),
}

LOGGER = logging.getLogger("metadata_figures")

# --------------------------------------------------------------------------- #
# Constants (study vocabulary; checked against the data, never assumed)
# --------------------------------------------------------------------------- #

CONDITION_ORDER: tuple[str, ...] = ("control", "raloxifene-d0")
# Second (redundant) channel for the run-layout markers so the figure also reads in
# grayscale / for color-vision-deficient viewers.
CONDITION_MARKERS: Mapping[str, str] = {"control": "o", "raloxifene-d0": "s"}
RUN_HALF_ORDER: tuple[str, ...] = ("early", "late")
REQUIRED_SAMPLE_COLUMNS: tuple[str, ...] = (
    "sample_id",
    "condition",
    "batch",
    "seq_number",
    "run_position_within_batch",
    "run_half",
    "candidate_pair",
)
STRATIFIED_SCOPE = "stratified_across_batches"
EDGE_COLOR = "#000000"
BAR_EDGE_WIDTH = 0.8
COUNT_LABEL_SIZE = 10

STEM_COUNTS = "cohort-counts"
STEM_BATCH = "0002-condition-by-batch"
STEM_RUN_HALF = "0001-condition-by-run-half"
STEM_RUN_LAYOUT = "0001-run-layout"

INPUT_FILES: tuple[str, ...] = (
    "samples.tsv",
    "crosstab_condition_batch.tsv",
    "crosstab_condition_run_half.tsv",
    "run_layout.tsv",
    "hypotheses.tsv",
    "data_version.json",
)


class MetadataConsistencyError(ValueError):
    """A precomputed table disagrees with ``samples.tsv`` (or is malformed)."""


class UpstreamRunFailedError(RuntimeError):
    """``metadata_dir`` carries a ``FAILED.json`` marker from a failed upstream run.

    Anything else in ``metadata_dir`` may be stale/inconsistent leftovers from an
    earlier, different (successful) run — rendering from it would silently show
    figures that do not correspond to the latest attempt. Refuse outright rather
    than render from possibly-stale tables.
    """


# --------------------------------------------------------------------------- #
# Loading + validation
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class MetadataTables:
    """The precomputed Stage-1 tables this script renders from.

    Attributes
    ----------
    samples:
        One row per run (``samples.tsv``).
    condition_by_batch, condition_by_run_half:
        Count crosstabs indexed by condition (``crosstab_*.tsv``).
    run_layout:
        Batch x run-position table of ``"<condition>/<sample_id>"`` cells.
    data_version:
        The ``data_version`` string from ``data_version.json``.
    stratified_p_two_sided:
        Two-sided p of the stratified exact permutation test (the direction was not
        pre-specified, so the two-sided value is the one shown) (``hypotheses.tsv``, H4,
        scope ``stratified_across_batches``).
    input_hashes:
        sha256 of every input file read, keyed by file name.
    """

    samples: pd.DataFrame
    condition_by_batch: pd.DataFrame
    condition_by_run_half: pd.DataFrame
    run_layout: pd.DataFrame
    data_version: str
    stratified_p_two_sided: float
    input_hashes: dict[str, str]


def _read_crosstab(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype={"condition": str})
    if "condition" not in table.columns:
        raise MetadataConsistencyError(f"{path} has no 'condition' column.")
    table = table.set_index("condition")
    table.columns = [str(c) for c in table.columns]
    return table.astype(int)


def _read_stratified_p(path: Path) -> float:
    hyp = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    rows = hyp[(hyp["hypothesis_id"] == "H4") & (hyp["scope"] == STRATIFIED_SCOPE)]
    if len(rows) != 1:
        raise MetadataConsistencyError(
            f"{path}: expected exactly one H4 row with scope {STRATIFIED_SCOPE!r}, "
            f"found {len(rows)}."
        )
    raw_p = rows["p_two_sided"].iloc[0]
    if not str(raw_p).strip():
        raise MetadataConsistencyError(
            f"{path}: the stratified H4 row's p_two_sided is blank. This happens "
            "when every batch is single-condition (u_max=0, no within-batch "
            "comparison is possible anywhere) — a degenerate design this figure "
            "cannot meaningfully summarize with a single p-value; it is not "
            "rendered rather than shown with a fabricated number."
        )
    p = float(raw_p)
    if not 0.0 <= p <= 1.0:
        raise MetadataConsistencyError(f"{path}: stratified p = {p} is not in [0, 1].")
    return p


def load_tables(metadata_dir: Path) -> MetadataTables:
    """Read and structurally validate the precomputed tables in ``metadata_dir``."""
    marker = metadata_dir / FAILURE_MARKER_NAME
    if marker.is_file():
        raise UpstreamRunFailedError(
            f"{marker} is present: the upstream metadata_characterize.py run "
            f"that should have populated {metadata_dir} failed. Refusing to "
            f"render figures from a directory that may hold stale/inconsistent "
            f"tables from an earlier run. Fix the upstream failure (see "
            f"{marker} for the recorded error) and re-run it before rendering."
        )
    for name in INPUT_FILES:
        if not (metadata_dir / name).is_file():
            raise FileNotFoundError(f"Required input {metadata_dir / name} is missing.")
    samples = pd.read_csv(metadata_dir / "samples.tsv", sep="\t", dtype=str)
    if samples.empty:
        raise MetadataConsistencyError(
            f"{metadata_dir / 'samples.tsv'} has no rows; there is no cohort to plot."
        )
    missing = [c for c in REQUIRED_SAMPLE_COLUMNS if c not in samples.columns]
    if missing:
        raise MetadataConsistencyError(f"samples.tsv lacks columns {missing}.")
    if samples[list(REQUIRED_SAMPLE_COLUMNS)].isna().any().any():
        raise MetadataConsistencyError("samples.tsv has blank required cells.")
    if samples["sample_id"].duplicated().any():
        raise MetadataConsistencyError("samples.tsv has duplicate sample_id values.")
    samples["seq_number"] = samples["seq_number"].astype(int)
    samples["run_position_within_batch"] = samples["run_position_within_batch"].astype(
        int
    )
    unknown = sorted(set(samples["condition"]) - set(CONDITION_ORDER))
    if unknown:
        raise MetadataConsistencyError(f"Unexpected condition levels {unknown}.")
    bad_half = sorted(set(samples["run_half"]) - set(RUN_HALF_ORDER))
    if bad_half:
        raise MetadataConsistencyError(f"Unexpected run_half levels {bad_half}.")

    version_raw: Any = json.loads(
        (metadata_dir / "data_version.json").read_text(encoding="utf-8")
    )
    if not isinstance(version_raw, dict) or not isinstance(
        version_raw.get("data_version"), str
    ):
        raise MetadataConsistencyError(
            "data_version.json must be an object with a string 'data_version'."
        )

    run_layout = pd.read_csv(metadata_dir / "run_layout.tsv", sep="\t", dtype=str)
    if "batch" not in run_layout.columns:
        raise MetadataConsistencyError("run_layout.tsv has no 'batch' column.")
    run_layout = run_layout.set_index("batch")

    return MetadataTables(
        samples=samples,
        condition_by_batch=_read_crosstab(
            metadata_dir / "crosstab_condition_batch.tsv"
        ),
        condition_by_run_half=_read_crosstab(
            metadata_dir / "crosstab_condition_run_half.tsv"
        ),
        run_layout=run_layout,
        data_version=str(version_raw["data_version"]),
        stratified_p_two_sided=_read_stratified_p(metadata_dir / "hypotheses.tsv"),
        input_hashes={n: sha256_of_file(metadata_dir / n) for n in INPUT_FILES},
    )


def _crosstab_from_samples(samples: pd.DataFrame, column: str) -> pd.DataFrame:
    counts = pd.crosstab(samples["condition"], samples[column]).astype(int)
    counts.columns = [str(c) for c in counts.columns]
    return counts


def _assert_same_counts(expected: pd.DataFrame, found: pd.DataFrame, name: str) -> None:
    """Compare two count tables, treating absent rows/columns as zero counts."""
    rows = sorted(set(expected.index) | set(found.index))
    cols = sorted(set(expected.columns) | set(found.columns))
    exp = expected.reindex(index=rows, columns=cols, fill_value=0)
    got = found.reindex(index=rows, columns=cols, fill_value=0)
    if not exp.equals(got):
        raise MetadataConsistencyError(
            f"{name} disagrees with samples.tsv:\nfrom samples.tsv:\n{exp}\n"
            f"precomputed:\n{got}"
        )


def verify_consistency(tables: MetadataTables) -> None:
    """Fail loud unless every precomputed table agrees with ``samples.tsv``."""
    s = tables.samples
    _assert_same_counts(
        _crosstab_from_samples(s, "batch"),
        tables.condition_by_batch,
        "crosstab_condition_batch.tsv",
    )
    _assert_same_counts(
        _crosstab_from_samples(s, "run_half"),
        tables.condition_by_run_half,
        "crosstab_condition_run_half.tsv",
    )
    # Run position must be the 1..n rank of seq_number within each batch, and the
    # run_half the first / second half of those positions.
    for batch, grp in s.groupby("batch"):
        ordered = grp.sort_values("seq_number")
        positions = ordered["run_position_within_batch"].tolist()
        if positions != list(range(1, len(grp) + 1)):
            raise MetadataConsistencyError(
                f"Batch {batch}: run positions {positions} are not the seq_number "
                f"rank 1..{len(grp)}."
            )
        halves = [
            run_half_label(int(pos), len(grp))
            for pos in ordered["run_position_within_batch"]
        ]
        if halves != ordered["run_half"].tolist():
            raise MetadataConsistencyError(
                f"Batch {batch}: run_half {ordered['run_half'].tolist()} does not "
                f"match the first/second half of run positions ({halves})."
            )
    # run_layout.tsv cells must be exactly "<condition>/<sample_id>" at each position.
    expected_cells: dict[tuple[str, str], str] = {
        (str(r.batch), str(r.run_position_within_batch)): f"{r.condition}/{r.sample_id}"
        for r in s.itertuples(index=False)
    }
    found_cells: dict[tuple[str, str], str] = {}
    for batch, row in tables.run_layout.iterrows():
        for pos, cell in row.items():
            if isinstance(cell, str) and cell.strip():
                found_cells[(str(batch), str(pos))] = cell.strip()
    if expected_cells != found_cells:
        raise MetadataConsistencyError(
            f"run_layout.tsv disagrees with samples.tsv: expected {expected_cells}, "
            f"found {found_cells}."
        )


# --------------------------------------------------------------------------- #
# Colors (registry, read-only)
# --------------------------------------------------------------------------- #


def registry_colors(
    category: str, levels: Sequence[str], registry_path: Path
) -> dict[str, str]:
    """Return ``{level: hex}`` for ``levels`` from the registry, never adding colors.

    Raises if the category or any level is absent: colors for this project's categories
    are curated in the registry, so a missing level is an error, not a new assignment.
    """
    registry = load_registry(registry_path)
    entry = registry.get(category)
    values = entry.get("values") if isinstance(entry, dict) else None
    if not isinstance(values, dict):
        raise KeyError(f"Color registry {registry_path} has no category {category!r}.")
    absent = [lv for lv in levels if lv not in values]
    if absent:
        raise KeyError(
            f"Color registry {registry_path} category {category!r} lacks levels "
            f"{absent}; add them to the registry (no ad-hoc colors)."
        )
    # Route through assign_colors so the >8-category guard applies; persist=False and
    # every level already present means nothing is (or could be) written.
    return assign_colors(
        category, list(levels), registry_path=registry_path, persist=False
    )


# --------------------------------------------------------------------------- #
# Legends (separate images)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class LegendSection:
    """One titled block of a legend image: ``(label, handle)`` pairs."""

    title: str
    entries: tuple[tuple[str, Artist], ...]


def bar_handle(color: str) -> Patch:
    """A swatch matching the bars drawn by this script."""
    return Patch(facecolor=color, edgecolor=EDGE_COLOR, linewidth=BAR_EDGE_WIDTH)


def marker_handle(color: str, marker: str) -> Line2D:
    """A marker swatch matching the run-layout markers."""
    return Line2D(
        [],
        [],
        linestyle="none",
        marker=marker,
        markersize=10,
        markerfacecolor=color,
        markeredgecolor=EDGE_COLOR,
        markeredgewidth=BAR_EDGE_WIDTH,
    )


def build_legend_figure(sections: Sequence[LegendSection]) -> Figure:
    """Render a standalone legend: one column per section, bold section titles."""
    if not sections:
        raise ValueError("A legend needs at least one section.")
    height = max(len(sec.entries) for sec in sections) + 1
    blank = Patch(visible=False)
    handles: list[Artist] = []
    labels: list[str] = []
    title_rows: list[int] = []
    for sec in sections:
        title_rows.append(len(labels))
        handles.append(blank)
        labels.append(sec.title)
        for label, handle in sec.entries:
            handles.append(handle)
            labels.append(label)
        pad = height - 1 - len(sec.entries)
        handles.extend([blank] * pad)
        labels.extend([""] * pad)
    fig = plt.figure(figsize=(0.1, 0.1))
    legend = fig.legend(
        handles,
        labels,
        loc="center",
        ncol=len(sections),
        frameon=False,
        handlelength=1.4,
        columnspacing=2.0,
    )
    texts = legend.get_texts()
    for idx in title_rows:
        texts[idx].set_fontweight("bold")
    return fig


# --------------------------------------------------------------------------- #
# Plot helpers
# --------------------------------------------------------------------------- #


def _annotate_bar(ax: Axes, x: float, count: int) -> None:
    ax.annotate(
        str(count),
        (x, count),
        xytext=(0, 2),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=COUNT_LABEL_SIZE,
    )


def _integer_count_axis(ax: Axes, max_count: int) -> None:
    ax.set_ylim(0, max_count + 1)
    ax.set_yticks(list(range(0, max_count + 1, 2 if max_count > 6 else 1)))


def _level_counts(
    samples: pd.DataFrame, column: str, levels: Sequence[str]
) -> list[int]:
    counts = samples[column].value_counts()
    return [int(counts.get(lv, 0)) for lv in levels]


# --------------------------------------------------------------------------- #
# Figure 1 - level counts (small multiples)
# --------------------------------------------------------------------------- #


def plot_cohort_counts(
    samples: pd.DataFrame,
    panels: Sequence[tuple[str, str, Sequence[str], Mapping[str, str]]],
) -> tuple[Figure, Figure]:
    """Bar counts per level, one panel per variable.

    ``panels`` holds ``(column, panel_title, levels, colors)`` per variable.
    """
    n = len(samples)
    widths = [max(len(levels), 2) for _, _, levels, _ in panels]
    fig, axes = plt.subplots(
        1,
        len(panels),
        figsize=(2.0 + 1.05 * sum(widths), 3.8),
        sharey=True,
        gridspec_kw={"width_ratios": widths},
    )
    axes_list: list[Axes] = list(axes) if len(panels) > 1 else [axes]
    max_count = 0
    sections: list[LegendSection] = []
    for ax, (column, title, levels, colors) in zip(axes_list, panels, strict=True):
        counts = _level_counts(samples, column, levels)
        if sum(counts) != n:
            raise MetadataConsistencyError(
                f"{column}: level counts sum to {sum(counts)}, expected n = {n}."
            )
        max_count = max(max_count, *counts)
        xs = list(range(len(levels)))
        ax.bar(
            xs,
            counts,
            width=0.7,
            color=[colors[lv] for lv in levels],
            edgecolor=EDGE_COLOR,
            linewidth=BAR_EDGE_WIDTH,
        )
        for x, c in zip(xs, counts, strict=True):
            _annotate_bar(ax, x, c)
        ax.set_xticks(xs)
        ax.set_xticklabels(
            list(levels), rotation=30, ha="right", rotation_mode="anchor"
        )
        ax.set_xlim(-0.6, max(len(levels), 2) - 0.4)
        ax.set_title(title, fontsize=12)
        sections.append(
            LegendSection(title, tuple((lv, bar_handle(colors[lv])) for lv in levels))
        )
    _integer_count_axis(axes_list[0], max_count)
    axes_list[0].set_ylabel("Samples (n)")
    fig.suptitle(f"Cohort composition (n = {n})", fontweight="bold", fontsize=14)
    fig.tight_layout()
    return fig, build_legend_figure(sections)


# --------------------------------------------------------------------------- #
# Figures 2 + 3 - condition x <variable> grouped bars
# --------------------------------------------------------------------------- #


def plot_condition_crosstab(
    crosstab: pd.DataFrame,
    column_levels: Sequence[str],
    condition_colors: Mapping[str, str],
    *,
    title: str,
    xlabel: str,
) -> tuple[Figure, Figure]:
    """Grouped bars: x = the crosstab's column levels, one bar per condition."""
    table = crosstab.reindex(
        index=list(CONDITION_ORDER), columns=list(column_levels), fill_value=0
    )
    n = int(table.to_numpy().sum())
    width = 0.36
    offsets = [
        (i - (len(CONDITION_ORDER) - 1) / 2) * width
        for i in range(len(CONDITION_ORDER))
    ]
    fig, ax = plt.subplots(figsize=(1.6 + 1.6 * len(column_levels), 3.8))
    for cond, off in zip(CONDITION_ORDER, offsets, strict=True):
        counts = [int(table.loc[cond, lv]) for lv in column_levels]
        xs = [i + off for i in range(len(column_levels))]
        ax.bar(
            xs,
            counts,
            width=width * 0.95,
            color=condition_colors[cond],
            edgecolor=EDGE_COLOR,
            linewidth=BAR_EDGE_WIDTH,
        )
        for x, c in zip(xs, counts, strict=True):
            _annotate_bar(ax, x, c)
    col_totals = [int(table[lv].sum()) for lv in column_levels]
    ax.set_xticks(list(range(len(column_levels))))
    ax.set_xticklabels(
        [f"{lv}\n(n = {t})" for lv, t in zip(column_levels, col_totals, strict=True)]
    )
    ax.set_xlim(-0.6, len(column_levels) - 0.4)
    _integer_count_axis(ax, int(table.to_numpy().max()))
    ax.set_ylabel("Samples (n)")
    ax.set_xlabel(xlabel)
    ax.set_title(f"{title} (n = {n})")
    fig.tight_layout()
    legend = build_legend_figure(
        [
            LegendSection(
                "Condition",
                tuple((c, bar_handle(condition_colors[c])) for c in CONDITION_ORDER),
            )
        ]
    )
    return fig, legend


# --------------------------------------------------------------------------- #
# Figure 4 - run layout
# --------------------------------------------------------------------------- #


def plot_run_layout(
    samples: pd.DataFrame,
    batch_levels: Sequence[str],
    condition_colors: Mapping[str, str],
    p_two_sided: float,
) -> tuple[Figure, Figure]:
    """One row per batch; x = within-batch run position; markers by condition."""
    n = len(samples)
    max_pos = int(samples["run_position_within_batch"].max())
    fig, ax = plt.subplots(
        figsize=(1.8 + 1.15 * max_pos, 1.2 + 1.1 * len(batch_levels))
    )
    y_of = {b: len(batch_levels) - 1 - i for i, b in enumerate(batch_levels)}
    for batch in batch_levels:
        grp = samples[samples["batch"] == batch]
        y = y_of[batch]
        # Guide line spans only this batch's runs (no implied empty run slots).
        positions = grp["run_position_within_batch"]
        ax.plot(
            [int(positions.min()) - 0.35, int(positions.max()) + 0.35],
            [y, y],
            color="#bbbbbb",
            linewidth=0.8,
            zorder=0,
        )
        for r in grp.itertuples(index=False):
            x = int(r.run_position_within_batch)
            cond = str(r.condition)
            ax.plot(
                x,
                y,
                linestyle="none",
                marker=CONDITION_MARKERS[cond],
                markersize=15,
                markerfacecolor=condition_colors[cond],
                markeredgecolor=EDGE_COLOR,
                markeredgewidth=BAR_EDGE_WIDTH,
                zorder=3,
            )
            ax.annotate(
                str(r.sample_id),
                (x, y),
                xytext=(0, 12),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=10,
            )
            ax.annotate(
                f"seq {int(r.seq_number)}",
                (x, y),
                xytext=(0, -12),
                textcoords="offset points",
                ha="center",
                va="top",
                fontsize=8,
                color="#555555",
            )
    ax.set_yticks([y_of[b] for b in batch_levels])
    ax.set_yticklabels(
        [f"{b}\n(n = {int((samples['batch'] == b).sum())})" for b in batch_levels]
    )
    ax.set_ylim(-0.6, len(batch_levels) - 0.4)
    ax.set_xticks(list(range(1, max_pos + 1)))
    ax.set_xlim(0.4, max_pos + 0.6)
    ax.set_xlabel(
        "Within-batch run position\n"
        "(rank of file sequence number; presumed acquisition order)"
    )
    ax.set_ylabel("Batch (acquisition date)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title(f"Run order by condition within batch (n = {n})", pad=30)
    ax.text(
        0.5,
        1.03,
        f"Stratified-by-batch exact permutation p = {p_two_sided:.3g} (two-sided)",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=11,
    )
    fig.tight_layout()
    legend = build_legend_figure(
        [
            LegendSection(
                "Condition",
                tuple(
                    (c, marker_handle(condition_colors[c], CONDITION_MARKERS[c]))
                    for c in CONDITION_ORDER
                ),
            )
        ]
    )
    return fig, legend


# --------------------------------------------------------------------------- #
# Provenance
# --------------------------------------------------------------------------- #


def _git_commit(project_root: Path) -> str | None:
    """HEAD commit of the project repo, or ``None`` when it is not a git repository."""
    try:
        out = subprocess.run(
            ["git", "-C", str(project_root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def _rel(path: Path, root: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(root.resolve()))
    except ValueError:
        return str(resolved)


def provenance_entry(
    artifacts: FigureArtifacts,
    *,
    description: str,
    tables: MetadataTables,
    inputs_used: Sequence[str],
    code_hashes: Mapping[str, str],
    params: Mapping[str, object],
    project_root: Path,
) -> dict[str, object]:
    """Per-figure provenance record (JSON-serializable).

    ``inputs`` always includes ``code_hashes`` (the color registry + the two
    figure-code modules that shaped every render — colors.py picks the hex
    values, figure_io.py controls the actual bytes written) IN ADDITION TO the
    ``inputs_used`` data tables for this specific figure, since the registry and
    code apply to every figure, not just some.
    """
    script = Path(__file__).resolve()
    return {
        "description": description,
        "svg": _rel(artifacts.svg, project_root),
        "png": _rel(artifacts.png, project_root),
        "legend_svg": _rel(artifacts.legend_svg, project_root)
        if artifacts.legend_svg
        else None,
        "legend_png": _rel(artifacts.legend_png, project_root)
        if artifacts.legend_png
        else None,
        "script": {
            "path": _rel(script, project_root),
            "sha256": sha256_of_file(script),
            # HEAD commit of the project repo AT RENDER TIME (i.e. when this
            # script ran) — not necessarily the commit that produced the
            # results/metadata tables being rendered (that is a separate,
            # earlier script run); None if the project is not a git repository.
            "git_commit": _git_commit(project_root),
            "seeded_from": __script_meta__["seeded_from"],
        },
        "data_version": tables.data_version,
        "inputs": {
            **{name: tables.input_hashes[name] for name in inputs_used},
            **dict(code_hashes),
        },
        "params": dict(params),
        "processing_state": None,
    }


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class OutputDirs:
    """Target directories under the structured ``figures/`` layout."""

    distributions: Path
    crosstabs: Path
    run_layout: Path


def render_all(
    metadata_dir: Path,
    registry_path: Path,
    out: OutputDirs,
    provenance_file: Path,
    *,
    dpi: int = 300,
    project_root: Path | None = None,
) -> dict[str, dict[str, object]]:
    """Render the four figures, write provenance, and return the provenance map."""
    root = project_root if project_root is not None else Path.cwd()
    tables = load_tables(metadata_dir)
    verify_consistency(tables)
    s = tables.samples
    LOGGER.info("Loaded %d samples; data_version=%s", len(s), tables.data_version)

    batch_levels = sorted(s["batch"].unique().tolist())
    pair_levels = sorted(s["candidate_pair"].unique().tolist())
    cond_colors = registry_colors("condition", CONDITION_ORDER, registry_path)
    batch_colors = registry_colors("batch", batch_levels, registry_path)
    pair_colors = registry_colors("candidate_pair", pair_levels, registry_path)
    common_params: dict[str, object] = {
        "metadata_dir": _rel(metadata_dir, root),
        "registry_path": _rel(registry_path, root),
        "dpi": dpi,
    }
    # Code/config that shaped EVERY render (not figure-specific data): the color
    # registry (which hex each level got) and the two figure-code modules
    # (colors.py picks the values; figure_io.py controls the bytes actually
    # written). Located via the imported modules' own __file__ so this is
    # correct regardless of whether this script is running from scripts/scratch
    # or (post-promotion) scripts/promoted.
    code_hashes: dict[str, str] = {
        _rel(registry_path, root): sha256_of_file(registry_path),
        _rel(Path(colors_module.__file__).resolve(), root): sha256_of_file(
            Path(colors_module.__file__).resolve()
        ),
        _rel(Path(figure_io_module.__file__).resolve(), root): sha256_of_file(
            Path(figure_io_module.__file__).resolve()
        ),
    }
    records: dict[str, dict[str, object]] = {}

    # save_figure() runs INSIDE publication_style()'s `with` block (not after
    # it): several PUBLICATION_RCPARAMS entries (svg.fonttype, savefig.dpi,
    # savefig.bbox) only take effect at the moment matplotlib actually draws
    # and serializes the figure (i.e. at savefig time), not at figure/Artist
    # construction time — calling save_figure() after the block exited would
    # silently apply matplotlib's global defaults instead.
    with publication_style():
        fig, leg = plot_cohort_counts(
            s,
            [
                ("condition", "Condition", CONDITION_ORDER, cond_colors),
                ("batch", "Batch (acquisition date)", batch_levels, batch_colors),
                ("candidate_pair", "Candidate pair", pair_levels, pair_colors),
            ],
        )
        art = save_figure(fig, out.distributions, STEM_COUNTS, legend_fig=leg, dpi=dpi)
    records[STEM_COUNTS] = provenance_entry(
        art,
        description="Sample counts per level: condition, batch, candidate_pair.",
        tables=tables,
        inputs_used=["samples.tsv", "data_version.json"],
        code_hashes=code_hashes,
        params={**common_params, "output_dir": _rel(out.distributions, root)},
        project_root=root,
    )

    with publication_style():
        fig, leg = plot_condition_crosstab(
            tables.condition_by_batch,
            batch_levels,
            cond_colors,
            title="Condition \u00d7 batch",
            xlabel="Batch (acquisition date)",
        )
        art = save_figure(fig, out.crosstabs, STEM_BATCH, legend_fig=leg, dpi=dpi)
    records[STEM_BATCH] = provenance_entry(
        art,
        description="Condition x batch sample counts (grouped bars).",
        tables=tables,
        inputs_used=[
            "crosstab_condition_batch.tsv",
            "samples.tsv",
            "data_version.json",
        ],
        code_hashes=code_hashes,
        params={**common_params, "output_dir": _rel(out.crosstabs, root)},
        project_root=root,
    )

    with publication_style():
        fig, leg = plot_condition_crosstab(
            tables.condition_by_run_half,
            RUN_HALF_ORDER,
            cond_colors,
            title="Condition \u00d7 within-batch run half",
            xlabel=(
                "Within-batch run half (by file sequence number, presumed "
                "acquisition order)"
            ),
        )
        art = save_figure(fig, out.crosstabs, STEM_RUN_HALF, legend_fig=leg, dpi=dpi)
    records[STEM_RUN_HALF] = provenance_entry(
        art,
        description="Condition x within-batch run half (early/late) sample counts.",
        tables=tables,
        inputs_used=[
            "crosstab_condition_run_half.tsv",
            "samples.tsv",
            "data_version.json",
        ],
        code_hashes=code_hashes,
        params={**common_params, "output_dir": _rel(out.crosstabs, root)},
        project_root=root,
    )

    with publication_style():
        fig, leg = plot_run_layout(
            s, batch_levels, cond_colors, tables.stratified_p_two_sided
        )
        art = save_figure(fig, out.run_layout, STEM_RUN_LAYOUT, legend_fig=leg, dpi=dpi)
    records[STEM_RUN_LAYOUT] = provenance_entry(
        art,
        description=(
            "Within-batch run layout: one row per batch, runs at their within-batch "
            "position, colored/shaped by condition, labeled with sample id."
        ),
        tables=tables,
        inputs_used=[
            "samples.tsv",
            "run_layout.tsv",
            "hypotheses.tsv",
            "data_version.json",
        ],
        code_hashes=code_hashes,
        params={
            **common_params,
            "output_dir": _rel(out.run_layout, root),
            "p_two_sided_source": (
                f"hypotheses.tsv H4 scope={STRATIFIED_SCOPE} column p_two_sided"
            ),
        },
        project_root=root,
    )

    provenance_file.parent.mkdir(parents=True, exist_ok=True)
    provenance_file.write_text(
        json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    LOGGER.info("Wrote provenance for %d figures to %s", len(records), provenance_file)
    return records


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--metadata-dir", type=Path, default=Path("results/metadata"))
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY_PATH)
    parser.add_argument(
        "--distributions-dir",
        type=Path,
        default=Path("figures/metadata/distributions"),
    )
    parser.add_argument(
        "--crosstabs-dir", type=Path, default=Path("figures/metadata/crosstabs")
    )
    parser.add_argument(
        "--run-layout-dir", type=Path, default=Path("figures/metadata/run-layout")
    )
    parser.add_argument(
        "--provenance-file",
        type=Path,
        default=Path("results/metadata/figure_provenance.json"),
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--log-level", type=str.upper, default="INFO")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=args.log_level, format="%(levelname)s %(name)s: %(message)s"
    )
    LOGGER.info("Parameters: %s", vars(args))
    records = render_all(
        args.metadata_dir,
        args.registry,
        OutputDirs(args.distributions_dir, args.crosstabs_dir, args.run_layout_dir),
        args.provenance_file,
        dpi=args.dpi,
    )
    for stem, rec in records.items():
        LOGGER.info("%s -> %s", stem, rec["png"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
