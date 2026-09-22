"""Save figures the way the workflow requires: dual export + a separate legend image.

Project module seeded from the plugin ``lib/figures`` template ``figure-io`` v0.3.
Held to the correctness charter (conventions/correctness.md): **assume nothing,
verify everything, fail loud.**

DEVIATION FROM THE TEMPLATE (beyond the header/``__script_meta__``): :func:`save_figure`
now pins ``svg.hashsalt`` to a fixed value and passes ``metadata={"Date": None}`` to
every SVG ``savefig`` call. Upstream matplotlib's SVG writer embeds a random
per-process hash salt (used to namespace element ids like clip paths) and the
current wall-clock time in a ``<dc:date>`` tag; both make two SVG renders of the
*exact same figure* byte-different, which breaks any downstream "is this
reproducible" diff (e.g. the two-runs-are-byte-identical check this project's
scripts are held to). Neither field carries information the project needs, so
both are pinned/suppressed here. PNG output is unaffected (Agg's PNG writer does
not embed either field).

Encodes three visualization conventions mechanically (conventions/visualization.md), so
every figure starts compliant instead of relying on each script to remember:

  * **Dual export.** Every figure is written as an **SVG** (vector master, for editing /
    publication) *and* a **PNG at 300 DPI** (the raster the figure-reviewer inspects and
    findings embed). :func:`save_figure` writes both from one call.
  * **Legend as a separate image.** A matplotlib legend baked into the plot routinely
    overlaps the data. So the legend is rendered as its **own figure** and dual-exported
    to ``<base>.legend.svg`` / ``<base>.legend.png`` beside the main figure, which stays
    clean. Pass it as ``legend_fig``; the plot template builds the swatch/colorbar key.
  * **Publication-ready defaults.** :func:`publication_style` is a shared matplotlib
    style (legible fonts at print scale, no chartjunk, clean spines). A figure script
    renders *inside* this context so its defaults are set centrally, not re-specified.

The set each plot produces — ``<base>.svg`` + ``<base>.png`` (+ ``<base>.legend.svg`` /
``<base>.legend.png`` when a legend figure is supplied) in the figure's directory under
``figures/`` (``figures/<phase>/<family>[/<label>]/`` — conventions/visualization.md,
*Where figures live*) — is what the figure-reviewer checks for.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": [
        "PUBLICATION_RCPARAMS",
        "FigureArtifacts",
        "publication_style",
        "save_figure",
    ],
    "uses": [],
    "seeded_from": {"template": "figure-io", "version": "0.3"},
    "description": (
        "Figure save helpers enforcing the visualization conventions: dual export "
        "(SVG vector + 300-DPI PNG) plus an optional companion legend figure exported "
        "as <base>.legend.{svg,png} (kept out of the plot so it cannot overlap the "
        "data), and a shared publication matplotlib style. Deviation from "
        "figure-io@0.3: SVG output is byte-deterministic (fixed svg.hashsalt, no Date "
        "metadata), and save_figure must be called inside publication_style(). "
        "Study-agnostic; fail-loud."
    ),
}

# Default raster resolution (DPI) for the PNG export — the convention's 300 DPI.
DEFAULT_DPI = 300

# Fixed (not None/random) so repeated renders of the *same* figure produce a
# byte-identical SVG. matplotlib's default ``svg.hashsalt`` is a fresh random
# value per process, used only to namespace internal element ids (e.g. clip
# paths) — it carries no information, but its randomness alone makes two SVGs
# of the same figure byte-different across runs. See the module docstring.
_SVG_HASHSALT = "findings-workflow-metadata-figures"

# Shared publication style. Sizes are tuned for legibility at print scale; spines are
# trimmed to reduce chartjunk. A script applies it via ``with publication_style():`` so
# the defaults live here, not scattered across figure scripts. Individual scripts may
# still override any rcParam locally for a specific figure.
PUBLICATION_RCPARAMS: dict[str, Any] = {
    "figure.dpi": 100,
    "savefig.dpi": DEFAULT_DPI,
    "savefig.bbox": "tight",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "figure.autolayout": False,
    "svg.fonttype": "none",  # keep text as text in the SVG (editable, smaller files)
}


@dataclass(frozen=True)
class FigureArtifacts:
    """The files written for one figure.

    Attributes
    ----------
    svg:
        Vector master ``<base>.svg``.
    png:
        Raster ``<base>.png`` at the requested DPI (the review/embed target).
    legend_svg, legend_png:
        The companion legend figure ``<base>.legend.{svg,png}``, or ``None`` when no
        legend figure was supplied.
    """

    svg: Path
    png: Path
    legend_svg: Path | None
    legend_png: Path | None


@contextmanager
def publication_style() -> Iterator[None]:
    """Apply :data:`PUBLICATION_RCPARAMS` for the duration of the ``with`` block.

    Uses :func:`matplotlib.rc_context`, so the changes are scoped and the caller's
    global rcParams are restored on exit (even if the body raises). Build AND save the
    figure *inside* the block: the style must be active while the artists are created
    and at ``savefig`` time (``svg.fonttype``, the SVG hash salt and metadata are read
    when saving)::

        with publication_style():
            fig = make_figure(...)
            save_figure(fig, "figures/qc/pca", "pca_by_genotype", legend_fig=legend)

    (Deviation from the figure-io@0.3 template, whose example saved after the block.)
    """
    # rc_context() with no argument snapshots the current rcParams and restores them on
    # exit; we mutate inside it so the override is scoped and exception-safe.
    with mpl.rc_context():
        # matplotlib types rcParams keys as a giant Literal; our style table is a plain
        # str-keyed dict by design (it is data, not code), so update() needs the escape.
        mpl.rcParams.update(PUBLICATION_RCPARAMS)  # type: ignore[arg-type]
        yield


def save_figure(
    fig: Figure,
    output_dir: str | Path,
    base_name: str,
    *,
    legend_fig: Figure | None = None,
    dpi: int = DEFAULT_DPI,
    close: bool = True,
) -> FigureArtifacts:
    """Dual-export ``fig`` (SVG + PNG) and any ``legend_fig`` as ``<base>.legend.*``.

    Parameters
    ----------
    fig:
        The main figure to save.
    output_dir:
        Destination directory (created if absent). In a project this is the figure's
        directory under the structured layout, e.g. ``figures/qc/pca`` or
        ``figures/analysis/differential-abundance/genotype-vs-wt``.
    base_name:
        Filename stem (no extension), e.g. ``"pca_experimental_by_genotype"``. Must be a
        bare name, not a path — a stem with a path separator is rejected so the figure
        and its legend cannot land in different places.
    legend_fig:
        Optional companion legend figure (a standalone swatch/colorbar key). When given
        it is dual-exported to ``<base>.legend.svg`` / ``<base>.legend.png`` so the
        legend never overlaps the plot. ``None`` for figures needing no separate legend.
    dpi:
        PNG resolution (default 300, the convention's floor). The SVG is scale-free.
    close:
        If ``True`` (default), close the figure(s) on **every** exit — success or
        exception — so a batch loop that catches errors never leaks open figures. Pass
        ``False`` to keep inspecting them (e.g. in tests); the caller then owns closing.

    Returns
    -------
    FigureArtifacts
        Paths to ``(svg, png, legend_svg, legend_png)``; the legend paths are ``None``
        when no ``legend_fig`` was supplied.
    """
    # With close=True, save_figure OWNS the figures and closes them on every exit —
    # success, a savefig error, or a bad-arg raise — so callers (e.g. a batch loop that
    # catches errors) never leak open figures. With close=False the caller keeps them.
    # Validation lives INSIDE the try so even a rejected base_name closes the figures.
    try:
        # Reject any path separator (both `/` and `\`, regardless of platform — these
        # seeds may run on Windows) and the `.`/`..` specials, so the figure and its
        # legend can't land in different directories.
        separators = {"/", "\\", os.sep} | ({os.altsep} if os.altsep else set())
        if (
            not base_name
            or base_name in {".", ".."}
            or any(s in base_name for s in separators)
        ):
            raise ValueError(
                f"base_name must be a bare filename stem, not a path; got "
                f"{base_name!r}."
            )
        if dpi <= 0:
            raise ValueError(f"dpi must be positive; got {dpi}.")

        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        svg_path = out / f"{base_name}.svg"
        png_path = out / f"{base_name}.png"
        # Pin the hash salt and suppress the embedded render date for the SVG
        # exports only (see module docstring) so re-rendering the same figure
        # is byte-identical; PNG is unaffected by either setting.
        with mpl.rc_context({"svg.hashsalt": _SVG_HASHSALT}):
            fig.savefig(
                svg_path, format="svg", bbox_inches="tight", metadata={"Date": None}
            )
            fig.savefig(png_path, format="png", dpi=dpi, bbox_inches="tight")

            legend_svg: Path | None = None
            legend_png: Path | None = None
            if legend_fig is not None:
                legend_svg = out / f"{base_name}.legend.svg"
                legend_png = out / f"{base_name}.legend.png"
                legend_fig.savefig(
                    legend_svg,
                    format="svg",
                    bbox_inches="tight",
                    metadata={"Date": None},
                )
                legend_fig.savefig(
                    legend_png, format="png", dpi=dpi, bbox_inches="tight"
                )

        return FigureArtifacts(
            svg=svg_path, png=png_path, legend_svg=legend_svg, legend_png=legend_png
        )
    finally:
        if close:
            plt.close(fig)
            if legend_fig is not None:
                plt.close(legend_fig)
