"""Stage-1 design rules shared by ``metadata_characterize.py`` and
``metadata_figures.py`` — anything both scripts must independently derive from
the same rule has to be defined exactly once here, not re-derived (possibly
inconsistently) in each script.

Currently this is:

* the run-half split (:func:`run_half_label`) — ``metadata_characterize.py``
  uses it to *write* the ``run_half`` column of ``samples.tsv``;
  ``metadata_figures.py`` uses the identical function to *re-derive* the
  expected value and cross-check it against what was written, so the two can
  never silently disagree on the tie-break for an odd-sized batch.
* :data:`FAILURE_MARKER_NAME` — the filename ``metadata_characterize.py`` writes
  into its output directory when a run fails, and that
  ``metadata_figures.py`` refuses to render from if present (stale/known-bad
  upstream tables).
"""

from __future__ import annotations

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": ["run_half_label", "FAILURE_MARKER_NAME"],
    "uses": [],
    "seeded_from": None,
    "description": (
        "The single run-half (ceil(n/2), middle sample -> early) rule and the "
        "upstream-failure-marker filename shared by metadata_characterize.py "
        "and metadata_figures.py."
    ),
}


# Filename metadata_characterize.py writes into its output directory when a run
# fails (cleared at the start of every run); metadata_figures.py refuses to
# render from a metadata_dir where this marker is present.
FAILURE_MARKER_NAME = "FAILED.json"


def run_half_label(position: int, batch_size: int) -> str:
    """ "early" or "late" for a 1-indexed run ``position`` within a batch of
    ``batch_size`` samples.

    The split is ``ceil(batch_size / 2)`` positions in "early", the rest
    "late" — a deliberate, documented tie-break for an ODD-sized batch: the
    extra (middle) sample is assigned to "early", not "late". This is the ONE
    place that rule is written down; every other computation of run_half must
    call this function rather than re-deriving the arithmetic.
    """
    if batch_size <= 0:
        raise ValueError(f"batch_size must be positive; got {batch_size}.")
    if not 1 <= position <= batch_size:
        raise ValueError(
            f"position {position} is out of range for batch_size {batch_size} "
            f"(expected 1..{batch_size})."
        )
    half_size = -(-batch_size // 2)  # ceil(batch_size / 2)
    return "early" if position <= half_size else "late"
