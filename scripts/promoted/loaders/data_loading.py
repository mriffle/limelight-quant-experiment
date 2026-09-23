"""The ``Dataset`` contract shared by every loader/QC module in this project.

PROJECT COPY seeded from the plugin template ``lib/common/data_loading.py``
(``wide-data-loader`` v0.5). **Trimmed, not verbatim**: this project's four raw files
(protein/peptide FlashLFQ quant matrices + Limelight table dumps) do not fit the
template's generic ``id_columns`` + single sample-column-family shape — they need
NaN/0 missing-token disambiguation, contaminant-regex feature metadata, a per-cell
Detection Type matrix (peptide), and value-correspondence label resolution
(Limelight). Per the template's own guidance ("a study whose files don't fit this
shape should get its own loader ... that returns the same Dataset structure"), the
project loaders (``protein_loader.py``, ``peptide_loader.py``,
``limelight_loader.py``) are custom, and only the template's **contract** —
``Dataset``/``Scale``/``LOG_SCALES`` — is carried over verbatim here, since it is
exactly what lets those custom loaders compose unchanged with the seeded
``normalize.py`` / ``batch_correct.py`` / ``dataset_io.py`` / ``missing_values.py``
modules in this package. ``ReplicateCollapse``/``load_wide_data``/
``load_precursor_data`` are deliberately NOT carried over: nothing in this project
calls them, and the maximal-testing rule (conventions/coding.md) would otherwise
require testing unused code, which is disproportionate (the scientist's
proportionality request in this task).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

# The scale an abundance matrix is on. Lets normalize/batch-correct refuse
# scale-incorrect operations (no double-log; ComBat only on a log-ish scale) instead
# of trusting a convention. See the lib template for the full rationale per member.
Scale = Literal["linear", "log2", "log10", "ln", "glog2", "zscore", "ratio"]

# The subset of scales on which ComBat and log-domain median centering are meaningful.
LOG_SCALES: frozenset[Scale] = frozenset({"log2", "log10", "ln", "glog2", "zscore"})

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": ["Dataset", "Scale", "LOG_SCALES"],
    "uses": [],
    "seeded_from": {"template": "wide-data-loader", "version": "0.5"},
    "description": (
        "The Dataset contract (abundances/feature_names/feature_metadata/metadata/"
        "scale) carried over from the wide-data-loader template so the project's "
        "custom protein/peptide/Limelight loaders and the seeded normalize/"
        "batch_correct/dataset_io/missing_values modules compose unchanged. Trimmed "
        "of the template's generic id_columns loaders (load_wide_data/"
        "load_precursor_data), which no project script calls."
    ),
}


@dataclass
class Dataset:
    """A wide omics dataset aligned to its sample metadata.

    Attributes
    ----------
    abundances:
        ``(n_samples, n_features)`` float array. Row ``i`` corresponds to
        ``metadata.iloc[i]``. On the file's scale (this project: linear intensity;
        missing values are ``NaN`` — see each loader's module docstring for the
        raw-token -> NaN policy).
    feature_names:
        ``(n_features,)`` str array — the chosen feature id column, in the same
        order as the columns of ``abundances``.
    feature_metadata:
        ``(n_features, ...)`` the data file's non-sample id/annotation columns, in
        feature order.
    metadata:
        ``(n_samples, ...)`` the sample-metadata rows, in ``abundances`` row order.
    scale:
        Which scale ``abundances`` is on (see :data:`Scale`). Stamped by the loader
        and updated by each transform (normalization, log2). Downstream code reads
        it to refuse scale-incorrect operations.
    """

    abundances: np.ndarray
    feature_names: np.ndarray
    feature_metadata: pd.DataFrame
    metadata: pd.DataFrame
    scale: Scale = "linear"
