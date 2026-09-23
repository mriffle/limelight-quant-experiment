"""Verified missing-value handling for a linear-scale ``Dataset``.

PROJECT COPY seeded from the plugin template ``lib/common/missing_values.py``
(``missing-values`` v0.1), carried over verbatim except the package-layout
adaptation: ``from loaders.data_loading import Dataset`` instead of
``from common.data_loading import Dataset``. No other logic differs from the
template. This project uses only the feature-filtering half
(``max_missing_fraction=0.0``, ``impute=None``) to build the "complete
non-contaminant features" QC state (`qc_prep.py`) — the scientist-confirmed
preprocessing decision is complete-case only, no imputation (`state/
DATA_DESCRIPTION.md`), so the ``impute`` methods are carried over unused for
now (kept, not deleted, because a left-censored-imputation sensitivity run is a
named candidate Stage-4 follow-up) and are exercised by tests regardless.

WHAT IT DOES. Two optional, ordered steps that turn a matrix with not-detected
values into a resolved one the downstream templates can consume:

1. **Feature filtering** (``max_missing_fraction``) — drop features detected in too
   few samples (missing in *more than* the given fraction).
2. **Imputation** (``impute``) — fill the remaining not-detected entries by one
   method: ``"zero"``, ``"mean"``/``"median"``, or ``"knn"``.

Either step is optional: ``handle_missing(ds)`` with both ``None`` validates and
passes the matrix through unchanged.

THE DETECTION PREDICATE. A value is **not detected** (missing) when it is
non-finite **or** ``<= min_intensity`` (default ``0.0``). This project's loaders
already convert both the ``0`` and literal ``"NaN"`` raw tokens to ``NaN``
(see ``protein_loader.py`` / ``peptide_loader.py``), so "non-finite" alone covers
this project's missingness; the ``<= min_intensity`` half of the predicate is inert
here (no positive-but-below-threshold values are ever produced) and is kept only
because it is part of the verified template contract.

SCALE IS A HARD REFUSE: the loader's raw-linear Dataset (``scale="linear"``) is the
only valid input — filtering/zero-fill only mean what they should on the raw
linear matrix.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace
from typing import Literal

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

from loaders.data_loading import Dataset

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": [
        "ImputeMethod",
        "MissingValueScaleError",
        "MissingValueResult",
        "handle_missing",
    ],
    "uses": ["loaders.data_loading"],
    "seeded_from": {"template": "missing-values", "version": "0.1"},
    "description": (
        "Verified Dataset missing-value handling: optional feature filtering by "
        "maximum per-feature missing fraction, then optional imputation (zero / "
        "mean / median / KNN) of the remaining not-detected entries. Fail-loud, "
        "returns an independent linear Dataset plus a drop/impute report. Unchanged "
        "from the lib template except the loaders.data_loading import path."
    ),
}

ImputeMethod = Literal["zero", "mean", "median", "knn"]

_PER_FEATURE_STAT: dict[ImputeMethod, Callable[[np.ndarray], np.ndarray]] = {
    "mean": lambda masked: np.asarray(np.nanmean(masked, axis=0), dtype=float),
    "median": lambda masked: np.asarray(np.nanmedian(masked, axis=0), dtype=float),
}


class MissingValueScaleError(ValueError):
    """Raised when missing-value handling is requested on a non-linear scale."""


@dataclass(frozen=True)
class MissingValueResult:
    """A resolved :class:`Dataset` plus a record of what handling was applied.

    Attributes
    ----------
    dataset:
        The new linear-scale :class:`Dataset` after filtering + imputation.
        Independent of the input. Finite everywhere when ``impute`` was set; may
        still carry not-detected values (NaN) when ``impute`` is ``None``.
    n_features_in, n_features_out:
        Feature counts before and after the ``max_missing_fraction`` filter.
    dropped_features:
        ``(n_features_in - n_features_out,)`` str array of the dropped feature names,
        in input order.
    n_values_imputed:
        Number of not-detected entries filled among the kept features (``0`` when
        ``impute`` is ``None``).
    max_missing_fraction, impute, min_intensity:
        The parameters used, echoed for the provenance record.
    """

    dataset: Dataset
    n_features_in: int
    n_features_out: int
    dropped_features: np.ndarray
    n_values_imputed: int
    max_missing_fraction: float | None
    impute: ImputeMethod | None
    min_intensity: float


def handle_missing(
    dataset: Dataset,
    *,
    max_missing_fraction: float | None = None,
    impute: ImputeMethod | None = None,
    min_intensity: float = 0.0,
    knn_neighbors: int = 5,
) -> MissingValueResult:
    """Filter high-missing features and/or impute the rest; return a resolved Dataset.

    The two steps run in order — **filter, then impute**.

    Parameters
    ----------
    dataset:
        Input dataset on the ``"linear"`` scale (raised otherwise).
    max_missing_fraction:
        Drop any feature whose missing fraction is **greater than** this. ``None``
        (default) keeps every feature. Must be in ``[0, 1]``.
    impute:
        How to fill the not-detected entries of the kept features. ``None`` (default)
        leaves them untouched.
    min_intensity:
        Detection threshold; an entry counts as detected when finite and
        ``> min_intensity`` (default ``0.0``).
    knn_neighbors:
        Neighbour count for ``impute="knn"``; default ``5``. Must be ``>= 1``.

    Returns
    -------
    MissingValueResult

    Raises
    ------
    MissingValueScaleError
        If ``dataset`` is not on the ``"linear"`` scale.
    ValueError
        On a non-2D / empty / non-row-aligned Dataset, an out-of-range parameter, a
        filter that removes every feature, or a kept fully-missing feature under a
        mean/median/KNN impute.
    """
    abundances = _validate(dataset, max_missing_fraction, impute, knn_neighbors)
    missing = _missing_mask(abundances, min_intensity)

    keep = _keep_mask(missing, max_missing_fraction)
    if not bool(keep.any()):
        raise ValueError(
            f"max_missing_fraction={max_missing_fraction!r} drops every feature "
            f"({abundances.shape[1]} total); loosen the threshold or check the data."
        )
    feature_names = np.asarray(dataset.feature_names, dtype=str)
    dropped_features = feature_names[~keep].copy()

    kept_abundances = abundances[:, keep]
    kept_missing = missing[:, keep]
    filled, n_imputed = _impute(kept_abundances, kept_missing, impute, knn_neighbors)

    new_dataset = _independent(
        dataset,
        filled,
        feature_names[keep].copy(),
        dataset.feature_metadata.loc[keep].reset_index(drop=True),
    )
    return MissingValueResult(
        dataset=new_dataset,
        n_features_in=int(abundances.shape[1]),
        n_features_out=int(keep.sum()),
        dropped_features=dropped_features,
        n_values_imputed=n_imputed,
        max_missing_fraction=max_missing_fraction,
        impute=impute,
        min_intensity=min_intensity,
    )


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #


def _validate(
    dataset: Dataset,
    max_missing_fraction: float | None,
    impute: ImputeMethod | None,
    knn_neighbors: int,
) -> np.ndarray:
    """Refuse a bad scale / shape / parameter; return the float abundances."""
    if dataset.scale != "linear":
        raise MissingValueScaleError(
            f"Missing-value handling requires linear-scale abundances but the "
            f"Dataset is on scale {dataset.scale!r}. On log/centred scales a 0 is an "
            f"ordinary value, so 'missing' is ill-defined and 'zero' fill is wrong. "
            f"Run it on the raw data, before normalization."
        )
    abundances = np.asarray(dataset.abundances, dtype=float)
    if abundances.ndim != 2:
        raise ValueError(
            f"abundances must be 2D (n_samples, n_features); got {abundances.shape}."
        )
    n_samples, n_features = abundances.shape
    if n_samples < 1 or n_features < 1:
        raise ValueError(
            f"Dataset is empty ({n_samples} samples x {n_features} features); "
            f"need >= 1 of each."
        )
    if len(dataset.feature_names) != n_features:
        raise ValueError(
            f"feature_names has {len(dataset.feature_names)} entries but there are "
            f"{n_features} features; the Dataset is not feature-aligned."
        )
    if len(dataset.feature_metadata) != n_features:
        raise ValueError(
            f"feature_metadata has {len(dataset.feature_metadata)} rows but there "
            f"are {n_features} features; the Dataset is not feature-aligned."
        )
    if max_missing_fraction is not None and not 0.0 <= max_missing_fraction <= 1.0:
        raise ValueError(
            f"max_missing_fraction must be in [0, 1]; got {max_missing_fraction}."
        )
    if impute is not None and impute not in ("zero", "mean", "median", "knn"):
        raise ValueError(
            f"Unknown impute method {impute!r}; use one of "
            f"'zero'/'mean'/'median'/'knn' (or None)."
        )
    if impute == "knn" and knn_neighbors < 1:
        raise ValueError(f"knn_neighbors must be >= 1; got {knn_neighbors}.")
    return abundances


def _missing_mask(abundances: np.ndarray, min_intensity: float) -> np.ndarray:
    """Boolean (n_samples, n_features): True where NOT detected."""
    detected = np.isfinite(abundances) & (abundances > min_intensity)
    return np.asarray(~detected, dtype=bool)


def _keep_mask(missing: np.ndarray, max_missing_fraction: float | None) -> np.ndarray:
    """Per-feature keep mask: keep where missing fraction <= max_missing_fraction."""
    if max_missing_fraction is None:
        return np.ones(missing.shape[1], dtype=bool)
    fraction = np.asarray(missing.mean(axis=0), dtype=float)
    return np.asarray(fraction <= max_missing_fraction, dtype=bool)


# --------------------------------------------------------------------------- #
# Imputation
# --------------------------------------------------------------------------- #


def _impute(
    abundances: np.ndarray,
    missing: np.ndarray,
    impute: ImputeMethod | None,
    knn_neighbors: int,
) -> tuple[np.ndarray, int]:
    """Fill the missing entries of ``abundances`` by ``impute``; return (filled, n)."""
    if impute is None:
        return abundances.copy(), 0

    n_to_impute = int(missing.sum())
    if impute == "zero":
        filled = np.where(missing, 0.0, abundances)
    elif impute == "knn":
        filled = _impute_knn(abundances, missing, knn_neighbors)
    else:
        filled = _impute_per_feature(abundances, missing, impute)

    if not np.isfinite(filled).all():
        n_bad = int((~np.isfinite(filled)).sum())
        raise ValueError(
            f"impute={impute!r} left {n_bad} non-finite value(s); this is a backstop "
            f"failure — check for fully-missing features or non-finite detected values."
        )
    if filled.shape != abundances.shape:
        raise ValueError(
            f"impute={impute!r} changed the matrix shape {abundances.shape} -> "
            f"{filled.shape}; expected it preserved."
        )
    return filled, n_to_impute


def _impute_per_feature(
    abundances: np.ndarray, missing: np.ndarray, impute: ImputeMethod
) -> np.ndarray:
    """Fill each feature's missing entries with its detected mean/median."""
    _require_no_fully_missing(missing, impute)
    detected = ~missing
    masked = np.where(detected, abundances, np.nan)
    stats = _PER_FEATURE_STAT[impute](masked)
    return np.where(missing, stats[np.newaxis, :], abundances)


def _impute_knn(
    abundances: np.ndarray, missing: np.ndarray, knn_neighbors: int
) -> np.ndarray:
    """Fill missing entries via sklearn ``KNNImputer`` (neighbours over samples)."""
    _require_no_fully_missing(missing, "knn")
    to_impute = np.where(missing, np.nan, abundances)
    imputer = KNNImputer(n_neighbors=knn_neighbors)
    # sklearn is untyped, so fit_transform is Any; np.asarray pins it to an ndarray.
    return np.asarray(imputer.fit_transform(to_impute), dtype=float)


def _require_no_fully_missing(missing: np.ndarray, impute: ImputeMethod) -> None:
    """Refuse if any kept feature is entirely missing (mean/median/KNN can't fill)."""
    fully_missing = missing.all(axis=0)
    n_bad = int(fully_missing.sum())
    if n_bad:
        raise ValueError(
            f"impute={impute!r} cannot fill {n_bad} feature(s) with no detected "
            f"value in any sample. Drop them first with a stricter "
            f"max_missing_fraction (e.g. < 1.0), or use impute='zero'."
        )


# --------------------------------------------------------------------------- #
# Dataset assembly
# --------------------------------------------------------------------------- #


def _independent(
    dataset: Dataset,
    abundances: np.ndarray,
    feature_names: np.ndarray,
    feature_metadata: pd.DataFrame,
) -> Dataset:
    """Build a new linear Dataset sharing NO mutable state with ``dataset``."""
    return replace(
        dataset,
        abundances=abundances,
        feature_names=feature_names,
        feature_metadata=feature_metadata,
        metadata=dataset.metadata.copy(),
        scale="linear",
    )
