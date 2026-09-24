"""Peptide-modification helpers: mass decomposition, dump localization, discordance.

Building blocks for ``scripts/scratch/peptide_mod_analysis.py`` (Stage 4, peptide-level
modification analysis, raloxifene-d0 vs control). Not seeded from a ``lib/`` template
-- nothing in ``lib/`` decomposes modification masses, parses localized peptide
strings, or counts paired presence/absence discordance. The moderated one-sample t
(:func:`moderated_one_sample`) copies Smyth's ``fitFDist`` prior fit from the project
copy of the ``differential-abundance`` template (``analysis/differential_abundance.py``,
``_fit_f_distribution_prior``, incl. the limma zero-variance floor) rather than
importing it, because that module is being modified concurrently; on a balanced
paired design it reproduces the template's paired-design contrast exactly (checked in
the runner and in the tests).

Sections
--------
* **Masses.** The raloxifene adduct delta is *derived* from the elemental formula
  (C28H27NO4S monoisotopic minus 2 H = a quinone-methide-type electrophile added
  whole), not typed in.
* **Mass decomposition.** A FlashLFQ peptide id carries ONE appended total mod mass;
  :func:`decompose_mass` enumerates non-negative integer combinations of the known
  deltas (carbamidomethyl, oxidation, raloxifene adduct) within a tolerance.
* **Localized sequences.** The Limelight peptide dump writes each modification after
  its residue with a 2-decimal display mass (``ELTNC[57.02]TR``).
  :func:`parse_localized_sequence` returns the base sequence plus 1-based sites.
* **Join.** Dump rows (positional isoforms) join many-to-one to quant features on
  (base sequence, n_CAM, n_Ox, n_adduct).
* **Discordance.** Per matched pair, a peptide is *treated-only* / *control-only* /
  *both* / *neither* detected under a detection definition.
* **Small statistics helpers.** BH, exact binomial + Clopper-Pearson, Fisher 2x2
  (conditional OR + exact CI), Mantel-Haenszel over strata, Welch, Mann-Whitney with a
  Hodges-Lehmann shift + CI, paired t, and the limma-style moderated one-sample t.
"""

from __future__ import annotations

import itertools
import math
import re
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special, stats
from scipy.stats import contingency
from statsmodels.stats.contingency_tables import StratifiedTable

__script_meta__: dict[str, object] = {
    "task": "peptide-mod-analysis",
    "kind": "module",
    "provides": [
        "RALOXIFENE_ADDUCT_DELTA",
        "MOD_DELTAS",
        "MOD_CLASSES",
        "ModComposition",
        "decompose_mass",
        "decompose_feature_masses",
        "parse_localized_sequence",
        "read_peptide_dump",
        "join_dump_to_features",
        "pair_layout",
        "pair_states",
        "discordance_counts",
        "moderated_one_sample",
        "bh_adjust",
    ],
    "uses": [],
    "seeded_from": (
        "no lib/ template for mass decomposition/discordance; moderated one-sample "
        "prior fit copied from analysis/differential_abundance.py "
        "(differential-abundance@0.1 project copy)"
    ),
    "description": (
        "Peptide modification decomposition (CAM / Ox / raloxifene adduct = "
        "C28H27NO4S - 2H), Limelight-dump localized-sequence parsing and dump<->quant "
        "join on (base, composition), matched-pair presence/absence discordance, and "
        "small exact/moderated statistics helpers (BH, binomial + Clopper-Pearson, "
        "Fisher OR + exact CI, Mantel-Haenszel, Welch, Mann-Whitney + Hodges-Lehmann, "
        "paired t, limma-style moderated one-sample t)."
    ),
}

# --------------------------------------------------------------------------- #
# Masses
# --------------------------------------------------------------------------- #
#: Monoisotopic element masses (u), IUPAC / NIST.
ELEMENT_MONOISOTOPIC: dict[str, float] = {
    "C": 12.0,
    "H": 1.00782503207,
    "N": 14.0030740048,
    "O": 15.99491461956,
    "S": 31.97207100,
}

#: Raloxifene, C28H27NO4S.
RALOXIFENE_FORMULA: dict[str, int] = {"C": 28, "H": 27, "N": 1, "O": 4, "S": 1}


def monoisotopic_mass(formula: Mapping[str, int]) -> float:
    """Monoisotopic mass of an elemental formula; fails loud on an unknown element."""
    total = 0.0
    for element, count in formula.items():
        if element not in ELEMENT_MONOISOTOPIC:
            raise ValueError(f"Unknown element {element!r} in formula {formula}.")
        if count < 0:
            raise ValueError(f"Negative element count in formula {formula}.")
        total += ELEMENT_MONOISOTOPIC[element] * count
    return total


RALOXIFENE_MONOISOTOPIC: float = monoisotopic_mass(RALOXIFENE_FORMULA)
#: Raloxifene - 2H: the (di)quinone-methide-type electrophile added whole (Michael
#: addition moves the nucleophile's H onto the adduct, so the net delta is M - 2H).
RALOXIFENE_ADDUCT_DELTA: float = RALOXIFENE_MONOISOTOPIC - 2 * ELEMENT_MONOISOTOPIC["H"]
CARBAMIDOMETHYL_DELTA: float = 57.021464
OXIDATION_DELTA: float = 15.994915

MOD_DELTAS: dict[str, float] = {
    "CAM": CARBAMIDOMETHYL_DELTA,
    "Ox": OXIDATION_DELTA,
    "Ralox": RALOXIFENE_ADDUCT_DELTA,
}

UNMODIFIED = "unmodified"
CAM_ONLY = "CAM-only"
OXIDIZED = "oxidized"
ADDUCT = "raloxifene-adduct"
UNEXPLAINED = "unexplained"
#: Class precedence: adduct > oxidized > CAM-only > unmodified.
MOD_CLASSES: tuple[str, ...] = (UNMODIFIED, CAM_ONLY, OXIDIZED, ADDUCT, UNEXPLAINED)

DEFAULT_MASS_TOLERANCE = 0.001


@dataclass(frozen=True, order=True)
class ModComposition:
    """Counts of each known modification on one peptide form."""

    n_cam: int
    n_ox: int
    n_ralox: int

    def __post_init__(self) -> None:
        if min(self.n_cam, self.n_ox, self.n_ralox) < 0:
            raise ValueError(f"Negative modification count: {self}.")

    @property
    def mass(self) -> float:
        """Summed delta mass of the composition."""
        return (
            self.n_cam * CARBAMIDOMETHYL_DELTA
            + self.n_ox * OXIDATION_DELTA
            + self.n_ralox * RALOXIFENE_ADDUCT_DELTA
        )

    @property
    def label(self) -> str:
        """Human-readable label, e.g. ``"2CAM+1Ox"``; ``"none"`` when unmodified."""
        parts = [
            f"{n}{name}"
            for n, name in (
                (self.n_cam, "CAM"),
                (self.n_ox, "Ox"),
                (self.n_ralox, "Ralox"),
            )
            if n
        ]
        return "+".join(parts) if parts else "none"

    @property
    def mod_class(self) -> str:
        """Class by precedence: adduct > oxidized > CAM-only > unmodified."""
        if self.n_ralox > 0:
            return ADDUCT
        if self.n_ox > 0:
            return OXIDIZED
        if self.n_cam > 0:
            return CAM_ONLY
        return UNMODIFIED


def decompose_mass(
    total: float,
    *,
    tolerance: float = DEFAULT_MASS_TOLERANCE,
    max_cam: int = 10,
    max_ox: int = 6,
    max_ralox: int = 4,
) -> list[tuple[ModComposition, float]]:
    """All compositions whose summed delta is within ``tolerance`` of ``total``.

    Returns ``(composition, error = total - composition.mass)`` sorted by ``|error|``
    (then composition). An empty list means the mass is unexplained by the known
    deltas. Fails loud on a non-finite or negative mass or a non-positive tolerance.
    """
    if not math.isfinite(total):
        raise ValueError(f"Modification mass must be finite; got {total!r}.")
    if total < -tolerance:
        raise ValueError(f"Modification mass must be >= 0; got {total!r}.")
    if tolerance <= 0:
        raise ValueError(f"tolerance must be > 0; got {tolerance!r}.")
    solutions: list[tuple[ModComposition, float]] = []
    for n_cam, n_ox, n_ralox in itertools.product(
        range(max_cam + 1), range(max_ox + 1), range(max_ralox + 1)
    ):
        comp = ModComposition(n_cam, n_ox, n_ralox)
        error = total - comp.mass
        if abs(error) <= tolerance:
            solutions.append((comp, error))
    solutions.sort(key=lambda item: (abs(item[1]), item[0]))
    return solutions


def decompose_feature_masses(
    masses: Sequence[float],
    base_sequences: Sequence[str],
    *,
    tolerance: float = DEFAULT_MASS_TOLERANCE,
) -> pd.DataFrame:
    """Decompose each feature's total mass; one row per feature, input order.

    Columns: ``total_mod_mass``, ``n_cam``/``n_ox``/``n_ralox`` (nullable ``Int64``,
    ``<NA>`` when unexplained or ambiguous), ``composition`` label, ``mass_error``,
    ``n_solutions``, ``mod_class`` (``unexplained`` when no or more than one
    composition fits), residue counts ``n_cys``/``n_met``/``n_tyr``, and
    ``cam_exceeds_cys`` (more CAM than Cys -- a decomposition inconsistent with the
    sequence).
    """
    if len(masses) != len(base_sequences):
        raise ValueError(
            f"masses ({len(masses)}) and base_sequences ({len(base_sequences)}) differ "
            "in length."
        )
    cache: dict[float, list[tuple[ModComposition, float]]] = {}
    rows: list[dict[str, object]] = []
    for mass, base in zip(masses, base_sequences, strict=True):
        key = float(mass)
        if key not in cache:
            cache[key] = decompose_mass(key, tolerance=tolerance)
        sols = cache[key]
        n_cys = base.count("C")
        row: dict[str, object] = {
            "total_mod_mass": key,
            "n_solutions": len(sols),
            "n_cys": n_cys,
            "n_met": base.count("M"),
            "n_tyr": base.count("Y"),
        }
        if len(sols) == 1:
            comp, err = sols[0]
            row.update(
                n_cam=comp.n_cam,
                n_ox=comp.n_ox,
                n_ralox=comp.n_ralox,
                composition=comp.label,
                mass_error=err,
                mod_class=comp.mod_class,
                cam_exceeds_cys=comp.n_cam > n_cys,
            )
        else:
            row.update(
                n_cam=pd.NA,
                n_ox=pd.NA,
                n_ralox=pd.NA,
                composition="unexplained" if not sols else "ambiguous",
                mass_error=float("nan"),
                mod_class=UNEXPLAINED,
                cam_exceeds_cys=False,
            )
        rows.append(row)
    out = pd.DataFrame(rows)
    for col in ("n_cam", "n_ox", "n_ralox"):
        out[col] = out[col].astype("Int64")
    return out


# --------------------------------------------------------------------------- #
# Localized sequences (Limelight peptide dump)
# --------------------------------------------------------------------------- #
#: Display masses are rounded to 2 decimals, so |display - true| <= 0.005.
DUMP_DISPLAY_TOLERANCE = 0.006
_LOCALIZED_FULL_RE = re.compile(r"^(?:[A-Z](?:\[[0-9]+(?:\.[0-9]+)?\])?)+$")
_LOCALIZED_TOKEN_RE = re.compile(r"([A-Z])(?:\[([0-9]+(?:\.[0-9]+)?)\])?")


@dataclass(frozen=True)
class ModSite:
    """One localized modification: 1-based position in the peptide, residue, name."""

    position: int
    residue: str
    mod: str

    @property
    def label(self) -> str:
        """``"C3:Ralox"``-style label."""
        return f"{self.residue}{self.position}:{self.mod}"


@dataclass(frozen=True)
class LocalizedPeptide:
    """A dump peptide string parsed into base sequence + localized sites."""

    base_sequence: str
    sites: tuple[ModSite, ...]

    @property
    def composition(self) -> ModComposition:
        """Modification counts implied by the sites."""
        mods = [site.mod for site in self.sites]
        return ModComposition(mods.count("CAM"), mods.count("Ox"), mods.count("Ralox"))

    def sites_of(self, mod: str) -> tuple[ModSite, ...]:
        """The sites carrying modification ``mod``."""
        return tuple(site for site in self.sites if site.mod == mod)


def _mod_name_for_display_mass(value: float) -> str:
    """Map a 2-decimal display mass to a known modification name, fail loud."""
    matches = [
        name
        for name, delta in MOD_DELTAS.items()
        if abs(value - delta) <= DUMP_DISPLAY_TOLERANCE
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Display mass {value!r} matches {len(matches)} known modification(s) "
            f"({matches}); expected exactly one of {sorted(MOD_DELTAS)}."
        )
    return matches[0]


def parse_localized_sequence(sequence: str) -> LocalizedPeptide:
    """Parse ``"ELTNC[57.02]TR"`` into base ``"ELTNCTR"`` + ``(C5:CAM,)``.

    Fails loud on any character outside upper-case residues + bracketed masses, or on a
    bracketed mass that is not one of the known modifications.
    """
    if not sequence or _LOCALIZED_FULL_RE.match(sequence) is None:
        raise ValueError(f"Unparseable localized peptide sequence {sequence!r}.")
    residues: list[str] = []
    sites: list[ModSite] = []
    for match in _LOCALIZED_TOKEN_RE.finditer(sequence):
        residues.append(match.group(1))
        if match.group(2) is not None:
            sites.append(
                ModSite(
                    position=len(residues),
                    residue=match.group(1),
                    mod=_mod_name_for_display_mass(float(match.group(2))),
                )
            )
    return LocalizedPeptide(base_sequence="".join(residues), sites=tuple(sites))


def composition_key(base_sequence: str, composition: ModComposition) -> str:
    """Join key ``"<BASE>|<nCAM>|<nOx>|<nRalox>"``."""
    return (
        f"{base_sequence}|{composition.n_cam}|{composition.n_ox}|{composition.n_ralox}"
    )


# --------------------------------------------------------------------------- #
# Limelight peptide dump
# --------------------------------------------------------------------------- #
_PSMS_COL_RE = re.compile(r"^PSMs \((?P<label>[^()]+)\)$")
_THOUSANDS_INT_RE = re.compile(r"^\d{1,3}(?:,\d{3})*$|^\d+$")
MBR_FLAG = "(MBR)"


def parse_psm_count(text: str) -> int:
    """Dump PSM cell -> int: ``""`` -> 0, ``"1,234"`` -> 1234; fail loud otherwise."""
    value = text.strip()
    if value == "":
        return 0
    if _THOUSANDS_INT_RE.match(value) is None:
        raise ValueError(f"Unparseable PSM count {text!r}.")
    return int(value.replace(",", ""))


@dataclass(frozen=True)
class PeptideDump:
    """The Limelight peptide dump, one row per localized peptide string.

    Attributes
    ----------
    rows:
        ``localized_sequence``, ``base_sequence``, ``n_cam``/``n_ox``/``n_ralox``,
        ``sites`` / ``ralox_sites`` / ``ox_sites`` (``;``-joined labels), ``proteins``
        (dump ``Protein(s)``), ``unique`` (dump ``Unique == "*"``), ``key``.
    psms:
        ``(n_rows, n_samples)`` int PSM counts (blank -> 0).
    mbr:
        ``(n_rows, n_samples)`` bool, the dump's ``(MBR)`` flag.
    sample_ids:
        Column order of ``psms`` / ``mbr``.
    """

    rows: pd.DataFrame
    psms: np.ndarray
    mbr: np.ndarray
    sample_ids: tuple[str, ...]


def read_peptide_dump(
    path: str | Path,
    label_to_sample: Mapping[str, str],
    sample_order: Sequence[str],
) -> PeptideDump:
    """Read the Limelight peptide dump, mapping per-run labels to ``sample_order``.

    ``label_to_sample`` is the verified Stage-2 label map. Fails loud on: a missing
    required column, a label with no mapping or a mapping that is not a bijection onto
    ``sample_order``, a duplicated localized sequence, an unparseable PSM count or MBR
    cell, or an unparseable localized sequence.
    """
    frame = pd.read_csv(Path(path), sep="\t", dtype=str, keep_default_na=False)
    for required in ("Peptide Sequence", "Protein(s)", "Unique"):
        if required not in frame.columns:
            raise ValueError(f"{path} has no column {required!r}.")
    labels = [
        m.group("label")
        for col in frame.columns
        if (m := _PSMS_COL_RE.match(str(col))) is not None
    ]
    unmapped = sorted(set(labels) - set(label_to_sample))
    if unmapped:
        raise ValueError(f"Dump labels without a sample mapping: {unmapped}.")
    mapped = {label_to_sample[label]: label for label in labels}
    if len(mapped) != len(labels) or set(mapped) != set(sample_order):
        raise ValueError(
            f"Dump labels {sorted(labels)} do not map 1:1 onto samples "
            f"{sorted(sample_order)}."
        )
    seqs = frame["Peptide Sequence"]
    if bool(seqs.duplicated().any()):
        dup = seqs[seqs.duplicated(keep=False)].unique().tolist()[:5]
        raise ValueError(f"Duplicated dump peptide sequences: {dup}.")

    n_rows, n_samples = len(frame), len(sample_order)
    psms = np.zeros((n_rows, n_samples), dtype=np.int64)
    mbr = np.zeros((n_rows, n_samples), dtype=bool)
    for j, sample in enumerate(sample_order):
        label = mapped[sample]
        psms[:, j] = [parse_psm_count(v) for v in frame[f"PSMs ({label})"]]
        mbr_col = frame[f"MBR ({label})"].str.strip()
        bad = ~mbr_col.isin(["", MBR_FLAG])
        if bool(bad.any()):
            raise ValueError(
                f"Unexpected MBR cell(s) for {label}: {mbr_col[bad].unique()[:5]}."
            )
        mbr[:, j] = (mbr_col == MBR_FLAG).to_numpy()

    parsed = [parse_localized_sequence(s) for s in seqs]
    rows = pd.DataFrame(
        {
            "localized_sequence": seqs.to_numpy(dtype=str),
            "base_sequence": [p.base_sequence for p in parsed],
            "n_cam": [p.composition.n_cam for p in parsed],
            "n_ox": [p.composition.n_ox for p in parsed],
            "n_ralox": [p.composition.n_ralox for p in parsed],
            "sites": [";".join(s.label for s in p.sites) for p in parsed],
            "ralox_sites": [
                ";".join(s.label for s in p.sites_of("Ralox")) for p in parsed
            ],
            "ox_sites": [";".join(s.label for s in p.sites_of("Ox")) for p in parsed],
            "ralox_residues": [
                "".join(s.residue for s in p.sites_of("Ralox")) for p in parsed
            ],
            "proteins": frame["Protein(s)"].to_numpy(dtype=str),
            "unique": (frame["Unique"].str.strip() == "*").to_numpy(),
            "key": [composition_key(p.base_sequence, p.composition) for p in parsed],
        }
    )
    return PeptideDump(
        rows=rows, psms=psms, mbr=mbr, sample_ids=tuple(str(s) for s in sample_order)
    )


@dataclass(frozen=True)
class DumpJoin:
    """Many-to-one dump-row -> quant-feature join on the composition key.

    ``dump_to_feature[i]`` is the quant feature index of dump row ``i`` (``-1`` when
    unmatched); ``n_isomers[f]`` counts dump rows joined to feature ``f``.
    """

    dump_to_feature: np.ndarray
    n_isomers: np.ndarray
    unmatched_dump_keys: tuple[str, ...]
    unmatched_feature_keys: tuple[str, ...]


def join_dump_to_features(
    dump_keys: Sequence[str], feature_keys: Sequence[str]
) -> DumpJoin:
    """Join dump keys to (unique) quant feature keys; fail loud on duplicate keys."""
    index: dict[str, int] = {}
    for i, key in enumerate(feature_keys):
        if key in index:
            raise ValueError(f"Duplicate quant feature key {key!r}.")
        index[key] = i
    dump_to_feature = np.array([index.get(k, -1) for k in dump_keys], dtype=np.int64)
    n_isomers = np.bincount(
        dump_to_feature[dump_to_feature >= 0], minlength=len(feature_keys)
    )
    matched_features = set(dump_to_feature[dump_to_feature >= 0].tolist())
    return DumpJoin(
        dump_to_feature=dump_to_feature,
        n_isomers=n_isomers,
        unmatched_dump_keys=tuple(
            k for k, f in zip(dump_keys, dump_to_feature, strict=True) if f < 0
        ),
        unmatched_feature_keys=tuple(
            k for i, k in enumerate(feature_keys) if i not in matched_features
        ),
    )


def aggregate_dump_to_features(
    join: DumpJoin, dump: PeptideDump, n_features: int
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Sum PSMs / OR the MBR flag over each feature's positional isoforms.

    Returns ``(psms (n_features, n_samples) float with NaN where no dump row joined,
    mbr_any (n_features, n_samples) bool, text (n_features rows: localized_sequences,
    sites, ralox_sites, ox_sites, ralox_residues, dump_unique))``.
    """
    n_samples = dump.psms.shape[1]
    psms = np.full((n_features, n_samples), np.nan)
    mbr_any = np.zeros((n_features, n_samples), dtype=bool)
    matched = join.dump_to_feature >= 0
    feats = join.dump_to_feature[matched]
    sums = np.zeros((n_features, n_samples))
    np.add.at(sums, feats, dump.psms[matched].astype(float))
    has = join.n_isomers > 0
    psms[has] = sums[has]
    np.logical_or.at(mbr_any, feats, dump.mbr[matched])

    text_cols = ("localized_sequence", "sites", "ralox_sites", "ox_sites")
    collected: dict[str, list[list[str]]] = {
        c: [[] for _ in range(n_features)] for c in text_cols
    }
    residues: list[set[str]] = [set() for _ in range(n_features)]
    unique: list[set[bool]] = [set() for _ in range(n_features)]
    sub = dump.rows.loc[matched]
    for f, (_, row) in zip(feats.tolist(), sub.iterrows(), strict=True):
        for c in text_cols:
            if row[c]:
                collected[c][f].append(str(row[c]))
        residues[f].update(str(row["ralox_residues"]))
        unique[f].add(bool(row["unique"]))
    text = pd.DataFrame(
        {
            "dump_localized_sequences": [
                " | ".join(v) for v in collected[text_cols[0]]
            ],
            "dump_sites": [" | ".join(v) for v in collected["sites"]],
            "dump_ralox_sites": [
                " | ".join(sorted(set(v))) for v in collected["ralox_sites"]
            ],
            "dump_ox_sites": [
                " | ".join(sorted(set(v))) for v in collected["ox_sites"]
            ],
            "dump_ralox_residues": ["".join(sorted(r)) for r in residues],
            "dump_unique": [
                (next(iter(u)) if len(u) == 1 else pd.NA) if u else pd.NA
                for u in unique
            ],
        }
    )
    return psms, mbr_any, text


# --------------------------------------------------------------------------- #
# Proteins
# --------------------------------------------------------------------------- #
_MEMBER_RE = re.compile(
    r"^psvid_(?P<psvid>\d+)_sp\|(?P<accession>[^|]*)\|(?P<entry>.*)$"
)


def parse_protein_members(groups: str) -> list[tuple[str, str, str]]:
    """``;``-joined ``psvid_<n>_sp|<acc>|<ENTRY>`` members -> ``(id, acc, entry)``.

    A no-accession contaminant-form member gives ``acc == ""`` and ``entry`` from the
    middle field. Fails loud on an unparseable member.
    """
    out: list[tuple[str, str, str]] = []
    for member in groups.split(";"):
        match = _MEMBER_RE.match(member)
        if match is None:
            raise ValueError(f"Unparseable protein-group member {member!r}.")
        acc, entry = match.group("accession"), match.group("entry")
        if entry == "" and acc:
            acc, entry = "", acc
        out.append((member, acc, entry))
    return out


# --------------------------------------------------------------------------- #
# Pairs and discordance
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PairLayout:
    """Row indices of each pair's control and treated sample (same pair order)."""

    pair_ids: tuple[str, ...]
    control_idx: np.ndarray
    treated_idx: np.ndarray


def pair_layout(
    metadata: pd.DataFrame,
    *,
    pair_col: str,
    condition_col: str,
    control: str,
    treated: str,
) -> PairLayout:
    """Build the pair layout from sample metadata rows; each pair must be 1 + 1."""
    pair_ids = sorted(str(p) for p in metadata[pair_col].unique())
    cond = metadata[condition_col].to_numpy(dtype=str)
    pairs = metadata[pair_col].to_numpy(dtype=str)
    ctrl_idx: list[int] = []
    trt_idx: list[int] = []
    for pair in pair_ids:
        c = np.flatnonzero((pairs == pair) & (cond == control))
        t = np.flatnonzero((pairs == pair) & (cond == treated))
        n_in_pair = int((pairs == pair).sum())
        if len(c) != 1 or len(t) != 1 or n_in_pair != 2:
            raise ValueError(
                f"Pair {pair!r} must hold exactly one {control!r} and one {treated!r} "
                f"sample; found {len(c)} / {len(t)} of {n_in_pair}."
            )
        ctrl_idx.append(int(c[0]))
        trt_idx.append(int(t[0]))
    return PairLayout(
        pair_ids=tuple(pair_ids),
        control_idx=np.array(ctrl_idx, dtype=np.int64),
        treated_idx=np.array(trt_idx, dtype=np.int64),
    )


def pair_states(
    detected: np.ndarray, layout: PairLayout
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per pair x feature booleans ``(treated_only, control_only, both, neither)``.

    ``detected`` is ``(n_samples, n_features)`` bool (samples in metadata row order).
    """
    det = np.asarray(detected)
    if det.dtype != bool or det.ndim != 2:
        raise ValueError("detected must be a 2-D boolean array (samples x features).")
    ctrl = det[layout.control_idx]
    trt = det[layout.treated_idx]
    return trt & ~ctrl, ctrl & ~trt, ctrl & trt, ~ctrl & ~trt


def discordance_counts(detected: np.ndarray, layout: PairLayout) -> pd.DataFrame:
    """Per feature: pairs treated-only / control-only / both / neither (sum n_pairs)."""
    t_only, c_only, both, neither = pair_states(detected, layout)
    return pd.DataFrame(
        {
            "n_treated_only": t_only.sum(axis=0).astype(np.int64),
            "n_control_only": c_only.sum(axis=0).astype(np.int64),
            "n_both": both.sum(axis=0).astype(np.int64),
            "n_neither": neither.sum(axis=0).astype(np.int64),
        }
    )


# --------------------------------------------------------------------------- #
# Statistics helpers
# --------------------------------------------------------------------------- #
def bh_adjust(pvalues: Sequence[float] | np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN entries pass through and leave the family."""
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if np.any((p[ok] < 0) | (p[ok] > 1)):
        raise ValueError("p-values must lie in [0, 1].")
    m = int(ok.sum())
    if m == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def binomial_proportion(k: int, n: int, *, alpha: float = 0.05) -> dict[str, float]:
    """Proportion ``k/n`` + Clopper-Pearson CI + exact two-sided binomial p vs 0.5."""
    if n < 0 or k < 0 or k > n:
        raise ValueError(f"Need 0 <= k <= n; got k={k}, n={n}.")
    if n == 0:
        nan = float("nan")
        return {"k": 0, "n": 0, "prop": nan, "ci_low": nan, "ci_high": nan, "p": nan}
    res = stats.binomtest(k, n, 0.5)
    ci = res.proportion_ci(confidence_level=1 - alpha, method="exact")
    return {
        "k": float(k),
        "n": float(n),
        "prop": k / n,
        "ci_low": float(ci.low),
        "ci_high": float(ci.high),
        "p": float(res.pvalue),
    }


def fisher_odds_ratio(table: np.ndarray, *, alpha: float = 0.05) -> dict[str, float]:
    """Fisher exact 2x2: conditional-MLE OR + exact CI + two-sided p.

    ``table = [[a, b], [c, d]]`` with rows = (group of interest, baseline) and columns =
    (treated-only, control-only); OR > 1 = the group leans treated-only more than the
    baseline.
    """
    tab = np.asarray(table, dtype=np.int64)
    if tab.shape != (2, 2) or np.any(tab < 0):
        raise ValueError(f"Need a non-negative 2x2 table; got {tab!r}.")
    if np.any(tab.sum(axis=1) == 0) or np.any(tab.sum(axis=0) == 0):
        nan = float("nan")
        return {"odds_ratio": nan, "ci_low": nan, "ci_high": nan, "p": nan}
    res = contingency.odds_ratio(tab, kind="conditional")
    ci = res.confidence_interval(confidence_level=1 - alpha)
    p = float(stats.fisher_exact(tab, alternative="two-sided").pvalue)
    return {
        "odds_ratio": float(res.statistic),
        "ci_low": float(ci.low),
        "ci_high": float(ci.high),
        "p": p,
    }


def mantel_haenszel(tables: np.ndarray, *, alpha: float = 0.05) -> dict[str, float]:
    """Mantel-Haenszel pooled OR over strata ``(K, 2, 2)`` + CI + CMH test p.

    Strata with an empty row or column carry no information and are dropped; with
    fewer than one informative stratum the result is NaN.
    """
    tabs = np.asarray(tables, dtype=float)
    if tabs.ndim != 3 or tabs.shape[1:] != (2, 2):
        raise ValueError(f"Need a (K, 2, 2) array; got shape {tabs.shape}.")
    keep = (tabs.sum(axis=2) > 0).all(axis=1) & (tabs.sum(axis=1) > 0).all(axis=1)
    nan = float("nan")
    if int(keep.sum()) == 0:
        return {"odds_ratio": nan, "ci_low": nan, "ci_high": nan, "p": nan}
    st = StratifiedTable(np.moveaxis(tabs[keep], 0, -1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        or_mh = float(st.oddsratio_pooled)
        lo, hi = st.oddsratio_pooled_confint(alpha=alpha)
        p = float(st.test_null_odds(correction=True).pvalue)
    return {"odds_ratio": or_mh, "ci_low": float(lo), "ci_high": float(hi), "p": p}


def welch_difference(
    x: np.ndarray, y: np.ndarray, *, alpha: float = 0.05
) -> dict[str, float]:
    """Welch two-sample: ``mean(x) - mean(y)`` + Welch CI + two-sided p."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or len(y) < 2:
        nan = float("nan")
        return {"diff": nan, "ci_low": nan, "ci_high": nan, "p": nan, "df": nan}
    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y)
    se = math.sqrt(vx + vy)
    df = (vx + vy) ** 2 / (vx**2 / (len(x) - 1) + vy**2 / (len(y) - 1))
    diff = float(x.mean() - y.mean())
    crit = float(stats.t.ppf(1 - alpha / 2, df))
    p = float(2 * stats.t.sf(abs(diff / se), df))
    return {
        "diff": diff,
        "ci_low": diff - crit * se,
        "ci_high": diff + crit * se,
        "p": p,
        "df": float(df),
    }


def mann_whitney_shift(
    x: np.ndarray,
    y: np.ndarray,
    *,
    alpha: float = 0.05,
    max_pairs: int = 60_000_000,
) -> dict[str, float]:
    """Mann-Whitney two-sided p + Hodges-Lehmann shift (x - y) + distribution-free CI.

    The CI is the classic order-statistic interval of the pairwise differences
    (normal approximation to the U distribution for the rank cut-off). Also returns the
    rank-biserial correlation ``2U/(n_x n_y) - 1`` (positive = x tends higher).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nx, ny = len(x), len(y)
    nan = float("nan")
    if nx == 0 or ny == 0:
        return {
            "hl_shift": nan,
            "ci_low": nan,
            "ci_high": nan,
            "rank_biserial": nan,
            "p": nan,
        }
    if nx * ny > max_pairs:
        raise ValueError(
            f"{nx * ny} pairwise differences exceed max_pairs={max_pairs}."
        )
    res = stats.mannwhitneyu(x, y, alternative="two-sided")
    diffs = np.subtract.outer(x, y).ravel()
    diffs.sort()
    n_pairs = diffs.size
    z = float(stats.norm.ppf(1 - alpha / 2))
    k = math.floor(n_pairs / 2 - z * math.sqrt(nx * ny * (nx + ny + 1) / 12))
    k = max(k, 0)
    lo = float(diffs[k]) if k < n_pairs else nan
    hi = float(diffs[n_pairs - 1 - k]) if k < n_pairs else nan
    return {
        "hl_shift": float(np.median(diffs)),
        "ci_low": lo,
        "ci_high": hi,
        "rank_biserial": float(2 * res.statistic / (nx * ny) - 1),
        "p": float(res.pvalue),
    }


def paired_t(diffs: np.ndarray, *, alpha: float = 0.05) -> dict[str, float]:
    """One-sample t on paired differences: mean + t CI + two-sided p (df = n - 1)."""
    d = np.asarray(diffs, dtype=float)
    if d.ndim != 1 or not np.all(np.isfinite(d)):
        raise ValueError("diffs must be a finite 1-D array.")
    n = len(d)
    nan = float("nan")
    if n < 2:
        return {"mean": nan, "ci_low": nan, "ci_high": nan, "p": nan, "df": nan}
    mean = float(d.mean())
    se = float(d.std(ddof=1) / math.sqrt(n))
    crit = float(stats.t.ppf(1 - alpha / 2, n - 1))
    p = float(2 * stats.t.sf(abs(mean / se), n - 1)) if se > 0 else nan
    return {
        "mean": mean,
        "ci_low": mean - crit * se,
        "ci_high": mean + crit * se,
        "p": p,
        "df": float(n - 1),
    }


def sign_test(n_positive: int, n_negative: int) -> float:
    """Exact two-sided sign test p (ties excluded by the caller)."""
    n = n_positive + n_negative
    if n == 0:
        return float("nan")
    return float(stats.binomtest(n_positive, n, 0.5).pvalue)


_ZERO_VARIANCE_FLOOR_FRACTION = 1e-5
_MIN_FEATURES_FOR_PRIOR = 3


def fit_f_prior(sigma2: np.ndarray, df: int) -> tuple[float, float]:
    """Smyth's ``fitFDist``: scaled-inverse-chi2 prior ``(s0^2, d0)`` by moments.

    Copied from ``analysis/differential_abundance.py::_fit_f_distribution_prior``
    (limma zero-variance floor at 1e-5 x median; ``d0 = inf`` when the variances are
    no more spread than sampling alone predicts).
    """
    s2 = np.asarray(sigma2, dtype=float)
    if s2.size < _MIN_FEATURES_FOR_PRIOR or not np.all(np.isfinite(s2)):
        raise ValueError(
            f"Need >= {_MIN_FEATURES_FOR_PRIOR} finite variances to fit the prior."
        )
    if df <= 0:
        raise ValueError(f"df must be positive; got {df}.")
    x = np.maximum(s2, 0.0)
    median = float(np.median(x))
    if median <= 0:
        raise ValueError("More than half of the variances are zero; prior unfittable.")
    x = np.maximum(x, _ZERO_VARIANCE_FLOOR_FRACTION * median)
    z = np.log(x)
    z_mean = float(z.mean())
    z_var = float(z.var(ddof=1))
    dig = float(special.digamma(df / 2))
    target = z_var - float(special.polygamma(1, df / 2))

    def s0_inf() -> float:
        return float(np.exp(z_mean - dig + np.log(df / 2)))

    if target <= 0:
        return s0_inf(), float("inf")

    def f(v: float) -> float:
        return float(special.polygamma(1, v / 2)) - target

    if f(1e6) > 0:
        return s0_inf(), float("inf")
    d0 = float(optimize.brentq(f, 1e-6, 1e6))
    log_s0 = (
        z_mean
        - dig
        + float(np.log(df / 2))
        + float(special.digamma(d0 / 2))
        - float(np.log(d0 / 2))
    )
    return float(np.exp(log_s0)), d0


@dataclass(frozen=True)
class ModeratedOneSample:
    """Per-feature moderated one-sample t on ``n`` paired differences."""

    mean: np.ndarray
    se: np.ndarray
    t: np.ndarray
    p: np.ndarray
    q: np.ndarray
    ci_low: np.ndarray
    ci_high: np.ndarray
    s0_sq: float
    d0: float
    df_residual: int
    n: int


def moderated_one_sample(
    diffs: np.ndarray, *, alpha: float = 0.05
) -> ModeratedOneSample:
    """limma-style moderated t for mean(diff) != 0, per column of ``(n, n_features)``.

    Equivalent (for a balanced two-condition paired design) to the ``condition +
    pair`` moderated contrast: ``sigma^2`` is shrunk toward the fitted prior,
    ``se = sqrt(s_post^2 / n)``, ``t`` on ``df_residual + d0`` degrees of freedom, BH
    over the features. Fails loud on any non-finite difference.
    """
    d = np.asarray(diffs, dtype=float)
    if d.ndim != 2 or d.shape[0] < 2:
        raise ValueError("diffs must be (n >= 2, n_features).")
    if not np.all(np.isfinite(d)):
        raise ValueError("diffs contain non-finite values; subset first.")
    n = d.shape[0]
    df = n - 1
    mean = d.mean(axis=0)
    s2 = d.var(axis=0, ddof=1)
    s0_sq, d0 = fit_f_prior(s2, df)
    if math.isinf(d0):
        post = np.full_like(s2, s0_sq)
        df_total = float("inf")
    else:
        post = (d0 * s0_sq + df * s2) / (d0 + df)
        df_total = d0 + df
    se = np.sqrt(post / n)
    t = mean / se
    if math.isinf(df_total):
        p = 2 * stats.norm.sf(np.abs(t))
        crit = float(stats.norm.ppf(1 - alpha / 2))
    else:
        p = 2 * stats.t.sf(np.abs(t), df_total)
        crit = float(stats.t.ppf(1 - alpha / 2, df_total))
    p = np.asarray(p, dtype=float)
    return ModeratedOneSample(
        mean=mean,
        se=se,
        t=t,
        p=p,
        q=bh_adjust(p),
        ci_low=mean - crit * se,
        ci_high=mean + crit * se,
        s0_sq=s0_sq,
        d0=d0,
        df_residual=df,
        n=n,
    )


def pair_level_log_odds_ratios(
    tables: np.ndarray, *, haldane: float = 0.5
) -> np.ndarray:
    """Per-stratum log OR of ``(K, 2, 2)`` tables, Haldane-corrected when a cell is 0.

    Row 0 = group of interest, row 1 = baseline; column 0 = treated-only, column 1 =
    control-only. A stratum with an empty row gives NaN.
    """
    tabs = np.asarray(tables, dtype=float)
    if tabs.ndim != 3 or tabs.shape[1:] != (2, 2):
        raise ValueError(f"Need a (K, 2, 2) array; got shape {tabs.shape}.")
    out = np.full(tabs.shape[0], np.nan)
    for k, tab in enumerate(tabs):
        if np.any(tab.sum(axis=1) == 0):
            continue
        t = tab + haldane if np.any(tab == 0) else tab
        out[k] = math.log((t[0, 0] * t[1, 1]) / (t[0, 1] * t[1, 0]))
    return out


def pair_level_class_shift(
    diffs: np.ndarray,
    classes: np.ndarray,
    baseline: str,
    *,
    alpha: float = 0.05,
) -> dict[str, dict[str, object]]:
    """Per class: per-pair shift ``mean(class) - mean(baseline)`` of the paired
    differences ``(n_pairs, n_features)``, then a one-sample t across pairs.

    The pair is the unit of replication, so this is the replication-honest class test
    (df = n_pairs - 1); peptide-level tests treat features as replicates of what is a
    run-level effect.
    """
    d = np.asarray(diffs, dtype=float)
    cls = np.asarray(classes, dtype=str)
    if d.ndim != 2 or d.shape[1] != len(cls):
        raise ValueError("diffs must be (n_pairs, n_features) matching classes.")
    base = cls == baseline
    if not base.any():
        raise ValueError(f"No features of baseline class {baseline!r}.")
    base_mean = np.nanmean(d[:, base], axis=1)
    out: dict[str, dict[str, object]] = {}
    for c in sorted(set(cls.tolist()) - {baseline}):
        sel = cls == c
        shift = np.nanmean(d[:, sel], axis=1) - base_mean
        res = paired_t(shift, alpha=alpha)
        out[c] = {
            "n_features": int(sel.sum()),
            "per_pair_shift": shift.tolist(),
            "n_pairs_positive": int((shift > 0).sum()),
            **res,
        }
    return out
