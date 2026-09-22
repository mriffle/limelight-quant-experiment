"""Stage-2 exploration of the quant matrices and Limelight dumps (scratch, read-only).

Prints a structural profile of each data file and resolves the Limelight dump run
labels by value correspondence with the FlashLFQ protein intensities. Writes
``results/stage2/`` summary tables. Run from the project root:

    ./.venv/bin/python scripts/scratch/stage2_explore.py
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("data")
OUT = Path("results/stage2")
SAMPLES = pd.read_csv("results/metadata/samples.tsv", sep="\t", dtype=str)
ID2SAMPLE = dict(zip(SAMPLES["search_scan_file_id"], SAMPLES["sample_id"], strict=True))


def read_text_table(path: Path) -> pd.DataFrame:
    """Read every cell as text; keep empty cells as '' (never type-infer ids)."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, quoting=3)


def intensity_cols(df: pd.DataFrame, prefix: str) -> dict[str, str]:
    """Map sample_id -> column for ``<prefix>_search_scan_file_id_<id>`` columns."""
    out = {}
    for col in df.columns:
        m = re.fullmatch(rf"{re.escape(prefix)}_search_scan_file_id_(\d+)", col)
        if m:
            out[ID2SAMPLE[m.group(1)]] = col
    return out


def parse_num(s: pd.Series) -> pd.Series:
    """Parse display-formatted numbers ('1,745', '1.19e+9'); non-numeric -> NaN."""
    return pd.to_numeric(s.str.replace(",", "", regex=False), errors="coerce")


def sig3(x: np.ndarray) -> np.ndarray:
    """Round to 3 significant figures (Limelight's display precision)."""
    out = np.zeros_like(x, dtype=float)
    nz = x > 0
    mag = np.floor(np.log10(x[nz]))
    out[nz] = np.round(x[nz] / 10**mag, 2) * 10**mag
    return out


def profile_values(name: str, mat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample, col in mat.items():
        v = col.to_numpy(dtype=float)
        pos = v[v > 0]
        rows.append(
            {
                "file": name,
                "sample_id": sample,
                "n": len(v),
                "n_nan": int(np.isnan(v).sum()),
                "n_zero": int((v == 0).sum()),
                "n_positive": int((v > 0).sum()),
                "n_negative": int((v < 0).sum()),
                "frac_zero": round(float((v == 0).mean()), 4),
                "min_pos": float(pos.min()),
                "median_pos": float(np.median(pos)),
                "max": float(pos.max()),
                "log10_range": round(float(np.log10(pos.max() / pos.min())), 2),
                "median_log2_pos": round(float(np.median(np.log2(pos))), 3),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    order = SAMPLES.sort_values(["batch", "run_position_within_batch"])["sample_id"]

    # ---------------- protein quants ----------------
    prot = read_text_table(DATA / "protein-quants.tsv")
    pcols = intensity_cols(prot, "Intensity")
    pmat = prot[[pcols[s] for s in order]].replace("NaN", np.nan).apply(pd.to_numeric, errors="raise")
    pmat.columns = list(order)
    print("PROTEIN quants shape", prot.shape, "| unique ids", prot["Protein Groups"].nunique())
    pid = prot["Protein Groups"]
    print("  id pattern psvid_<n>_<db>|<acc>|<entry>:",
          pid.str.fullmatch(r"psvid_\d+_\w+\|[^|]+\|\S+").sum(), "of", len(pid))
    print("  ids containing ';' (multi-member groups):", pid.str.contains(";").sum())
    for tok in ("DECOY", "REV", "rev_", "CON", "contam", "cRAP", "TREMBL", "tr|"):
        print(f"  ids containing {tok!r}:", pid.str.contains(tok, regex=False).sum())
    print("  Gene Name non-empty:", (prot["Gene Name"] != "").sum(),
          "| Organism non-empty:", (prot["Organism"] != "").sum())
    print("  duplicate id rows:", pid.duplicated().sum(),
          "| duplicate value rows:", pmat.duplicated().sum())
    print("  all-zero rows:", int((pmat == 0).all(axis=1).sum()),
          "| complete (no zero) rows:", int((pmat > 0).all(axis=1).sum()))
    print("  non-integer-looking values present:", (pmat % 1 != 0).any().any())
    pprof = profile_values("protein", pmat)
    print(pprof.to_string(index=False))
    det_counts = (pmat > 0).sum(axis=1).value_counts().sort_index()
    print("  proteins by #samples quantified:", det_counts.to_dict())

    # ---------------- peptide quants ----------------
    pep = read_text_table(DATA / "peptide-quants.tsv")
    icol = intensity_cols(pep, "Intensity")
    dcol = intensity_cols(pep, "Detection Type")
    emat = pep[[icol[s] for s in order]].apply(pd.to_numeric, errors="raise")
    emat.columns = list(order)
    dmat = pep[[dcol[s] for s in order]]
    dmat.columns = list(order)
    print("\nPEPTIDE quants shape", pep.shape, "| unique Sequence", pep["Sequence"].nunique(),
          "| unique Base Sequence", pep["Base Sequence"].nunique())
    mods = pep["Sequence"].str.findall(r"\[[^\]]*\]").explode().value_counts()
    print("  modification tokens:", mods.head(15).to_dict())
    print("  Protein Groups with ';':", pep["Protein Groups"].str.contains(";").sum(),
          "| Gene Names non-empty:", (pep["Gene Names"] != "").sum())
    print("  duplicate Sequence rows:", pep["Sequence"].duplicated().sum())
    print("  all-zero rows:", int((emat == 0).all(axis=1).sum()),
          "| complete rows:", int((emat > 0).all(axis=1).sum()))
    stacked = pd.DataFrame({"type": dmat.to_numpy().ravel(), "intensity": emat.to_numpy().ravel()})
    xt = stacked.assign(positive=stacked["intensity"] > 0).groupby("type")["positive"].agg(
        n="size", n_positive="sum"
    )
    print("  detection type vs intensity>0:\n", xt.to_string())
    eprof = profile_values("peptide", emat)
    print(eprof.to_string(index=False))
    pd.concat([pprof, eprof]).to_csv(OUT / "value_profile.tsv", sep="\t", index=False)
    xt.to_csv(OUT / "peptide_detection_type_vs_intensity.tsv", sep="\t")

    # ---------------- Limelight protein dump: resolve run labels ----------------
    dump = read_text_table(DATA / "protein-limelight-table-dump.txt")
    labels = [re.search(r"\(([^()]*)\)$", c).group(1) for c in dump.columns if c.startswith("PSMs (")]
    print("\nPROTEIN DUMP shape", dump.shape, "| labels", labels)
    print("  multi-protein rows (',' or ';' in Protein(s)):",
          dump["Protein(s)"].str.contains("[,;]").sum())
    lim_tokens = pd.Series(
        dump[[f"Quant (Limelight) ({lab})" for lab in labels]].to_numpy().ravel()
    )
    print("  non-numeric Quant (Limelight) tokens:",
          lim_tokens[parse_num(lim_tokens).isna()].value_counts().to_dict())
    ff_tokens = pd.Series(dump[[f"Quant (FlashLFQ) ({lab})" for lab in labels]].to_numpy().ravel())
    print("  non-numeric Quant (FlashLFQ) tokens:",
          ff_tokens[parse_num(ff_tokens).isna()].value_counts().to_dict())
    acc = prot["Protein Groups"].str.replace(r"^psvid_\d+_", "", regex=True)
    pq = pmat.set_axis(acc.to_numpy())
    common = pd.Index(dump["Protein(s)"]).intersection(pq.index)
    print("  dump proteins matched to protein-quants by accession:", len(common),
          "of", len(dump), "(quants:", len(pq), ")")
    dq = dump.set_index("Protein(s)")
    score = pd.DataFrame(index=labels, columns=list(order), dtype=float)
    for lab in labels:
        d = parse_num(dq.loc[common, f"Quant (FlashLFQ) ({lab})"]).to_numpy(dtype=float)
        for s in order:
            q = sig3(pq.loc[common, s].to_numpy(dtype=float))
            ok = ~np.isnan(d)
            score.loc[lab, s] = float(np.isclose(d[ok], q[ok], rtol=1e-6).mean())
    print("  fraction of proteins whose dump FlashLFQ == sig3(quant), label x sample:")
    print(score.round(3).to_string())
    best = score.idxmax(axis=1)
    mapping = pd.DataFrame({"label": labels, "sample_id": best.to_numpy(),
                            "match_fraction": score.max(axis=1).to_numpy(),
                            "second_best": score.apply(lambda r: sorted(r)[-2], axis=1).to_numpy()})
    mapping = mapping.merge(SAMPLES[["sample_id", "condition", "batch", "search_scan_file_id"]])
    print(mapping.to_string(index=False))
    print("  bijective:", mapping["sample_id"].is_unique)
    mapping.to_csv(OUT / "limelight_label_map.tsv", sep="\t", index=False)
    score.to_csv(OUT / "limelight_label_match_scores.tsv", sep="\t")

    # PSM / NSAF columns: sanity
    psm = dq[[f"PSMs ({lab})" for lab in labels]].apply(parse_num)
    nsaf = dq[[f"NSAF ({lab})" for lab in labels]].apply(parse_num)
    print("  PSM NaN:", int(psm.isna().sum().sum()), "| PSM zero:", int((psm == 0).sum().sum()),
          "| PSM non-integer:", int((psm % 1 != 0).sum().sum()))
    print("  NSAF col sums:", nsaf.sum().round(3).to_dict())
    print("  NSAF min positive:", float(nsaf[nsaf > 0].min().min()),
          "| NSAF zero:", int((nsaf == 0).sum().sum()))
    print("  Protein Group Number unique:", dump["Protein Group Number"].nunique())

    # ---------------- Limelight peptide dump: check same label map ----------------
    pdump = read_text_table(DATA / "peptide-limelight-table-dump.txt")
    print("\nPEPTIDE DUMP shape", pdump.shape)
    print("  Unique col values:", pdump["Unique"].value_counts().to_dict())
    mbr_vals = pd.Series(pdump[[f"MBR ({lab})" for lab in labels]].to_numpy().ravel())
    print("  MBR values:", mbr_vals.value_counts().head(5).to_dict())
    sg = [c for c in pdump.columns if c.startswith("Shared group")]
    print("  Shared-group values:", pd.Series(pdump[sg].to_numpy().ravel()).value_counts().head(5).to_dict())
    pdump["base"] = pdump["Peptide Sequence"].str.replace(r"\[[^\]]*\]", "", regex=True)
    pep_sum = emat.groupby(pep["Base Sequence"].to_numpy()).sum()
    pd_sum = pdump.groupby("base")[[f"Quant ({lab})" for lab in labels]].agg(
        lambda s: parse_num(s).sum(min_count=1)
    )
    common_p = pd_sum.index.intersection(pep_sum.index)
    print("  base sequences in both:", len(common_p), "| dump:", len(pd_sum), "| quants:", len(pep_sum))
    for _, r in mapping.iterrows():
        d = pd_sum.loc[common_p, f"Quant ({r['label']})"].to_numpy(dtype=float)
        q = pep_sum.loc[common_p, r["sample_id"]].to_numpy(dtype=float)
        ok = ~np.isnan(d) & (q > 0) & (d > 0)
        rho = pd.Series(np.log(d[ok])).corr(pd.Series(np.log(q[ok])))
        rel = np.median(np.abs(d[ok] / q[ok] - 1))
        print(f"  {r['label']:>9} -> {r['sample_id']}: n={ok.sum()} log-r={rho:.4f} median|rel diff|={rel:.4f}")

    # ---------------- cross-level consistency ----------------
    prot_in_pep = pep["Protein Groups"].str.split(";").explode().unique()
    print("\nprotein-quants ids appearing in peptide Protein Groups:",
          int(prot["Protein Groups"].isin(prot_in_pep).sum()), "of", len(prot))
    _ = math  # (kept for quick interactive use)


if __name__ == "__main__":
    main()
