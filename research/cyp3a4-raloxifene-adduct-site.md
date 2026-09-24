---
# Research-finding document — external knowledge, structured like a finding but for the literature/tools.
# Spec: doc 04.4. Stored in research/<slug>.md. A research finding without verified references is NOT accepted.

topic: "CYP3A4 raloxifene/mechanism-based-inactivator adduction site (Cys239) and mapping of the tryptic peptide GFCMFDMECHK"
type: protein
created: 2026-09-24
updated: 2026-09-24

status: reviewed
reviewed_by: research-reviewer
reviewed_date: 2026-09-24

entities:
  - { db: uniprot, id: "P08684", label: "CYP3A4 (CP3A4_HUMAN)" }
  - { db: pdb, id: "5VCC", label: "Crystal structure of human CYP3A4 bound to glycerol" }

references:
  - id: "uniprot:P08684"
    type: url
    claim: "Canonical human CYP3A4 sequence (503 aa, SV4); peptide GFCMFDMECHK maps to residues 56-66 with Cys58 and Cys64; residue 239 is Cys; heme axial (thiolate) ligand is Cys442; initiator Met1 is removed but chain/numbering runs 2-503 with no cleaved signal peptide (N-terminal residues 2-22 annotated as a single-pass transmembrane/signal-anchor helix, ECO:0000255 predicted); PDB 5VCC-derived secondary structure places residue 239 in a non-helical loop between Helix(230-235) and Helix(243-260)."
    verified: true
    verified_by: research-reviewer
  - id: "doi:10.1021/tx700037e"
    type: doi
    claim: "Baer BR, Wienkers LC, Rock DA (2007) 'Time-dependent inactivation of P450 3A4 by raloxifene: identification of Cys239 as the site of apoprotein alkylation,' Chem Res Toxicol 20(6):954-964. Mass spectrometry of proteinase-K-digested P450 3A4 apoprotein/peptides showed a +471 Da mass shift localized to Cys239, consistent with covalent adduction by the raloxifene diquinone methide; pre-alkylation of Cys239 with iodoacetamide or N-(1-pyrene)iodoacetamide blocked raloxifene-mediated time-dependent inactivation."
    verified: true
    verified_by: research-reviewer
  - id: "pmid:17497897"
    type: pmid
    claim: "PubMed record confirms citation (authors Baer BR, Wienkers LC, Rock DA; title; Chem Res Toxicol 2007;20(6):954-64; DOI 10.1021/tx700037e; epub 2007 May 12) and provides the full abstract used above."
    verified: true
    verified_by: research-reviewer
  - id: "doi:10.1021/tx700207u"
    type: doi
    claim: "Pearson JT, Wahlstrom JL, Dickmann LJ, Kumar S, Halpert JR, Wienkers LC, Foti RS, Rock DA (2007) 'Differential time-dependent inactivation of P450 3A4 and P450 3A5 by raloxifene: a key role for C239 in quenching reactive intermediates,' Chem Res Toxicol 20(12):1778-1786. CYP3A5 carries Ser at the position corresponding to CYP3A4 Cys239 and shows only reversible (not time-dependent) inhibition by raloxifene; a CYP3A4 C239A mutant likewise loses time-dependent inactivation, confirming Cys239 as the nucleophile responsible for raloxifene-mediated TDI."
    verified: true
    verified_by: research-reviewer
  - id: "pmid:18001057"
    type: pmid
    claim: "PubMed record confirms citation and abstract for the Pearson et al. 2007 CYP3A4/CYP3A5 C239 follow-up study."
    verified: true
    verified_by: research-reviewer
  - id: "doi:10.1124/mol.112.080739"
    type: doi
    claim: "VandenBrink BM, Davis JA, Pearson JT, Foti RS, Wienkers LC, Rock DA (2012) 'Cytochrome P450 architecture and cysteine nucleophile placement impact raloxifene-mediated mechanism-based inactivation,' Mol Pharmacol 82(5):835-842. Proteolytic digests of recombinant CYP2C8 and CYP3A4 Supersomes after raloxifene time-dependent inactivation showed adducts localized to Cys225 (CYP2C8) and Cys239 (CYP3A4); only CYP2C8 and CYP3A4 among the P450s tested (1A2, 2C8, 2C9, 2C19, 2D6, 2E1, 3A4, 3A5) possess an accessible active-site/channel cysteine and were the only ones showing raloxifene time-dependent inhibition."
    verified: true
    verified_by: research-reviewer
  - id: "pmid:22859722"
    type: pmid
    claim: "PubMed record confirms citation and abstract for VandenBrink et al. 2012 Mol Pharmacol paper."
    verified: true
    verified_by: research-reviewer
  - id: "doi:10.1021/acs.biochem.7b00334"
    type: doi
    claim: "Sevrioukova IF (2017) 'High-Level Production and Properties of the Cysteine-Depleted Cytochrome P450 3A4,' Biochemistry 56(24):3058-3067 (PMCID PMC5858725). CYP3A4 has six non-heme-ligating cysteines (Cys58, Cys64, Cys98, Cys239, Cys377, Cys468); in the crystal structure (e.g. PDB 5VCC and related entries) Cys239 is located in the G'-G connecting loop, part of the F-G structural fragment; the engineered C239T substitution promotes new H-bonding to the Phe241 amide nitrogen and a recruited water molecule, and is proposed (speculatively, in the source's own hedged wording) to lower the motional freedom of the B-C and F-G fragments."
    verified: true
    verified_by: research-reviewer
  - id: "pmid:28590129"
    type: pmid
    claim: "PubMed record confirms citation and abstract for Sevrioukova 2017 Biochemistry paper on cysteine-depleted CYP3A4."
    verified: true
    verified_by: research-reviewer
  - id: "pdb:5VCC"
    type: url
    claim: "RCSB PDB entry 5VCC, 'Crystal structure of human CYP3A4 bound to glycerol' (Sevrioukova IF), initial release 2017-05-31; this is the structure UniProt cites as the evidence source for the Helix(230-235)/Helix(243-260) secondary-structure features flanking residue 239 in P08684."
    verified: true
    verified_by: research-reviewer
---

# Research: CYP3A4 raloxifene/mechanism-based-inactivator adduction site (Cys239) and mapping of the tryptic peptide GFCMFDMECHK

## Summary

Using the canonical UniProt P08684 (CP3A4_HUMAN, SV4, 503 aa) sequence, the tryptic peptide `GFCMFDMECHK` maps unambiguously to residues **56–66**, with its two cysteines at **Cys58** (peptide position 3) and **Cys64** (peptide position 9 — not position 10 as assumed in the task; H is position 10). **This peptide does NOT contain Cys239** and its first cysteine does **not** correspond to Cys239; Cys239 lies in a different, non-overlapping tryptic peptide (residues 213–243, `FDFLDPFFLSITVFPFLIPILEVLNICVFPR`) far downstream in the sequence [uniprot:P08684]. The premise in the task that "the first Cys of that peptide corresponds to Cys239" is therefore **not correct** under standard trypsin digestion of the UniProt canonical sequence.

Separately, the literature citation is verified as real and accurately described: Baer BR, Wienkers LC, Rock DA, "Time-dependent inactivation of P450 3A4 by raloxifene: identification of Cys239 as the site of apoprotein alkylation," *Chem Res Toxicol* 2007;20(6):954–964, DOI 10.1021/tx700037e, PMID 17497897 [doi:10.1021/tx700037e; pmid:17497897]. That paper used **proteinase K** (not trypsin) digestion and mass spectrometry to localize a **+471 Da** covalent mass shift to **Cys239** of P450 3A4, attributed to alkylation by the raloxifene diquinone methide reactive intermediate; the exact proteinase-K-derived peptide sequence and whether recombinant enzyme or human liver microsomes (HLM) were used are **not stated in the retrievable abstract** and could not be confirmed from full text (paywalled) — this specific detail is marked unverified below. Independent follow-up work from the same group (recombinant CYP3A4 Supersomes) reproduces the Cys239 finding [doi:10.1124/mol.112.080739] and a related paper shows CYP3A5 (which has Ser instead of Cys at the homologous position) lacks time-dependent inactivation by raloxifene, and a CYP3A4 C239A mutant likewise loses it [doi:10.1021/tx700207u]. Note that Baer 2007, Pearson 2007 and VandenBrink 2012 share authors (Wienkers, Rock) and originate from the same group, so these are replications rather than fully independent confirmations; genuinely external corroboration comes from Sevrioukova 2017, which cites the Baer and Pearson work and states that Cys239 "is solely alkylated by the raloxifene metabolite" [doi:10.1021/acs.biochem.7b00334].

UniProt numbers CYP3A4 as the full 503-residue translated precursor; Met1 is proteolytically removed (experimentally confirmed) but the mature chain retains its original numbering (Chain 2–503) — CYP3A4 has **no cleaved signal peptide**, only an N-terminal single-pass transmembrane/signal-anchor helix (residues 2–22, computationally predicted, ECO:0000255) that remains part of the mature protein. Consequently, literature residue numbers such as "Cys239" map directly onto UniProt numbering with **no offset** — confirmed directly, since UniProt P08684 residue 239 is indeed Cys [uniprot:P08684].

Structurally, per a crystallographic study of cysteine-depleted CYP3A4 (Sevrioukova 2017), Cys239 is one of CYP3A4's six non-heme-ligating cysteines and sits in the **G′–G connecting loop, part of the F–G structural fragment** [doi:10.1021/acs.biochem.7b00334]; it is distinct in both sequence and (per UniProt's PDB-derived secondary-structure annotation) local fold from the heme-thiolate axial ligand **Cys442** [uniprot:P08684], which lies near the C-terminal "Cys pocket" heme-binding motif. No source found in this search reports a non-raloxifene mechanism-based inactivator alkylating CYP3A4 apoprotein at Cys239; all located Cys239-apoprotein-adduct literature concerns raloxifene specifically.

## Detailed findings

### 1. Peptide mapping against UniProt P08684

- UniProt accession **P08684** (CP3A4_HUMAN, cytochrome P450 3A4, *Homo sapiens*), reviewed/Swiss-Prot entry, **sequence version 4** (integrated 23-JAN-2007), **entry version 254**, last annotation update 2026-09-02, retrieved 2026-09-24 via `https://rest.uniprot.org/uniprotkb/P08684` (FASTA and JSON) [uniprot:P08684]. Canonical sequence length: **503 aa**.
- The exact substring `GFCMFDMECHK` occurs **exactly once** in the canonical sequence, at 1-based residues **56–66** (context: `...FLGNILSYHK | GFCMFDMECHK | KYGKVWGFYD...`) [uniprot:P08684].
- Within the peptide (G¹F²C³M⁴F⁵D⁶M⁷E⁸C⁹H¹⁰K¹¹), the cysteines fall at peptide positions **3 and 9**:
  - Peptide C3 → absolute residue **Cys58**
  - Peptide C9 → absolute residue **Cys64**
  - Note: the task described these as "C3 and C10" — position 10 in the peptide is **His**, not Cys; the second cysteine is at position 9. This is a simple off-by-one; the absolute residue numbers (58 and 64) are unaffected by this indexing slip.
- **UniProt numbering convention**: the entry's `Chain` feature spans residues 2–503 ("Cytochrome P450 3A4"), and the `Initiator methionine` at position 1 is annotated as experimentally **removed** (evidence: PMID 3243766, PMID 3898085, per UniProt) [uniprot:P08684]. Critically, UniProt does **not** renumber the mature protein after removing Met1 — position 2 remains "2," not "1." There is **no `Signal peptide` feature** annotated for P08684; instead, residues 2–22 are annotated as a single-pass **Transmembrane** ("Helical") region (evidence code ECO:0000255, i.e., computationally predicted, not experimentally confirmed) [uniprot:P08684], consistent with the well-established membrane-anchored (not signal-peptide-cleaved) topology of microsomal P450s. Practically, this means literature residue numbers for CYP3A4 (e.g., "Cys239," "Thr309," "Cys442") map **directly and without offset** onto UniProt P08684 numbering — there is no reconciliation needed for this entry, and the task's caution about "older literature using a different convention" does not apply to the sources checked here (all consulted papers use the same 503-aa numbering as UniProt, corroborated directly by UniProt residue 239 = Cys, residue 442 = Cys — see below).

### 2. Is the peptide's first Cys = Cys239? What does Baer et al. 2007 actually say?

- **No.** Cys58 (the peptide's first cysteine, at position 56–66 in the sequence) is not Cys239. An in-silico full tryptic digest (cleave C-terminal to K/R, not before P; no missed cleavages) of the P08684 canonical sequence places **Cys239 in the peptide spanning residues 213–243**, sequence `FDFLDPFFLSITVFPFLIPILEVLNICVFPR` — a completely different, non-overlapping peptide from `GFCMFDMECHK` (56–66) [uniprot:P08684, in-silico digest performed for this research]. The peptide `GFCMFDMECHK` instead contains Cys58 and Cys64, near the N-terminal, membrane-proximal region of the protein (immediately following the substrate-recognition/catalytic-domain start, downstream of the transmembrane anchor 2–22), not in the region containing Cys239.
- **Citation verification**: The paper exists exactly as the user described: **Baer BR, Wienkers LC, Rock DA**, "**Time-dependent inactivation of P450 3A4 by raloxifene: identification of Cys239 as the site of apoprotein alkylation**," *Chemical Research in Toxicology*, **2007**, **20**(6), **954–964**; Epub 2007 May 12; DOI **10.1021/tx700037e**; PMID **17497897** [doi:10.1021/tx700037e; pmid:17497897]. Confirmed via PubMed ESummary/EFetch (NCBI E-utilities) against PMID 17497897.
- **What the paper (per its PubMed abstract) actually claims**:
  - Mass spectrometry showed **a single equivalent of raloxifene bound to the intact P450 apoprotein**.
  - Digestion was performed with **proteinase K** (not trypsin); mass analysis of the resulting peptides localized the covalent drug adduct to residue **Cys239**.
  - A **mass shift of +471 Da** was observed on both the intact protein and the modified peptide, relative to unmodified controls.
  - The chemistry proposed: raloxifene is bioactivated to a **diquinone methide** reactive intermediate, which is then attacked by the **sulfur nucleophile of Cys239**, forming the covalent apoprotein adduct responsible for time-dependent (mechanism-based) inactivation.
  - Functional confirmation: pre-treating P450 3A4 with cysteine-alkylating reagents **iodoacetamide** or **N-(1-pyrene)iodoacetamide** — which the authors state modify Cys239 "exclusively" — **prevented** subsequent raloxifene-mediated time-dependent inactivation, supporting Cys239 as the functionally relevant site [doi:10.1021/tx700037e; pmid:17497897].
  - **Not verifiable from the abstract alone (marked unverified)**: (a) the exact proteinase-K-derived peptide sequence carrying the +471 Da adduct — the abstract does not give the sequence, and the full text is paywalled (ACS; a WebFetch attempt returned HTTP 403); (b) whether the enzyme source was recombinant, baculovirus-expressed CYP3A4 ("Supersomes") or human liver microsomes (HLM) — the abstract only says "P450 3A4" / "P450 apoprotein" generically. Indirect (not conclusive) support that this lab's standard practice is recombinant Supersomes comes from the same group's related 2012 paper, which explicitly used "proteolytic digests of CYP2C8 and CYP3A4 Supersomes" to map the analogous raloxifene-Cys adducts [doi:10.1124/mol.112.080739], but this is a different publication and does not by itself confirm the 2007 paper's exact system.
  - Because `GFCMFDMECHK` is not the Cys239-containing tryptic peptide, and because Baer et al. used proteinase K (a non-specific protease) rather than trypsin, `GFCMFDMECHK` cannot be (and is not claimed by this research to be) the peptide reported in Baer et al. 2007. If the user's dataset observed `GFCMFDMECHK`, it reflects the **Cys58/Cys64** region of CYP3A4, not the raloxifene-Cys239 adduct site.
- **Independent corroboration that Cys239 (UniProt numbering) is correct**:
  - Pearson JT et al. (2007), *Chem Res Toxicol* 20(12):1778–1786, DOI 10.1021/tx700207u, PMID 18001057 [doi:10.1021/tx700207u; pmid:18001057]: shows the CYP3A4 paralog **CYP3A5 has Ser at the position homologous to CYP3A4 Cys239** ("S239") and exhibits only reversible (not time-dependent) inhibition by raloxifene; a CYP3A4 **C239A** point mutant likewise loses time-dependent inactivation by raloxifene, directly implicating residue 239 functionally.
  - VandenBrink BM et al. (2012), *Mol Pharmacol* 82(5):835–842, DOI 10.1124/mol.112.080739, PMID 22859722 [doi:10.1124/mol.112.080739; pmid:22859722]: proteolytic digests of recombinant **CYP3A4 Supersomes** after raloxifene time-dependent inactivation localized the adduct to **Cys239** (and, in CYP2C8, to the analogous Cys225), reproducing the Baer et al. finding in an independent experiment/paper from the same lab.

### 3. Structural context of Cys239

- UniProt's `Binding site` annotation for P08684 identifies **Cys442** as the "axial binding residue" for heme (ChEBI:CHEBI:30413) — i.e., **Cys442 is the heme-thiolate proximal ligand** of CYP3A4 [uniprot:P08684]. This resides within the canonical P450 heme-binding "Cys pocket" motif region near the C-terminus of the protein (UniProt sequence context around residue 442: `...FGSGPRNCIGMRFALMNMK...`, consistent with the FxxGxxxCxG P450 signature).
- Cys239 is sequence- and structure-distal from Cys442: per UniProt's PDB-derived secondary-structure features (evidence source PDB 5VCC [pdb:5VCC]), residue 239 falls in a non-helical loop connecting `Helix(230–235)` and `Helix(243–260)` [uniprot:P08684].
- A dedicated crystallographic/biochemical study of a cysteine-depleted CYP3A4 construct (Sevrioukova 2017, *Biochemistry* 56(24):3058–3067, DOI 10.1021/acs.biochem.7b00334, PMID 28590129, PMCID PMC5858725) explicitly states **"Cys239 is located in the G′–G connecting loop, part of the F–G [structural] fragment"** and identifies Cys239 as one of CYP3A4's six non-heme-ligating cysteines (Cys58, Cys64, Cys98, Cys239, Cys377, Cys468) [doi:10.1021/acs.biochem.7b00334; pmid:28590129]. In that study, engineering a C239T substitution (part of a Cys-depleted CYP3A4 variant used for crystallography) produced new hydrogen-bonding contacts to the Phe241 backbone amide and a recruited water molecule, and was associated with altered conformational flexibility of the B–C and F–G structural fragments — consistent with Cys239 sitting in a mobile loop region rather than a rigid core helix, though this study does not itself make a claim about "substrate access channel" membership, and no such specific claim is cited here beyond what the source states.
- The F–G region of mammalian microsomal P450s (including CYP3A4) is broadly known in the structural literature to participate in substrate access/egress; however, **no specific, retrievable source was verified in this research that explicitly states Cys239 lines the substrate access channel** (an earlier AI-generated web-search summary claimed this citing a review, PMC3787833, but on direct inspection of that review's text, no mention of "Cys239," "F-G loop," or residues 202–258 could be found — that specific claim is therefore **not verified** and is deliberately **excluded** from the sourced findings above rather than asserted).

### Other mechanism-based inactivators and Cys239

- All Cys239-apoprotein-adduct literature located in this research (Baer et al. 2007 [doi:10.1021/tx700037e]; Pearson et al. 2007 [doi:10.1021/tx700207u]; VandenBrink et al. 2012 [doi:10.1124/mol.112.080739]) concerns **raloxifene** specifically. No source was found in this research reporting a different (non-raloxifene) mechanism-based inactivator alkylating CYP3A4 apoprotein at Cys239; this absence-of-evidence is based on the searches performed (PubMed E-utilities queries for "Cys239 CYP3A4," "C239 AND CYP3A4," and related terms — see Methods) and should not be read as a proof that no such report exists elsewhere.

## Methods / sources consulted

- **UniProt**: retrieved P08684 in FASTA (`https://rest.uniprot.org/uniprotkb/P08684.fasta`) and full JSON record (`https://rest.uniprot.org/uniprotkb/P08684.json`) on 2026-09-24; confirmed accession, entry type (Swiss-Prot reviewed), sequence version (4), entry version (254), last annotation update (2026-09-02). Located `GFCMFDMECHK` in the canonical sequence programmatically (Python substring search); confirmed uniqueness (single occurrence) and 1-based coordinates. Extracted `Binding site` (heme, Cys442), `Initiator methionine`, `Chain`, `Transmembrane`, and `Helix` (PDB-derived secondary structure) features from the JSON record.
- **In-silico tryptic digest**: performed a standard full in-silico trypsin digest (cleave C-terminal to K or R, except when followed by P; zero missed cleavages) of the P08684 canonical sequence in Python to determine which tryptic peptide contains residue 239, and to confirm it is distinct from the peptide containing residues 56–66.
- **NCBI PubMed E-utilities** (`esearch`, `esummary`, `efetch`, `rettype=abstract`): used to confirm the exact citation (authors, title, journal, year, volume, issue, pages, DOI) and to retrieve the full abstract text for PMID 17497897 (Baer et al. 2007), PMID 18001057 (Pearson et al. 2007), PMID 22859722 (VandenBrink et al. 2012), and PMID 28590129 (Sevrioukova 2017).
- **WebFetch** of PMC5858725 (open-access full text of Sevrioukova 2017) to extract the specific structural statement about Cys239's location (G′–G loop / F-G fragment). A WebFetch attempt on the paywalled ACS full text of Baer et al. 2007 (`https://pubs.acs.org/doi/10.1021/tx700037e`) returned HTTP 403 (no access) — full-text-only details (exact proteinase K peptide sequence; recombinant vs. HLM enzyme source) could not be verified and are explicitly flagged as unverified above rather than asserted.
- **WebFetch** of PMC3787833 (a general CYP3A4 mechanism review) was attempted to check a web-search-generated claim about a "F-G loop, residues 202–258" location for Cys239; on direct text inspection, no such content was found in that article, so the claim was **discarded** as unverifiable rather than included.
- **RCSB PDB REST API** (`https://data.rcsb.org/rest/v1/core/entry/5VCC`) to confirm PDB entry 5VCC's title and initial release date (2017-05-31), the structure UniProt cites as evidence for the Helix(230–235)/Helix(243–260) secondary-structure annotation flanking residue 239.

## References
- Baer BR, Wienkers LC, Rock DA. Time-dependent inactivation of P450 3A4 by raloxifene: identification of Cys239 as the site of apoprotein alkylation. Chem Res Toxicol. 2007;20(6):954-964. doi:10.1021/tx700037e. PMID:17497897.
- Pearson JT, Wahlstrom JL, Dickmann LJ, Kumar S, Halpert JR, Wienkers LC, Foti RS, Rock DA. Differential time-dependent inactivation of P450 3A4 and P450 3A5 by raloxifene: a key role for C239 in quenching reactive intermediates. Chem Res Toxicol. 2007;20(12):1778-1786. doi:10.1021/tx700207u. PMID:18001057.
- VandenBrink BM, Davis JA, Pearson JT, Foti RS, Wienkers LC, Rock DA. Cytochrome P450 architecture and cysteine nucleophile placement impact raloxifene-mediated mechanism-based inactivation. Mol Pharmacol. 2012;82(5):835-842. doi:10.1124/mol.112.080739. PMID:22859722.
- Sevrioukova IF. High-Level Production and Properties of the Cysteine-Depleted Cytochrome P450 3A4. Biochemistry. 2017;56(24):3058-3067. doi:10.1021/acs.biochem.7b00334. PMID:28590129. PMCID:PMC5858725.
- UniProt Consortium. UniProtKB entry P08684 (CP3A4_HUMAN), sequence version 4, entry version 254, last annotation update 2026-09-02. https://www.uniprot.org/uniprotkb/P08684 (accessed 2026-09-24).
- RCSB Protein Data Bank entry 5VCC, "Crystal structure of human CYP3A4 bound to glycerol" (Sevrioukova IF), initial release 2017-05-31. https://www.rcsb.org/structure/5VCC (accessed 2026-09-24).
