---
# Research-finding document — external knowledge, structured like a finding but for the literature/tools.
# Spec: doc 04.4. Stored in research/<slug>.md. A research finding without verified references is NOT accepted.

topic: "Bioactivation chemistry of raloxifene to reactive metabolites and the resulting covalent adduct mass on protein cysteine/nucleophiles"
type: general                            # protein | gene | disease | pathway | publication | software | general
created: 2026-09-24
updated: 2026-09-24

status: reviewed                         # draft | reviewed  (reviewed = passed the research-reviewer)
reviewed_by: research-reviewer           # set by the librarian on the reviewer's ACCEPT (review provenance)
reviewed_date: 2026-09-24                # YYYY-MM-DD, stamped alongside reviewed_by

# Normalized domain-entity references this research is about (same scheme as findings; canonical IDs).
entities:
  - { db: pubchem, id: "CID5035", label: "raloxifene" }
  - { db: chebi, id: "CHEBI:8772", label: "raloxifene" }
  - { db: uniprot, id: "P08684", label: "CYP3A4 (human)" }

# (required, non-empty) Every external claim traces to a reference here. Each is fact-checked by the
# research-reviewer: it must EXIST and SUPPORT the claim attributed to it. Hallucinated citations are
# the known failure mode this section + review exist to stop.
references:
  - id: "url:https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/5035/property/MolecularFormula,MonoisotopicMass,IUPACName/JSON"
    type: url
    claim: "Raloxifene (PubChem CID 5035) has molecular formula C28H27NO4S and monoisotopic/exact mass 473.16607952 Da, confirming the element-wise calculation and the observed value of 473.1661 Da for the intact (unmodified) free base. The same record gives the IUPAC name [6-hydroxy-2-(4-hydroxyphenyl)-1-benzothiophen-3-yl]-[4-(2-piperidin-1-ylethoxy)phenyl]methanone, establishing that raloxifene's two free hydroxyls are the benzothiophene 6-OH and the 4-hydroxyphenyl ring attached at benzothiophene C2, while the ring bearing the piperidinylethoxy substituent is an ether with no free hydroxyl."
    verified: true
    verified_by: research-reviewer
  - id: "url:https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:8772"
    type: url
    claim: "Canonical ChEBI identifier for raloxifene is CHEBI:8772 (ChEBI name 'raloxifene')."
    verified: true
    verified_by: research-reviewer
  - id: "url:https://www.uniprot.org/uniprotkb/P08684/entry"
    type: url
    claim: "CYP3A4 (human) canonical UniProt accession is P08684 (entry name CP3A4_HUMAN), gene name CYP3A4, organism Homo sapiens. The entry's only mention of raloxifene is a DrugBank cross-reference (DB00481); it does not itself state a hepatic-abundance figure or a bioactivation/mechanism-based-inactivation mechanism for raloxifene (those claims are sourced to Chen2002/Baer2007 below, not to this UniProt entry)."
    verified: true
    verified_by: research-reviewer
  - id: "doi:10.1021/tx0200109"
    type: doi
    claim: "Chen Q, Ngui JS, Doss GA, Wang RW, Cai X, DiNinno FP, Blizzard TA, Hammond ML, Stearns RA, Evans DC, Baillie TA, Tang W. 'Cytochrome P450 3A4-Mediated Bioactivation of Raloxifene: Irreversible Enzyme Inhibition and Thiol Adduct Formation.' Chem Res Toxicol. 2002 Jul;15(7):907-14. PMID 12119000. Raloxifene irreversibly (NADPH- and preincubation-time-dependent) inhibits CYP3A4 in human liver microsomes (K_I = 9.9 microM, k_inact = 0.16 min^-1); loss of activity is partially attenuated by glutathione, implicating a reactive metabolite. GSH adducts of raloxifene form by substitution at the 5- or 7-position of the benzothiophene ring or the 3'-position of the phenol ring (7-glutathionyl derivative most abundant), postulated to arise either from GSH addition to a raloxifene arene oxide followed by dehydration/aromatization, or from GSH trapping of an extended quinone intermediate. GSH-adduct formation was almost abolished by ketoconazole or anti-CYP3A4 IgG, implicating CYP3A4 as the responsible enzyme. GSH adducts were also detected in rat/human hepatocytes, and the corresponding N-acetylcysteine (mercapturate) conjugates were found in rat bile and urine after oral dosing."
    verified: true
    verified_by: research-reviewer
  - id: "doi:10.1021/tx0342722"
    type: doi
    claim: "Yu L, Liu H, Li W, Zhang F, Luckie C, van Breemen RB, Thatcher GRJ, Bolton JL. 'Oxidation of Raloxifene to Quinoids: Potential Toxic Pathways via a Diquinone Methide and o-Quinones.' Chem Res Toxicol. 2004 Jul;17(7):879-88. PMID 15257612. Raloxifene incubated with GSH and rat/human liver microsomes is converted to raloxifene diquinone-methide GSH conjugates, raloxifene o-quinone GSH conjugates, and raloxifene catechols (7-hydroxyraloxifene oxidizes further to the 6,7-o-quinone). Raloxifene diquinone methide has a half-life <1 s in phosphate buffer, whereas raloxifene 6,7-o-quinone is far more stable (t1/2 = 69 +/- 2.5 min), implying the o-quinone may be the more persistent/toxic species; 7-hydroxyraloxifene was more cytotoxic than raloxifene in S30 and MDA-MB-231 human breast cancer cell lines."
    verified: true
    verified_by: research-reviewer
  - id: "doi:10.1021/tx700037e"
    type: doi
    claim: "Baer BR, Wienkers LC, Rock DA. 'Time-Dependent Inactivation of P450 3A4 by Raloxifene: Identification of Cys239 as the Site of Apoprotein Alkylation.' Chem Res Toxicol. 2007 Jun;20(6):954-64. PMID 17497897. Mass spectrometry of intact CYP3A4 protein and of proteinase-K-generated peptides showed a mass shift of 471 Da on both the intact apoprotein and the modified peptide, localized to residue Cys239, consistent with covalent attack of the Cys239 thiolate on the raloxifene diquinone-methide intermediate (not heme adduction). Pre-alkylation of Cys239 with iodoacetamide or N-(1-pyrene)iodoacetamide, which modified Cys239 exclusively, prevented raloxifene time-dependent inactivation of CYP3A4, providing functional confirmation that Cys239 apoprotein adduction (mass change +471 Da, i.e., raloxifene minus 2 H) accounts for the mechanism-based inactivation."
    verified: true
    verified_by: research-reviewer
  - id: "doi:10.1021/bi902213r"
    type: doi
    claim: "Moore CD, Reilly CA, Yost GS. 'CYP3A4-Mediated Oxygenation versus Dehydrogenation of Raloxifene.' Biochemistry. 2010 Jun 1;49(21):4466-75. PMID 20405834; PMCID PMC2887288. Using 18O-incorporation experiments, the authors showed 3'-hydroxyraloxifene is produced exclusively via CYP3A4-mediated oxygenation (the arene-oxide/hydroxylation pathway), whereas CYP3A4-mediated dehydrogenation of raloxifene to the reactive diquinone methide (not the arene-oxide/oxygenation pathway) is the mechanism responsible for the protein-binding, enzyme-inactivating intermediate. 7-Hydroxyraloxifene, previously assumed to be an ordinary oxygenation product, instead arises via hydrolysis of a putative ester formed by conjugation of the raloxifene diquinone methide with a carboxylic acid moiety on CYP3A4 or other proteins."
    verified: true
    verified_by: research-reviewer
  - id: "doi:10.1021/bi101139q"
    type: doi
    claim: "Moore CD, Shahrokh K, Sontum SF, Cheatham TE 3rd, Yost GS. 'Improved Cytochrome P450 3A4 Molecular Models Accurately Predict the Phe215 Requirement for Raloxifene Dehydrogenation Selectivity.' Biochemistry. 2010 Oct 19;49(41):9011-9. PMID 20812728; PMCID PMC2958526. Molecular modeling (CYP3A4 crystal structure 1W0E) plus site-directed mutagenesis identified active-site residue Phe215 as required for orienting raloxifene for CYP3A4-mediated dehydrogenation to the electrophilic diquinone-methide intermediate that has been linked to CYP3A4 inactivation; Phe215-to-Gly/Gln substitutions decreased dehydrogenation rates without changing oxygenation rates, supporting dehydrogenation and oxygenation as mechanistically distinct CYP3A4 reaction pathways on raloxifene."
    verified: true
    verified_by: research-reviewer
  - id: "doi:10.1021/tx7001367"
    type: doi
    claim: "Liu H, Qin Z, Thatcher GRJ, Bolton JL. 'Uterine Peroxidase-Catalyzed Formation of Diquinone Methides from the Selective Estrogen Receptor Modulators Raloxifene and Desmethylated Arzoxifene.' Chem Res Toxicol. 2007 Nov;20(11):1676-84. PMID 17630709; PMCID PMC2507766. Rat uterine microsomes oxidize raloxifene to an electrophilic diquinone methide (trapped as GSH conjugates) via an NADPH-independent, cyanide-inhibitable, H2O2-enhanced (i.e., peroxidase-catalyzed, not CYP-catalyzed) pathway; a raloxifene Cys-Gly conjugate was identified, formed by gamma-glutamyl-transpeptidase hydrolysis of 7-glutathionyl raloxifene. Horseradish peroxidase directly oxidized raloxifene to phenoxyl-radical-derived dimers without added H2O2. Biotinylated raloxifene 'COATag' probes covalently modified several uterine-microsomal proteins, detected by Western blot, with modification enhanced by H2O2 and decreased by NADPH — i.e., a peroxidase-mediated, extrahepatic bioactivation route distinct from the hepatic CYP3A4 pathway."
    verified: true
    verified_by: research-reviewer
---

# Research: Bioactivation chemistry of raloxifene to reactive metabolites and the resulting covalent adduct mass on protein cysteine/nucleophiles

## Summary

Raloxifene free base (C28H27NO4S, PubChem CID 5035, ChEBI:8772) has a monoisotopic mass of 473.16608 Da, confirmed independently by element-wise summation and by PubChem's computed exact mass (473.16607952 Da) — matching the user-supplied value of 473.1661 Da [pubchem]. CYP3A4 (UniProt P08684) bioactivates raloxifene by two distinct oxidative reaction types: (1) dehydrogenation (formal loss of 2 H) to an electrophilic diquinone-methide, and (2) oxygenation, which is established (via 18O-labeling) to produce 3'-hydroxyraloxifene [Moore2010a]. An arene-oxide chemical intermediate has been postulated by one study to explain certain GSH adducts [Chen2002], but is explicitly excluded by another as the source of the CYP3A4-inactivating species [Moore2010a]; 7-hydroxyraloxifene (which further oxidizes to the 6,7-o-quinone) arises not from oxygenation but from hydrolysis of a putative ester formed between the diquinone-methide and a carboxylic acid moiety of CYP3A4 or other proteins in the reconstituted system [Moore2010a][Yu2004]. For the specific, directly measured case of covalent modification of CYP3A4's own apoprotein at Cys239 — the mechanistically best-characterized raloxifene-protein-thiol adduct in the literature — mass spectrometry of the intact protein and of the modified peptide showed a mass shift of **471 Da**, matching **M − 2H (471.1504 Da)**, i.e., addition of the diquinone-methide-derived species to the cysteine thiol, not the +489 Da (M + O, arene-oxide) or +473 Da (M + 0, no covalent change) alternatives [Baer2007]. 18O-labeling experiments confirmed that dehydrogenation to the diquinone methide (not the oxygenation route, which instead yields the stable 3'-hydroxyraloxifene metabolite) is the pathway responsible for the enzyme-inactivating, protein-adducting species, and explicitly excluded an arene-oxide intermediate as the source of that species [Moore2010a]. The arene-oxide route is, at best, postulated (via GSH trapping in liver microsomes/hepatocytes) as giving adducts at the 5-/7-positions of the benzothiophene ring and the 3'-position of the phenol ring, formed via addition-dehydration-aromatization rather than as a simple, un-rearranged epoxide-opening adduct — with an extended-quinone mechanism offered by the same source as an explicit alternative [Chen2002]. Raloxifene is a well-established mechanism-based (time-dependent, NADPH- and preincubation-time-dependent, irreversible) inactivator of CYP3A4, and this inactivation is attributed to apoprotein (Cys239) adduction rather than heme adduction, based on direct intact-protein/peptide mass spectrometry and on protection from inactivation by pre-alkylating Cys239 [Chen2002][Baer2007]. A separate, non-CYP (uterine peroxidase) route to the same diquinone-methide chemistry has also been investigated in extrahepatic (uterine) tissue [Liu2007].

## Detailed findings

### 1. Monoisotopic mass of raloxifene and the load-bearing delta-mass arithmetic

Raloxifene free base: **C28H27NO4S**, PubChem CID 5035, ChEBI:8772 [pubchem][chebi].

Element-wise monoisotopic mass calculation (standard monoisotopic atomic masses: ¹²C = 12.0000000, ¹H = 1.0078250319, ¹⁴N = 14.0030740052, ¹⁶O = 15.9949146221, ³²S = 31.97207069):

| Element | Count | Mass each (Da) | Subtotal (Da) |
|---|---|---|---|
| C | 28 | 12.0000000 | 336.0000000 |
| H | 27 | 1.0078250319 | 27.2112759 |
| N | 1 | 14.0030740052 | 14.0030740 |
| O | 4 | 15.9949146221 | 63.9796585 |
| S | 1 | 31.9720706900 | 31.9720707 |
| **Total (M)** | | | **473.1660790 Da** |

This matches PubChem's computed monoisotopic/exact mass for CID 5035 (473.16607952 Da) [pubchem] to 5 decimal places (both round to 473.16608 Da), and matches the value the user supplied (473.1661 Da).

Two chemically distinct delta masses follow directly from arithmetic on M, corresponding to two distinct bioactivation chemistries:

- **M − 2H = 473.1660790 − 2×1.0078250 = 471.1504290 Da.** This is the delta expected for a **quinone-methide-type Michael addition**: CYP-mediated dehydrogenation of raloxifene's phenol(s) forms an electrophilic (di)quinone-methide (formal loss of 2 H relative to raloxifene); a thiol nucleophile (Cys-SH, GSH-SH) then adds across the exocyclic methylene in a conjugate (Michael) addition, and the ring re-aromatizes/re-phenolizes by internal proton transfer from the thiol's own SH to the developing phenolate oxygen — no atoms are lost to solvent, so the net mass added to the protein/peptide equals the mass of the diquinone-methide species itself, i.e., M − 2H.
- **M + O = 473.1660790 + 15.9949146 = 489.1609937 Da.** This is the delta expected for a simple, un-rearranged **arene-oxide/epoxide ring-opening** adduct: CYP-mediated monooxygenation inserts one oxygen atom to form an arene oxide, and if a thiol nucleophile opens the epoxide directly (SN2-type addition across the C–O bond, retaining the oxygen as a new hydroxyl) without subsequent dehydration, the adduct's added mass equals M + O.
- **M + 0 = 473.1660790 Da** — no mass change — is not a covalent-adduct chemistry at all for either mechanism above; both quinone-methide and epoxide chemistries necessarily change the mass (by −2H or, before any dehydration, by +O). An unshifted mass would correspond only to non-covalent, reversible binding, which is not the mechanism reported for the reactive metabolites discussed here.

**Which is reported in the literature for the raloxifene-protein-cysteine adduct:** direct mass spectrometry of the intact CYP3A4 apoprotein, and of the modified proteolytic peptide, showed a **mass shift of 471 Da** relative to unmodified control protein, localized to Cys239 — i.e., **M − 2H (471.1504 Da), not M + O (489.16 Da)** [Baer2007]. This was corroborated mechanistically by 18O-incorporation experiments showing that CYP3A4-mediated **dehydrogenation** (not oxygenation/arene-oxide formation) produces the diquinone-methide responsible for enzyme inactivation, while the oxygenation/arene-oxide route instead yields the isolable, non-inactivating metabolite 3'-hydroxyraloxifene [Moore2010a]. I did not find any primary study reporting a directly measured +489 Da (M + O) mass shift on a protein or peptide for raloxifene; the arene-oxide route's characterized GSH conjugates are reported as arising via epoxide addition **followed by dehydration and aromatization** [Chen2002] — i.e., the isolated/characterized arene-oxide-derived adducts are not simple, un-rearranged M + O epoxide adducts either (their net mass addition relative to raloxifene, after dehydration, is GSH mass minus 2H, not GSH mass plus O — this composition is for the small-molecule GSH conjugate, not for a direct protein-Cys delta mass, and I did not find it characterized as such on a protein/peptide in the literature I verified). **This is an unverified gap**: no paper I located reports a +489 Da (M + O) covalent adduct mass measured directly on a protein/peptide for raloxifene; only the +471 Da (M − 2H) value is directly measured on protein.

### 2. Bioactivation pathways and the CYPs mediating them

CYP3A4 acts on raloxifene via two distinct, established oxidative reaction types:

- **(a) Dehydrogenation to a (di)quinone methide.** Raloxifene has two free phenolic hydroxyls: the benzothiophene 6-OH and the 4-hydroxyphenyl ring attached at benzothiophene C2 (the third aromatic ring — the 4-(2-piperidin-1-ylethoxy)phenyl/benzoyl ring — bears the piperidinylethoxy ether side chain and has no free hydroxyl) [pubchem][chebi]. CYP3A4-mediated dehydrogenation across these two phenols, conjugated through the benzothiophene C2-aryl linkage, generates an electrophilic quinone-methide at each ring; because raloxifene can be doubly dehydrogenated in this way the reactive species is termed a "diquinone methide" [Yu2004][Chen2002]. This diquinone methide is extremely short-lived (t1/2 < 1 s in phosphate buffer) and is trapped efficiently by GSH to give diquinone-methide GSH conjugates [Yu2004]. The dehydrogenation route (rather than oxygenation) was directly implicated as the CYP3A4-inactivating pathway by 18O-labeling, which also explicitly excluded an arene-oxide intermediate as the source of the inactivating species [Moore2010a], and active-site residue Phe215 was shown by molecular modeling and site-directed mutagenesis to orient raloxifene for this dehydrogenation selectively (Phe215Gly/Gln mutants reduced dehydrogenation without affecting oxygenation) [Moore2010b].
- **(b) Oxygenation to 3'-hydroxyraloxifene.** CYP3A4-mediated oxygenation (monooxygenation) is established, via 18O-incorporation, to produce 3'-hydroxyraloxifene, confirmed to derive specifically from the oxygenation (not dehydrogenation) pathway [Moore2010a]. An arene-oxide chemical intermediate is sometimes invoked to explain this oxygenation and separately reported GSH adducts at the 5- or 7-position of the benzothiophene ring or the 3'-position of the phenol ring; Chen2002 states these adducts are "postulated to derive from addition of GSH to raloxifene arene oxides" (with subsequent dehydration/aromatization), while also proposing, as an explicit alternative, GSH trapping of an extended quinone intermediate [Chen2002]. Moore2010a's 18O-labeling data explicitly exclude the arene-oxide pathway as the source of the CYP3A4-inactivating, protein-adducting intermediate — only the dehydrogenation/diquinone-methide route accounts for that species [Moore2010a]. 7-Hydroxyraloxifene, previously assumed to be a standard oxygenation product, was instead shown to arise via hydrolysis of a putative ester formed between the diquinone methide and a carboxylic acid moiety on CYP3A4 or other proteins — i.e., it is a downstream product of the dehydrogenation pathway, not of oxygenation [Moore2010a]. Separately, 7-hydroxyraloxifene itself oxidizes further to the more stable 6,7-o-quinone (t1/2 = 69 ± 2.5 min), a redox-active quinoid that is trapped by GSH as o-quinone conjugates and that may be more persistent (and thus potentially more toxic) than the diquinone methide itself [Yu2004].

**CYP mediating bioactivation:** CYP3A4 is the enzyme directly implicated by chemical-inhibitor (ketoconazole) and antibody (anti-CYP3A4 IgG) experiments, which almost abolished GSH-adduct formation in human liver microsomes [Chen2002], and by all of the mechanistic dehydrogenation/oxygenation work reviewed above, which used purified/reconstituted CYP3A4 [Moore2010a][Moore2010b][Baer2007]. I did not find, in the sources verified for this finding, primary literature establishing other CYP isoforms (e.g., CYP2C8, CYP2D6) as mediators of raloxifene bioactivation to these reactive quinoids — **that claim is unverified** and is not asserted here. A distinct, non-CYP bioactivation pathway was also investigated: oxidation of raloxifene to the same diquinone-methide chemistry in rat uterine microsomes was NADPH-independent, cyanide-inhibitable, and H2O2-enhanced, suggesting this extrahepatic bioactivation could be mediated by uterine peroxidases rather than by P450s [Liu2007].

### 3. Verification of the specific primary-literature citations named by the requester

- **"Chen Q et al. 2002 (Chem Res Toxicol, raloxifene bioactivation / GSH conjugates)"** — verified. Chen Q, Ngui JS, Doss GA, Wang RW, Cai X, DiNinno FP, Blizzard TA, Hammond ML, Stearns RA, Evans DC, Baillie TA, Tang W. "Cytochrome P450 3A4-Mediated Bioactivation of Raloxifene: Irreversible Enzyme Inhibition and Thiol Adduct Formation." Chem Res Toxicol. 2002 Jul;15(7):907-14. doi:10.1021/tx0200109, PMID 12119000 [Chen2002]. Confirmed via PubMed/NCBI E-utilities full abstract retrieval; content matches the description above (irreversible, NADPH/time-dependent CYP3A4 inhibition; GSH adducts at benzothiophene 5-/7- and phenol 3'-positions; CYP3A4 implicated by ketoconazole/anti-CYP3A4-IgG blockade; hepatocyte/bile/urine adducts).
- **"work by Yu L / Yu et al."** — verified as Yu L, Liu H, Li W, Zhang F, Luckie C, van Breemen RB, Thatcher GRJ, Bolton JL. "Oxidation of Raloxifene to Quinoids: Potential Toxic Pathways via a Diquinone Methide and o-Quinones." Chem Res Toxicol. 2004 Jul;17(7):879-88. doi:10.1021/tx0342722, PMID 15257612 [Yu2004]. Confirmed via PubMed/NCBI E-utilities; content matches (diquinone methide and o-quinone GSH conjugates; half-lives; cytotoxicity of 7-hydroxyraloxifene). Note: this is the paper that reports the diquinone-methide/o-quinone chemistry and half-lives; if a different Yu L paper was intended, it was not located and is not claimed here.
- **"Moore CD et al. 2010"** — verified, and in fact there are **two** 2010 Moore CD papers on raloxifene/CYP3A4, both confirmed via PubMed/NCBI E-utilities:
  1. Moore CD, Reilly CA, Yost GS. "CYP3A4-Mediated Oxygenation versus Dehydrogenation of Raloxifene." Biochemistry. 2010 Jun 1;49(21):4466-75. doi:10.1021/bi902213r, PMID 20405834, PMCID PMC2887288 [Moore2010a]. This is the mechanistically load-bearing one for the M−2H vs M+O question: 18O-labeling distinguishes the oxygenation (→3'-hydroxyraloxifene) pathway from the dehydrogenation (→diquinone methide) pathway, explicitly excludes an arene-oxide intermediate as the source of the CYP3A4-inactivating species, and assigns dehydrogenation as responsible for CYP3A4 inactivation.
  2. Moore CD, Shahrokh K, Sontum SF, Cheatham TE 3rd, Yost GS. "Improved Cytochrome P450 3A4 Molecular Models Accurately Predict the Phe215 Requirement for Raloxifene Dehydrogenation Selectivity." Biochemistry. 2010 Oct 19;49(41):9011-9. doi:10.1021/bi101139q, PMID 20812728, PMCID PMC2958526 [Moore2010b]. Structural/mutagenesis follow-up identifying Phe215 as controlling dehydrogenation selectivity.
  Both are legitimate, verified 2010 Moore CD papers on this exact topic; neither reports content inconsistent with what is summarized above.
- Additionally verified (not explicitly named by the requester but load-bearing for the delta-mass question): Baer BR, Wienkers LC, Rock DA. "Time-Dependent Inactivation of P450 3A4 by Raloxifene: Identification of Cys239 as the Site of Apoprotein Alkylation." Chem Res Toxicol. 2007 Jun;20(6):954-64. doi:10.1021/tx700037e, PMID 17497897 [Baer2007]. This is the paper that directly measures the 471 Da (M − 2H) mass shift on CYP3A4 protein/peptide and localizes it to Cys239 — it is the primary evidentiary source for the delta-mass answer in Finding 1 and should be treated as at least as important as the three requester-named papers for that specific question.

No misattributions were found: all four requester-referenced works (Chen 2002, Yu et al., and both Moore CD 2010 papers) exist as described and report the content attributed to them.

### 4. Glutathione/cysteine adducts: positions and stability

Characterized GSH conjugates of raloxifene include:
- **7-glutathionyl raloxifene** (most abundant GSH adduct in human liver microsomes) and **5-glutathionyl raloxifene**, both substitutions on the benzothiophene ring, plus a **3'-glutathionyl** conjugate on the phenol ring of the 4-hydroxyphenyl ring attached at benzothiophene C2 [Chen2002]. These are postulated to derive from GSH addition to a raloxifene arene oxide followed by dehydration/aromatization, or alternatively from GSH trapping of an extended quinone intermediate [Chen2002].
- **Raloxifene diquinone-methide GSH conjugates** and **raloxifene o-quinone (6,7-o-quinone) GSH conjugates**, characterized by LC-MS/MS and (for synthesized catechol standards) NMR [Yu2004].
- **Raloxifene Cys-Gly**, a dipeptide conjugate formed by gamma-glutamyl-transpeptidase-mediated hydrolysis of 7-glutathionyl raloxifene, characterized in rat uterine microsomal incubations [Liu2007].
- Downstream **N-acetylcysteine (mercapturic acid) conjugates**, the expected further metabolic processing (via the mercapturic acid pathway) of the initial GSH adducts, detected in rat bile and urine after oral raloxifene dosing [Chen2002].

The site of thiol attack is thus the electrophilic carbon centers generated on the benzothiophene ring (positions 5 and 7, flanking the ring's phenolic oxygen) and on the phenolic ring of the 4-hydroxyphenyl ring attached at benzothiophene C2 (3'-position) — i.e., both of raloxifene's two free-hydroxyl-bearing aromatic rings are sites of quinone-methide formation (and, per Chen2002's postulated mechanism, arene-oxide formation) and subsequent thiol conjugation [Chen2002][Yu2004].

**Stability/reversibility:** the GSH/Cys-Gly/NAC conjugates recovered from microsomes, hepatocytes, bile, and urine are stable, isolable, characterizable species (by LC-MS/MS and NMR), consistent with irreversible (covalent) thioether bond formation, not a reversible/dissociable adduct [Chen2002][Yu2004][Liu2007]. For the specific protein-cysteine case (CYP3A4 Cys239), the adduct survived intact-protein and proteolytic-peptide sample preparation for mass spectrometry and was associated with permanent (mechanism-based) loss of enzyme activity that was prevented by blocking Cys239 in advance (the one intervention tested; the study does not establish this as the only possible way to prevent inactivation) — direct evidence the protein-Cys adduct is stable/irreversible under the conditions studied [Baer2007].

### 5. Raloxifene as a mechanism-based (time-dependent) CYP3A4 inactivator; apoprotein vs. heme adduction

Raloxifene is an established mechanism-based (time-dependent) inactivator of CYP3A4: inhibition in human liver microsomes was irreversible, and required both NADPH and a preincubation period, with kinetic constants K_I = 9.9 microM and k_inact = 0.16 min^-1 [Chen2002]. The inactivation is attributed to **apoprotein adduction, not heme adduction**: direct mass spectrometry of intact CYP3A4 protein and of a proteinase-K-generated peptide showed a single equivalent of raloxifene-derived mass (471 Da) bound to the apoprotein, localized specifically to residue **Cys239**; pre-blocking Cys239 with cysteine-reactive alkylating reagents (iodoacetamide, N-(1-pyrene)iodoacetamide) prevented the time-dependent inactivation, directly linking Cys239 apoprotein alkylation (rather than any heme modification) to loss of catalytic activity [Baer2007]. The chemical identity of the adducting species — the raloxifene diquinone methide reacting via Michael addition with the Cys239 thiolate — was established by the 471 Da (M − 2H) mass shift and corroborated by 18O-exclusion experiments showing the inactivating pathway is dehydrogenation, not oxygenation [Baer2007][Moore2010a].

## Methods / sources consulted

- **PubChem PUG REST API** (`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/5035/property/...`) queried directly for CID 5035 (raloxifene): confirmed molecular formula C28H27NO4S and monoisotopic/exact mass 473.16607952 Da.
- **ChEBI** (EBI) entry lookup confirmed CHEBI:8772 as the canonical raloxifene identifier.
- **UniProtKB** entry P08684 confirmed as the canonical human CYP3A4 accession.
- **NCBI PubMed via E-utilities `efetch`** (`db=pubmed`, `rettype=abstract`) used to retrieve full, verbatim abstract text (not just search-snippet paraphrase) for all cited primary papers, to confirm exact author lists, titles, journal, year/volume/pages, DOI, and PMID, and to confirm that each paper's actual reported content matches the claims attributed to it: PMID 12119000 (Chen Q et al. 2002), PMID 15257612 (Yu L et al. 2004), PMID 17497897 (Baer BR et al. 2007), PMID 20405834 (Moore CD et al. 2010, Biochemistry 49:4466), PMID 20812728 (Moore CD et al. 2010, Biochemistry 49:9011), PMID 17630709 (Liu H et al. 2007).
- **NCBI PubMed via E-utilities `esearch`** used to confirm there are exactly two 2010 "Moore CD" + "raloxifene" papers in PubMed (disambiguating the requester's "Moore CD et al. 2010" reference), and to locate the uterine-peroxidase paper PMID.
- **Manual element-wise monoisotopic mass calculation** for C28H27NO4S performed and cross-checked against PubChem's computed exact mass (agreement to the fifth decimal place; both round to 473.16608 Da).
- Web search (general) used only to locate candidate PMIDs/DOIs for follow-up verification via the primary-source methods above; no claim in this document rests on a search-engine AI summary alone — every claim was checked against the retrieved primary abstract text.

## References

- Chen Q, Ngui JS, Doss GA, Wang RW, Cai X, DiNinno FP, Blizzard TA, Hammond ML, Stearns RA, Evans DC, Baillie TA, Tang W. "Cytochrome P450 3A4-Mediated Bioactivation of Raloxifene: Irreversible Enzyme Inhibition and Thiol Adduct Formation." *Chem Res Toxicol.* 2002 Jul;15(7):907-14. doi:10.1021/tx0200109. PMID: 12119000. [Chen2002]
- Yu L, Liu H, Li W, Zhang F, Luckie C, van Breemen RB, Thatcher GRJ, Bolton JL. "Oxidation of Raloxifene to Quinoids: Potential Toxic Pathways via a Diquinone Methide and o-Quinones." *Chem Res Toxicol.* 2004 Jul;17(7):879-88. doi:10.1021/tx0342722. PMID: 15257612. [Yu2004]
- Baer BR, Wienkers LC, Rock DA. "Time-Dependent Inactivation of P450 3A4 by Raloxifene: Identification of Cys239 as the Site of Apoprotein Alkylation." *Chem Res Toxicol.* 2007 Jun;20(6):954-64. doi:10.1021/tx700037e. PMID: 17497897. [Baer2007]
- Moore CD, Reilly CA, Yost GS. "CYP3A4-Mediated Oxygenation versus Dehydrogenation of Raloxifene." *Biochemistry.* 2010 Jun 1;49(21):4466-75. doi:10.1021/bi902213r. PMID: 20405834. PMCID: PMC2887288. [Moore2010a]
- Moore CD, Shahrokh K, Sontum SF, Cheatham TE 3rd, Yost GS. "Improved Cytochrome P450 3A4 Molecular Models Accurately Predict the Phe215 Requirement for Raloxifene Dehydrogenation Selectivity." *Biochemistry.* 2010 Oct 19;49(41):9011-9. doi:10.1021/bi101139q. PMID: 20812728. PMCID: PMC2958526. [Moore2010b]
- Liu H, Qin Z, Thatcher GRJ, Bolton JL. "Uterine Peroxidase-Catalyzed Formation of Diquinone Methides from the Selective Estrogen Receptor Modulators Raloxifene and Desmethylated Arzoxifene." *Chem Res Toxicol.* 2007 Nov;20(11):1676-84. doi:10.1021/tx7001367. PMID: 17630709. PMCID: PMC2507766. [Liu2007]
- PubChem Compound Summary for CID 5035, Raloxifene. National Center for Biotechnology Information, PubChem PUG REST API (`MolecularFormula`, `MonoisotopicMass`, `IUPACName` properties retrieved directly). [pubchem]
- ChEBI entry CHEBI:8772, raloxifene. European Bioinformatics Institute. [chebi]
- UniProtKB entry P08684 (CP3A4_HUMAN), Cytochrome P450 3A4, *Homo sapiens*. [uniprot]
