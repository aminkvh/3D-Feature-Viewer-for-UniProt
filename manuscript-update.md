# 3D Feature Viewer for UniProt — manuscript update

All **bold** text marks an addition or edit relative to the tentative (v1.5.1) manuscript, made after comparing
the text against the current code base. Unchanged text is reproduced in plain type for context; sections not
shown (3.1–3.6 prioritization algorithms; Discussion; Funding; Conflict of Interest) were reviewed and are
unchanged in substance.

---

## Abstract

**Summary:** 3D Feature Viewer for UniProt is a browser-native extension that converts any UniProt entry into
an interactive 3D workspace. It integrates functional annotations, including post-translational modifications,
disease-associated variants, AlphaMissense scores, ProtVar predictions, and other residue-level features, with
available protein structures. **For each residue, a per-substitution prediction table reports AlphaMissense
together with EVE, ESM-1b, CADD, and FoldX ΔΔG, plus per-residue evolutionary conservation and a
Missense3D structural-effect call.** The tool allows users to examine these features in a structural context,
assess their spatial co-occurrence, and explore potential functional crosstalk. **In heteromeric complexes,
residues of partner chains — distinct UniProt entries resolved through SIFTS — are identified, and their disease
variants, PTMs, and sites can be overlaid and inspected together with the protein of interest.** Binding pocket
modules provide a rapid view of pocket residues, associated ligands, **pocket descriptors (predicted-pocket
count, buriedness, druggability, mean pocket confidence, and residue-class composition),** and proteins with
similar pockets to the protein of interest **(through an RCSB Structure Motif search seeded from the pocket
residues, supporting both experimental PDB entries and AlphaFold models)**. **Ligands carry chemical
descriptors with links to PubChem 2D and 3D similarity searches and an in-structure fingerprint comparison.**
Built-in prioritization routines support fast residue-level triage, while export options enable seamless
downstream analysis.

**Availability and Implementation:** The extension is implemented in JavaScript as a browser extension
compatible with Chrome and Firefox. It has no server-side component and retrieves data from public resources,
including UniProt, PDBe/SIFTS, PDB/RCSB **(including the RCSB Structure Motif Search service)**, AlphaFold
DB, 3D-Beacons, AlphaFill, AlphaMissense, ProtVar, PubChem, ChEMBL, and EBI variation resources where
available. The current viewer uses **Mol\* (Sehnal et al. 2021)** inside a sandboxed iframe **(migrated from the
3Dmol.js renderer used in earlier versions)**. Source code is freely available at
https://github.com/aminkvh/3D-Feature-Viewer-for-UniProt. Archived version: https://zenodo.org/records/20574851.

---

## 2 System overview and user workflow

(Paragraph 1 — entry points and structure retrieval — unchanged.)

(Paragraph 2 — annotation mapping. Edited:)

The workspace maps UniProt annotations onto the selected structure, including PTMs, disease-associated
variants with ClinVar (Landrum et al. 2018) and gnomAD (Gudmundsson et al. 2021) metadata, functional sites,
domain and region features, membrane topology, and predictive scores such as AlphaMissense and FoldX ΔΔG
(Valanciute et al. 2023). Additional residue-level scores and metadata are provided through ProtVar and related
resources (Frazer et al. 2021; Schubach et al. 2024; Rives et al. 2021; Hanna et al. 2024). **For a selected
residue, these are consolidated into a per-substitution prediction table that places AlphaMissense alongside the
ProtVar-served EVE, ESM-1b, CADD, and FoldX ΔΔG scores, together with a per-residue conservation value and
a Missense3D structural-effect call; substitutions that correspond to observed variants are marked and colored
by their clinical consequence. Because CADD is a genomic score, it is shown for the substitutions reachable by a
single nucleotide change.** These features can be displayed individually or layered together. Because all
annotations are shown in the same structural coordinate frame, the workspace makes spatial relationships central
to interpretation. When a residue is selected, the view centers on its local environment and reports nearby
annotations. This allows a variant, nearby PTMs, functional sites, domains, ligands, and other relevant features
to be examined together. **In multi-chain complexes, neighbouring residues that belong to partner chains are
colored by the protein of interest's disease scheme and listed on a separate line in the nearby-residue report;
selecting a partner-chain residue opens a minimal panel that links to that protein's own UniProt entry, and an
"other chains" menu offers the same links.** The same view supports one of the main goals of the tool: identifying
where disease-associated variation intersects with other functional features in three-dimensional space.

(Paragraph 3 — ligands and pockets. Edited:)

The workspace also integrates what the structure itself contains. Ligands and cofactors present in the model are
shown with basic chemical information **(molecular formula, molecular weight, hydrogen-bond donor and
acceptor counts, SMILES, and InChIKey, with external links to PubChem and DrugBank)**, can be compared by
CACTVS fingerprint (Kim et al. 2023) **and through PubChem 2D and 3D structure-similarity searches**, and are
viewed in their local structural context alongside the surrounding annotations. **For a ligand-binding site, a
pocket-evidence panel reports the closest protein–ligand contact, the mean pocket confidence (pLDDT for
predicted models, crystallographic B-factor for experimental structures), overlap with annotated functional
sites, and the residue-class composition of the pocket (hydrophobic, aromatic, acidic, basic, and polar), shown
as a compact radar summary. A structure-based similar-pocket search submits the pocket residues to the RCSB
Structure Motif Search service to return entries with a comparable three-dimensional residue arrangement,
working from experimental PDB entries and from AlphaFold computed-structure-model identifiers.** Beyond
direct annotation display, exploratory overlays highlight residues that may merit attention, turning the same
integrated view into a means of triage.

(New paragraph — multichain partner display and chimeric flag:)

**Where a loaded structure is a complex, partner chains — chains that map to a different UniProt entry than the
protein of interest — are detected through SIFTS and handled distinctly from the target protein. Their disease
variants, PTMs, and functional sites are fetched on demand and can be toggled on as spheres for minimal
visualization, while the partner residues themselves are not treated as the protein of interest (consistent with
the prioritization routines, which use partner residues only as contextual signal). The structure selector also
flags constructs in which the protein's own chain is a cross-species or multi-entry fusion (a chimera), since in
those cases the UniProt mapping covers only part of the chain and neighbouring residues may fall in unmapped
regions.**

(Paragraph 4 — export — unchanged.)

---

## 3 Residue-prioritization algorithms

Unchanged. (Sub-sections 3.1 PTM–variant proximity, 3.2 pathogenic-variant enrichment hotspots, 3.3
long-range contact hubs, 3.4 mutation/phenotype burden, 3.5 AlphaMissense residue aggregation, and 3.6
constraint-pocket prioritization match the current implementation.)

---

## Reference additions

- **Sehnal D, Bittrich S, Deshpande M, et al. Mol\* Viewer: modern web app for 3D visualization and analysis of
  large biomolecular structures. Nucleic Acids Res 2021;49:W431–7.** (Cited for the current viewer; the
  existing 3Dmol.js citation — Rego & Koes 2015 — can be retained to credit the earlier renderer or noted as
  superseded.)

---

## Notes for the authors (not for the manuscript)

- The per-substitution predictor panel, the pocket-evidence/similar-pocket modules, the partner-chain
  (multichain) handling, the ligand 2D/3D similarity links, and the chimeric flag are present in the current code
  but were absent or under-described in the v1.5.1 text; the edits above bring the manuscript in line with the code.
- The text already states the viewer uses Mol\*; only the **reference list** lacked a Mol\* citation (it cited
  3Dmol.js, the superseded renderer). The reference addition fixes that.
- No new external resources beyond those already listed are required, except naming the **RCSB Structure Motif
  Search** service used for the similar-pocket feature (same RCSB provider already cited).
