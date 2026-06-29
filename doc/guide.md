# 3D Feature Viewer for UniProt: User Guide

This guide explains what the extension does and how to use it.

---

## Table of Contents

1. [What This Extension Does](#what-this-extension-does)
2. [Opening the Viewer](#opening-the-viewer)
3. [Choosing a Structure](#choosing-a-structure)
4. [Navigating the 3D View](#navigating-the-3d-view)
5. [Annotation Panels: What You Can Show on the Structure](#annotation-panels)
   - [PTMs](#ptms-post-translational-modifications)
   - [Disease Variants](#disease-variants)
   - [Functional Sites](#functional-sites)
   - [Domains](#domains)
   - [Ligands](#ligands)
6. [Color Modes](#color-modes)
7. [Clicking a Residue: The Details Panel](#clicking-a-residue-the-details-panel)
8. [Clicking a Ligand: The Ligand Panel](#clicking-a-ligand-the-ligand-panel)
9. [Exporting Your Results](#exporting-your-results)
10. [Settings](#settings)

---

## What This Extension Does

UniProt annotates proteins richly: post-translational modifications, disease-associated variants, functional sites, domain boundaries, and more. These annotations are presented as flat text lists organized by category. What is hard to see from a list is how these annotations relate to each other in three-dimensional space, and what their structural context is.

This extension maps all of those annotation layers onto the protein's 3D structure simultaneously, directly inside the UniProt page. You can see a disease-associated variant at position 52, a phosphorylation site at position 50, and a known active-site residue at position 48 all at once, on the structure, with their spatial relationships visible. Switching between annotation types, color modes, and structures happens without leaving the page or re-entering residue numbers.

The viewer also provides, for each residue you click, aggregated information from several external databases: clinical variant classifications, per-substitution effect predictions from multiple models, binding pocket evidence, and drug tractability data.

---

## Opening the Viewer

When you visit a UniProt protein entry (for example `https://www.uniprot.org/uniprotkb/P14867/entry`), the extension adds **View in 3D** buttons directly next to relevant sections on the page. Click the button closest to the annotation type you want to start with.

| Button | Where it appears | What it opens the viewer with |
|---|---|---|
| **View PTMs in 3D** | PTM / Processing section | PTM sites shown as colored spheres |
| **View Variants in 3D** | Disease & Variants section | Disease-associated variants shown as colored spheres |
| **View Sites in 3D** | Function section | Active sites and binding residues shown |
| **View Domains in 3D** | Family & Domains section | Domain boundaries shown and colored |
| **View in 3D** | Structure section | Full structure, all annotation layers available |
| **View in 3D** | Subcellular Location | Structure colored by membrane topology |

The viewer opens on the right side of the page. The UniProt page remains on the left so you can read annotations and see them on the structure at the same time.

---

## Choosing a Structure

The **structure selector** is at the top of the viewer. Use the **← →** arrows or the dropdown menu to switch between available structures.

**What loads and in what order:**
1. **AlphaFold model:** A computationally predicted structure from AlphaFold (Jumper et al., 2021, *Nature*) is available for most reviewed UniProt entries and loads first. AlphaFold models cover the full sequence but are predicted; confidence varies by region (see pLDDT in Color Modes).
2. **Experimental PDB structures:** Crystal structures, cryo-EM maps, and NMR ensembles from the PDB are discovered in the background. These are determined from physical experiments but may cover only part of the sequence, may be from a truncated or mutant construct, and represent the protein under specific experimental conditions.
3. **Isoform models:** AlphaFold models for reviewed isoforms, where available.
4. **Other computed models:** From SWISS-MODEL, ModelArchive, and similar repositories.

Each entry in the dropdown shows the sequence coverage (as a percentage) and, for experimental structures, the method and resolution. A **⚛** icon means the structure contains chains from more than one organism.

For multi-chain structures, the **Other chains** button lists partner proteins with links to their UniProt pages. The **Partners** button overlays annotations from partner chains onto the same view.

Which structure type loads by default can be set in [Settings](#settings).

---

## Navigating the 3D View

**Mouse controls:**

| Action | What it does |
|---|---|
| Click a residue | Opens the details panel for that position |
| Click a ligand | Opens the ligand panel |
| Double-click | Closes the details panel and resets the view |
| Scroll wheel | Zoom in and out |
| Click and drag | Rotate the structure |
| Right-click and drag | Pan without rotating |

**Buttons in the top bar:**

- **Reset:** Returns to the default structure, default coloring, and no focused residue.
- **☀ / ☾:** Switch between light and dark background.
- **📷 Screenshot:** Save an image of the current view.
- **Tractability:** Drug tractability assessment for this protein from Open Targets (Ochoa et al., 2021, *Nat Genet*), including approved drugs, clinical candidates, and small-molecule or antibody tractability assessments.
- **Download:** All export options (see [Exporting Your Results](#exporting-your-results)).
- **✕:** Close the viewer.

**The sequence track** below the 3D canvas shows the full protein sequence as a scrollable strip. Positions are colored to show annotation status. Grey means the residue is not present in the loaded structure; blue means it is modeled. Orange underlines mark PTM sites; red dots mark variant positions. Click any modeled residue to select it and center the 3D view on it.

---

## Annotation Panels

The right-side panel controls which annotation layers are visible on the structure. Every section is collapsible. **All** and **None** toggle all items in a section. The **C** button changes the coloring scheme for spheres in that layer.

The key capability here is layering: you can show PTMs, disease-associated variants, functional sites, and domains simultaneously on the same structure, with different sphere colors per layer, so their spatial relationships become visible.

### PTMs (Post-Translational Modifications)

Post-translational modifications are covalent chemical changes made to a protein after translation, such as phosphorylation, glycosylation, ubiquitination, and acetylation. These are annotated in UniProt from the literature and curated databases.

The PTM panel groups modifications by type. Expand any category to see individual sites. Each row has a zoom button that centers the 3D view on that residue.

Other annotation layers (variants, functional sites, domains, ligands) can be overlaid in the same view using the additional sections in the panel.

### Disease Variants

These are amino acid positions where variants have been reported in the context of disease or where computational tools predict a change would be deleterious. The annotations are drawn from UniProt, which integrates ClinVar, the published literature, and computational predictors.

**Filtering options:**
- **By disease:** Restrict to variants associated with a specific disease.
- **By consequence:** Filter by annotation confidence. Categories include Likely pathogenic, Predicted deleterious, Uncertain significance, and Likely benign, following ClinVar's classification framework.
- **By provenance:** Filter by annotation source (ClinVar submissions, published literature, or computational predictors).

### Functional Sites

Residues annotated in UniProt as directly involved in the protein's biochemical activity: catalytic residues, ligand-binding residues, metal-coordinating residues, and similar functional annotations derived from the literature and curated resources such as the Catalytic Site Atlas.

### Domains

Annotated domains, regions, coiled coils, transmembrane segments, and signal peptides from UniProt. The **C (color)** button assigns a distinct color to each domain along the protein backbone, making it easy to see which part of the 3D structure corresponds to which annotated region.

### Ligands

Small molecules present in the loaded structure. For AlphaFold models, the extension integrates **AlphaFill** (Hekkelman et al., 2023, *Nat Methods*), which transplants ligands from experimentally determined structures of homologous proteins into the AlphaFold model at structurally equivalent positions.

- Each ligand is listed by its CCD code (e.g. ATP, HEM, ZN).
- For AlphaFill transplants, the donor protein's sequence identity to the query and the source PDB entry are shown. Transplant accuracy is expected to increase with sequence identity, but the placement is a structural inference, not an experimentally determined position for this protein.
- The **exclude ions** toggle hides monoatomic ions (Na+, Cl-, Mg2+, etc.) that are often crystallographic rather than functional.
- Click the zoom button or ligand name to focus the view and open the [Ligand Panel](#clicking-a-ligand-the-ligand-panel).

---

## Color Modes

The **color mode** dropdown changes how the protein backbone is colored. Annotation spheres keep their own colors regardless of the backbone color mode.

| Mode | What it shows |
|---|---|
| **Default** | Uniform backbone color; annotation spheres are most visible here |
| **pLDDT score** | Per-residue pLDDT from the AlphaFold model, a measure of local structural confidence used by AlphaFold authors (Jumper et al., 2021): blue (pLDDT > 90, high local confidence) to red (pLDDT < 50, low confidence, often disordered). Not applicable to experimental structures. |
| **Experimental B-factor** | Crystallographic or NMR B-factor per residue, reflecting atomic displacement in the experimental data: blue (low B-factor, more ordered) to red (high B-factor, more disordered or flexible). Applicable to experimental structures only. |
| **Membrane topology** | Colors transmembrane, cytoplasmic, and extracellular segments when topology annotations are present in UniProt. |
| **AlphaMissense summary** | Mean AlphaMissense pathogenicity score (Cheng et al., 2023, *Science*) across all 19 possible amino acid substitutions at each position. A high mean score indicates that most substitutions at that position are predicted to be damaging, suggesting the position is functionally constrained. |

**Exploratory modes** are research-oriented overlays for hypothesis generation. They appear only when enabled in [Settings](#settings). None of these constitute clinical evidence.

| Mode | What it shows |
|---|---|
| **Pathogenic variant hotspots** | Residues where disease-associated variants are more spatially concentrated than expected by permutation, using an approach analogous to published spatial clustering methods (e.g. Kamburov et al., 2015, *PNAS*). Tiered by statistical support (strong, moderate, weak). |
| **Contact-network centrality** | Residues with high betweenness centrality in the Ca-Ca contact graph (8 Å cutoff), computed with Brandes' exact algorithm. High-centrality residues lie on many shortest paths through the contact network and may be relevant to structural communication or stability. |
| **Recurrent phenotype residues** | Positions that accumulate multiple distinct disease or phenotype labels across independent variant reports, scored by a composite of variant count and phenotype diversity. |
| **Constraint pocket clusters** | Residue groups forming geometrically buried, evolutionarily constrained cavities, identified using a heuristic analogous to pocket detection methods. A sensitivity slider controls the FDR threshold. |

---

## Clicking a Residue: The Details Panel

Click any modeled residue in the 3D canvas or the sequence strip to open the details panel. It aggregates annotations and predictions for that specific position from multiple sources.

### Residue header

Shows the three-letter amino acid code and UniProt sequence position (e.g. **ALA 421**). Colored flags indicate which exploratory overlays have flagged this residue:

- Red: Pathogenic variant hotspot
- Orange: Recurrent phenotype residue
- Purple: Contact-network hub
- Teal: Constraint pocket cluster

### Nearby residues

A distance slider (default 5 Å, range 2 to 30 Å) lists all residues within that radius in the current structure. Each entry is color-coded by its annotation type. Click any residue in the list to move the focus there.

This is one of the core uses of the extension: selecting a disease-associated variant and immediately seeing which PTMs, functional sites, and other variants are spatially close to it on the structure.

### PTMs, Sites, and Mutagenesis

PTM annotations, functional site annotations, and experimental mutagenesis results at this position are shown as labeled chips.

### Variants

All variants annotated at this position in UniProt. For each variant:
- The amino acid change (e.g. W123M: tryptophan to methionine at position 123)
- Associated disease or phenotype labels
- Clinical significance classification (Pathogenic, Likely pathogenic, Variant of uncertain significance, Likely benign, Benign) following ClinVar's five-tier framework

Click **Show evidence** to expand each variant and see: the ClinVar review status (number of independent submitting laboratories), the dbSNP identifier, the gnomAD population allele frequency (Karczewski et al., 2020, *Nature*), and the genomic coordinates. A variant that is very rare in large population cohorts is more consistent with pathogenicity than one that is common, though rarity alone does not establish disease association.

### PTM–Variant Proximity

Two distance sliders search for PTMs or disease-associated variants within a configurable radius of the selected residue. Results are shown as colored chips. This highlights spatial co-occurrence of different annotation types, such as a variant landing adjacent to a phosphorylation site.

### Binding & Pockets

**Predicted pockets** (from ProtVar/AutoSite): Geometric cavities identified computationally from the structure. For each pocket near this residue:
- The residues lining the pocket (click to highlight them in 3D)
- Buriedness (0 = surface-exposed, 1 = fully buried interior)
- pLDDT or B-factor at pocket residues, as a local confidence indicator
- Radius of gyration (a measure of the cavity's spatial extent)
- Amino acid composition of the pocket lining (fractions hydrophobic, aromatic, acidic, basic, polar)
- **Find Similar Motifs:** Links to RCSB's structural motif search for this pocket geometry across all PDB structures

These are predicted cavities. Their relevance as binding sites depends on additional evidence.

**Experimental binding** (from PDBe-KB): Ligand contacts and protein-protein interface residues observed at this position in experimentally determined structures in the PDB. This is direct structural evidence that the position contacts a ligand or an interacting protein.

### Predictions

Per-position and per-substitution computational predictions from ProtVar (Nightingale et al., 2023, *Nat Biotechnol*).

**Per-position:**
- **Conservation** (0 to 1): Evolutionary conservation score. A score near 1 indicates the amino acid at this position is nearly invariant across homologous sequences, consistent with functional or structural constraint.
- **M3D:** A structure-based predictor of mutational impact at this position, including the structural feature identified as the main driver (e.g. buried position, active site proximity).

**Per-substitution table** (one row per possible amino acid change):
- **AlphaMissense** (Cheng et al., 2023, *Science*): A missense pathogenicity score from 0 to 1. Values above 0.564 are classified as likely pathogenic; below 0.34 as likely benign; the interval between is ambiguous. Trained on evolutionary and structural features.
- **EVE** (Fraternali et al. / Fraternali group): An unsupervised evolutionary model scoring variant fitness from deep sequence alignments. More negative values indicate a substitution less consistent with the evolutionary record.
- **CADD** (Kircher et al., 2014, *Nat Genet*): A combined annotation score integrating many genomic and functional features. Phred-scaled; values above 20 correspond to the top 1% of deleteriousness across the human genome.
- **ESM-1b** (Rives et al., 2021, *PNAS*): A protein language model likelihood ratio score. More negative values indicate a substitution less expected by the model trained on millions of protein sequences.
- **FoldX DDG** (Schymkowitz et al., 2005, *Nucleic Acids Res*): Predicted change in folding free energy in kcal/mol. Positive values indicate the substitution is predicted to destabilize the folded structure.

Variants already documented in UniProt are shown in bold. Click any row to open the full ProtVar record.

---

## Clicking a Ligand: The Ligand Panel

Click any ligand in the 3D view or in the Ligands section to open the ligand panel.

If the same ligand appears in more than one binding site in the structure, use the left and right arrows to navigate between copies.

### Nearby residues

Same distance slider as the residue panel. Lists protein residues within range of the ligand, annotated by their type.

### AlphaFill transplant evidence

For ligands placed by AlphaFill (Hekkelman et al., 2023, *Nat Methods*):
- **Sequence identity:** The pairwise sequence identity between the donor protein and this protein. Higher identity supports a more reliable structural transfer.
- **Donor PDB entry:** The experimentally determined structure the ligand was taken from.
- **Clash score:** An estimate of how well the transplanted ligand fits the AlphaFold model geometry (low = few steric conflicts, high = significant overlap with model atoms).

These are structural inferences. Binding of this ligand to this protein has not necessarily been demonstrated experimentally.

### Chemical information

Loaded from RCSB: the ligand's IUPAC name, molecular formula, molecular weight, hydrogen-bond donor and acceptor counts. Copy buttons are provided for the SMILES string and InChIKey.

### Pocket evidence

Same predicted pocket and experimental binding data as the residue panel, restricted to the region surrounding this ligand.

### External links

- **PubChem:** Compound record with bioassay and literature data
- **DrugBank:** Drug and pharmacological target information where applicable
- **Similarity search:** PubChem 2D or 3D similarity search for structurally related compounds

### Ligand similarity

Other ligands in the same structure ranked by Tanimoto similarity of their CACTVS 881-bit fingerprints (a standard 2D structural fingerprint). A score of 1.0 indicates identical fingerprints; 0 indicates no shared bits. Click any entry to navigate to that ligand.

---

## Exporting Your Results

Click **Download** in the top bar for all export options.

### CSV annotation table

A tabular file with one row per UniProt residue position. Columns include: amino acid identity, PTM annotations (one column per category), disease-associated variant counts and labels, functional site flags, gnomAD population frequencies, AlphaMissense per-position statistics, exploratory algorithm outputs, and nearby ligand CCD codes.

**Download CSV + ProtVar predictions** extends this with per-residue conservation scores and per-position EVE, ESM-1b, FoldX, and CADD summaries from ProtVar. This requires one additional network request to the ProtVar API.

### PyMOL session

A `.pml` script that reproduces the current view in PyMOL (Schrodinger LLC). The protein cartoon, annotation spheres, ligand representations, and (if a residue is focused) the zoomed selection are all included with the same coloring as the extension. Useful for figure preparation or further analysis.

### VMD session

A `.vmd` script that reproduces the current view in VMD (Humphrey et al., 1996, *J Mol Graph*).

### Copy to clipboard

A residue selection string in PyMOL (`resi` syntax) or VMD (`resid` syntax), listing the currently visible annotated residues. The format is set in [Settings](#settings).

---

## Settings

Open the settings page from your browser's extension management menu. Changes apply the next time you open the viewer.

| Setting | What it controls |
|---|---|
| **Default structure** | Which structure type loads first: AlphaFold model, experimental (PDB), or highest sequence coverage |
| **Default color mode** | The backbone coloring applied when the viewer opens |
| **Clipboard format** | Whether copied selections use PyMOL or VMD syntax |
| **PTM search radius** | Default radius for nearby PTM search in the details panel |
| **Variant search radius** | Default radius for nearby variant search in the details panel |
| **Font size** | Text size in the annotation panels |
| **Show exploratory algorithms** | When off, the four investigational color modes are hidden from the color dropdown |
