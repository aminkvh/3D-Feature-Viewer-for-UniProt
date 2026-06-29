# 3D Feature Viewer for UniProt — User Guide

This extension lets you look at a protein's 3D structure and its biological annotations — mutations, modification sites, binding pockets, and more — all without leaving the UniProt page. This guide explains what everything does and when you would use it.

---

## Table of Contents

1. [What This Extension Does](#what-this-extension-does)
2. [Opening the Viewer](#opening-the-viewer)
3. [Choosing a Structure](#choosing-a-structure)
4. [Navigating the 3D View](#navigating-the-3d-view)
5. [Annotation Panels — What You Can Show on the Structure](#annotation-panels)
   - [PTMs](#ptms-post-translational-modifications)
   - [Disease Variants](#disease-variants)
   - [Functional Sites](#functional-sites)
   - [Domains](#domains)
   - [Ligands](#ligands)
6. [Color Modes](#color-modes)
7. [Clicking a Residue — The Details Panel](#clicking-a-residue--the-details-panel)
8. [Clicking a Ligand — The Ligand Panel](#clicking-a-ligand--the-ligand-panel)
9. [Exporting Your Results](#exporting-your-results)
10. [Settings](#settings)

---

## What This Extension Does

UniProt is a great starting point for learning about a protein, but understanding *where* an annotation sits in 3D space usually means opening a separate viewer, finding the right structure, and manually re-entering position numbers. This extension handles all of that automatically — it reads the UniProt page you are already on, fetches the best available structure, maps all the annotations onto it, and shows everything in one window.

You can explore:
- Where disease-causing mutations cluster in the structure
- Which residues are near a binding pocket
- What effect a substitution is predicted to have
- Which ligands fit into the structure and where

---

## Opening the Viewer

When you visit a UniProt protein entry (for example `https://www.uniprot.org/uniprotkb/P14867/entry`), the extension adds **View in 3D** buttons directly next to relevant sections on the page. You do not need to do anything else — just click the button closest to the annotations you are interested in.

| Button | Where it appears | What it opens the viewer with |
|---|---|---|
| **View PTMs in 3D** | PTM / Processing section | Chemical modifications highlighted |
| **View Variants in 3D** | Disease & Variants section | Disease-associated mutations highlighted |
| **View Sites in 3D** | Function section | Active sites and binding residues highlighted |
| **View Domains in 3D** | Family & Domains section | Domain boundaries shown and colored |
| **View in 3D** | Structure section | Full structure, all layers available |
| **View in 3D** | Subcellular Location | Structure colored by membrane regions |

The viewer opens as a panel on the right side of the page. You can keep reading UniProt on the left while interacting with the 3D viewer on the right.

---

## Choosing a Structure

The **structure selector** is at the top of the viewer. Use the **← →** arrows or the dropdown menu to switch between available structures for this protein.

**What loads and in what order:**
1. **AlphaFold model** — An AI-predicted structure is available for almost every UniProt entry and loads first so you are never waiting.
2. **Experimental PDB structures** — Crystal structures, cryo-EM maps, and NMR models are discovered in the background and added to the list. These show where the protein has actually been observed.
3. **Isoform models** — If the protein has reviewed isoforms, their AlphaFold models appear here.
4. **Other computed models** — From SWISS-MODEL, ModelArchive, and similar sources.

Each entry in the dropdown shows how much of the protein it covers (as a percentage) and, for experimental structures, the resolution and method. A **⚛** icon means the structure contains chains from more than one organism — useful context if you are working with a complex.

If the protein is part of a multi-chain complex, the **Other chains** button lists the partner proteins and links to their UniProt pages. The **Partners** button lets you overlay annotations from those partner proteins onto the same 3D view.

You can set which type of structure loads by default in [Settings](#settings).

---

## Navigating the 3D View

**Mouse controls:**

| Action | What it does |
|---|---|
| Click a residue or atom | Opens the details panel for that position |
| Click a ligand | Opens the ligand panel |
| Double-click | Closes the details panel and resets the view |
| Scroll wheel | Zoom in and out |
| Click and drag | Rotate the structure |
| Right-click and drag | Move the view without rotating |

**Buttons in the top bar:**

- **Reset** — Returns everything to its starting state: default structure, default colors, no focused residue.
- **☀ / ☾** — Switch between a light and dark background. Useful depending on what you are trying to see.
- **📷 Screenshot** — Saves an image of the current view.
- **Tractability** — Shows drug development data for this protein from Open Targets, including whether approved drugs or clinical candidates exist, and whether the protein is considered a good small-molecule or antibody target.
- **Download** — All export options (see [Exporting Your Results](#exporting-your-results)).
- **✕** — Close the viewer.

**The sequence track** appears below the 3D canvas as a scrollable strip of amino acid letters. Each position is colored to show what is known about it — grey means it is not resolved in the loaded structure, blue means it is present. Orange underlines mark PTM sites; red dots mark variant positions. Click any resolved residue in the sequence strip to select it and focus the 3D view on it.

---

## Annotation Panels

The right-side panel shows the biological annotations you can display on the structure. Every section is collapsible — click the header to expand or collapse it. **All** and **None** buttons toggle everything in a section at once. The **C** button switches the coloring style for spheres in that layer.

### PTMs (Post-Translational Modifications)

Post-translational modifications are chemical changes that happen to a protein after it is made — phosphorylation, glycosylation, acetylation, and many others. These can switch a protein's activity on or off, change where it goes in the cell, or affect how it interacts with other molecules.

The PTM panel groups modifications by type. Expand any category to see individual sites, each with a zoom button that centers the 3D view on that residue.

You can also overlay disease variants, functional sites, domains, and ligands in the same view using the additional sections lower in the panel.

### Disease Variants

These are positions in the protein where a mutation has been linked to a disease, or where a computational tool predicts the change would be damaging.

**Filtering options:**
- **By disease** — Show only the variants associated with a disease you care about.
- **By consequence** — Filter by how confident the annotation is: Likely pathogenic, Predicted deleterious, Uncertain significance, Likely benign. These categories come from ClinVar (a clinical database) and computational predictions.
- **By provenance** — Filter by where the annotation came from: ClinVar submissions, published literature, or computational predictors.

Use these filters together to focus on, for example, only the high-confidence pathogenic variants in a specific disease.

### Functional Sites

These are residues that are directly involved in the protein's biological function — for example, the amino acids in an enzyme's active site that catalyze a reaction, or a metal-binding residue that coordinates a zinc ion.

Each site has a zoom button. Overlapping PTMs, variants, and domains can be shown at the same time.

### Domains

Domains are stretches of the protein that fold into a recognizable structural unit and often carry a specific function. The **C (color)** button colors the protein backbone by domain, giving each a different color so you can immediately see which part of the structure corresponds to which functional region.

Each domain entry has a zoom button to bring that region into view.

### Ligands

Ligands are small molecules bound to the protein in the loaded structure. For AlphaFold models, the extension uses **AlphaFill**, which transplants ligands from experimentally solved structures of similar proteins into the AlphaFold model to show where they would likely bind.

- Each ligand is listed by its chemical code (e.g. ATP, HEM).
- For transplanted ligands, you can see the **sequence identity** of the donor protein — how similar it was to yours — and the **donor PDB entry** it came from. Higher identity means the transplant is more likely to be accurate.
- Use the **exclude ions** toggle to hide small ions (Na⁺, Cl⁻, Mg²⁺, etc.) that are often just structural and not the main ligand of interest.
- Click the zoom button or ligand name to focus the view and open the [Ligand Panel](#clicking-a-ligand--the-ligand-panel).

---

## Color Modes

The **color mode** dropdown (top-right of the viewer) changes how the protein backbone is colored. It does not affect the annotation spheres — those keep their own colors.

| Mode | What it shows | Best used when |
|---|---|---|
| **Default** | Uniform color | You want annotation spheres to stand out clearly |
| **pLDDT confidence** | AlphaFold's confidence in each residue: blue = high confidence, red = low confidence | You are working with an AlphaFold model and want to know which regions to trust |
| **Experimental B-factor** | How much each residue moved in the crystal or solution: blue = rigid, red = flexible | You are working with an experimental structure |
| **Membrane topology** | Colors different parts of a membrane protein by location (inside cell, membrane, outside) | You are studying a transmembrane protein |
| **AlphaMissense summary** | Average predicted harmfulness of all possible mutations at each position | You want a quick overview of which regions are most sensitive to change |

**Exploratory modes** — These four are investigational tools to help generate hypotheses. They are shown only when enabled in [Settings](#settings).

| Mode | What it shows |
|---|---|
| **Pathogenic variant hotspots** | Positions where disease-causing mutations are more spatially clustered than expected by chance — potential functional or structurally important regions |
| **Contact-network centrality** | Residues that act as structural hubs, with many of the shortest communication paths in the protein passing through them — candidates for allosteric or stability-critical positions |
| **Recurrent phenotype residues** | Positions that accumulate multiple different disease labels across independent reports — suggesting broad functional importance |
| **Constraint pocket clusters** | Groups of residues that form buried, evolutionarily constrained cavities — candidate binding sites even without known ligands. A sensitivity slider lets you make this more or less stringent. |

These four modes highlight regions worth investigating further. They are not clinical predictions.

---

## Clicking a Residue — The Details Panel

Click any residue in the 3D canvas or the sequence strip to open the details panel. It pulls together everything known about that specific position.

### Residue header

Shows the amino acid name and position (e.g. **ALA 421**). Small colored flags appear if any exploratory algorithm flagged this residue:

- 🔴 Pathogenic variant hotspot
- 🟠 Recurrent phenotype residue
- 🟣 Contact-network hub
- 🟢 Constraint pocket cluster

### Nearby residues

A distance slider (default ~5 Å, adjustable from 2–30 Å) lists all residues within that distance in the structure. Each nearby residue is labeled with its annotation colors. Click any of them to refocus the view there.

This is useful for understanding local structural context — for example, whether a disease mutation is sitting next to an active site residue.

### PTMs, Sites, and Mutagenesis

If this position has a post-translational modification, a functional site annotation, or an experimentally tested mutagenesis result, it appears here as labeled chips.

### Variants

All mutations annotated at this position are listed. For each one you can see:
- The amino acid change (e.g. W123M — tryptophan to methionine at position 123)
- The associated disease(s)
- The clinical significance (Pathogenic, Likely pathogenic, Uncertain, etc.)

Click **Show evidence** to expand each variant and see: the ClinVar review level (how many independent submissions support this classification), the dbSNP identifier, how common the variant is in the general population (gnomAD frequency — a very rare variant in healthy people is more likely to be disease-causing), and the genomic coordinates.

### PTM–Variant Proximity

Two sliders let you search for PTMs or pathogenic variants near this position within a chosen radius. Results appear as colored chips you can click to jump to that position. This helps you spot cases where a disease mutation is close to a modification site — which can suggest they interfere with each other.

### Binding & Pockets

Two sources of information about whether this residue is near a binding site:

**Predicted pockets** — Computed from the structure. For each pocket at or near this residue:
- How many residues line the pocket, with a button to highlight them in 3D
- How buried the pocket is (0 = fully exposed on the surface, 1 = deeply buried)
- How confident the structure is in this region (pLDDT for AlphaFold, or B-factor for experimental)
- The pocket's physical size (radius of gyration in Ångströms — roughly how big the cavity is)
- The amino acid composition of the pocket: what fraction is hydrophobic, aromatic, positively charged, negatively charged, or polar. This tells you what kind of molecule might bind there.
- **Find Similar Motifs** — Searches the entire Protein Data Bank for structures with a similar pocket geometry. Useful for finding proteins with known ligands that might bind here too.

**Experimental binding** (from PDBe-KB) — If any experimental PDB structure has a ligand or protein-protein interface at this position, it is listed here with links to those structures. This is the strongest evidence that a pocket is real and relevant.

### Predictions

Computational predictions of how harmful a mutation at this position would be, sourced from ProtVar.

**Per-position scores (apply to the whole position):**
- **Conservation** (0–1): How consistent this amino acid is across evolution. A score close to 1 means this position is almost never mutated across species — suggesting it is important.
- **M3D**: A structural predictor that tells you whether a change here is likely to be damaging, and which structural feature (e.g. buried hydrophobic core, active site proximity) makes it sensitive.

**Per-substitution table** (one row per possible amino acid swap): For each of the 19 possible changes you could make at this position:
- **AlphaMissense** — Google DeepMind's model trained on evolutionary and structural data. Scores range 0–1; above 0.564 is considered likely pathogenic.
- **EVE** — An evolutionary model that looks at how often this change appears (or not) across thousands of related sequences.
- **CADD** — A score combining many different signals. Higher numbers mean a more likely damaging change (typically > 20 is considered notable).
- **ESM-1b** — A language model trained on protein sequences, predicting whether this substitution is unusual relative to what evolution has allowed.
- **FoldX ΔΔG** — Predicted change in folding stability in kcal/mol. Positive values mean the mutation destabilizes the protein.

Variants already documented in UniProt are shown in bold with a marker. Click any row to open the full ProtVar record.

---

## Clicking a Ligand — The Ligand Panel

Click any ligand in the 3D view or in the Ligands section to open the ligand panel.

If the same ligand appears in more than one location in the structure (for example, an enzyme with two active sites), use the **◄ ►** arrows to move between copies.

### Nearby residues

Same distance slider as the residue panel. Lists protein residues within range of the ligand — these are the residues most likely to be important for binding.

### AlphaFill transplant evidence

For ligands placed by AlphaFill (into AlphaFold models):
- **Sequence identity** — How similar the donor protein was to this one. Higher identity = more reliable transplant.
- **Donor PDB entry** — The experimental structure the ligand came from.
- **Clash score** — Whether the transplanted ligand fits cleanly (low clash) or overlaps with the model atoms (high clash, less reliable).

### Chemical information

The ligand's name, molecular formula, molecular weight, and how many hydrogen-bond donors and acceptors it has. Copy buttons for the SMILES string (a text representation of the molecule's structure) and InChIKey (a unique chemical identifier useful for searching databases).

### Pocket evidence

Same predicted pocket and experimental binding data as the residue panel, focusing on the residues surrounding this ligand.

### External links

- **PubChem** — Detailed chemistry and bioactivity data for this compound
- **DrugBank** — Drug and target information if this is a known pharmaceutical
- **Similarity search** — Find structurally similar molecules in PubChem (2D or 3D similarity)

### Ligand similarity

Other ligands in the same structure ranked by how chemically similar they are to this one (using a standard molecular fingerprint comparison called Tanimoto similarity — a score of 1.0 means identical, 0 means completely different). Click any entry to focus on that ligand.

---

## Exporting Your Results

Click **Download** in the top bar to access all export options.

### CSV — Annotation table

A spreadsheet with one row per residue and one column per annotation type. Useful for further analysis in Python, R, or Excel. Columns include: UniProt position, amino acid, PTMs (one column per category), diseases (one column per disease), variant counts, how common variants are in the healthy population (gnomAD frequency), AlphaMissense summary scores, exploratory algorithm results, and which ligands are nearby.

**Download CSV + ProtVar predictions** — Same file but also includes the per-residue conservation, stability, and effect prediction scores from ProtVar for every residue in the protein. This takes longer because it queries ProtVar's API for the whole sequence.

### PyMOL session

Downloads a script file that you can open in PyMOL to reproduce the exact view from the extension. The protein will appear with the same coloring, the same annotation spheres, and — if you had a residue selected — the same zoomed focus. Useful for making publication-quality figures or continuing analysis in PyMOL.

### VMD session

Same as the PyMOL session but for VMD, a molecular dynamics visualization program widely used in computational biophysics.

### Copy to clipboard

Copies a selection string (in PyMOL or VMD format) listing the currently visible annotated residues. You can paste this directly into PyMOL or VMD to select those residues in a structure you already have open. The format (PyMOL or VMD) is set in [Settings](#settings).

---

## Settings

Open the settings page from your browser's extension management (the puzzle-piece icon in Chrome, or the extensions menu in Firefox). Changes take effect next time you open the viewer.

| Setting | What it controls |
|---|---|
| **Default structure** | Whether AlphaFold, the best experimental structure, or the highest-coverage structure loads first |
| **Default color mode** | The coloring applied when the viewer first opens |
| **Clipboard format** | Whether copied selections use PyMOL or VMD syntax |
| **PTM search radius** | The default distance used when searching for nearby PTMs in the details panel |
| **Variant search radius** | The default distance used when searching for nearby variants |
| **Font size** | Makes the text in the side panels smaller or larger |
| **Show optional tracks** | Enables additional conservation and population frequency displays |
| **Show exploratory algorithms** | When turned off, the four investigational color modes (hotspots, hubs, burden, pockets) are hidden from the color dropdown — useful if you want a simpler interface focused on direct annotations |
