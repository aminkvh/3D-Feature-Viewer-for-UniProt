# 3D Feature Viewer for UniProt — User Guide

This guide covers everything the extension can do. If you are new, start with [Getting Started](#getting-started). Everything else is reference.

---

## Table of Contents

1. [Getting Started](#getting-started)
2. [Opening the Viewer](#opening-the-viewer)
3. [Choosing a Structure](#choosing-a-structure)
4. [The 3D Viewer](#the-3d-viewer)
5. [Annotation Panels](#annotation-panels)
   - [PTMs](#ptm-panel)
   - [Disease Variants](#disease--variants-panel)
   - [Functional Sites](#functional-features-panel)
   - [Family & Domains](#family--domains-panel)
   - [Ligands](#ligands)
6. [Color Modes](#color-modes)
7. [Residue Details Panel](#residue-details-panel)
8. [Ligand Details Panel](#ligand-details-panel)
9. [Export](#export)
10. [Settings](#settings)

---

## Getting Started

Install the extension from the [Chrome Web Store](https://chromewebstore.google.com/detail/uniprot-3d-feature-viewer/fplpkbigppbpbcdmpkefdoilfgcicaof) or [Mozilla Add-ons](https://addons.mozilla.org/en-US/firefox/addon/3d-feature-viewer-for-uniprot/), then open any UniProt protein entry page (e.g. `https://www.uniprot.org/uniprotkb/P14867/entry`). You will see **View in 3D** buttons appear next to relevant sections on the page.

---

## Opening the Viewer

The extension adds buttons directly into UniProt pages. Each button opens the viewer pre-filtered to that annotation type:

| Button | Where it appears | What it shows |
|---|---|---|
| **View PTMs in 3D** | PTM / Processing section | Post-translational modification sites |
| **View Variants in 3D** | Disease & Variants section | Disease-associated variants |
| **View Sites in 3D** | Function section | Active sites, binding sites, metal sites |
| **View Domains in 3D** | Family & Domains section | Domain and region boundaries |
| **View in 3D** | Structure section | Full structure with all layers available |
| **View in 3D** | Subcellular Location section | Structure colored by membrane topology |

Clicking any button opens the 3D viewer window on the right side of the page.

---

## Choosing a Structure

The **structure selector** appears at the top of the viewer. Use the **← →** arrows or the dropdown to switch between available structures.

Structures are loaded in this order:
1. **AlphaFold model** — loads immediately (predicted structure for every UniProt entry)
2. **Experimental PDB structures** — discovered and listed after the initial load; ordered by coverage
3. **Isoform AlphaFold models** — for proteins with reviewed isoforms
4. **Computed models** — from SWISS-MODEL, ModelArchive, and similar providers

Each structure entry shows:
- The source (AlphaFold, PDB accession, or model provider)
- Chain identifiers and residue coverage percentage
- Experimental method (X-ray, cryo-EM, NMR) and resolution where applicable
- A **⚛** flag if the structure contains chains from more than one organism

Your default structure preference (AlphaFold first, Experimental first, or Best coverage first) can be set in [Settings](#settings).

**Multi-chain structures:** For PDB structures with multiple proteins, the **Other chains** button lists all partner proteins with links to their UniProt entries. Use the **Partners** button to toggle annotation layers (variants, PTMs, sites) from partner chains onto the 3D view.

---

## The 3D Viewer

### Controls

| Action | What it does |
|---|---|
| **Click a residue** | Opens the details panel for that residue |
| **Click a ligand** | Opens the ligand details panel |
| **Double-click** | Closes the details panel and resets the focus |
| **Scroll / pinch** | Zoom in and out |
| **Click and drag** | Rotate the structure |
| **Right-click drag** | Pan |

### Header Buttons

- **Reset** — Returns to the default view: original structure, default colors, no focused residue.
- **Theme** (☀/☾) — Switch between light and dark background.
- **Screenshot** (📷) — Save an image of the current 3D view.
- **Tractability** — Open Targets drug tractability data for this protein: approved drugs, clinical candidates, and pocket-based small-molecule or antibody assessments.
- **Download** — Export options (see [Export](#export)).
- **✕** — Close the viewer.

### Sequence Track

Below the 3D canvas, a scrollable sequence ribbon shows the full protein sequence. Residues are color-coded:

- **Light grey** — not resolved in the loaded structure
- **Blue** — resolved in the structure
- **Orange underline** — has a PTM annotation
- **Red dot** — has a disease variant
- **Highlighted** — currently visible in the active annotation layer
- **Bold / outlined** — currently selected residue

Click any resolved residue in the sequence track to open its details panel and focus the 3D view on it.

---

## Annotation Panels

The right-side panel shows annotation layers relevant to the current view mode. Each section is collapsible. Use **All / None** to toggle all items at once, and the **C** button to switch sphere coloring within that layer.

### PTM Panel

Shows all post-translational modification sites annotated in UniProt for this protein.

- **PTM types** — Categories (e.g. Phosphoserine, N-linked glycosylation) with individual site rows. Each row shows the position and modification description. Click the zoom (🔍) icon to focus the 3D view on that residue.
- **Functional sites** — Active sites, binding sites, metal-binding sites. Shown as a separate sub-section.
- **Ligands** — Ligands in the loaded structure (see [Ligands](#ligands) below).
- **Disease variants** — Optionally overlay disease variant spheres alongside PTMs.
- **Family & Domains** — Optionally overlay domain annotation.

### Disease & Variants Panel

Shows disease-associated and pathogenicity-annotated variants from UniProt and ClinVar.

- **Disease filter** — Group variants by disease. Toggle individual diseases or use All/None. Expand a disease to see each variant with zoom buttons.
- **Provenance filter** — Filter by annotation source (ClinVar, literature, predicted).
- **Consequence filter** — Filter by variant effect category: Likely pathogenic, Predicted deleterious, Uncertain significance, Likely benign, etc.
- **PTM co-display** — Optionally show PTMs alongside variants.
- **Functional sites / Ligands / Domains** — Same overlay options as other panels.

### Functional Features Panel

Shows catalytic and binding sites from UniProt's functional annotation.

- Active sites, binding sites, metal-binding sites, and other functional regions.
- Zoom buttons per site.
- Overlay PTMs, variants, and domains as additional layers.

### Family & Domains Panel

Shows annotated domains, coiled coils, transmembrane regions, signal peptides, and other linear features.

- **C (color)** — Colors the protein cartoon backbone by domain, assigning each domain a distinct color.
- Zoom buttons bring each domain region into focus.
- Overlay PTMs, variants, and sites.

### Ligands

Ligands appear as a collapsible section in all panels (when the loaded structure contains ligands). Most commonly appears with AlphaFill models.

- Each ligand is grouped by its chemical ID (CCD code), with individual entries for each copy in the structure.
- **AlphaFill badge** — Shows the sequence identity (%) of the donor protein from which the ligand was transplanted, and the donor PDB ID.
- **Ion exclusion toggle** — Hide small ions (Na⁺, Cl⁻, etc.) to reduce clutter.
- Click the zoom button or the ligand name to focus the 3D view on that ligand and open the [Ligand Details Panel](#ligand-details-panel).

---

## Color Modes

The color mode dropdown (top-right of the viewer) changes how the structure backbone is colored. It does not affect annotation spheres.

| Mode | What it shows |
|---|---|
| **Default** | Uniform color; annotation spheres stand out clearly |
| **pLDDT confidence** | AlphaFold model confidence per residue: blue (high) → red (low). Not available for experimental structures. |
| **Experimental B-factor** | Structural flexibility from experimental data: blue (rigid) → red (flexible). PDB structures only. |
| **Membrane topology** | Colors transmembrane, cytoplasmic, and extracellular regions when topology annotation is present. |
| **AlphaMissense summary** | Mean predicted pathogenicity across all substitutions at each position. |

**Exploratory modes** (shown only when enabled in [Settings](#settings)):

| Mode | What it shows |
|---|---|
| **Pathogenic variant hotspots** | Residues where pathogenic variants cluster in 3D space, beyond what would be expected by chance. Tiered: strong, moderate, weak. |
| **Contact-network centrality** | Residues that act as structural bridges in the protein's residue contact network — positions where many shortest communication paths pass through. |
| **Recurrent phenotype residues** | Positions that accumulate multiple distinct disease or phenotype labels across different variant records. |
| **Constraint pocket clusters** | Residue groups that form geometrically buried, sequence-constrained cavities — candidate binding sites. A sensitivity slider (FDR threshold) controls how stringent the call is. |

These four modes are investigational heuristics for hypothesis generation. They are not clinical classifiers.

---

## Residue Details Panel

Click any residue in the 3D viewer or the sequence track to open the details panel on the right side.

### Header

Shows the three-letter amino acid code and UniProt position (e.g. **ALA 421**). Colored flags indicate which exploratory algorithm tiers this residue has reached:

- 🔴 Pathogenic variant hotspot
- 🟠 Recurrent phenotype residue
- 🟣 Contact-network hub
- 🟢 Constraint pocket cluster

### Nearby

An interactive distance slider (2–30 Å). Residues within the selected radius of the clicked residue are listed with their distance. Each nearby residue is color-coded by its annotations. Click any nearby residue to refocus the view.

### PTMs / Sites / Mutagenesis

Any PTM, functional site, or experimental mutagenesis annotation at this position is shown as labeled chips.

### Variants

Lists all variants annotated at this position. Each variant shows:
- The amino acid change (e.g. W123M), colored by ClinVar significance
- Associated diseases
- **Show evidence** — expands ClinVar review status, dbSNP IDs, gnomAD allele frequency, and genomic coordinates

### PTM–Variant Proximity

Two configurable sliders let you search for PTMs or pathogenic variants within a chosen radius of the selected residue. Results appear as clickable chips that refocus the view.

### Binding & Pockets

Two data sources combined:

**Predicted pockets** — Computed directly from the structure. For each pocket at or near this residue:
- Number of lining residues (click to highlight in 3D)
- Buriedness (0 = surface, 1 = fully buried)
- pLDDT or B-factor confidence of the pocket
- Radius of gyration (cavity size)
- Amino acid composition (hydrophobic / aromatic / acidic / basic / polar)
- **Find Similar Motifs** — launches an RCSB Structure Motif search for this pocket geometry

**Experimental binding** (from PDBe-KB) — Known ligands and protein-protein interfaces observed at this position in experimental structures, with links to the relevant PDB entries.

### Predictions

Per-residue and per-substitution effect predictions sourced from ProtVar:

- **Conservation** (0–1): evolutionary conservation across homologues
- **M3D** prediction: whether the substitution is predicted damaging, and the structural feature driving that call
- **Per-substitution table**: for each of the 19 possible amino acid changes at this position —
  - **AlphaMissense** score (0–1; ≥ 0.564 = likely pathogenic)
  - **EVE** score (evolutionary model)
  - **CADD** phred score (combined annotation)
  - **ESM-1b** score (language model)
  - **FoldX ΔΔG** (predicted folding stability change in kcal/mol)
  - Variants observed in UniProt are shown in bold with a dot marker
  - Click any row to view full details on ProtVar

---

## Ligand Details Panel

Click any ligand in the 3D viewer or in the Ligands section to open the ligand panel.

### Navigation

Use **◄ ►** to cycle between copies of the same ligand (for structures where a ligand appears in multiple binding sites). The counter shows which copy you are viewing.

### Nearby

Same distance slider as the residue panel. Shows protein residues near the ligand, color-coded by annotation.

### AlphaFill Evidence

For transplanted ligands from AlphaFill:
- Donor protein's sequence identity to the query
- Donor PDB ID
- Transplant clash assessment (low / moderate / high steric clash with the model)

### Chemical Identity

Loaded asynchronously from RCSB:
- Full chemical name
- Molecular formula
- Molecular weight
- H-bond donor and acceptor counts
- SMILES string (with copy button)
- InChIKey (with copy button)

### Pocket Evidence

Same predicted pocket and experimental binding data as the residue panel, focused on residues in contact with this ligand.

### External Links

- **PubChem** — Compound page for this ligand
- **DrugBank** — If the ligand is a known drug
- **2D / 3D similarity** — PubChem similarity search by structure

### Ligand Similarity

Tanimoto similarity scores against other ligands present in the structure (CACTVS 881-bit fingerprints). Listed by score, descending. Click any entry to refocus on that ligand.

---

## Export

Click the **Download** button in the viewer header to access export options.

### CSV — Annotation table

A spreadsheet with one row per residue. Columns include:
- UniProt position and amino acid
- PDB residue number (mapped per chain for multi-chain structures)
- One column per PTM category (presence/absence)
- One column per disease label (presence/absence)
- Variant counts (total and pathogenic)
- gnomAD maximum allele frequency
- AlphaMissense per-residue statistics (mean, max, number of scored substitutions)
- Exploratory algorithm results (hotspot tier, contact-hub tier, burden flag, constraint-pocket score)
- PTM–variant proximity statistics
- Nearby ligand CCD codes

**Download CSV + ProtVar predictions** — Same as above but also fetches per-residue ProtVar data (conservation, EVE, ESM-1b, FoldX, CADD) for every residue. Slower because it queries the ProtVar API for the whole protein.

### PyMOL Session

Downloads a `.pml` script. Open it in PyMOL (`File → Run Script`) to reproduce the exact on-screen view:
- Protein cartoon with the same color mode active in the viewer
- Annotation spheres (VDW) with the same colors and filter state
- Ligands shown as colored sticks
- If a residue was focused: that residue and its 5 Å neighborhood shown in ball-and-stick

### VMD Session

Downloads a `.vmd` script. Source it in VMD (`source session.vmd` in the Tk Console) for the same result:
- NewCartoon representation with per-residue colors
- Annotation VDW spheres
- Ligands as Licorice
- Focus state reproduced

### Copy to Clipboard

Copies a residue selection string (for the currently visible/focused residues) in either PyMOL or VMD syntax. The format is set in [Settings](#settings).

---

## Settings

Open the settings page from your browser's extension management menu.

| Setting | Options | Default |
|---|---|---|
| Default structure | AlphaFold first / Experimental first / Best coverage | AlphaFold first |
| Default color mode | Default / pLDDT / B-factor / AlphaMissense | Default |
| Copy-to-clipboard format | PyMOL / VMD | PyMOL |
| PTM search radius | 2–30 Å | 8 Å |
| Pathogenic variant search radius | 2–30 Å | 12 Å |
| Font size | Smaller / Default / Larger | Default |
| Show optional conservation / frequency tracks | On / Off | Off |
| Show exploratory algorithms | On / Off | On |

**Show exploratory algorithms** — When turned off, the four investigational color modes (hotspots, contact hubs, burden residues, constraint pockets) are hidden from the color mode dropdown. Turn this off if you want a cleaner interface focused on direct UniProt annotations.
