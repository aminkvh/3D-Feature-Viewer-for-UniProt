# 3D Feature Viewer for UniProt [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20274183.svg)](https://doi.org/10.5281/zenodo.20274183)

3D Feature Viewer for UniProt is a browser extension that brings structural protein interpretation directly into UniProt. It adds an interactive 3D workspace to UniProt entry and variant pages so users can inspect post-translational modifications, variants, clinical annotations, prediction scores, and residue-prioritization overlays without leaving the UniProt workflow.

![UniProt 3D Viewer Screenshot](/icons/Screenshot.png)

## What It Helps You Do

- View UniProt PTMs, disease-associated variants, ClinVar annotations, and AlphaMissense scores on available 3D structures.
- Switch between AlphaFold models and mapped experimental PDB structures when available.
- Inspect residues in local 3D context, including nearby annotated residues.
- Use color overlays for pLDDT, beta-factor, AlphaMissense, enrichment hotspots, long-range contact hubs, and mutation/phenotype burden.
- Export residue sets and results for follow-up work in PyMOL, VMD, PDB-based workflows, spreadsheets, or machine-learning pipelines.

## Why Use It

UniProt is often the first stop for protein annotation, but structural interpretation usually requires moving between multiple sites, choosing a structure, matching residue numbers, and manually transferring features into a molecular viewer. This extension keeps that process on the UniProt page and handles residue mapping for AlphaFold and chain-specific PDB structures.

The built-in residue-prioritization views are intended for exploration and hypothesis generation. They highlight regions that may deserve closer inspection; they are not clinical classifiers or validated predictors of pathogenicity.

## Install

### Chrome Web Store

[3D Feature Viewer for UniProt](https://chromewebstore.google.com/detail/uniprot-3d-feature-viewer/fplpkbigppbpbcdmpkefdoilfgcicaof?authuser=1)

### Manual Install

1. Download or clone this repository.
2. Open your browser's extensions page.
   - Chrome: `chrome://extensions/`
   - Firefox: `about:debugging#/runtime/this-firefox`
3. Enable developer mode or temporary add-on loading.
4. Load this project folder as an unpacked extension.

## Quick Start

1. Open a UniProt entry page, such as `https://www.uniprot.org/uniprotkb/P14867/entry`.
2. Go to a section with PTMs, variants, or related residue annotations.
3. Click **View in 3D**.
4. Choose a structure, filter annotations, select residues, and export any results you want to reuse.

## What's new in v1.5.0

* You can now visualize calculated features directly on AlphaFold or PDB structures.
* New structure-level annotations include beta values, pLDDT, hotspot enrichment, long-range contact hubs, mutation/phenotype burden scores, and AlphaMissense scores.
* Calculated features can now be downloaded as a one-hot encoded CSV file.
* A new interactive sequence ribbon shows annotations directly on the protein sequence.
* Structure visualization is smoother and more stable, with improved support for different chains and complex structures.
* The extension has been renamed to **3D Feature Viewer for UniProt**.

## Data and Privacy

The extension runs in the browser and has no server-side component. It retrieves protein, structure, and annotation data from public resources including UniProt, PDBe/PDB, AlphaFold DB, and related public endpoints.

## License

Released under the MIT License. See [LICENSE](LICENSE).
