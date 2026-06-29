# 3D Feature Viewer for UniProt [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20274183.svg)](https://doi.org/10.5281/zenodo.20274183)

A browser extension that brings structural protein interpretation directly into UniProt. It adds an interactive 3D workspace to UniProt entry pages so you can inspect post-translational modifications, variants, clinical annotations, prediction scores, and residue-prioritization overlays without leaving UniProt.

![UniProt 3D Viewer Screenshot](/icons/Screenshot.png)

## Features

- View UniProt PTMs, disease-associated variants, ClinVar annotations, and AlphaMissense scores on available 3D structures.
- Switch between AlphaFold models and mapped experimental PDB structures.
- Per-residue predictor table: EVE, ESM1b, FoldX ΔΔG, conservation, and CADD scores via ProtVar.
- Binding pocket analysis: PDBe-KB known sites, constraint pocket detection, and pocket confidence scoring.
- Open Targets tractability and drug evidence per residue.
- Color overlays for pLDDT, AlphaMissense, enrichment hotspots, long-range contact hubs, and mutation burden.
- Ligand similarity by CACTVS/Tanimoto fingerprint against AlphaFill transplants.
- Export residue sets and sessions for PyMOL, VMD, spreadsheets, or ML pipelines.

## Install

### From a GitHub Release (recommended)

1. Go to the [Releases page](../../releases) and download the latest release.
   - `chrome-extension-v2.0.0.zip` for Chrome
   - `firefox-extension-v2.0.0.zip` for Firefox
2. Unzip the downloaded file.

**Chrome:**
1. Open `chrome://extensions/`
2. Enable **Developer mode** (top-right toggle)
3. Click **Load unpacked** and select the unzipped `chrome-extension-v2.0.0/` folder

**Firefox:**
1. Open `about:debugging#/runtime/this-firefox`
2. Click **Load Temporary Add-on**
3. Select the `manifest.json` file inside the unzipped `firefox-extension-v2.0.0/` folder

> Firefox temporary add-ons are removed on browser restart. For persistent installation, submit to [Mozilla Add-ons](https://addons.mozilla.org) or use a signed `.xpi`.

### Build from source

```powershell
git clone https://github.com/aminkvh/3D-Feature-Viewer-for-UniProt.git
cd 3D-Feature-Viewer-for-UniProt
pwsh ./build-all.ps1
```

This produces `chrome-build/` and `firefox-build/` ready to load as unpacked extensions, plus `.zip` files for distribution.

## Quick Start

1. Open any UniProt entry, e.g. `https://www.uniprot.org/uniprotkb/P14867/entry`
2. Go to a section with PTMs, variants, or related annotations.
3. Click **View in 3D**.
4. Choose a structure, filter annotations, click residues or the Nearby panel to inspect.
5. Export results via the Download menu.

## Changelog

### v2.0.0
- **3D viewer migrated to Mol\*** — WebGL2 rendering via a sandboxed iframe; smoother interaction on large and multi-chain structures
- **Per-residue predictor table** — EVE, ESM1b, FoldX ΔΔG, conservation, and CADD scores sourced from ProtVar at the per-substitution level
- **Binding pocket analysis** — PDBe-KB known binding sites, constraint pocket heuristic (UFVPocket), and mean pocket pLDDT confidence
- **Open Targets tractability** — drug tractability and clinical evidence per protein, shown in the residue report
- **Ligand panel improvements** — AlphaFill transplant metadata, Tanimoto similarity display, and direct PubChem links
- **PyMOL and VMD session export** with full feature parity (all layers, color modes, hotspots, pockets)
- **Chrome (MV3) + Firefox (MV2) build** — separate manifests with a single `build-all.ps1` producing both packages

### v1.7.3
- Functional Features and Family & Domains windows
- Streaming structure load: AlphaFold model loads immediately, experimental/isoform/computed discovery streams behind it
- Full isoform AlphaFold support with correct annotation mapping
- PyMOL/VMD session export added to Download menu
- Ligand similarity by CACTVS fingerprint (Tanimoto vs PubChem)

## Residue-Prioritization Algorithms

See [METHODS.md](METHODS.md) for a detailed description of all four algorithms: pathogenic-variant enrichment hotspots, long-range contact hubs, mutation/phenotype burden, and AlphaMissense residue scores.

## Data and Privacy

The extension runs entirely in the browser with no server-side component. It retrieves data from public resources (UniProt, PDBe, AlphaFold DB, ProtVar, Open Targets, PubChem, AlphaFill). No personal data is collected. See [PRIVACY_POLICY.md](PRIVACY_POLICY.md) for details.

## License

Released under the MIT License. See [LICENSE](LICENSE).
