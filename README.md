# UniProt 3D Feature Viewer

A Chrome extension that injects an interactive 3D protein structure viewer directly into UniProt entry and variant-viewer pages.

![UniProt 3D Viewer Screenshot](icon/Screenshot.png) 

## Motivation

When working with proteins, understanding the structural context of Post-Translational Modifications (PTMs) and genetic variants is crucial. A simple 1D sequence position often doesn't tell the whole story.

*   *Is this mutation buried in the hydrophobic core, or exposed on the surface?*
*   *Is this phosphorylation site accessible to kinases?*
*   *Does this disease variant cluster with other known mutations in 3D space?*

Currently, researchers have to manually cross-reference 1D sequence positions from UniProt with a separate 3D viewer (like PyMOL or Chimera). This extension solves that problem by providing a zero-setup, on-the-fly 3D viewer right where you need it—directly on the UniProt page.

## Features

- **Seamless Integration**: Automatically injects a "View in 3D" button into the "PTM / Processing" and "Disease & Variants" sections of UniProt.
- **SPA Aware**: Works perfectly with UniProt's single-page application navigation.
- **Large-Scale Data**: Pulls in high-throughput proteomics data (PRIDE / PTMeXchange) to visualize modifications that aren't in the standard UniProt features track.
- **Interactive Filtering**: Filter variants by disease, consequence, or provenance using an interactive side panel.
- **Visual Focus**: Click on any variant or PTM to automatically zoom in, highlight the residue in ball-and-stick mode, and display local interactions (within 5Å).
- **Copy IDs**: Quickly copy all currently visualized residue positions to your clipboard for downstream analysis.
- **Secure**: Built with strict XSS prevention and size guards to handle large AlphaFold models safely.

## Installation

### For Users (Chrome Web Store)
*(Link coming soon once published to the Chrome Web Store!)*

### For Developers (Manual Install)
1. Clone this repository or download the source code.
2. Open Google Chrome and navigate to `chrome://extensions/`.
3. Enable **Developer mode** in the top right corner.
4. Click **Load unpacked** and select the folder containing the extension files.

## Usage
1. Navigate to any UniProt protein entry (e.g., `https://www.uniprot.org/uniprotkb/P14867/entry`).
2. Scroll down to the **Disease & Variants** or **PTM / Processing** sections.
3. Click the **"View in 3D"** button.
4. Use the side panel to filter the displayed features. Click on any colored sphere on the 3D model to see detailed information.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
