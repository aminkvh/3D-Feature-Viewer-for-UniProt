# UniProt 3D Feature Viewer

A Chrome extension that injects an interactive 3D protein structure viewer directly into UniProt entry and variant-viewer pages.

![UniProt 3D Viewer Screenshot](/icons/Screenshot.png) 

## Motivation

When starting a new protein project, I usually want a quick way to see where mutations and PTMs are located in 3D space. UniProt already has great annotations, but the default feature viewer can feel limited and inconsistent for structural exploration. I found myself repeatedly jumping between UniProt and external visualization tools, so I built a small extension to make that workflow easier and faster directly inside the UniProt page.

## Features

- **Seamless Integration**: Automatically injects a "View in 3D" button into the "PTM / Processing" and "Disease & Variants" sections of UniProt.
- **Large-Scale Data**: Pulls in high-throughput proteomics data (PRIDE / PTMeXchange) to visualize modifications that aren't in the standard UniProt features track.
- **Interactive Filtering**: Filter variants by disease, consequence, or provenance using an interactive side panel.
- **Copy IDs**: Quickly copy all currently visualized residue positions to your clipboard for downstream analysis.
- 
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
