# Privacy Policy

> Last updated: June 2026

This extension does **not** collect, store, or share any personal data of any kind.

---

## How It Works

The extension operates entirely within your browser. When you open a UniProt entry page, it reads the page URL to identify the protein accession, then fetches publicly available scientific data from external APIs to populate the 3D viewer and annotation panels. No information about you, your searches, or your browser is included in any of these requests.

---

## External Services Accessed

All network requests are read-only calls to public scientific APIs. The requests contain only the protein accession or structure identifier needed to retrieve data. No user identifiers or browsing history are transmitted beyond the network metadata inherent to any web request.

| Service | Data Retrieved | Privacy Information |
|---|---|---|
| [UniProt](https://www.uniprot.org) | Protein sequence, variants, PTMs, functional features | [uniprot.org/help/privacy](https://www.uniprot.org/help/privacy) |
| [PDBe / EBI](https://www.ebi.ac.uk) | Experimental structures, residue mapping, binding sites, PDBe-KB annotations, 3D-Beacons structure index | [ebi.ac.uk/data-protection](https://www.ebi.ac.uk/data-protection) |
| [AlphaFold DB](https://alphafold.ebi.ac.uk) | Predicted protein structures and AlphaMissense scores | [ebi.ac.uk/data-protection](https://www.ebi.ac.uk/data-protection) |
| [AlphaFill](https://alphafill.eu) | Transplanted ligand metadata for AlphaFold models | [alphafill.eu](https://alphafill.eu) |
| [ProtVar](https://www.ebi.ac.uk/ProtVar) | Per-substitution effect predictions (EVE, ESM1b, FoldX, CADD, conservation) | [ebi.ac.uk/data-protection](https://www.ebi.ac.uk/data-protection) |
| [Open Targets Platform](https://platform.opentargets.org) | Drug tractability and clinical evidence | [platform.opentargets.org/privacy-policy](https://platform.opentargets.org/privacy-policy) |
| [RCSB PDB](https://www.rcsb.org) | Structure coordinate files and structural motif search | [rcsb.org/pages/privacy-policy](https://www.rcsb.org/pages/privacy-policy) |
| [PubChem](https://pubchem.ncbi.nlm.nih.gov) | Ligand chemical descriptors and fingerprints | [ncbi.nlm.nih.gov/home/about/policies](https://www.ncbi.nlm.nih.gov/home/about/policies) |
| [SWISS-MODEL](https://swissmodel.expasy.org) | Comparative homology models | [expasy.org/terms-of-use](https://www.expasy.org/terms-of-use) |
| [ModelArchive](https://www.modelarchive.org) | Deposited computed models | [modelarchive.org](https://www.modelarchive.org) |

---

## Local Data Storage

The extension uses `chrome.storage.local` (or `localStorage` as a fallback) exclusively to store your **settings** (default structure preference, color mode, search radii, font scale). This data:

- Never leaves your device
- Contains no personal information
- Can be cleared at any time via your browser's extension management page

No session data, viewed proteins, residue selections, or analysis results are persisted between sessions.

---

## What This Extension Does NOT Do

- Collect any personal or identifying information
- Track which proteins you view or which residues you click
- Store data beyond your settings preferences
- Share, sell, or transmit any data to third parties
- Use cookies, analytics, or telemetry of any kind
- Make requests to any server operated by this extension (there is no backend)

---

## Browser Permissions

The extension requests the following permissions:

| Permission | Why |
|---|---|
| `storage` | Save your settings (structure preference, color mode, etc.) locally |
| Host permissions for listed scientific APIs | Fetch protein and structure data from those services |

No other permissions are requested or used.

---

## Contact

If you have questions or concerns, please open an issue on [GitHub](https://github.com/aminkvh/3D-Feature-Viewer-for-UniProt/issues).
