# 3D Feature Viewer for UniProt — PyMOL & VMD plugins

Project UniProt residue annotations — PTMs, disease variants, AlphaMissense, functional
sites, domains, and membrane topology — onto a structure in PyMOL or VMD. Loads AlphaFold,
experimental PDB, or computed models directly, or annotates anything you already have open.

Same public data as the [browser extension](../README.md) (UniProtKB, EBI Proteins API,
AlphaFold DB, PDBe/SIFTS, EMBL-EBI ProtVar). No API keys.

| File | Role |
|---|---|
| `ufv_pymol.py` | PyMOL plugin **and** the shared data backend (standard library only) |
| `ufv_vmd.tcl`  | VMD plugin (calls `ufv_pymol.py` to fetch — needs **python3 on PATH**) |

## Install

**PyMOL** — *Plugin → Plugin Manager → Install New Plugin → Choose file…* → `ufv_pymol.py`.

**VMD** — put both files in one folder, then `source ufv_vmd.tcl` in the VMD console.
If your interpreter isn't `python`, set it once: `ufv_python python3`.

## Use

Open the panel and drive everything from there:

```
ufv_gui
```

Enter an accession → **Fetch** → pick a structure → **Load**. Then toggle annotation layers,
colour the cartoon (domains, AlphaMissense, hotspots, …), and click any residue to zoom in
and read its report. **Export CSV** writes a per-residue annotation matrix.

Or script it:

```
ufv_load P35498            # download AlphaFold model + annotate
ufv_ptms                   # PTM spheres        (also: ufv_variants, ufv_sites, ufv_ligands)
ufv_alphamissense          # colour cartoon     (also: ufv_domains, ufv_hotspots, ufv_pockets…)
ufv_report 1378            # residue report + zoom in
ufv_csv out.csv            # per-residue matrix
ufv_clear
```

## Residue numbering

Loaded residue numbers rarely equal UniProt positions. Pick one (panel → *Numbering*, or
`ufv_map` / `ufv_chain`):

- **Identity** — `resi == UniProt position` (AlphaFold / UniProt-numbered models).
- **SIFTS** — map through PDBe/SIFTS for a PDB id (handles gaps and author numbering).
- **Manual** — anchor each chain by hand; ideal for trajectories with custom numbering.

---

`ufv_pymol.py` also runs as a CLI used internally by the VMD plugin
(`python ufv_pymol.py --emit-tcl <uid>`, `--download <uid>`, …).
