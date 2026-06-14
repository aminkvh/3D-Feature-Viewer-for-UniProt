# 3D Feature Viewer for UniProt — desktop plugins (PyMOL & VMD)

Desktop counterparts of the browser extension. They bring UniProt residue-level
annotations onto a structure **or trajectory** already open in PyMOL/VMD — or download
the AlphaFold model directly — and project PTMs, disease variants, ClinVar /
AlphaMissense, functional sites, domains/regions, and membrane topology onto it.

Same public data sources as the extension (UniProtKB, EBI Proteins API, AlphaFold DB,
PDBe/SIFTS). No API keys.

| File | Role |
|---|---|
| `ufv_pymol.py` | PyMOL plugin **and** the shared data backend (standard library only) |
| `ufv_vmd.tcl`  | VMD plugin (native `atomselect`/representations; calls `ufv_pymol.py` to fetch) |

The VMD plugin shells out to `ufv_pymol.py` for fetching/mapping (HTTPS + JSON), so it
needs **python3 on PATH**. PyMOL uses its own bundled Python and needs nothing extra.

---

## PyMOL

**Install:** *Plugin* → *Plugin Manager* → *Install New Plugin* → *Choose file…* →
`ufv_pymol.py`. Or `run /path/to/ufv_pymol.py`.

```
ufv_load P35498                  # download AlphaFold model + annotate
ufv_gui                          # graphical panel

# annotate something you loaded yourself:
load mytraj.pdb, traj
ufv_fetch P35498
ufv_map traj, P35498, identity            # resi == UniProt position
ufv_map traj, P35498, sifts, 7dtd         # map via PDBe/SIFTS for PDB 7DTD
ufv_chain traj, A, 200, 5, 480            # chain A resi 5 == UniProt 200, valid 5..480
ufv_ptms traj ; ufv_variants traj ; ufv_sites traj
ufv_domains traj ; ufv_topology traj ; ufv_alphamissense traj
ufv_clear traj
```

## VMD

**Install:** put both files in a folder and `source ufv_vmd.tcl` from the VMD console
(or drop into a VMD plugin directory). If your Python isn't called `python`, set it once:
`ufv_python python3`. If `ufv_pymol.py` lives elsewhere: `ufv_backend /path/ufv_pymol.py`.

```
ufv_load P35498                  ;# download AlphaFold model + annotate
ufv_gui                          ;# graphical panel

# annotate something you loaded yourself (molid defaults to top):
mol new mytraj.pdb
ufv_fetch P35498
ufv_map identity                 ;# resid == UniProt position
ufv_map sifts 7dtd               ;# map via PDBe/SIFTS for PDB 7DTD
ufv_chain A 200 5 480            ;# chain A resid 5 == UniProt 200, valid 5..480
ufv_ptms ; ufv_variants ; ufv_variants pathogenic ; ufv_sites
ufv_domains ; ufv_topology ; ufv_alphamissense
ufv_clear
```

---

## Residue numbering (both)

A loaded object's residue numbers rarely equal UniProt positions. Choose one:

- **Identity** — `resi == UniProt position` (AlphaFold and UniProt-numbered models).
- **SIFTS** — map through PDBe/SIFTS for a given PDB id (handles gaps, author numbering,
  cleaved/renumbered chains, e.g. insulin 3i40: chain A 90–110→1–21, chain B 25–53→1–29).
- **Manual** — anchor each chain by hand: residue `resiStart` corresponds to UniProt
  `uStart`, linearly, optionally bounded by `resiEnd`. Ideal for **trajectories** with
  custom numbering or to define each chain's reference point / start–finish.

## Layers (both)

`ptms` · `variants` (filter: pathogenic/benign/uncertain/deleterious) · `sites` ·
`domains` · `topology` · `alphamissense`. Point features render as Cα spheres; ranges
colour the cartoon. Colours and category logic mirror the browser extension exactly.

## Notes / limitations

- **Trajectories:** representations/colours apply across all frames, so a projected
  annotation persists through the trajectory.
- **SIFTS author numbering:** segment endpoints use SIFTS `author_residue_number` when
  present, else SEQRES `residue_number`. For the rare entries where neither matches the
  ATOM-record numbering, anchor the chains manually (`ufv_chain`).
- `ufv_pymol.py` doubles as a CLI: `python ufv_pymol.py --emit-tcl <uid>` /
  `--emit-sifts <pdb> <uid>` / `--download <uid>` (used internally by the VMD plugin).
