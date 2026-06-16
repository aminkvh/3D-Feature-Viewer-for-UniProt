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
ufv_load P35498                  # download AlphaFold model + fetch (loads a plain grey cartoon)
ufv_gui                          # Qt control panel — pick layers, filter variants, read counts
```

`ufv_load` no longer dumps every layer at once; you choose what to show, from the panel or
commands, so the view stays legible. Each layer is independently toggleable:

```
# annotate something you loaded yourself:
load mytraj.pdb, traj
ufv_fetch P35498
ufv_map traj, P35498, identity            # resi == UniProt position
ufv_map traj, P35498, sifts, 7dtd         # map via PDBe/SIFTS for PDB 7DTD
ufv_chain traj, A, 200, 5, 480            # chain A resi 5 == UniProt 200, valid 5..480

ufv_ptms traj
ufv_variants traj                          # defaults to pathogenic only (avoids flooding)
ufv_variants traj, P35498, "pathogenic deleterious"   # combine consequences, or "all"
ufv_sites traj
# cartoon colouring (one at a time): structure + analysis modes
ufv_domains traj ; ufv_topology traj ; ufv_alphamissense traj
ufv_plddt traj ; ufv_bfactor traj ; ufv_burden traj
ufv_hotspots traj ; ufv_contacthubs traj
ufv_hide traj, variants                    # toggle a single layer off
ufv_clear traj                             # remove everything
ufv_csv P35498_residue_annotations.csv     # per-residue one-hot annotation matrix

# residue report + zoom-in
ufv_report 1378, traj                      # print nearby annotations + distances + AM, and zoom in
ufv_focus traj, 1378 ; ufv_resetview traj
```

The **Qt panel** (`ufv_gui`) is the easiest way in. Each sphere layer (PTMs / variants / sites)
is drawn on its **own lightweight object** (`ufv_<obj>_<tag>` — a copy of just the annotated Cα
atoms, no cartoon), so showing thousands of variant points never rebuilds or repaints the main
structure's cartoon and stays fast. Network work runs off the UI thread (the panel shows a status
and stays responsive).

1. Enter an accession → **Fetch** (reports counts). **Open in UniProt ↗** opens the entry page;
   **AI name (ProtNLM)** shows UniProt's ProtNLM/ProtNLM2 AI-predicted protein name when the entry
   has one (predicted from sequence for otherwise-uncharacterised proteins; curated entries say so).
2. Pick a **structure** — AlphaFold, any experimental PDB chain (best-coverage first), or a
   computed model — and **Load** (or **Load all** to pull every listed structure at once).
   Numbering is automatic (AlphaFold/computed → identity, PDB → SIFTS). `ufv_structures` /
   `ufv_use <key>` do the same from the command line.
3. **Layers**: PTMs / Sites / **Ligands** checkboxes; Variants with per-consequence chips +
   **reviewed** only (off = the full variant-viewer set).
4. **Cartoon colouring** (single selector): Domains, Topology, pLDDT (predicted only), B-factor
   (experimental only), AlphaMissense, Burden, **Hotspots**, **Contact hubs**, **Constraint
   pocket**. Colours match the extension exactly; the structure-dependent analyses run on the
   loaded coordinates off the UI thread (numpy-accelerated — hotspots/pockets ≈0.5 s, contact
   hubs sampled for large chains). Applying a colouring is a single `alter` pass over the atoms
   (no per-residue selection parsing), so it stays fast even on large/AlphaFill models.
5. **Report**: pick an annotation type — PTMs / Variants / Sites / Domains / **Ligands** — and
   browse the **colour-coded list** (filter box — large sets are capped to the first 600 rows for
   responsiveness; narrow with the filter to see the rest). **Click any row to zoom in**: the view
   focuses the residue + its 5 Å neighbourhood as annotation-coloured sticks (cartoon dimmed) while a
   **sectioned, colour-coded detail panel** shows variants (with **collapsible evidence** — ClinVar /
   review status / dbSNP / disease behind a *show evidence* toggle), features, AlphaMissense mean +
   the full substitution grid, nearby ligands, and **clickable** nearby residues split into **≤ 8 Å**
   and **8–12 Å** bands (click one to jump the report there). Tick **ProtVar** to add EMBL-EBI
   ProtVar predictors — per-position **conservation** and **M3D** structural effect, plus
   **EVE / ESM1b / FoldX ΔΔG** resolved to each of the residue's actual variants (each links out to
   ProtVar) — complementing AlphaMissense. Click a **ligand** row to zoom to it and
   see its chemistry (name / formula / SMILES / DrugBank links out) and the most **Tanimoto-similar**
   ligands; when a component occurs several times, a **◂ i/N ▸** stepper circulates the copies (each
   highlighted in 3D with a label). The button row under the list is **Show · Hide** (the
   filtered-row markers), **Pick 3D**, and **Reset view**. Tick **Pick 3D** then click any atom in the
   PyMOL viewer to zoom straight to its residue/ligand. **Align** superposes multiple loaded structures.

Each sphere layer is a separate lightweight object (`ufv_<obj>_<tag>`) so large variant sets stay
fast, and the whole panel is scrollable so it never outgrows the screen. **Export CSV** (or
`ufv_csv`) writes the same per-residue one-hot annotation matrix the browser extension produces
(position/aa, one-hot PTM-category & disease columns, site/mutagenesis flags, variant counts,
gnomAD AF, AlphaMissense summary, burden, and hotspot/contact-hub/constraint-pocket tiers for the
loaded structure). **+ ProtVar** (or `ufv_csv …, protvar=1`) additionally fetches ProtVar
conservation / EVE / ESM1b / M3D for every residue — one request per residue, so it's slow
(progress is shown in the status line).

For a structure/trajectory you loaded yourself, expand **Advanced numbering** to set identity /
SIFTS / manual-per-chain (also `ufv_map` / `ufv_chain` on the command line).

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

`ptms` · `variants` (filter: pathogenic/benign/uncertain/deleterious, or several at once) ·
`sites` · `domains` · `topology` · `alphamissense`. Point features render as Cα spheres
(independent, stackable); the three cartoon colourings are mutually exclusive. Each layer is
toggleable — `ufv_hide <obj>, <layer>` removes one without disturbing the others. Colours and
category logic mirror the browser extension exactly.

## Notes / limitations

- **Trajectories:** representations/colours apply across all frames, so a projected
  annotation persists through the trajectory.
- **SIFTS author numbering:** segment endpoints use SIFTS `author_residue_number` when
  present, else SEQRES `residue_number`. For the rare entries where neither matches the
  ATOM-record numbering, anchor the chains manually (`ufv_chain`).
- `ufv_pymol.py` doubles as a CLI: `python ufv_pymol.py --emit-tcl <uid>` /
  `--emit-sifts <pdb> <uid>` / `--download <uid>` (used internally by the VMD plugin).
