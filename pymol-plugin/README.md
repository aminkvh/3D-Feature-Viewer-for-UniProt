# 3D Feature Viewer for UniProt — PyMOL plugin

The desktop counterpart of the browser extension. It brings UniProt residue-level
annotations onto a structure **or trajectory** already open in PyMOL — or downloads
the AlphaFold model directly — and projects PTMs, disease variants, ClinVar /
AlphaMissense, functional sites, domains/regions, and membrane topology onto it.

Same public data sources as the extension (UniProtKB, EBI Proteins API, AlphaFold DB,
PDBe/SIFTS). No API keys, no third-party Python packages — standard library only.

## Install

**Plugin Manager (recommended):** PyMOL → *Plugin* → *Plugin Manager* → *Install New
Plugin* → *Choose file…* → select `ufv_pymol.py`. A menu item *3D Feature Viewer for
UniProt* appears under *Plugin*.

**Or just run it:**

```
run /path/to/ufv_pymol.py
```

## Usage

### Download a model and annotate it

```
ufv_load P35498        # download the AlphaFold model, then project PTMs + variants + sites
ufv_domains P35498     # add domain/region colouring
ufv_alphamissense P35498
```

### Annotate a structure / trajectory you already loaded

```
load mytraj.pdb, traj
ufv_fetch P35498                       # fetch annotations once (cached)
```

Tell the plugin how `traj`'s residue numbers relate to UniProt positions — pick one:

```
ufv_map traj, P35498, identity         # resi == UniProt position (AlphaFold-style)
ufv_map traj, P35498, sifts, 7dtd      # map through PDBe/SIFTS for PDB 7DTD
```

…or define each chain by hand (ideal for trajectories with custom numbering, or to set
each chain's reference point / start–finish):

```
ufv_chain traj, A, 1                   # chain A resi 1   == UniProt 1
ufv_chain traj, A, 200, 5, 480         # chain A resi 5   == UniProt 200, valid resi 5..480
ufv_chain traj, B, 1, 1, 350           # chain B resi 1   == UniProt 1, valid resi 1..350
```

Then project any layer (object name optional — defaults to the active object):

```
ufv_ptms traj
ufv_variants traj                      # all variants, coloured by clinical consequence
ufv_variants traj, P35498, pathogenic  # filter: pathogenic | benign | uncertain | deleterious
ufv_sites traj
ufv_domains traj
ufv_topology traj
ufv_alphamissense traj
ufv_clear traj
```

### GUI

```
ufv_gui
```

A panel for accession, target object, numbering (Identity / SIFTS / Manual chain), and
one-click layer buttons.

## Commands

| Command | Purpose |
|---|---|
| `ufv_load uid [, name]` | Download the AlphaFold model and project PTMs + variants + sites |
| `ufv_fetch uid` | Fetch & cache all annotation layers for an accession |
| `ufv_map obj, uid, mode [, pdb]` | Numbering: `identity` or `sifts, <pdbid>` |
| `ufv_chain obj, chain, uStart [, resiStart [, resiEnd]]` | Manually anchor one chain's numbering |
| `ufv_ptms` / `ufv_variants` / `ufv_sites` | Point annotations as Cα spheres |
| `ufv_domains` / `ufv_topology` | Range annotations as cartoon colouring |
| `ufv_alphamissense` | Colour by mean AlphaMissense pathogenicity |
| `ufv_clear [obj]` | Remove UFV overlays |
| `ufv_info [uid]` | Print an annotation summary |

## Notes / limitations

- **Trajectories:** colouring and representations apply across all states, so a projected
  annotation persists through every frame.
- **SIFTS author numbering:** segment endpoints use SIFTS `author_residue_number` when
  present, otherwise the SEQRES `residue_number`. For the rare entries where neither matches
  the ATOM-record numbering, anchor the chains manually with `ufv_chain`.
- Colours and category logic mirror the browser extension exactly.
