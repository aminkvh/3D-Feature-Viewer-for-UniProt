# Migration plan: 3Dmol.js → Mol\* (browser extension viewer)

Status: **proposal — not started.** This documents what a switch would involve so we can decide
before touching code.

## Why consider it

- **Large structures.** Mol\* (molstar) uses GPU-instanced rendering and streaming; it handles
  100k+ atom assemblies and multi-chain complexes far better than 3Dmol.js, which rebuilds
  geometry on the main thread and stutters on big AlphaFill / experimental complexes.
- **Built-in representations** we currently hand-roll (cartoon by property, surfaces, measurement
  labels, sequence panel) are first-class in Mol\*.
- **Maintained & standard.** It's what PDBe and RCSB ship; long-term it's the safer bet.

## Why it's not a drop-in

The extension is a **Manifest V3 content script injected into every UniProt entry page**
(`manifest.json` → `content_scripts`). That makes the cost different from a standalone web app:

- **Bundle size / load.** `lib/3Dmol-min.js` is ~1 MB and loads on every entry page. Mol\* is
  several MB (WebGL2 + workers). Injecting that into UniProt's page on `document_idle` will hurt
  page load unless we **lazy-load it only when the viewer modal opens** (dynamic import from
  `web_accessible_resources`, not a static content-script entry).
- **CSP / workers.** Mol\* spins up web workers and WASM; both need entries in
  `web_accessible_resources` and may collide with UniProt's page CSP. Needs a spike to confirm it
  runs injected (vs. in an iframe we control).
- **Shadow DOM / styling.** The modal renders into the page; Mol\*'s canvas + UI must be scoped so
  UniProt's CSS doesn't leak in (and vice-versa).

## What changes in our code

All viewer access goes through `viewer.js` (`StructureViewer`, ~1300 lines). That's the boundary —
keep its public API identical and rewrite the internals. Consumers that must keep working unchanged:

- **modal.js** — the bulk of the calls (layers, focus, click/hover, screenshot).
- **analysis.js** — `residueGeometry()`, `mappedResidues()` (Cα coords feed hotspot/hub/pocket).
- **export.js** — `currentPdbText`, `getSceneState()` (B-factor PDB rewrite, VMD/PyMOL session export).

### API surface to re-implement (3Dmol → Mol\*)

| `StructureViewer` method | 3Dmol today | Mol\* equivalent |
|---|---|---|
| `init` / `loadStructure` | `createViewer`, `addModel(pdb/cif)` | `createPluginUI` / `builders.structure` from URL or string |
| `applyCartoonColoring` | `setStyle({cartoon:{colorfunc}})` | custom `ColorTheme` keyed by residue → value |
| `showPTMs/Variants/sites` (spheres) | `addStyle({sphere})` per residue | overpaint representation or shape primitives |
| `showLigands` / `focusLigand` | `addStyle`, `zoomTo` | `Structure` selection + `setSubtreeVisibility`, `camera.focus` |
| `focusResidue` | `zoomTo(sel)` + sticks | `camera.focus` on `StructureSelection` |
| `showProximityLines` | `addLine` shapes | `Shape`/`Mesh` representation |
| `residueGeometry` / `mappedResidues` | iterate model atoms | `StructureProperties` over `StructureElement` loop |
| `currentPdbText` | cached source text | keep the fetched text as-is (independent of viewer) |
| `screenshot` | `pngURI()` | `plugin.helpers.viewportScreenshot` |
| click / hover / dblclick | `setClickable/Hoverable` | `plugin.behaviors.interaction` event subscriptions |

The hardest items are **per-residue colouring** (we map UniProt position → colour; Mol\* wants a
`ColorTheme` provider) and **annotation spheres / proximity lines** (custom shapes), because those
are where our extension differs most from a stock structure viewer.

## Phased approach (each phase shippable / revertable)

1. **Spike (½–1 day):** prove Mol\* loads lazily inside the injected modal on a real UniProt page —
   CSP, workers, Shadow-DOM scoping, bundle from `web_accessible_resources`. *Go / no-go gate.*
2. **Parallel viewer behind a flag:** new `viewer-molstar.js` implementing the same `StructureViewer`
   contract; an options toggle picks the backend. 3Dmol stays default.
3. **Port representations:** cartoon colouring → custom `ColorTheme`; spheres/lines → shapes; ligands.
4. **Port interaction + geometry:** click/hover/focus, `residueGeometry`, `getSceneState`, screenshot.
5. **Validate** against the existing structures we test (AlphaFold, big PDB complex, AlphaFill,
   computed model) for parity, then flip the default and remove 3Dmol.

## Effort & risk

- **Effort:** ~3–5 focused days *after* a successful spike; most of it is phase 3–4.
- **Risk:** medium-high. The spike (phase 1) is the real risk — if Mol\* can't be injected cleanly
  into UniProt's page, we'd need an iframe-hosted viewer, which is a bigger architectural change.
- **Reversible:** yes, as long as we keep the `StructureViewer` contract and ship behind a flag.

## Recommendation

Do **phase 1 (the spike) only**, then reassess. Don't start the full port until we've confirmed
Mol\* runs injected on a live UniProt page within MV3's constraints. If the spike fails, the cheaper
win is to optimize the current 3Dmol path (lazy geometry, fewer full re-renders) rather than migrate.
