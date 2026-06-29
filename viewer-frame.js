/*
 * viewer-frame.js — Mol* driver, runs inside the sandboxed viewer-frame.html iframe.
 *
 * Driven by viewer-molstar.js (content script) over postMessage. Sandboxed (opaque origin, no
 * chrome.*) because Mol*'s WASM glue needs 'unsafe-eval' (only allowed in a sandbox CSP).
 *
 * Protocol
 *   parent → frame:  { ns:'ufv', cmd, reqId, ...args }
 *   frame  → parent: { ns:'ufv', evt:'ready' }                      once Mol* is created
 *                    { ns:'ufv', evt:'atoms', atoms:[...] }         after each structure load
 *                    { ns:'ufv', evt:'result', reqId, ok, value }   per-command completion
 *
 * Commands: loadData, loadUrl, clear, reset, resize, screenshot, background, setCartoon.
 * (M2b/M2c add: spheres, focus, picking.)
 */
(function () {
  'use strict';

  const NS = 'ufv';
  const L = () => window.molstar.lib;
  const statusEl = document.getElementById('ufv-spike-status');
  const setStatus = (msg, isError) => {
    if (!statusEl) return;
    statusEl.style.display = msg ? 'block' : 'none';
    statusEl.textContent = msg || '';
    statusEl.style.color = isError ? '#ff8a80' : '#cfe8ff';
  };
  const post = (msg) => { try { parent.postMessage(Object.assign({ ns: NS }, msg), '*'); } catch (_) {} };
  const reply = (reqId, ok, value) => post({ evt: 'result', reqId, ok, value });
  const fail = (reqId, message) => { if (reqId != null) reply(reqId, false, { error: String(message) }); else post({ evt: 'error', message: String(message) }); };

  const AF_API = (acc) => `https://alphafold.ebi.ac.uk/api/prediction/${acc}`;

  let viewer = null;     // molstar Viewer wrapper
  let ready = false;
  let overpaintRefs = []; // state refs of the overpaint nodes we add, so we can replace them
  let ligandTransparencyRefs = []; // transparency nodes hiding individually-toggled ligands
  let _lastLigandHidden = [];      // last hidden-ligand list, re-applied after a focus rep is created
  let markerByColor = new Map(); // colorInt -> {sel, rep} refs; persistent so toggles update in place
  let markerRepRefs = new Set(); // rep refs of all spacefill markers — excluded from cartoon overpaint
  let initialRepRefs = new Set(); // reps present right after structure load (cartoon backbone) — focus sticks appear later and are NOT in this set
  let focusSelRef = null, focusRepRef = null; // our own focus stick component + rep (replaces Mol* focus manager)
  let _focusGen = 0; // bumped per focus so a superseded focus's deferred work (interactions rep) is skipped
  let _pendingInter = null; // { ref, gen } — interactions rep to add AFTER the camera animation settles

  // Add the non-covalent interactions rep to a focus selection, unless a newer focus superseded it.
  function addFocusInteractions(ref, gen) {
    if (gen !== _focusGen) return;
    try { viewer.plugin.builders.structure.representation.addRepresentation(ref, { type: 'interactions' }).catch(() => {}); } catch (_) {}
  }
  let _focusedResidue = null;    // { chain, resi } of selected residue in focus — excluded from annotation in sticks
  let _focusAnnotations = null;  // [{ chain, resi, color }] from sphere buffer, used to color focus sticks
  let labelMap = {};      // 'chain|resi' -> annotation text, shown in Mol*'s native hover label
  let _lastCartoonPayload = null; // cached so focus representation inherits annotation overpaint
  let _bgColor = 0x0c111b; // stored so it can be re-applied after each structure load

  // ---- colour helpers ------------------------------------------------------------------------
  function colorToInt(c) {
    if (typeof c === 'number') return c;
    if (!c) return 0xd0d0d0;
    let m = /^#?([0-9a-f]{6})$/i.exec(c);
    if (m) return parseInt(m[1], 16);
    m = /rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)/i.exec(c);
    if (m) return (parseInt(m[1]) << 16) | (parseInt(m[2]) << 8) | parseInt(m[3]);
    return 0xd0d0d0;
  }

  // ---- structure access ----------------------------------------------------------------------
  function currentStructureCell() {
    const cur = viewer.plugin.managers.structure.hierarchy.current;
    return cur && cur.structures && cur.structures[0];
  }
  function currentStructureData() {
    const s = currentStructureCell();
    return s && s.cell && s.cell.obj && s.cell.obj.data;
  }
  function allRepresentationRefs() {
    const refs = [];
    const cur = viewer.plugin.managers.structure.hierarchy.current;
    (cur && cur.structures || []).forEach(st => (st.components || []).forEach(c =>
      (c.representations || []).forEach(r => { if (r.cell && r.cell.transform) refs.push(r.cell.transform.ref); })));
    return refs;
  }
  // Structure-wide element bundles (all / polymer / polymer-carbons / het-carbons / empty) depend ONLY on
  // the loaded structure (group_PDB + element), but setCartoon / setLigandVisibility used to rebuild them
  // on EVERY focus — each a full pass over ~25k atoms (AlphaFill). Cache them per structure so a zoom
  // in/out doesn't recompute them. PERF.
  let _sBundles = null;
  function structureBundles() {
    const structure = currentStructureData();
    const cell = currentStructureCell();
    const ref = cell && cell.cell && cell.cell.transform && cell.cell.transform.ref;
    if (_sBundles && _sBundles.ref === ref && structure) return _sBundles;
    const Q = L().structure.Queries, SP = L().structure.StructureProperties;
    const QC = L().structure.QueryContext, Bundle = L().structure.StructureElement.Bundle;
    const isPolymer = (ctx) => SP.residue.group_PDB(ctx.element) !== 'HETATM';
    const isC = (ctx) => SP.atom.type_symbol(ctx.element) === 'C';
    const mk = (q) => Bundle.fromSelection(q(new QC(structure)));
    _sBundles = {
      ref,
      all: mk(Q.generators.all),
      polymer: mk(Q.generators.atoms({ residueTest: isPolymer })),
      polymerCarbon: mk(Q.generators.atoms({ atomTest: isC, residueTest: isPolymer })),
      hetCarbon: mk(Q.generators.atoms({ atomTest: isC, residueTest: (ctx) => !isPolymer(ctx) })),
      empty: mk(Q.generators.atoms({ residueTest: () => false })),
    };
    return _sBundles;
  }
  // Refs of representations that must keep their OWN colouring and never be overpainted/transparency'd:
  // SNFG carbohydrate symbols (glycans like NAG stay blue, not green) and non-covalent `interactions`
  // (dashed H-bond / salt-bridge lines keep their interaction-type colours).
  function carbohydrateRepRefs() {
    const set = new Set();
    const cur = viewer.plugin.managers.structure.hierarchy.current;
    (cur && cur.structures || []).forEach(st => (st.components || []).forEach(c =>
      (c.representations || []).forEach(r => {
        const t = r.cell && r.cell.transform && r.cell.transform.params && r.cell.transform.params.type;
        if (t && (t.name === 'carbohydrate' || t.name === 'interactions') && r.cell.transform) set.add(r.cell.transform.ref);
      })));
    return set;
  }
  // Map each representation ref to its type name ('cartoon', 'ball-and-stick', …) so colouring can be
  // decided by what the rep IS, not by when it was created.
  function repTypeByRef() {
    const m = new Map();
    const cur = viewer.plugin.managers.structure.hierarchy.current;
    (cur && cur.structures || []).forEach(st => (st.components || []).forEach(c =>
      (c.representations || []).forEach(r => {
        const t = r.cell && r.cell.transform && r.cell.transform.params && r.cell.transform.params.type;
        if (r.cell && r.cell.transform) m.set(r.cell.transform.ref, t && t.name);
      })));
    return m;
  }

  // ---- atom extraction (feeds analysis/mapping on the parent) ---------------------------------
  function extractAtoms(structure) {
    const SP = L().structure.StructureProperties;
    const StructureElement = L().structure.StructureElement;
    const loc = StructureElement.Location.create(structure);
    const atoms = [];
    for (const unit of structure.units) {
      loc.unit = unit;
      const els = unit.elements;
      for (let i = 0, n = els.length; i < n; i++) {
        loc.element = els[i];
        atoms.push({
          x: SP.atom.x(loc), y: SP.atom.y(loc), z: SP.atom.z(loc),
          resi: SP.residue.auth_seq_id(loc),
          chain: SP.chain.auth_asym_id(loc),
          resn: SP.atom.auth_comp_id(loc),
          atom: SP.atom.auth_atom_id(loc) || SP.atom.label_atom_id(loc),
          b: SP.atom.B_iso_or_equiv(loc),
          elem: SP.atom.type_symbol(loc),
          hetflag: SP.residue.group_PDB(loc) === 'HETATM',
        });
      }
    }
    return atoms;
  }

  // ---- cartoon colouring via overpaint --------------------------------------------------------
  // groups: [{ color, residues: [[chain, resi], ...] }]  + base colour for everything else.
  async function setCartoon(base, groups) {
    const structure = currentStructureData();
    if (!structure) return false;
    const Q = L().structure.Queries;
    const SP = L().structure.StructureProperties;
    const QueryContext = L().structure.QueryContext;
    const Bundle = L().structure.StructureElement.Bundle;

    // Polymer (protein) vs ligand: ATOM vs HETATM. This is reliable for experimental + most model
    // ligands. (entity.type proved unreliable here — it reported 'polymer' for ligands too, which made
    // ligands inherit annotation/grey instead of green. Branched-glycan edge cases are handled per-rep.)
    const isPolymer = (ctx) => SP.residue.group_PDB(ctx.element) !== 'HETATM';

    // Annotation bundles are POLYMER-only: a ligand whose residue number happens to match an
    // annotated protein position must never inherit the annotation colour.
    const bundleForSet = (resSet) => Bundle.fromSelection(Q.generators.atoms({
      residueTest: (ctx) => isPolymer(ctx)
        && resSet.has(SP.chain.auth_asym_id(ctx.element) + '|' + SP.residue.auth_seq_id(ctx.element)),
    })(new QueryContext(structure)));
    const sb = structureBundles(); // cached per structure (PERF)
    const allBundle = sb.all;
    const polymerBundle = sb.polymer;

    // Cartoon layers (applied to the polymer CARTOON rep): grey base + per-residue annotation groups.
    const layers = [{ bundle: polymerBundle, color: colorToInt(base), clear: false }];
    for (const g of (groups || [])) {
      if (!g.residues || !g.residues.length) continue;
      const set = new Set(g.residues.map(([c, r]) => (c == null ? '' : c) + '|' + r));
      const anyChain = g.residues.some(([c]) => c == null);
      const bundle = anyChain
        ? Bundle.fromSelection(Q.generators.atoms({
            residueTest: (ctx) => {
              if (!isPolymer(ctx)) return false;
              const ch = SP.chain.auth_asym_id(ctx.element), ri = SP.residue.auth_seq_id(ctx.element);
              return set.has(ch + '|' + ri) || set.has('|' + ri);
            },
          })(new QueryContext(structure)))
        : bundleForSet(set);
      layers.push({ bundle, color: colorToInt(g.color), clear: false });
    }

    // Focus/stick layers (applied to BALL-AND-STICK reps — both the ligand rep and the focus sticks):
    //   1) clear → element theme (N=blue, O=red, S=yellow; carbons are chain-id-coloured),
    //   2) polymer carbons → grey CPK baseline,
    //   3) annotation → recolour annotated polymer residues' carbons,
    //   4) ligand (non-polymer) carbons → consistent green.
    // The selected residue is excluded from (3) so its carbons stay grey = plain Mol* default look.
    const focusKey     = _focusedResidue ? (_focusedResidue.chain || '') + '|' + _focusedResidue.resi : null;
    const focusKeyAny  = _focusedResidue ? '|' + _focusedResidue.resi : null;
    const carbonBundle = sb.polymerCarbon;
    const LIGAND_CARBON = 0x33aa33;
    const hetCarbonBundle = sb.hetCarbon;
    const focusLayers = [
      { bundle: allBundle, color: 0, clear: true },           // revert to element theme
      { bundle: carbonBundle, color: 0x909090, clear: false }, // polymer carbons → CPK grey
    ];
    const focusGroups = _focusAnnotations
      ? (() => {
          const byColor = new Map();
          for (const a of _focusAnnotations) {
            if (!byColor.has(a.color)) byColor.set(a.color, []);
            byColor.get(a.color).push([a.chain, a.resi]);
          }
          return [...byColor.entries()].map(([color, residues]) => ({ color, residues }));
        })()
      : (groups || []);
    for (const g of focusGroups) {
      if (!g.residues || !g.residues.length) continue;
      // Exclude exactly the selected residue; keep all annotated neighbors
      const focusResidues = focusKey
        ? g.residues.filter(([c, r]) => {
            const k = (c == null ? '' : c) + '|' + r;
            return k !== focusKey && k !== focusKeyAny;
          })
        : g.residues;
      if (focusResidues.length === 0) continue;
      const focusSet = new Set(focusResidues.map(([c, r]) => (c == null ? '' : c) + '|' + r));
      // Paint ONLY the carbon atoms of annotated residues with the annotation colour, so N/O/S keep
      // their element colours (CPK convention: carbons = category colour, heteroatoms = by element).
      const focusBundle = Bundle.fromSelection(Q.generators.atoms({
        atomTest: (ctx) => SP.atom.type_symbol(ctx.element) === 'C',
        residueTest: (ctx) => {
          if (!isPolymer(ctx)) return false; // never colour ligand carbons with annotation
          const ch = SP.chain.auth_asym_id(ctx.element), ri = SP.residue.auth_seq_id(ctx.element);
          return focusSet.has(ch + '|' + ri) || focusSet.has('|' + ri);
        },
      })(new QueryContext(structure)));
      focusLayers.push({ bundle: focusBundle, color: colorToInt(g.color), clear: false });
    }
    // Ligand carbons → consistent green (after the clear→element revert).
    focusLayers.push({ bundle: hetCarbonBundle, color: LIGAND_CARBON, clear: false });

    const OT = L().plugin.StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle;
    const carbRefs = carbohydrateRepRefs(); // SNFG glycan symbols — keep their default colours
    const typeByRef = repTypeByRef();
    const b = viewer.plugin.build();
    overpaintRefs.forEach(ref => { try { b.delete(ref); } catch (_) {} });
    overpaintRefs = [];
    for (const ref of allRepresentationRefs().filter(r => !markerRepRefs.has(r) && !carbRefs.has(r))) {
      // Decide by REPRESENTATION TYPE, not load-time timing: the polymer cartoon gets the grey base +
      // annotation; every ball-and-stick rep (ligand AND focus sticks) gets the element/grey/green
      // stick layers. This is timing-independent, so colouring is consistent without reloads.
      const isCartoonRep = typeByRef.get(ref) === 'cartoon';
      const node = b.to(ref).apply(OT, { layers: isCartoonRep ? layers : focusLayers });
      overpaintRefs.push(node.ref);
    }
    await b.commit();
    return true;
  }

  // Refs of the static-ligand component's representations.
  function ligandComponentRepRefs() {
    const set = new Set();
    const cur = viewer.plugin.managers.structure.hierarchy.current;
    (cur && cur.structures || []).forEach(st => (st.components || []).forEach(c => {
      if (c.key !== 'structure-component-static-ligand') return;
      (c.representations || []).forEach(r => { if (r.cell && r.cell.transform) set.add(r.cell.transform.ref); });
    }));
    return set;
  }

  // ---- per-ligand visibility -----------------------------------------------------------------
  // Two things hidden via full transparency (value:1), as ONE transparency node per rep so nothing
  // overrides anything else:
  //   1. `hidden` ligand residues — applied to EVERY rep (focus: hide non-focused ligands).
  //   2. AlphaFill donor PEPTIDE residues — applied to the LIGAND rep ONLY. AlphaFill drags the donor's
  //      binding-site amino acids (group_PDB='ATOM') into the static-ligand component; they'd render as
  //      grey/white "faraway side chains". Hiding them only on the ligand rep leaves the real protein
  //      cartoon and the focus pocket (also ATOM) untouched. Re-applied on every call so it always holds.
  async function setLigandVisibility(hidden) {
    _lastLigandHidden = hidden || []; // remembered so focus can re-hide on its newly-created rep
    const structure = currentStructureData();
    if (!structure) return false;
    const Q = L().structure.Queries;
    const SP = L().structure.StructureProperties;
    const QueryContext = L().structure.QueryContext;
    const Bundle = L().structure.StructureElement.Bundle;
    const TT = L().plugin.StateTransforms.Representation.TransparencyStructureRepresentation3DFromBundle;
    const b = viewer.plugin.build();
    ligandTransparencyRefs.forEach(ref => { try { b.delete(ref); } catch (_) {} });
    ligandTransparencyRefs = [];
    let hiddenBundle = null;
    if (hidden && hidden.length) {
      const set = new Set(hidden.map(([c, r]) => (c == null ? '' : c) + '|' + r));
      const anyChain = hidden.some(([c]) => c == null);
      hiddenBundle = Bundle.fromSelection(Q.generators.atoms({
        residueTest: (ctx) => {
          const ch = SP.chain.auth_asym_id(ctx.element), ri = SP.residue.auth_seq_id(ctx.element);
          return set.has(ch + '|' + ri) || (anyChain && set.has('|' + ri));
        },
      })(new QueryContext(structure)));
    }
    // Donor-peptide bundle = all ATOM (polymer) residues — only meaningful where it intersects the
    // ligand rep. Reuse the cached polymer bundle (PERF: don't rebuild it every focus).
    const peptideBundle = structureBundles().polymer;
    const ligRepRefs = ligandComponentRepRefs();
    const carbRefs = carbohydrateRepRefs(); // never make glycan SNFG / interaction lines transparent
    for (const ref of allRepresentationRefs().filter(r => !markerRepRefs.has(r) && !carbRefs.has(r))) {
      const layers = [];
      if (hiddenBundle) layers.push({ bundle: hiddenBundle, value: 1, clear: false });
      if (ligRepRefs.has(ref)) layers.push({ bundle: peptideBundle, value: 1, clear: false }); // peptides: ligand rep only
      if (!layers.length) continue;
      const node = b.to(ref).apply(TT, { layers });
      ligandTransparencyRefs.push(node.ref);
    }
    await b.commit();
    return true;
  }

  // ---- picking: resolve a Mol* loci to a plain atom descriptor for the parent ----------------
  // Resolve any Mol* loci to its first StructureElement.Location. A click/hover on a STICK lands on a
  // Bond.Loci, not an element loci — without this it reads as nothing (NO-ATOM, empty hover). Handle it
  // by kind (works even if the Bond namespace isn't on lib.structure): prefer Bond.toStructureElementLoci,
  // else read the bond's first atom.
  function firstLocation(loci) {
    const SE = L().structure.StructureElement;
    const Bond = L().structure.Bond;
    if (!loci) return null;
    if (loci.kind === 'bond-loci') {
      try { if (Bond && Bond.toStructureElementLoci) { const el = Bond.toStructureElementLoci(loci); if (el && !SE.Loci.isEmpty(el)) return SE.Loci.getFirstLocation(el); } } catch (_) {}
      if (loci.bonds && loci.bonds.length) { const b = loci.bonds[0]; try { return SE.Location.create(loci.structure, b.aUnit, b.aUnit.elements[b.aIndex]); } catch (_) {} }
      return null;
    }
    if (SE.Loci.is(loci) && !SE.Loci.isEmpty(loci)) return SE.Loci.getFirstLocation(loci);
    return null;
  }
  function lociToAtom(loci) {
    const SP = L().structure.StructureProperties;
    const loc = firstLocation(loci);
    if (!loc) return null;
    return {
      chain: SP.chain.auth_asym_id(loc), resi: SP.residue.auth_seq_id(loc),
      resn: SP.atom.auth_comp_id(loc), hetflag: SP.residue.group_PDB(loc) === 'HETATM',
      x: SP.atom.x(loc), y: SP.atom.y(loc), z: SP.atom.z(loc),
    };
  }

  // ---- sphere markers (PTM / variant / site Cα dots) via spacefill components -----------------
  // groups: [{ color, radius, residues: [[chain, resi], ...] }]
  // skipCamRestore: when true (focus/unfocus in progress) skip snapshot+restore so the camera
  // animation is not interrupted. When false (sphere toggle) restore the camera to prevent Mol*'s
  // auto-zoom from shifting the view when components are added/removed.
  let _setMarkersRunning = false;
  let _setMarkersPending = null;
  async function setMarkers(groups, skipCamRestore) {
    // Serialize: if a setMarkers is already running (dataTransaction in flight), queue only the
    // latest call so rapid focus-switches don't race and corrupt markerByColor.
    if (_setMarkersRunning) { _setMarkersPending = { groups, skipCamRestore }; return true; }
    _setMarkersRunning = true;
    try {
      await _doSetMarkers(groups, skipCamRestore);
      while (_setMarkersPending) {
        const { groups: pg, skipCamRestore: ps } = _setMarkersPending;
        _setMarkersPending = null;
        await _doSetMarkers(pg, ps);
      }
    } finally { _setMarkersRunning = false; }
    return true;
  }
  async function _doSetMarkers(groups, skipCamRestore) {
    const structure = currentStructureData();
    if (!structure) return false;
    const Q = L().structure.Queries;
    const SP = L().structure.StructureProperties;
    const QueryContext = L().structure.QueryContext;
    const Bundle = L().structure.StructureElement.Bundle;
    const StateTransforms = L().plugin.StateTransforms;
    const structureRef = currentStructureCell().cell.transform.ref;

    const SelT = StateTransforms.Model.StructureSelectionFromBundle;
    const buildSel = (g) => {
      const set = new Set(g.residues.map(([c, r]) => (c == null ? '' : c) + '|' + r));
      const anyChain = g.residues.some(([c]) => c == null);
      const query = Q.generators.atoms({
        atomTest: (ctx) => SP.atom.label_atom_id(ctx.element) === 'CA',
        residueTest: (ctx) => {
          const ch = SP.chain.auth_asym_id(ctx.element), ri = SP.residue.auth_seq_id(ctx.element);
          return set.has(ch + '|' + ri) || (anyChain && set.has('|' + ri));
        },
      });
      return Bundle.fromSelection(query(new QueryContext(structure)));
    };
    // Group incoming residues by colour. Persistent marker components keyed by colour mean a toggle
    // only UPDATES the affected colour's selection in place (no delete+recreate → no blink).
    const incoming = new Map();
    for (const g of (groups || [])) {
      if (!g.residues || !g.residues.length) continue;
      incoming.set(colorToInt(g.color), { bundle: buildSel(g), radius: g.radius || 1.6 });
    }
    // A bundle that matches nothing — colours that drop out are updated to this (component kept alive).
    const emptyBundle = structureBundles().empty; // cached (PERF)
    // Camera snapshot (clone Vec3s — getSnapshot returns the live state) so add/remove can't shift view.
    // Skipped when skipCamRestore is true (a focus/unfocus animation is in flight — let it complete).
    let camSnap = null;
    if (!skipCamRestore) {
      try {
        const s = viewer.plugin.canvas3d.camera.getSnapshot();
        camSnap = { ...s, position: Array.from(s.position), up: Array.from(s.up), target: Array.from(s.target) };
      } catch (_) {}
    }
    await viewer.plugin.dataTransaction(async () => {
      const b = viewer.plugin.build();
      const creates = [];
      for (const [color, info] of incoming) {
        const existing = markerByColor.get(color);
        if (existing) b.to(existing.sel).update(SelT, (old) => ({ ...old, bundle: info.bundle }));   // in place
        else creates.push({ color, radius: info.radius,
          to: b.to(structureRef).apply(SelT, { bundle: info.bundle, label: 'ufv-marker' }, { tags: 'ufv-marker' }) });
      }
      // Colours no longer present → update their selection to an EMPTY bundle instead of DELETING the
      // component. The component + its spacefill rep persist (rendering nothing), so re-toggling just
      // refills the bundle in place — no delete+recreate, no addRepresentation-after-commit → no flicker,
      // and the component count stays stable so Mol* doesn't auto-zoom (no camera bounce).
      for (const [color, m] of markerByColor) {
        if (!incoming.has(color)) b.to(m.sel).update(SelT, (old) => ({ ...old, bundle: emptyBundle }));
      }
      await b.commit();
      for (const c of creates) {
        const rep = await viewer.plugin.builders.structure.representation.addRepresentation(c.to.ref, {
          type: 'spacefill', color: 'uniform', colorParams: { value: c.color },
          size: 'uniform', sizeParams: { value: c.radius },
        });
        markerByColor.set(c.color, { sel: c.to.ref, rep: rep && rep.ref });
      }
      // Keep markerRepRefs in sync so setCartoon never overpaints spacefill spheres
      markerRepRefs = new Set([...markerByColor.values()].map(m => m.rep).filter(Boolean));
    });
    if (camSnap) {
      // Re-assert the pre-edit camera at spaced intervals. Creating a NEW marker colour component (e.g.
      // toggling Glycosylation on for the first time) makes Mol* auto-fit to the changed bounds — sometimes
      // instantly, sometimes as a short animation that starts several frames AFTER the transaction commits.
      // Snapping back at 0…400 ms (duration 0 cancels any in-flight camera animation) catches it whenever it
      // fires, without continuously fighting a single frame. Skipped during focus/unfocus (camSnap is null).
      const restoreCam = () => { try { viewer.plugin.canvas3d.camera.setState(camSnap, 0); } catch (_) {} };
      restoreCam();                          // SYNCHRONOUS — before any auto-fit render frame can paint
      try { requestAnimationFrame(restoreCam); } catch (_) {}  // and on the very next frame
      [16, 33, 66, 120, 250, 400].forEach(ms => setTimeout(restoreCam, ms));
    }
    return true;
  }

  // Hide / show ALL marker spheres via a full-transparency layer (value:1) on each marker rep. The sphere
  // GEOMETRY stays in the scene, so the bounding sphere is unchanged → no camera auto-fit and no flash;
  // only the alpha flips. Used by the "show other spheres" toggle while focused.
  let _markerHideRefs = [];
  async function setMarkersHidden(hidden) {
    const structure = currentStructureData();
    if (!structure) return true;
    const Q = L().structure.Queries;
    const QueryContext = L().structure.QueryContext;
    const Bundle = L().structure.StructureElement.Bundle;
    const TT = L().plugin.StateTransforms.Representation.TransparencyStructureRepresentation3DFromBundle;
    const b = viewer.plugin.build();
    _markerHideRefs.forEach(ref => { try { b.delete(ref); } catch (_) {} });
    _markerHideRefs = [];
    if (hidden) {
      const allBundle = Bundle.fromSelection(Q.generators.atoms({})(new QueryContext(structure)));
      for (const m of markerByColor.values()) {
        if (!m.rep) continue;
        const node = b.to(m.rep).apply(TT, { layers: [{ bundle: allBundle, value: 1, clear: false }] });
        _markerHideRefs.push(node.ref);
      }
    }
    await b.commit();
    return true;
  }

  // ---- focus (Mol* native): residue/ligand + surroundings + zoom -----------------------------
  function lociForResidue(chain, resi) {
    const structure = currentStructureData();
    const Q = L().structure.Queries, SP = L().structure.StructureProperties;
    const QueryContext = L().structure.QueryContext, StructureSelection = L().structure.StructureSelection;
    const query = Q.generators.atoms({
      residueTest: (ctx) => (chain == null || SP.chain.auth_asym_id(ctx.element) === chain) && SP.residue.auth_seq_id(ctx.element) === resi,
    });
    return StructureSelection.toLociWithSourceUnits(query(new QueryContext(structure)));
  }
  // Build a loci covering all residues in the neighbors list (array of {chain, resi}).
  // Used to set the sticks representation to exactly the nearby-distance set.
  function lociForResidueList(neighbors) {
    const structure = currentStructureData();
    if (!structure || !neighbors.length) return null;
    const Q = L().structure.Queries, SP = L().structure.StructureProperties;
    const QueryContext = L().structure.QueryContext, StructureSelection = L().structure.StructureSelection;
    const keys = new Set(neighbors.map(r => (r.chain == null ? '' : r.chain) + '|' + r.resi));
    const query = Q.generators.atoms({
      residueTest: ctx => {
        const ri = SP.residue.auth_seq_id(ctx.element);
        return keys.has(SP.chain.auth_asym_id(ctx.element) + '|' + ri) || keys.has('|' + ri);
      },
    });
    return StructureSelection.toLociWithSourceUnits(query(new QueryContext(structure)));
  }
  // Bundle covering a list of [chain, resi] residues (focused residue/ligand + nearby protein neighbours).
  function bundleForResidueList(residues) {
    const structure = currentStructureData();
    const Q = L().structure.Queries, SP = L().structure.StructureProperties;
    const QueryContext = L().structure.QueryContext, Bundle = L().structure.StructureElement.Bundle;
    const set = new Set(residues.map(([c, r]) => (c == null ? '' : c) + '|' + r));
    const anyChain = residues.some(([c]) => c == null);
    const query = Q.generators.atoms({
      residueTest: (ctx) => {
        const ch = SP.chain.auth_asym_id(ctx.element), ri = SP.residue.auth_seq_id(ctx.element);
        return set.has(ch + '|' + ri) || (anyChain && set.has('|' + ri));
      },
    });
    return Bundle.fromSelection(query(new QueryContext(structure)));
  }
  async function removeFocusRep() {
    focusSelRef = null; focusRepRef = null;
    // Delete EVERY 'ufv-focus'-tagged component (not just the tracked ref) so a leftover focus rep from
    // a rapid re-click can't linger as stray faraway sticks.
    try {
      const refs = [];
      viewer.plugin.state.data.cells.forEach((cell, ref) => {
        const tags = cell.transform && cell.transform.tags;
        const has = Array.isArray(tags) ? tags.indexOf('ufv-focus') >= 0 : tags === 'ufv-focus';
        if (has) refs.push(ref);
      });
      if (refs.length) {
        const b = viewer.plugin.build();
        refs.forEach(r => { try { b.delete(r); } catch (_) {} });
        await b.commit();
      }
    } catch (_) {}
  }
  // ---- predicted-pocket highlight: a translucent molecular surface over the pocket's residues, so the
  // user can SEE the predicted cavity (clicking a pocket in the side panel). Tagged 'ufv-pocket' so it's
  // removed wholesale; independent of the focus rep. residues = [{chain, resi}].
  async function doClearPocket() {
    try {
      const refs = [];
      viewer.plugin.state.data.cells.forEach((cell, ref) => {
        const tags = cell.transform && cell.transform.tags;
        const has = Array.isArray(tags) ? tags.indexOf('ufv-pocket') >= 0 : tags === 'ufv-pocket';
        if (has) refs.push(ref);
      });
      if (refs.length) {
        const b = viewer.plugin.build();
        refs.forEach(r => { try { b.delete(r); } catch (_) {} });
        await b.commit();
      }
    } catch (_) {}
  }
  async function doShowPocket(residues, surfColor) {
    if (!currentStructureData() || !residues || !residues.length) return false;
    await removeFocusRep();   // the pocket view supersedes the single-residue focus sticks
    await doClearPocket();
    try {
      const SelT = L().plugin.StateTransforms.Model.StructureSelectionFromBundle;
      const structureRef = currentStructureCell().cell.transform.ref;
      // 1) translucent molecular surface over all pocket residues (the cavity shape).
      const b = viewer.plugin.build();
      const surfSel = b.to(structureRef).apply(SelT, { bundle: bundleForResidueList(residues.map(r => [r.chain, r.resi])), label: 'ufv-pocket' }, { tags: 'ufv-pocket' });
      await b.commit();
      await viewer.plugin.builders.structure.representation.addRepresentation(surfSel.ref, {
        type: 'molecular-surface', typeParams: { alpha: 0.4, ignoreHydrogens: true },
        color: 'uniform', colorParams: { value: colorToInt(surfColor || '#26c6da') },
      });
      // 2) ball-and-stick for the residues, coloured by annotation (per-colour groups; uncoloured = element).
      const byColor = new Map(); const noColor = [];
      for (const r of residues) { if (r.color) { if (!byColor.has(r.color)) byColor.set(r.color, []); byColor.get(r.color).push([r.chain, r.resi]); } else noColor.push([r.chain, r.resi]); }
      const addSticks = async (resList, colorHex) => {
        const bb = viewer.plugin.build();
        const sel = bb.to(structureRef).apply(SelT, { bundle: bundleForResidueList(resList), label: 'ufv-pocket' }, { tags: 'ufv-pocket' });
        await bb.commit();
        await viewer.plugin.builders.structure.representation.addRepresentation(sel.ref, {
          type: 'ball-and-stick', typeParams: { visuals: ['element-sphere', 'intra-bond'], sizeFactor: 0.22, ignoreHydrogens: true },
          ...(colorHex ? { color: 'uniform', colorParams: { value: colorToInt(colorHex) } } : { color: 'element-symbol' }),
        });
      };
      for (const [colorHex, resList] of byColor) await addSticks(resList, colorHex);
      if (noColor.length) await addSticks(noColor, null);
      try { viewer.plugin.managers.camera.focusLoci(lociForResidueList(residues), { extraRadius: 4 }); } catch (_) {}
    } catch (_) { return false; }
    return true;
  }
  // Focus = OUR OWN ball-and-stick rep of exactly the focused residue/ligand + its nearby set. We do not
  // use Mol*'s focus manager: it only dims (doesn't remove) the rest, adds its own ~5 Å surroundings, and
  // builds the rep across several async transactions (the colour/blob race). Building our own rep is
  // synchronous and exact — other ligands are genuinely hidden, only this set shows, colour is right
  // immediately. visuals exclude 'inter-bond' so packed ligands never get spurious bonds.
  // neighbors: [{chain, resi}] (nearby protein). skipCam: don't move the camera (filter/colour re-focus).
  async function doFocus(chain, resi, neighbors, skipCam, annotations) {
    if (!currentStructureData()) return false;
    const myGen = ++_focusGen; // this focus's generation (deferred work below bails if superseded)
    // (No camera snapshot/restore needed here: canvas3d.manualReset=true disables Mol*'s auto-fit, so the
    // rep rebuild below can't move the camera. The parent's camFocus is the sole camera mover.)
    await removeFocusRep();
    _focusedResidue = { chain, resi };
    _focusAnnotations = annotations || null;
    const residues = [[chain, resi], ...((neighbors || []).map(n => [n.chain, n.resi]))];
    try {
      const SelT = L().plugin.StateTransforms.Model.StructureSelectionFromBundle;
      const structureRef = currentStructureCell().cell.transform.ref;
      const b = viewer.plugin.build();
      const sel = b.to(structureRef).apply(SelT, { bundle: bundleForResidueList(residues), label: 'ufv-focus' }, { tags: 'ufv-focus' });
      await b.commit();
      focusSelRef = sel.ref;
      const rep = await viewer.plugin.builders.structure.representation.addRepresentation(sel.ref, {
        type: 'ball-and-stick',
        typeParams: { visuals: ['element-sphere', 'intra-bond'], sizeFactor: 0.22 },
        color: 'element-symbol',
      }, { initialState: { isHidden: true } }); // created HIDDEN — revealed only after setCartoon colours it
      focusRepRef = rep && rep.ref;
      // Non-covalent interactions (H-bonds, salt bridges, …) are the heaviest step, so they're deferred OFF
      // the focus path. For a real focus we hand the request to doCamFocus, which fires it AFTER the camera
      // animation finishes (so the compute never hiccups the zoom). For a skipCam re-colour/re-filter there's
      // no animation, so schedule it directly. Either way it's skipped if a newer focus supersedes this one.
      if (skipCam) {
        const ref = sel.ref;
        setTimeout(() => addFocusInteractions(ref, myGen), 260);
      } else {
        _pendingInter = { ref: sel.ref, gen: myGen };
      }
    } catch (_) {}
    // Force inter-bond OFF on every ball-and-stick rep (incl. our focus rep): otherwise the focused
    // ligand renders tiny stub bonds to atoms of nearby (now-hidden) ligands/residues, since Mol* bonds
    // very-close atoms across chains.
    try { await dropLigandInterBonds(); } catch (_) {}
    // Colour the new focus rep (ball-and-stick → grey polymer carbons / annotation / green ligand) and
    // keep the overview colouring, via the same type-based overpaint. Synchronous → no flash.
    if (_lastCartoonPayload) { try { await setCartoon(_lastCartoonPayload.base, _lastCartoonPayload.groups); } catch (_) {} }
    // Reveal the now-coloured focus sticks (created hidden so the element-symbol intermediate — the
    // "default Mol* colouring flash" when switching residues — is never shown).
    if (focusRepRef) { try { await viewer.plugin.state.data.updateCellState(focusRepRef, { isHidden: false }); } catch (_) {} }
    // The CAMERA is not moved here — the parent sends `camFocus` LAST so the zoom animation is uncontended
    // (and manualReset=true means the rep rebuild above never auto-fit the camera).
    return true;
  }
  // Camera focus on a residue/ligand by chain + auth residue — sent LAST so the animation is uncontended.
  // panOnly: keep the CURRENT zoom (camera radius) and just recentre on the new residue — used when moving
  // residue→residue while already zoomed in, so there's no extra zoom-out/zoom-in between residues.
  function doCamFocus(chain, resi, radius, panOnly) {
    try {
      let r = radius;
      if (panOnly) {
        const cur = viewer.plugin.canvas3d?.camera?.getSnapshot?.();
        if (cur && cur.radius > 0) r = cur.radius; // re-frame at the existing zoom (pure pan)
      }
      // Centre the camera on the CLICKED entity (residue or ligand) — focusLoci centres on the loci's own
      // sphere, so no "jump" to a neighbour centroid. Apply r as an ABSOLUTE radius via minRadius
      // (extraRadius 0). (focusLoci radius = max(sphere+extra, min).)
      const opts = r ? { minRadius: r, extraRadius: 0 } : undefined;
      viewer.plugin.managers.camera.focusLoci(lociForResidue(chain, resi), opts);
    } catch (_) {}
    // Now that the zoom animation has STARTED, schedule this focus's interactions rep to be added AFTER it
    // finishes (~250 ms) — so the heavy compute never hiccups the camera motion. Cleared so a later
    // rezoom (slider) camFocus, which has no pending focus, doesn't re-add it.
    const pend = _pendingInter; _pendingInter = null;
    if (pend) setTimeout(() => addFocusInteractions(pend.ref, pend.gen), 320);
  }
  async function doUnfocus() {
    _focusedResidue = null;
    _focusAnnotations = null;
    _focusGen++; // cancel any pending deferred interactions from the focus we're leaving
    await removeFocusRep();
    resetCamera(300);
    return true;
  }

  // ---- core operations -----------------------------------------------------------------------
  function applyBg() {
    try { viewer.plugin.canvas3d?.setProps({ renderer: { backgroundColor: _bgColor } }); } catch (_) {}
  }
  async function afterLoad() {
    overpaintRefs = []; ligandTransparencyRefs = []; markerByColor = new Map(); labelMap = {}; _markerHideRefs = [];
    _sBundles = null; // structure changed → invalidate the cached element bundles
    _focusedResidue = null; _focusAnnotations = null; focusSelRef = null; focusRepRef = null;
    applyBg();
    const structure = currentStructureData();
    if (structure) post({ evt: 'atoms', atoms: extractAtoms(structure) });
    // Override Mol*'s default preset colouring (pLDDT orange/yellow for AlphaFold) with our neutral base
    // ASAP, before the modal's first real setCartoon — otherwise the protein flashes orange on load.
    if (structure && _lastCartoonPayload) { try { await setCartoon(_lastCartoonPayload.base, _lastCartoonPayload.groups); } catch (_) {} }
    // manualReset=true disables Mol*'s auto-fit, so frame the freshly-loaded structure explicitly here
    // (requestCameraReset is an explicit command and works regardless of manualReset).
    if (structure) resetCamera(0);
    await new Promise(r => setTimeout(r, 0));
    // The preset renders the protein as cartoon-only (verified), so there are no protein side-chain
    // sticks to hide here — any protein sticks seen in the overview are the leftover focus rep.
    // Snapshot the reps now present. Focus sticks created later are NOT in this set.
    initialRepRefs = new Set(allRepresentationRefs());
    await dropLigandInterBonds();
    // Late-rep safety net: ligand / SNFG reps can finish creating just after the preset settles, so the
    // first pass can miss them (the "needs a couple of reloads" bug). Re-apply once.
    setTimeout(() => {
      if (_lastCartoonPayload) setCartoon(_lastCartoonPayload.base, _lastCartoonPayload.groups).catch(() => {});
      dropLigandInterBonds().catch(() => {});
      setLigandVisibility(_lastLigandHidden).catch(() => {}); // applies the donor-peptide hide on the ligand rep
    }, 450);
  }

  // AlphaFill (and any model) packs distinct ligands close together; Mol*'s ball-and-stick draws
  // spurious INTER-unit bonds between them (each transplant is its own chain/unit), producing the
  // "blob". The data is correct — only the rendering is wrong — so we drop the 'inter-bond' visual
  // from ligand ball-and-stick reps, keeping each ligand's own (intra-unit) bonds intact.
  async function dropLigandInterBonds() {
    const refs = [];
    const cur = viewer.plugin.managers.structure.hierarchy.current;
    (cur && cur.structures || []).forEach(st => (st.components || []).forEach(c =>
      (c.representations || []).forEach(r => {
        const t = r.cell && r.cell.transform && r.cell.transform.params && r.cell.transform.params.type;
        if (t && t.name === 'ball-and-stick' && r.cell.transform) refs.push(r.cell.transform.ref);
      })));
    if (!refs.length) return;
    try {
      const b = viewer.plugin.build();
      for (const ref of refs) {
        b.to(ref).update(o => {
          if (o.type && o.type.name === 'ball-and-stick' && o.type.params) {
            o.type.params.visuals = ['element-sphere', 'intra-bond'];
          }
        });
      }
      await b.commit();
    } catch (_) { /* representation params shape changed — leave default */ }
  }
  // Neutral grey until the modal sends this structure's real colouring — so neither the PREVIOUS
  // structure's overpaint nor Mol*'s default AlphaFold pLDDT theme (orange/yellow) flashes through.
  const NEUTRAL_CARTOON = { base: '#d0d0d0', groups: [] };
  async function loadData(data, format, isBinary) {
    _lastCartoonPayload = NEUTRAL_CARTOON;
    await viewer.plugin.clear();
    await viewer.loadStructureFromData(data, format || 'mmcif', !!isBinary);
    await afterLoad();
  }
  async function loadUrl(url, format, isBinary) {
    _lastCartoonPayload = NEUTRAL_CARTOON;
    await viewer.plugin.clear();
    await viewer.loadStructureFromUrl(url, format || 'mmcif', !!isBinary);
    await afterLoad();
  }
  async function clearModel() { overpaintRefs = []; markerByColor = new Map(); labelMap = {}; await viewer.plugin.clear(); }
  function resetCamera(durationMs) { try { viewer.plugin.canvas3d?.requestCameraReset({ durationMs: durationMs ?? 0 }); } catch (_) {} }
  function handleResize() { try { viewer.handleResize(); } catch (_) {} }
  async function screenshot() {
    const vs = viewer.plugin.helpers?.viewportScreenshot;
    if (!vs) throw new Error('screenshot helper unavailable');
    return await vs.getImageDataUri();
  }
  function setBackground(color) {
    try { viewer.plugin.canvas3d?.setProps({ renderer: { backgroundColor: colorToInt(color) } }); } catch (_) {}
  }

  // ---- command dispatch ----------------------------------------------------------------------
  const HANDLERS = {
    async loadData(m) { await loadData(m.data, m.format, m.isBinary); return true; },
    async loadUrl(m) { await loadUrl(m.url, m.format, m.isBinary); return true; },
    async clear() { await clearModel(); return true; },
    reset(m) { resetCamera(m.duration); return true; },
    resize() { handleResize(); return true; },
    background(m) { setBackground(m.color); return true; },
    async screenshot() { return await screenshot(); },
    async setCartoon(m) { _lastCartoonPayload = m; return await setCartoon(m.base, m.groups); },
    async setLigandVisibility(m) { return await setLigandVisibility(m.hidden); },
    async setMarkers(m) { return await setMarkers(m.groups, m.skipCamRestore); },
    markersHidden(m) { return setMarkersHidden(m.hidden); },
    setLabels(m) { labelMap = m.map || {}; return true; },
    focus(m) { return doFocus(m.chain, m.resi, m.neighbors, false, m.annotations); },
    refocus(m) { return doFocus(m.chain, m.resi, m.neighbors, true, m.annotations); },
    camFocus(m) { doCamFocus(m.chain, m.resi, m.radius, m.panOnly); return true; },
    showPocket(m) { return doShowPocket(m.residues, m.color); },
    clearPocket() { return doClearPocket(); },
    frameResidues(m) { try { if (m.residues && m.residues.length) viewer.plugin.managers.camera.focusLoci(lociForResidueList(m.residues), { extraRadius: 4 }); } catch (_) {} return true; },
    unfocus() { return doUnfocus(); },
    background(m) { _bgColor = parseInt((m.color || '#0c111b').replace('#', ''), 16); applyBg(); return true; },
  };

  // Commands are SERIALISED through this queue. onMessage is async, so without a queue two messages
  // (e.g. setLigandVisibility then focus) interleave their awaited state transactions and race on the
  // same Mol* reps — the first click on an AlphaFill ligand then leaves the other transplants visible
  // (blob, chain-id blue/cyan) and the focused ligand un-recoloured (blue). The queue guarantees each
  // command's state edits fully commit before the next command starts.
  let _cmdQueue = Promise.resolve();
  window.addEventListener('message', (ev) => {
    const m = ev.data;
    if (!m || m.ns !== NS || !m.cmd) return;
    _cmdQueue = _cmdQueue.then(async () => {
      const h = HANDLERS[m.cmd];
      if (!h) { fail(m.reqId, 'unknown command: ' + m.cmd); return; }
      if (!ready) { fail(m.reqId, 'viewer not ready'); return; }
      try {
        const value = await h(m);
        if (m.reqId != null) reply(m.reqId, true, value);
      } catch (e) { fail(m.reqId, (e && e.message) || e); }
    });
  });

  // ---- init + self-test ----------------------------------------------------------------------
  async function init() {
    if (!window.molstar || !molstar.Viewer) {
      setStatus('Mol* bundle failed to load (window.molstar missing).', true);
      post({ evt: 'error', message: 'molstar-bundle-missing' });
      return;
    }
    try {
      viewer = await molstar.Viewer.create('app', {
        layoutIsExpanded: false, layoutShowControls: false, layoutShowSequence: false,
        layoutShowLog: false, layoutShowLeftPanel: false, viewportShowExpand: false,
        viewportShowSelectionMode: false, viewportShowAnimation: false,
        pdbProvider: 'rcsb', emdbProvider: 'rcsb',
      });
    } catch (e) {
      setStatus('Mol* Viewer.create failed: ' + (e && e.message || e), true);
      post({ evt: 'error', message: 'viewer-create-failed' });
      return;
    }
    // Orthographic projection + no bottom-left axis gizmo + initial background colour.
    // pickPadding widens the pick radius so thin ligand sticks aren't a pixel-perfect target (was the
    // "10 clicks to select a ligand / NO-ATOM" problem).
    // manualReset:true disables Mol*'s automatic camera fit on scene/bounds changes — so adding/removing the
    // focus rep, the green highlight sphere, or marker components never auto-zooms. OUR camFocus is then the
    // ONLY camera mover, which removes the "double zoom" (auto-fit + camFocus) on a residue/ligand focus.
    try { viewer.plugin.canvas3d?.setProps({ camera: { mode: 'orthographic', manualReset: true, helper: { axes: { name: 'off', params: {} } } }, renderer: { backgroundColor: _bgColor }, pickPadding: 5 }); } catch (_) {}
    // SUPPRESS *all* of Mol*'s native click reactions so a click on the 3D structure does ONLY our
    // focusResidue — byte-for-byte the same code path as clicking a Nearby chip in the side panel (which
    // never sends a click into Mol* at all). Without this, a single primary click on the structure fires
    // THREE separate camera/geometry reactions that stack on top of (and stutter) our own camFocus:
    //   • `camera-focus-loci`        — `clickCenterFocus` binding ⇒ camera.focusLoci(clicked)   (native zoom #1)
    //   • `representation-focus-loci`— `clickFocus` binding ⇒ sets focus + camera.focusLoci(focus) (native zoom #2)
    //   • `structure-focus-representation` — draws the target + EVERY residue within `expandRadius` (5 Å) as
    //     ball-and-stick (~900 atoms) plus a non-covalent `interactions` calc (the async geometry build)
    // That stack is the "three levels of zoom" / "not as smooth as Nearby" the user reported.
    //
    // Fix, applied to the live behaviour cells (matched by their params, not the minified transformer name):
    //   1. Empty the click-focus bindings on the two camera behaviours so a click NEVER moves the camera
    //      natively — our parent-driven camFocus becomes the sole, single zoom (identical to the Nearby path).
    //   2. Zero `expandRadius` on the focus-representation so even if some other path sets the focus, the
    //      expensive surroundings shell + interactions are never built.
    const EMPTY_BINDING = () => ({ triggers: [], action: '', description: '' });
    try {
      const bState = viewer.plugin.state.behaviors;
      bState.cells.forEach((cell) => {
        const vals = cell.params && cell.params.values;
        const tr = cell.transform && cell.transform.transformer;
        if (!tr || !vals) return;
        if (typeof vals.expandRadius === 'number') {
          try { viewer.plugin.state.updateBehavior(tr, (p) => { p.expandRadius = 0; }); } catch (_) {}
        }
        const b = vals.bindings;
        if (b && (b.clickCenterFocus || b.clickFocus)) {
          try {
            viewer.plugin.state.updateBehavior(tr, (p) => {
              for (const k in p.bindings) {
                if (/^click(Center)?Focus/.test(k)) p.bindings[k] = EMPTY_BINDING();
              }
            });
          } catch (_) {}
        }
      });
    } catch (_) {}
    // Belt-and-suspenders: if anything still sets the structure focus (hover-driven paths, future bindings),
    // clear it immediately so it never draws over our own `ufv-focus`. Clearing re-emits `undefined`, which
    // the `if (e)` guard ignores (no loop).
    try {
      viewer.plugin.managers.structure.focus.behaviors.current.subscribe((e) => {
        if (e) Promise.resolve().then(() => { try { viewer.plugin.managers.structure.focus.clear(); } catch (_) {} });
      });
    } catch (_) {}
    setupInteraction();
    ready = true;
    setStatus('');
    post({ evt: 'ready' });
    runSelfTestFromQuery();
  }

  // Our own hover tooltip (Mol*'s default label is hidden via CSS): on hover over an ANNOTATED
  // residue, show our annotation text at the cursor. Clicks are forwarded to the parent (report).
  let mouseX = 0, mouseY = 0;
  const tipEl = document.getElementById('ufv-tip');
  function tipReposition() {
    if (!tipEl || tipEl.style.display !== 'block') return;
    const pad = 14, w = tipEl.offsetWidth, h = tipEl.offsetHeight;
    let x = mouseX + pad, y = mouseY + pad;
    if (x + w > window.innerWidth) x = mouseX - w - pad;
    if (y + h > window.innerHeight) y = mouseY - h - pad;
    tipEl.style.left = Math.max(2, x) + 'px';
    tipEl.style.top = Math.max(2, y) + 'px';
  }
  function tipShow(resn, ri, bodyHtml) {
    if (!tipEl) return;
    const body = bodyHtml ? `<div class="ufv-tip-body">${bodyHtml}</div>` : '';
    tipEl.innerHTML = `<div class="ufv-tip-hdr">${resn} ${ri}</div>${body}`;
    tipEl.style.display = 'block';
    tipReposition();
  }
  function tipHide() { if (tipEl) tipEl.style.display = 'none'; }
  function setupInteraction() {
    document.addEventListener('mousemove', (e) => { mouseX = e.clientX; mouseY = e.clientY; tipReposition(); });
    document.addEventListener('dblclick', () => post({ evt: 'pick', kind: 'dblclick' }));
    try {
      const SP = L().structure.StructureProperties;
      const SE = L().structure.StructureElement;
      viewer.plugin.behaviors.interaction.click.subscribe((e) => {
        const loci = e && e.current && e.current.loci;
        const atom = lociToAtom(loci);
        if (atom) post({ evt: 'pick', kind: 'click', atom });
        else if (loci && loci.kind === 'empty-loci') post({ evt: 'pick', kind: 'background' });
      });
      viewer.plugin.behaviors.interaction.hover.subscribe((e) => {
        const loc = firstLocation(e && e.current && e.current.loci); // handles stick (bond) hover too
        if (!loc) { tipHide(); return; }
        const ri = SP.residue.auth_seq_id(loc);
        const resn = SP.atom.auth_comp_id(loc);
        // Ligands have no per-residue annotation — show just the component name (no "Ligand" label).
        if (SP.residue.group_PDB(loc) === 'HETATM') { tipShow(resn, ri, ''); return; }
        const body = labelMap[SP.chain.auth_asym_id(loc) + '|' + ri] || labelMap['|' + ri] || null;
        if (!body) { tipHide(); return; }
        tipShow(resn, ri, body);
      });
    } catch (_) {}
  }

  async function loadAlphaFold(acc) {
    setStatus('Resolving AlphaFold model…');
    try {
      const r = await fetch(AF_API(acc));
      if (!r.ok) throw new Error('AlphaFold API ' + r.status);
      const arr = await r.json();
      const cifUrl = arr && arr[0] && (arr[0].cifUrl || arr[0].pdbUrl);
      if (!cifUrl) throw new Error('no model URL in AlphaFold API response');
      setStatus('Loading structure…');
      await loadUrl(cifUrl, cifUrl.endsWith('.pdb') ? 'pdb' : 'mmcif', false);
      setStatus('');
    } catch (e) { setStatus('Self-test load failed: ' + (e && e.message || e), true); }
  }
  function runSelfTestFromQuery() {
    const q = new URLSearchParams(location.search);
    const af = q.get('af'), url = q.get('url');
    if (af) { loadAlphaFold(af.toUpperCase()); return; }
    if (url) {
      setStatus('Loading structure…');
      loadUrl(url, q.get('format') || 'mmcif', q.get('binary') === '1')
        .then(() => setStatus('')).catch(e => setStatus('load failed: ' + (e && e.message || e), true));
    }
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
