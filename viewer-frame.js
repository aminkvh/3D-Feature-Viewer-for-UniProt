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
  let markerByColor = new Map(); // colorInt -> {sel, rep} refs; persistent so toggles update in place
  let markerRepRefs = new Set(); // rep refs of all spacefill markers — excluded from cartoon overpaint
  let initialRepRefs = new Set(); // reps present right after structure load (cartoon backbone) — focus sticks appear later and are NOT in this set
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

    const bundleForSet = (resSet) => {
      const query = Q.generators.atoms({
        residueTest: (ctx) => resSet.has(SP.chain.auth_asym_id(ctx.element) + '|' + SP.residue.auth_seq_id(ctx.element)),
      });
      return Bundle.fromSelection(query(new QueryContext(structure)));
    };
    const allBundle = Bundle.fromSelection(Q.generators.all(new QueryContext(structure)));

    // Cartoon layers: base paints everything grey, then annotation groups override per residue.
    const layers = [{ bundle: allBundle, color: colorToInt(base), clear: false }];
    for (const g of (groups || [])) {
      if (!g.residues || !g.residues.length) continue;
      const set = new Set(g.residues.map(([c, r]) => (c == null ? '' : c) + '|' + r));
      const anyChain = g.residues.some(([c]) => c == null);
      const bundle = anyChain
        ? Bundle.fromSelection(Q.generators.atoms({
            residueTest: (ctx) => {
              const ch = SP.chain.auth_asym_id(ctx.element), ri = SP.residue.auth_seq_id(ctx.element);
              return set.has(ch + '|' + ri) || set.has('|' + ri);
            },
          })(new QueryContext(structure)))
        : bundleForSet(set);
      layers.push({ bundle, color: colorToInt(g.color), clear: false });
    }

    // Focus/stick layers, built bottom-up (CPK convention: heteroatoms by element, carbons by category):
    //   1) clear → revert to base theme (gives correct N=blue, O=red, S=yellow; but carbons are
    //      teal because Mol*'s element-symbol theme defaults carbonColor to chain-id),
    //   2) paint ALL carbons grey → CPK carbon baseline (non-annotated residues & selected residue),
    //   3) annotation layers → recolour ONLY the carbons of each annotated neighbour by its colour.
    // The selected residue is excluded from (3) so its carbons stay grey = plain Mol* default look.
    // To exclude it, match against both the chain-qualified key ('A|169') AND the chain-agnostic
    // key ('|169') — annotations may store chain:null while focusedResidue may have 'A'.
    const focusKey     = _focusedResidue ? (_focusedResidue.chain || '') + '|' + _focusedResidue.resi : null;
    const focusKeyAny  = _focusedResidue ? '|' + _focusedResidue.resi : null;
    const carbonBundle = Bundle.fromSelection(Q.generators.atoms({
      atomTest: (ctx) => SP.atom.type_symbol(ctx.element) === 'C',
    })(new QueryContext(structure)));
    const focusLayers = [
      { bundle: allBundle, color: 0, clear: true },        // revert to element theme
      { bundle: carbonBundle, color: 0x909090, clear: false }, // carbons → CPK grey
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
          const ch = SP.chain.auth_asym_id(ctx.element), ri = SP.residue.auth_seq_id(ctx.element);
          return focusSet.has(ch + '|' + ri) || focusSet.has('|' + ri);
        },
      })(new QueryContext(structure)));
      focusLayers.push({ bundle: focusBundle, color: colorToInt(g.color), clear: false });
    }

    const OT = L().plugin.StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle;
    const b = viewer.plugin.build();
    overpaintRefs.forEach(ref => { try { b.delete(ref); } catch (_) {} });
    overpaintRefs = [];
    for (const ref of allRepresentationRefs().filter(r => !markerRepRefs.has(r))) {
      // Reps present at load time are cartoon — use grey base + annotation groups.
      // Reps added later (focus sticks) use clear+annotation only to preserve element colors.
      const isCartoon = initialRepRefs.size > 0 && initialRepRefs.has(ref);
      const node = b.to(ref).apply(OT, { layers: isCartoon ? layers : focusLayers });
      overpaintRefs.push(node.ref);
    }
    await b.commit();
    return true;
  }

  // ---- picking: resolve a Mol* loci to a plain atom descriptor for the parent ----------------
  function lociToAtom(loci) {
    const SE = L().structure.StructureElement;
    if (!loci || !SE.Loci.is(loci) || SE.Loci.isEmpty(loci)) return null;
    const loc = SE.Loci.getFirstLocation(loci);
    if (!loc) return null;
    const SP = L().structure.StructureProperties;
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
      for (const [color, m] of [...markerByColor]) {                 // remove colours no longer present
        if (!incoming.has(color)) { try { b.delete(m.sel); } catch (_) {} markerByColor.delete(color); }
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
      const restoreCam = () => { try { viewer.plugin.canvas3d.camera.setState(camSnap, 0); } catch (_) {} };
      restoreCam(); requestAnimationFrame(restoreCam);
    }
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
  // neighbors: array of {chain, resi} covering the main residue + all nearby ones.
  // skipCam: true means update sticks only (no camera zoom).
  function doFocus(chain, resi, neighbors, skipCam, annotations) {
    if (!currentStructureData()) return false;
    const mainLoci = lociForResidue(chain, resi);
    let sticksLoci = mainLoci;
    if (neighbors && neighbors.length > 0) {
      try { const nb = lociForResidueList(neighbors); if (nb) sticksLoci = nb; } catch (_) {}
    }
    // Record selected residue and sphere annotation colors for stick coloring
    _focusedResidue = { chain, resi };
    _focusAnnotations = annotations || null;
    try { viewer.plugin.managers.structure.focus.setFromLoci(sticksLoci); } catch (_) {}
    if (!skipCam) {
      try { viewer.plugin.managers.camera.focusLoci(mainLoci); } catch (_) {}
    }
    // Re-apply annotation overpaint once Mol* finishes creating the focus sticks.
    // Use two-phase subscription: wait for isUpdating→true then →false (avoids firing immediately
    // on subscribe since BehaviorSubject emits the current value right away).
    if (_lastCartoonPayload) {
      const payload = _lastCartoonPayload;
      let fired = false;
      const run = () => {
        if (fired) return; fired = true;
        setCartoon(payload.base, payload.groups).catch(() => {});
      };
      try {
        let sawTrue = false;
        const sub = viewer.plugin.behaviors.state.isUpdating.subscribe(isUpdating => {
          if (isUpdating) { sawTrue = true; }
          else if (sawTrue) { try { sub.unsubscribe(); } catch (_) {} run(); }
        });
        // Fallback: if update already completed before we subscribed, run after a frame
        setTimeout(() => { try { sub.unsubscribe(); } catch (_) {} run(); }, 300);
      } catch (_) {
        setTimeout(run, 150);
      }
    }
    return true;
  }
  function doUnfocus() {
    _focusedResidue = null;
    _focusAnnotations = null;
    try { viewer.plugin.managers.structure.focus.clear(); } catch (_) {}
    resetCamera(300);
    return true;
  }

  // ---- core operations -----------------------------------------------------------------------
  function applyBg() {
    try { viewer.plugin.canvas3d?.setProps({ renderer: { backgroundColor: _bgColor } }); } catch (_) {}
  }
  async function afterLoad() {
    overpaintRefs = []; markerByColor = new Map(); labelMap = {};
    _focusedResidue = null; _focusAnnotations = null;
    applyBg();
    const structure = currentStructureData();
    if (structure) post({ evt: 'atoms', atoms: extractAtoms(structure) });
    // Snapshot the reps present after load (cartoon backbone).
    // Focus sticks created later will NOT be in this set → get element coloring, not grey base.
    await new Promise(r => setTimeout(r, 0));
    initialRepRefs = new Set(allRepresentationRefs());
  }
  async function loadData(data, format, isBinary) {
    await viewer.plugin.clear();
    await viewer.loadStructureFromData(data, format || 'mmcif', !!isBinary);
    await afterLoad();
  }
  async function loadUrl(url, format, isBinary) {
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
    async setMarkers(m) { return await setMarkers(m.groups, m.skipCamRestore); },
    setLabels(m) { labelMap = m.map || {}; return true; },
    focus(m) { return doFocus(m.chain, m.resi, m.neighbors, false, m.annotations); },
    refocus(m) { return doFocus(m.chain, m.resi, m.neighbors, true, m.annotations); },
    unfocus() { return doUnfocus(); },
    background(m) { _bgColor = parseInt((m.color || '#0c111b').replace('#', ''), 16); applyBg(); return true; },
  };

  window.addEventListener('message', async (ev) => {
    const m = ev.data;
    if (!m || m.ns !== NS || !m.cmd) return;
    const h = HANDLERS[m.cmd];
    if (!h) { fail(m.reqId, 'unknown command: ' + m.cmd); return; }
    if (!ready) { fail(m.reqId, 'viewer not ready'); return; }
    try {
      const value = await h(m);
      if (m.reqId != null) reply(m.reqId, true, value);
    } catch (e) { fail(m.reqId, (e && e.message) || e); }
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
    try { viewer.plugin.canvas3d?.setProps({ camera: { mode: 'orthographic', helper: { axes: { name: 'off', params: {} } } }, renderer: { backgroundColor: _bgColor } }); } catch (_) {}
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
    tipEl.innerHTML = `<div class="ufv-tip-hdr">${resn} ${ri}</div><div class="ufv-tip-body">${bodyHtml}</div>`;
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
        const atom = lociToAtom(e && e.current && e.current.loci);
        if (atom) post({ evt: 'pick', kind: 'click', atom });
      });
      viewer.plugin.behaviors.interaction.hover.subscribe((e) => {
        const loci = e && e.current && e.current.loci;
        if (!loci || !SE.Loci.is(loci) || SE.Loci.isEmpty(loci)) { tipHide(); return; }
        const loc = SE.Loci.getFirstLocation(loci); if (!loc) { tipHide(); return; }
        const ri = SP.residue.auth_seq_id(loc);
        const body = labelMap[SP.chain.auth_asym_id(loc) + '|' + ri] || labelMap['|' + ri] || null;
        if (!body) { tipHide(); return; }
        tipShow(SP.atom.auth_comp_id(loc), ri, body);
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
