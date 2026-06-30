/*
 * viewer-molstar.js — content-script proxy mirroring the StructureViewer (3Dmol) API but rendering
 * via Mol* inside a sandboxed iframe (viewer-frame.html). Real rendering happens in the iframe; this
 * object owns the UniProt↔PDB residue mapping and posts draw commands.
 *
 * Flag-gated: replaces the global StructureViewer only when Mol* is selected
 * (localStorage 'ufv-viewer-engine' === 'molstar', or window.__UFV_VIEWER__ === 'molstar').
 *
 * The residue-mapping + colouring methods (_resiToPdb, residueSelector, allChainsAddStyle,
 * applyCartoonColoring, _applyModeStyles, residueGeometry, …) are copied from viewer.js so behaviour
 * is identical. The only change: `this.viewer` is a CAPTURING SHIM — its addStyle/setStyle record
 * cartoon colours into a buffer and its getModel().selectedAtoms() reads the atom list the frame
 * extracted, so analysis.js works unchanged. After a colour pass the buffer is flushed to the frame
 * as overpaint groups. (Spheres / focus / picking land in M2b/M2c.)
 */
const MolstarViewer = (() => {
  'use strict';
  const NS = 'ufv';
  const FRAME_URL = (typeof chrome !== 'undefined' && chrome.runtime && chrome.runtime.getURL)
    ? chrome.runtime.getURL('viewer-frame.html') : 'viewer-frame.html';

  let iframe = null, readyResolve = null, readyPromise = null, reqSeq = 0;
  const pending = new Map();
  let atoms = [];                 // atom list extracted by the frame (for getModel emulation)
  let cartoonBase = null, cartoonOverrides = [];  // capture buffer for cartoon colouring
  let sphereBuffer = [];          // capture buffer for marker spheres: {chain, resi, color, radius}
  let _partnerSpheres = [];       // multichain: partner-protein spheres by PDB author residue {chain, resi, color, radius, label}
  let _partnerCartoon = [];       // multichain: partner-protein cartoon overrides by PDB author residue {chain, resi, color}
  let _chimericCartoon = [];      // chimeric-partner highlighting overrides {chain, resi, color}
  let cartoonDirty = false, markersDirty = false;
  let lastCartoonJSON = '', lastMarkersJSON = '', lastLabelsJSON = '';  // dedup: skip unchanged sends
  let _clickFn = null, _hoverFn = null, _hoverOutFn = null;  // native-pick handlers (from _bindHover)
  let _focusRegion = null;        // Set of 'chain|pdbResi' shown as sticks in focus → their spheres are excluded
  let _focusHighlightResi = null; // {chain, resi} PDB coords of the clicked residue — rendered as a green CA sphere
  let _showOtherSpheres = true;   // when false, hide all annotation spheres outside the focus region
  // The live object: after Object.assign(StructureViewer, MolstarViewer) the modal writes state onto
  // StructureViewer, not onto `api`. Module-level callbacks (render→_flush, atoms→cache) must read
  // that same object, so init() captures it here.
  let liveRef = null;

  function ensureReady() { if (!readyPromise) readyPromise = new Promise(res => { readyResolve = res; }); return readyPromise; }

  // Lighten an annotation colour so it's legible as TEXT on the native hover tooltip's fixed-dark
  // background. The disease palette (purples/magentas) and the "benign" blue become consequenceColor in
  // disease/variant colouring and were rendered as dark text on dark — unreadable. Colours that are
  // already bright pass through unchanged; dark ones are mixed toward white proportionally.
  function legibleOnDark(hex) {
    const m = /^#?([0-9a-f]{6})$/i.exec(String(hex || ''));
    if (!m) return hex || '#cfd6e0';
    let r = parseInt(m[1].slice(0, 2), 16), g = parseInt(m[1].slice(2, 4), 16), b = parseInt(m[1].slice(4, 6), 16);
    const lum = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255; // perceived luminance 0..1
    if (lum >= 0.5) return '#' + m[1];
    const mix = Math.min(0.72, (0.5 - lum) * 1.5); // how far toward white
    const h = x => Math.round(x + (255 - x) * mix).toString(16).padStart(2, '0');
    return '#' + h(r) + h(g) + h(b);
  }

  // Robust structure-file fetch. AlphaFold/PDBe endpoints occasionally hang, return a transient 5xx, or
  // drop the connection mid-download — leaving "the rest of the structure" missing. We retry with backoff,
  // bound each attempt with a timeout (AbortController), and verify the downloaded body matches
  // Content-Length so a truncated download is treated as a failure (and retried) rather than parsed as a
  // half-built model. Retries bypass the HTTP cache so a corrupt cached entry can't be re-served.
  async function fetchStructureData(url, isBinary, { retries = 3, timeoutMs = 30000 } = {}) {
    let lastErr = null;
    for (let attempt = 0; attempt <= retries; attempt++) {
      const ctrl = (typeof AbortController !== 'undefined') ? new AbortController() : null;
      const timer = ctrl ? setTimeout(() => ctrl.abort(), timeoutMs) : null;
      try {
        const opts = {};
        if (ctrl) opts.signal = ctrl.signal;
        if (attempt > 0) opts.cache = 'reload'; // a retry must not re-read a possibly-truncated cache entry
        const resp = await fetch(url, opts);
        if (!resp.ok) throw new Error('structure fetch ' + resp.status);
        const expected = parseInt(resp.headers.get('content-length') || '0', 10);
        const buf = await resp.arrayBuffer();          // rejects if the connection drops mid-body
        if (expected && buf.byteLength < expected) throw new Error(`truncated structure (${buf.byteLength}/${expected})`);
        if (timer) clearTimeout(timer);
        return isBinary ? new Uint8Array(buf) : new TextDecoder().decode(buf);
      } catch (e) {
        if (timer) clearTimeout(timer);
        lastErr = (e && e.name === 'AbortError') ? new Error('structure fetch timed out') : e;
        if (attempt < retries) await new Promise(r => setTimeout(r, 400 * (attempt + 1))); // 0.4s, 0.8s, 1.2s backoff
      }
    }
    throw lastErr || new Error('structure fetch failed');
  }

  function onMessage(ev) {
    if (!iframe || ev.source !== iframe.contentWindow) return;
    const m = ev.data;
    if (!m || m.ns !== NS) return;
    if (m.evt === 'ready') { if (readyResolve) readyResolve(); return; }
    if (m.evt === 'atoms') { atoms = m.atoms || []; const a = (liveRef || api); a._buildObservedResiCache(); if (a.observedResiCb) { try { a.observedResiCb(); } catch (_) {} } return; }
    if (m.evt === 'pick') {
      if (m.kind === 'click') { if (_clickFn && m.atom) _clickFn(m.atom); }
      else if (m.kind === 'hover') { if (m.atom) { if (_hoverFn) _hoverFn(m.atom, null, {}); } else if (_hoverOutFn) _hoverOutFn(); }
      else if (m.kind === 'dblclick' || m.kind === 'background') { const cb = (liveRef || api).dblClickCb; if (cb) cb(); }
      return;
    }
    if (m.evt === 'result' && m.reqId != null) {
      const p = pending.get(m.reqId); if (!p) return; pending.delete(m.reqId);
      m.ok ? p.resolve(m.value) : p.reject(new Error(m.value && m.value.error || 'command failed'));
    }
  }
  async function send(cmd, args) {
    await ensureReady();
    const reqId = ++reqSeq;
    const p = new Promise((resolve, reject) => pending.set(reqId, { resolve, reject }));
    iframe.contentWindow.postMessage(Object.assign({ ns: NS, cmd, reqId }, args), '*');
    return p;
  }

  // ---- atom filtering: emulate 3Dmol model.selectedAtoms(selector) over the extracted list ----
  function filterAtoms(sel) {
    sel = sel || {};
    const chains = sel.chain != null ? (Array.isArray(sel.chain) ? new Set(sel.chain) : new Set([sel.chain])) : null;
    const resns = sel.resn != null ? new Set(Array.isArray(sel.resn) ? sel.resn : [sel.resn]) : null;
    const out = [];
    for (const a of atoms) {
      if (sel.atom === 'CA' && a.atom !== 'CA') continue;
      if (sel.hetflag === true && !a.hetflag) continue;
      if (sel.hetflag === false && a.hetflag) continue;
      if (chains && !chains.has(a.chain)) continue;
      if (sel.resi != null && a.resi !== sel.resi) continue;
      if (resns && !resns.has(a.resn)) continue;
      out.push(a);
    }
    return out;
  }

  function inferFormat(url, structure) {
    const fmt = structure && (structure.format || structure.fileFormat);
    if (fmt) return fmt === 'cif' || fmt === 'mmcif' ? 'mmcif' : fmt;
    const u = (url || '').toLowerCase();
    if (u.endsWith('.bcif')) return { format: 'bcif', isBinary: true };
    if (u.endsWith('.pdb') || u.endsWith('.ent')) return 'pdb';
    return 'mmcif';
  }

  const api = {
    // ---- state mirrored from StructureViewer ----
    viewer: null,
    currentStructure: null, currentPdbText: null, currentFormat: null,
    nearbyDistance: 5,       // Å — focus neighbourhood radius (configurable via residue-report slider)
    activeColoringMode: 'default', _lastColoringContext: {},
    _observedResi: null, _observedResiByChain: null,
    _selectedResi: null, _inFocusMode: false, _focusState: null,
    excludeIons: false, showLigands: true,
    hiddenLigands: new Set(),   // 'chain|resi' keys of individually-hidden ligands
    ION_CODES: new Set(['NA','K','LI','RB','CS','MG','CA','SR','BA','ZN','FE','FE2','MN','MN3','CU','CU1','CO','NI','CD','HG','CL','BR','IOD','FLO','F','AL','PT','AU','AG','PB','SO4','PO4']),
    clickCb: null, hoverCb: null, dblClickCb: null, ligandClickCb: null, partnerClickCb: null,

    // ---- lifecycle ----
    init(container) {
      liveRef = this;                 // the object the modal actually drives (StructureViewer post-swap)
      if (iframe) return;
      ensureReady();
      window.addEventListener('message', onMessage);
      iframe = document.createElement('iframe');
      iframe.src = FRAME_URL;
      iframe.setAttribute('allow', 'xr-spatial-tracking; fullscreen; clipboard-write');
      iframe.style.cssText = 'border:0;width:100%;height:100%;display:block;background:#0c111b;';
      container.appendChild(iframe);
      // Capturing shim: cartoon → overpaint buffer; spheres → marker buffer; render() flushes both.
      // getModel().selectedAtoms() reads the frame-extracted atom list so analysis.js works unchanged.
      this.viewer = {
        getModel: () => ({ selectedAtoms: (sel) => filterAtoms(sel) }),
        addStyle: (selv, style) => {
          if (!style || !selv || selv.resi == null) return;
          if (style.cartoon && style.cartoon.color != null && selv.atom == null) {
            cartoonOverrides.push({ chain: selv.chain != null ? selv.chain : null, resi: selv.resi, color: style.cartoon.color }); cartoonDirty = true;
          } else if (style.sphere && style.sphere.color != null) {
            sphereBuffer.push({ chain: selv.chain != null ? selv.chain : null, resi: selv.resi, color: style.sphere.color, radius: style.sphere.radius || 1.6 }); markersDirty = true;
          }
        },
        setStyle: (selv, style) => {
          if (selv && selv.resn != null) return;                 // ligand hide — M2b ligand layer
          if (!selv || Object.keys(selv).length === 0) {
            if (!style || Object.keys(style).length === 0) { cartoonBase = null; cartoonOverrides = []; sphereBuffer = []; cartoonDirty = true; markersDirty = true; return; }
            if (style.cartoon) {
              cartoonDirty = true;
              if (typeof style.cartoon.colorfunc === 'function') {
                for (const a of filterAtoms({ atom: 'CA' })) cartoonOverrides.push({ chain: a.chain, resi: a.resi, color: style.cartoon.colorfunc(a) });
              } else if (style.cartoon.color != null) { cartoonBase = style.cartoon.color; }
            }
          }
        },
        setHoverable: (_sel, _on, fnHover, fnOut) => { _hoverFn = fnHover; _hoverOutFn = fnOut; },
        setClickable: (_sel, _on, fnClick) => { _clickFn = fnClick; },
        removeAllShapes: () => {}, removeAllLabels: () => {}, render: () => (liveRef || api)._flush(),
        addLabel: () => null, removeLabel: () => {}, addLine: () => {}, removeShape: () => {},
        zoomTo: (_s, dur) => send('reset', { duration: dur ?? 0 }).catch(() => {}),
        zoom: (_f, dur) => send('reset', { duration: dur ?? 0 }).catch(() => {}),
        resize: () => send('resize').catch(() => {}),
        pngURI: () => null,
        setBackgroundColor: (color) => {
          if (iframe) iframe.style.background = color;
          send('background', { color }).catch(() => {});
        },
      };
    },

    async loadStructure(url, structure = null) {
      await ensureReady();
      this.currentStructure = structure;
      this._observedResi = this._observedResiByChain = null; this._waterKeysCache = null;
      atoms = []; _partnerSpheres = []; _partnerCartoon = []; _chimericCartoon = []; _focusRegion = null; _focusHighlightResi = null; this._lastCamTarget = null; this._lastCamRadius = null; this._pocketShown = false; lastCartoonJSON = lastMarkersJSON = lastLabelsJSON = '';
      // Reset all focus / ligand-visibility state BEFORE the new structure loads. Otherwise a focus left
      // active from the previous structure carries over — in particular a hide-all _focusNearbyLigands
      // (residue focus on a transplant-crowded model) would hide every ligand on the freshly loaded one.
      this._selectedResi = null; this._inFocusMode = false; this._focusState = null;
      this._focusNearbyLigands = null; this.hiddenLigands = new Set(); this._lastLigandJSON = '';
      _showOtherSpheres = true;
      const f = inferFormat(url, structure);
      const format = typeof f === 'string' ? f : f.format;
      const isBinary = typeof f === 'string' ? false : f.isBinary;
      // Retry/timeout/truncation-guarded fetch — a transient drop no longer fails the whole load.
      const data = await fetchStructureData(url, isBinary);
      this.currentPdbText = isBinary ? null : data;
      this.currentFormat = format === 'pdb' ? 'pdb' : 'mmcif';
      await send('loadData', { data, format, isBinary });   // frame replies after load; 'atoms' event follows
      return true;
    },

    async clearModel() {
      this.currentStructure = this.currentPdbText = this.currentFormat = null;
      this._observedResi = this._observedResiByChain = null; this._waterKeysCache = null; atoms = [];
      this._selectedResi = null; this._inFocusMode = false; this._focusState = null;
      this.hiddenLigands = new Set(); this._lastLigandJSON = '';
      this._focusNearbyLigands = null;
      if (iframe) { try { await send('clear'); } catch (_) {} }
    },

    resize() { send('resize').catch(() => {}); },
    resetView() { this.clearPocket(); this._selectedResi = null; this._inFocusMode = false; this._focusState = null; this._focusLabels = null; _focusRegion = null; _focusHighlightResi = null; this._lastCamTarget = null; this._lastCamRadius = null; _showOtherSpheres = true; this._focusNearbyLigands = null; send('unfocus').catch(() => {}); markersDirty = true; this._flush(true); this._drawLigands(); this._applyMarkerVisibility(); /* clear the toggle-off marker transparency → overview shows all spheres */ },
    // Re-issue ONLY the camera zoom for the active focus target at the current nearbyDistance. Used by the
    // Nearby slider so dragging the radius re-frames the pocket without rebuilding the focus sticks.
    rezoomFocus() { const t = this._lastCamTarget; if (t) send('camFocus', { chain: t.chain, resi: t.resi, radius: this._lastCamRadius || (this.nearbyDistance || 5) + 3 }).catch(() => {}); },
    // Highlight a predicted pocket as a translucent molecular surface over its residues (UniProt positions
    // → PDB) and frame it. clearPocket() removes it; any subsequent focus/reset also clears it.
    showPocket(positions, chain = null, annotatedResidues = null) {
      const residues = []; const seen = new Set();
      for (const p of (positions || [])) {
        const sel = this.residueSelectorForChain(p, chain);
        if (!sel || sel.resi === -999999) continue;
        const ch = sel.chain != null ? sel.chain : (chain != null ? chain : null);
        const k = (ch == null ? '' : ch) + '|' + sel.resi; if (seen.has(k)) continue; seen.add(k);
        const meta = annotatedResidues && annotatedResidues.get(p); // per-residue annotation colour for sticks
        residues.push({ chain: ch, resi: sel.resi, color: (meta && meta.color) || null });
      }
      if (!residues.length) return false;
      this._pocketShown = true;
      send('showPocket', { residues, color: '#26c6da' }).catch(() => {});
      return true;
    },
    clearPocket() { if (this._pocketShown) { this._pocketShown = false; send('clearPocket').catch(() => {}); } },
    // Frame the camera on a set of residues (e.g. a whole domain range) — UniProt positions → PDB — without
    // entering focus mode or drawing sticks. Used by the range-feature magnifier ("zoom to a–b").
    frameResidues(positions, chain = null) {
      const residues = []; const seen = new Set();
      for (const p of (positions || [])) {
        const sel = this.residueSelectorForChain(p, chain);
        if (!sel || sel.resi === -999999) continue;
        const ch = sel.chain != null ? sel.chain : (chain != null ? chain : null);
        const k = (ch == null ? '' : ch) + '|' + sel.resi; if (seen.has(k)) continue; seen.add(k);
        residues.push({ chain: ch, resi: sel.resi });
      }
      if (residues.length) send('frameResidues', { residues }).catch(() => {});
    },
    // Toggle annotation spheres OUTSIDE the focus region. The non-focus spheres are ALWAYS drawn (stable
    // marker geometry); the toggle only flips a full-transparency layer on the marker reps — no geometry
    // change, so no flash and no camera auto-fit bounce. No re-focus, no _flush.
    setOtherSpheresVisible(visible) { _showOtherSpheres = visible; this._applyMarkerVisibility(); },
    // Hide the non-focus marker spheres (transparency) only while focused with the toggle off.
    _applyMarkerVisibility() { send('markersHidden', { hidden: !!_focusRegion && !_showOtherSpheres }).catch(() => {}); },
    async screenshot() {
      try {
        const uri = await send('screenshot'); if (!uri) return;
        const a = document.createElement('a');
        a.href = uri; a.download = (this.currentStructure && (this.currentStructure.id || this.currentStructure.pdbId) || 'structure') + '.png';
        a.click();
      } catch (_) {}
    },

    // ---- residue mapping (copied verbatim from viewer.js) ----
    _resiToPdb(uniprotResi, r, chainId = null) {
      const map = (chainId != null && this.currentStructure?.chainSeqresToAuthor?.[chainId]) || this.currentStructure?.seqresToAuthor;
      if (map) { const seqresPos = (r.seqresStart ?? r.pdbStart) + (uniprotResi - r.uniprotStart); return map.get(seqresPos) ?? seqresPos; }
      return r.pdbStart + (uniprotResi - r.uniprotStart);
    },
    _buildObservedResiCache() {
      this._observedResi = null; this._observedResiByChain = null;
      const structure = this.currentStructure;
      if (!structure || (structure.source === 'AlphaFold' && !structure.isoform) || !structure.mappedRanges?.length) return;
      const model = this.viewer?.getModel(); if (!model) return;
      if (structure.chainIds?.length > 1) {
        this._observedResiByChain = {};
        structure.chainIds.forEach(chain => { const a = model.selectedAtoms({ chain, atom: 'CA' }); this._observedResiByChain[chain] = a?.length ? new Set(a.map(x => x.resi)) : null; });
        this._observedResi = this._observedResiByChain[structure.chainId] || null;
      } else {
        const chain = structure.chainId;
        const a = model.selectedAtoms({ ...(chain ? { chain } : {}), atom: 'CA' });
        if (a?.length) this._observedResi = new Set(a.map(x => x.resi));
      }
    },
    // Is a UniProt position actually resolved (has coordinates) in the loaded structure? Full-coverage
    // models (canonical AlphaFold) resolve everything → true. Experimental / mapped structures only
    // resolve the observed residues; an annotation on an unresolved residue can't be drawn, so the UI
    // greys it out. Returns true when coverage is unknown (no observed cache yet) to avoid false greys.
    isResolved(uniPos, chain = null) {
      const structure = this.currentStructure;
      if (!structure || (structure.source === 'AlphaFold' && !structure.isoform) || !structure.mappedRanges?.length) return true;
      if (!this._observedResi && !this._observedResiByChain) return true; // cache not built yet
      const sel = this.residueSelectorForChain(uniPos, chain);
      return sel.resi !== -999999;
    },
    residueSelector(resi) {
      const structure = this.currentStructure;
      if (!structure || (structure.source === 'AlphaFold' && !structure.isoform) || !structure.mappedRanges?.length) return { resi };
      for (const r of structure.mappedRanges) {
        if (resi >= r.uniprotStart && resi <= r.uniprotEnd) {
          const pdbResi = this._resiToPdb(resi, r);
          if (this._observedResi && !this._observedResi.has(pdbResi)) return { resi: -999999 };
          return r.chainId ? { chain: r.chainId, resi: pdbResi } : { resi: pdbResi };
        }
      }
      return { resi: -999999 };
    },
    residueSelectorForChain(resi, chain) {
      const structure = this.currentStructure;
      if (chain == null) return this.residueSelector(resi);
      if (!structure || (structure.source === 'AlphaFold' && !structure.isoform)) return { resi };
      const ranges = structure.chainMappings?.[chain] || structure.mappedRanges || [];
      const obs = this._observedResiByChain?.[chain] ?? this._observedResi;
      for (const r of ranges) {
        if (resi >= r.uniprotStart && resi <= r.uniprotEnd) {
          const pdbResi = this._resiToPdb(resi, r, chain);
          if (obs && !obs.has(pdbResi)) return { chain, resi: -999999 };
          return { chain, resi: pdbResi };
        }
      }
      return { chain, resi: -999999 };
    },
    chainAddStyle(chain, resi, style, atomSel = {}) {
      const structure = this.currentStructure;
      if (chain == null) return this.allChainsAddStyle(resi, style, atomSel);
      const ranges = structure?.chainMappings?.[chain] || structure?.mappedRanges || [];
      const obs = this._observedResiByChain?.[chain] ?? this._observedResi;
      const r = ranges.find(mr => resi >= mr.uniprotStart && resi <= mr.uniprotEnd);
      if (!r) return false;
      const pdbResi = this._resiToPdb(resi, r, chain);
      if (obs && !obs.has(pdbResi)) return false;
      this.viewer.addStyle({ chain, resi: pdbResi, ...atomSel }, style);
      return true;
    },
    allChainsAddStyle(resi, style, atomSel = {}) {
      const structure = this.currentStructure;
      if (!structure || (structure.source === 'AlphaFold' && !structure.isoform) || !structure.mappedRanges?.length) { this.viewer.addStyle({ resi, ...atomSel }, style); return true; }
      if (structure.chainIds?.length > 1) {
        let applied = false;
        structure.chainIds.forEach(chain => {
          const ranges = structure.chainMappings?.[chain] || structure.mappedRanges;
          const r = ranges.find(mr => resi >= mr.uniprotStart && resi <= mr.uniprotEnd);
          if (!r) return;
          const chainPdbResi = this._resiToPdb(resi, r, chain);
          const chainObs = this._observedResiByChain?.[chain] ?? this._observedResi;
          if (chainObs && !chainObs.has(chainPdbResi)) return;
          this.viewer.addStyle({ chain, resi: chainPdbResi, ...atomSel }, style); applied = true;
        });
        return applied;
      }
      for (const r of structure.mappedRanges) {
        if (resi >= r.uniprotStart && resi <= r.uniprotEnd) {
          const pdbResi = this._resiToPdb(resi, r);
          if (this._observedResi && !this._observedResi.has(pdbResi)) return false;
          const chainSel = r.chainId ? { chain: r.chainId } : {};
          this.viewer.addStyle({ ...chainSel, resi: pdbResi, ...atomSel }, style);
          return true;
        }
      }
      return false;
    },
    mappedResidues() { return this.mappedResiduesForChain(null); },
    mappedResiduesForChain(chain) {
      const s = this.currentStructure;
      const ranges = (chain != null && s?.chainMappings?.[chain]) || s?.mappedRanges;
      if (!ranges?.length) return null;
      const obs = (chain != null ? (this._observedResiByChain?.[chain] ?? null) : this._observedResi);
      const residues = [];
      ranges.forEach(r => { for (let i = r.uniprotStart; i <= r.uniprotEnd; i++) { if (!obs) residues.push(i); else { const pdbResi = this._resiToPdb(i, r, chain); if (obs.has(pdbResi)) residues.push(i); } } });
      return residues.length > 0 ? residues : null;
    },
    _reverseResidueMapForChain(chain = null) {
      const out = new Map(); const s = this.currentStructure; if (!s) return out;
      // A partner chain (a DIFFERENT protein in a complex) must NOT inherit our protein's mapping — doing
      // so MIRRORS our annotations onto it (wrong data on cartoon/sticks/nearby/clicks). Only our chains
      // (chainIds) get the mappedRanges fallback; everything else returns an empty map (no UniProt mapping).
      if (chain != null && Array.isArray(s.chainIds) && s.chainIds.length && !s.chainIds.includes(chain)) return out;
      const ranges = (chain != null && s.chainMappings?.[chain]) || s.mappedRanges || [];
      const map = (chain != null && s.chainSeqresToAuthor?.[chain]) || s.seqresToAuthor || null;
      const obs = (chain != null ? (this._observedResiByChain?.[chain] ?? this._observedResi) : this._observedResi);
      ranges.forEach(r => {
        if (map && r.seqresStart != null) {
          const seqresEnd = r.seqresStart + (r.uniprotEnd - r.uniprotStart);
          map.forEach((author, seqres) => { if (seqres < r.seqresStart || seqres > seqresEnd) return; if (obs && !obs.has(author)) return; out.set(author, r.uniprotStart + (seqres - r.seqresStart)); });
        } else { for (let pdb = r.pdbStart; pdb <= r.pdbEnd; pdb++) { if (obs && !obs.has(pdb)) continue; out.set(pdb, r.uniprotStart + (pdb - r.pdbStart)); } }
      });
      return out;
    },
    residueGeometry() {
      const model = this.viewer?.getModel?.(); if (!model) return [];
      const s = this.currentStructure;
      const ourChains = s?.chainIds?.length ? new Set(s.chainIds) : (s?.chainId ? new Set([s.chainId]) : null);
      const reverseCache = new Map();
      const toUni = (chain, resi) => { if (!reverseCache.has(chain)) reverseCache.set(chain, this._reverseResidueMapForChain(chain)); return reverseCache.get(chain).get(resi); };
      const seen = new Set(); const out = [];
      (model.selectedAtoms({ atom: 'CA' }) || []).forEach(a => {
        const key = `${a.chain}|${a.resi}`; if (seen.has(key)) return; seen.add(key);
        const isOurs = !ourChains || ourChains.has(a.chain);
        const uni = isOurs ? toUni(a.chain, a.resi) : null;
        out.push({ uniPos: uni != null ? uni : null, chain: a.chain ?? null, resi: a.resi, ca: { x: a.x, y: a.y, z: a.z }, b: a.b });
      });
      return out;
    },
    computeActualCoverage(structure) {
      if (!this.viewer || !structure?.mappedRanges?.length) return null;
      const model = this.viewer.getModel(); if (!model) return null;
      const chain = structure.chainId;
      const a = model.selectedAtoms(chain ? { chain } : {}); if (!a || !a.length) return null;
      const observed = new Set(a.map(x => x.resi)); let count = 0;
      structure.mappedRanges.forEach(r => { for (let uni = r.uniprotStart; uni <= r.uniprotEnd; uni++) { const pdb = this._resiToPdb(uni, r); if (observed.has(pdb)) count++; } });
      return count;
    },

    // ---- cartoon colouring (copied from viewer.js; flushes the capture buffer to the frame) ----
    applyCartoonColoring(mode = 'default', context = {}, defer = false) {
      if (!this.viewer) return;
      this.activeColoringMode = mode; this._lastColoringContext = context;
      this._inFocusMode = false; this._focusState = null;
      this.viewer.setStyle({}, {});                       // reset capture buffer
      const base = { opacity: 0.82, thickness: 0.25, ribbonWidth: 0.6 };
      this._applyModeStyles(mode, context, base);
      this._drawLigands();
      if (!defer) this.viewer.render();
    },
    _applyModeStyles(mode, context, base) {
      if (mode === 'plddt') {
        this.viewer.setStyle({}, { cartoon: { ...base, colorfunc: atom => { const b = atom.b; if (b >= 90) return '#0053d6'; if (b >= 70) return '#65cbf3'; if (b >= 50) return '#ffdb13'; return '#ff7d45'; } } });
      } else if (mode === 'bfactor') {
        this.viewer.setStyle({}, { cartoon: { ...base, colorfunc: atom => {
          const v = Math.max(0, Math.min(100, atom.b));
          if (v <= 50) { const t = v / 50; return `rgb(${Math.round(49 + 198 * t)},${Math.round(54 + 193 * t)},${Math.round(149 + 98 * t)})`; }
          const t = (v - 50) / 50; return `rgb(${Math.round(247 - 32 * t)},${Math.round(247 - 199 * t)},${Math.round(247 - 208 * t)})`;
        } } });
      } else if (mode === 'hotspots') {
        this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
        this._applyTierColoring(context.hotspots, context.hotspotsByChain, { strong: '#b71c1c', moderate: '#e64a19', weak: '#ffa726' }, base);
      } else if (mode === 'distantContacts') {
        this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
        this._applyTierColoring(context.distantContacts, context.distantContactsByChain, { strong: '#6a1b9a', moderate: '#ab47bc' }, base);
      } else if (mode === 'alphaMissense') {
        this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
        (context.alphaMissense || new Map()).forEach((d, pos) => { const avg = d.avg; const color = avg >= 0.78 ? '#b71c1c' : avg >= 0.564 ? '#e06666' : avg >= 0.34 ? '#b9c2cf' : '#3d85c8'; this.allChainsAddStyle(pos, { cartoon: { ...base, color } }); });
      } else if (mode === 'residueBurden') {
        this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
        (context.residueBurden || new Set()).forEach(pos => this.allChainsAddStyle(pos, { cartoon: { ...base, color: '#e65100' } }));
      } else if (mode === 'topology' || mode === 'domains') {
        const posColor = context.topologyByPos instanceof Map ? context.topologyByPos : context.domainByPos instanceof Map ? context.domainByPos : new Map();
        // Domains use a near-white base (off-white) for non-domain residues so the coloured domains stand
        // out without the blue-grey tint; topology keeps the blue-grey base.
        const baseCol = mode === 'domains' ? '#ededed' : '#b9c2cf';
        const s = this.currentStructure;
        const isAF = !s || (s.source === 'AlphaFold' && !s.isoform) || !s.mappedRanges?.length;
        const ourChains = s?.chainIds?.length ? new Set(s.chainIds) : (s?.chainId ? new Set([s.chainId]) : null);
        const reverseCache = new Map();
        const toUni = (chain, resi) => { if (!reverseCache.has(chain)) reverseCache.set(chain, this._reverseResidueMapForChain(chain)); return reverseCache.get(chain).get(resi); };
        this.viewer.setStyle({}, { cartoon: { ...base, colorfunc: atom => { if (!isAF && ourChains && atom.chain != null && !ourChains.has(atom.chain)) return baseCol; const uni = isAF ? atom.resi : toUni(atom.chain, atom.resi); return posColor.get(uni) || baseCol; } } });
      } else if (mode === 'prism') {
        this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
        const catColors = { pocket: '#00897b', exposed: '#8e24aa' };
        if (context.pocketByPos instanceof Map) context.pocketByPos.forEach((info, pos) => this.allChainsAddStyle(pos, { cartoon: { ...base, color: catColors[info.cat] || '#00897b' } }));
      } else {
        this.viewer.setStyle({}, { cartoon: { ...base, color: '#d0d0d0' } });
      }
    },
    _applyTierColoring(merged, byChain, tierColors, base) {
      const structure = this.currentStructure;
      if (structure?.chainIds?.length > 1 && byChain) { byChain.forEach((tierMap, chain) => { tierMap.forEach((tier, pos) => { this.chainAddStyle(chain, pos, { cartoon: { ...base, color: tierColors[tier] || '#b9c2cf' } }); }); }); return; }
      (merged || new Map()).forEach((tier, pos) => { this.allChainsAddStyle(pos, { cartoon: { ...base, color: tierColors[tier] || '#b9c2cf' } }); });
    },
    // Flush both capture channels to the frame (cartoon overpaint + marker spheres + native labels).
    // Each channel is sent only when its payload actually changed, so toggling a marker doesn't
    // rebuild the cartoon (which otherwise flickers) and vice-versa.
    // skipCamRestore: pass true when a camera animation was just requested (focus/unfocus) so the
    // frame's setMarkers doesn't snapshot+restore the pre-animation camera position, which would
    // cancel the zoom partway through.
    _flush(skipCamRestore = false) {
      if (cartoonDirty) {
        // Resolve to ONE colour per residue (last write wins) so a partner-protein override reliably beats
        // the mode colourfunc's neutral colour for that same partner residue, regardless of group order.
        const resColor = new Map();
        for (const o of cartoonOverrides) resColor.set((o.chain == null ? '' : o.chain) + '|' + o.resi, o);
        // Multichain: partner-protein cartoon overrides (disease colours on partner chains), applied LAST.
        for (const pc of _partnerCartoon) resColor.set((pc.chain == null ? '' : pc.chain) + '|' + pc.resi, pc);
        // Chimeric-partner highlighting: overrides the mode colour for unmapped residues in our chain.
        for (const cc of _chimericCartoon) resColor.set((cc.chain == null ? '' : cc.chain) + '|' + cc.resi, cc);
        const byColor = new Map();
        for (const o of resColor.values()) { if (!byColor.has(o.color)) byColor.set(o.color, []); byColor.get(o.color).push([o.chain, o.resi]); }
        const payload = { base: cartoonBase || '#d0d0d0', groups: [...byColor.entries()].map(([color, residues]) => ({ color, residues })) };
        const j = JSON.stringify(payload);
        if (j !== lastCartoonJSON) { lastCartoonJSON = j; send('setCartoon', payload).catch(() => {}); }
        cartoonDirty = false;
      }
      if (markersDirty) {
        // Dedup overlapping annotations to ONE sphere per residue (last push wins: site → PTM → variant),
        // so a residue carrying two annotations shows a single deterministic colour and toggling the
        // winning annotation re-reveals the next one's colour instead of leaving a stale stacked sphere.
        const byPos = new Map();
        for (const s of sphereBuffer) byPos.set((s.chain == null ? '' : s.chain) + '|' + s.resi, s);
        const byColor = new Map();
        for (const s of byPos.values()) {
          const key = (s.chain == null ? '' : s.chain) + '|' + s.resi;
          if (_focusRegion && _focusRegion.has(key)) continue;          // shown as sticks in focus
          // NB: the _showOtherSpheres toggle is NOT applied here. The non-focus spheres are always part
          // of the marker geometry (stable bounds → no camera auto-fit). The toggle hides them via a
          // transparency command (markersHidden) instead, which doesn't change geometry → no flash/bounce.
          if (!byColor.has(s.color)) byColor.set(s.color, { radius: s.radius, residues: [] }); byColor.get(s.color).residues.push([s.chain, s.resi]);
        }
        // Multichain: partner-protein spheres are placed by PDB author residue on the partner chain
        // directly (they don't pass through our UniProt mapping). Same colour grouping.
        for (const ps of _partnerSpheres) {
          const pk = (ps.chain == null ? '' : ps.chain) + '|' + ps.resi;
          if (_focusRegion && _focusRegion.has(pk)) continue;            // partner neighbour shown as a focus stick
          if (!byColor.has(ps.color)) byColor.set(ps.color, { radius: ps.radius || 1.6, residues: [] });
          byColor.get(ps.color).residues.push([ps.chain, ps.resi]);
        }
        // Focused-residue highlight: a green CA sphere marking the clicked residue so it stands out
        // from the surrounding neighbour sticks. Not filtered by _focusRegion — intentionally visible
        // on top of the ball-and-stick rep. Radius 0.7 Å sits just outside the B&S atom spheres.
        if (_focusHighlightResi) {
          const fh = _focusHighlightResi;
          const FH_COLOR = '#48c78e'; // Mol*-ish teal green
          if (!byColor.has(FH_COLOR)) byColor.set(FH_COLOR, { radius: 0.7, residues: [] });
          byColor.get(FH_COLOR).residues.push([fh.chain, fh.resi]);
        }
        const payload = { groups: [...byColor.entries()].map(([color, g]) => ({ color, radius: g.radius, residues: g.residues })) };
        const j = JSON.stringify(payload);
        if (j !== lastMarkersJSON) { lastMarkersJSON = j; send('setMarkers', { ...payload, skipCamRestore }).catch(() => {}); }
        // Native hover label map: chain|authResi → annotation text (reverse-map each sphere to its uniPos).
        const labels = {};
        for (const s of sphereBuffer) {
          const uni = (this._reverseResidueMapForChain(s.chain).get(s.resi)) || s.resi;
          const html = this._richLabel(this._hoverMap && this._hoverMap.get(uni));
          if (html) labels[(s.chain == null ? '' : s.chain) + '|' + s.resi] = html;
        }
        // Partner spheres carry their own pre-rendered label (different protein, not in our hover map).
        for (const ps of _partnerSpheres) { if (ps.label) labels[(ps.chain == null ? '' : ps.chain) + '|' + ps.resi] = ps.label; }
        // Focus mode: basic residue labels for every focus side chain, UNDER any richer annotation label.
        if (this._inFocusMode && this._focusLabels) for (const k in this._focusLabels) { if (!labels[k]) labels[k] = this._focusLabels[k]; }
        const lj = JSON.stringify(labels);
        if (lj !== lastLabelsJSON) { lastLabelsJSON = lj; send('setLabels', { map: labels }).catch(() => {}); }
        markersDirty = false;
      }
    },
    // Rich tooltip body HTML, matching the modal's residueSummary(): variants as bold WT-pos-MUT +
    // consequence (coloured), PTMs/sites as category/description. Sent to the frame's tooltip.
    _richLabel(d) {
      if (!d) return null;
      const esc = (x) => String(x == null ? '' : x).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
      if (d.variants && d.variants.length) {
        return d.variants.slice(0, 4).map(v =>
          `<strong>${esc(v.wildType)}${v.position}${esc(v.mutant)}</strong> <span style="color:${esc(legibleOnDark(v.consequenceColor || '#9e9e9e'))}">${esc(v.consequence)}</span>`
        ).join(' &nbsp;|&nbsp; ');
      }
      if (d.isSite) return esc(d.description || 'Site');
      if (d.category) return esc(d.description ? `${d.category}: ${d.description}` : d.category);
      return esc(d.description || d.label || '') || null;
    },
    // Push current ligand visibility to the frame. Ligands render via Mol*'s default preset, so we
    // "hide" them with a transparency layer: all of them when showLigands is off, otherwise the
    // individually-toggled hiddenLigands set. Deduped so repeated applyMode calls don't re-send.
    _lastLigandJSON: '',
    _focusNearbyLigands: null,   // while focused: the ONLY ligand keys to keep visible (hide the rest)
    _drawLigands() {
      if (!this.viewer) return;
      let hidden;
      if (!this.showLigands) {
        hidden = this.enumerateLigands().map(l => [l.chain, l.resi]);
      } else {
        const keys = new Set(this.hiddenLigands || []);
        // Exclude water & ions toggle: hide every ion ligand AND every water (HOH/WAT). Water isn't in
        // the ligand list, so its residue keys are enumerated separately (cached per structure).
        if (this.excludeIons) {
          for (const l of this.enumerateLigands()) {
            if (this.ION_CODES.has(l.resn)) keys.add((l.chain == null ? '' : l.chain) + '|' + l.resi);
          }
          for (const k of this._waterKeys()) keys.add(k);
        }
        // Focus mode: hide every ligand NOT in the keep-set (so distant AlphaFill ligands don't
        // clutter the zoom). Kept ligands are shown even if the user had hidden them.
        if (this._focusNearbyLigands) {
          for (const l of this.enumerateLigands()) {
            const k = (l.chain == null ? '' : l.chain) + '|' + l.resi;
            if (!this._focusNearbyLigands.has(k)) keys.add(k);
          }
          for (const k of this._focusNearbyLigands) keys.delete(k);
        }
        hidden = [...keys].map(k => { const i = k.indexOf('|'); const c = k.slice(0, i); return [c === '' ? null : c, Number(k.slice(i + 1))]; });
      }
      const j = JSON.stringify(hidden);
      if (j === this._lastLigandJSON) return;
      this._lastLigandJSON = j;
      send('setLigandVisibility', { hidden }).catch(() => {});
    },
    // Multichain: set the partner-protein annotation spheres (placed by PDB author residue on partner
    // chains). Each: { chain, resi, color, radius?, label? }. Re-flushes the marker channel.
    setPartnerSpheres(spheres) {
      _partnerSpheres = Array.isArray(spheres) ? spheres : [];
      markersDirty = true;
      this._flush(true);
    },
    // Multichain: tint ENTIRE partner chains a single identity colour (the "Partners" toggle) so the user
    // can see which subunits are other proteins. NOT disease/PTM colours. Empty list / no colour clears it.
    tintPartnerChains(chainIds, color) {
      _partnerCartoon = [];
      if (chainIds && chainIds.length && color) {
        const set = new Set(chainIds); const seen = new Set();
        for (const a of atoms) {
          if (a.hetflag || a.chain == null || !set.has(a.chain)) continue;
          const k = a.chain + '|' + a.resi; if (seen.has(k)) continue; seen.add(k);
          _partnerCartoon.push({ chain: a.chain, resi: a.resi, color });
        }
      }
      cartoonDirty = true; lastCartoonJSON = '';
      this._flush(true);
    },
    // Chimeric-partner highlighting: override the backbone colour of unmapped residues in our chain
    // so they read as a distinct group from the annotation-coloured residues. Persists through
    // colour-mode changes (applied after partner overrides in _flush). Pass [] to clear.
    setChimericHighlight(residues, color = '#8fa5c0') {
      _chimericCartoon = residues.map(r => ({ chain: r.chain, resi: r.resi, color }));
      cartoonDirty = true;
      this._flush();
    },
    clearChimericHighlight() {
      if (!_chimericCartoon.length) return;
      _chimericCartoon = [];
      cartoonDirty = true;
      this._flush();
    },
    // Multichain: focus a PARTNER residue with its surroundings (normal Mol* behaviour) — zoom + sticks for
    // the residue and its neighbourhood, by PDB author residue (bypassing our UniProt mapping).
    focusPartnerResidue(chain, resi) {
      this.clearPocket();
      const wasFocused = !!this._inFocusMode;
      this._selectedResi = null; this._inFocusMode = true;
      const targets = atoms.filter(a => !a.hetflag && a.resi === resi && a.chain === chain);
      let neighbors = [], region = null;
      try { if (targets.length) { const fn = this._focusNeighbourhood(targets, chain, resi); neighbors = fn.nearPdb || []; region = fn.region; } } catch (_) {}
      // Set the exclusion region + flush markers so the partner/our spheres at the focus residues are HIDDEN
      // (shown as sticks instead) — without this the partner spheres overlaid the focus sticks.
      _focusRegion = region; _focusHighlightResi = { chain, resi }; this._lastCamTarget = { chain, resi };
      this._lastCamRadius = this._pocketCamRadius(targets, neighbors);
      // Colour the partner focus sticks with OUR disease colours (via _partnerColorMap) — no our-protein map.
      const annotations = this._stickAnnotations(neighbors, null);
      send('focus', { chain, resi, neighbors, annotations }).catch(() => {});
      markersDirty = true; this._flush(true); this._applyMarkerVisibility();
      send('camFocus', { chain, resi, radius: this._lastCamRadius, panOnly: wasFocused }).catch(() => {}); // camera LAST (uncontended)
    },
    // Water residue keys ('chain|resi') for the loaded structure (cached — there can be thousands).
    _waterKeysCache: null,
    _waterKeys() {
      if (this._waterKeysCache) return this._waterKeysCache;
      const out = []; const seen = new Set();
      for (const a of atoms) {
        if (!a.hetflag || (a.resn !== 'HOH' && a.resn !== 'WAT')) continue;
        const k = (a.chain == null ? '' : a.chain) + '|' + a.resi;
        if (!seen.has(k)) { seen.add(k); out.push(k); }
      }
      this._waterKeysCache = out;
      return out;
    },
    // Toggle a single ligand's visibility (modal ligand list). visible=false hides just this copy.
    setLigandVisible(chain, resi, visible) {
      const key = (chain == null ? '' : chain) + '|' + resi;
      if (visible) this.hiddenLigands.delete(key); else this.hiddenLigands.add(key);
      this._drawLigands();
    },
    isLigandVisible(chain, resi) {
      if (!this.showLigands) return false;
      return !this.hiddenLigands.has((chain == null ? '' : chain) + '|' + resi);
    },

    // ---- annotation sphere layers (copied from viewer.js; render() flushes to setMarkers) ----
    _drawSiteSpheres(sites, hoverMap, active) {
      if (!sites || !sites.length) return;
      sites.forEach(site => {
        const positions = [site.position, site.endPosition].filter((p, i, a) => p && a.indexOf(p) === i);
        positions.forEach(pos => {
          if (hoverMap.has(pos)) return;
          if (!this.allChainsAddStyle(pos, { sphere: { radius: 1.8, color: site.color, opacity: 0.92 } }, { atom: 'CA' })) return;
          hoverMap.set(pos, { position: pos, color: site.color, isSite: true, description: site.description, category: 'Site' });
          active.set(pos, site.color);
        });
      });
    },
    showPTMs(ptms, ptmGroups, sites = [], extras = []) {
      this._selectedResi = null; this._lastRender = () => this.showPTMs(ptms, ptmGroups, sites, extras);
      const hoverMap = new Map(); const active = new Map(); let count = 0;
      (ptms || []).forEach(ptm => {
        // item.visible is authoritative (group toggles sync all items). This lets an INDIVIDUAL PTM be
        // shown even when its group's master checkbox is off — e.g. in the Structure/Subcellular windows.
        const g = ptmGroups[ptm.category]; if (!g || ptm.visible === false) return;
        if (!this.allChainsAddStyle(ptm.position, { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }, { atom: 'CA' })) return;
        hoverMap.set(ptm.position, ptm); active.set(ptm.position, ptm.color); count++;
        if (ptm.endPosition && ptm.endPosition !== ptm.position) {
          if (this.allChainsAddStyle(ptm.endPosition, { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }, { atom: 'CA' })) active.set(ptm.endPosition, ptm.color);
          if (!hoverMap.has(ptm.endPosition)) hoverMap.set(ptm.endPosition, { ...ptm, position: ptm.endPosition, description: ptm.description + ` (bonded to ${ptm.position})` });
        }
      });
      this._drawSiteSpheres(sites, hoverMap, active);
      (extras || []).forEach(sp => { if (sp.position == null) return; if (this.allChainsAddStyle(sp.position, { sphere: { radius: 1.8, color: sp.color, opacity: 0.92 } }, { atom: 'CA' })) { active.set(sp.position, sp.color); if (sp.hover && !hoverMap.has(sp.position)) hoverMap.set(sp.position, sp.hover); } });
      this._activeSpheres = active; this._bindHover(hoverMap, 'ptm'); this.viewer.render();
      return count;
    },
    showAnnotationSpheres(spheres) {
      if (!this.viewer) return 0;
      this._selectedResi = null; this._lastRender = () => this.showAnnotationSpheres(spheres);
      const hoverMap = new Map(); const active = new Map(); let count = 0;
      (spheres || []).forEach(sp => {
        if (sp.position == null) return;
        if (!this.allChainsAddStyle(sp.position, { sphere: { radius: 1.8, color: sp.color, opacity: 0.92 } }, { atom: 'CA' })) return;
        active.set(sp.position, sp.color); if (sp.hover) hoverMap.set(sp.position, sp.hover); count++;
        if (sp.endPosition && sp.endPosition !== sp.position) {
          if (this.allChainsAddStyle(sp.endPosition, { sphere: { radius: 1.8, color: sp.color, opacity: 0.92 } }, { atom: 'CA' })) active.set(sp.endPosition, sp.color);
          if (sp.hover && !hoverMap.has(sp.endPosition)) hoverMap.set(sp.endPosition, { ...sp.hover, position: sp.endPosition });
        }
      });
      this._activeSpheres = active; this._bindHover(hoverMap, 'feature'); this.viewer.render();
      return count;
    },
    showVariants(filtered, coPtms = [], sites = []) {
      this._selectedResi = null; this._lastRender = () => this.showVariants(filtered, coPtms, sites);
      const severity = ['Likely pathogenic or pathogenic', 'Predicted deleterious', 'Uncertain significance', 'Likely benign or benign'];
      const posMap = new Map();
      (filtered || []).forEach(v => {
        const ex = posMap.get(v.position);
        if (!ex) posMap.set(v.position, { position: v.position, color: v.consequenceColor, topConsequence: v.consequence, variants: [v] });
        else { ex.variants.push(v); const ei = severity.indexOf(ex.topConsequence), ni = severity.indexOf(v.consequence); if (ni >= 0 && (ei < 0 || ni < ei)) { ex.topConsequence = v.consequence; ex.color = v.consequenceColor; } }
      });
      const hoverMap = new Map(); const active = new Map(); let posCount = 0, varCount = 0;
      posMap.forEach((d, pos) => { if (!this.allChainsAddStyle(pos, { sphere: { radius: 1.8, color: d.color, opacity: 0.92 } }, { atom: 'CA' })) return; hoverMap.set(pos, d); active.set(pos, d.color); posCount++; varCount += d.variants.length; });
      let ptmCount = 0;
      (coPtms || []).forEach(ptm => {
        if (posMap.has(ptm.position)) return;
        if (!this.allChainsAddStyle(ptm.position, { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }, { atom: 'CA' })) return;
        if (!hoverMap.has(ptm.position)) hoverMap.set(ptm.position, { position: ptm.position, color: ptm.color, category: ptm.category });
        active.set(ptm.position, ptm.color); ptmCount++;
        if (ptm.endPosition && ptm.endPosition !== ptm.position && !posMap.has(ptm.endPosition)) { if (this.allChainsAddStyle(ptm.endPosition, { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }, { atom: 'CA' })) active.set(ptm.endPosition, ptm.color); }
      });
      this._drawSiteSpheres(sites, hoverMap, active);
      this._activeSpheres = active; this._bindHover(hoverMap, 'variant'); this.viewer.render();
      return { posCount, varCount, ptmCount };
    },
    refreshPTMDisplay(ptms, ptmGroups, sites = [], extras = []) {
      if (!this.viewer) return false;   // markers are their own channel — updatable even while focused
      this.viewer.setStyle({}, {});
      const base = { opacity: 0.82, thickness: 0.25, ribbonWidth: 0.6 };
      this._applyModeStyles(this.activeColoringMode, this._lastColoringContext || {}, base);
      this._drawLigands();
      return this.showPTMs(ptms, ptmGroups, sites, extras);
    },
    _bindHover(map, mode) {
      const self = this; this._hoverMap = map;
      const resolveUni = atom => self._reverseResidueMapForChain(atom.chain).get(atom.resi) || atom.resi;
      // A protein residue on a chain that isn't ours (a partner subunit in a complex) — we have no
      // annotations for it; never resolve it through our mapping (would mirror wrong data / jump camera).
      const isPartner = atom => { const our = self.currentStructure?.chainIds; return Array.isArray(our) && our.length && atom.chain != null && !our.includes(atom.chain); };
      this.viewer.setHoverable({}, true,
        (atom) => {
          if (!atom || !atom.resi) return;
          if (atom.hetflag && atom.resn !== 'HOH' && atom.resn !== 'WAT') { self._hoverLigand(atom); if (self.hoverCb) self.hoverCb(null, mode, null); return; }
          if (!atom.hetflag && isPartner(atom)) { if (self.hoverCb) self.hoverCb(null, mode, null); return; }
          const d = map.get(resolveUni(atom)); if (d && self.hoverCb) self.hoverCb({ ...d, pdbResi: atom.resi, chain: atom.chain }, mode, {}, atom.chain);
        },
        () => { self._clearLigandHover(); if (self.hoverCb) self.hoverCb(null, mode, null); });
      this.viewer.setClickable({}, true,
        (atom) => {
          if (!atom || !atom.resi) return;
          if (atom.hetflag) {
            // HETATMs never fall through to the protein-residue click. Water is inert (clicking it did
            // nothing useful and wrongly opened a residue panel); other het → the ligand panel.
            if (atom.resn === 'HOH' || atom.resn === 'WAT') return;
            if (self.ligandClickCb) self.ligandClickCb({ resn: atom.resn, resi: atom.resi, chain: atom.chain ?? null });
            return;
          }
          // Partner-subunit residue (different protein): route to the partner handler (multichain) with the
          // PDB author residue — NOT our clickCb, which would mis-map it.
          if (isPartner(atom)) { if (self.partnerClickCb) self.partnerClickCb({ chain: atom.chain, resi: atom.resi, resn: atom.resn }); return; }
          const uni = resolveUni(atom); const d = map.get(uni); if (self.clickCb) self.clickCb(d || { position: uni, pdbResi: atom.resi }, mode, atom.chain);
        });
    },
    _hoverLigand() {},        // Mol* handles native ligand hover highlight
    _clearLigandHover() {},

    // ---- focus: zoom + Mol* native focus rep. Returns nearby UniProt positions (≤5 Å) for the panel,
    // and sets _focusRegion (chain|pdbResi shown as sticks) so those residues' spheres are hidden. ----
    // region holds BOTH 'chain|resi' AND chain-agnostic '|resi' so it matches markers whether they were
    // stored with a chain (PDB) or without (AlphaFold, where allChainsAddStyle omits the chain).
    // Stick colours for the focus view. Mol*'s focus rep shows our nearby loci PLUS its own ~5 Å
    // surrounding residues, so colouring only the nearby set leaves Mol*-added pathogenic residues
    // grey. Instead convert EVERY annotated residue (Map<uniProtPos,{color}>) to PDB coords across
    // all chains; overpaint only affects atoms actually rendered, so whatever sticks Mol* shows, any
    // annotated one is coloured. The selected residue is excluded frame-side.
    _stickAnnotations(_nearPdb, annotatedResidues) {
      const out = [];
      const s = this.currentStructure;
      if (annotatedResidues && annotatedResidues.size) {
        const chains = s?.chainIds?.length ? s.chainIds : [s?.chainId ?? null];
        for (const ch of chains) {
          const rev = this._reverseResidueMapForChain(ch); // pdbResi -> uniProtPos
          for (const [pdbResi, uni] of rev) {
            const meta = annotatedResidues.get(uni);
            if (meta && meta.color) out.push({ chain: ch, resi: pdbResi, color: meta.color });
          }
        }
      }
      // Multichain: partner-chain neighbours in the focus set get OUR disease colour for that partner
      // residue (chain|resi -> colour), so focusing one of our interface residues colours the partner
      // contacts the same way our own neighbours are coloured.
      if (this._partnerColorMap && this._partnerColorMap.size && Array.isArray(_nearPdb)) {
        for (const nb of _nearPdb) {
          const color = this._partnerColorMap.get((nb.chain == null ? '' : nb.chain) + '|' + nb.resi);
          if (color) out.push({ chain: nb.chain, resi: nb.resi, color });
        }
      }
      return out;
    },
    // Basic hover labels (resn + resi + chain + UniProt pos) for EVERY focus-stick residue, so hovering an
    // unannotated side chain still identifies it. Merged in _flush UNDER the richer annotation labels.
    _buildFocusLabels(nearPdb) {
      const out = {};
      if (Array.isArray(nearPdb) && nearPdb.length) {
        for (const nb of nearPdb) {
          const k = (nb.chain == null ? '' : nb.chain) + '|' + nb.resi;
          // The tooltip already shows "resn resi" as its header (frame side) — only add the NEW info here
          // (chain + UniProt position) so the hover isn't redundant.
          const uni = this._reverseResidueMapForChain(nb.chain).get(nb.resi);
          const parts = [];
          if (nb.chain != null) parts.push(`chain ${nb.chain}`);
          if (uni != null) parts.push(`UniProt ${uni}`);
          out[k] = parts.length ? parts.join(' · ') : `residue ${nb.resi}`;
        }
      }
      this._focusLabels = out;
    },
    // Ligand residue keys ('chain|resi') with any atom within `radius` Å of the focus targets — used
    // to keep only pocket ligands visible while focused (distant AlphaFill ligands are hidden).
    // Camera radius (Å) that frames the WHOLE pocket the same way for a residue OR a ligand: centred on the
    // focused entity, reaching out to its farthest nearby residue. Without this the zoom was entity-size +
    // shell, so a tiny residue zoomed in hard while a large ligand stayed far out — the "ligands act
    // differently" feel. Computed from the same atoms/nearby set both paths already use.
    _pocketCamRadius(targets, nearPdb) {
      const nd = this.nearbyDistance || 5;
      if (!targets || !targets.length) return nd + 3;
      let cx = 0, cy = 0, cz = 0;
      for (const t of targets) { cx += t.x; cy += t.y; cz += t.z; }
      cx /= targets.length; cy /= targets.length; cz /= targets.length;
      const keys = new Set((nearPdb || []).map(n => (n.chain == null ? '' : n.chain) + '|' + n.resi));
      let max2 = 0;
      const consider = (x, y, z) => { const dx = x - cx, dy = y - cy, dz = z - cz, d2 = dx * dx + dy * dy + dz * dz; if (d2 > max2) max2 = d2; };
      for (const t of targets) consider(t.x, t.y, t.z);                       // the entity's own extent
      for (const a of atoms) {                                                // + every nearby residue atom
        if (a.hetflag) continue;
        if (keys.has((a.chain == null ? '' : a.chain) + '|' + a.resi)) consider(a.x, a.y, a.z);
      }
      return Math.sqrt(max2) + 3; // small constant padding so the outermost sticks aren't flush to the edge
    },
    _nearbyLigandKeys(targets, radius) {
      const r2 = radius * radius; const keep = new Set();
      for (const a of atoms) {
        if (!a.hetflag || a.resn === 'HOH' || a.resn === 'WAT') continue;
        const k = (a.chain == null ? '' : a.chain) + '|' + a.resi;
        if (keep.has(k)) continue;
        for (const t of targets) {
          const dx = a.x - t.x, dy = a.y - t.y, dz = a.z - t.z;
          if (dx * dx + dy * dy + dz * dz <= r2) { keep.add(k); break; }
        }
      }
      return keep;
    },
    _focusNeighbourhood(targets, selChain, pdbResi) {
      const nd = (liveRef || api).nearbyDistance || 5;
      const near = new Set(); const region = new Set(); const d2 = nd * nd;
      const nearPdb = new Map(); // "chain|resi" -> { chain, resi } — PDB coords for sticks loci
      // AlphaFill drags donor binding-site PEPTIDE residues (real amino acids, group_PDB=ATOM) into the
      // transplant chains alongside the ligand. They'd be collected as "protein" here and shown as grey
      // pocket sticks. Identify transplant chains = chains carrying a HETATM ligand that are NOT the main
      // (largest) protein chain, and skip them. (Experimental structures keep protein+ligand in one big
      // chain, which stays the main chain, so nothing real is excluded there.)
      const resByChain = new Map(); const hetChains = new Set();
      for (const a of atoms) {
        const c = a.chain == null ? '' : a.chain;
        if (!resByChain.has(c)) resByChain.set(c, new Set());
        resByChain.get(c).add(a.resi);
        if (a.hetflag && a.resn !== 'HOH' && a.resn !== 'WAT') hetChains.add(c);
      }
      let mainChain = null, mainCount = -1;
      for (const [c, s] of resByChain) if (s.size > mainCount) { mainCount = s.size; mainChain = c; }
      // The donor-peptide transplant artefact is an AlphaFill (computed-model) thing. In a multi-chain
      // EXPERIMENTAL complex (e.g. a glycosylated receptor) EVERY subunit carries HETATM glycans, so the old
      // "non-largest HETATM chain = transplant" rule wrongly excluded real partner chains — and even the
      // FOCUSED chain when it wasn't the largest, emptying Nearby and showing no neighbour sticks. So only
      // apply it to AlphaFold-based models; for PDB/experimental there are no transplant peptides to skip.
      const _st = this.currentStructure;
      const _isModel = !!_st && (_st.source === 'AlphaFold' || /alphafill/i.test(_st.provider || _st.source || _st.label || ''));
      const transplantChains = _isModel ? new Set([...hetChains].filter(c => c !== mainChain)) : new Set();
      // Reverse pdb→UniProt maps are cached per chain — rebuilding one per neighbour atom was O(N·ranges)
      // and a big chunk of the focus lag (and risked an empty "Nearby" list when it timed out mid-loop).
      const revCache = new Map();
      const toUni = (ch, ri) => { const k = ch == null ? '' : ch; if (!revCache.has(k)) revCache.set(k, this._reverseResidueMapForChain(ch)); return revCache.get(k).get(ri); };
      const add = (ch, ri) => {
        region.add((ch == null ? '' : ch) + '|' + ri); region.add('|' + ri);
        nearPdb.set((ch == null ? '' : ch) + '|' + ri, { chain: ch == null ? null : ch, resi: ri });
      };
      // Multichain: neighbours on a PARTNER chain (a different protein) don't map to our UniProt — collect
      // them separately so the "Nearby" list can still show them (with the partner's annotations).
      const ourChains = Array.isArray(this.currentStructure?.chainIds) ? new Set(this.currentStructure.chainIds) : null;
      const partnerNear = new Map(); // "chain|resi" -> { chain, resi }
      add(selChain, pdbResi);
      for (const a of atoms) {
        if (a.hetflag) continue;
        if (transplantChains.has(a.chain == null ? '' : a.chain)) continue; // skip AlphaFill donor peptides
        for (const t of targets) {
          const dx = a.x - t.x, dy = a.y - t.y, dz = a.z - t.z;
          if (dx * dx + dy * dy + dz * dz <= d2) {
            add(a.chain, a.resi);
            const uni = toUni(a.chain, a.resi);
            if (uni != null) near.add(uni);
            else if (ourChains && ourChains.size && a.chain != null && !ourChains.has(a.chain)) partnerNear.set((a.chain == null ? '' : a.chain) + '|' + a.resi, { chain: a.chain, resi: a.resi });
            break;
          }
        }
      }
      // Mol*'s focus representation ALSO renders residues within ~5 Å of the stick set as sticks
      // (its expandRadius default). Grow the sphere-exclusion region (only) to cover that shell so a
      // residue never shows BOTH a stick and a sphere. near (the list) + nearPdb (sticks loci sent to
      // Mol*) stay at the user's nearby distance.
      const SURROUND2 = 5 * 5;
      const stickAtoms = atoms.filter(a => !a.hetflag && nearPdb.has((a.chain == null ? '' : a.chain) + '|' + a.resi));
      for (const a of atoms) {
        if (a.hetflag) continue;
        const key = (a.chain == null ? '' : a.chain) + '|' + a.resi;
        if (region.has(key)) continue;
        for (const t of stickAtoms) {
          const dx = a.x - t.x, dy = a.y - t.y, dz = a.z - t.z;
          if (dx * dx + dy * dy + dz * dz <= SURROUND2) { region.add(key); region.add('|' + a.resi); break; }
        }
      }
      return { near, region, nearPdb: Array.from(nearPdb.values()), partnerNear: Array.from(partnerNear.values()) };
    },
    // _annotations and opts are passed by modal.js but were previously ignored.
    // opts.rezoom === false means "update sticks/exclusion region but don't move camera" —
    // used by applyMode() when re-applying focus after a filter/colour toggle.
    focusResidue(resi, chain = null, _annotations = null, opts = {}) {
      this.clearPocket(); // a residue focus supersedes a pocket-surface highlight
      const wasFocused = !!this._inFocusMode; // moving residue→residue (already zoomed in) ⇒ pan, don't re-zoom
      this._selectedResi = resi; this._inFocusMode = true;
      const sel = this.residueSelectorForChain(resi, chain);
      const pdbResi = sel.resi, selChain = sel.chain != null ? sel.chain : (chain != null ? chain : null);
      if (pdbResi === -999999) return new Set([resi]);
      const targets = atoms.filter(a => a.resi === pdbResi && (selChain == null || a.chain === selChain));
      const { near, region, nearPdb, partnerNear } = this._focusNeighbourhood(targets, selChain, pdbResi);
      this._nearbyPartners = partnerNear || []; // multichain: partner-chain neighbours for the Nearby list
      // Ligands while focused on a RESIDUE: in a normal structure (few ligands) keep the ones within the
      // pocket so you can see what the residue contacts. But an AlphaFill model packs hundreds of
      // overlapping transplant copies — there, "nearby" is a blob (and long lipids trail off-screen), so
      // hide ALL ligands for a clean residue pocket. Threshold: >20 ligands ⇒ transplant-crowded.
      const transplantCrowded = this.enumerateLigands().length > 20;
      this._focusNearbyLigands = transplantCrowded
        ? new Set()                                                          // hide every ligand
        : this._nearbyLigandKeys(targets, (this.nearbyDistance || 5) + 3);   // keep pocket ligands
      near.add(resi); _focusRegion = region; _focusHighlightResi = { chain: selChain, resi: pdbResi }; this._lastCamTarget = { chain: selChain, resi: pdbResi };
      if (opts.showOtherSpheres !== undefined) _showOtherSpheres = opts.showOtherSpheres;
      const rezoom = opts.rezoom !== false;
      // Colour focus sticks by each nearby residue's annotation (variant pathogenicity / PTM / site),
      // independent of which spheres are shown. annotatedResidues is Map<uniProtPos, {color}>.
      const annotations = this._stickAnnotations(nearPdb, _annotations && _annotations.annotatedResidues);
      this._buildFocusLabels(nearPdb); // basic hover labels for every focus side chain
      if (rezoom) send('focus',   { chain: selChain, resi: pdbResi, neighbors: nearPdb, annotations }).catch(() => {});
      else        send('refocus', { chain: selChain, resi: pdbResi, neighbors: nearPdb, annotations }).catch(() => {});
      // Send the ligand-hide AFTER the focus command so the camera zoom (in doFocus) is queued first and
      // isn't gated behind this potentially-expensive hide on AlphaFill models (hundreds of ligands).
      this._drawLigands();
      // Record the focus state for scene export (PyMOL/VMD): selected residue + neighbour sticks with
      // their annotation colours, in PDB numbering. Sticks without a colour export as element/Jmol.
      const selKey = (selChain == null ? '' : selChain) + '|' + pdbResi;
      const annoByKey = new Map(annotations.map(a => [(a.chain == null ? '' : a.chain) + '|' + a.resi, a.color]));
      this._focusState = {
        selChain: selChain == null ? '' : selChain, pdbResi,
        sticks: nearPdb
          .filter(nb => ((nb.chain == null ? '' : nb.chain) + '|' + nb.resi) !== selKey)
          .map(nb => { const k = (nb.chain == null ? '' : nb.chain) + '|' + nb.resi; return { chain: nb.chain == null ? '' : nb.chain, resi: nb.resi, het: false, color: annoByKey.get(k) }; }),
      };
      // skipCamRestore=true when zooming: the camera is animating toward the residue — don't
      // let setMarkers snapshot+restore the pre-animation position and cancel the zoom.
      markersDirty = true; this._flush(rezoom); this._applyMarkerVisibility();
      // Camera LAST — after the cartoon/marker/ligand commits above — so the zoom animation isn't stuttered
      // by concurrent state updates. radius frames the whole pocket (same logic as the ligand path).
      // panOnly: when already zoomed in, moving to another residue PANS at the current zoom (no re-zoom),
      // so the camera doesn't do an extra zoom-out/zoom-in between residues.
      this._lastCamRadius = this._pocketCamRadius(targets, nearPdb);
      if (rezoom) send('camFocus', { chain: selChain, resi: pdbResi, radius: this._lastCamRadius, panOnly: wasFocused }).catch(() => {});
      return near;
    },
    focusLigand(resn, resi, chain = null, opts = {}) {
      this.clearPocket();
      this._inFocusMode = true;
      this._focusState = null; // residue-focus export state doesn't apply to a ligand focus
      const targets = atoms.filter(a => a.resn === resn && a.resi === resi && (chain == null || a.chain === chain));
      const { near, region, nearPdb } = this._focusNeighbourhood(targets, chain, resi);
      _focusRegion = region; _focusHighlightResi = null; this._lastCamTarget = { chain, resi }; // ligand focus has no single "clicked residue"
      if (opts.showOtherSpheres !== undefined) _showOtherSpheres = opts.showOtherSpheres;
      const rezoom = opts.rezoom !== false;
      // Show only the focused ligand (AlphaFill models pack many close together).
      this._focusNearbyLigands = new Set([(chain == null ? '' : chain) + '|' + resi]);
      // focusLigand receives annotatedResidues via opts (4th arg), not a separate _annotations param.
      const annotations = this._stickAnnotations(nearPdb, opts.annotatedResidues);
      this._buildFocusLabels(nearPdb); // basic hover labels for every focus side chain
      if (rezoom) send('focus',   { chain, resi, neighbors: nearPdb, annotations }).catch(() => {});
      else        send('refocus', { chain, resi, neighbors: nearPdb, annotations }).catch(() => {});
      this._drawLigands(); // after the focus command, so the camera zoom isn't gated behind the ligand-hide
      markersDirty = true; this._flush(rezoom); this._applyMarkerVisibility();
      // radius frames the whole pocket centred on the ligand — same computation as the residue path, so the
      // zoom <-> nearby relationship is identical whether a residue or a ligand is selected.
      this._lastCamRadius = this._pocketCamRadius(targets, nearPdb);
      if (rezoom) send('camFocus', { chain, resi, radius: this._lastCamRadius }).catch(() => {}); // camera LAST → uncontended zoom animation
      return near;
    },
    enumerateLigands() { return filterAtoms({ hetflag: true }).reduce((acc, a) => { const k = (a.chain ?? '') + '|' + a.resi; if (a.resn !== 'HOH' && a.resn !== 'WAT' && !acc._s.has(k)) { acc._s.add(k); acc.push({ resn: a.resn, resi: a.resi, chain: a.chain ?? null }); } return acc; }, Object.assign([], { _s: new Set() })); },
    ligandContactsByResidue() { return new Map(); },
    // Serialise the on-screen scene for PyMOL/VMD export, rebuilt from the capture buffers (already in
    // PDB author numbering). Reproduces the annotated overview: per-residue cartoon overrides, every
    // visible annotation Cα sphere, ligand groups, and — when zoomed in — the focus stick set.
    getSceneState() {
      if (!this.currentPdbText) return null;
      const base = cartoonBase || '#d0d0d0';
      const baseLc = base.toLowerCase();
      // Cartoon: only residues whose colour differs from the base (keeps the script compact).
      const cartoon = [];
      for (const o of cartoonOverrides) {
        if (o.color && o.color.toLowerCase() !== baseLc) cartoon.push({ chain: o.chain || '', resi: o.resi, color: o.color });
      }
      // Annotation Cα spheres: current sphere buffer minus any residue shown as a focus stick.
      const spheres = [];
      for (const sp of sphereBuffer) {
        const key = (sp.chain == null ? '' : sp.chain) + '|' + sp.resi;
        if (_focusRegion && (_focusRegion.has(key) || _focusRegion.has('|' + sp.resi))) continue;
        spheres.push({ chain: sp.chain || '', resi: sp.resi, color: sp.color });
      }
      // Ligands / cofactors (non-water HETATM), unique per chain+resi+resn.
      const ligands = []; const seen = new Set();
      for (const a of atoms) {
        if (!a.hetflag || a.resn === 'HOH' || a.resn === 'WAT') continue;
        const key = (a.chain || '') + '|' + a.resi + '|' + a.resn;
        if (seen.has(key)) continue; seen.add(key);
        ligands.push({ chain: a.chain || '', resi: a.resi, resn: a.resn, ion: this.ION_CODES.has(a.resn) });
      }
      const focus = (this._inFocusMode && this._focusState) ? this._focusState : null;
      return {
        format: this.currentFormat === 'mmcif' ? 'cif' : 'pdb',
        coordinates: this.currentPdbText,
        coloringMode: this.activeColoringMode || 'default',
        cartoonBase: base,
        cartoonOpacity: this._inFocusMode ? 0.42 : 0.82,
        sphereRadius: 1.8, sphereOpacity: 0.92,
        camera: null,
        cartoon, spheres, ligands, focus,
      };
    },
  };

  return api;
})();

// Mol* is the only viewer engine; the former 3Dmol path (viewer.js, lib/3Dmol-min.js) has been
// removed. viewer-molstar.js defines StructureViewer directly so modal.js picks it up.
try {
  if (typeof StructureViewer !== 'undefined' && StructureViewer) { Object.assign(StructureViewer, MolstarViewer); }
  else { window.StructureViewer = MolstarViewer; }
  console.info('[UFV] viewer engine: Mol* (sandboxed iframe)');
} catch (e) { console.warn('[UFV] Mol* engine setup failed:', e); }
