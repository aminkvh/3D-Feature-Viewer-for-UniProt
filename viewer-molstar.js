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
  let cartoonDirty = false, markersDirty = false;
  let lastCartoonJSON = '', lastMarkersJSON = '', lastLabelsJSON = '';  // dedup: skip unchanged sends
  let _clickFn = null, _hoverFn = null, _hoverOutFn = null;  // native-pick handlers (from _bindHover)
  let _focusRegion = null;        // Set of 'chain|pdbResi' shown as sticks in focus → their spheres are excluded
  let _showOtherSpheres = true;   // when false, hide all annotation spheres outside the focus region
  // The live object: after Object.assign(StructureViewer, MolstarViewer) the modal writes state onto
  // StructureViewer, not onto `api`. Module-level callbacks (render→_flush, atoms→cache) must read
  // that same object, so init() captures it here.
  let liveRef = null;

  function ensureReady() { if (!readyPromise) readyPromise = new Promise(res => { readyResolve = res; }); return readyPromise; }

  function onMessage(ev) {
    if (!iframe || ev.source !== iframe.contentWindow) return;
    const m = ev.data;
    if (!m || m.ns !== NS) return;
    if (m.evt === 'ready') { if (readyResolve) readyResolve(); return; }
    if (m.evt === 'atoms') { atoms = m.atoms || []; (liveRef || api)._buildObservedResiCache(); return; }
    if (m.evt === 'pick') {
      if (m.kind === 'click') { if (_clickFn && m.atom) _clickFn(m.atom); }
      else if (m.kind === 'hover') { if (m.atom) { if (_hoverFn) _hoverFn(m.atom, null, {}); } else if (_hoverOutFn) _hoverOutFn(); }
      else if (m.kind === 'dblclick') { const cb = (liveRef || api).dblClickCb; if (cb) cb(); }
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
    ION_CODES: new Set(['NA','K','LI','RB','CS','MG','CA','SR','BA','ZN','FE','FE2','MN','MN3','CU','CU1','CO','NI','CD','HG','CL','BR','IOD','FLO','F','AL','PT','AU','AG','PB','SO4','PO4']),
    clickCb: null, hoverCb: null, dblClickCb: null, ligandClickCb: null,

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
        setBackgroundColor: (color) => send('background', { color }).catch(() => {}),
      };
    },

    async loadStructure(url, structure = null) {
      await ensureReady();
      this.currentStructure = structure;
      this._observedResi = this._observedResiByChain = null;
      atoms = []; _focusRegion = null; lastCartoonJSON = lastMarkersJSON = lastLabelsJSON = '';
      const f = inferFormat(url, structure);
      const format = typeof f === 'string' ? f : f.format;
      const isBinary = typeof f === 'string' ? false : f.isBinary;
      const resp = await fetch(url);
      if (!resp.ok) throw new Error('structure fetch ' + resp.status);
      let data;
      if (isBinary) { data = new Uint8Array(await resp.arrayBuffer()); this.currentPdbText = null; }
      else { data = await resp.text(); this.currentPdbText = data; }
      this.currentFormat = format === 'pdb' ? 'pdb' : 'mmcif';
      await send('loadData', { data, format, isBinary });   // frame replies after load; 'atoms' event follows
      return true;
    },

    async clearModel() {
      this.currentStructure = this.currentPdbText = this.currentFormat = null;
      this._observedResi = this._observedResiByChain = null; atoms = [];
      this._selectedResi = null; this._inFocusMode = false; this._focusState = null;
      if (iframe) { try { await send('clear'); } catch (_) {} }
    },

    resize() { send('resize').catch(() => {}); },
    resetView() { this._selectedResi = null; this._inFocusMode = false; this._focusState = null; _focusRegion = null; _showOtherSpheres = true; send('unfocus').catch(() => {}); markersDirty = true; this._flush(true); },
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
        out.push({ uniPos: uni != null ? uni : null, chain: a.chain ?? null, resi: a.resi, ca: { x: a.x, y: a.y, z: a.z } });
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
        const s = this.currentStructure;
        const isAF = !s || (s.source === 'AlphaFold' && !s.isoform) || !s.mappedRanges?.length;
        const ourChains = s?.chainIds?.length ? new Set(s.chainIds) : (s?.chainId ? new Set([s.chainId]) : null);
        const reverseCache = new Map();
        const toUni = (chain, resi) => { if (!reverseCache.has(chain)) reverseCache.set(chain, this._reverseResidueMapForChain(chain)); return reverseCache.get(chain).get(resi); };
        this.viewer.setStyle({}, { cartoon: { ...base, colorfunc: atom => { if (!isAF && ourChains && atom.chain != null && !ourChains.has(atom.chain)) return '#b9c2cf'; const uni = isAF ? atom.resi : toUni(atom.chain, atom.resi); return posColor.get(uni) || '#b9c2cf'; } } });
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
        const byColor = new Map();
        for (const o of cartoonOverrides) { if (!byColor.has(o.color)) byColor.set(o.color, []); byColor.get(o.color).push([o.chain, o.resi]); }
        const payload = { base: cartoonBase || '#d0d0d0', groups: [...byColor.entries()].map(([color, residues]) => ({ color, residues })) };
        const j = JSON.stringify(payload);
        if (j !== lastCartoonJSON) { lastCartoonJSON = j; send('setCartoon', payload).catch(() => {}); }
        cartoonDirty = false;
      }
      if (markersDirty) {
        const byColor = new Map();
        for (const s of sphereBuffer) {
          const key = (s.chain == null ? '' : s.chain) + '|' + s.resi;
          if (_focusRegion && _focusRegion.has(key)) continue;          // shown as sticks in focus
          if (_focusRegion && !_showOtherSpheres) continue;             // hide all non-focus when checkbox off
          if (!byColor.has(s.color)) byColor.set(s.color, { radius: s.radius, residues: [] }); byColor.get(s.color).residues.push([s.chain, s.resi]);
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
          `<strong>${esc(v.wildType)}${v.position}${esc(v.mutant)}</strong> <span style="color:${esc(v.consequenceColor || '#9e9e9e')}">${esc(v.consequence)}</span>`
        ).join(' &nbsp;|&nbsp; ');
      }
      if (d.isSite) return esc(d.description || 'Site');
      if (d.category) return esc(d.description ? `${d.category}: ${d.description}` : d.category);
      return esc(d.description || d.label || '') || null;
    },
    _drawLigands() {},                                      // M2b ligand layer — wired in M2c

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
        const g = ptmGroups[ptm.category]; if (!g || !g.visible || ptm.visible === false) return;
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
      this.viewer.setHoverable({}, true,
        (atom) => {
          if (!atom || !atom.resi) return;
          if (atom.hetflag && atom.resn !== 'HOH' && atom.resn !== 'WAT') { self._hoverLigand(atom); if (self.hoverCb) self.hoverCb(null, mode, null); return; }
          const d = map.get(resolveUni(atom)); if (d && self.hoverCb) self.hoverCb(d, mode, {}, atom.chain);
        },
        () => { self._clearLigandHover(); if (self.hoverCb) self.hoverCb(null, mode, null); });
      this.viewer.setClickable({}, true,
        (atom) => {
          if (!atom || !atom.resi) return;
          if (atom.hetflag && atom.resn !== 'HOH' && self.ligandClickCb) { self.ligandClickCb({ resn: atom.resn, resi: atom.resi, chain: atom.chain ?? null }); return; }
          const d = map.get(resolveUni(atom)); if (self.clickCb) self.clickCb(d || { position: resolveUni(atom) }, mode, atom.chain);
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
      if (!annotatedResidues || !annotatedResidues.size) return [];
      const s = this.currentStructure;
      const chains = s?.chainIds?.length ? s.chainIds : [s?.chainId ?? null];
      const out = [];
      for (const ch of chains) {
        const rev = this._reverseResidueMapForChain(ch); // pdbResi -> uniProtPos
        for (const [pdbResi, uni] of rev) {
          const meta = annotatedResidues.get(uni);
          if (meta && meta.color) out.push({ chain: ch, resi: pdbResi, color: meta.color });
        }
      }
      return out;
    },
    _focusNeighbourhood(targets, selChain, pdbResi) {
      const nd = (liveRef || api).nearbyDistance || 5;
      const near = new Set(); const region = new Set(); const d2 = nd * nd;
      const nearPdb = new Map(); // "chain|resi" -> { chain, resi } — PDB coords for sticks loci
      const add = (ch, ri) => {
        region.add((ch == null ? '' : ch) + '|' + ri); region.add('|' + ri);
        nearPdb.set((ch == null ? '' : ch) + '|' + ri, { chain: ch == null ? null : ch, resi: ri });
      };
      add(selChain, pdbResi);
      for (const a of atoms) {
        if (a.hetflag) continue;
        for (const t of targets) {
          const dx = a.x - t.x, dy = a.y - t.y, dz = a.z - t.z;
          if (dx * dx + dy * dy + dz * dz <= d2) { add(a.chain, a.resi); const uni = this._reverseResidueMapForChain(a.chain).get(a.resi); if (uni != null) near.add(uni); break; }
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
      return { near, region, nearPdb: Array.from(nearPdb.values()) };
    },
    // _annotations and opts are passed by modal.js but were previously ignored.
    // opts.rezoom === false means "update sticks/exclusion region but don't move camera" —
    // used by applyMode() when re-applying focus after a filter/colour toggle.
    focusResidue(resi, chain = null, _annotations = null, opts = {}) {
      this._selectedResi = resi; this._inFocusMode = true;
      const sel = this.residueSelectorForChain(resi, chain);
      const pdbResi = sel.resi, selChain = sel.chain != null ? sel.chain : (chain != null ? chain : null);
      if (pdbResi === -999999) return new Set([resi]);
      const targets = atoms.filter(a => a.resi === pdbResi && (selChain == null || a.chain === selChain));
      const { near, region, nearPdb } = this._focusNeighbourhood(targets, selChain, pdbResi);
      near.add(resi); _focusRegion = region;
      if (opts.showOtherSpheres !== undefined) _showOtherSpheres = opts.showOtherSpheres;
      const rezoom = opts.rezoom !== false;
      // Colour focus sticks by each nearby residue's annotation (variant pathogenicity / PTM / site),
      // independent of which spheres are shown. annotatedResidues is Map<uniProtPos, {color}>.
      const annotations = this._stickAnnotations(nearPdb, _annotations && _annotations.annotatedResidues);
      if (rezoom) send('focus',   { chain: selChain, resi: pdbResi, neighbors: nearPdb, annotations }).catch(() => {});
      else        send('refocus', { chain: selChain, resi: pdbResi, neighbors: nearPdb, annotations }).catch(() => {});
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
      markersDirty = true; this._flush(rezoom);
      return near;
    },
    focusLigand(resn, resi, chain = null, opts = {}) {
      this._inFocusMode = true;
      this._focusState = null; // residue-focus export state doesn't apply to a ligand focus
      const targets = atoms.filter(a => a.resn === resn && a.resi === resi && (chain == null || a.chain === chain));
      const { near, region, nearPdb } = this._focusNeighbourhood(targets, chain, resi);
      _focusRegion = region;
      if (opts.showOtherSpheres !== undefined) _showOtherSpheres = opts.showOtherSpheres;
      const rezoom = opts.rezoom !== false;
      // focusLigand receives annotatedResidues via opts (4th arg), not a separate _annotations param.
      const annotations = this._stickAnnotations(nearPdb, opts.annotatedResidues);
      if (rezoom) send('focus',   { chain, resi, neighbors: nearPdb, annotations }).catch(() => {});
      else        send('refocus', { chain, resi, neighbors: nearPdb, annotations }).catch(() => {});
      markersDirty = true; this._flush(rezoom);
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
