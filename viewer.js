/* ============================================
   3Dmol.js Viewer Wrapper — Chrome Extension
   ============================================ */

const StructureViewer = {
    viewer: null,
    hoverCb: null,
    clickCb: null,
    ligandClickCb: null,
    dblClickCb: null,
    _selectedResi: null,
    currentPdbText: '',
    currentFormat: 'pdb',
    currentStructure: null,
    activeColoringMode: 'default',
    _lastColoringContext: {},
    _observedResi: null,
    _observedResiByChain: null,
    _inFocusMode: false,
    _resizeBound: false,
    _lastRender: null, // closure that re-draws the active PTM/variant spheres (for WebGL restore)
    _activeSpheres: null, // Map<uniProtPos, color> of the currently-drawn annotation spheres
                          // (PTM and/or variant) — used to keep them visible while zoomed in

    init(container) {
        const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
        const bgColor = prefersDark ? '#0c111b' : '#f0f2f5';
        this.viewer = $3Dmol.createViewer(container, {
            backgroundColor: bgColor,
            antialias: true,
            cartoonQuality: 8,
        });
        // Fire hover near-instantly (3Dmol default is 500 ms) so ligand highlight feels immediate.
        try { this.viewer.setHoverDuration(5); } catch (_) {}
        if (!this._resizeBound) {
            window.addEventListener('resize', () => { if (this.viewer) this.viewer.resize(); });
            this._resizeBound = true;
        }

        // Double-click on background = reset view
        container.addEventListener('dblclick', () => {
            if (this.dblClickCb) this.dblClickCb();
        });

        // Recover from WebGL context loss (browsers may discard contexts in background tabs).
        container.addEventListener('webglcontextlost', e => { e.preventDefault(); }, false);
        container.addEventListener('webglcontextrestored', () => {
            if (this.currentPdbText && this.currentStructure) {
                // Clear any half-restored state first, then re-load the already-fetched PDB
                // text to rebuild the WebGL model from scratch.
                this.viewer.removeAllModels();
                this.viewer.removeAllLabels();
                this.viewer.removeAllShapes();
                this.viewer.addModel(this.currentPdbText, this.currentFormat);
                this._buildObservedResiCache();
                this.applyCartoonColoring(this.activeColoringMode, this._lastColoringContext);
                // Re-draw the active PTM/variant spheres and re-bind hover/click so the
                // recovered context is fully interactive again (not just a static cartoon).
                if (this._lastRender) this._lastRender();
                this.viewer.zoomTo();
                this.viewer.render();
            }
        }, false);

        // Invert scroll-wheel zoom direction so scroll-out = zoom in
        let _reWheel = false;
        container.addEventListener('wheel', e => {
            if (_reWheel) return;
            e.preventDefault();
            e.stopImmediatePropagation();
            _reWheel = true;
            e.target.dispatchEvent(new WheelEvent('wheel', {
                bubbles: true, cancelable: true,
                deltaX: e.deltaX, deltaY: e.deltaY, deltaZ: e.deltaZ,
                deltaMode: e.deltaMode, view: e.view,
                clientX: e.clientX, clientY: e.clientY,
                ctrlKey: e.ctrlKey, shiftKey: e.shiftKey,
                altKey: e.altKey, metaKey: e.metaKey,
            }));
            _reWheel = false;
        }, { capture: true, passive: false });
    },

    async loadStructure(pdbUrl, structure = null) {
        // Validate URL — only allow HTTPS from trusted structure hosts (experimental + the
        // computed-model providers surfaced via 3D-Beacons).
        const ALLOWED_HOSTS = ['alphafold.ebi.ac.uk', 'www.ebi.ac.uk', 'files.rcsb.org',
            'swissmodel.expasy.org', 'www.modelarchive.org', 'modelarchive.org',
            'proteinensemble.org', 'pdb-ihm.org', 'alphafill.eu'];
        let parsed;
        try {
            parsed = new URL(pdbUrl);
        } catch (_) {
            throw new Error('Invalid PDB URL');
        }
        if (parsed.protocol !== 'https:' || !ALLOWED_HOSTS.includes(parsed.hostname)) {
            throw new Error(`Untrusted PDB source: ${parsed.hostname}`);
        }

        const res = await fetch(pdbUrl);
        if (!res.ok) throw new Error(`Structure fetch failed: ${res.status}`);

        // Content-type sanity check
        const ct = (res.headers.get('content-type') || '').toLowerCase();
        if (ct && !ct.includes('text/') && !ct.includes('chemical/') && !ct.includes('application/octet-stream')) {
            console.warn('[UniProt 3D] Unexpected content-type for PDB:', ct);
        }

        // Size guard — reject responses over 50 MB to prevent memory pressure
        const contentLength = parseInt(res.headers.get('content-length') || '0');
        if (contentLength > 50 * 1024 * 1024) {
            throw new Error(`PDB file too large: ${contentLength} bytes`);
        }

        const pdb = await res.text();
        if (pdb.length > 50 * 1024 * 1024) {
            throw new Error(`PDB data too large: ${pdb.length} chars`);
        }
        this.currentPdbText = pdb;
        // Prefer an explicit format from the structure record (computed models carry one); else
        // infer from the URL. Use includes('.cif') so query strings on model URLs don't break it.
        this.currentFormat = structure?.format || (pdbUrl.toLowerCase().includes('.cif') ? 'mmcif' : 'pdb');
        this.currentStructure = structure;

        this.viewer.removeAllModels();
        this.viewer.removeAllLabels();
        this.viewer.removeAllShapes();

        this.viewer.addModel(pdb, this.currentFormat);
        // Build SIFTS-validated residue index: find which PDB residues are actually modeled
        this._buildObservedResiCache();
        this._sanitizeLigandBonds();
        this.applyCartoonColoring('default');
        this.viewer.zoomTo();
        this.viewer.zoom(1.15); // slightly tighter default framing than zoomTo's exact fit
        this.viewer.render();
    },

    /**
     * Build a Set of PDB residue numbers (resi) that are physically present
     * in the loaded model for this structure's chain. Used to skip unmodeled
     * residues (flexible loops, disordered termini) when placing PTM/variant spheres.
     * This cross-references the SIFTS-derived mappedRanges with actual atom coordinates.
     */
    /**
     * Convert a UniProt position to the PDB author residue number for a given mapped range.
     * For most structures this is a linear offset from pdbStart.  For chimeric / engineered
     * structures where SIFTS has no entry for this protein, api.js attaches a
     * `seqresToAuthor` Map and a `seqresStart` field per range so we can do a direct lookup.
     */
    _resiToPdb(uniprotResi, r, chainId = null) {
        // Prefer the per-chain seqresToAuthor map if available (chimeric multi-chain
        // structures where each chain was decorated independently), then fall back to
        // the structure-level map stored on the primary chain.
        const map = (chainId != null && this.currentStructure?.chainSeqresToAuthor?.[chainId])
            || this.currentStructure?.seqresToAuthor;
        if (map) {
            const seqresPos = (r.seqresStart ?? r.pdbStart) + (uniprotResi - r.uniprotStart);
            return map.get(seqresPos) ?? seqresPos;
        }
        return r.pdbStart + (uniprotResi - r.uniprotStart);
    },

    _buildObservedResiCache() {
        this._observedResi = null;
        this._observedResiByChain = null;
        const structure = this.currentStructure;
        if (!structure || (structure.source === 'AlphaFold' && !structure.isoform) || !structure.mappedRanges?.length) return;
        const model = this.viewer?.getModel();
        if (!model) return;
        if (structure.chainIds?.length > 1) {
            // Multi-chain: build a separate observed-resi Set for each chain so that
            // allChainsAddStyle can correctly validate each chain's author residue numbers
            // (each chain may have independent author numbering, e.g. chimeric structures).
            this._observedResiByChain = {};
            structure.chainIds.forEach(chain => {
                const atoms = model.selectedAtoms({ chain, atom: 'CA' });
                this._observedResiByChain[chain] = atoms?.length ? new Set(atoms.map(a => a.resi)) : null;
            });
            // Also set the primary-chain set for single-chain code paths
            this._observedResi = this._observedResiByChain[structure.chainId] || null;
        } else {
            const chain = structure.chainId;
            const atoms = model.selectedAtoms({ ...(chain ? { chain } : {}), atom: 'CA' });
            if (atoms?.length) this._observedResi = new Set(atoms.map(a => a.resi));
        }
    },

    residueSelector(resi) {
        const structure = this.currentStructure;
        if (!structure || (structure.source === 'AlphaFold' && !structure.isoform) || !structure.mappedRanges?.length) return { resi };
        for (const r of structure.mappedRanges) {
            if (resi >= r.uniprotStart && resi <= r.uniprotEnd) {
                const pdbResi = this._resiToPdb(resi, r);
                // Only return a valid selector if the residue is actually modeled
                if (this._observedResi && !this._observedResi.has(pdbResi)) {
                    return { resi: -999999 }; // In mapped range but unmodeled (flexible loop etc.)
                }
                return r.chainId ? { chain: r.chainId, resi: pdbResi } : { resi: pdbResi };
            }
        }
        return { resi: -999999 }; // Outside all SIFTS segments
    },

    /**
     * Resolve a UniProt position to a {chain, resi} selector for a SPECIFIC chain of a
     * multi-chain structure.  Each subunit may have its own author residue numbering
     * (chainMappings / chainSeqresToAuthor), so clicking a residue on chain D must address
     * chain D's atoms — not the primary chain's copy.  Falls back to residueSelector() when
     * no chain is given.
     */
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

    /**
     * Apply a style to one chain's copy of a UniProt position (per-chain coloring for
     * structure-dependent modes like hotspots / contact hubs, where each subunit may carry
     * a different tier).  Returns true when the style was applied.
     */
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

    /**
     * Applies an addStyle to the correct PDB residue(s) for a UniProt position.
     * Returns true if the style was applied, false if the position is outside the
     * mapped range or maps to an unmodelled (unobserved) residue.
     */
    allChainsAddStyle(resi, style, atomSel = {}) {
        const structure = this.currentStructure;
        if (!structure || (structure.source === 'AlphaFold' && !structure.isoform) || !structure.mappedRanges?.length) {
            this.viewer.addStyle({ resi, ...atomSel }, style);
            return true;
        }
        if (structure.chainIds?.length > 1) {
            // Multi-chain (homodimers / homo-oligomers): use per-chain mappings so that
            // chains with different residue offsets are addressed at the correct PDB resi.
            let applied = false;
            structure.chainIds.forEach(chain => {
                const ranges = structure.chainMappings?.[chain] || structure.mappedRanges;
                const r = ranges.find(mr => resi >= mr.uniprotStart && resi <= mr.uniprotEnd);
                if (!r) return;
                const chainPdbResi = this._resiToPdb(resi, r, chain);
                // Use per-chain observed set so chimeric chains with independent author
                // numbering are validated correctly (not against the primary chain's atoms).
                const chainObs = this._observedResiByChain?.[chain] ?? this._observedResi;
                if (chainObs && !chainObs.has(chainPdbResi)) return;
                this.viewer.addStyle({ chain, resi: chainPdbResi, ...atomSel }, style);
                applied = true;
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
        return false; // resi outside all mapped ranges
    },

    mappedResidues() {
        return this.mappedResiduesForChain(null);
    },

    /**
     * UniProt positions that are actually modelled in a given chain (or the primary chain
     * when chain is null).  Used to render one sequence-ribbon track per chain, each showing
     * which residues are resolved in that specific subunit.
     */
    mappedResiduesForChain(chain) {
        const s = this.currentStructure;
        const ranges = (chain != null && s?.chainMappings?.[chain]) || s?.mappedRanges;
        if (!ranges?.length) return null;
        const obs = (chain != null ? (this._observedResiByChain?.[chain] ?? null) : this._observedResi);
        const residues = [];
        ranges.forEach(r => {
            for (let i = r.uniprotStart; i <= r.uniprotEnd; i++) {
                if (!obs) {
                    residues.push(i); // no cache yet: include all mapped
                } else {
                    const pdbResi = this._resiToPdb(i, r, chain);
                    if (obs.has(pdbResi)) residues.push(i);
                }
            }
        });
        return residues.length > 0 ? residues : null;
    },

    applyCartoonColoring(mode = 'default', context = {}, defer = false) {
        if (!this.viewer) return;
        this.activeColoringMode = mode;
        this._lastColoringContext = context;
        this._inFocusMode = false; // exit focus mode on any full re-color
        this._focusState = null;
        this.viewer.removeAllShapes();
        this.viewer.removeAllLabels();
        this.viewer.setStyle({}, {}); // wipe accumulated sphere/stick addStyle records
        const base = { opacity: 0.82, thickness: 0.25, ribbonWidth: 0.6 };
        this._applyModeStyles(mode, context, base);
        this._drawLigands();
        if (!defer) this.viewer.render();
    },

    /**
     * For AlphaFill models, draw the transplanted cofactors/ligands (HETATM) as licorice —
     * the ligands are the whole point of an AlphaFill model, so showing only the protein cartoon
     * would hide them. setStyle() in _applyModeStyles replaces all per-atom styles, so this must
     * be called AFTER it on every (re)style.
     */
    showLigands: true, // global "show all ligands" toggle (Ligands list All/None)
    excludeIons: false, // when true, hide monatomic ions (and water) — Ligands header toggle

    // Common monatomic-ion CCD codes (so they can be drawn as spheres / optionally hidden).
    ION_CODES: new Set(['NA', 'K', 'LI', 'RB', 'CS', 'MG', 'CA', 'SR', 'BA', 'ZN', 'FE', 'FE2',
        'MN', 'MN3', 'CU', 'CU1', 'CO', 'NI', 'CD', 'HG', 'CL', 'BR', 'IOD', 'FLO', 'F', 'AL',
        'PT', 'AU', 'AG', 'PB', 'SO4', 'PO4']),

    /**
     * Remove spurious bonds that 3Dmol's distance-based bond assignment creates BETWEEN distinct
     * hetero groups (and between a ligand and protein). AlphaFill transplants many alternative
     * ligands into the same site, so their atoms overlap and get cross-bonded — making everything
     * look covalently joined (e.g. a sodium "bonded" to a ligand). We keep only intra-group bonds.
     */
    _sanitizeLigandBonds() {
        const model = this.viewer?.getModel?.();
        const atoms = model?.selectedAtoms?.({}) || [];
        if (!atoms.length) return;
        const byIdx = new Map(atoms.map(a => [a.index, a]));
        for (const a of atoms) {
            if (!a.bonds || !a.bonds.length) continue;
            const keepB = [], keepO = [];
            for (let k = 0; k < a.bonds.length; k++) {
                const b = byIdx.get(a.bonds[k]);
                if (!b) continue;
                const sameGroup = a.chain === b.chain && a.resi === b.resi && a.resn === b.resn;
                // Drop any bond involving a hetero atom that crosses groups (ligand↔ligand,
                // ligand↔protein, ligand↔ion); keep protein backbone and intra-ligand bonds.
                if (((a.hetflag || b.hetflag) && !sameGroup)) continue;
                keepB.push(a.bonds[k]); keepO.push(a.bondOrder ? a.bondOrder[k] : 1);
            }
            a.bonds = keepB; a.bondOrder = keepO;
        }
    },

    /**
     * Show ligands/cofactors (HETATM) as licorice for ANY structure that has them — AlphaFill
     * transplanted cofactors, experimental-structure ligands, etc. Bulk water is always hidden;
     * monatomic ions are hidden when excludeIons is set. Called after every cartoon (re)style.
     */
    _drawLigands() {
        if (this.showLigands === false || !this.viewer) return;
        this.viewer.addStyle({ hetflag: true }, { stick: { radius: 0.18, colorscheme: 'Jmol' }, sphere: { radius: 0.35, colorscheme: 'Jmol' } });
        const hide = ['HOH', 'WAT', 'DOD'];
        if (this.excludeIons) hide.push(...this.ION_CODES);
        this.viewer.setStyle({ resn: hide }, {});
    },

    // Brighten the hovered ligand and label it with its CCD code (cleared on hover-out / move).
    _hoverLigand(atom) {
        if (this.excludeIons && this.ION_CODES.has(atom.resn)) return; // hidden ion
        const key = `${atom.chain ?? ''}|${atom.resi}|${atom.resn}`;
        if (this._hoverLigKey === key) return;
        this._clearLigandHover();
        this._hoverLigKey = key;
        const sel = { resn: atom.resn, resi: atom.resi, ...(atom.chain != null ? { chain: atom.chain } : {}) };
        this.viewer.setStyle(sel, { stick: { radius: 0.32, color: '#ffeb3b' }, sphere: { radius: 0.5, color: '#ffeb3b' } });
        try {
            this._hoverLabel = this.viewer.addLabel(atom.resn, {
                position: { x: atom.x, y: atom.y, z: atom.z },
                backgroundColor: '#11161f', backgroundOpacity: 0.9, fontColor: '#ffeb3b', fontSize: 13, borderThickness: 0,
                screenOffset: { x: 16, y: -16 }, // nudge off the cursor so it isn't hidden by it
                inFront: true,
            });
        } catch (_) {}
        this.viewer.render();
    },

    _clearLigandHover() {
        if (!this._hoverLigKey) return;
        const [chain, resi, resn] = this._hoverLigKey.split('|');
        this.viewer.setStyle({ resn, resi: parseInt(resi, 10), ...(chain ? { chain } : {}) },
            { stick: { radius: 0.18, colorscheme: 'Jmol' }, sphere: { radius: 0.35, colorscheme: 'Jmol' } });
        if (this._hoverLabel) { try { this.viewer.removeLabel(this._hoverLabel); } catch (_) {} this._hoverLabel = null; }
        this._hoverLigKey = null;
        this.viewer.render();
    },

    /** Shared helper: apply cartoon color styles for a given mode (used by applyCartoonColoring and focusResidue). */
    _applyModeStyles(mode, context, base) {
        if (mode === 'plddt') {
            // Match AlphaFold DB pLDDT color scheme exactly
            this.viewer.setStyle({}, { cartoon: { ...base, colorfunc: atom => {
                const b = atom.b;
                if (b >= 90) return '#0053d6';
                if (b >= 70) return '#65cbf3';
                if (b >= 50) return '#ffdb13';
                return '#ff7d45';
            }}});
        } else if (mode === 'bfactor') {
            // Blue (low/rigid) → white (midpoint) → red (high/flexible)
            // Matches legend: #313695 → #f7f7f7 → #d73027
            this.viewer.setStyle({}, { cartoon: { ...base, colorfunc: atom => {
                const v = Math.max(0, Math.min(100, atom.b));
                if (v <= 50) {
                    const t = v / 50;
                    const r = Math.round(49  + (247 - 49)  * t);
                    const g = Math.round(54  + (247 - 54)  * t);
                    const b = Math.round(149 + (247 - 149) * t);
                    return `rgb(${r},${g},${b})`;
                } else {
                    const t = (v - 50) / 50;
                    const r = Math.round(247 + (215 - 247) * t);
                    const g = Math.round(247 + (48  - 247) * t);
                    const b = Math.round(247 + (39  - 247) * t);
                    return `rgb(${r},${g},${b})`;
                }
            }}});
        } else if (mode === 'hotspots') {
            this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
            const tierColors = { strong: '#b71c1c', moderate: '#e64a19', weak: '#ffa726' };
            this._applyTierColoring(context.hotspots, context.hotspotsByChain, tierColors, base);
        } else if (mode === 'distantContacts') {
            this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
            const dcColors = { strong: '#6a1b9a', moderate: '#ab47bc' };
            this._applyTierColoring(context.distantContacts, context.distantContactsByChain, dcColors, base);
        } else if (mode === 'alphaMissense') {
            this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
            (context.alphaMissense || new Map()).forEach((d, pos) => {
                const avg = d.avg;
                const color = avg >= 0.78 ? '#b71c1c'
                            : avg >= 0.564 ? '#e06666'
                            : avg >= 0.34 ? '#b9c2cf'
                            : '#3d85c8';
                this.allChainsAddStyle(pos, { cartoon: { ...base, color } });
            });
        } else if (mode === 'residueBurden') {
            this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
            (context.residueBurden || new Set()).forEach(pos => this.allChainsAddStyle(pos, { cartoon: { ...base, color: '#e65100' } }));
        } else if (mode === 'topology' || mode === 'domains') {
            // Colour the cartoon by per-position segment map: membrane-topology segments, or (in the
            // Family & Domains window) the domain/region/repeat ranges. Same machinery either way.
            const posColor = context.topologyByPos instanceof Map ? context.topologyByPos
                : context.domainByPos instanceof Map ? context.domainByPos : new Map();
            const s = this.currentStructure;
            const isAF = !s || (s.source === 'AlphaFold' && !s.isoform) || !s.mappedRanges?.length;
            // Restrict colouring to chains that actually belong to THIS protein. In a hetero-complex
            // the partner chains aren't in chainMappings, so _reverseResidueMapForChain would fall
            // back to our primary mappedRanges and mis-map their residues onto our topology — paint
            // them neutral instead. (Homo-oligomers list every copy in chainIds, so they're kept.)
            const ourChains = s?.chainIds?.length ? new Set(s.chainIds)
                            : (s?.chainId ? new Set([s.chainId]) : null);
            const reverseCache = new Map();
            const toUni = (chain, resi) => {
                if (!reverseCache.has(chain)) reverseCache.set(chain, this._reverseResidueMapForChain(chain));
                return reverseCache.get(chain).get(resi);
            };
            this.viewer.setStyle({}, { cartoon: { ...base, colorfunc: atom => {
                if (!isAF && ourChains && atom.chain != null && !ourChains.has(atom.chain)) return '#b9c2cf';
                const uni = isAF ? atom.resi : toUni(atom.chain, atom.resi);
                return posColor.get(uni) || '#b9c2cf';
            } } });
        } else if (mode === 'prism') {
            // Constraint-pocket mode: colour candidate residues by geometric class.
            this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
            const catColors = { pocket: '#00897b', exposed: '#8e24aa' };
            if (context.pocketByPos instanceof Map) {
                context.pocketByPos.forEach((info, pos) => this.allChainsAddStyle(pos, { cartoon: { ...base, color: catColors[info.cat] || '#00897b' } }));
            }
        } else {
            this.viewer.setStyle({}, { cartoon: { ...base, color: '#00bcd4' } });
        }
    },

    /**
     * Color a structure-dependent tier map.  For multi-chain structures each chain's copy is
     * colored from ITS OWN per-chain tier map (byChain) so subunits that differ structurally
     * show different tiers; single-chain/AlphaFold uses the merged map across all chains.
     */
    _applyTierColoring(merged, byChain, tierColors, base) {
        const structure = this.currentStructure;
        if (structure?.chainIds?.length > 1 && byChain) {
            byChain.forEach((tierMap, chain) => {
                tierMap.forEach((tier, pos) => {
                    this.chainAddStyle(chain, pos, { cartoon: { ...base, color: tierColors[tier] || '#b9c2cf' } });
                });
            });
            return;
        }
        (merged || new Map()).forEach((tier, pos) => {
            this.allChainsAddStyle(pos, { cartoon: { ...base, color: tierColors[tier] || '#b9c2cf' } });
        });
    },

    /**
     * Render PTM residues as Cα VDW spheres.
     * For disulfide bonds, both begin and end positions are rendered.
     */
    /**
     * Draw "Site" annotation spheres (amber) for the given sites at positions not already
     * occupied by a PTM/variant sphere. Adds them to the active + hover maps so they persist in
     * focus mode and show their annotation on hover/click.
     */
    _drawSiteSpheres(sites, hoverMap, active) {
        if (!sites || !sites.length) return;
        sites.forEach(site => {
            const positions = [site.position, site.endPosition].filter((p, i, a) => p && a.indexOf(p) === i);
            positions.forEach(pos => {
                if (hoverMap.has(pos)) return; // PTM/variant sphere already here
                const placed = this.allChainsAddStyle(pos, { sphere: { radius: 1.8, color: site.color, opacity: 0.92 } }, { atom: 'CA' });
                if (!placed) return;
                hoverMap.set(pos, { position: pos, color: site.color, isSite: true, description: site.description, category: 'Site' });
                active.set(pos, site.color);
            });
        });
    },

    showPTMs(ptms, ptmGroups, sites = [], extras = []) {
        this._selectedResi = null;
        this._lastRender = () => this.showPTMs(ptms, ptmGroups, sites, extras);
        const hoverMap = new Map();
        const active = new Map(); // uniProtPos → color, for zoom-mode sphere persistence
        let count = 0;

        ptms.forEach(ptm => {
            const g = ptmGroups[ptm.category];
            if (!g || !g.visible) return;
            if (ptm.visible === false) return;

            // Only count + hover this PTM when the sphere is actually placed in the 3D model
            // (allChainsAddStyle returns false for positions outside the mapped range or
            // in unmodelled/disordered regions so they don't show in the current structure).
            const placed = this.allChainsAddStyle(ptm.position, { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }, { atom: 'CA' });
            if (!placed) return;

            hoverMap.set(ptm.position, ptm);
            active.set(ptm.position, ptm.color);
            count++;

            // For disulfide bonds, also render the end position
            if (ptm.endPosition && ptm.endPosition !== ptm.position) {
                if (this.allChainsAddStyle(ptm.endPosition, { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }, { atom: 'CA' })) {
                    active.set(ptm.endPosition, ptm.color);
                }
                // Create a hover entry for the end position too
                if (!hoverMap.has(ptm.endPosition)) {
                    hoverMap.set(ptm.endPosition, {
                        ...ptm,
                        position: ptm.endPosition,
                        description: ptm.description + ` (bonded to ${ptm.position})`
                    });
                }
            }
        });

        this._drawSiteSpheres(sites, hoverMap, active);
        // Optional extra spheres (e.g. the secondary "Disease variants" group toggled on in the
        // PTM window). Drawn after PTMs/sites so they share the same hover + persistence machinery.
        (extras || []).forEach(sp => {
            if (sp.position == null) return;
            if (this.allChainsAddStyle(sp.position, { sphere: { radius: 1.8, color: sp.color, opacity: 0.92 } }, { atom: 'CA' })) {
                active.set(sp.position, sp.color);
                if (sp.hover && !hoverMap.has(sp.position)) hoverMap.set(sp.position, sp.hover);
            }
        });
        this._activeSpheres = active;
        this._bindHover(hoverMap, 'ptm');
        this.viewer.render();
        return count;
    },

    /**
     * Generic CA-sphere overlay used by the Functional-features and Family & Domains windows.
     * Assumes the cartoon has already been (re)coloured via applyCartoonColoring(..., defer=true);
     * adds one CA sphere per annotation (plus the end position for ranged ones) and binds hover.
     * `spheres`: [{ position, endPosition?, color, hover }]. Returns the number actually placed
     * (positions outside the mapped/modelled range are skipped, like showPTMs).
     */
    showAnnotationSpheres(spheres) {
        if (!this.viewer) return 0;
        this._selectedResi = null;
        this._lastRender = () => this.showAnnotationSpheres(spheres);
        const hoverMap = new Map();
        const active = new Map();
        let count = 0;
        (spheres || []).forEach(sp => {
            if (sp.position == null) return;
            const placed = this.allChainsAddStyle(sp.position, { sphere: { radius: 1.8, color: sp.color, opacity: 0.92 } }, { atom: 'CA' });
            if (!placed) return;
            active.set(sp.position, sp.color);
            if (sp.hover) hoverMap.set(sp.position, sp.hover);
            count++;
            if (sp.endPosition && sp.endPosition !== sp.position) {
                if (this.allChainsAddStyle(sp.endPosition, { sphere: { radius: 1.8, color: sp.color, opacity: 0.92 } }, { atom: 'CA' }))
                    active.set(sp.endPosition, sp.color);
                if (sp.hover && !hoverMap.has(sp.endPosition)) hoverMap.set(sp.endPosition, { ...sp.hover, position: sp.endPosition });
            }
        });
        this._activeSpheres = active;
        this._bindHover(hoverMap, 'feature');
        this.viewer.render();
        return count;
    },

    /**
     * PTM sphere refresh — re-applies base cartoon (clearing accumulated addStyle layers)
     * then re-renders only the currently visible PTM spheres.
     * Returns false when in focus mode so the caller can do a full applyMode() instead.
     */
    refreshPTMDisplay(ptms, ptmGroups, sites = [], extras = []) {
        if (!this.viewer || this._inFocusMode) return false;
        this.viewer.removeAllLabels();
        // Hard-clear ALL per-atom style records first.  setStyle({}, {cartoon}) replaces the
        // cartoon property but may leave sphere properties from previous addStyle calls intact
        // in 3Dmol's internal style list, causing spheres to accumulate brightness each call.
        // The empty setStyle({}, {}) wipes every style record from every atom before we
        // re-apply the cartoon and then add spheres cleanly.
        this.viewer.setStyle({}, {});
        const base = { opacity: 0.82, thickness: 0.25, ribbonWidth: 0.6 };
        this._applyModeStyles(this.activeColoringMode, this._lastColoringContext || {}, base);
        this._drawLigands();
        return this.showPTMs(ptms, ptmGroups, sites, extras);
    },

    /**
     * Render variant spheres, plus (optionally) a set of co-displayed PTMs so the user can
     * see post-translational modifications alongside variants in the Disease & Variants view.
     * Where a residue carries both, the variant sphere wins (its severity colour is kept) but
     * the PTM still shows in that residue's hover/details.
     */
    showVariants(filtered, coPtms = [], sites = []) {
        this._selectedResi = null;
        this._lastRender = () => this.showVariants(filtered, coPtms, sites);
        const severity = [
            'Likely pathogenic or pathogenic',
            'Predicted deleterious',
            'Uncertain significance',
            'Likely benign or benign'
        ];
        const posMap = new Map();

        filtered.forEach(v => {
            const ex = posMap.get(v.position);
            if (!ex) {
                posMap.set(v.position, { position: v.position, color: v.consequenceColor, topConsequence: v.consequence, variants: [v] });
            } else {
                ex.variants.push(v);
                const ei = severity.indexOf(ex.topConsequence);
                const ni = severity.indexOf(v.consequence);
                if (ni >= 0 && (ei < 0 || ni < ei)) {
                    ex.topConsequence = v.consequence;
                    ex.color = v.consequenceColor;
                }
            }
        });

        const hoverMap = new Map();
        const active = new Map(); // uniProtPos → color, for zoom-mode sphere persistence
        let placedPosCount = 0;
        let placedVarCount = 0;
        posMap.forEach((d, pos) => {
            const placed = this.allChainsAddStyle(pos, { sphere: { radius: 1.8, color: d.color, opacity: 0.92 } }, { atom: 'CA' });
            if (!placed) return;
            hoverMap.set(pos, d);
            active.set(pos, d.color);
            placedPosCount++;
            placedVarCount += d.variants.length;
        });

        // Co-displayed PTM spheres (Disease & Variants "show PTMs too" option).
        let placedPtmCount = 0;
        coPtms.forEach(ptm => {
            if (posMap.has(ptm.position)) return; // variant sphere already occupies this Cα
            const placed = this.allChainsAddStyle(ptm.position, { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }, { atom: 'CA' });
            if (!placed) return;
            if (!hoverMap.has(ptm.position)) hoverMap.set(ptm.position, { position: ptm.position, color: ptm.color, category: ptm.category });
            active.set(ptm.position, ptm.color);
            placedPtmCount++;
            if (ptm.endPosition && ptm.endPosition !== ptm.position && !posMap.has(ptm.endPosition)) {
                if (this.allChainsAddStyle(ptm.endPosition, { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }, { atom: 'CA' })) {
                    active.set(ptm.endPosition, ptm.color);
                }
            }
        });

        this._drawSiteSpheres(sites, hoverMap, active);
        this._activeSpheres = active;
        this._bindHover(hoverMap, 'variant');
        this.viewer.render();
        return { posCount: placedPosCount, varCount: placedVarCount, ptmCount: placedPtmCount };
    },

    _bindHover(map, mode) {
        const self = this;
        // Resolve the atom's author residue number to a UniProt position using the mapping of
        // the chain the atom actually belongs to, so clicks on any subunit are interpreted
        // (and later focused) on that subunit rather than the primary chain.
        const resolveUni = atom => self._reverseResidueMapForChain(atom.chain).get(atom.resi) || atom.resi;
        // Hover stays annotated-only: it fires on every mousemove, so rebuilding the tooltip
        // (which scans all PTMs/variants) for bare backbone residues would jank large proteins.
        this.viewer.setHoverable({}, true,
            (atom, _v, ev) => {
                if (!atom?.resi) return;
                // Ligand hover: brighten that ligand and label it with its CCD code.
                if (atom.hetflag && atom.resn !== 'HOH' && atom.resn !== 'WAT') {
                    self._hoverLigand(atom);
                    if (self.hoverCb) self.hoverCb(null, mode, null);
                    return;
                }
                if (self._hoverLigKey) self._clearLigandHover();
                const uniResi = resolveUni(atom);
                const d = map.get(uniResi);
                if (d && self.hoverCb) self.hoverCb(d, mode, ev, atom.chain);
            },
            () => { self._clearLigandHover(); if (self.hoverCb) self.hoverCb(null, mode, null); }
        );
        // Click is bound on the WHOLE model ({}) so EVERY residue is clickable — not just
        // annotated ones.  Unannotated residues fall back to a bare { position } payload so the
        // cartoon is always focusable/zoomable regardless of whether it carries a PTM/variant.
        this.viewer.setClickable({}, true,
            (atom) => {
                if (!atom?.resi) return;
                // Ligand/cofactor (HETATM) — route to the ligand handler instead of treating it
                // as a protein residue (which would map it to a bogus UniProt position and show
                // variants / AlphaMissense / hotspot flags for it).
                if (atom.hetflag && atom.resn !== 'HOH' && self.ligandClickCb) {
                    self.ligandClickCb({ resn: atom.resn, resi: atom.resi, chain: atom.chain ?? null });
                    return;
                }
                const uniResi = resolveUni(atom);
                const d = map.get(uniResi);
                if (self.clickCb) self.clickCb(d || { position: uniResi }, mode, atom.chain);
            }
        );
    },

    /**
     * Focus on a residue: keep cartoon visible, show selected residue
     * and its 5Å neighbors as ball-and-stick overlay.
     * No interaction lines. No residue label.
     */
    focusResidue(resi, chain = null, annotations = {}, opts = {}) {
        if (!this.viewer) return;
        this.viewer.removeAllShapes();
        this.viewer.removeAllLabels();
        this._selectedResi = resi;

        // Resolve to the clicked chain's copy so focus lands on the right subunit.
        const selector = this.residueSelectorForChain(resi, chain);
        const selChain = selector.chain ?? null;
        const pdbResi = selector.resi;
        const model = this.viewer.getModel();
        // Search the WHOLE model (all chains) so interface neighbours from other subunits
        // are still captured, then restrict the "selected" atoms to the clicked chain.
        const allAtoms = model.selectedAtoms({});
        const selAtoms = allAtoms.filter(a => a.resi === pdbResi && (selChain == null || a.chain === selChain));
        if (selAtoms.length === 0) return;

        const CONTACT_DIST = 5.0;
        // Residue numbers repeat across chains, so track neighbours by chain+resi.
        const keyOf = (c, r) => `${c ?? ''}|${r}`;
        const nearby = new Map(); // key → { chain, resi, het }
        nearby.set(keyOf(selChain, pdbResi), { chain: selChain, resi: pdbResi, het: false });
        allAtoms.forEach(a => {
            if (a.resi === pdbResi && a.chain === selChain) return;
            for (const sa of selAtoms) {
                const dx = a.x - sa.x, dy = a.y - sa.y, dz = a.z - sa.z;
                if (Math.sqrt(dx*dx + dy*dy + dz*dz) < CONTACT_DIST) {
                    nearby.set(keyOf(a.chain, a.resi), { chain: a.chain, resi: a.resi, het: !!a.hetflag });
                    break;
                }
            }
        });

        const annotated = annotations.annotatedResidues || new Map();
        const reverseCache = new Map(); // chain → reverse map (built lazily, reused)
        const toUni = (c, r) => {
            if (!reverseCache.has(c)) reverseCache.set(c, this._reverseResidueMapForChain(c));
            return reverseCache.get(c).get(r) || r;
        };

        // Focused neighbourhood as UniProt positions — these become sticks below, so the rest
        // of the annotation spheres should stay as spheres (the user wants other PTMs/variants
        // to remain visible while zoomed into one of them).
        const nearbyUni = new Set(Array.from(nearby.values()).filter(n => !n.het).map(n => toUni(n.chain, n.resi)));

        // Dim the cartoon to background opacity.  Re-applied on every focus (not only the
        // first) because setStyle replaces the per-atom style layers — this clears the
        // previous click's sticks/context spheres so they don't accumulate, and lets us
        // re-draw the surrounding annotation spheres cleanly each time.
        this._inFocusMode = true;
        const focusBase = { opacity: 0.42, thickness: 0.2, ribbonWidth: 0.5 };
        this._applyModeStyles(this.activeColoringMode, this._lastColoringContext, focusBase);
        // NB: do NOT draw all ligands here — in a residue zoom view, showing every ligand clutters
        // the pocket. Ligands that fall within the focused neighbourhood are already drawn as
        // sticks below (via the `nearby` set).

        // Keep every OTHER annotation sphere visible — controlled by opts.showOtherSpheres (default true).
        if (opts.showOtherSpheres !== false) {
            (this._activeSpheres || new Map()).forEach((color, uniPos) => {
                if (nearbyUni.has(uniPos)) return; // shown as a stick instead
                this.allChainsAddStyle(uniPos, { sphere: { radius: 1.8, color, opacity: 0.92 } }, { atom: 'CA' });
            });
        }

        // Show nearby residues as thin sticks (each on its own chain's copy).
        // Capture them for the session export (PyMOL/VMD) so the zoom-in is reproduced faithfully.
        const focusSticks = [];
        nearby.forEach(({ chain: nc, resi: nr, het }) => {
            if (nc === selChain && nr === pdbResi) return;
            if (het) {
                // Ligand/cofactor near the residue — keep default element colours (grey C, red O,
                // …), never an annotation colour, and add small spheres so single-atom ions show.
                this.viewer.addStyle(nc != null ? { chain: nc, resi: nr, hetflag: true } : { resi: nr, hetflag: true },
                    { stick: { radius: 0.15, colorscheme: 'Jmol' }, sphere: { radius: 0.3, colorscheme: 'Jmol' } });
                focusSticks.push({ chain: nc || '', resi: nr, het: true, color: null });
                return;
            }
            const uniResi = toUni(nc, nr);
            const color = annotated.get(uniResi)?.color;
            this.viewer.addStyle(
                nc != null ? { chain: nc, resi: nr } : { resi: nr },
                { stick: color ? { radius: 0.15, color, opacity: 0.8 } : { radius: 0.12, colorscheme: 'Jmol', opacity: 0.6 } }
            );
            focusSticks.push({ chain: nc || '', resi: nr, het: false, color: color || null });
        });
        this._focusState = { selChain: selChain || '', pdbResi, sticks: focusSticks };

        // Show selected residue as ball-and-stick (prominent)
        this.viewer.addStyle(
            selector,
            {
                stick: { radius: 0.2, colorscheme: 'Jmol' },
                sphere: { radius: 0.4, colorscheme: 'Jmol' }
            }
        );

        // Bind hover/click BEFORE render so callbacks survive geometry rebuild
        const focusHoverMap = new Map();
        nearby.forEach(({ chain: nc, resi: nr }) => {
            const uni = toUni(nc, nr);
            focusHoverMap.set(uni, { position: uni, color: annotated.get(uni)?.color });
        });
        this._bindHover(focusHoverMap, 'focus');

        // Zoom to the local pocket on the selected chain (frame the clicked subunit).
        // opts.rezoom === false keeps the current camera (e.g. when only toggling sphere
        // visibility) so the view doesn't perform a distracting zoom animation.
        if (opts.rezoom !== false) {
            const zoomResis = Array.from(nearby.values())
                .filter(n => n.chain === selChain)
                .map(n => n.resi);
            this._zoomToWithMargin(selChain != null ? { chain: selChain, resi: zoomResis } : { resi: zoomResis });
        }
        this.viewer.render();

        return nearbyUni;
    },

    /**
     * List the distinct ligand/cofactor instances (HETATM groups, excluding water) present in
     * the loaded model. Returns [{ resn, resi, chain }] — one per (chain, resi).
     */
    enumerateLigands() {
        const model = this.viewer?.getModel?.();
        if (!model) return [];
        const seen = new Map();
        (model.selectedAtoms({ hetflag: true }) || []).forEach(a => {
            if (a.resn === 'HOH' || a.resn === 'WAT') return;
            const key = `${a.chain ?? ''}|${a.resi}`;
            if (!seen.has(key)) seen.set(key, { resn: a.resn, resi: a.resi, chain: a.chain ?? null });
        });
        return [...seen.values()];
    },

    /**
     * Map each protein residue (UniProt position) to the set of ligand CCD codes whose atoms lie
     * within maxDist Å of it. Used to add a "nearby_ligands" column to the CSV export.
     */
    ligandContactsByResidue(maxDist = 5) {
        const model = this.viewer?.getModel?.();
        if (!model) return new Map();
        const ligAtoms = (model.selectedAtoms({ hetflag: true }) || []).filter(a => a.resn !== 'HOH' && a.resn !== 'WAT');
        if (!ligAtoms.length) return new Map();
        const d2 = maxDist * maxDist;
        const reverseCache = new Map();
        const toUni = (c, r) => { if (!reverseCache.has(c)) reverseCache.set(c, this._reverseResidueMapForChain(c)); return reverseCache.get(c).get(r); };
        const out = new Map(); // uniPos → Set<CCD>
        (model.selectedAtoms({}) || []).forEach(a => {
            if (a.hetflag) return;
            let uni, computed = false;
            for (const la of ligAtoms) {
                const dx = a.x - la.x, dy = a.y - la.y, dz = a.z - la.z;
                if (dx * dx + dy * dy + dz * dz <= d2) {
                    if (!computed) { uni = toUni(a.chain ?? null, a.resi); computed = true; }
                    if (uni == null) break;
                    if (!out.has(uni)) out.set(uni, new Set());
                    out.get(uni).add(la.resn);
                }
            }
        });
        return out;
    },

    /**
     * Focus a ligand: dim the cartoon, HIDE all other hetero groups, show the selected ligand as
     * prominent licorice, draw its nearby protein residues as sticks, and zoom to it. Returns the
     * set of nearby protein residues as UniProt positions (for the detail panel's "Nearby" list).
     */
    focusLigand(resn, resi, chain = null, opts = {}) {
        if (!this.viewer) return new Set();
        this.viewer.removeAllShapes();
        this.viewer.removeAllLabels();
        this.clearProximityLines();
        this._selectedResi = null;
        this._focusState = null; // residue-focus export state doesn't apply to a ligand focus
        const prevActive = this._activeSpheres || new Map(); // PTM/variant/site spheres to optionally keep
        const model = this.viewer.getModel();
        const ligSel = { resn, resi, ...(chain != null ? { chain } : {}) };
        const ligAtoms = (model.selectedAtoms(ligSel) || []).filter(a => a.hetflag);
        if (!ligAtoms.length) return new Set();

        // Nearby PROTEIN residues within 5 Å of any ligand atom.
        const CONTACT2 = 25; // 5 Å squared
        const nearby = new Map(); // chain|resi → {chain, resi}
        (model.selectedAtoms({}) || []).forEach(a => {
            if (a.hetflag) return; // protein atoms only
            for (const la of ligAtoms) {
                const dx = a.x - la.x, dy = a.y - la.y, dz = a.z - la.z;
                if (dx * dx + dy * dy + dz * dz < CONTACT2) { nearby.set(`${a.chain ?? ''}|${a.resi}`, { chain: a.chain ?? null, resi: a.resi }); break; }
            }
        });

        // Dim cartoon; setStyle({}, {}) wipes prior styles so OTHER ligands disappear (they get no
        // style here), satisfying "selecting a ligand hides the others".
        this._inFocusMode = true;
        this.viewer.setStyle({}, {});
        const focusBase = { opacity: 0.42, thickness: 0.2, ribbonWidth: 0.5 };
        this._applyModeStyles(this.activeColoringMode, this._lastColoringContext, focusBase);
        // Selected ligand — prominent ball-and-stick (default element colours).
        this.viewer.addStyle(ligSel, { stick: { radius: 0.25, colorscheme: 'Jmol' }, sphere: { radius: 0.45, colorscheme: 'Jmol' } });

        // Map nearby protein residues to UniProt positions (for colouring + the "Nearby" list).
        const reverseCache = new Map();
        const toUni = (c, r) => { if (!reverseCache.has(c)) reverseCache.set(c, this._reverseResidueMapForChain(c)); return reverseCache.get(c).get(r) || r; };
        const nearbyUni = new Set([...nearby.values()].map(n => toUni(n.chain, n.resi)));

        // Nearby protein residues — thin sticks, coloured by their annotation (variant/PTM) when
        // one applies, so the pocket's disease residues stand out (same as residue focus).
        const annotated = opts.annotatedResidues || new Map();
        nearby.forEach(({ chain: nc, resi: nr }) => {
            const color = annotated.get(toUni(nc, nr))?.color;
            this.viewer.addStyle(nc != null ? { chain: nc, resi: nr } : { resi: nr },
                { stick: color ? { radius: 0.15, color, opacity: 0.85 } : { radius: 0.15, colorscheme: 'Jmol', opacity: 0.8 } });
        });

        // Keep the other PTM/variant/site annotation spheres visible (toggleable), exactly like
        // residue focus — so the header sphere toggle does something during ligand focus too.
        if (opts.showOtherSpheres !== false) {
            prevActive.forEach((color, uniPos) => {
                if (nearbyUni.has(uniPos)) return;
                this.allChainsAddStyle(uniPos, { sphere: { radius: 1.8, color, opacity: 0.92 } }, { atom: 'CA' });
            });
        }
        this._activeSpheres = prevActive; // retain so toggling re-draws them

        const focusHoverMap = new Map();
        nearby.forEach(({ chain: nc, resi: nr }) => { const u = toUni(nc, nr); focusHoverMap.set(u, { position: u }); });
        this._bindHover(focusHoverMap, 'focus');

        if (opts.rezoom !== false) this._zoomToWithMargin(ligSel);
        this.viewer.render();
        return nearbyUni;
    },

    /**
     * author-residue-number → UniProt-position map for ONE chain (or the primary/whole
     * structure when chain is null).  Honours per-chain mappings and chimeric seqres→author
     * maps so a click on any subunit resolves to the correct UniProt position.
     */
    _reverseResidueMapForChain(chain = null) {
        const out = new Map();
        const s = this.currentStructure;
        if (!s) return out;
        const ranges = (chain != null && s.chainMappings?.[chain]) || s.mappedRanges || [];
        const map = (chain != null && s.chainSeqresToAuthor?.[chain]) || s.seqresToAuthor || null;
        const obs = (chain != null ? (this._observedResiByChain?.[chain] ?? this._observedResi) : this._observedResi);
        ranges.forEach(r => {
            if (map && r.seqresStart != null) {
                // Non-linear mapping (e.g. chimeric structures like 6CDU):
                // invert seqresToAuthor over this range to get correct UniProt positions.
                const seqresEnd = r.seqresStart + (r.uniprotEnd - r.uniprotStart);
                map.forEach((author, seqres) => {
                    if (seqres < r.seqresStart || seqres > seqresEnd) return;
                    if (obs && !obs.has(author)) return;
                    out.set(author, r.uniprotStart + (seqres - r.seqresStart));
                });
            } else {
                for (let pdb = r.pdbStart; pdb <= r.pdbEnd; pdb++) {
                    if (obs && !obs.has(pdb)) continue; // skip unmodeled residues
                    out.set(pdb, r.uniprotStart + (pdb - r.pdbStart));
                }
            }
        });
        return out;
    },

    /**
     * Extract one record per modelled Cα for the constraint-pocket analysis:
     *   { uniPos, chain, resi, ca:{x,y,z} }
     * EVERY Cα in the model is returned — including partner-protein chains and other copies of
     * the entry protein — so the analysis can account for burial at subunit/partner interfaces.
     * `uniPos` is the UniProt position for residues of the entry protein, and null for
     * partner-protein / unmapped atoms, which serve as burial context only (never scored).
     * `resi` is the PDB author number (used for the sequential pseudo-Cβ direction).
     */
    residueGeometry() {
        const model = this.viewer?.getModel?.();
        if (!model) return [];
        const s = this.currentStructure;
        // Chains belonging to OUR protein (the entry).  Everything else is partner context.
        // null ⇒ no chain info (AlphaFold / single unmapped) ⇒ treat all chains as ours.
        const ourChains = s?.chainIds?.length ? new Set(s.chainIds)
                        : (s?.chainId ? new Set([s.chainId]) : null);
        const reverseCache = new Map(); // chain → (authorResi → uniPos)
        const toUni = (chain, resi) => {
            if (!reverseCache.has(chain)) reverseCache.set(chain, this._reverseResidueMapForChain(chain));
            return reverseCache.get(chain).get(resi);
        };
        const seen = new Set(); // chain|resi — one CA per residue
        const out = [];
        (model.selectedAtoms({ atom: 'CA' }) || []).forEach(a => {
            const key = `${a.chain}|${a.resi}`;
            if (seen.has(key)) return;
            seen.add(key);
            const isOurs = !ourChains || ourChains.has(a.chain);
            const uni = isOurs ? toUni(a.chain, a.resi) : null;
            out.push({ uniPos: uni != null ? uni : null, chain: a.chain ?? null, resi: a.resi, ca: { x: a.x, y: a.y, z: a.z } });
        });
        return out;
    },

    computeActualCoverage(structure) {
        if (!this.viewer || !structure?.mappedRanges?.length) return null;
        const model = this.viewer.getModel();
        if (!model) return null;
        const chain = structure.chainId;
        const atoms = model.selectedAtoms(chain ? { chain } : {});
        if (!atoms || atoms.length === 0) return null;
        const observedPdbResi = new Set(atoms.map(a => a.resi));
        let count = 0;
        structure.mappedRanges.forEach(r => {
            for (let uni = r.uniprotStart; uni <= r.uniprotEnd; uni++) {
                const pdb = this._resiToPdb(uni, r);
                if (observedPdbResi.has(pdb)) count++;
            }
        });
        return count;
    },

    // ── Session export ─────────────────────────────────────────────────────────────────────────
    // The cartoon background colour shared by every "highlight a subset" mode.
    _BG_COLOR: '#b9c2cf',

    // Tier colour for a residue, matching _applyTierColoring: per-chain map for multi-chain
    // structures, merged map otherwise.
    _tierForAtom(merged, byChain, chain, uni) {
        const s = this.currentStructure;
        if (s?.chainIds?.length > 1 && byChain) {
            const m = byChain.get(chain);
            return m ? m.get(uni) : undefined;
        }
        return (merged || new Map()).get(uni);
    },

    // The exact cartoon colour the viewer paints a given CA atom, mirroring _applyModeStyles.
    _cartoonColorForAtom(atom, mode, ctx, isAF, ourChains, toUni) {
        const BG = this._BG_COLOR;
        if (mode === 'plddt') {
            const b = atom.b;
            return b >= 90 ? '#0053d6' : b >= 70 ? '#65cbf3' : b >= 50 ? '#ffdb13' : '#ff7d45';
        }
        if (mode === 'bfactor') {
            const v = Math.max(0, Math.min(100, atom.b));
            if (v <= 50) { const t = v / 50; return `rgb(${Math.round(49 + (247 - 49) * t)},${Math.round(54 + (247 - 54) * t)},${Math.round(149 + (247 - 149) * t)})`; }
            const t = (v - 50) / 50; return `rgb(${Math.round(247 + (215 - 247) * t)},${Math.round(247 + (48 - 247) * t)},${Math.round(247 + (39 - 247) * t)})`;
        }
        const uni = isAF ? atom.resi : toUni(atom.chain, atom.resi);
        if (mode === 'hotspots') {
            const tc = { strong: '#b71c1c', moderate: '#e64a19', weak: '#ffa726' };
            return tc[this._tierForAtom(ctx.hotspots, ctx.hotspotsByChain, atom.chain, uni)] || BG;
        }
        if (mode === 'distantContacts') {
            const dc = { strong: '#6a1b9a', moderate: '#ab47bc' };
            return dc[this._tierForAtom(ctx.distantContacts, ctx.distantContactsByChain, atom.chain, uni)] || BG;
        }
        if (mode === 'alphaMissense') {
            const d = (ctx.alphaMissense || new Map()).get(uni);
            if (!d) return BG;
            return d.avg >= 0.78 ? '#b71c1c' : d.avg >= 0.564 ? '#e06666' : d.avg >= 0.34 ? BG : '#3d85c8';
        }
        if (mode === 'residueBurden') return (ctx.residueBurden instanceof Set && ctx.residueBurden.has(uni)) ? '#e65100' : BG;
        if (mode === 'prism') {
            const cc = { pocket: '#00897b', exposed: '#8e24aa' };
            const info = ctx.pocketByPos instanceof Map ? ctx.pocketByPos.get(uni) : null;
            return info ? (cc[info.cat] || '#00897b') : BG;
        }
        if (mode === 'topology' || mode === 'domains') {
            if (!isAF && ourChains && atom.chain != null && !ourChains.has(atom.chain)) return BG;
            const posColor = ctx.topologyByPos instanceof Map ? ctx.topologyByPos
                : ctx.domainByPos instanceof Map ? ctx.domainByPos : new Map();
            return posColor.get(uni) || BG;
        }
        return '#00bcd4';
    },

    /**
     * Capture the on-screen scene as a plain, serialisable object so it can be reproduced in
     * PyMOL/VMD: per-residue cartoon colours, annotation Cα spheres, and ligand groups — all in
     * author (PDB) chain+residue numbering, which is what PyMOL/VMD selections use. Reproduces the
     * annotated overview (cartoon + every visible PTM/variant/site/domain sphere + ligands); when
     * zoomed into a residue the same full overview is exported.
     */
    getSceneState() {
        if (!this.viewer || !this.currentPdbText) return null;
        const model = this.viewer.getModel?.();
        if (!model) return null;
        const s = this.currentStructure;
        const isAF = !s || (s.source === 'AlphaFold' && !s.isoform) || !s.mappedRanges?.length;
        const ourChains = s?.chainIds?.length ? new Set(s.chainIds) : (s?.chainId ? new Set([s.chainId]) : null);
        const reverseCache = new Map();
        const toUni = (chain, resi) => {
            if (isAF) return resi;
            if (!reverseCache.has(chain)) reverseCache.set(chain, this._reverseResidueMapForChain(chain));
            return reverseCache.get(chain).get(resi);
        };

        const mode = this.activeColoringMode || 'default';
        const ctx = this._lastColoringContext || {};
        const base = mode === 'default' ? '#00bcd4' : this._BG_COLOR;

        const caAtoms = (model.selectedAtoms({ atom: 'CA' }) || []).filter(a => !a.hetflag);
        // Cartoon: only residues whose colour differs from `base` (keeps the script compact).
        const cartoon = [];
        caAtoms.forEach(a => {
            const color = this._cartoonColorForAtom(a, mode, ctx, isAF, ourChains, toUni);
            if (color && color.toLowerCase() !== base.toLowerCase()) cartoon.push({ chain: a.chain || '', resi: a.resi, color });
        });

        // Focus (zoom-in) state: selected residue + 5 Å neighbourhood drawn as sticks. Those
        // residues are NOT spheres on screen, so exclude them from the sphere list below.
        const focus = (this._inFocusMode && this._focusState) ? this._focusState : null;
        const stickKeys = new Set();
        if (focus) {
            focus.sticks.filter(x => !x.het).forEach(x => stickKeys.add(`${x.chain}|${x.resi}`));
            stickKeys.add(`${focus.selChain}|${focus.pdbResi}`);
        }

        // Annotation Cα spheres: _activeSpheres is uniProt→colour, drawn on every chain.
        const spheres = [];
        const activeSph = this._activeSpheres || new Map();
        if (activeSph.size) {
            caAtoms.forEach(a => {
                if (stickKeys.has(`${a.chain || ''}|${a.resi}`)) return; // shown as a stick in focus mode
                const uni = toUni(a.chain, a.resi);
                if (uni != null && activeSph.has(uni)) spheres.push({ chain: a.chain || '', resi: a.resi, color: activeSph.get(uni) });
            });
        }

        // Ligands / cofactors (non-water HETATM), one entry per chain+resi+resn.
        const ligands = [];
        const seen = new Set();
        (model.selectedAtoms({ hetflag: true }) || []).forEach(a => {
            if (a.resn === 'HOH' || a.resn === 'WAT') return;
            const key = `${a.chain || ''}|${a.resi}|${a.resn}`;
            if (seen.has(key)) return;
            seen.add(key);
            ligands.push({ chain: a.chain || '', resi: a.resi, resn: a.resn, ion: this.ION_CODES.has(a.resn) });
        });

        return {
            format: this.currentFormat === 'mmcif' ? 'cif' : 'pdb',
            coordinates: this.currentPdbText,
            coloringMode: mode,
            cartoonBase: base,
            cartoonOpacity: this._inFocusMode ? 0.42 : 0.82,
            sphereRadius: 1.8,
            sphereOpacity: 0.92,
            cartoon, spheres, ligands, focus,
        };
    },

    /**
     * Tear down the currently-loaded model and per-structure caches without destroying the
     * 3Dmol viewer/WebGL context (recreating the context risks hitting the browser's live
     * context limit).  Called when the page navigates to a different protein so a stale model
     * can't be mistaken for an already-loaded structure on the next open.
     */
    clearModel() {
        this._selectedResi = null;
        this._inFocusMode = false;
        this._focusState = null;
        this._observedResi = null;
        this._observedResiByChain = null;
        this._lastRender = null;
        this._activeSpheres = null;
        this.currentStructure = null;
        this.currentPdbText = '';
        this._lastColoringContext = {};
        try {
            this.viewer?.removeAllModels();
            this.viewer?.removeAllLabels();
            this.viewer?.removeAllShapes();
            this.viewer?.render();
        } catch (_) { /* viewer may not be initialised yet */ }
    },

    // Zoom to a selection with a little breathing room, in a SINGLE smooth motion. zoomTo fits
    // the selection tightly; the same-duration zoom() blends in so it eases to ~12% wider without
    // the distracting "zoom in then pull back" of a delayed second animation.
    _zoomToWithMargin(sel) {
        this.viewer.zoomTo(sel, 450);
        this.viewer.zoom(0.95, 450); // small breathing room, a touch tighter than before
    },

    _proximityShapes: [],  // shape handles for PTM–variant dashed lines

    /**
     * Draw dashed lines from a PTM Cα to the Cα of each nearby variant residue.
     * geometry: from residueGeometry() — {uniPos, ca:{x,y,z}}
     * pairs: [{variantPos, dist, tier}] — what to draw lines to
     * tierColors: {1:'#...', 2:'#...', 3:'#...'}
     */
    showProximityLines(ptmPos, pairs, geometry) {
        this.clearProximityLines();
        if (!this.viewer || !pairs.length || !geometry.length) return;
        const caByUni = new Map();
        geometry.forEach(g => { if (g.uniPos != null && !caByUni.has(g.uniPos)) caByUni.set(g.uniPos, g.ca); });
        const ptmCa = caByUni.get(ptmPos);
        if (!ptmCa) return;
        const tierColors = { 1: '#ef5350', 2: '#ff7043', 3: '#ffa726' };
        pairs.forEach(({ variantPos, tier }) => {
            const ca = caByUni.get(variantPos);
            if (!ca) return;
            const color = tierColors[tier] || '#ffa726';
            const shape = this.viewer.addLine({
                start: { x: ptmCa.x, y: ptmCa.y, z: ptmCa.z },
                end: { x: ca.x, y: ca.y, z: ca.z },
                color, dashed: true, dashLength: 0.4, gapLength: 0.3, lineWidth: 2,
            });
            this._proximityShapes.push(shape);
        });
        this.viewer.render();
    },

    clearProximityLines() {
        if (!this._proximityShapes.length) return;
        this._proximityShapes.forEach(s => { try { this.viewer?.removeShape(s); } catch (_) {} });
        this._proximityShapes = [];
        this.viewer?.render(); // re-draw so the removed lines actually disappear
    },

    resetView(animated = false) {
        this._selectedResi = null;
        this._inFocusMode = false;
        this.clearProximityLines();
        this.viewer?.removeAllShapes();
        this.viewer?.removeAllLabels();
        this.viewer?.zoomTo({}, animated ? 600 : 0);
        this.viewer?.render();
    },

    screenshot() {
        const png = this.viewer?.pngURI();
        if (!png) return;
        const a = document.createElement('a');
        a.href = png; a.download = 'protein_3d.png'; a.click();
    },

    resize() { this.viewer?.resize(); }
};
