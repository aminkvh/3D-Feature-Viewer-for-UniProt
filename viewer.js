/* ============================================
   3Dmol.js Viewer Wrapper — Chrome Extension
   ============================================ */

const StructureViewer = {
    viewer: null,
    hoverCb: null,
    clickCb: null,
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

    init(container) {
        const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
        const bgColor = prefersDark ? '#0c111b' : '#f0f2f5';
        this.viewer = $3Dmol.createViewer(container, {
            backgroundColor: bgColor,
            antialias: true,
            cartoonQuality: 8,
        });
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
                // Re-load the already-fetched PDB text to rebuild the WebGL state.
                this.viewer.addModel(this.currentPdbText, this.currentFormat);
                this._buildObservedResiCache();
                this.applyCartoonColoring(this.activeColoringMode, this._lastColoringContext);
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
        // Validate URL — only allow HTTPS from trusted structure hosts.
        const ALLOWED_HOSTS = ['alphafold.ebi.ac.uk', 'www.ebi.ac.uk', 'files.rcsb.org'];
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
        this.currentFormat = pdbUrl.toLowerCase().endsWith('.cif') ? 'mmcif' : 'pdb';
        this.currentStructure = structure;

        this.viewer.removeAllModels();
        this.viewer.removeAllLabels();
        this.viewer.removeAllShapes();

        this.viewer.addModel(pdb, this.currentFormat);
        // Build SIFTS-validated residue index: find which PDB residues are actually modeled
        this._buildObservedResiCache();
        this.applyCartoonColoring('default');
        this.viewer.zoomTo();
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
        if (!structure || structure.source === 'AlphaFold' || !structure.mappedRanges?.length) return;
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
        if (!structure || structure.source === 'AlphaFold' || !structure.mappedRanges?.length) return { resi };
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
     * Applies an addStyle to the correct PDB residue(s) for a UniProt position.
     * Returns true if the style was applied, false if the position is outside the
     * mapped range or maps to an unmodelled (unobserved) residue.
     */
    allChainsAddStyle(resi, style, atomSel = {}) {
        const structure = this.currentStructure;
        if (!structure || structure.source === 'AlphaFold' || !structure.mappedRanges?.length) {
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
        const s = this.currentStructure;
        if (!s?.mappedRanges?.length) return null;
        const obs = this._observedResi; // Set of observed PDB resi with CA atoms
        const residues = [];
        s.mappedRanges.forEach(r => {
            for (let i = r.uniprotStart; i <= r.uniprotEnd; i++) {
                if (!obs) {
                    residues.push(i); // no cache yet: include all mapped
                } else {
                    const pdbResi = this._resiToPdb(i, r);
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
        this.viewer.removeAllShapes();
        this.viewer.removeAllLabels();
        const base = { opacity: 0.82, thickness: 0.25, ribbonWidth: 0.6 };
        this._applyModeStyles(mode, context, base);
        if (!defer) this.viewer.render();
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
            (context.hotspots || new Map()).forEach((tier, pos) => {
                const color = tierColors[tier] || '#b9c2cf';
                this.allChainsAddStyle(pos, { cartoon: { ...base, color } });
            });
        } else if (mode === 'distantContacts') {
            this.viewer.setStyle({}, { cartoon: { ...base, color: '#b9c2cf' } });
            const dcColors = { strong: '#6a1b9a', moderate: '#ab47bc' };
            (context.distantContacts || new Map()).forEach((tier, pos) => {
                const color = dcColors[tier] || '#b9c2cf';
                this.allChainsAddStyle(pos, { cartoon: { ...base, color } });
            });
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
        } else {
            this.viewer.setStyle({}, { cartoon: { ...base, color: '#00bcd4' } });
        }
    },

    /**
     * Render PTM residues as Cα VDW spheres.
     * For disulfide bonds, both begin and end positions are rendered.
     */
    showPTMs(ptms, ptmGroups) {
        this._selectedResi = null;
        const hoverMap = new Map();
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
            count++;

            // For disulfide bonds, also render the end position
            if (ptm.endPosition && ptm.endPosition !== ptm.position) {
                this.allChainsAddStyle(ptm.endPosition, { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }, { atom: 'CA' });
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

        this._bindHover(hoverMap, 'ptm');
        this.viewer.render();
        return count;
    },

    /**
     * PTM sphere refresh — re-applies base cartoon (clearing accumulated addStyle layers)
     * then re-renders only the currently visible PTM spheres.
     * Returns false when in focus mode so the caller can do a full applyMode() instead.
     */
    refreshPTMDisplay(ptms, ptmGroups) {
        if (!this.viewer || this._inFocusMode) return false;
        this.viewer.removeAllLabels();
        // Re-apply the base cartoon styles to clear any accumulated addStyle sphere layers
        // from previous showPTMs() calls — 3Dmol's removeAllShapes() only removes geometry
        // objects (addSphere/addCylinder), not addStyle layers, so they must be reset here.
        const base = { opacity: 0.82, thickness: 0.25, ribbonWidth: 0.6 };
        this._applyModeStyles(this.activeColoringMode, this._lastColoringContext || {}, base);
        return this.showPTMs(ptms, ptmGroups);
    },

    showVariants(filtered) {
        this._selectedResi = null;
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
        let placedPosCount = 0;
        let placedVarCount = 0;
        posMap.forEach((d, pos) => {
            const placed = this.allChainsAddStyle(pos, { sphere: { radius: 1.8, color: d.color, opacity: 0.92 } }, { atom: 'CA' });
            if (!placed) return;
            hoverMap.set(pos, d);
            placedPosCount++;
            placedVarCount += d.variants.length;
        });

        this._bindHover(hoverMap, 'variant');
        this.viewer.render();
        return { posCount: placedPosCount, varCount: placedVarCount };
    },

    _bindHover(map, mode) {
        const self = this;
        const reverseMap = this._reverseResidueMap();
        this.viewer.setHoverable({}, true,
            (atom, _v, ev) => {
                if (!atom?.resi) return;
                const uniResi = reverseMap.get(atom.resi) || atom.resi;
                const d = map.get(uniResi);
                if (d && self.hoverCb) self.hoverCb(d, mode, ev);
            },
            () => { if (self.hoverCb) self.hoverCb(null, mode, null); }
        );
        this.viewer.setClickable({}, true,
            (atom) => {
                if (!atom?.resi) return;
                const uniResi = reverseMap.get(atom.resi) || atom.resi;
                const d = map.get(uniResi);
                if (d && self.clickCb) self.clickCb(d, mode);
            }
        );
    },

    /**
     * Focus on a residue: keep cartoon visible, show selected residue
     * and its 5Å neighbors as ball-and-stick overlay.
     * No interaction lines. No residue label.
     */
    focusResidue(resi, annotations = {}) {
        if (!this.viewer) return;
        this.viewer.removeAllShapes();
        this.viewer.removeAllLabels();
        this._selectedResi = resi;

        const selector = this.residueSelector(resi);
        const pdbResi = selector.resi;
        const allAtoms = this.viewer.getModel().selectedAtoms(selector.chain ? { chain: selector.chain } : {});
        const selAtoms = allAtoms.filter(a => a.resi === pdbResi);
        if (selAtoms.length === 0) return;

        const CONTACT_DIST = 5.0;
        const nearbyResis = new Set();
        nearbyResis.add(pdbResi);

        // Find nearby residues: any atom within 5Å of any atom in selected residue
        allAtoms.forEach(a => {
            if (a.resi === pdbResi) return;
            for (const sa of selAtoms) {
                const dx = a.x - sa.x, dy = a.y - sa.y, dz = a.z - sa.z;
                if (Math.sqrt(dx*dx + dy*dy + dz*dz) < CONTACT_DIST) {
                    nearbyResis.add(a.resi);
                    break;
                }
            }
        });

        // Dim the cartoon to background opacity only on first focus (avoids expensive
        // full geometry rebuild on every residue click while already in focus mode)
        if (!this._inFocusMode) {
            this._inFocusMode = true;
            const focusBase = { opacity: 0.42, thickness: 0.2, ribbonWidth: 0.5 };
            this._applyModeStyles(this.activeColoringMode, this._lastColoringContext, focusBase);
        }

        // Show nearby residues as thin sticks
        const reverseMap = this._reverseResidueMap();
        const annotated = annotations.annotatedResidues || new Map();
        nearbyResis.forEach(nResi => {
            if (nResi === pdbResi) return;
            const uniResi = reverseMap.get(nResi) || nResi;
            const color = annotated.get(uniResi)?.color;
            this.viewer.addStyle(
                selector.chain ? { chain: selector.chain, resi: nResi } : { resi: nResi },
                { stick: color ? { radius: 0.15, color, opacity: 0.8 } : { radius: 0.12, colorscheme: 'Jmol', opacity: 0.6 } }
            );
        });

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
        nearbyResis.forEach(nResi => {
            const uni = reverseMap.get(nResi) || nResi;
            focusHoverMap.set(uni, { position: uni, color: annotated.get(uni)?.color });
        });
        this._bindHover(focusHoverMap, 'focus');

        // Zoom to the neighborhood with animation and render
        this.viewer.zoomTo(selector.chain ? { chain: selector.chain, resi: Array.from(nearbyResis) } : { resi: Array.from(nearbyResis) }, 600);
        this.viewer.render();

        return new Set(Array.from(nearbyResis).map(p => reverseMap.get(p) || p));
    },

    _reverseResidueMap() {
        const out = new Map();
        const s = this.currentStructure;
        const obs = this._observedResi; // only map residues that are physically present
        (s?.mappedRanges || []).forEach(r => {
            if (s.seqresToAuthor && r.seqresStart != null) {
                // Non-linear mapping (e.g. chimeric structures like 6CDU):
                // invert seqresToAuthor over this range to get correct UniProt positions.
                const seqresEnd = r.seqresStart + (r.uniprotEnd - r.uniprotStart);
                s.seqresToAuthor.forEach((author, seqres) => {
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

    resetView(animated = false) {
        this._selectedResi = null;
        this._inFocusMode = false;
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
