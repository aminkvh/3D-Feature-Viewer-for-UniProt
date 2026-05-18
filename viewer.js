/* ============================================
   3Dmol.js Viewer Wrapper — Chrome Extension
   ============================================ */

const StructureViewer = {
    viewer: null,
    hoverCb: null,
    clickCb: null,
    dblClickCb: null,
    _selectedResi: null,

    init(container) {
        const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
        const bgColor = prefersDark ? '#0c111b' : '#f0f2f5';
        this.viewer = $3Dmol.createViewer(container, {
            backgroundColor: bgColor,
            antialias: true,
            cartoonQuality: 8,
        });
        window.addEventListener('resize', () => { if (this.viewer) this.viewer.resize(); });

        // Double-click on background = reset
        container.addEventListener('dblclick', () => {
            if (this._selectedResi !== null && this.dblClickCb) {
                this.dblClickCb();
            }
        });
    },

    async loadStructure(pdbUrl) {
        // Validate URL — only allow HTTPS from trusted AlphaFold host
        const ALLOWED_HOSTS = ['alphafold.ebi.ac.uk'];
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

        this.viewer.removeAllModels();
        this.viewer.removeAllLabels();
        this.viewer.removeAllShapes();

        this.viewer.addModel(pdb, 'pdb');
        this.viewer.setStyle({}, {
            cartoon: { color: '#00bcd4', opacity: 0.88, thickness: 0.25, ribbonWidth: 0.6 }
        });
        this.viewer.zoomTo();
        this.viewer.render();
    },

    /**
     * Render PTM residues as Cα VDW spheres.
     * For disulfide bonds, both begin and end positions are rendered.
     */
    showPTMs(ptms, ptmGroups) {
        this._selectedResi = null;
        this.viewer.setStyle({}, {
            cartoon: { color: '#00bcd4', opacity: 0.82, thickness: 0.25, ribbonWidth: 0.6 }
        });
        this.viewer.removeAllLabels();
        this.viewer.removeAllShapes();

        const hoverMap = new Map();
        let count = 0;

        ptms.forEach(ptm => {
            const g = ptmGroups[ptm.category];
            if (!g || !g.visible) return;

            // Render begin position
            this.viewer.addStyle(
                { resi: ptm.position, atom: 'CA' },
                { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }
            );
            hoverMap.set(ptm.position, ptm);

            // For disulfide bonds, also render the end position
            if (ptm.endPosition && ptm.endPosition !== ptm.position) {
                this.viewer.addStyle(
                    { resi: ptm.endPosition, atom: 'CA' },
                    { sphere: { radius: 1.8, color: ptm.color, opacity: 0.92 } }
                );
                // Create a hover entry for the end position too
                if (!hoverMap.has(ptm.endPosition)) {
                    hoverMap.set(ptm.endPosition, {
                        ...ptm,
                        position: ptm.endPosition,
                        description: ptm.description + ` (bonded to ${ptm.position})`
                    });
                }
            }
            count++;
        });

        this._bindHover(hoverMap, 'ptm');
        this.viewer.render();
        return count;
    },

    showVariants(filtered) {
        this._selectedResi = null;
        this.viewer.setStyle({}, {
            cartoon: { color: '#00bcd4', opacity: 0.82, thickness: 0.25, ribbonWidth: 0.6 }
        });
        this.viewer.removeAllLabels();
        this.viewer.removeAllShapes();

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
        posMap.forEach((d, pos) => {
            this.viewer.addStyle(
                { resi: pos, atom: 'CA' },
                { sphere: { radius: 1.8, color: d.color, opacity: 0.92 } }
            );
            hoverMap.set(pos, d);
        });

        this._bindHover(hoverMap, 'variant');
        this.viewer.render();
        return { posCount: posMap.size, varCount: filtered.length };
    },

    _bindHover(map, mode) {
        const self = this;
        this.viewer.setHoverable({}, true,
            (atom, _v, ev) => {
                if (!atom?.resi) return;
                const d = map.get(atom.resi);
                if (d && self.hoverCb) self.hoverCb(d, mode, ev);
            },
            () => { if (self.hoverCb) self.hoverCb(null, mode, null); }
        );
        this.viewer.setClickable({}, true,
            (atom) => {
                if (!atom?.resi) return;
                const d = map.get(atom.resi);
                if (d && self.clickCb) self.clickCb(d, mode);
            }
        );
    },

    /**
     * Focus on a residue: keep cartoon visible, show selected residue
     * and its 5Å neighbors as ball-and-stick overlay.
     * No interaction lines. No residue label.
     */
    focusResidue(resi) {
        if (!this.viewer) return;
        this.viewer.removeAllShapes();
        this.viewer.removeAllLabels();
        this._selectedResi = resi;

        const allAtoms = this.viewer.getModel().selectedAtoms({});
        const selAtoms = allAtoms.filter(a => a.resi === resi);
        if (selAtoms.length === 0) return;

        const CONTACT_DIST = 5.0;
        const nearbyResis = new Set();
        nearbyResis.add(resi);

        // Find nearby residues: any atom within 5Å of any atom in selected residue
        allAtoms.forEach(a => {
            if (a.resi === resi) return;
            for (const sa of selAtoms) {
                const dx = a.x - sa.x, dy = a.y - sa.y, dz = a.z - sa.z;
                if (Math.sqrt(dx*dx + dy*dy + dz*dz) < CONTACT_DIST) {
                    nearbyResis.add(a.resi);
                    break;
                }
            }
        });

        // Keep the full cartoon backbone visible
        this.viewer.setStyle({}, {
            cartoon: { color: '#00bcd4', opacity: 0.55, thickness: 0.2, ribbonWidth: 0.5 }
        });

        // Show nearby residues as thin sticks
        nearbyResis.forEach(nResi => {
            if (nResi === resi) return;
            this.viewer.addStyle(
                { resi: nResi },
                { stick: { radius: 0.12, colorscheme: 'Jmol', opacity: 0.6 } }
            );
        });

        // Show selected residue as ball-and-stick (prominent)
        this.viewer.addStyle(
            { resi: resi },
            {
                stick: { radius: 0.2, colorscheme: 'Jmol' },
                sphere: { radius: 0.4, colorscheme: 'Jmol' }
            }
        );

        // Zoom to the neighborhood
        this.viewer.zoomTo({ resi: Array.from(nearbyResis) });
        this.viewer.render();
    },

    resetView() {
        this._selectedResi = null;
        this.viewer?.removeAllShapes();
        this.viewer?.removeAllLabels();
        this.viewer?.zoomTo();
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
