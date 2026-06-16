const UFVExport = (() => {
    'use strict';

    // Map a UniProt position to the author PDB residue number for the given structure.
    // Returns the original position for AlphaFold (1:1 mapping) or when no mapping exists.
    function uniprotToPdbResi(pos, structure) {
        if (!structure || structure.source === 'AlphaFold' || !structure.mappedRanges?.length) return pos;
        for (const r of structure.mappedRanges) {
            if (pos >= r.uniprotStart && pos <= r.uniprotEnd) {
                if (structure.seqresToAuthor && r.seqresStart != null) {
                    const seqresPos = r.seqresStart + (pos - r.uniprotStart);
                    return structure.seqresToAuthor.get(seqresPos) ?? null;
                }
                return r.pdbStart + (pos - r.uniprotStart);
            }
        }
        return null;
    }

    function formatSelection(positions, format = 'pymol', structure = null) {
        const mapped = structure
            ? positions.map(p => uniprotToPdbResi(p, structure)).filter(p => p != null)
            : positions;
        const values = Array.from(new Set(mapped.map(Number).filter(Boolean))).sort((a, b) => a - b);
        if (!values.length) return '';
        if (format === 'vmd') return `resid ${values.join(' ')}`;
        return `resi ${values.join('+')}`;
    }

    function setBFactor(line, value) {
        if (!/^(ATOM  |HETATM)/.test(line)) return line;
        const beta = value.toFixed(2).padStart(6, ' ');
        if (line.length < 66) return line.padEnd(60, ' ') + beta;
        return line.slice(0, 60) + beta + line.slice(66);
    }

    function rewritePdbBeta(pdbText, positions, structure, colorMode, colorContext) {
        const selected = new Set(positions.map(Number));
        const ranges = structure?.mappedRanges || [];
        const chainId = structure?.chainId || null;
        function mapsToUniProt(pdbResi) {
            if (!ranges.length) return pdbResi;
            for (const r of ranges) {
                if (chainId && r.chainId && r.chainId !== chainId) continue;
                if (pdbResi >= r.pdbStart && pdbResi <= r.pdbEnd) {
                    return r.uniprotStart + (pdbResi - r.pdbStart);
                }
            }
            return null;
        }
        return pdbText.split(/\r?\n/).map(line => {
            if (!/^(ATOM  |HETATM)/.test(line)) return line;
            const atomChain = line.slice(21, 22).trim();
            if (chainId && atomChain && atomChain !== chainId) return setBFactor(line, 0);
            const pdbResi = parseInt(line.slice(22, 26), 10);
            const uniprotResi = mapsToUniProt(pdbResi);
            if (!uniprotResi) return setBFactor(line, 0);
            // pLDDT and experimental B-factor: keep original file values unchanged
            if (colorMode === 'plddt' || colorMode === 'bfactor') return line;
            // AlphaMissense: use per-residue average score (0–1)
            if (colorMode === 'alphaMissense') {
                const d = colorContext?.get(uniprotResi);
                return setBFactor(line, d != null ? d.avg : 0);
            }
            // Hotspots: encode tier as fractional B-factor (strong=1.0, moderate=0.67, weak=0.33, none=0)
            if (colorMode === 'hotspots') {
                const tierVal = { strong: 1.0, moderate: 0.67, weak: 0.33 };
                const tier = colorContext instanceof Map ? colorContext.get(uniprotResi) : null;
                return setBFactor(line, tierVal[tier] ?? 0);
            }
            // Distant contacts: encode tier (strong=1.0, moderate=0.5)
            if (colorMode === 'distantContacts') {
                const tierVal = { strong: 1.0, moderate: 0.5 };
                const tier = colorContext instanceof Map ? colorContext.get(uniprotResi) : null;
                return setBFactor(line, tierVal[tier] ?? 0);
            }
            // Residue burden: 1 if flagged, 0 otherwise
            if (colorMode === 'residueBurden') {
                return setBFactor(line, colorContext instanceof Set && colorContext.has(uniprotResi) ? 1 : 0);
            }
            // Constraint pocket: encode geometric class as fractional B-factor
            // (buried pocket = 1.0, exposed = 0.66, none = 0).
            if (colorMode === 'prism') {
                const catVal = { pocket: 1.0, exposed: 0.66 };
                const info = colorContext instanceof Map ? colorContext.get(uniprotResi) : null;
                return setBFactor(line, info ? (catVal[info.cat] ?? 1.0) : 0);
            }
            // Default (cyan): 1 for currently displayed (selected) residues, 0 elsewhere
            return setBFactor(line, selected.has(uniprotResi) ? 1 : 0);
        }).join('\n');
    }

    function downloadText(filename, text, type = 'chemical/x-pdb') {
        const blob = new Blob([text], { type });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        a.click();
        setTimeout(() => URL.revokeObjectURL(url), 500);
    }

    async function copyText(text) {
        if (!text) return false;
        try {
            await navigator.clipboard.writeText(text);
            return true;
        } catch (_) {
            const ta = document.createElement('textarea');
            ta.value = text;
            ta.style.cssText = 'position:fixed;opacity:0';
            document.body.appendChild(ta);
            ta.select();
            document.execCommand('copy');
            document.body.removeChild(ta);
            return true;
        }
    }

    /**
     * Build a CSV annotation matrix: one row per residue, one-hot columns for
     * each PTM category and each disease, plus the per-residue AlphaMissense
     * average score when the CSV map is available.
     */
    // Resolve UniProt position → author residue number for a SPECIFIC chain of a structure.
    function uniprotToPdbResiForChain(pos, structure, chain) {
        if (!structure || structure.source === 'AlphaFold') return pos;
        const ranges = (structure.chainMappings?.[chain]) || structure.mappedRanges || [];
        const map = (structure.chainSeqresToAuthor?.[chain]) || structure.seqresToAuthor || null;
        if (!ranges.length) return pos;
        for (const r of ranges) {
            if (pos >= r.uniprotStart && pos <= r.uniprotEnd) {
                if (map && r.seqresStart != null) {
                    const seqresPos = r.seqresStart + (pos - r.uniprotStart);
                    return map.get(seqresPos) ?? null;
                }
                return r.pdbStart + (pos - r.uniprotStart);
            }
        }
        return null;
    }

    function buildResidueMatrix(sequence, ptms, ptmGroups, variants, amMap, analysis = {}, structure = null, sites = [], mutagenesis = []) {
        const AM_AAS = 'ACDEFGHIKLMNPQRSTVWY';
        const ptmCats = Object.keys(ptmGroups).sort();
        const diseaseSet = new Set();
        variants.forEach(v => (v.diseasePairs || []).forEach(p => { if (p.label) diseaseSet.add(p.label); }));
        const diseases = [...diseaseSet].sort();

        // Structure-dependent columns (pdb_residue, hotspot_tier, contact_hub_tier) are repeated
        // per chain for multi-chain structures because each subunit resolves and scores residues
        // slightly differently.  Single-chain output keeps the original column layout unchanged.
        const chains = structure?.chainIds?.length > 1 ? structure.chainIds : null;
        const hotspotTierNum = { strong: 3, moderate: 2, weak: 1 };
        const hubTierNum = { strong: 2, moderate: 1 };
        const tierFor = (map, pos, lut) => (map instanceof Map ? (lut[map.get(pos)] ?? 0) : 0);

        // Constraint-pocket values (UniProt-position-keyed, structure-dependent): per-residue
        // BH-FDR q-value and geometric class for candidates (positions with a positive residual).
        const pocketByPos = analysis.prism?.byPos instanceof Map ? analysis.prism.byPos : null;
        const pocketCols = pos => {
            const info = pocketByPos?.get(pos);
            return [info ? info.q.toFixed(4) : '', info ? info.cat : ''];
        };

        // PTM–Variant Proximity columns (structure-dependent, PTM positions only).
        // For residue-centric CSV, nearby variants are serialised as a semicolon-delimited list.
        // Per-residue ligand contacts (CCD codes within 5 Å), structure-dependent.
        const ligandContacts = analysis.ligandContacts instanceof Map ? analysis.ligandContacts : null;
        const ligandCol = pos => {
            const set = ligandContacts?.get(pos);
            return set && set.size ? [...set].sort().join(';') : '';
        };

        const proxMap = analysis.ptmVariantProximity instanceof Map ? analysis.ptmVariantProximity : null;
        const proxCols = pos => {
            const p = proxMap?.get(pos);
            if (!p) return ['', '', '', '', ''];
            const allVars = [
                ...p.tier1.map(v => `${v.wildType||''}${v.position}${v.mutant||''}(T1)`),
                ...p.tier2.map(({variant:v, dist}) => `${v.wildType||''}${v.position}${v.mutant||''}(T2,${dist.toFixed(1)}Å)`),
                ...p.tier3.map(({variant:v, dist}) => `${v.wildType||''}${v.position}${v.mutant||''}(T3,${dist.toFixed(1)}Å)`),
            ];
            return [
                p.tier ?? '',
                p.nearestVariant ?? '',
                p.nearestDist !== null ? p.nearestDist.toFixed(1) : '',
                p.nearbyCount8A,
                p.pathCount8A,
            ];
        };

        // Build per-position lookups
        const ptmByPos = new Map(); // pos → Set<category>
        ptms.forEach(p => {
            [p.position, p.endPosition].filter(Boolean).forEach(pos => {
                if (!ptmByPos.has(pos)) ptmByPos.set(pos, new Set());
                ptmByPos.get(pos).add(p.category);
            });
        });
        const diseaseByPos = new Map(); // pos → Set<label>
        variants.forEach(v => {
            if (!diseaseByPos.has(v.position)) diseaseByPos.set(v.position, new Set());
            (v.diseasePairs || []).forEach(p => { if (p.label) diseaseByPos.get(v.position).add(p.label); });
        });

        // Sites (one-hot by category), mutagenesis flag, gnomAD AF and variant counts — per position.
        const siteCat = d => { const t = (d || '').toLowerCase(); return t.includes('active') ? 'active' : t.includes('metal') ? 'metal' : t.includes('binding') ? 'binding' : 'other'; };
        const siteByPos = new Map(); // pos → Set<category>
        (sites || []).forEach(x => {
            for (let p = x.position; p <= (x.endPosition || x.position); p++) {
                if (!siteByPos.has(p)) siteByPos.set(p, new Set());
                siteByPos.get(p).add(siteCat(x.description));
            }
        });
        const mutByPos = new Set();
        (mutagenesis || []).forEach(mtg => { for (let p = mtg.position; p <= (mtg.endPosition || mtg.position); p++) mutByPos.add(p); });
        const nVarByPos = new Map(), nPathByPos = new Map(), gnomadByPos = new Map();
        variants.forEach(v => {
            nVarByPos.set(v.position, (nVarByPos.get(v.position) || 0) + 1);
            const path = `${v.consequence || ''} ${v.clinVarSignificance || ''}`.toLowerCase().includes('pathogenic');
            if (path) nPathByPos.set(v.position, (nPathByPos.get(v.position) || 0) + 1);
            if (Number.isFinite(v.gnomadAf)) gnomadByPos.set(v.position, Math.max(gnomadByPos.get(v.position) ?? 0, v.gnomadAf));
        });
        const featCols = pos => {
            const sc = siteByPos.get(pos);
            const af = gnomadByPos.get(pos);
            return [
                sc ? 1 : 0, sc?.has('active') ? 1 : 0, sc?.has('binding') ? 1 : 0, sc?.has('metal') ? 1 : 0,
                mutByPos.has(pos) ? 1 : 0,
                nVarByPos.get(pos) || 0, nPathByPos.get(pos) || 0,
                af != null ? af.toExponential(3) : '',
            ];
        };

        // Sanitise column names for CSV
        const safe = s => s.replace(/[^A-Za-z0-9_]/g, '_').replace(/_+/g, '_');
        // Single-chain: original layout (pdb_residue before the one-hot columns; tiers after
        // am_avg_score).  Multi-chain: structure-dependent columns repeated per chain at the end.
        const proxHeaders = ['ptm_variant_tier', 'nearest_variant_to_ptm', 'ptm_variant_distance_A', 'nearby_variant_count_8A', 'nearby_pathogenic_count_8A'];
        const featHeaders = ['site', 'site_active', 'site_binding', 'site_metal', 'mutagenesis', 'n_variants', 'n_pathogenic', 'gnomad_max_af'];
        const headers = chains
            ? [
                'position', 'aa',
                ...ptmCats.map(c => `ptm_${safe(c)}`),
                ...diseases.map(d => `disease_${safe(d)}`),
                ...featHeaders,
                'am_avg_score',
                'am_max_score',
                'am_n_subs',
                'residue_burden',
                'constraint_pocket_q',
                'constraint_pocket_class',
                ...proxHeaders,
                'nearby_ligands',
                ...chains.flatMap(c => [`pdb_residue_${c}`, `hotspot_tier_${c}`, `contact_hub_tier_${c}`]),
            ]
            : [
                'position', 'aa', 'pdb_residue',
                ...ptmCats.map(c => `ptm_${safe(c)}`),
                ...diseases.map(d => `disease_${safe(d)}`),
                ...featHeaders,
                'am_avg_score',
                'am_max_score',
                'am_n_subs',
                'hotspot_tier',
                'contact_hub_tier',
                'residue_burden',
                'constraint_pocket_q',
                'constraint_pocket_class',
                ...proxHeaders,
                'nearby_ligands',
            ];
        const rows = [headers.join(',')];
        sequence.split('').forEach((aa, i) => {
            const pos = i + 1;
            const ptmFlags = ptmCats.map(c => ptmByPos.get(pos)?.has(c) ? 1 : 0);
            const diseaseFlags = diseases.map(d => diseaseByPos.get(pos)?.has(d) ? 1 : 0);
            let amAvg = '', amMax = '', amN = '';
            if (amMap && amMap.size > 0) {
                const scores = [];
                for (const mut of AM_AAS) {
                    if (mut === aa) continue;
                    const sc = amMap.get(`${aa}${pos}${mut}`);
                    if (Number.isFinite(sc)) scores.push(sc);
                }
                if (scores.length > 0) {
                    amAvg = (scores.reduce((a, b) => a + b, 0) / scores.length).toFixed(4);
                    amMax = Math.max(...scores).toFixed(4);
                    amN = scores.length;
                }
            }
            const burden = analysis.residueBurden instanceof Set && analysis.residueBurden.has(pos) ? 1 : 0;
            const [pocketQ, pocketClass] = pocketCols(pos);
            const proxVals = proxCols(pos);
            if (chains) {
                const structVals = chains.flatMap(c => [
                    uniprotToPdbResiForChain(pos, structure, c) ?? '',
                    tierFor(analysis.hotspotsByChain?.get(c), pos, hotspotTierNum),
                    tierFor(analysis.distantContactsByChain?.get(c), pos, hubTierNum),
                ]);
                rows.push([pos, aa, ...ptmFlags, ...diseaseFlags, ...featCols(pos), amAvg, amMax, amN, burden, pocketQ, pocketClass, ...proxVals, ligandCol(pos), ...structVals].join(','));
            } else {
                const pdbResi = uniprotToPdbResi(pos, structure) ?? '';
                const hotspotTier = tierFor(analysis.hotspots, pos, hotspotTierNum);
                const hubTier = tierFor(analysis.distantContacts, pos, hubTierNum);
                rows.push([pos, aa, pdbResi, ...ptmFlags, ...diseaseFlags, ...featCols(pos), amAvg, amMax, amN, hotspotTier, hubTier, burden, pocketQ, pocketClass, ...proxVals, ligandCol(pos)].join(','));
            }
        });
        return rows.join('\n');
    }

    // ── Session export (PyMOL / VMD) ─────────────────────────────────────────────────────────────
    // Parse '#rrggbb' or 'rgb(r,g,b)' into [r,g,b] floats in 0–1 (the form PyMOL/VMD want).
    function colorToUnit(c) {
        if (!c) return [0.5, 0.5, 0.5];
        const m = /^rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)$/i.exec(c.trim());
        if (m) return [(+m[1]) / 255, (+m[2]) / 255, (+m[3]) / 255];
        const h = c.replace('#', '');
        const n = parseInt(h.length === 3 ? h.split('').map(x => x + x).join('') : h, 16);
        return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255];
    }
    const fix = n => n.toFixed(3);

    // [{chain,resi}] → Map<chain, resi[]>.
    function itemsByChain(items) {
        const m = new Map();
        items.forEach(({ chain, resi }) => { if (!m.has(chain)) m.set(chain, []); m.get(chain).push(resi); });
        return m;
    }

    // Group [{chain,resi,color}] by colour, then by chain → sorted resi list.
    function groupByColor(items) {
        const byColor = new Map(); // color → Map<chain, number[]>
        items.forEach(({ chain, resi, color }) => {
            if (!byColor.has(color)) byColor.set(color, new Map());
            const byChain = byColor.get(color);
            if (!byChain.has(chain)) byChain.set(chain, []);
            byChain.get(chain).push(resi);
        });
        return byColor;
    }

    /**
     * Build a self-contained PyMOL script (.pml) that reproduces the current view: it embeds the
     * structure coordinates, loads them, then re-applies the cartoon colours, annotation spheres,
     * and ligand sticks exactly as shown. Camera is auto-framed (orient) — orientation may differ.
     */
    function buildPymolSession(scene, objName = 'structure') {
        if (!scene) return '';
        const L = [];
        L.push('# 3D Feature Viewer for UniProt — PyMOL session');
        L.push(`# Coloring mode: ${scene.coloringMode}`);
        L.push('from pymol import cmd, util');
        L.push('python');
        L.push('import tempfile, os');
        L.push('_data = r"""');
        L.push(scene.coordinates.replace(/\r\n/g, '\n').replace(/\s+$/,''));
        L.push('"""');
        L.push(`_fd, _p = tempfile.mkstemp(suffix=".${scene.format}")`);
        L.push('with os.fdopen(_fd, "w") as _fh: _fh.write(_data)');
        L.push(`cmd.load(_p, "${objName}")`);
        L.push('os.remove(_p)');
        L.push('python end');
        L.push('');
        L.push('hide everything');
        L.push('bg_color white');
        L.push('show cartoon, polymer');
        L.push(`set cartoon_transparency, ${fix(1 - scene.cartoonOpacity)}`);
        // Base cartoon colour.
        const [br, bg, bb] = colorToUnit(scene.cartoonBase);
        L.push(`set_color ufv_base, [${fix(br)}, ${fix(bg)}, ${fix(bb)}]`);
        L.push('color ufv_base, polymer');
        // Per-residue cartoon overrides.
        let ci = 0;
        groupByColor(scene.cartoon).forEach((byChain, color) => {
            const name = `ufv_c${ci++}`;
            const [r, g, b] = colorToUnit(color);
            L.push(`set_color ${name}, [${fix(r)}, ${fix(g)}, ${fix(b)}]`);
            L.push(`color ${name}, (${pymolSel(byChain)})`);
        });
        // Annotation Cα spheres.
        if (scene.spheres.length) {
            const sphereSels = [];
            let si = 0;
            groupByColor(scene.spheres).forEach((byChain, color) => {
                const name = `ufv_s${si++}`;
                const [r, g, b] = colorToUnit(color);
                const sel = `(${pymolSel(byChain)}) and name CA`;
                sphereSels.push(sel);
                L.push(`set_color ${name}, [${fix(r)}, ${fix(g)}, ${fix(b)}]`);
                L.push(`show spheres, ${sel}`);
                L.push(`color ${name}, ${sel}`);
            });
            L.push(`select ufv_spheres, (${sphereSels.join(') or (')})`);
            L.push(`alter ufv_spheres, vdw=${fix(scene.sphereRadius)}`);
            L.push('set sphere_scale, 1.0, ufv_spheres');
            L.push(`set sphere_transparency, ${fix(1 - scene.sphereOpacity)}, ufv_spheres`);
            L.push('rebuild');
            L.push('deselect');
        }
        // Ligands / cofactors: sticks coloured by element; ions as small spheres.
        if (scene.ligands.length) {
            L.push('show sticks, (not polymer and not solvent)');
            L.push('util.cnc (not polymer and not solvent)');
            const ions = scene.ligands.filter(l => l.ion);
            if (ions.length) {
                L.push(`show spheres, (${pymolSel(itemsByChain(ions))})`);
            }
        }
        // Zoom-in (focus) state: selected residue + 5 Å neighbourhood as sticks.
        if (scene.focus) {
            const f = scene.focus;
            const allItems = [{ chain: f.selChain, resi: f.pdbResi }, ...f.sticks];
            L.push('# Zoom-in: focused residue + 5 A neighbourhood (ball-and-stick)');
            L.push(`show sticks, (${pymolSel(itemsByChain(allItems))})`);
            L.push(`select ufv_focus_sel, (${pymolSel(itemsByChain([{ chain: f.selChain, resi: f.pdbResi }]))})`);
            L.push('show spheres, ufv_focus_sel');
            L.push(`alter ufv_focus_sel, vdw=${fix(0.4)}`);
            // Element (Jmol-style) colouring for the selected residue + uncoloured/ligand neighbours.
            const jmolItems = [{ chain: f.selChain, resi: f.pdbResi }, ...f.sticks.filter(x => !x.color)];
            const jmolSel = pymolSel(itemsByChain(jmolItems));
            L.push(`color grey80, (${jmolSel})`);
            L.push(`util.cnc (${jmolSel})`);
            // Annotation-coloured neighbours: whole stick painted the annotation colour.
            let fi = 0;
            groupByColor(f.sticks.filter(x => x.color)).forEach((byChain, color) => {
                const name = `ufv_f${fi++}`;
                const [r, g, b] = colorToUnit(color);
                L.push(`set_color ${name}, [${fix(r)}, ${fix(g)}, ${fix(b)}]`);
                L.push(`color ${name}, (${pymolSel(byChain)})`);
            });
            L.push('rebuild');
            L.push('deselect');
        }
        // Frame the view: if zoom-in is active, zoom to the selected residue; otherwise orient to fit.
        if (scene.focus) {
            const selSel = pymolSel(itemsByChain([{ chain: scene.focus.selChain, resi: scene.focus.pdbResi }]));
            L.push(`zoom (${selSel})`);
        } else {
            L.push('orient');
        }
        return L.join('\n') + '\n';
    }

    // PyMOL selection from Map<chain, resi[]>: '(chain A and resi 1+2+3) or (resi 9+10)'.
    function pymolSel(byChain) {
        const parts = [];
        byChain.forEach((resis, chain) => {
            const list = [...new Set(resis)].sort((a, b) => a - b).join('+');
            parts.push(chain ? `chain ${chain} and resi ${list}` : `resi ${list}`);
        });
        return parts.join(') or (');
    }

    // VMD resid list (space-separated) per chain.
    function vmdSelParts(byChain) {
        const parts = [];
        byChain.forEach((resis, chain) => {
            const list = [...new Set(resis)].sort((a, b) => a - b).join(' ');
            parts.push(chain ? `(chain ${chain} and resid ${list})` : `(resid ${list})`);
        });
        return parts.join(' or ');
    }

    /**
     * Build a self-contained VMD script (.tcl/.vmd): writes the embedded coordinates to a file in
     * the current directory, loads it, then adds one representation per colour group (cartoon,
     * VDW spheres, licorice ligands). Custom colours are registered from ColorID 33 upward.
     */
    function buildVmdSession(scene, objName = 'structure') {
        if (!scene) return '';
        const L = [];
        L.push('# 3D Feature Viewer for UniProt — VMD session');
        L.push(`# Coloring mode: ${scene.coloringMode}`);
        L.push(`set _data {${scene.coordinates.replace(/\r\n/g, '\n')}}`);
        L.push(`set _p [file join [pwd] "ufv_${objName}.${scene.format}"]`);
        L.push('set _fh [open $_p w]');
        L.push('puts -nonewline $_fh $_data');
        L.push('close $_fh');
        L.push('mol new $_p waitfor all');
        L.push('mol delrep 0 top');
        let cid = 33; // first user-definable ColorID
        const defColor = c => { const [r, g, b] = colorToUnit(c); L.push(`color change rgb ${cid} ${fix(r)} ${fix(g)} ${fix(b)}`); return cid++; };
        // Base cartoon.
        const baseId = defColor(scene.cartoonBase);
        L.push('mol representation NewCartoon');
        L.push(`mol color ColorID ${baseId}`);
        L.push('mol selection {protein}');
        L.push('mol material Opaque');
        L.push('mol addrep top');
        // Cartoon overrides.
        groupByColor(scene.cartoon).forEach((byChain, color) => {
            const id = defColor(color);
            L.push('mol representation NewCartoon');
            L.push(`mol color ColorID ${id}`);
            L.push(`mol selection {protein and (${vmdSelParts(byChain)})}`);
            L.push('mol addrep top');
        });
        // Annotation spheres (VDW on CA).
        groupByColor(scene.spheres).forEach((byChain, color) => {
            const id = defColor(color);
            L.push('mol representation VDW 1.0 16');
            L.push(`mol color ColorID ${id}`);
            L.push(`mol selection {name CA and (${vmdSelParts(byChain)})}`);
            L.push('mol addrep top');
        });
        // Ligands: licorice, coloured by element name.
        if (scene.ligands.length) {
            L.push('mol representation Licorice 0.3 12');
            L.push('mol color Name');
            L.push('mol selection {not protein and not water}');
            L.push('mol addrep top');
        }
        // Zoom-in (focus) state: selected residue + 5 Å neighbourhood as licorice.
        if (scene.focus) {
            const f = scene.focus;
            // Element-coloured licorice: selected residue + uncoloured/ligand neighbours.
            const jmolItems = [{ chain: f.selChain, resi: f.pdbResi }, ...f.sticks.filter(x => !x.color)];
            L.push('mol representation Licorice 0.2 12');
            L.push('mol color Name');
            L.push(`mol selection {${vmdSelParts(itemsByChain(jmolItems))}}`);
            L.push('mol addrep top');
            // Annotation-coloured neighbours.
            groupByColor(f.sticks.filter(x => x.color)).forEach((byChain, color) => {
                const id = defColor(color);
                L.push('mol representation Licorice 0.2 12');
                L.push(`mol color ColorID ${id}`);
                L.push(`mol selection {${vmdSelParts(byChain)}}`);
                L.push('mol addrep top');
            });
            // Selected residue: small spheres.
            L.push('mol representation VDW 0.4 16');
            L.push('mol color Name');
            L.push(`mol selection {${vmdSelParts(itemsByChain([{ chain: f.selChain, resi: f.pdbResi }]))}}`);
            L.push('mol addrep top');
            // Zoom to the selected residue pocket
            L.push(`molinfo top set center_matrix [molinfo top get center_matrix]`);
            L.push(`set sel [atomselect top {${vmdSelParts(itemsByChain([{ chain: f.selChain, resi: f.pdbResi }]))}}]`);
            L.push('set center [measure center $sel]');
            L.push('set matrix [trans center $center]');
            L.push('molinfo top set center_matrix $matrix');
            L.push('$sel delete');
            L.push('display resetview');
        } else {
            L.push('display resetview');
        }
        return L.join('\n') + '\n';
    }

    return { formatSelection, rewritePdbBeta, buildResidueMatrix, buildPymolSession, buildVmdSession, downloadText, copyText };
})();
