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

    function buildResidueMatrix(sequence, ptms, ptmGroups, variants, amMap, analysis = {}, structure = null) {
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

        // Sanitise column names for CSV
        const safe = s => s.replace(/[^A-Za-z0-9_]/g, '_').replace(/_+/g, '_');
        // Single-chain: original layout (pdb_residue before the one-hot columns; tiers after
        // am_avg_score).  Multi-chain: structure-dependent columns repeated per chain at the end.
        const headers = chains
            ? [
                'position', 'aa',
                ...ptmCats.map(c => `ptm_${safe(c)}`),
                ...diseases.map(d => `disease_${safe(d)}`),
                'am_avg_score',
                'residue_burden',
                ...chains.flatMap(c => [`pdb_residue_${c}`, `hotspot_tier_${c}`, `contact_hub_tier_${c}`]),
            ]
            : [
                'position', 'aa', 'pdb_residue',
                ...ptmCats.map(c => `ptm_${safe(c)}`),
                ...diseases.map(d => `disease_${safe(d)}`),
                'am_avg_score',
                'hotspot_tier',
                'contact_hub_tier',
                'residue_burden',
            ];
        const rows = [headers.join(',')];
        sequence.split('').forEach((aa, i) => {
            const pos = i + 1;
            const ptmFlags = ptmCats.map(c => ptmByPos.get(pos)?.has(c) ? 1 : 0);
            const diseaseFlags = diseases.map(d => diseaseByPos.get(pos)?.has(d) ? 1 : 0);
            let amAvg = '';
            if (amMap && amMap.size > 0) {
                const scores = [];
                for (const mut of AM_AAS) {
                    if (mut === aa) continue;
                    const sc = amMap.get(`${aa}${pos}${mut}`);
                    if (Number.isFinite(sc)) scores.push(sc);
                }
                if (scores.length > 0) amAvg = (scores.reduce((a, b) => a + b, 0) / scores.length).toFixed(4);
            }
            const burden = analysis.residueBurden instanceof Set && analysis.residueBurden.has(pos) ? 1 : 0;
            if (chains) {
                const structVals = chains.flatMap(c => [
                    uniprotToPdbResiForChain(pos, structure, c) ?? '',
                    tierFor(analysis.hotspotsByChain?.get(c), pos, hotspotTierNum),
                    tierFor(analysis.distantContactsByChain?.get(c), pos, hubTierNum),
                ]);
                rows.push([pos, aa, ...ptmFlags, ...diseaseFlags, amAvg, burden, ...structVals].join(','));
            } else {
                const pdbResi = uniprotToPdbResi(pos, structure) ?? '';
                const hotspotTier = tierFor(analysis.hotspots, pos, hotspotTierNum);
                const hubTier = tierFor(analysis.distantContacts, pos, hubTierNum);
                rows.push([pos, aa, pdbResi, ...ptmFlags, ...diseaseFlags, amAvg, hotspotTier, hubTier, burden].join(','));
            }
        });
        return rows.join('\n');
    }

    return { formatSelection, rewritePdbBeta, buildResidueMatrix, downloadText, copyText };
})();
