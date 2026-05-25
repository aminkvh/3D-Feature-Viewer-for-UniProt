const UFVAnalysis = (() => {
    'use strict';

    function residueCaAtoms(viewer) {
        const atoms = viewer?.getModel?.()?.selectedAtoms?.({ atom: 'CA' }) || [];
        const out = new Map();
        atoms.forEach(a => {
            if (!out.has(a.resi)) out.set(a.resi, a);
        });
        return out;
    }

    /**
     * Return a Map keyed by chain → (Map of pdbResi → atom).
     * Used for per-subunit hotspot calculations.
     */
    function residueCaAtomsPerChain(viewer) {
        const atoms = viewer?.getModel?.()?.selectedAtoms?.({ atom: 'CA' }) || [];
        const out = new Map();
        atoms.forEach(a => {
            if (!out.has(a.chain)) out.set(a.chain, new Map());
            if (!out.get(a.chain).has(a.resi)) out.get(a.chain).set(a.resi, a);
        });
        return out;
    }

    function distance(a, b) {
        const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    function residueNeighborhood(viewer, resi, threshold = 5) {
        const model = viewer?.getModel?.();
        const atoms = model?.selectedAtoms?.({}) || [];
        const selected = atoms.filter(a => a.resi === resi);
        const nearby = new Set([resi]);
        atoms.forEach(atom => {
            if (atom.resi === resi) return;
            for (const s of selected) {
                if (distance(atom, s) <= threshold) {
                    nearby.add(atom.resi);
                    break;
                }
            }
        });
        return nearby;
    }

    function mapUniToPdb(pos, structure) {
        if (!structure?.mappedRanges?.length) return pos;
        for (const r of structure.mappedRanges) {
            if (pos >= r.uniprotStart && pos <= r.uniprotEnd) {
                // Chimeric/non-SIFTS structures store a seqresToAuthor map on the structure
                // object and save the original SEQRES start as r.seqresStart.  Use it here
                // so that 3-D analyses find the correct ATOM record author residue number.
                if (structure.seqresToAuthor && r.seqresStart != null) {
                    const seqresPos = r.seqresStart + (pos - r.uniprotStart);
                    return structure.seqresToAuthor.get(seqresPos) ?? null;
                }
                return r.pdbStart + (pos - r.uniprotStart);
            }
        }
        return null;
    }

    function mapPdbToUni(pos, structure) {
        if (!structure?.mappedRanges?.length) return pos;
        for (const r of structure.mappedRanges) {
            if (structure.seqresToAuthor && r.seqresStart != null) {
                // For chimeric structures author residue numbers are non-linear; reverse the
                // seqresToAuthor map by scanning through the range to find the UniProt offset.
                const rangeLen = r.uniprotEnd - r.uniprotStart;
                for (let i = 0; i <= rangeLen; i++) {
                    if (structure.seqresToAuthor.get(r.seqresStart + i) === pos) {
                        return r.uniprotStart + i;
                    }
                }
            } else if (pos >= r.pdbStart && pos <= r.pdbEnd) {
                return r.uniprotStart + (pos - r.pdbStart);
            }
        }
        return null;
    }

    function logChoose(n, k) {
        if (k < 0 || k > n) return -Infinity;
        k = Math.min(k, n - k);
        let r = 0;
        for (let i = 0; i < k; i++) r += Math.log(n - i) - Math.log(i + 1);
        return r;
    }

    function fisherExactOneSided(a, b, c, d) {
        // One-sided Fisher's exact test: P(X >= a) where X ~ Hypergeometric
        // Table: [[a,b],[c,d]]; testing over-representation of 'a' in row 1
        const r1 = a + b, r2 = c + d, c1 = a + c, n = r1 + r2;
        if (n === 0 || c1 === 0 || r1 === 0) return 1;
        const logDenom = logChoose(n, c1);
        const kMax = Math.min(r1, c1);
        let pVal = 0;
        for (let k = a; k <= kMax; k++) {
            const logP = logChoose(r1, k) + logChoose(r2, c1 - k) - logDenom;
            if (Number.isFinite(logP)) pVal += Math.exp(logP);
        }
        return Math.min(1, pVal);
    }

    function computeHotspots(viewer, variants, structure = null, threshold = 8) {
        const pathogenicPositions = new Set();
        variants.forEach(v => {
            if (/pathogenic|deleterious/i.test(v.consequence || '')) pathogenicPositions.add(v.position);
        });
        const benignOnlyPositions = new Set();
        variants.forEach(v => {
            if (!pathogenicPositions.has(v.position) && /benign/i.test(v.consequence || '')) benignOnlyPositions.add(v.position);
        });

        const globalPath = pathogenicPositions.size;
        const globalBenign = benignOnlyPositions.size;
        const globalTotal = globalPath + globalBenign;
        if (globalTotal === 0) return new Map();
        const globalPathFraction = globalPath / globalTotal;

        /**
         * Compute hotspot stats for a single chain's CA atom map.
         * For multi-chain structures each subunit is analysed independently using its own
         * spatial context so that per-subunit differences are captured; the caller then
         * promotes each UniProt position to the highest tier seen across all chains.
         * @param {Map<number,object>} ca  pdbResi → CA atom for the chain being analysed
         * @param {function(number):number|null} uniToPdb  UniProt pos → PDB resi for this chain
         */
        function computeChainStats(ca, uniToPdb) {
            const stats = [];
            pathogenicPositions.forEach(pos => {
                const pdbPos = uniToPdb(pos);
                const a = ca.get(pdbPos);
                if (!a) return;
                let localPath = 1; // includes self
                let localBenign = 0;
                pathogenicPositions.forEach(q => {
                    if (q === pos) return;
                    const pdbQ = uniToPdb(q);
                    if (ca.has(pdbQ) && distance(a, ca.get(pdbQ)) <= threshold) localPath++;
                });
                benignOnlyPositions.forEach(q => {
                    const pdbQ = uniToPdb(q);
                    if (ca.has(pdbQ) && distance(a, ca.get(pdbQ)) <= threshold) localBenign++;
                });
                const localTotal = localPath + localBenign;
                const localPathFraction = localTotal > 0 ? localPath / localTotal : 0;
                const enrichmentRatio = localPathFraction / Math.max(globalPathFraction, 1e-9);
                stats.push({ pos, localPath, localBenign, localTotal, localPathFraction, enrichmentRatio });
            });
            return stats;
        }

        // Build per-chain CA maps when the structure has multiple chains so each
        // subunit's 3-D context is evaluated independently.
        const chainIds = structure?.chainIds?.length > 1 ? structure.chainIds : null;
        const caPerChain = chainIds ? residueCaAtomsPerChain(viewer) : null;

        // Fallback: single chain or no structure — behave as before.
        const fallbackCa = caPerChain ? null : residueCaAtoms(viewer);

        // Collect stats across all chains; a UniProt position may appear multiple times
        // (once per chain) — we take the entry with the highest localPath to represent
        // the most pathogenic spatial context for that residue.
        const bestByPos = new Map(); // UniProt pos → best stats entry
        const accumulateStats = stats => {
            stats.forEach(s => {
                const prev = bestByPos.get(s.pos);
                if (!prev || s.localPath > prev.localPath ||
                    (s.localPath === prev.localPath && s.enrichmentRatio > prev.enrichmentRatio)) {
                    bestByPos.set(s.pos, s);
                }
            });
        };

        if (chainIds) {
            chainIds.forEach(chain => {
                const ca = caPerChain.get(chain);
                if (!ca?.size) return;
                const chainRanges = structure.chainMappings?.[chain] || structure.mappedRanges;
                const chainMap = structure.chainSeqresToAuthor?.[chain] || structure.seqresToAuthor || null;
                const uniToPdb = pos => {
                    for (const r of chainRanges) {
                        if (pos >= r.uniprotStart && pos <= r.uniprotEnd) {
                            if (chainMap && r.seqresStart != null) {
                                const seqresPos = r.seqresStart + (pos - r.uniprotStart);
                                return chainMap.get(seqresPos) ?? null;
                            }
                            return r.pdbStart + (pos - r.uniprotStart);
                        }
                    }
                    return null;
                };
                accumulateStats(computeChainStats(ca, uniToPdb));
            });
        } else {
            const uniToPdb = pos => mapUniToPdb(pos, structure);
            accumulateStats(computeChainStats(fallbackCa, uniToPdb));
        }

        const stats = Array.from(bestByPos.values());

        // BH-corrected Fisher p-values for strong-tier candidates only
        const strongCandidates = stats.filter(s => s.localPath >= 3 && s.enrichmentRatio >= 2.0);
        const strongAdjP = new Map();
        if (strongCandidates.length > 0) {
            strongCandidates.forEach(s => {
                const outsidePath = Math.max(0, globalPath - s.localPath);
                const outsideBenign = Math.max(0, globalBenign - s.localBenign);
                s.p = fisherExactOneSided(s.localPath, s.localBenign, outsidePath, outsideBenign);
            });
            strongCandidates.sort((a, b) => a.p - b.p);
            const N = strongCandidates.length;
            strongCandidates.forEach((c, i) => { c.adjP = Math.min(c.p * N / (i + 1), 1); });
            for (let i = N - 2; i >= 0; i--) strongCandidates[i].adjP = Math.min(strongCandidates[i].adjP, strongCandidates[i + 1].adjP);
            strongCandidates.forEach(c => strongAdjP.set(c.pos, c.adjP));
        }

        const result = new Map();
        stats.forEach(s => {
            // Strong: localPath ≥ 3, enrichmentRatio ≥ 2.0, Fisher q ≤ 0.10
            if (s.localPath >= 3 && s.enrichmentRatio >= 2.0 && (strongAdjP.get(s.pos) ?? 1) <= 0.10) {
                result.set(s.pos, 'strong');
                return;
            }
            // Moderate: localPath ≥ 2, enrichmentRatio ≥ 1.5, localPathFraction ≥ 0.50
            if (s.localPath >= 2 && s.enrichmentRatio >= 1.5 && s.localPathFraction >= 0.50) {
                result.set(s.pos, 'moderate');
                return;
            }
            // Weak: localPath ≥ 2, localPathFraction ≥ 0.40
            if (s.localPath >= 2 && s.localPathFraction >= 0.40) {
                result.set(s.pos, 'weak');
            }
        });
        return result;
    }

    function computeDistantContacts(viewer, structure, variants = [], distanceThreshold = 8, sequenceThreshold = 40) {
        const caMap = residueCaAtoms(viewer);
        if (caMap.size === 0) return new Map();
        const ca = Array.from(caMap.entries()); // [pdbResi, atom]

        // Build pathogenic residue set in PDB coordinates
        const pathogenicPdb = new Set();
        variants.forEach(v => {
            if (/pathogenic|deleterious/i.test(v.consequence || '')) {
                const pdb = mapUniToPdb(v.position, structure);
                if (pdb != null) pathogenicPdb.add(pdb);
            }
        });

        // Score each residue's long-range contacts:
        // contactScore = longRangeCount + 2*pathogenicPartners + 2*cross-chain partners
        const scores = new Map(); // pdbResi → { count, total }
        for (let i = 0; i < ca.length; i++) {
            for (let j = i + 1; j < ca.length; j++) {
                const [aResi, a] = ca[i];
                const [bResi, b] = ca[j];
                if (Math.abs(aResi - bResi) < sequenceThreshold) continue;
                if (distance(a, b) > distanceThreshold) continue;
                const isInterface = a.chain && b.chain && a.chain !== b.chain;
                const addScore = (resi, partnerIsPath) => {
                    const s = scores.get(resi) || { count: 0, total: 0 };
                    s.count++;
                    s.total += 1 + (partnerIsPath ? 2 : 0) + (isInterface ? 2 : 0);
                    scores.set(resi, s);
                };
                addScore(aResi, pathogenicPdb.has(bResi));
                addScore(bResi, pathogenicPdb.has(aResi));
            }
        }

        // Require at least 3 long-range contacts
        const qualified = Array.from(scores.entries())
            .filter(([, s]) => s.count >= 3)
            .map(([pdbResi, s]) => ({ pdbResi, score: s.total }));
        if (qualified.length === 0) return new Map();

        qualified.sort((a, b) => b.score - a.score);
        const top5  = Math.max(1, Math.ceil(0.05 * qualified.length));
        const top15 = Math.max(1, Math.ceil(0.15 * qualified.length));
        const result = new Map();
        qualified.forEach(({ pdbResi }, idx) => {
            const tier = idx < top5 ? 'strong' : idx < top15 ? 'moderate' : null;
            if (!tier) return;
            const uniPos = mapPdbToUni(pdbResi, structure) ?? pdbResi;
            result.set(uniPos, tier);
        });
        return result;
    }

    function aggregateAlphaMissense(variants, amMap = null) {
        // If the full AlphaMissense CSV map is available use it — it covers every residue
        // with all 19 possible substitutions, avoiding 'No score' on unannotated positions.
        if (amMap && amMap.size > 0) {
            const byPos = new Map();
            amMap.forEach((score, key) => {
                const m = key.match(/^[A-Z](-?\d+)[A-Z]$/);
                if (!m) return;
                const pos = Number(m[1]);
                if (!byPos.has(pos)) byPos.set(pos, []);
                byPos.get(pos).push(score);
            });
            const out = new Map();
            byPos.forEach((scores, pos) => {
                out.set(pos, {
                    avg: scores.reduce((a, b) => a + b, 0) / scores.length,
                    max: Math.max(...scores),
                    count: scores.length,
                });
            });
            return out;
        }
        // Fallback: aggregate from annotated variants only (used when CSV is unavailable)
        const byPos = new Map();
        variants.forEach(v => {
            const score = Number(v.alphaMissenseScore);
            if (!Number.isFinite(score)) return;
            if (!byPos.has(v.position)) byPos.set(v.position, []);
            byPos.get(v.position).push(score);
        });
        const out = new Map();
        byPos.forEach((scores, pos) => {
            out.set(pos, {
                avg: scores.reduce((a, b) => a + b, 0) / scores.length,
                max: Math.max(...scores),
                count: scores.length,
            });
        });
        return out;
    }

    function computeResidueBurden(variants) {
        const mutCountByPos = new Map();
        const phenByPos = new Map();
        variants.forEach(v => {
            const pos = v.position;
            mutCountByPos.set(pos, (mutCountByPos.get(pos) || 0) + 1);
            if (!phenByPos.has(pos)) phenByPos.set(pos, new Set());
            (v.diseasePairs || []).forEach(p => {
                const key = p.id || p.label;
                if (key) phenByPos.get(pos).add(key);
            });
        });
        if (mutCountByPos.size === 0) return new Set();
        const counts = [...mutCountByPos.values()].sort((a, b) => a - b);
        const p90Idx = Math.ceil(0.9 * counts.length) - 1;
        const p90 = counts[Math.max(0, p90Idx)];
        const flagged = new Set();
        mutCountByPos.forEach((count, pos) => {
            const phenCount = phenByPos.get(pos)?.size || 0;
            if (count >= 3 && phenCount >= 2 && count >= p90) flagged.add(pos);
        });
        return flagged;
    }

    return { residueNeighborhood, computeHotspots, computeDistantContacts, aggregateAlphaMissense, computeResidueBurden };
})();
