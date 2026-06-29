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

    /**
     * Build a UniProt→PDB-author-residue mapper for ONE chain of a structure.
     * For multi-chain structures each chain has its own mappedRanges (chainMappings) and,
     * for chimeric chains, its own seqresToAuthor map (chainSeqresToAuthor), so distance
     * calculations resolve to the correct author residue number for that subunit.
     * chain === null falls back to the structure-level (primary chain) mapping.
     */
    function uniToPdbForChain(structure, chain) {
        const ranges = (chain != null && structure?.chainMappings?.[chain]) || structure?.mappedRanges || [];
        const map = (chain != null && structure?.chainSeqresToAuthor?.[chain]) || structure?.seqresToAuthor || null;
        return pos => {
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
        };
    }

    /** Inverse of uniToPdbForChain: PDB author residue number → UniProt position for one chain. */
    function pdbToUniForChain(structure, chain) {
        const ranges = (chain != null && structure?.chainMappings?.[chain]) || structure?.mappedRanges || [];
        const map = (chain != null && structure?.chainSeqresToAuthor?.[chain]) || structure?.seqresToAuthor || null;
        return pdb => {
            if (!ranges.length) return pdb;
            for (const r of ranges) {
                if (map && r.seqresStart != null) {
                    const rangeLen = r.uniprotEnd - r.uniprotStart;
                    for (let i = 0; i <= rangeLen; i++) {
                        if (map.get(r.seqresStart + i) === pdb) return r.uniprotStart + i;
                    }
                } else if (pdb >= r.pdbStart && pdb <= r.pdbEnd) {
                    return r.uniprotStart + (pdb - r.pdbStart);
                }
            }
            return null;
        };
    }

    /** Combine per-chain tier Maps into one Map keeping the highest tier seen per position. */
    function mergeByChain(byChain, rank) {
        const merged = new Map();
        byChain.forEach(m => m.forEach((tier, pos) => {
            if (!merged.has(pos) || (rank[tier] || 0) > (rank[merged.get(pos)] || 0)) merged.set(pos, tier);
        }));
        return merged;
    }

    // Deterministic, seedable PRNG (mulberry32) so the permutation-based hotspot test gives
    // reproducible results across runs/sessions instead of flickering on each recompute.
    function mulberry32(seed) {
        let a = seed >>> 0;
        return function () {
            a |= 0; a = (a + 0x6D2B79F5) | 0;
            let t = Math.imul(a ^ (a >>> 15), 1 | a);
            t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
            return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
        };
    }

    // Number of label permutations and PRNG seed for the spatial null (see computeHotspots).
    const HOTSPOT_PERMUTATIONS = 1000;
    const HOTSPOT_SEED = 0x9e3779b9;
    // Minimum classified-variant counts for the permutation test to be meaningful: with too few
    // pathogenic or benign labels there is no spatial contrast to detect and the empirical null
    // is degenerate, so we stay silent rather than emit noise.
    const HOTSPOT_MIN_PATH = 3;
    const HOTSPOT_MIN_BENIGN = 2;

    /**
     * 3-D enrichment hotspots via a spatial label-permutation null.
     *
     * For each chain we fix the modelled positions of the classified variants (pathogenic vs
     * benign) and their 3-D coordinates, then build a local pathogenic-density statistic for
     * every classified position (number of pathogenic labels within `threshold` Å, including
     * self).  Significance is assessed by randomly permuting the pathogenic/benign labels over
     * the fixed positions many times and comparing each observed local density to its empirical
     * null — this respects the spatial autocorrelation and non-independence of overlapping
     * neighbourhoods that a per-residue Fisher test ignores.
     *
     * This Monte-Carlo / permutation approach to spatial clustering of (germline pathogenic)
     * missense variants on protein structures follows established methods:
     *   • Kamburov et al., PNAS 2015 (doi:10.1073/pnas.1516373112) — empirical null from
     *     permuting mutated positions across the structure.
     *   • Tokheim et al. "HotMAPS", Cancer Res 2016 (doi:10.1158/0008-5472.CAN-15-3190) —
     *     significantly increased local mutation density vs. an empirical null.
     *   • Sivley et al., Am J Hum Genet 2018 (doi:10.1016/j.ajhg.2018.01.017) — spatial
     *     clustering of germline pathogenic vs. neutral variants in human protein structures.
     * Empirical p-values use the (b+1)/(m+1) estimator (Phipson & Smyth, SAGMB 2010) so a
     * p-value is never exactly zero.
     */
    function computeHotspots(viewer, variants, structure = null, threshold = 8, partnerPoints = []) {
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
        if (globalPath < HOTSPOT_MIN_PATH) return { merged: new Map(), byChain: new Map() };
        // With enough benign controls, use the case-control label-permutation test (it also
        // controls for which residues were studied/sequenced — the stronger test).  Otherwise
        // fall back to a spatial-placement null (pathogenic clustering vs. random placement
        // over the whole chain), which needs NO benign controls — the same approach as HotMAPS
        // / Kamburov — so benign-poor disease genes still get hotspot results.
        const method = globalBenign >= HOTSPOT_MIN_BENIGN ? 'casecontrol' : 'spatial';

        // Resolve neighbouring partner-protein classified residues to CA atoms (other chains
        // of a complex, e.g. other GABA-A subunits).  These are folded into the test so a
        // residue of OUR protein sitting next to a partner's pathogenic cluster is scored
        // accordingly — but they are never reported/visualised (report:false).
        const partnerItems = [];
        if (partnerPoints && partnerPoints.length) {
            const caAll = residueCaAtomsPerChain(viewer);
            partnerPoints.forEach(pp => {
                const atom = caAll.get(pp.chainId)?.get(pp.pdbResi);
                if (atom) partnerItems.push({ pos: null, path: !!pp.path, atom, report: false });
            });
        }
        const partnerPathAtoms = partnerItems.filter(p => p.path).map(p => p.atom);

        /**
         * Tier the classified positions of one chain using the spatial label-permutation null.
         * @param {Map<number,object>} ca  pdbResi → CA atom for the chain being analysed
         * @param {function(number):number|null} uniToPdb  UniProt pos → PDB resi for this chain
         * @returns {Map<number,string>} UniProt position → tier ('strong'|'moderate'|'weak')
         */
        function computeChainTiers(ca, uniToPdb) {
            // Collect classified positions of OUR protein that are modelled in this chain.
            const items = []; // { pos, path, atom, report }
            const collect = (set, isPath) => set.forEach(pos => {
                const atom = ca.get(uniToPdb(pos));
                if (atom) items.push({ pos, path: isPath, atom, report: true });
            });
            collect(pathogenicPositions, true);
            collect(benignOnlyPositions, false);
            const ourCount = items.length;
            if (ourCount === 0) return new Map();

            // Add partner residues that are spatially close to this chain to the test pool.
            partnerItems.forEach(pi => {
                for (let q = 0; q < ourCount; q++) {
                    if (distance(pi.atom, items[q].atom) <= threshold) { items.push(pi); break; }
                }
            });

            const k = items.length;
            if (k < HOTSPOT_MIN_PATH + 1) return new Map();

            // Fixed neighbourhood structure (indices within threshold; includes self).
            const neigh = Array.from({ length: k }, () => []);
            for (let i = 0; i < k; i++) {
                neigh[i].push(i);
                for (let j = i + 1; j < k; j++) {
                    if (distance(items[i].atom, items[j].atom) <= threshold) { neigh[i].push(j); neigh[j].push(i); }
                }
            }

            const labels = items.map(it => it.path);
            const localPath = new Array(k).fill(0);
            for (let i = 0; i < k; i++) for (const j of neigh[i]) if (labels[j]) localPath[i]++;
            // Effect-size baseline = pathogenic fraction of the whole pool (incl. partners).
            const baseFraction = labels.reduce((a, b) => a + (b ? 1 : 0), 0) / k;

            // Permute labels over the fixed positions; count how often the permuted local
            // pathogenic density meets/exceeds the observed value at each position.
            const ge = new Array(k).fill(0);
            const perm = labels.slice();
            const rand = mulberry32(HOTSPOT_SEED);
            for (let p = 0; p < HOTSPOT_PERMUTATIONS; p++) {
                for (let i = k - 1; i > 0; i--) { const r = Math.floor(rand() * (i + 1)); const t = perm[i]; perm[i] = perm[r]; perm[r] = t; }
                for (let i = 0; i < k; i++) {
                    let c = 0; for (const j of neigh[i]) if (perm[j]) c++;
                    if (c >= localPath[i]) ge[i]++;
                }
            }

            // Empirical p (Phipson & Smyth), BH-FDR across OUR pathogenic centres only
            // (partner residues contribute to the null/density but are never reported).
            const centres = [];
            for (let i = 0; i < k; i++) {
                if (!items[i].path || !items[i].report) continue;
                const localTotal = neigh[i].length;
                const localPathFraction = localTotal > 0 ? localPath[i] / localTotal : 0;
                centres.push({
                    pos: items[i].pos,
                    p: (ge[i] + 1) / (HOTSPOT_PERMUTATIONS + 1),
                    localPath: localPath[i],
                    localPathFraction,
                    enrichmentRatio: localPathFraction / Math.max(baseFraction, 1e-9),
                });
            }
            centres.sort((a, b) => a.p - b.p);
            const N = centres.length;
            centres.forEach((c, i) => { c.q = Math.min(c.p * N / (i + 1), 1); });
            for (let i = N - 2; i >= 0; i--) centres[i].q = Math.min(centres[i].q, centres[i + 1].q);

            const result = new Map();
            centres.forEach(c => {
                // Significance (empirical FDR) gates the tier; effect size separates strong/moderate.
                if (c.q <= 0.10 && c.localPath >= 3 && c.enrichmentRatio >= 2.0) result.set(c.pos, 'strong');
                else if (c.q <= 0.25 && c.localPath >= 2 && c.enrichmentRatio >= 1.5) result.set(c.pos, 'moderate');
                else if (c.q <= 0.25 && c.localPath >= 2 && c.localPathFraction >= 0.40) result.set(c.pos, 'weak');
            });
            return result;
        }

        /**
         * Spatial-placement null (HotMAPS / Kamburov) for benign-poor proteins: score each
         * pathogenic residue by its local pathogenic density (neighbours within threshold, plus
         * nearby partner-protein pathogenic residues) against the density expected if the same
         * number of pathogenic residues were placed at random over ALL modelled residues of the
         * chain.  No benign controls needed.
         * @param {Map<number,object>} ca  pdbResi → CA atom for the chain
         * @param {function(number):number|null} uniToPdb  UniProt pos → PDB resi
         * @param {function(number):number|null} toUni     PDB resi → UniProt pos (for reporting)
         */
        function computeChainTiersSpatial(ca, uniToPdb, toUni) {
            const universe = [];
            ca.forEach((atom, resi) => universe.push({ resi, atom }));
            const N = universe.length;
            if (N < 10 || N > 6000) return new Map(); // too small to test / too large to be cheap
            const resiToIdx = new Map(universe.map((u, i) => [u.resi, i]));
            const pathIdx = [];
            pathogenicPositions.forEach(pos => { const idx = resiToIdx.get(uniToPdb(pos)); if (idx != null) pathIdx.push(idx); });
            const m = pathIdx.length;
            if (m < HOTSPOT_MIN_PATH) return new Map();

            // Neighbourhoods among the chain's residues (self excluded) + fixed partner density.
            const neigh = Array.from({ length: N }, () => []);
            for (let i = 0; i < N; i++) {
                for (let j = i + 1; j < N; j++) {
                    if (distance(universe[i].atom, universe[j].atom) <= threshold) { neigh[i].push(j); neigh[j].push(i); }
                }
            }
            const partnerDen = new Array(N).fill(0);
            partnerPathAtoms.forEach(pa => {
                for (let i = 0; i < N; i++) if (distance(pa, universe[i].atom) <= threshold) partnerDen[i]++;
            });

            const isPath = new Uint8Array(N);
            pathIdx.forEach(i => { isPath[i] = 1; });
            const obsDen = new Array(N);
            for (let i = 0; i < N; i++) { let c = partnerDen[i]; for (const j of neigh[i]) if (isPath[j]) c++; obsDen[i] = c; }

            // Empirical null: place m pathogenic residues uniformly at random over the chain.
            const ge = new Array(N).fill(0);
            const rand = mulberry32(HOTSPOT_SEED);
            const pool = [...Array(N).keys()];
            const simPath = new Uint8Array(N);
            for (let p = 0; p < HOTSPOT_PERMUTATIONS; p++) {
                simPath.fill(0);
                for (let i = 0; i < m; i++) { const r = i + Math.floor(rand() * (N - i)); const t = pool[i]; pool[i] = pool[r]; pool[r] = t; simPath[pool[i]] = 1; }
                for (const ci of pathIdx) {
                    let c = partnerDen[ci]; for (const j of neigh[ci]) if (simPath[j]) c++;
                    if (c >= obsDen[ci]) ge[ci]++;
                }
            }

            const centres = pathIdx.map(ci => ({ ci, p: (ge[ci] + 1) / (HOTSPOT_PERMUTATIONS + 1), den: obsDen[ci] }));
            centres.sort((a, b) => a.p - b.p);
            const Nn = centres.length;
            centres.forEach((c, i) => { c.q = Math.min(c.p * Nn / (i + 1), 1); });
            for (let i = Nn - 2; i >= 0; i--) centres[i].q = Math.min(centres[i].q, centres[i + 1].q);

            const result = new Map();
            centres.forEach(c => {
                const uni = toUni(universe[c.ci].resi);
                if (uni == null) return;
                if (c.q <= 0.10 && c.den >= 3) result.set(uni, 'strong');
                else if (c.q <= 0.25 && c.den >= 2) result.set(uni, 'moderate');
            });
            return result;
        }

        // Tier each chain independently and keep the per-chain results (the user wants the
        // structure-dependent hotspot calls separated, not pooled across subunits).
        const byChain = new Map(); // chainId → Map<uniProtPos, tier>
        const multi = structure?.chainIds?.length > 1;
        const tiersFor = (ca, uniToPdb, toUni) => method === 'casecontrol'
            ? computeChainTiers(ca, uniToPdb)
            : computeChainTiersSpatial(ca, uniToPdb, toUni);
        if (multi) {
            const caPerChain = residueCaAtomsPerChain(viewer);
            structure.chainIds.forEach(chain => {
                const ca = caPerChain.get(chain);
                if (!ca?.size) return;
                byChain.set(chain, tiersFor(ca, uniToPdbForChain(structure, chain), pdbToUniForChain(structure, chain)));
            });
        } else {
            const key = structure?.chainId ?? null;
            // Resolve our chain's CA atoms by chain id (avoids residue-number collisions with
            // partner chains in a complex); AlphaFold has no chain id, so use the flat map.
            const ca = key != null ? (residueCaAtomsPerChain(viewer).get(key) || new Map()) : residueCaAtoms(viewer);
            byChain.set(key, tiersFor(ca, pos => mapUniToPdb(pos, structure), pdb => mapPdbToUni(pdb, structure)));
        }

        const merged = mergeByChain(byChain, { strong: 3, moderate: 2, weak: 1 });
        return { merged, byChain, method };
    }

    // Betweenness z-score cut-offs for contact-hub tiers (absolute statistical thresholds on
    // the observed distribution, not a fixed percentile — so a featureless chain emits nothing).
    const HUB_Z_STRONG = 3.0;
    const HUB_Z_MODERATE = 2.0;
    // Brandes' betweenness is O(V·E); for an 8 Å contact graph |E| ≈ 4V, so cost ≈ O(4V²):
    // ~20 ms at 2000 residues, ~80 ms at 4000, ~0.6 s at 10000. The analysis is run off the main
    // thread (requestIdleCallback) so the only thing this cap protects against is a pathologically
    // huge chain (e.g. a giant complex parsed as one chain) stalling for seconds. It is set well
    // above any real single protein chain — a 2000-residue protein is computed normally — so it
    // never nerfs ordinary entries; above it we skip the hub overlay for that chain.
    const HUB_MAX_CA = 6000;

    /**
     * Identify "long-range contact hubs" as residues with high BETWEENNESS CENTRALITY in the
     * residue contact network (nodes = Cα atoms, edges = Cα–Cα ≤ threshold Å).
     *
     * Why betweenness instead of raw contact count: a plain contact-degree score mostly
     * rediscovers buried core residues (degree ≈ burial), which is a confound.  Betweenness
     * instead rewards residues that lie on many shortest paths — i.e. that bridge otherwise
     * distant regions/domains.  Residue-interaction networks and short-path communication have
     * established precedent in protein-structure analysis:
     *   • del Sol et al., Mol Syst Biol 2006 (doi:10.1038/msb4100063)
     *   • Amitai et al., J Mol Biol 2004 (doi:10.1016/j.jmb.2004.10.055)
     * Centrality is computed with Brandes' algorithm (J Math Sociol 2001) on the unweighted
     * graph; residues are flagged by an absolute z-score on the betweenness distribution
     * rather than always taking a top percentile.
     *
     * @param {object[]} caAtoms  Cα atoms for one chain
     * @param {number} threshold  contact distance cut-off (Å)
     * @param {function(number):number|null} toUni  author resi → UniProt position
     * @returns {Map<number,string>} UniProt position → tier ('strong'|'moderate')
     */
    function betweennessHubs(caAtoms, threshold, toUni) {
        const n = caAtoms.length;
        if (n < 8) return new Map(); // too small for a meaningful network
        if (n > HUB_MAX_CA) { // too large — skip rather than freeze the main thread (see HUB_MAX_CA)
            console.warn(`[UniProt 3D] Skipping contact-hub analysis for a ${n}-residue chain (cap ${HUB_MAX_CA}).`);
            return new Map();
        }
        // Build adjacency by Cα–Cα distance.
        const th2 = threshold * threshold;
        const adj = Array.from({ length: n }, () => []);
        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                const a = caAtoms[i], b = caAtoms[j];
                const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
                if (dx * dx + dy * dy + dz * dz <= th2) { adj[i].push(j); adj[j].push(i); }
            }
        }
        // Brandes' betweenness centrality (unweighted, undirected).
        const bc = new Array(n).fill(0);
        for (let s = 0; s < n; s++) {
            const stack = [];
            const pred = Array.from({ length: n }, () => []);
            const sigma = new Array(n).fill(0); sigma[s] = 1;
            const dist = new Array(n).fill(-1); dist[s] = 0;
            const queue = [s]; let qi = 0;
            while (qi < queue.length) {
                const v = queue[qi++]; stack.push(v);
                for (const w of adj[v]) {
                    if (dist[w] < 0) { dist[w] = dist[v] + 1; queue.push(w); }
                    if (dist[w] === dist[v] + 1) { sigma[w] += sigma[v]; pred[w].push(v); }
                }
            }
            const delta = new Array(n).fill(0);
            for (let idx = stack.length - 1; idx >= 0; idx--) {
                const w = stack[idx];
                for (const v of pred[w]) delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w]);
                if (w !== s) bc[w] += delta[w];
            }
        }
        // Normalise (undirected): each pair counted twice, then scale by max possible pairs.
        const scale = 1 / ((n - 1) * (n - 2));
        for (let i = 0; i < n; i++) bc[i] *= scale;
        // Absolute cut-off via z-score on the betweenness distribution.
        const mean = bc.reduce((a, b) => a + b, 0) / n;
        const std = Math.sqrt(bc.reduce((a, b) => a + (b - mean) * (b - mean), 0) / n);
        if (!(std > 0)) return new Map();
        const result = new Map();
        for (let i = 0; i < n; i++) {
            if (bc[i] <= 0) continue;
            const z = (bc[i] - mean) / std;
            const tier = z >= HUB_Z_STRONG ? 'strong' : z >= HUB_Z_MODERATE ? 'moderate' : null;
            if (!tier) continue;
            const uni = toUni(caAtoms[i].resi);
            if (uni == null) continue; // residue outside the UniProt mapping (tags, ligands)
            if (result.get(uni) !== 'strong') result.set(uni, tier);
        }
        return result;
    }

    function computeDistantContacts(viewer, structure, variants = [], distanceThreshold = 8) {
        const atoms = (viewer?.getModel?.()?.selectedAtoms?.({ atom: 'CA' }) || []);
        if (atoms.length === 0) return { merged: new Map(), byChain: new Map() };
        const multi = structure?.chainIds?.length > 1;

        // Per-chain betweenness so each subunit's own packing/topology is scored independently.
        const byChain = new Map(); // chainId → Map<uniProtPos, tier>
        if (multi) {
            structure.chainIds.forEach(chain => {
                const chainAtoms = atoms.filter(a => a.chain === chain);
                const tiers = betweennessHubs(chainAtoms, distanceThreshold, pdbToUniForChain(structure, chain));
                if (tiers.size) byChain.set(chain, tiers);
            });
        } else {
            const cid = structure?.chainId || null;
            const chainAtoms = cid ? atoms.filter(a => a.chain === cid) : atoms;
            const tiers = betweennessHubs(chainAtoms, distanceThreshold, pdb => mapPdbToUni(pdb, structure));
            if (tiers.size) byChain.set(cid, tiers);
        }

        const merged = mergeByChain(byChain, { strong: 2, moderate: 1 });
        return { merged, byChain };
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

    /**
     * Residue-burden hotspots via a percentile-rank composite of two per-position metrics:
     *   (1) number of distinct phenotypes (pleiotropy)   — the reporting-bias-resistant signal
     *   (2) number of reported mutations (recurrence)
     * This mirrors the multi-metric hotspot prioritisation in Akbari Ahangar et al., iScience
     * 2024 (PMID 39286500), which ranks positions by phenotype count, mutation count, and the
     * number of proteins mutated at equivalent positions.  The third (cross-paralog) metric
     * needs family/MSA data not available from a single UniProt entry, so it is omitted here.
     *
     * Each metric is converted to a percentile rank across positions (robust to the absolute
     * scale and to the small, skewed count vectors typical of variant data) and averaged.
     * A position is flagged when its composite rank is high AND it shows genuine recurrence +
     * pleiotropy (absolute floors), so sparse proteins don't trigger on singletons.
     */
    function computeResidueBurden(variants) {
        const mutCountByPos = new Map();
        const phenByPos = new Map();
        const identifiedLabelsByPos = new Map();
        variants.forEach(v => {
            const pos = v.position;
            mutCountByPos.set(pos, (mutCountByPos.get(pos) || 0) + 1);
            if (!phenByPos.has(pos)) phenByPos.set(pos, new Set());
            if (!identifiedLabelsByPos.has(pos)) identifiedLabelsByPos.set(pos, new Set());
            const identifiedLabels = identifiedLabelsByPos.get(pos);
            (v.diseasePairs || []).forEach(p => {
                const key = p.id || p.label;
                if (p.id && p.label) {
                    identifiedLabels.add(p.label);
                    phenByPos.get(pos).delete(p.label);
                }
                if (key) phenByPos.get(pos).add(key);
            });
            // Structured disease IDs are preferred. Description-derived labels are the
            // fallback when the API does not provide a paired identifier.
            (v.diseases || []).forEach(label => {
                if (label && !identifiedLabels.has(label)) phenByPos.get(pos).add(label);
            });
        });
        if (mutCountByPos.size === 0) return new Set();

        const positions = [...mutCountByPos.keys()];
        const phenCountByPos = new Map(positions.map(pos => [pos, phenByPos.get(pos)?.size || 0]));

        // Percentile rank of a value within a value list: fraction of entries ≤ value (max → 1).
        const rankerFor = values => {
            const sorted = [...values].sort((a, b) => a - b);
            const N = sorted.length;
            return v => {
                // count of entries ≤ v via binary search for the upper bound
                let lo = 0, hi = N;
                while (lo < hi) { const mid = (lo + hi) >> 1; if (sorted[mid] <= v) lo = mid + 1; else hi = mid; }
                return lo / N;
            };
        };
        const rankMut = rankerFor([...mutCountByPos.values()]);
        const rankPhen = rankerFor([...phenCountByPos.values()]);

        const COMPOSITE_MIN = 0.90; // top-decile composite rank
        const flagged = new Set();
        positions.forEach(pos => {
            const count = mutCountByPos.get(pos);
            const phenCount = phenCountByPos.get(pos);
            // Absolute floors: real recurrence + pleiotropy (keeps the signal interpretable).
            if (count < 2 || phenCount < 2) return;
            const composite = 0.5 * rankMut(count) + 0.5 * rankPhen(phenCount);
            if (composite >= COMPOSITE_MIN) flagged.add(pos);
        });
        return flagged;
    }

    /**
     * PTM–Variant Proximity analysis.
     * For each PTM site that is modelled in the current structure, finds variant-containing
     * residues within 12 Å Cα distance and assigns a tier:
     *   Tier 1 — variant on the same residue as the PTM
     *   Tier 2 — pathogenic/deleterious variant within 8 Å
     *   Tier 3 — any variant within 8–12 Å (or same residue, non-pathogenic)
     *
     * Returns Map<ptmUniPos, {
     *   tier:          number (1|2|3|null),
     *   tier1:         [{variant}],
     *   tier2:         [{variant, dist}],
     *   tier3:         [{variant, dist}],
     *   nearbyCount8A: number,
     *   pathCount8A:   number,
     *   nearestDist:   number|null,
     *   nearestVariant: string|null,
     * }>
     *
     * geometry: array from StructureViewer.residueGeometry() — {uniPos, chain, resi, ca:{x,y,z}}
     */
    function computePtmVariantProximity(ptms, variants, geometry, d8 = 8, d12 = 12) {
        // Build uniPos → Cα position map (first modelled copy per UniProt position)
        const caByUni = new Map();
        for (const g of geometry) {
            if (g.uniPos != null && !caByUni.has(g.uniPos)) caByUni.set(g.uniPos, g.ca);
        }

        // Index variants by position for fast lookup
        const varsByPos = new Map();
        variants.forEach(v => {
            if (!varsByPos.has(v.position)) varsByPos.set(v.position, []);
            varsByPos.get(v.position).push(v);
        });

        const isPathogenic = v => /pathogenic|deleterious/i.test(v.consequence || v.clinVarSignificance || '');

        const result = new Map();
        const seenPtmPos = new Set();

        ptms.forEach(ptm => {
            const ptmPos = ptm.position;
            if (seenPtmPos.has(ptmPos)) return;
            seenPtmPos.add(ptmPos);
            const ptmCa = caByUni.get(ptmPos);
            if (!ptmCa) return; // not modelled in this structure

            const tier1 = [], tier2 = [], tier3 = [];
            let nearestDist = Infinity, nearestVariant = null;

            // Tier 1 — same residue
            (varsByPos.get(ptmPos) || []).forEach(v => tier1.push(v));

            // Tiers 2 & 3 — spatial neighbours
            caByUni.forEach((ca, uniPos) => {
                if (uniPos === ptmPos) return;
                const dx = ca.x - ptmCa.x, dy = ca.y - ptmCa.y, dz = ca.z - ptmCa.z;
                const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
                if (dist > d12) return;
                const varsHere = varsByPos.get(uniPos) || [];
                if (!varsHere.length) return;
                varsHere.forEach(v => {
                    if (!isPathogenic(v)) return;
                    const label = `${v.wildType || ''}${v.position}${v.mutant || ''}`;
                    (dist <= d8 ? tier2 : tier3).push({ variant: v, dist });
                    if (dist < nearestDist) { nearestDist = dist; nearestVariant = label; }
                });
            });
            // Nearest for tier1 is 0 Å
            if (tier1.length && nearestDist > 0) { nearestDist = 0; nearestVariant = `${tier1[0].wildType || ''}${ptmPos}${tier1[0].mutant || ''}`; }

            const pathCount8A = tier2.length;
            const nearbyCount8A = tier1.length + tier2.length;
            const tier = tier1.length ? 1 : tier2.length ? 2 : tier3.length ? 3 : null;

            if (tier !== null) {
                result.set(ptmPos, {
                    tier, tier1, tier2, tier3,
                    nearbyCount8A, pathCount8A, d8, d12,
                    nearestDist: nearestDist === Infinity ? null : nearestDist,
                    nearestVariant,
                });
            }
        });
        return result;
    }

    /**
     * Per-residue proximity analysis (residue-centric, not PTM-centric).
     * All distances are measured FROM `pos`, not from any PTM site.
     *
     * Returns {
     *   nearbyPtms:      [{ptm, dist}]     — PTMs within ptmRadius Å, sorted by dist
     *   nearbyVariants:  [{variant, dist}]  — pathogenic variants within varRadius Å, sorted by dist
     * }
     */
    function computeResidueProximity(pos, ptms, variants, geometry, ptmRadius = 8, varRadius = 12) {
        const caByUni = new Map();
        for (const g of geometry) {
            if (g.uniPos != null && !caByUni.has(g.uniPos)) caByUni.set(g.uniPos, g.ca);
        }
        const rCa = caByUni.get(pos);
        if (!rCa) return { nearbyPtms: [], nearbyVariants: [] };

        const seenPtmPos = new Set();
        const nearbyPtms = [];
        for (const ptm of ptms) {
            const p = ptm.position;
            if (p === pos) continue; // exclude the selected residue itself (distance 0)
            if (seenPtmPos.has(p)) continue;
            seenPtmPos.add(p);
            const ca = caByUni.get(p);
            if (!ca) continue;
            const dx = ca.x - rCa.x, dy = ca.y - rCa.y, dz = ca.z - rCa.z;
            const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (dist <= ptmRadius) nearbyPtms.push({ ptm, dist });
        }
        nearbyPtms.sort((a, b) => a.dist - b.dist);

        const isPathogenic = v => /pathogenic|deleterious/i.test(v.consequence || v.clinVarSignificance || '');
        const nearbyVariants = [];
        for (const v of variants) {
            if (v.position === pos) continue; // exclude the selected residue itself (distance 0)
            if (!isPathogenic(v)) continue;
            const ca = caByUni.get(v.position);
            if (!ca) continue;
            const dx = ca.x - rCa.x, dy = ca.y - rCa.y, dz = ca.z - rCa.z;
            const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (dist <= varRadius) nearbyVariants.push({ variant: v, dist });
        }
        nearbyVariants.sort((a, b) => a.dist - b.dist);

        return { nearbyPtms, nearbyVariants };
    }

    /**
     * Aggregate pocket-level evidence for a set of pocket residues (UniProt positions).
     *   meanB    — mean CA B-factor over pocket residues (= mean pLDDT for AlphaFold/AlphaFill models)
     *   count    — number of pocket residues that had a CA B-factor
     *   siteHits — UniProt functional sites (ACT_SITE/BINDING/METAL/SITE) overlapping the pocket
     * @param {Iterable<number>} pocketPositions  UniProt positions forming the pocket
     * @param {Array} geometry  residueGeometry() output ({ uniPos, ca, b })
     * @param {Array} sites      UniProt site features ({ position, endPosition, description })
     */
    function pocketEvidence(pocketPositions, geometry, sites, sequence) {
        const posSet = pocketPositions instanceof Set ? pocketPositions : new Set(pocketPositions || []);
        let sum = 0, count = 0;
        for (const g of (geometry || [])) {
            if (g.uniPos != null && posSet.has(g.uniPos) && typeof g.b === 'number') { sum += g.b; count++; }
        }
        const meanB = count ? sum / count : null;
        const siteHits = (sites || []).filter(site => {
            const end = site.endPosition || site.position;
            for (let p = site.position; p <= end; p++) if (posSet.has(p)) return true;
            return false;
        });
        // Pocket-residue fingerprint: amino-acid CLASS composition of the lining residues (from the
        // canonical sequence), a quick read on whether the pocket is hydrophobic, polar or charged.
        const composition = pocketComposition(posSet, sequence);
        return { meanB, count, siteHits, composition };
    }

    const _AA_CLASS = {};
    'AVLIMFWPGC'.split('').forEach(a => (_AA_CLASS[a] = 'hydrophobic'));
    'STNQY'.split('').forEach(a => (_AA_CLASS[a] = 'polar'));
    'DE'.split('').forEach(a => (_AA_CLASS[a] = 'acidic'));
    'KRH'.split('').forEach(a => (_AA_CLASS[a] = 'basic'));
    function pocketComposition(posSet, sequence) {
        if (!sequence) return null;
        const comp = { hydrophobic: 0, polar: 0, acidic: 0, basic: 0, aromatic: 0, total: 0 };
        posSet.forEach(p => {
            const aa = sequence[p - 1]; if (!aa) return;
            comp.total++;
            const c = _AA_CLASS[aa]; if (c) comp[c]++;
            if ('FWYH'.includes(aa)) comp.aromatic++;
        });
        return comp.total ? comp : null;
    }

    /**
     * Closest heavy-atom contact between a ligand and the protein — a pure-JS clash check (works for any
     * ligand, not just AlphaFill transplants). Returns the min distance and a band: severe (<2.0 Å,
     * atoms interpenetrate), suspicious (<2.6 Å, tighter than a typical non-bonded contact), or normal.
     * @param {Array} atoms  full atom list (getModel().selectedAtoms({}))
     */
    function ligandClash(ligResn, ligResi, ligChain, atoms) {
        if (!atoms || !atoms.length) return null;
        const isHeavy = a => a.elem !== 'H' && a.elem !== 'D';
        const lig = atoms.filter(a => a.hetflag && a.resn === ligResn && a.resi === ligResi && (ligChain == null || a.chain === ligChain) && isHeavy(a));
        if (!lig.length) return null;
        const prot = atoms.filter(a => !a.hetflag && isHeavy(a));
        let min2 = Infinity;
        for (const l of lig) for (const p of prot) {
            const dx = l.x - p.x, dy = l.y - p.y, dz = l.z - p.z;
            const d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < min2) min2 = d2;
        }
        if (!isFinite(min2)) return null;
        const d = Math.sqrt(min2);
        return { minDist: d, band: d < 2.0 ? 'severe' : d < 2.6 ? 'suspicious' : 'normal' };
    }

    return { residueNeighborhood, computeHotspots, computeDistantContacts, aggregateAlphaMissense, computeResidueBurden, computePtmVariantProximity, computeResidueProximity, pocketEvidence, ligandClash };
})();
