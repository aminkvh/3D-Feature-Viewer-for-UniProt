/* ================================================================================
   Constraint-pocket residue prioritization
   --------------------------------------------------------------------------------
   Highlights buried, evolutionarily-constrained residues that cluster together in 3-D —
   candidate catalytic or binding-pocket regions.  For each residue the mean AlphaMissense
   pathogenicity is taken as a constraint signal and regressed against structural burial
   (coordination number) with a local LOESS fit; the residual is the constraint that is NOT
   explained by burial.  A Getis-Ord Gi* spatial-autocorrelation statistic on these residuals,
   assessed with a within-protein permutation null and Benjamini-Hochberg FDR control, flags
   residues whose neighbourhoods are unexpectedly constrained for their depth.

   Runs in pure JavaScript on data already in memory:
     • AlphaMissense substitution matrix (amMap)   — constraint signal
     • Cα coordinates from the loaded model         — burial + spatial weights
     • AlphaFold Predicted Aligned Error (optional) — gates the spatial weights
   ================================================================================ */
const UFVPocket = (() => {
    'use strict';

    const AM_AAS = 'ACDEFGHIKLMNPQRSTVWY';
    const AM_PATH = 0.564;     // AlphaMissense "likely pathogenic" threshold
    const HSE_RADIUS = 13;     // Å — half-sphere-exposure / coordination sphere radius
    const WEIGHT_CUTOFF = 13;  // Å — spatial-weight neighbourhood radius
    const SIGMA = 5;           // Å — Gaussian spatial-weight bandwidth
    const SEED = 0x9e3779b9;
    const MAX_N = 3000;        // residue cap to keep the in-browser compute responsive

    // Permutation count scales down with size to bound compute time; the FDR threshold
    // stays resolvable (min p = 1/(PERM+1) ≤ 0.0025 even at the smallest count).
    function permCount(n) { return n <= 600 ? 999 : n <= 1200 ? 599 : 399; }

    // Deterministic PRNG (mulberry32) so the permutation null is reproducible run-to-run.
    function mulberry32(seed) {
        let a = seed >>> 0;
        return function () {
            a |= 0; a = (a + 0x6D2B79F5) | 0;
            let t = Math.imul(a ^ (a >>> 15), 1 | a);
            t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
            return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
        };
    }

    function shuffleInPlace(arr, rand) {
        for (let i = arr.length - 1; i > 0; i--) {
            const r = Math.floor(rand() * (i + 1));
            const t = arr[i]; arr[i] = arr[r]; arr[r] = t;
        }
    }

    // Per-position AlphaMissense mean pathogenicity (over the 19 possible substitutions).
    function amMean(uniPos, wt, amMap) {
        if (!wt) return null;
        let sum = 0, cnt = 0;
        for (const mut of AM_AAS) {
            if (mut === wt) continue;
            const sc = amMap.get(`${wt}${uniPos}${mut}`);
            if (Number.isFinite(sc)) { sum += sc; cnt++; }
        }
        return cnt === 0 ? null : sum / cnt;
    }

    /**
     * Build per-residue geometry, burial and the spatial-weight graph.
     *
     * Multi-chain handling (matters for complexes such as the GABA-A receptor):
     *   • BURIAL CONTEXT = every modelled Cα in the file — all copies of the entry protein AND
     *     the partner-protein chains around it — so a residue buried at a subunit interface is
     *     correctly seen as buried.  Geometry records with uniPos == null are partner-protein (or
     *     unmapped) atoms; they contribute to burial but are never scored or reported.
     *   • SCORED RESIDUES = the entry protein's residues that carry AlphaMissense scores.  Every
     *     physical copy is scored (the copies of a homo-oligomer sit in different environments)
     *     and collapsed back to one call per UniProt position downstream.
     */
    function prepare(geometry, amMap, sequence, pae, maxN) {
        if (!amMap || amMap.size === 0) return { reason: 'No AlphaMissense data for this protein.' };
        if (!geometry || geometry.length === 0) return { reason: 'No modelled residues.' };

        const nCtx = geometry.length;
        const cx = new Float64Array(nCtx), cy = new Float64Array(nCtx), cz = new Float64Array(nCtx);
        for (let i = 0; i < nCtx; i++) { cx[i] = geometry[i].ca.x; cy[i] = geometry[i].ca.y; cz[i] = geometry[i].ca.z; }

        // `ctxIdx` ties each scored residue to its atom in the context arrays, so burial can
        // count every nearby atom while skipping the residue's own atom.
        const recs = [];
        for (let gi = 0; gi < nCtx; gi++) {
            const g = geometry[gi];
            if (g.uniPos == null) continue;                  // partner/unmapped → context only
            const wt = sequence ? sequence[g.uniPos - 1] : null;
            const am = amMean(g.uniPos, wt, amMap);
            if (am == null) continue;
            recs.push({ uniPos: g.uniPos, chain: g.chain, resi: g.resi, ca: g.ca, am, ctxIdx: gi });
        }
        if (recs.length < 12) return { reason: 'Too few residues with AlphaMissense scores.' };
        if (recs.length > maxN) return { reason: `Protein too large for in-browser analysis (${recs.length} > ${maxN} residues).` };

        const n = recs.length;
        const xs = new Float64Array(n), ys = new Float64Array(n), zs = new Float64Array(n);
        for (let i = 0; i < n; i++) { xs[i] = recs[i].ca.x; ys[i] = recs[i].ca.y; zs[i] = recs[i].ca.z; }

        // Sequential Cα-neighbour lookup per chain for the half-sphere pseudo-Cβ direction.
        const byChainResi = new Map(); // chain → (resi → scored index)
        recs.forEach((r, i) => {
            if (!byChainResi.has(r.chain)) byChainResi.set(r.chain, new Map());
            byChainResi.get(r.chain).set(r.resi, i);
        });
        const norm = v => { const L = Math.hypot(v.x, v.y, v.z) || 1; return { x: v.x / L, y: v.y / L, z: v.z / L }; };
        const dirs = new Array(n).fill(null);
        for (let i = 0; i < n; i++) {
            const cm = byChainResi.get(recs[i].chain);
            const prev = cm.get(recs[i].resi - 1);
            const next = cm.get(recs[i].resi + 1);
            // Direction pointing away from sequential neighbours ≈ Cα→Cβ (side-chain side).
            if (prev != null && next != null) dirs[i] = norm({ x: 2 * xs[i] - xs[prev] - xs[next], y: 2 * ys[i] - ys[prev] - ys[next], z: 2 * zs[i] - zs[prev] - zs[next] });
            else if (prev != null) dirs[i] = norm({ x: xs[i] - xs[prev], y: ys[i] - ys[prev], z: zs[i] - zs[prev] });
            else if (next != null) dirs[i] = norm({ x: xs[i] - xs[next], y: ys[i] - ys[next], z: zs[i] - zs[next] });
        }

        // PAE gate e^(-PAE/τ).  Index the matrix by the residue's absolute model-order index
        // (ctxIdx), not its UniProt position — uniPos repeats across the chains of a homo-
        // oligomer.  For a monomer ctxIdx == uniPos-1; an inter-chain lookup outside a monomeric
        // PAE matrix returns gate = 1 (pure Euclidean weight).
        const TAU = 10;
        const paeGate = (a, b) => {
            if (!pae) return 1;
            const i = recs[a].ctxIdx, j = recs[b].ctxIdx;
            if (i >= pae.n || j >= pae.n) return 1;
            const e = 0.5 * (pae.data[i * pae.n + j] + pae.data[j * pae.n + i]);
            return Math.exp(-e / TAU);
        };

        const twoSigma2 = 2 * SIGMA * SIGMA;
        const cut2 = WEIGHT_CUTOFF * WEIGHT_CUTOFF, hse2 = HSE_RADIUS * HSE_RADIUS;

        // Burial (coordination number, upward half-sphere count) of each scored residue against
        // the FULL context, so interface burial by partner subunits is counted. O(n·nCtx).
        const hseUp = new Float64Array(n), hseDown = new Float64Array(n), cn = new Float64Array(n);
        for (let i = 0; i < n; i++) {
            const ci = recs[i].ctxIdx, xi = xs[i], yi = ys[i], zi = zs[i], di = dirs[i];
            for (let k = 0; k < nCtx; k++) {
                if (k === ci) continue;
                const dx = cx[k] - xi, dy = cy[k] - yi, dz = cz[k] - zi;
                if (dx * dx + dy * dy + dz * dz <= hse2) {
                    cn[i]++;
                    if (di) (dx * di.x + dy * di.y + dz * di.z) > 0 ? hseUp[i]++ : hseDown[i]++;
                }
            }
        }

        // Gaussian, PAE-gated spatial-weight graph among scored residues. O(n²).
        const neigh = Array.from({ length: n }, () => []);
        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                const dx = xs[i] - xs[j], dy = ys[i] - ys[j], dz = zs[i] - zs[j];
                const d2 = dx * dx + dy * dy + dz * dz;
                if (d2 <= cut2) {
                    const w = Math.exp(-d2 / twoSigma2) * paeGate(i, j);
                    if (w > 1e-6) { neigh[i].push({ j, w }); neigh[j].push({ j: i, w }); }
                }
            }
        }

        return { recs, n, neigh, hseUp, hseDown, cn, hasPae: !!pae };
    }

    // 1-D LOESS (local weighted linear regression, tricube kernel) of Y over X.
    // Returns the smoothed prediction at each X_i; bandwidth = fraction f of points.
    function loess(X, Y, f = 0.3) {
        const n = X.length;
        const k = Math.max(2, Math.min(n, Math.round(f * n)));
        const order = [...X.keys()].sort((a, b) => X[a] - X[b]);
        const pred = new Float64Array(n);
        for (let oi = 0; oi < n; oi++) {
            const i = order[oi];
            const dists = order.map(j => Math.abs(X[j] - X[i])).sort((a, b) => a - b);
            const h = dists[Math.min(k - 1, n - 1)] || 1e-9;
            let sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
            for (let j = 0; j < n; j++) {
                const u = Math.abs(X[j] - X[i]) / h;
                if (u >= 1) continue;
                const w = Math.pow(1 - u * u * u, 3);        // tricube
                sw += w; swx += w * X[j]; swy += w * Y[j]; swxx += w * X[j] * X[j]; swxy += w * X[j] * Y[j];
            }
            const denom = sw * swxx - swx * swx;
            if (Math.abs(denom) < 1e-12) { pred[i] = sw ? swy / sw : Y[i]; continue; }
            const b = (sw * swxy - swx * swy) / denom;
            pred[i] = (swy - b * swx) / sw + b * X[i];
        }
        return pred;
    }

    // Geometric label for a flagged residue: buried/concave → 'pocket', otherwise 'exposed'.
    // A visualization aid, not a validated site-type call.
    function categorise(hseUp, hseDown, cn, medCn) {
        const concavity = hseUp / (hseDown + 1);
        return (cn >= medCn || concavity >= 1.3) ? 'pocket' : 'exposed';
    }

    /**
     * Constraint pockets: Getis-Ord Gi* on the residual pathogenicity (AlphaMissense mean minus
     * its LOESS fit against coordination number).  Returns byPos = Map(uniPos → {cat, p, q,
     * score}) with EVERY positive-residual candidate and its BH-FDR q-value; the consumer applies
     * the display threshold so a UI control can re-threshold without recomputing.
     */
    function computePockets(geometry, amMap, sequence, pae = null, options = {}) {
        const prep = prepare(geometry, amMap, sequence, pae, options.maxN || MAX_N);
        if (prep.reason) return { byPos: new Map(), reason: prep.reason };
        const { recs, n, neigh, hseUp, hseDown, cn, hasPae } = prep;

        // Residual pathogenicity vs a per-protein burial (coordination-number) baseline.
        const am = recs.map(r => r.am);
        const cnArr = Array.from(cn);
        const expected = loess(cnArr, am, 0.3);
        const rp = am.map((v, i) => v - expected[i]);

        // Getis-Ord Gi*_i = Σ_j w_ij x_j (self included via a w_ii ≈ 1 term).
        const giObs = new Float64Array(n);
        for (let i = 0; i < n; i++) {
            let acc = rp[i];
            for (const { j, w } of neigh[i]) acc += w * rp[j];
            giObs[i] = acc;
        }

        // Permutation null: shuffle residuals, recompute Gi*, count one-sided exceedances.
        const PERM = permCount(n);
        const ge = new Int32Array(n);
        const rpPerm = rp.slice();
        const rand = mulberry32(SEED);
        for (let p = 0; p < PERM; p++) {
            shuffleInPlace(rpPerm, rand);
            for (let i = 0; i < n; i++) {
                let acc = rpPerm[i];
                for (const { j, w } of neigh[i]) acc += w * rpPerm[j];
                if (acc >= giObs[i]) ge[i]++;
            }
        }

        const cnSorted = cnArr.slice().sort((a, b) => a - b);
        const medCn = cnSorted[cnSorted.length >> 1] || 0;

        // Candidates = positive residual (unexpectedly constrained).  Collapse homo-oligomer
        // copies to the most-significant per UniProt position, then BH-FDR over unique positions.
        const bestByPos = new Map(); // uniPos → { i, p }
        for (let i = 0; i < n; i++) {
            if (rp[i] <= 0) continue;
            const p = (ge[i] + 1) / (PERM + 1);
            const u = recs[i].uniPos;
            const ex = bestByPos.get(u);
            if (!ex || p < ex.p) bestByPos.set(u, { i, p });
        }
        const cands = [...bestByPos.values()];
        cands.sort((a, b) => a.p - b.p);
        const Nc = cands.length;
        cands.forEach((c, k) => { c.q = Math.min(c.p * Nc / (k + 1), 1); });
        for (let k = Nc - 2; k >= 0; k--) cands[k].q = Math.min(cands[k].q, cands[k + 1].q);

        const byPos = new Map();
        cands.forEach(c => byPos.set(recs[c.i].uniPos, {
            cat: categorise(hseUp[c.i], hseDown[c.i], cn[c.i], medCn), p: c.p, q: c.q, score: rp[c.i],
        }));
        return { byPos, reason: null, hasPae, n };
    }

    return { computePockets };
})();

// When loaded inside a Web Worker (no document) via importScripts(), expose on the worker
// global so the worker entry script can reach it; in the content-script context this is skipped.
if (typeof self !== 'undefined' && typeof document === 'undefined') self.UFVPocket = UFVPocket;
