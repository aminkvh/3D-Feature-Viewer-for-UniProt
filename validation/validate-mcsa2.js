/* ============================================================================
   M-CSA extended characterization — three analyses on the same enzyme set:
     1. Baseline comparison: raw AM (no burial normalization) vs PRISM
     2. Enrichment distribution (percentiles, fraction beating chance)
     3. Per-protein predictors of success (size, catalytic-site count, flag rate)
   Run:  node validation/validate-mcsa2.js [N]   (N = target, default 40)
   ============================================================================ */
const fs = require('fs');
const path = require('path');
const code = fs.readFileSync(path.join(__dirname, '..', 'algorithms.js'), 'utf8');
let UFVPocket; eval(code.replace('const UFVPocket', 'UFVPocket'));

const TARGET = parseInt(process.argv[2]) || 40;
const Q = 0.10;

async function txt(u) { const r = await fetch(u); if (!r.ok) throw new Error(r.status); return r.text(); }
async function jsn(u) { const r = await fetch(u, { headers: { Accept: 'application/json' } }); if (!r.ok) throw new Error(r.status); return r.json(); }
async function head(u) { try { return (await fetch(u, { method: 'HEAD' })).ok; } catch (_) { return false; } }
const dist = (a, b) => Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);

function parseAmCsv(t) { const m = new Map(); for (const l of t.split('\n').slice(1)) { const s = l.trim(); if (!s) continue; const c = s.indexOf(','); if (c < 0) continue; const v = s.slice(0, c), sc = Number(s.slice(c + 1).split(',')[0]); if (v && Number.isFinite(sc)) m.set(v, sc); } return m; }
function seqFromAm(am) { const wt = new Map(); for (const k of am.keys()) { const m = k.match(/^([A-Z])(\d+)[A-Z]$/); if (m) wt.set(+m[2], m[1]); } const L = Math.max(...wt.keys()); let s = ''; for (let i = 1; i <= L; i++) s += wt.get(i) || 'X'; return s; }
function pdbCA(p) { const o = [], seen = new Set(); for (const l of p.split('\n')) { if (!l.startsWith('ATOM') || l.slice(12, 16).trim() !== 'CA') continue; const resi = parseInt(l.slice(22, 26)); if (seen.has(resi)) continue; seen.add(resi); o.push({ resi, ca: { x: +l.slice(30, 38), y: +l.slice(38, 46), z: +l.slice(46, 54) } }); } return o; }

// ---- Raw-AM baseline: Gi* on raw pathogenic fraction, no burial regression ----
// This is the control condition — the same spatial clustering applied to raw AM scores
// rather than LOESS residuals. Measures whether the burial normalization adds any value.
function computeRawAM(geometry, amMap, sequence) {
    const AAs = 'ACDEFGHIKLMNPQRSTVWY';
    const AM_PATH = 0.564;
    const WEIGHT_CUTOFF = 13, SIGMA = 5, twoSigma2 = 2 * SIGMA * SIGMA;
    const cut2 = WEIGHT_CUTOFF * WEIGHT_CUTOFF;
    const SEED = 0x9e3779b9;
    function mulberry32(seed) { let a = seed >>> 0; return () => { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
    function shuffle(arr, rand) { for (let i = arr.length - 1; i > 0; i--) { const r = Math.floor(rand() * (i + 1)); [arr[i], arr[r]] = [arr[r], arr[i]]; } }

    // compute raw PF per residue
    const recs = [];
    for (const g of geometry) {
        const wt = sequence[g.uniPos - 1];
        if (!wt) continue;
        let cnt = 0, path = 0;
        for (const m of AAs) { if (m === wt) continue; const sc = amMap.get(`${wt}${g.uniPos}${m}`); if (Number.isFinite(sc)) { cnt++; if (sc > AM_PATH) path++; } }
        if (cnt === 0) continue;
        recs.push({ uniPos: g.uniPos, pf: path / cnt, ca: g.ca });
    }
    if (recs.length < 12) return { byPos: new Map(), reason: 'too few' };
    const n = recs.length;
    const xs = recs.map(r => r.ca.x), ys = recs.map(r => r.ca.y), zs = recs.map(r => r.ca.z);

    // spatial weights (Gaussian)
    const neigh = Array.from({ length: n }, () => []);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
        const dx = xs[i]-xs[j], dy = ys[i]-ys[j], dz = zs[i]-zs[j], d2 = dx*dx+dy*dy+dz*dz;
        if (d2 <= cut2) { const w = Math.exp(-d2 / twoSigma2); if (w > 1e-6) { neigh[i].push({j,w}); neigh[j].push({j:i,w}); } }
    }

    const pf = recs.map(r => r.pf);
    const giObs = new Float64Array(n);
    for (let i = 0; i < n; i++) { let acc = pf[i]; for (const {j,w} of neigh[i]) acc += w * pf[j]; giObs[i] = acc; }

    const PERM = n <= 600 ? 999 : n <= 1200 ? 599 : 399;
    const ge = new Int32Array(n);
    const pfP = pf.slice();
    const rand = mulberry32(SEED);
    for (let p = 0; p < PERM; p++) {
        shuffle(pfP, rand);
        for (let i = 0; i < n; i++) { let acc = pfP[i]; for (const {j,w} of neigh[i]) acc += w * pfP[j]; if (acc >= giObs[i]) ge[i]++; }
    }

    // BH-FDR same as PRISM
    const cands = recs.map((r,i) => ({ i, p: (ge[i]+1)/(PERM+1), uniPos: r.uniPos })).sort((a,b) => a.p-b.p);
    cands.forEach((c,k) => { c.q = Math.min(c.p * cands.length / (k+1), 1); });
    for (let k = cands.length-2; k >= 0; k--) cands[k].q = Math.min(cands[k].q, cands[k+1].q);
    const byPos = new Map();
    cands.forEach(c => byPos.set(c.uniPos, { p: c.p, q: c.q }));
    return { byPos, reason: null };
}

async function gatherEnzymes(target) {
    const out = [];
    let url = 'https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&page_size=100';
    let scanned = 0;
    while (url && out.length < target && scanned < 600) {
        const page = await jsn(url);
        for (const e of page.results) {
            scanned++;
            const id = e.reference_uniprot_id; if (!id) continue;
            const cat = new Set();
            (e.residues||[]).forEach(r => (r.residue_sequences||[]).forEach(rs => { if (rs.uniprot_id===id && Number.isFinite(rs.resid)) cat.add(rs.resid); }));
            if (cat.size < 2) continue;
            if (out.some(o => o.id === id)) continue;
            if (!await head(`https://alphafold.ebi.ac.uk/files/AF-${id}-F1-aa-substitutions.csv`)) continue;
            out.push({ id, name: e.enzyme_name, cat, ec: e.all_ecs?.[0] || '' });
            if (out.length >= target) break;
            process.stdout.write(`\r  gathering… ${out.length}/${target} (scanned ${scanned})`);
        }
        url = page.next;
    }
    process.stdout.write('\n');
    return out;
}

const recallOf = (flags, gt, caByPos, r=5) => { if (!gt.size) return null; let h=0; for (const t of gt) { const tc=caByPos.get(t); if(!tc) continue; for (const f of flags) { const fc=caByPos.get(f); if(fc && dist(tc,fc)<=r){h++;break;} } } return h/gt.size; };
let rs = 42; const rng = () => { rs|=0; rs=(rs+0x6D2B79F5)|0; let t=Math.imul(rs^(rs>>>15),1|rs); t=(t+Math.imul(t^(t>>>7),61|t))^t; return((t^(t>>>14))>>>0)/4294967296; };
const expRecall = (k,gt,caByPos,allPos) => { if(!gt.size||!k) return 0; let s=0; for(let d=0;d<100;d++){const pool=allPos.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(pool.length-i));[pool[i],pool[j]]=[pool[j],pool[i]];}s+=recallOf(pool.slice(0,k),gt,caByPos);} return s/100; };
const flagsAt = (res, q) => [...res.byPos.entries()].filter(([,v]) => v.q <= q).map(([k]) => k);

const pct = (n, d, dec=1) => d ? (n/d*100).toFixed(dec)+'%' : 'n/a';
const pctile = (arr, p) => { const s = arr.slice().sort((a,b)=>a-b); return s[Math.floor(p*(s.length-1))]; };

(async () => {
    console.log(`M-CSA extended characterization — targeting ${TARGET} enzymes\n`);
    const enzymes = await gatherEnzymes(TARGET);
    console.log(`\nRunning ${enzymes.length} enzymes…\n`);

    // per-enzyme results
    const rows = [];
    let ok = 0;

    for (const e of enzymes) {
        try {
            const [amT, pdb] = await Promise.all([
                txt(`https://alphafold.ebi.ac.uk/files/AF-${e.id}-F1-aa-substitutions.csv`),
                txt(`https://alphafold.ebi.ac.uk/files/AF-${e.id}-F1-model_v6.pdb`).catch(() => txt(`https://alphafold.ebi.ac.uk/files/AF-${e.id}-F1-model_v4.pdb`)),
            ]);
            const am = parseAmCsv(amT), seq = seqFromAm(am), cas = pdbCA(pdb);
            if (cas.length > 1400) continue;
            const caByPos = new Map(cas.map(c => [c.resi, c.ca]));
            const geometry = cas.map(c => ({ uniPos: c.resi, chain: 'A', resi: c.resi, ca: c.ca }));
            const gt = new Set([...e.cat].filter(p => caByPos.has(p)));
            if (gt.size < 2) continue;
            const allPos = geometry.map(g => g.uniPos);
            const n = geometry.length;

            const P = UFVPocket.computePockets(geometry, am, seq, null);
            const R = computeRawAM(geometry, am, seq);
            if (P.reason) continue;

            const pf = flagsAt(P, Q), rf = flagsAt(R, Q);
            const pRecall = recallOf(pf, gt, caByPos);
            const rRecall = recallOf(rf, gt, caByPos);
            const pExp = expRecall(pf.length, gt, caByPos, allPos);
            const rExp = expRecall(rf.length, gt, caByPos, allPos);
            const pEnr = Number.isFinite(pRecall) ? pRecall / (pExp || 1) : null;
            const rEnr = Number.isFinite(rRecall) ? rRecall / (rExp || 1) : null;

            rows.push({
                id: e.id, name: e.name, ec: e.ec,
                n,                          // protein length
                nCat: gt.size,              // catalytic residue count
                flagRate: pf.length / n,    // PRISM flag rate
                prismEnr: pEnr,             // PRISM enrichment
                rawEnr: rEnr,               // raw-AM enrichment
                prismRecall: pRecall,
                rawRecall: rRecall,
            });
            ok++;
            process.stdout.write(`\r  ran ${ok} enzymes…`);
        } catch (_) { continue; }
    }
    process.stdout.write('\n\n');

    // ---- Analysis 1: Baseline comparison ----
    const validPrism = rows.filter(r => r.prismEnr != null);
    const validRaw   = rows.filter(r => r.rawEnr   != null);
    const prismEnrs  = validPrism.map(r => r.prismEnr);
    const rawEnrs    = validRaw.map(r => r.rawEnr);
    const prismBeat  = prismEnrs.filter(e => e > 1).length;
    const rawBeat    = rawEnrs.filter(e => e > 1).length;
    console.log('=== 1. BASELINE: PRISM vs Raw AlphaMissense (no burial normalization) ===\n');
    console.log(`  Enzymes with valid enrichment: PRISM ${validPrism.length}  Raw-AM ${validRaw.length}`);
    console.log(`  Median enrichment:             PRISM ${pctile(prismEnrs,0.5).toFixed(2)}×   Raw-AM ${pctile(rawEnrs,0.5).toFixed(2)}×`);
    console.log(`  Mean enrichment:               PRISM ${(prismEnrs.reduce((a,b)=>a+b,0)/prismEnrs.length).toFixed(2)}×   Raw-AM ${(rawEnrs.reduce((a,b)=>a+b,0)/rawEnrs.length).toFixed(2)}×`);
    console.log(`  Fraction beating chance (>1×): PRISM ${pct(prismBeat,prismEnrs.length)}   Raw-AM ${pct(rawBeat,rawEnrs.length)}`);
    console.log(`  Verdict: ${pctile(prismEnrs,0.5) > pctile(rawEnrs,0.5) + 0.1 ? 'LOESS normalization HELPS' : Math.abs(pctile(prismEnrs,0.5) - pctile(rawEnrs,0.5)) <= 0.1 ? 'normalization makes little difference' : 'Raw AM is better (unexpected)'}`);

    // ---- Analysis 2: Enrichment distribution ----
    console.log('\n=== 2. ENRICHMENT DISTRIBUTION (PRISM, q≤0.10) ===\n');
    console.log('  Percentiles (enrichment × chance):');
    for (const p of [10,25,50,75,90]) console.log(`    ${p}th: ${pctile(prismEnrs, p/100).toFixed(2)}×`);
    const hist = [0,0,0,0,0]; // <0.5, 0.5-1, 1-2, 2-4, >4
    prismEnrs.forEach(e => { if(e<0.5) hist[0]++; else if(e<1) hist[1]++; else if(e<2) hist[2]++; else if(e<4) hist[3]++; else hist[4]++; });
    console.log('\n  Distribution of per-enzyme enrichment:');
    const labels = ['<0.5× (harmful)', '0.5–1× (below chance)', '1–2× (weak)', '2–4× (moderate)', '>4× (strong)'];
    hist.forEach((n,i) => console.log(`    ${labels[i].padEnd(24)}: ${'█'.repeat(n)} ${n} (${pct(n,prismEnrs.length)})`));
    console.log(`\n  Fraction with >1× (beats chance): ${pct(prismBeat,prismEnrs.length)}`);
    console.log(`  Fraction with >2×: ${pct(prismEnrs.filter(e=>e>2).length,prismEnrs.length)}`);

    // ---- Analysis 3: Predictors of success ----
    console.log('\n=== 3. PREDICTORS OF SUCCESS ===\n');
    // split by median protein length
    const medLen = pctile(rows.map(r=>r.n), 0.5);
    const short = rows.filter(r => r.n <= medLen && r.prismEnr != null);
    const long  = rows.filter(r => r.n >  medLen && r.prismEnr != null);
    const medEnr = arr => pctile(arr.map(r=>r.prismEnr), 0.5).toFixed(2);
    console.log(`  Protein length split at median (${medLen} residues):`);
    console.log(`    Short (≤${medLen}): n=${short.length}, median enrichment ${medEnr(short)}×`);
    console.log(`    Long  (>${medLen}): n=${long.length},  median enrichment ${medEnr(long)}×`);
    // split by catalytic site count
    const fewCat = rows.filter(r => r.nCat <= 3 && r.prismEnr != null);
    const manyCat = rows.filter(r => r.nCat > 3 && r.prismEnr != null);
    console.log(`\n  Catalytic site count split at 3:`);
    console.log(`    Few  (≤3): n=${fewCat.length}, median enrichment ${fewCat.length?medEnr(fewCat):'n/a'}×`);
    console.log(`    Many (>3): n=${manyCat.length}, median enrichment ${manyCat.length?medEnr(manyCat):'n/a'}×`);
    // flag rate correlation (split at median flag rate)
    const medFlag = pctile(rows.map(r=>r.flagRate), 0.5);
    const lowFlag  = rows.filter(r => r.flagRate <= medFlag && r.prismEnr != null);
    const highFlag = rows.filter(r => r.flagRate >  medFlag && r.prismEnr != null);
    console.log(`\n  Flag-rate split at median (${(medFlag*100).toFixed(0)}%):`);
    console.log(`    Low flag rate:  n=${lowFlag.length}, median enrichment ${lowFlag.length?medEnr(lowFlag):'n/a'}×`);
    console.log(`    High flag rate: n=${highFlag.length}, median enrichment ${highFlag.length?medEnr(highFlag):'n/a'}×`);

    // top 5 and bottom 5
    const sorted = rows.filter(r=>r.prismEnr!=null).sort((a,b)=>b.prismEnr-a.prismEnr);
    console.log('\n  Best 5:');
    sorted.slice(0,5).forEach(r=>console.log(`    ${r.id} ${r.name.slice(0,30).padEnd(30)} n=${r.n} cat=${r.nCat} flag=${pct(r.flagRate,1,0)} enr=${r.prismEnr.toFixed(2)}× vs raw=${r.rawEnr?.toFixed(2)||'?'}×`));
    console.log('  Worst 5:');
    sorted.slice(-5).reverse().forEach(r=>console.log(`    ${r.id} ${r.name.slice(0,30).padEnd(30)} n=${r.n} cat=${r.nCat} flag=${pct(r.flagRate,1,0)} enr=${r.prismEnr.toFixed(2)}× vs raw=${r.rawEnr?.toFixed(2)||'?'}×`));
})();
