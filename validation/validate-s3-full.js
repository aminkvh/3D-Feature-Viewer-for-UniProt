/* ============================================================================
   S3 (mutation/phenotype burden) — full DMS validation.
   Aggregates the COMPLETE Findlay 2018 BRCA1 saturation-genome-editing map (all 13
   exon scoresets, 2 replicates each) for proper coverage, and re-runs TP53
   (Giacomelli 2018, whole-protein). Each scoreset is centred on its own median to
   align scales before pooling per-position medians. Burden-positive positions are
   tested for larger experimental |effect| than other variant-containing positions.
   Run:  node validation/validate-s3-full.js
   ============================================================================ */
const fs = require('fs'), path = require('path');
const dpCode = fs.readFileSync(path.join(__dirname, '..', 'data.js'), 'utf8');
let DataProcessor; eval(dpCode.replace('const DataProcessor', 'DataProcessor'));
const anCode = fs.readFileSync(path.join(__dirname, '..', 'analysis.js'), 'utf8');
let UFVAnalysis; eval(anCode.replace('const UFVAnalysis', 'UFVAnalysis'));

async function jsn(u) { const r = await fetch(u, { headers: { Accept: 'application/json' } }); if (!r.ok) throw new Error(u + ' -> ' + r.status); return r.json(); }
const median = a => { if (!a.length) return NaN; const s = a.slice().sort((x, y) => x - y); return s[s.length >> 1]; };
const parseHgvsPro = h => { if (!h || h === 'NA') return null; const m = h.match(/p\.[A-Za-z*]{1,3}(\d+)/); return m ? parseInt(m[1]) : null; };

// Per-scoreset: per-position median of missense scores, then centre on the scoreset median.
async function scoresetPositionMedians(urn) {
    const r = await fetch(`https://api.mavedb.org/api/v1/score-sets/${urn}/scores`);
    if (!r.ok) return null;
    const lines = (await r.text()).split(/\r?\n/);
    const header = lines[0].split(',').map(h => h.trim());
    const si = header.indexOf('score'), hi = header.indexOf('hgvs_pro');
    if (si < 0 || hi < 0) return null;
    const byPos = new Map();
    for (const line of lines.slice(1)) {
        if (!line.trim()) continue;
        const cols = line.split(',');
        const hgvs = (cols[hi] || '').trim();
        if (/=\s*$/.test(hgvs)) continue;                       // skip synonymous
        const score = parseFloat((cols[si] || '').trim());
        const pos = parseHgvsPro(hgvs);
        if (pos == null || !Number.isFinite(score)) continue;
        if (!byPos.has(pos)) byPos.set(pos, []);
        byPos.get(pos).push(score);
    }
    const med = new Map();
    byPos.forEach((scores, pos) => med.set(pos, median(scores)));
    // centre on this scoreset's own median so exon-to-exon scale offsets don't bias pooling
    const centre = median([...med.values()]);
    const out = new Map();
    med.forEach((v, pos) => out.set(pos, v - centre));
    return out;
}

// Pool position medians across many scoresets → one centred median per position.
async function aggregateDms(urns) {
    const pooled = new Map(); // pos -> [centred medians]
    for (const urn of urns) {
        const m = await scoresetPositionMedians(urn);
        if (!m) continue;
        m.forEach((v, pos) => { if (!pooled.has(pos)) pooled.set(pos, []); pooled.get(pos).push(v); });
        process.stdout.write(`\r  fetched ${urn.split(':')[1]}…   `);
    }
    process.stdout.write('\n');
    const out = new Map();
    pooled.forEach((vals, pos) => out.set(pos, median(vals)));
    return out;
}

async function runS3(id, label, urns, maveLabel) {
    console.log(`\n--- S3: ${label} (${id}) | DMS: ${maveLabel} (${urns.length} scoreset${urns.length > 1 ? 's' : ''}) ---`);
    const variation = await jsn(`https://www.ebi.ac.uk/proteins/api/variation/${id}`);
    const variants = DataProcessor.extractVariants(variation);
    const burden = UFVAnalysis.computeResidueBurden(variants);
    const dms = await aggregateDms(urns);
    if (burden.size === 0) { console.log('  no burden-positive residues'); return; }

    const burdenScores = [], nonScores = [];
    dms.forEach((score, pos) => { (burden.has(pos) ? burdenScores : nonScores).push(score); });
    const pct = (n, d) => d ? (n / d * 100).toFixed(1) + '%' : 'n/a';
    console.log(`  Burden-positive positions: ${burden.size} | with DMS coverage: ${burdenScores.length} | non-burden w/ DMS: ${nonScores.length}`);
    if (burdenScores.length === 0) { console.log('  no overlap between burden positions and DMS coverage'); return; }

    // Direction-agnostic effect: |score − neutral|, neutral = pooled median.
    const all = [...dms.values()];
    const neutral = median(all);
    const eff = arr => arr.map(s => Math.abs(s - neutral));
    const effB = eff(burdenScores), effN = eff(nonScores);
    console.log(`  DMS |effect| from neutral (median): burden+ ${median(effB).toFixed(3)}   non-burden ${median(effN).toFixed(3)}  → burden+ ${median(effB) > median(effN) ? 'LARGER (expected)' : 'not larger'}`);
    const sortedAbs = all.map(s => Math.abs(s - neutral)).sort((a, b) => b - a);
    const top20 = sortedAbs[Math.floor(0.20 * sortedAbs.length)];
    const fB = effB.filter(e => e >= top20).length / effB.length;
    const fN = effN.filter(e => e >= top20).length / effN.length;
    console.log(`  Fraction in top-20% strongest DMS effect: burden+ ${pct(effB.filter(e=>e>=top20).length,effB.length)}   non-burden ${pct(effN.filter(e=>e>=top20).length,effN.length)}`);
    console.log(`  Enrichment of strong-effect DMS at burden+ positions: ${(fB / (fN || 1e-9)).toFixed(2)}×`);
}

(async () => {
    const brca1 = 'abcdefghijklmnopqrstuvwxyz'.split('').map(L => `urn:mavedb:00000097-${L}-1`);
    console.log('=== S3 full deep-mutational-scanning validation ===');
    await runS3('P04637', 'TP53', ['urn:mavedb:00000068-b-1'], 'Giacomelli 2018 nutlin, whole protein');
    await runS3('P38398', 'BRCA1', brca1, 'Findlay 2018 SGE, all 13 exons');
})();
