/* ============================================================================
   Validation of the two published algorithms (SI S1 and S3).

   S1 — Pathogenic-variant enrichment hotspots
     Ground truth: ClinVar "reviewed by expert panel" or "multiple submitters,
     no conflicts" (3-4 ★) pathogenic variants = high-confidence true positives.
     Negative proxy: positions with only "no assertion criteria" benign variants
     = low-confidence benign.  Tests whether spatial clustering selects for
     higher-evidence positions among those with any pathogenic annotation.
     Also measures if expert-panel positions are enriched at hotspot residues
     versus non-hotspot pathogenic positions (partially non-circular because
     the algorithm uses binary pathogenic/benign labels, not review stars).

   S3 — Mutation/phenotype burden
     Ground truth: DMS (deep mutational scanning) functional scores from MaveDB.
     Burden-positive positions should show lower (more deleterious) functional
     scores in independent experimental assays.
     Proteins: TP53 (Giacomelli 2018, nutlin selection, urn:mavedb:00000068-b-1)
               BRCA1 (SGE Exon 2, urn:mavedb:00000097-a-1)

   Run:  node validation/validate-s1-s3.js
   ============================================================================ */
const fs = require('fs'), path = require('path');

// Load DataProcessor (for extractVariants)
const dpCode = fs.readFileSync(path.join(__dirname, '..', 'data.js'), 'utf8');
let DataProcessor; eval(dpCode.replace('const DataProcessor', 'DataProcessor'));

// Load UFVAnalysis (for computeHotspots and computeResidueBurden)
const anCode = fs.readFileSync(path.join(__dirname, '..', 'analysis.js'), 'utf8');
let UFVAnalysis; eval(anCode.replace('const UFVAnalysis', 'UFVAnalysis'));

async function txt(u) { const r = await fetch(u); if (!r.ok) throw new Error(u + ' -> ' + r.status); return r.text(); }
async function jsn(u) { const r = await fetch(u, { headers: { Accept: 'application/json' } }); if (!r.ok) throw new Error(u + ' -> ' + r.status); return r.json(); }

// ---- PDB CA parser ----
function pdbCA(pdb) {
    const out = [], seen = new Set();
    for (const line of pdb.split('\n')) {
        if (!line.startsWith('ATOM') || line.slice(12, 16).trim() !== 'CA') continue;
        const resi = parseInt(line.slice(22, 26));
        if (seen.has(resi)) continue; seen.add(resi);
        out.push({ resi, chain: 'A', atom: 'CA', x: +line.slice(30, 38), y: +line.slice(38, 46), z: +line.slice(46, 54) });
    }
    return out;
}

// Mock viewer that UFVAnalysis.computeHotspots needs to call viewer.getModel().selectedAtoms()
function mockViewer(caAtoms) {
    return { getModel: () => ({ selectedAtoms: () => caAtoms }) };
}

// AlphaFold-like structure stub for computeHotspots (resi == uniPos for AlphaFold)
const AF_STRUCTURE = { source: 'AlphaFold', chainIds: null, chainId: null, mappedRanges: [] };

// ClinVar review-star rating
function reviewStars(statusStr) {
    if (!statusStr) return 0;
    const s = statusStr.toLowerCase();
    if (s.includes('reviewed by expert panel') || s.includes('practice guideline')) return 4;
    if (s.includes('multiple submitters') && !s.includes('conflict')) return 3;
    if (s.includes('criteria provided')) return 2;
    return 1;
}

// Maximum star rating across all ClinVar significance entries for a variant
function maxStars(variant) {
    return Math.max(0, ...(variant.clinVarReviewStatus || '').split(',').map(s => reviewStars(s.trim())));
}

// ================================================================
// S1 VALIDATION
// ================================================================
async function validateS1(id, label) {
    console.log(`\n--- S1: ${label} (${id}) ---`);
    let variation, pdb;
    try {
        [variation, pdb] = await Promise.all([
            jsn(`https://www.ebi.ac.uk/proteins/api/variation/${id}`),
            txt(`https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v6.pdb`)
                .catch(() => txt(`https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v4.pdb`)),
        ]);
    } catch (e) { console.log('  SKIPPED:', e.message); return; }

    const variants = DataProcessor.extractVariants(variation);
    const cas = pdbCA(pdb);
    const viewer = mockViewer(cas);
    const caByResi = new Map(cas.map(a => [a.resi, a]));

    // Run S1 hotspot algorithm
    const hotspots = UFVAnalysis.computeHotspots(viewer, variants, AF_STRUCTURE);
    const hotPos = hotspots.merged;  // Map<uniPos, tier>

    if (!hotPos || hotPos.size === 0) { console.log('  No hotspots detected (too few classified variants)'); return; }

    // Classify positions by star rating and hotspot status
    const pathPositions = new Map(); // uniPos -> {maxStars, inHotspot}
    variants.forEach(v => {
        if (!caByResi.has(v.position)) return;
        const stars = maxStars(v);
        const isPat = /pathogenic|deleterious/i.test(v.consequence || '');
        if (!isPat) return;
        const ex = pathPositions.get(v.position);
        if (!ex || stars > ex.stars) pathPositions.set(v.position, { stars, inHotspot: hotPos.has(v.position) });
    });

    const all = [...pathPositions.values()];
    const inHot = all.filter(p => p.inHotspot);
    const notHot = all.filter(p => !p.inHotspot);
    const pct = (n, d) => d ? (n / d * 100).toFixed(1) + '%' : 'n/a';
    const fracHighStar = arr => arr.filter(p => p.stars >= 3).length / (arr.length || 1);

    console.log(`  Pathogenic positions in model: ${all.length} | in hotspot: ${inHot.length} | not in hotspot: ${notHot.length}`);
    console.log(`  Hotspot tiers: ${[...new Set([...hotPos.values()])].join(', ')}`);
    console.log(`  3★+ pathogenic at HOTSPOT positions: ${pct(inHot.filter(p=>p.stars>=3).length, inHot.length)}`);
    console.log(`  3★+ pathogenic at NON-HOTSPOT positions: ${pct(notHot.filter(p=>p.stars>=3).length, notHot.length)}`);
    console.log(`  3★+ pathogenic across ALL positions: ${pct(all.filter(p=>p.stars>=3).length, all.length)}`);
    const enr = fracHighStar(inHot) / (fracHighStar(all) || 1e-9);
    console.log(`  Enrichment of 3★+ at hotspot vs overall: ${enr.toFixed(2)}×`);

    // Also check star-0 positions (no assertion criteria) — should be LESS represented in hotspots
    const star0Hot = inHot.filter(p => p.stars <= 1).length / (inHot.length || 1);
    const star0All = all.filter(p => p.stars <= 1).length / (all.length || 1);
    console.log(`  Low-evidence (≤1★) fraction in hotspots: ${(star0Hot*100).toFixed(1)}%  vs overall: ${(star0All*100).toFixed(1)}%  → hotspots ${star0Hot < star0All ? 'deplete' : 'do NOT deplete'} low-evidence variants`);
}

// ================================================================
// S3 VALIDATION
// ================================================================
function parseHgvsPro(hgvs) {
    // Parse "p.Arg273His" → position 273, or "p.Glu32Val" → 32
    if (!hgvs || hgvs === 'NA') return null;
    const m = hgvs.match(/p\.[A-Za-z*]{1,3}(\d+)/);
    return m ? parseInt(m[1]) : null;
}

async function parseMavedbScores(urn) {
    const r = await fetch(`https://api.mavedb.org/api/v1/score-sets/${urn}/scores`);
    if (!r.ok) return null;
    const text = await r.text();
    const lines = text.split(/\r?\n/);                          // MaveDB uses CRLF
    const header = lines[0].split(',').map(h => h.trim());
    const scoreIdx = header.indexOf('score');
    const hgvsIdx = header.indexOf('hgvs_pro');
    if (scoreIdx < 0 || hgvsIdx < 0) return null;
    const byPos = new Map(); // pos -> [scores]
    for (const line of lines.slice(1)) {
        if (!line.trim()) continue;
        const cols = line.split(',');
        const hgvs = (cols[hgvsIdx] || '').trim();
        const score = parseFloat((cols[scoreIdx] || '').trim());
        const pos = parseHgvsPro(hgvs);
        if (pos == null || !Number.isFinite(score)) continue;
        // skip synonymous (p.Xxx123= ) — not informative for missense effect
        if (/=\s*$/.test(hgvs)) continue;
        if (!byPos.has(pos)) byPos.set(pos, []);
        byPos.get(pos).push(score);
    }
    // Median score per position
    const medByPos = new Map();
    byPos.forEach((scores, pos) => {
        const s = scores.slice().sort((a, b) => a - b);
        medByPos.set(pos, s[s.length >> 1]);
    });
    return medByPos;
}

function median(arr) {
    if (!arr.length) return NaN;
    const s = arr.slice().sort((a, b) => a - b);
    return s[s.length >> 1];
}
function mean(arr) { return arr.length ? arr.reduce((a, b) => a + b, 0) / arr.length : NaN; }

async function validateS3(id, label, maveUrn, maveLabel, maveOffset = 0) {
    // maveOffset: add to MaveDB position to get UniProt position (for partial constructs)
    console.log(`\n--- S3: ${label} (${id}) | DMS: ${maveLabel} ---`);
    let variation, dmsScores;
    try {
        [variation, dmsScores] = await Promise.all([
            jsn(`https://www.ebi.ac.uk/proteins/api/variation/${id}`),
            parseMavedbScores(maveUrn),
        ]);
    } catch (e) { console.log('  SKIPPED:', e.message); return; }
    if (!dmsScores || dmsScores.size === 0) { console.log('  MaveDB scores unavailable'); return; }

    const variants = DataProcessor.extractVariants(variation);
    const burden = UFVAnalysis.computeResidueBurden(variants);

    if (burden.size === 0) { console.log('  No burden-positive residues detected'); return; }

    // Map MaveDB scores to UniProt positions
    const dmsUni = new Map();
    dmsScores.forEach((score, pos) => dmsUni.set(pos + maveOffset, score));

    // Split positions into burden+ and burden-
    const burdenPos = new Set(burden);
    const burdenPosScores = [], nonBurdenScores = [];
    dmsUni.forEach((score, uniPos) => {
        if (burdenPos.has(uniPos)) burdenPosScores.push(score);
        else nonBurdenScores.push(score);
    });

    const pct = (n, d) => d ? (n / d * 100).toFixed(1) + '%' : 'n/a';
    console.log(`  Burden-positive positions: ${burdenPos.size} | with DMS coverage: ${burdenPosScores.length} | non-burden w/ DMS: ${nonBurdenScores.length}`);
    if (burdenPosScores.length === 0) { console.log('  No overlap between burden positions and DMS coverage'); return; }

    // Direction-agnostic: a functionally important position has a LARGE-MAGNITUDE DMS effect
    // (mutations there strongly change function, either gain or loss). Measure |score − neutral|
    // where neutral = the median DMS score across all positions (the bulk-tolerant baseline).
    const allScores = [...dmsUni.values()];
    const neutral = median(allScores);
    const effB = burdenPosScores.map(s => Math.abs(s - neutral));
    const effN = nonBurdenScores.map(s => Math.abs(s - neutral));
    console.log(`  DMS |effect| from neutral (median): burden+ ${median(effB).toFixed(3)}   non-burden ${median(effN).toFixed(3)}  → burden+ effects are ${median(effB) > median(effN) ? 'LARGER (expected)' : 'not larger'}`);

    // Tail enrichment: fraction of positions among the top-20% strongest |effect| genome-wide.
    const sortedEff = allScores.map(s => Math.abs(s - neutral)).sort((a, b) => b - a);
    const top20 = sortedEff[Math.floor(0.20 * sortedEff.length)];
    const fracB = effB.filter(e => e >= top20).length / effB.length;
    const fracN = effN.filter(e => e >= top20).length / effN.length;
    console.log(`  Fraction in top-20% strongest DMS effect: burden+ ${pct(effB.filter(e=>e>=top20).length,effB.length)}   non-burden ${pct(effN.filter(e=>e>=top20).length,effN.length)}`);
    console.log(`  Enrichment of strong-effect DMS at burden+ positions: ${(fracB / (fracN || 1e-9)).toFixed(2)}×`);
}

(async () => {
    console.log('=== S1: Pathogenic-variant enrichment hotspot validation ===');
    await validateS1('P04637', 'TP53 (disease mutations, DNA-binding domain)');
    await validateS1('P01116', 'KRAS (cancer driver, G12 hotspot)');
    await validateS1('P38398', 'BRCA1 (hereditary breast cancer)');

    console.log('\n\n=== S3: Mutation/phenotype burden validation vs DMS ===');
    // TP53 Giacomelli 2018 nutlin selection (negative score = deleterious to p53 function)
    await validateS3('P04637', 'TP53', 'urn:mavedb:00000068-b-1', 'Giacomelli 2018 nutlin-3', 0);
    // BRCA1 SGE Exon 2 (covers residues in RING domain, offset ~2 for construct start)
    await validateS3('P38398', 'BRCA1', 'urn:mavedb:00000097-a-1', 'Findlay SGE Exon2', 0);
})();
