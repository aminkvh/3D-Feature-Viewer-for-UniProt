/* ================================================================================
   Web Worker — ECFP/Morgan ligand similarity via RDKit (compiled to WebAssembly).
   Loaded from a web_accessible_resource by the content script. RDKit is initialised lazily on
   the first message (its WASM is ~7 MB) and then reused for the session. Computes the Morgan
   fingerprint (ECFP4: radius 2, 2048 bits) for each SMILES and the exact Tanimoto coefficient.
   Runs off the main thread so the WASM load / fingerprinting never freezes the UI.
   ================================================================================ */
let _rdkitReady = null;
function ensureRDKit() {
    if (_rdkitReady) return _rdkitReady;
    importScripts('lib/RDKit_minimal.js');
    _rdkitReady = self.initRDKitModule({
        locateFile: () => new URL('lib/RDKit_minimal.wasm', self.location.href).href,
    });
    return _rdkitReady;
}

const _fpCache = new Map(); // SMILES → Morgan fingerprint bit string ('0'/'1' chars) or null
function morganFp(RDKit, smiles) {
    if (_fpCache.has(smiles)) return _fpCache.get(smiles);
    let fp = null;
    try {
        const m = RDKit.get_mol(smiles || '');
        if (m) {
            if (m.is_valid()) fp = m.get_morgan_fp(JSON.stringify({ radius: 2, nBits: 2048 }));
            m.delete();
        }
    } catch (_) { /* invalid SMILES → no fingerprint */ }
    _fpCache.set(smiles, fp);
    return fp;
}

function tanimoto(a, b) {
    if (!a || !b) return 0;
    let inter = 0, union = 0;
    for (let k = 0; k < a.length; k++) {
        const x = a.charCodeAt(k) - 48, y = b.charCodeAt(k) - 48; // '0'→0, '1'→1
        inter += x & y; union += x | y;
    }
    return union ? inter / union : 0;
}

self.onmessage = async (e) => {
    const { id, focusSmiles, others } = e.data || {};
    try {
        const RDKit = await ensureRDKit();
        const fp0 = morganFp(RDKit, focusSmiles);
        const scores = (others || []).map(o => ({ ccd: o.ccd, score: fp0 ? tanimoto(fp0, morganFp(RDKit, o.smiles)) : 0 }));
        self.postMessage({ id, result: scores });
    } catch (err) {
        self.postMessage({ id, error: (err && err.message) ? err.message : String(err) });
    }
};
