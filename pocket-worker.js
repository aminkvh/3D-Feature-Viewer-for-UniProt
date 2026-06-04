/* ================================================================================
   Web Worker — off-main-thread compute for the constraint-pocket coloring mode.
   Keeps the permutation test / burial / spatial weighting from freezing the UI on large
   proteins.  Loaded from a web_accessible_resource by the content script; pulls in the shared
   algorithm code via importScripts (which sets self.UFVPocket).  The content script
   transparently falls back to a synchronous main-thread compute if this worker can't be created
   or errors out (e.g. a restrictive page CSP), so the feature never hard-fails.
   ================================================================================ */
/* global UFVPocket */
try {
    importScripts('algorithms.js');
} catch (_) {
    // Surfaced per-message below if UFVPocket ends up unavailable.
}

self.onmessage = (e) => {
    const { geometry, amEntries, sequence, pae } = e.data || {};
    try {
        const mod = self.UFVPocket || (typeof UFVPocket !== 'undefined' ? UFVPocket : null);
        if (!mod) throw new Error('algorithms module failed to load in worker');
        const amMap = new Map(amEntries);
        const res = mod.computePockets(geometry, amMap, sequence, pae, { maxN: 3000 });
        // Map isn't cloneable in a friendly way — send entries and rebuild on the main side.
        self.postMessage({ result: { byPos: [...res.byPos], reason: res.reason, hasPae: res.hasPae, n: res.n } });
    } catch (err) {
        self.postMessage({ error: (err && err.message) ? err.message : String(err) });
    }
};
