# Milestone 0 — Mol* spike (go / no-go)

Goal: prove Mol\* renders a structure as an **extension page** before we port the real viewer.
Because the page runs under the *extension's* origin and CSP (not UniProt's), hosting Mol\* in an
iframe sidesteps UniProt's page CSP entirely — this is why we chose the iframe approach.

## What was added
- `lib/molstar.js` + `lib/molstar.css` — Mol\* viewer bundle, pinned to **5.10.1** (~4.95 MB).
- `viewer-frame.html` — the extension page that hosts Mol\* (loads the bundle, full-window canvas).
- `viewer-frame.js` — bootstrap: creates the viewer, self-tests via a query param, and seeds the
  parent⇄frame `postMessage` protocol M1 will build on.
- `manifest.json` — `viewer-frame.html` is now a **sandboxed page** (`sandbox.pages`) with a
  `content_security_policy.sandbox` that permits `'unsafe-eval'`, plus a `web_accessible_resources`
  entry so it can be embedded as an iframe on `uniprot.org`.

## Why a sandboxed iframe (the key constraint)
Two CSP failures, in order:
1. `window.molstar missing` under the default CSP — Mol\* **instantiates WASM at load**.
2. After adding `'wasm-unsafe-eval'`: `EvalError … 'unsafe-eval' is not allowed`, thrown from
   `lib/molstar.js` Emscripten glue (`dynCall` builds call trampolines with `new Function`).

So Mol\*'s WASM needs **both** `wasm-unsafe-eval` *and* full `unsafe-eval`. MV3 **forbids
`unsafe-eval` on normal extension pages** but **allows it in a page's `sandbox` CSP**. Hence the
viewer must live in a **sandboxed iframe**. This is also why the whole architecture is iframe +
`postMessage`: the sandboxed frame can't call `chrome.*`, so the parent drives it by messages.

Notes:
- WASM is **inlined as base64** in the 4.95 MB bundle (no external `.wasm`), so enabling it adds **no
  extra download**. It can't be dropped without a custom Mol\* build (`DYNAMIC_EXECUTION=0`).
- Ligand-similarity/Tanimoto is unrelated — plain JS over PubChem fingerprints, no WASM.
- **Firefox caveat:** Firefox ignores the `sandbox` manifest key, so this exact path is Chrome-first.
  Firefox handling (separate build/manifest, or a Mol\* eval-free build) is deferred to M1.

## How to verify (you run this — I can't drive a browser)
1. Load the extension unpacked:
   - **Chrome:** `chrome://extensions` → enable *Developer mode* → *Load unpacked* → select the repo
     folder. Copy the extension **ID**.
   - **Firefox:** `about:debugging#/runtime/this-firefox` → *Load Temporary Add-on* → pick
     `manifest.json`. Note the internal UUID (Inspect → URL).
2. Open the **host** page (which embeds the sandboxed frame as an iframe — the only way Chrome
   applies the `sandbox` CSP):
   ```
   chrome-extension://<ID>/spike-host.html
   ```
   - The host embeds `viewer-frame.html?af=P35498`. Edit `spike-host.html` to try another accession
     or `?url=<cif-or-pdb>&format=mmcif`.
   - **Do NOT open `viewer-frame.html` directly** as a top-level tab — Chrome applies the *sandbox*
     CSP only to subframes; top-level it falls back to `extension_pages` and Mol*'s WASM/`new Function`
     glue is blocked (the `'wasm-unsafe-eval'` CSP error). It must be loaded inside an iframe.

After any manifest change, reload the extension. **For `sandbox` key changes the ↻ button is not
enough** — Chrome registers sandboxed pages only at load time, so **Remove the extension and
*Load unpacked* again**. Confirm the sandbox is active: in the frame's DevTools console,
`window.origin` should print `"null"` (opaque origin = sandboxed). If it prints `chrome-extension://…`,
the sandbox isn't applied and Mol*'s eval glue will be blocked.

### Success
The P35498 AlphaFold model renders as a cartoon you can rotate/zoom; the "Loading…" overlay
disappears and the console is free of CSP errors. → **GO.** The sandboxed extension page renders
Mol\*, so embedding it as an iframe in the modal (M1) will work.

### Failure to watch for (open the page's DevTools console)
- **Still an `unsafe-eval` CSP error** → the sandbox CSP didn't take; confirm the extension was
  reloaded and `manifest.json` has both `sandbox.pages` and `content_security_policy.sandbox`.
- **Blank canvas / WebGL error** → GPU/WebGL2 issue in your browser, unrelated to our wiring.
- **CORS / network error fetching the CIF** → the sandboxed frame fetches with a null/extension
  origin; AlphaFold/EBI send `Access-Control-Allow-Origin: *` so this should be fine. In M1 the
  parent will fetch (using host permissions) and post the text to the frame, removing all frame CORS.

## Recommendation
**Sandboxed iframe.** Forced by Mol\*'s Emscripten WASM glue (`new Function`), which needs
`'unsafe-eval'` — only allowed in a `sandbox` CSP under MV3. This confirms the architecture: a
sandboxed `viewer-frame.html` embedded in the modal, driven entirely over `postMessage` (load
structure, colour, focus, …). Firefox parity is a separate M1 task (it ignores `sandbox.pages`).

## Nothing else changed
The existing 3Dmol-based extension is untouched and still the active viewer. This spike only adds new
files + a `web_accessible_resources` entry + the sandbox declaration.
