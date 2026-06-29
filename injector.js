/* global UFVState, UFVModal, StructureViewer */
const UFVInjector = (() => {
    'use strict';

    const runtime = window.__UFV_RUNTIME__ || {
        booted: false,
        observer: null,
        entryPoll: null,
        variantPoll: null,
        featurePoll: null,
        urlPoll: null,
        scrollBound: false,
        visibilityBound: false,
        historyPatched: false,
        lastUrl: '',
    };
    window.__UFV_RUNTIME__ = runtime;
    if (!runtime.injectedIds) runtime.injectedIds = new Set(); // ids of buttons we've actually placed

    function getUniProtId() {
        const m = window.location.pathname.match(/\/uniprotkb\/([A-Za-z0-9_-]+)/);
        return m ? m[1].toUpperCase() : null;
    }

    function isVariantViewerPage() {
        return window.location.pathname.includes('/variant-viewer');
    }

    function isFeatureViewerPage() {
        return window.location.pathname.includes('/feature-viewer');
    }

    function scrapeDiseaseHeadings(sectionHeading) {
        const diseases = [];
        let section = document.getElementById('disease_variants') || sectionHeading.closest('section, .card, [class*="card"]') || document.body;
        section.querySelectorAll('h4').forEach(h4 => {
            const link = h4.querySelector('a[href*="/diseases/"]');
            if (!link) return;
            const id = (link.getAttribute('href') || '').match(/\/diseases\/([^/?#]+)/)?.[1] || '';
            const label = h4.textContent.trim().replace(/^["']|["']$/g, '');
            diseases.push({ id, label });
        });
        return diseases;
    }

    function findSectionHeading(sectionId) {
        const section = document.getElementById(sectionId);
        if (section) {
            // Only treat a real heading element as the anchor.  Returning the <section>
            // itself (when its <h2> hasn't rendered yet) would append the button to the very
            // bottom of the section AND mark injection as "done" — stranding the button far
            // from the title until a manual refresh.  Returning null instead lets the poll /
            // MutationObserver retry until the heading actually exists.  This was the main
            // cause of "sometimes I have to refresh to see the button".
            const h = section.querySelector('h1, h2, h3');
            if (h) return h;
        }
        const textMap = {
            ptm_processing: ['PTM/Processing', 'PTM / Processing'],
            disease_variants: ['Disease & Variants', 'Disease and Variants'],
            function: ['Function'],
            structure: ['Structure'],
            subcellular_location: ['Subcellular location', 'Subcellular Location'],
        };
        const main = document.querySelector('main') || document.querySelector('[role="main"]') || document.body;
        for (const h of main.querySelectorAll('h1, h2, h3')) {
            const text = h.textContent.trim().toLowerCase();
            if ((textMap[sectionId] || []).some(term => text.includes(term.toLowerCase()))) return h;
        }
        return null;
    }

    // Append the button after the heading text, wrapped in a ZERO-SIZE inline anchor. The button is
    // absolutely positioned relative to that anchor, so it adds NO height to the heading's line box and
    // never shifts the sections below it. That keeps UniProt's left-nav scroll-spy section offsets intact —
    // the real cause of the active-section indicator landing on the wrong item (no scroll/resize nudges
    // needed once the layout is genuinely untouched).
    function placeButton(heading, btn) {
        btn.classList.add('ufv-3d-btn--inline');
        const anchor = document.createElement('span');
        anchor.className = 'ufv-3d-btn-anchor';
        anchor.appendChild(btn);
        heading.appendChild(anchor);
        if (btn.id) runtime.injectedIds.add(btn.id); // remember it so the observer re-injects only if dropped
    }

    function tryInjectPTMButton() {
        const existing = document.getElementById('ufv-btn-ptm');
        if (existing?.isConnected) return;
        // Anchor to the <h3>Features</h3> inside PTM/Processing; fall back to the <h2> if absent.
        const h = findSectionSubHeading('ptm_processing', ['ptm/processing', 'ptm / processing'], ['features'])
               || findSectionHeading('ptm_processing');
        if (!h) return;
        placeButton(h, UFVModal.createButton('ufv-btn-ptm', 'View PTMs in 3D', () => UFVModal.open('ptm')));
    }

    function tryInjectVariantButton() {
        const existing = document.getElementById('ufv-btn-variant');
        if (existing?.isConnected) return;
        // Anchor to the <h3>Involvement in disease</h3>; fall back to the <h2> if absent.
        const h = findSectionSubHeading('disease_variants', ['disease & variants', 'disease and variants'], ['involvement in disease'])
               || findSectionHeading('disease_variants');
        if (!h) return;
        UFVState.state.scrapedDiseases = scrapeDiseaseHeadings(h);
        placeButton(h, UFVModal.createButton('ufv-btn-variant', 'View Variants in 3D', () => UFVModal.open('variant')));
    }

    function tryInjectFeaturesButton() {
        const existing = document.getElementById('ufv-btn-features');
        if (existing?.isConnected) return;
        const h = findSectionSubHeading('function', ['function'], ['features']);
        if (!h) return;
        placeButton(h, UFVModal.createButton('ufv-btn-features', 'View Sites in 3D', () => UFVModal.open('sites')));
    }

    function tryInjectDomainsButton() {
        const existing = document.getElementById('ufv-btn-domains');
        if (existing?.isConnected) return;
        const h = findSectionSubHeading('family_and_domains', ['family & domains'], ['features'])
               || findSectionSubHeading('family_&_domains', ['family & domains'], ['features']);
        if (!h) return;
        placeButton(h, UFVModal.createButton('ufv-btn-domains', 'View Domains in 3D', () => UFVModal.open('domains')));
    }

    function tryInjectStructureButton() {
        const existing = document.getElementById('ufv-btn-structure');
        if (existing?.isConnected) return;
        const h = findSectionHeading('structure');
        if (!h) return;
        placeButton(h, UFVModal.createButton('ufv-btn-structure', 'View in 3D', () => UFVModal.open('structure')));
    }

    function tryInjectSubcellularButton() {
        const existing = document.getElementById('ufv-btn-subcellular');
        if (existing?.isConnected) return;
        // Anchor to the "Features" sub-heading inside the Subcellular location section. If the entry has
        // no subcellular Features sub-heading, the button is intentionally NOT shown.
        const h = findSectionSubHeading('subcellular_location', ['subcellular location'], ['features']);
        if (!h) return;
        placeButton(h, UFVModal.createButton('ufv-btn-subcellular', 'View in 3D', () => UFVModal.open('subcellular')));
    }

    function injectAllEntryButtons() {
        tryInjectPTMButton();
        tryInjectVariantButton();
        tryInjectFeaturesButton();
        tryInjectDomainsButton();
        tryInjectStructureButton();
        tryInjectSubcellularButton();
    }

    // Find an <h3> sub-heading (matching one of h3Texts) inside a section. We anchor buttons to the
    // <h3> rather than the section's <h2> on purpose: UniProt's left-nav scroll-spy tracks the
    // <h2>, so injecting a child into the <h2> disturbs the active-section indicator. The <h3> is
    // safe. Returns null when the sub-heading isn't present (caller can fall back to the <h2>).
    function findSectionSubHeading(sectionId, sectionTexts, h3Texts) {
        let section = document.getElementById(sectionId);
        if (!section) {
            const main = document.querySelector('main') || document.body;
            section = [...main.querySelectorAll('h2')]
                .find(h => sectionTexts.includes(h.textContent.trim().toLowerCase()))?.closest('section') || null;
        }
        if (!section) return null;
        for (const h of section.querySelectorAll('h3')) {
            if (h3Texts.includes(h.textContent.trim().toLowerCase())) return h;
        }
        return null;
    }

    function tryInjectVariantViewerButton() {
        const existing = document.getElementById('ufv-btn-variant-page');
        if (existing && existing.isConnected) return;
        // Anchor inline next to the page's "Variants" heading (the poll re-injects if React drops it).
        const main = document.querySelector('main') || document.body;
        const h = [...main.querySelectorAll('h1, h2')].find(x => /^variants?\b/i.test((x.textContent || '').trim()));
        if (!h) return;
        placeButton(h, UFVModal.createButton('ufv-btn-variant-page', 'View Variants in 3D', () => UFVModal.open('variant')));
    }

    function injectEntryButtons() {
        // MutationObserver: only re-injects when a button is MISSING (React dropped it).
        // Short-circuits immediately when all buttons are present, so it never fires during
        // normal scroll / nav-link highlighting updates — that was what caused the sidebar lag.
        if (!runtime.observer) {
            let _t = null;
            runtime.observer = new MutationObserver(() => {
                // Re-inject ONLY when a button we previously PLACED was dropped by a React re-render.
                // Buttons whose section doesn't exist on this entry (e.g. PTM/Processing absent on some
                // entries) never enter injectedIds, so a missing section can no longer keep the observer
                // firing on every mutation — which churned the DOM during scroll/lazy-load and disturbed
                // UniProt's nav. (allPresent() never became true on those entries — the old bug.)
                let dropped = false;
                runtime.injectedIds.forEach(id => { if (!document.getElementById(id)?.isConnected) dropped = true; });
                if (!dropped) return;
                if (_t) return;
                _t = setTimeout(() => { _t = null; injectAllEntryButtons(); }, 200);
            });
            runtime.observer.observe(document.body, { childList: true, subtree: true });
        }
        if (!runtime.visibilityBound) {
            runtime.visibilityBound = true;
            document.addEventListener('visibilitychange', () => { if (!document.hidden) injectAllEntryButtons(); });
        }
        if (!runtime.entryPoll) {
            runtime.entryPoll = setInterval(injectAllEntryButtons, 2000);
        }
        // Buttons are layout-neutral (zero-size anchor + absolute button) and the observer no longer churns
        // (see injectedIds), so the extension does NOT perturb UniProt's left-nav scroll-spy. We deliberately
        // do NOT touch the page's scroll/indicator ourselves — earlier scroll/resize/scrollIntoView "revive"
        // hooks fought UniProt's own settle mechanism and caused random scroll jumps. Let UniProt own it.
        [0, 120, 300, 600, 1000, 1600, 2500, 4000].forEach(d => setTimeout(injectAllEntryButtons, d));
    }

    function injectVariantViewerButton() {
        // Persistent poll (cheap — returns early if present): the button is now anchored to the page
        // heading, which React can re-render away, so keep re-injecting. Cleared on page navigation.
        if (!runtime.variantPoll) runtime.variantPoll = setInterval(tryInjectVariantViewerButton, 800);
        [0, 120, 300, 600, 1000].forEach(d => setTimeout(tryInjectVariantViewerButton, d));
    }

    // The expanded feature-viewer page (/uniprotkb/<id>/feature-viewer) lists every feature. The button
    // opens the unified Structure window (all layers available, toggle on what you want), anchored inline
    // next to the "Feature viewer" heading.
    function tryInjectFeatureViewerButton() {
        const existing = document.getElementById('ufv-btn-feature-page');
        if (existing && existing.isConnected) return;
        const main = document.querySelector('main') || document.body;
        const heads = [...main.querySelectorAll('h1, h2, h3')].filter(x => x.offsetParent !== null); // visible headings
        // Prefer a "Feature viewer" / "Features" heading; otherwise fall back to the page's first visible
        // heading (the feature-viewer page sometimes titles itself by the protein name, not "Feature viewer").
        const h = heads.find(x => /feature\s*viewer|^features?\b/i.test((x.textContent || '').trim())) || heads[0];
        if (!h) return;
        placeButton(h, UFVModal.createButton('ufv-btn-feature-page', 'View Features in 3D', () => UFVModal.open('structure')));
    }

    function injectFeatureViewerButton() {
        // Persistent poll (anchored to the heading; cleared on page navigation).
        if (!runtime.featurePoll) runtime.featurePoll = setInterval(tryInjectFeatureViewerButton, 800);
        [0, 120, 300, 600, 1000].forEach(d => setTimeout(tryInjectFeatureViewerButton, d));
    }

    function handlePageType() {
        const id = getUniProtId();
        if (!id) return;
        if (UFVState.state.uniprotId !== id) {
            UFVState.resetForProtein(id);
            UFVModal.close();
            // Drop the previous protein's loaded model + caches so a stale structure can't be
            // mistaken for an already-loaded one (and so the next open starts clean).
            try { StructureViewer?.clearModel(); } catch (_) {}
            // Explicitly remove previous protein's injected buttons so tryInject* doesn't treat
            // them as already-placed before React unmounts the old page DOM.
            const BTNS = ['ufv-btn-ptm', 'ufv-btn-variant', 'ufv-btn-features', 'ufv-btn-domains', 'ufv-btn-structure', 'ufv-btn-subcellular'];
            BTNS.forEach(btnId => document.getElementById(btnId)?.remove());
            runtime.injectedIds.clear();
        }
        UFVState.state.pageContext = isVariantViewerPage() ? 'variant-viewer' : 'entry';
        // Kick off background annotation fetch so the modal opens instantly
        UFVModal.prefetchData();
        const ENTRY_BTNS = ['ufv-btn-ptm', 'ufv-btn-variant', 'ufv-btn-features', 'ufv-btn-domains', 'ufv-btn-structure', 'ufv-btn-subcellular'];
        // Stop the entry-page injection machinery — shared by the variant-viewer and feature-viewer sub-pages.
        const stopEntryMachinery = () => {
            if (runtime.observer) { runtime.observer.disconnect(); runtime.observer = null; }
            if (runtime.entryPoll) { clearInterval(runtime.entryPoll); runtime.entryPoll = null; }
            runtime.injectedIds.clear();
        };
        if (isVariantViewerPage()) {
            stopEntryMachinery();
            if (runtime.featurePoll) { clearInterval(runtime.featurePoll); runtime.featurePoll = null; }
            [...ENTRY_BTNS, 'ufv-btn-variant-page', 'ufv-btn-feature-page'].forEach(id => document.getElementById(id)?.remove());
            injectVariantViewerButton();
        } else if (isFeatureViewerPage()) {
            stopEntryMachinery();
            if (runtime.variantPoll) { clearInterval(runtime.variantPoll); runtime.variantPoll = null; }
            [...ENTRY_BTNS, 'ufv-btn-variant-page', 'ufv-btn-feature-page'].forEach(id => document.getElementById(id)?.remove());
            injectFeatureViewerButton();
        } else {
            // Entry page: cancel the sub-page polls so they can't re-inject their fixed buttons.
            if (runtime.variantPoll) { clearInterval(runtime.variantPoll); runtime.variantPoll = null; }
            if (runtime.featurePoll) { clearInterval(runtime.featurePoll); runtime.featurePoll = null; }
            document.getElementById('ufv-btn-variant-page')?.remove();
            document.getElementById('ufv-btn-feature-page')?.remove();
            injectEntryButtons();
        }
    }

    function patchHistory() {
        if (runtime.historyPatched) return;
        runtime.historyPatched = true;
        const origPush = history.pushState;
        const origReplace = history.replaceState;
        history.pushState = function() {
            origPush.apply(this, arguments);
            onUrlChange();
        };
        history.replaceState = function() {
            origReplace.apply(this, arguments);
            onUrlChange();
        };
        window.addEventListener('popstate', onUrlChange);
        runtime.urlPoll = runtime.urlPoll || setInterval(onUrlChange, 1000);
    }

    function onUrlChange() {
        if (window.location.href === runtime.lastUrl) return;
        const prevPath = runtime.lastUrl.split('#')[0];
        const newPath = window.location.href.split('#')[0];
        runtime.lastUrl = window.location.href;
        // Ignore hash-only changes (sidebar anchor links) — they don't switch pages
        // and firing handlePageType() would race with UniProt's sidebar indicator updates.
        if (prevPath === newPath) return;
        runtime.injectedIds.clear(); // new page/protein → forget which buttons we'd placed
        handlePageType();
    }

    function boot() {
        runtime.lastUrl = window.location.href;
        // Patch history FIRST and unconditionally. The content script now loads across /uniprotkb* (incl.
        // the search/listing pages and the bare /uniprotkb/<id> URL). The user often lands on one of those
        // and then client-side-navigates into an entry; without an early history patch that SPA navigation
        // is never observed and the button only appears after a manual page refresh.
        patchHistory();
        const id = getUniProtId();
        if (!id) return;
        if (!UFVState.state.uniprotId) UFVState.resetForProtein(id);
        handlePageType();
        // Restored from the back/forward (bfcache) cache: the content script doesn't re-run
        // and the injection timers may have been frozen, so re-arm injection on pageshow.
        if (!runtime.pageshowBound) {
            runtime.pageshowBound = true;
            window.addEventListener('pageshow', e => { if (e.persisted) handlePageType(); });
        }
    }

    return { boot, getUniProtId };
})();
