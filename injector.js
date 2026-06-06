/* global UFVState, UFVModal, StructureViewer */
const UFVInjector = (() => {
    'use strict';

    const runtime = window.__UFV_RUNTIME__ || {
        booted: false,
        observer: null,
        entryPoll: null,
        variantPoll: null,
        urlPoll: null,
        scrollBound: false,
        visibilityBound: false,
        historyPatched: false,
        lastUrl: '',
    };
    window.__UFV_RUNTIME__ = runtime;

    function getUniProtId() {
        const m = window.location.pathname.match(/\/uniprotkb\/([A-Za-z0-9_-]+)/);
        return m ? m[1].toUpperCase() : null;
    }

    function isVariantViewerPage() {
        return window.location.pathname.includes('/variant-viewer');
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
        };
        const main = document.querySelector('main') || document.querySelector('[role="main"]') || document.body;
        for (const h of main.querySelectorAll('h1, h2, h3')) {
            const text = h.textContent.trim().toLowerCase();
            if ((textMap[sectionId] || []).some(term => text.includes(term.toLowerCase()))) return h;
        }
        return null;
    }

    // Place the button inline, right after the section title text (inside the heading element),
    // so it sits "by the title" yet stays on the heading's existing line box — it does not add
    // height or reflow the sections below it, which would otherwise disturb UniProt's
    // scroll/IntersectionObserver-driven left-nav active slider.
    function placeAnchoredButton(heading, btn) {
        btn.classList.add('ufv-3d-btn--inline');
        heading.appendChild(btn);
    }

    function tryInjectPTMButton() {
        const existing = document.getElementById('ufv-btn-ptm');
        if (existing && document.body.contains(existing)) return;
        const heading = findSectionHeading('ptm_processing');
        if (!heading) return;
        const btn = UFVModal.createButton('ufv-btn-ptm', 'View PTMs in 3D', () => UFVModal.open('ptm'));
        placeAnchoredButton(heading, btn);
    }

    function tryInjectVariantButton() {
        const existing = document.getElementById('ufv-btn-variant');
        if (existing && document.body.contains(existing)) return;
        const heading = findSectionHeading('disease_variants');
        if (!heading) return;
        UFVState.state.scrapedDiseases = scrapeDiseaseHeadings(heading);
        const btn = UFVModal.createButton('ufv-btn-variant', 'View Variants in 3D', () => UFVModal.open('variant'));
        placeAnchoredButton(heading, btn);
    }

    // The "Features" sub-viewer inside the Function section is where active/binding/metal sites
    // are tabulated, so the button belongs next to that <h3>Features</h3> — not the section's
    // top-level <h2>Function</h2>.
    function findFunctionFeaturesHeading() {
        const section = document.getElementById('function')
            || [...(document.querySelector('main') || document.body).querySelectorAll('h2')]
                .find(h => h.textContent.trim().toLowerCase() === 'function')?.closest('section');
        if (!section) return null;
        for (const h of section.querySelectorAll('h3')) {
            if (h.textContent.trim().toLowerCase() === 'features') return h;
        }
        return null;
    }

    // Features sub-section under Function (holds active/binding/metal sites). Opens with sites shown.
    function tryInjectFeaturesButton() {
        const existing = document.getElementById('ufv-btn-features');
        if (existing && document.body.contains(existing)) return;
        const heading = findFunctionFeaturesHeading();
        if (!heading) return;
        const btn = UFVModal.createButton('ufv-btn-features', 'View Sites in 3D', () => UFVModal.open('ptm', { sitesOn: true }));
        placeAnchoredButton(heading, btn);
    }

    function tryInjectVariantViewerButton() {
        const existing = document.getElementById('ufv-btn-variant-page');
        if (existing && document.body.contains(existing)) return;
        const btn = UFVModal.createButton('ufv-btn-variant-page', 'View Variants in 3D', () => UFVModal.open('variant'));
        // Always attach to body with fixed positioning — anchoring to a DOM element
        // inside the variant-viewer causes the button to jump when React re-renders
        // on tab switches (each render may produce a different h2 target).
        btn.style.cssText = 'position:fixed;top:80px;right:20px;z-index:9999;box-shadow:0 4px 16px rgba(0,0,0,0.3);';        document.body.appendChild(btn);
    }

    function injectAllEntryButtons() {
        tryInjectPTMButton();
        tryInjectVariantButton();
        tryInjectFeaturesButton();
    }
    function allEntryButtonsPresent() {
        return document.getElementById('ufv-btn-ptm')?.isConnected &&
            document.getElementById('ufv-btn-variant')?.isConnected &&
            document.getElementById('ufv-btn-features')?.isConnected;
    }

    function injectEntryButtons() {
        if (!runtime.observer) {
            // Debounce: React batch-renders cause many rapid mutations; collapse them
            // into a single inject attempt so we never mutate the DOM mid-scroll.
            let _injectTimer = null;
            runtime.observer = new MutationObserver(() => {
                // Skip when all buttons are already present — prevents timer churn during
                // scroll-triggered React re-renders.
                if (allEntryButtonsPresent()) return;
                if (_injectTimer) return;
                _injectTimer = setTimeout(() => { _injectTimer = null; injectAllEntryButtons(); }, 150);
            });
            runtime.observer.observe(document.body, { childList: true, subtree: true });
        }
        // Re-inject when switching back to this tab (UniProt SPA re-renders on visibility)
        if (!runtime.visibilityBound) {
            runtime.visibilityBound = true;
            document.addEventListener('visibilitychange', () => { if (!document.hidden) injectAllEntryButtons(); });
        }
        if (!runtime.entryPoll) {
            // Permanent low-frequency keepalive. UniProt's React app re-renders sections on
            // scroll / tab changes and silently drops our buttons, and the Function "Features"
            // sub-view renders lazily only when scrolled into view — so a one-shot poll that
            // stops once buttons exist left them gone until a manual refresh (the "inconsistent
            // injection" bug). Re-checking every second is cheap (each tryInject* short-circuits
            // on getElementById when present) and guarantees the buttons always come back.
            // Cleared only when navigating to a non-entry page (handlePageType).
            runtime.entryPoll = setInterval(injectAllEntryButtons, 1000);
        }
        // NOTE: scroll listener removed — MutationObserver already covers lazy-loaded
        // section re-renders, and the scroll listener was racing with UniProt's
        // IntersectionObserver that drives the sidebar active-indicator.
        // Immediate burst of attempts for a snappy first paint (covers the common case where
        // the sections are already in the DOM at document_idle).
        [0, 120, 300, 600, 1000, 1600].forEach(d => setTimeout(injectAllEntryButtons, d));
    }

    function injectVariantViewerButton() {
        if (!runtime.variantPoll) {
            let attempts = 0;
            runtime.variantPoll = setInterval(() => {
                tryInjectVariantViewerButton();
                attempts++;
                if (document.getElementById('ufv-btn-variant-page') || attempts > 40) {
                    clearInterval(runtime.variantPoll);
                    runtime.variantPoll = null;
                }
            }, 500);
        }
        [0, 120, 300, 600, 1000].forEach(d => setTimeout(tryInjectVariantViewerButton, d));
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
        }
        UFVState.state.pageContext = isVariantViewerPage() ? 'variant-viewer' : 'entry';
        // Kick off background annotation fetch so the modal opens instantly
        UFVModal.prefetchData();
        if (isVariantViewerPage()) {
            // Cancel any entry-page injection machinery so it can't re-inject
            // entry buttons (or re-trigger the MutationObserver) on this page.
            if (runtime.observer) { runtime.observer.disconnect(); runtime.observer = null; }
            if (runtime.entryPoll) { clearInterval(runtime.entryPoll); runtime.entryPoll = null; }
            // Remove all extension buttons — no button on the variant-viewer page.
            ['ufv-btn-ptm', 'ufv-btn-variant', 'ufv-btn-features', 'ufv-btn-variant-page'].forEach(id => document.getElementById(id)?.remove());
            injectVariantViewerButton();
        } else {
            // Cancel the variant-viewer poll so it can't re-inject ufv-btn-variant-page
            // into whatever main h2 happens to be present on the new page.
            if (runtime.variantPoll) { clearInterval(runtime.variantPoll); runtime.variantPoll = null; }
            // Remove variant-viewer button when back on the entry page
            document.getElementById('ufv-btn-variant-page')?.remove();
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
        handlePageType();
    }

    function boot() {
        runtime.lastUrl = window.location.href;
        const id = getUniProtId();
        if (!id) return;
        if (!UFVState.state.uniprotId) UFVState.resetForProtein(id);
        handlePageType();
        patchHistory();
        // Restored from the back/forward (bfcache) cache: the content script doesn't re-run
        // and the injection timers may have been frozen, so re-arm injection on pageshow.
        if (!runtime.pageshowBound) {
            runtime.pageshowBound = true;
            window.addEventListener('pageshow', e => { if (e.persisted) handlePageType(); });
        }
    }

    return { boot, getUniProtId };
})();
