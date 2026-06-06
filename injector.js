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

    // Append the button inline into the heading (after the heading's own text content).
    // The button uses vertical-align:middle and zero top/bottom margin so it sits on the
    // heading's existing line without adding height — the sections below don't shift.
    function placeButton(heading, btn) {
        btn.classList.add('ufv-3d-btn--inline');
        heading.appendChild(btn);
    }

    function tryInjectPTMButton() {
        const existing = document.getElementById('ufv-btn-ptm');
        if (existing?.isConnected) return;
        const h = findSectionHeading('ptm_processing');
        if (!h) return;
        placeButton(h, UFVModal.createButton('ufv-btn-ptm', 'View PTMs in 3D', () => UFVModal.open('ptm')));
    }

    function tryInjectVariantButton() {
        const existing = document.getElementById('ufv-btn-variant');
        if (existing?.isConnected) return;
        const h = findSectionHeading('disease_variants');
        if (!h) return;
        UFVState.state.scrapedDiseases = scrapeDiseaseHeadings(h);
        placeButton(h, UFVModal.createButton('ufv-btn-variant', 'View Variants in 3D', () => UFVModal.open('variant')));
    }

    function tryInjectFeaturesButton() {
        const existing = document.getElementById('ufv-btn-features');
        if (existing?.isConnected) return;
        const h = findSectionFeaturesHeading('function', 'function');
        if (!h) return;
        placeButton(h, UFVModal.createButton('ufv-btn-features', 'View Sites in 3D', () => UFVModal.open('sites')));
    }

    function tryInjectDomainsButton() {
        const existing = document.getElementById('ufv-btn-domains');
        if (existing?.isConnected) return;
        const h = findSectionFeaturesHeading('family_and_domains', 'family & domains')
               || findSectionFeaturesHeading('family_&_domains', 'family & domains');
        if (!h) return;
        placeButton(h, UFVModal.createButton('ufv-btn-domains', 'View Domains in 3D', () => UFVModal.open('domains')));
    }

    function injectAllEntryButtons() {
        tryInjectPTMButton();
        tryInjectVariantButton();
        tryInjectFeaturesButton();
        tryInjectDomainsButton();
    }

    function allPresent() {
        return ['ufv-btn-ptm','ufv-btn-variant','ufv-btn-features','ufv-btn-domains']
            .every(id => document.getElementById(id)?.isConnected);
    }

    // The "Features" sub-viewer inside a section (Function, Family & Domains, …) is where that
    // section's features are tabulated, so the button belongs next to that <h3>Features</h3> —
    // not the section's top-level <h2>.
    function findSectionFeaturesHeading(sectionId, sectionText) {
        const section = document.getElementById(sectionId)
            || [...(document.querySelector('main') || document.body).querySelectorAll('h2')]
                .find(h => h.textContent.trim().toLowerCase() === sectionText)?.closest('section');
        if (!section) return null;
        for (const h of section.querySelectorAll('h3')) {
            if (h.textContent.trim().toLowerCase() === 'features') return h;
        }
        return null;
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

    function injectEntryButtons() {
        // MutationObserver: only re-injects when a button is MISSING (React dropped it).
        // Short-circuits immediately when all buttons are present, so it never fires during
        // normal scroll / nav-link highlighting updates — that was what caused the sidebar lag.
        if (!runtime.observer) {
            let _t = null;
            runtime.observer = new MutationObserver(() => {
                if (allPresent()) return;
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
        // PTM nav fix: dispatch scroll when PTM section enters the viewport, nudging UniProt's
        // scroll-spy to re-evaluate the active nav link.
        schedulePTMNavFix();
        [0, 120, 300, 600, 1000, 1600].forEach(d => setTimeout(injectAllEntryButtons, d));
    }

    function schedulePTMNavFix() {
        if (runtime.ptmNavObserver) return;
        const tryFind = () => {
            if (runtime.ptmNavObserver) return;
            const el = document.getElementById('ptm_processing')
                || [...document.querySelectorAll('section[id]')].find(s => s.id.toLowerCase().includes('ptm'));
            if (!el) { setTimeout(tryFind, 1500); return; }
            runtime.ptmNavObserver = new IntersectionObserver(entries => {
                if (entries.some(e => e.isIntersecting)) {
                    window.dispatchEvent(new Event('scroll'));
                    setTimeout(() => window.dispatchEvent(new Event('scroll')), 80);
                }
            }, { threshold: 0.05 });
            runtime.ptmNavObserver.observe(el);
        };
        tryFind();
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
            // Cancel entry-page injection machinery.
            if (runtime.observer) { runtime.observer.disconnect(); runtime.observer = null; }
            if (runtime.entryPoll) { clearInterval(runtime.entryPoll); runtime.entryPoll = null; }
            if (runtime.ptmNavObserver) { runtime.ptmNavObserver.disconnect(); runtime.ptmNavObserver = null; }
            // Remove all extension buttons — no button on the variant-viewer page.
            ['ufv-btn-ptm', 'ufv-btn-variant', 'ufv-btn-features', 'ufv-btn-domains', 'ufv-btn-variant-page'].forEach(id => document.getElementById(id)?.remove());
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
