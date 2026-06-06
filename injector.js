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

    // Entry-page buttons FLOAT on document.body, positioned over each section title — they are NOT
    // injected into UniProt's section headings. Putting a child into a heading made React re-render
    // that section, which both dropped our button (the "disappearing buttons" bug) and made the
    // section lose its left-nav scroll-spy registration (the PTM active-indicator bug). Keeping the
    // buttons out of the React-managed DOM avoids both: body children survive re-renders, and the
    // headings are never mutated. position:absolute in document coordinates means they scroll with
    // the page naturally; we only reposition when layout/structure changes.
    const FLOAT_BTN_SPECS = [
        { id: 'ufv-btn-ptm', label: 'View PTMs in 3D', open: () => UFVModal.open('ptm'),
          find: () => findSectionHeading('ptm_processing') },
        { id: 'ufv-btn-variant', label: 'View Variants in 3D', open: () => UFVModal.open('variant'),
          find: () => findSectionHeading('disease_variants'),
          onPosition: h => { try { UFVState.state.scrapedDiseases = scrapeDiseaseHeadings(h); } catch (_) {} } },
        { id: 'ufv-btn-features', label: 'View Sites in 3D', open: () => UFVModal.open('sites'),
          find: () => findSectionFeaturesHeading('function', 'function') },
        { id: 'ufv-btn-domains', label: 'View Domains in 3D', open: () => UFVModal.open('domains'),
          find: () => findSectionFeaturesHeading('family_and_domains', 'family & domains')
                     || findSectionFeaturesHeading('family_&_domains', 'family & domains') },
    ];

    // (Re)create each floating button and place it just after its section's title text. Cheap — a
    // getBoundingClientRect per button — so it's safe to call on mutations, resize, and a keepalive.
    function positionEntryButtons() {
        if (isVariantViewerPage()) return; // no entry buttons on the variant-viewer page
        FLOAT_BTN_SPECS.forEach(spec => {
            const heading = spec.find();
            let btn = document.getElementById(spec.id);
            if (!heading) { if (btn) btn.style.display = 'none'; return; }
            if (!btn) {
                btn = UFVModal.createButton(spec.id, spec.label, spec.open);
                btn.classList.add('ufv-3d-btn--float');
                document.body.appendChild(btn);
            }
            if (spec.onPosition) spec.onPosition(heading);
            // Measure where the title text ends so the button sits right after it (not at the far
            // right of a full-width heading). A Range over the heading's contents gives that.
            let tr;
            try { const range = document.createRange(); range.selectNodeContents(heading); tr = range.getBoundingClientRect(); }
            catch (_) { tr = heading.getBoundingClientRect(); }
            if (tr.width === 0 && tr.height === 0) { btn.style.display = 'none'; return; }
            btn.style.display = '';
            const top = window.scrollY + tr.top + tr.height / 2 - btn.offsetHeight / 2;
            let left = window.scrollX + tr.right + 12;
            const maxLeft = window.scrollX + document.documentElement.clientWidth - btn.offsetWidth - 10;
            if (left > maxLeft) left = Math.max(window.scrollX, maxLeft);
            btn.style.top = Math.max(0, top) + 'px';
            btn.style.left = Math.max(0, left) + 'px';
        });
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
        // Reposition (and create-if-missing) on any DOM change — covers sections rendering late and
        // React re-renders that move a heading. The buttons live on body, so this never removes
        // them; it only updates their position, and it never touches the headings.
        if (!runtime.observer) {
            let _t = null;
            runtime.observer = new MutationObserver(() => {
                if (_t) return;
                _t = setTimeout(() => { _t = null; positionEntryButtons(); }, 150);
            });
            runtime.observer.observe(document.body, { childList: true, subtree: true });
        }
        if (!runtime.visibilityBound) {
            runtime.visibilityBound = true;
            document.addEventListener('visibilitychange', () => { if (!document.hidden) positionEntryButtons(); });
        }
        if (!runtime.scrollBound) {
            runtime.scrollBound = true;
            // Reposition on resize (layout width changes move the title-text end). Scrolling needs
            // no handler: the buttons are absolute in document coordinates, so they scroll with the
            // page. rAF-coalesced to avoid redundant work.
            let raf = 0;
            const reflow = () => { if (raf) return; raf = requestAnimationFrame(() => { raf = 0; positionEntryButtons(); }); };
            window.addEventListener('resize', reflow, { passive: true });
        }
        if (!runtime.entryPoll) {
            // Keepalive: catches headings that render later and layout shifts (lazy images/fonts)
            // that move a heading without firing a childList mutation. Cheap and idempotent.
            runtime.entryPoll = setInterval(positionEntryButtons, 1000);
        }
        // Immediate burst for a snappy first paint.
        [0, 120, 300, 600, 1000, 1600].forEach(d => setTimeout(positionEntryButtons, d));
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
