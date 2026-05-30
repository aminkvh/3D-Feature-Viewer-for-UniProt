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
        if (section) return section.querySelector('h2') || section;
        const textMap = {
            ptm_processing: ['PTM/Processing', 'PTM / Processing'],
            disease_variants: ['Disease & Variants', 'Disease and Variants'],
        };
        const main = document.querySelector('main') || document.querySelector('[role="main"]') || document.body;
        for (const h of main.querySelectorAll('h2')) {
            const text = h.textContent.trim().toLowerCase();
            if ((textMap[sectionId] || []).some(term => text.includes(term.toLowerCase()))) return h;
        }
        return null;
    }

    // Attach a section-header button without changing the section's vertical layout: the
    // button is absolutely positioned (see .ufv-3d-btn--anchored), so it can't reflow the
    // sections below it and disturb UniProt's scroll/IntersectionObserver active-slider.
    function placeAnchoredButton(heading, btn) {
        const container = heading.closest('.card__header, [class*="card__header"]') || heading;
        try {
            if (getComputedStyle(container).position === 'static') container.style.position = 'relative';
        } catch (_) {}
        btn.classList.add('ufv-3d-btn--anchored');
        container.appendChild(btn);
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
        if (!runtime.observer) {
            // Debounce: React batch-renders cause many rapid mutations; collapse them
            // into a single inject attempt so we never mutate the DOM mid-scroll.
            let _injectTimer = null;
            runtime.observer = new MutationObserver(() => {
                // Skip when both buttons are already present — prevents timer churn
                // during scroll-triggered React re-renders which would otherwise
                // interfere with UniProt's IntersectionObserver sidebar indicator.
                if (document.getElementById('ufv-btn-ptm')?.isConnected &&
                    document.getElementById('ufv-btn-variant')?.isConnected) return;
                if (_injectTimer) return;
                _injectTimer = setTimeout(() => {
                    _injectTimer = null;
                    tryInjectPTMButton();
                    tryInjectVariantButton();
                }, 150);
            });
            runtime.observer.observe(document.body, { childList: true, subtree: true });
        }
        // Re-inject when switching back to this tab (UniProt SPA re-renders on visibility)
        if (!runtime.visibilityBound) {
            runtime.visibilityBound = true;
            document.addEventListener('visibilitychange', () => {
                if (!document.hidden) {
                    tryInjectPTMButton();
                    tryInjectVariantButton();
                }
            });
        }
        if (!runtime.entryPoll) {
            let attempts = 0;
            runtime.entryPoll = setInterval(() => {
                tryInjectPTMButton();
                tryInjectVariantButton();
                attempts++;
                if ((document.getElementById('ufv-btn-ptm') && document.getElementById('ufv-btn-variant')) || attempts > 80) {
                    clearInterval(runtime.entryPoll);
                    runtime.entryPoll = null;
                }
            }, 1500);
        }
        // NOTE: scroll listener removed — MutationObserver already covers lazy-loaded
        // section re-renders, and the scroll listener was racing with UniProt's
        // IntersectionObserver that drives the sidebar active-indicator.
        tryInjectPTMButton();
        tryInjectVariantButton();
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
            }, 1500);
        }
        tryInjectVariantViewerButton();
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
            ['ufv-btn-ptm', 'ufv-btn-variant', 'ufv-btn-variant-page'].forEach(id => document.getElementById(id)?.remove());
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
    }

    return { boot, getUniProtId };
})();
