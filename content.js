/* ============================================
   Content Script — UniProt 3D Feature Viewer
   Injected into UniProt entry pages.
   Adds "View in 3D" buttons next to PTM and
   Disease & Variants section headings.
   ============================================ */

(function () {
    'use strict';

    // ---- State ----
    let uniprotId = null;
    let currentMode = 'ptm'; // 'ptm' | 'variant'
    let featuresData = null;
    let variationData = null;
    let ptms = [];
    let ptmGroups = {};
    let variants = [];
    let activeConsequences = new Set();
    let activeProvenances = new Set();
    let activeDiseases = null; // null = no disease filter, Set = active diseases
    let scrapedDiseases = []; // disease names scraped from page h4 headings
    let modalEl = null;
    let overlayEl = null;
    let viewerInitialized = false;
    let currentTheme = 'auto'; // 'auto' | 'light' | 'dark'
    let lastPageContext = null; // 'entry' | 'variant-viewer' — tracks which page built the filters

    // ---- API URLs ----
    const ALPHAFOLD_API = (id) => `https://alphafold.ebi.ac.uk/api/prediction/${id}`;
    const FEATURES_URL  = (id) => `https://www.ebi.ac.uk/proteins/api/features/${id}`;
    const VARIATION_URL = (id) => `https://www.ebi.ac.uk/proteins/api/variation/${id}`;
    const PROTEOMICS_PTM_URL = (id) => `https://www.ebi.ac.uk/proteins/api/proteomics-ptm/${id}`;

    // ---- Icons ----
    const ICON_3D = `<svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"><path d="M12 2L2 7l10 5 10-5-10-5z"/><path d="M2 17l10 5 10-5"/><path d="M2 12l10 5 10-5"/></svg>`;
    const ICON_RESET = `<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M3 12a9 9 0 1 0 9-9 9.75 9.75 0 0 0-6.74 2.74L3 8"/><path d="M3 3v5h5"/></svg>`;

    const ICON_CAM   = `<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M23 19a2 2 0 0 1-2 2H3a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h4l2-3h6l2 3h4a2 2 0 0 1 2 2z"/><circle cx="12" cy="13" r="4"/></svg>`;
    const ICON_THEME = `<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/></svg>`;
    const ICON_COPY = `<svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="9" y="9" width="13" height="13" rx="2" ry="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg>`;

    // ---- State: currently displayed residue positions ----
    let displayedPositions = [];

    // ================================================================
    // 1) Extract UniProt ID from URL
    // ================================================================
    function getUniProtId() {
        const m = window.location.pathname.match(/\/uniprotkb\/([A-Za-z0-9_-]+)/);
        return m ? m[1].toUpperCase() : null;
    }

    // ================================================================
    // 2) Inject "View in 3D" buttons into the page
    // ================================================================
    function injectButtons() {
        // Use a MutationObserver to watch for lazy-loaded sections
        const observer = new MutationObserver(() => {
            tryInjectPTMButton();
            tryInjectVariantButton();
        });
        observer.observe(document.body, { childList: true, subtree: true });

        // Also use periodic polling as a fallback — UniProt's React lazy-loading
        // can sometimes evade MutationObserver callbacks
        const pollId = setInterval(() => {
            tryInjectPTMButton();
            tryInjectVariantButton();
            // Stop polling once both buttons are injected
            if (document.getElementById('ufv-btn-ptm') && document.getElementById('ufv-btn-variant')) {
                clearInterval(pollId);
            }
        }, 1500);

        // Stop polling after 2 minutes to avoid infinite loops
        setTimeout(() => clearInterval(pollId), 120000);

        // Try immediately
        tryInjectPTMButton();
        tryInjectVariantButton();

        // Also try on scroll — sections may lazy-load on scroll into view
        const scrollContainer = document.querySelector('.vJtX6') || window;
        scrollContainer.addEventListener('scroll', () => {
            tryInjectPTMButton();
            tryInjectVariantButton();
        }, { passive: true });
    }

    function tryInjectPTMButton() {
        // Check if button already exists and is still in the DOM
        const existing = document.getElementById('ufv-btn-ptm');
        if (existing && document.body.contains(existing)) return;

        // Look for the PTM / Processing section heading
        const heading = findSectionHeading('ptm_processing');
        if (!heading) return;

        const btn = createButton('ufv-btn-ptm', 'View PTMs in 3D', () => openModal('ptm'));
        // Insert next to the heading — prefer parent container for inline placement
        const container = heading.closest('.card__header, [class*="card__header"]') || heading;
        container.appendChild(btn);
        // console.log('\[UniProt 3D\] PTM button injected');
    }

    function tryInjectVariantButton() {
        const existing = document.getElementById('ufv-btn-variant');
        if (existing && document.body.contains(existing)) return;

        const heading = findSectionHeading('disease_variants');
        if (!heading) return;

        // Scrape disease names from h4 headings in the Disease & Variants section
        scrapedDiseases = scrapeDiseaseHeadings(heading);
        // console.log('\[UniProt 3D\] Scraped diseases from page:', scrapedDiseases);

        const btn = createButton('ufv-btn-variant', 'View Variants in 3D', () => openModal('variant'));
        const container = heading.closest('.card__header, [class*="card__header"]') || heading;
        container.appendChild(btn);
        // console.log('\[UniProt 3D\] Variant button injected');
    }

    /**
     * Scrape disease names from h4 elements within the Disease & Variants section.
     * UniProt renders diseases as <h4> with <a href="/diseases/DI-..."> children.
     * The h4s live inside card__content which is a sibling of card__header.
     */
    function scrapeDiseaseHeadings(sectionHeading) {
        const diseases = [];

        // Strategy 1: Find the section by its ID directly
        let section = document.getElementById('disease_variants');

        // Strategy 2: Walk up from the heading to find the broadest card container
        if (!section) {
            let el = sectionHeading;
            while (el && el !== document.body) {
                el = el.parentElement;
                if (el && (el.id === 'disease_variants' ||
                    el.tagName === 'SECTION' ||
                    (el.className && typeof el.className === 'string' &&
                     el.className.includes('card') && !el.className.includes('card__header')))) {
                    section = el;
                }
            }
        }

        if (!section) section = sectionHeading.parentElement?.parentElement || document.body;

        // Find all h4 elements within this section
        const h4s = section.querySelectorAll('h4');
        h4s.forEach(h4 => {
            const link = h4.querySelector('a[href*="/diseases/"]');
            if (link) {
                // Extract the full text, clean up surrounding quotes
                let name = h4.textContent.trim()
                    .replace(/^["'\u201C\u201D]|["'\u201C\u201D]$/g, '')
                    .trim();
                if (name) diseases.push(name);
            }
        });

        // console.log('\[UniProt 3D\] Disease scraping from section:', section?.id || section?.tagName, 'found:', diseases);
        return diseases;
    }

    /**
     * Find the section heading element to inject a button near.
     * Strategy: Look for the section by its ID, then find the h2 inside it.
     * Also tries text-based matching as a fallback, but only within
     * actual content sections (not sidebar navigation links).
     */
    function findSectionHeading(sectionId) {
        // Strategy 1: Find by section ID (most reliable)
        const section = document.getElementById(sectionId);
        if (section) {
            const h2 = section.querySelector('h2');
            if (h2) return h2;
            return section;
        }

        // Strategy 2: Find by heading text, but only in the main content area
        // (avoid matching sidebar navigation links)
        const textMap = {
            'ptm_processing': ['PTM/Processing', 'PTM / Processing'],
            'disease_variants': ['Disease & Variants', 'Disease and Variants']
        };
        const searchTerms = textMap[sectionId] || [];
        const mainContent = document.querySelector('main') || document.querySelector('[role="main"]') || document.body;
        const headings = mainContent.querySelectorAll('h2');
        for (const h of headings) {
            const text = h.textContent.trim();
            for (const term of searchTerms) {
                if (text.toLowerCase().includes(term.toLowerCase())) {
                    // Make sure this heading is inside a card/section, not the sidebar
                    const card = h.closest('section, .card, [class*="card"]');
                    if (card) return h;
                }
            }
        }

        return null;
    }

    function createButton(id, label, onClick) {
        const btn = document.createElement('button');
        btn.id = id;
        btn.className = 'ufv-3d-btn';
        btn.innerHTML = `${ICON_3D} ${label}`;
        btn.addEventListener('click', (e) => {
            e.preventDefault();
            e.stopPropagation();
            onClick();
        });
        return btn;
    }

    // ================================================================
    // 3) Modal — build once, reuse
    // ================================================================
    function buildModal() {
        if (overlayEl) return;

        overlayEl = document.createElement('div');
        overlayEl.className = 'ufv-overlay';
        overlayEl.style.display = 'none';

        overlayEl.innerHTML = `
        <div class="ufv-modal">
            <!-- Header -->
            <div class="ufv-modal-header">
                <div class="ufv-modal-title">
                    <span class="ufv-badge" id="ufv-id-badge">${uniprotId}</span>
                    <h2 id="ufv-modal-heading">3D Feature Viewer</h2>
                </div>
                <div class="ufv-modal-actions">
                    <button class="ufv-icon-btn" id="ufv-btn-theme" title="Toggle theme">${ICON_THEME}</button>
                    <button class="ufv-icon-btn" id="ufv-btn-reset" title="Reset view">${ICON_RESET}</button>
                    <button class="ufv-icon-btn" id="ufv-btn-screenshot" title="Screenshot">${ICON_CAM}</button>
                    <button class="ufv-close-btn" id="ufv-close" title="Close">×</button>
                </div>
            </div>

            <!-- Body -->
            <div class="ufv-body">
                <!-- 3D viewer -->
                <div class="ufv-viewer-wrap">
                    <div class="ufv-viewer" id="ufv-mol-viewer"></div>

                    <!-- Loading -->
                    <div class="ufv-loading" id="ufv-loading">
                        <div class="ufv-spinner"></div>
                        <div class="ufv-loading-text" id="ufv-loading-text">Loading structure...</div>
                    </div>

                    <!-- Tooltip -->
                    <div class="ufv-tooltip" id="ufv-tooltip">
                        <div class="ufv-tooltip-hdr" id="ufv-tooltip-hdr"></div>
                        <div class="ufv-tooltip-body" id="ufv-tooltip-body"></div>
                    </div>
                </div>

                <!-- Side panel -->
                <div class="ufv-side">
                    <!-- PTM filters -->
                    <div id="ufv-ptm-panel" class="ufv-filter-scroll">
                        <div class="ufv-panel-hdr">
                            <h3>PTM Types</h3>
                            <div class="ufv-panel-actions">
                                <button class="ufv-sm-btn" id="ufv-ptm-all">All</button>
                                <button class="ufv-sm-btn" id="ufv-ptm-none">None</button>
                            </div>
                        </div>
                        <div id="ufv-ptm-list"></div>
                    </div>

                    <!-- Variant filters -->
                    <div id="ufv-var-panel" class="ufv-filter-scroll ufv-hidden">
                        <div id="ufv-dis-section" class="ufv-hidden">
                            <div class="ufv-section-title">
                                <span>Disease</span>
                                <div class="ufv-section-actions">
                                    <button class="ufv-section-btn" id="ufv-dis-all">All</button>
                                    <button class="ufv-section-btn" id="ufv-dis-none">None</button>
                                </div>
                            </div>
                            <div id="ufv-dis-list"></div>
                        </div>
                        <div class="ufv-collapsible" id="ufv-prov-section">
                            <div class="ufv-collapsible-hdr" id="ufv-prov-toggle">
                                <span class="ufv-collapsible-chevron">&#9654;</span>
                                <span>Provenance</span>
                                <div class="ufv-section-actions">
                                    <button class="ufv-section-btn" id="ufv-prov-all">All</button>
                                    <button class="ufv-section-btn" id="ufv-prov-none">None</button>
                                </div>
                            </div>
                            <div class="ufv-collapsible-body ufv-collapsed" id="ufv-prov-body">
                                <div id="ufv-prov-list"></div>
                            </div>
                        </div>
                        <div class="ufv-collapsible" id="ufv-cons-section">
                            <div class="ufv-collapsible-hdr" id="ufv-cons-toggle">
                                <span class="ufv-collapsible-chevron">&#9654;</span>
                                <span>Consequence</span>
                                <div class="ufv-section-actions">
                                    <button class="ufv-section-btn" id="ufv-cons-all">All</button>
                                    <button class="ufv-section-btn" id="ufv-cons-none">None</button>
                                </div>
                            </div>
                            <div class="ufv-collapsible-body ufv-collapsed" id="ufv-cons-body">
                                <div id="ufv-cons-list"></div>
                            </div>
                        </div>
                    </div>

                    <!-- Footer -->
                    <div class="ufv-panel-footer">
                        <span class="ufv-count-text" id="ufv-count-text">—</span>
                        <button class="ufv-copy-btn" id="ufv-btn-copy" title="Copy residue IDs to clipboard">${ICON_COPY} Copy IDs</button>
                    </div>

                    <!-- Click details -->
                    <div class="ufv-details" id="ufv-details">
                        <div class="ufv-details-hdr">
                            <h4 id="ufv-details-title">Details</h4>
                            <button class="ufv-details-close" id="ufv-details-close">×</button>
                        </div>
                        <div class="ufv-details-body" id="ufv-details-body"></div>
                    </div>
                </div>
            </div>
        </div>`;

        document.body.appendChild(overlayEl);
        modalEl = overlayEl.querySelector('.ufv-modal');

        // --- Bind events ---
        document.getElementById('ufv-close').addEventListener('click', closeModal);
        overlayEl.addEventListener('click', (e) => { if (e.target === overlayEl) closeModal(); });
        document.addEventListener('keydown', (e) => { if (e.key === 'Escape') closeModal(); });

        // Theme toggle
        document.getElementById('ufv-btn-theme').addEventListener('click', cycleTheme);

        document.getElementById('ufv-btn-reset').addEventListener('click', () => {
            StructureViewer.resetView();
            if (viewerInitialized) applyMode();
            document.getElementById('ufv-details').classList.remove('show');
        });
        document.getElementById('ufv-btn-screenshot').addEventListener('click', () => StructureViewer.screenshot());

        // Copy residue IDs button
        document.getElementById('ufv-btn-copy').addEventListener('click', () => {
            if (!displayedPositions.length) return;
            const text = displayedPositions.join(', ');
            navigator.clipboard.writeText(text).then(() => {
                const btn = document.getElementById('ufv-btn-copy');
                const orig = btn.textContent;
                btn.textContent = '✓ Copied!';
                btn.classList.add('ufv-copied');
                setTimeout(() => {
                    btn.innerHTML = `${ICON_COPY} Copy IDs`;
                    btn.classList.remove('ufv-copied');
                }, 1500);
            }).catch(() => {
                // Fallback for older browsers
                const ta = document.createElement('textarea');
                ta.value = text;
                ta.style.cssText = 'position:fixed;opacity:0';
                document.body.appendChild(ta);
                ta.select();
                document.execCommand('copy');
                document.body.removeChild(ta);
            });
        });

        // PTM all/none
        document.getElementById('ufv-ptm-all').addEventListener('click', () => ptmSetAll(true));
        document.getElementById('ufv-ptm-none').addEventListener('click', () => ptmSetAll(false));

        // Variant consequence all/none
        document.getElementById('ufv-cons-all').addEventListener('click', () => varSectionSetAll('consequence', true));
        document.getElementById('ufv-cons-none').addEventListener('click', () => varSectionSetAll('consequence', false));

        // Variant provenance all/none
        document.getElementById('ufv-prov-all').addEventListener('click', () => varSectionSetAll('provenance', true));
        document.getElementById('ufv-prov-none').addEventListener('click', () => varSectionSetAll('provenance', false));

        // Disease all/none
        document.getElementById('ufv-dis-all').addEventListener('click', () => varSectionSetAll('disease', true));
        document.getElementById('ufv-dis-none').addEventListener('click', () => varSectionSetAll('disease', false));

        // Collapsible toggles for Provenance and Consequence
        document.getElementById('ufv-prov-toggle').addEventListener('click', (e) => {
            if (e.target.closest('.ufv-section-btn')) return; // Don't toggle when clicking All/None
            toggleCollapsible('ufv-prov-body', 'ufv-prov-toggle');
        });
        document.getElementById('ufv-cons-toggle').addEventListener('click', (e) => {
            if (e.target.closest('.ufv-section-btn')) return;
            toggleCollapsible('ufv-cons-body', 'ufv-cons-toggle');
        });

        // Details close
        document.getElementById('ufv-details-close').addEventListener('click', () => {
            document.getElementById('ufv-details').classList.remove('show');
        });
    }

    // ================================================================
    // 4) Open / Close modal
    // ================================================================
    async function openModal(mode) {
        buildModal();
        currentMode = mode;
        overlayEl.style.display = 'flex';
        document.body.style.overflow = 'hidden';

        // Update heading based on mode
        const heading = document.getElementById('ufv-modal-heading');
        if (heading) {
            heading.textContent = mode === 'ptm' ? 'PTM Viewer' : 'Variant Viewer';
        }

        // Keep the ID badge in sync with the currently viewed protein
        // (the modal is built once and reused across SPA navigations)
        const idBadge = document.getElementById('ufv-id-badge');
        if (idBadge) idBadge.textContent = uniprotId;

        // Show correct panel
        document.getElementById('ufv-ptm-panel').classList.toggle('ufv-hidden', mode !== 'ptm');
        document.getElementById('ufv-var-panel').classList.toggle('ufv-hidden', mode !== 'variant');
        document.getElementById('ufv-details').classList.remove('show');

        const currentPageContext = isVariantViewerPage() ? 'variant-viewer' : 'entry';

        if (!viewerInitialized) {
            lastPageContext = currentPageContext;
            await loadAllData();
        } else if (lastPageContext !== currentPageContext) {
            // Page context changed (entry ↔ variant-viewer)
            // Rebuild filters to pick up correct disease source and accordion state
            // console.log('\[UniProt 3D\] Page context changed:', lastPageContext, '→', currentPageContext);
            lastPageContext = currentPageContext;
            buildVariantFilters();
            applyMode();
        } else {
            applyMode();
        }
    }

    function closeModal() {
        if (!overlayEl) return;
        overlayEl.style.display = 'none';
        document.body.style.overflow = '';
    }

    // ================================================================
    // 5) Discover AlphaFold PDB URL via API
    // ================================================================
    async function getAlphaFoldPdbUrl(id) {
        try {
            const res = await fetch(ALPHAFOLD_API(id), {
                headers: { 'Accept': 'application/json' }
            });
            if (res.ok) {
                const data = await res.json();
                // API returns an array; pick the first canonical entry
                const entry = Array.isArray(data)
                    ? data.find(e => e.uniprotAccession === id) || data[0]
                    : data;
                if (entry && entry.pdbUrl) {
                    // console.log('\[UniProt 3D\] AlphaFold PDB URL from API:', entry.pdbUrl);
                    return entry.pdbUrl;
                }
            }
        } catch (e) {
            console.warn('[UniProt 3D] AlphaFold API lookup failed, trying fallback:', e);
        }

        // Fallback: try common versions in descending order
        for (const ver of [6, 5, 4, 3]) {
            const url = `https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v${ver}.pdb`;
            try {
                const probe = await fetch(url, { method: 'HEAD' });
                if (probe.ok) {
                    console.log(`[UniProt 3D] AlphaFold fallback hit v${ver}`);
                    return url;
                }
            } catch (_) { /* continue */ }
        }

        throw new Error('No AlphaFold structure found for ' + id);
    }

    // ================================================================
    // 6) Load data + structure
    // ================================================================
    async function loadAllData() {
        const loading = document.getElementById('ufv-loading');
        const loadText = document.getElementById('ufv-loading-text');
        loading.classList.remove('hidden');

        try {
            loadText.textContent = 'Fetching protein annotations...';

            const [featRes, varRes, ptmExRes] = await Promise.all([
                fetch(FEATURES_URL(uniprotId), { headers: { 'Accept': 'application/json' } }),
                fetch(VARIATION_URL(uniprotId), { headers: { 'Accept': 'application/json' } }),
                fetch(PROTEOMICS_PTM_URL(uniprotId), { headers: { 'Accept': 'application/json' } }).catch(() => null)
            ]);

            if (!featRes.ok) throw new Error(`Features API: ${featRes.status}`);
            featuresData = await featRes.json();
            variationData = varRes.ok ? await varRes.json() : { features: [] };
            const ptmExData = (ptmExRes && ptmExRes.ok) ? await ptmExRes.json() : null;

            // Process standard PTMs
            ptms = DataProcessor.extractPTMs(featuresData);

            // Merge large-scale PTMs from proteomics-ptm API
            if (ptmExData) {
                const lsPtms = DataProcessor.extractProteomicsPTMs(ptmExData);
                // Only add if not already present at that position with same category
                const existingKeys = new Set(ptms.map(p => `${p.position}-${p.category}`));
                lsPtms.forEach(lsp => {
                    if (!existingKeys.has(`${lsp.position}-${lsp.category}`)) {
                        ptms.push(lsp);
                        existingKeys.add(`${lsp.position}-${lsp.category}`);
                    }
                });
                console.log(`[UniProt 3D] Merged ${lsPtms.length} large-scale PTMs`);
            }

            ptmGroups = DataProcessor.groupPTMsByCategory(ptms);
            variants = DataProcessor.extractVariants(variationData);

            // Build filter UI (sets activeConsequences, activeProvenances, activeDiseases)
            buildPTMFilters();
            buildVariantFilters();

            // Load 3D structure — discover URL from AlphaFold API
            loadText.textContent = 'Loading AlphaFold structure...';
            StructureViewer.init(document.getElementById('ufv-mol-viewer'));
            StructureViewer.hoverCb = onHover;
            StructureViewer.clickCb = onClick;
            StructureViewer.dblClickCb = () => {
                StructureViewer.resetView();
                applyMode();
                document.getElementById('ufv-details').classList.remove('show');
            };

            const pdbUrl = await getAlphaFoldPdbUrl(uniprotId);
            await StructureViewer.loadStructure(pdbUrl);
            viewerInitialized = true;

            loading.classList.add('hidden');
            applyMode();

        } catch (err) {
            console.error('[UniProt 3D] Error:', err);
            loadText.textContent = `Error: ${err.message}`;
        }
    }

    // ================================================================
    // 6) Mode switching
    // ================================================================
    function switchMode(mode) {
        currentMode = mode;
        document.getElementById('ufv-ptm-panel').classList.toggle('ufv-hidden', mode !== 'ptm');
        document.getElementById('ufv-var-panel').classList.toggle('ufv-hidden', mode !== 'variant');
        document.getElementById('ufv-details').classList.remove('show');
        if (viewerInitialized) applyMode();
    }

    function applyMode() {
        if (currentMode === 'ptm') {
            const n = StructureViewer.showPTMs(ptms, ptmGroups);
            document.getElementById('ufv-count-text').textContent = `${n} PTM${n !== 1 ? 's' : ''} displayed`;
            // Collect displayed PTM positions
            const posSet = new Set();
            ptms.forEach(p => {
                const g = ptmGroups[p.category];
                if (g && g.visible) {
                    posSet.add(p.position);
                    if (p.endPosition && p.endPosition !== p.position) posSet.add(p.endPosition);
                }
            });
            displayedPositions = Array.from(posSet).sort((a, b) => a - b);
        } else {
            const filtered = DataProcessor.filterVariants(variants, activeConsequences, activeProvenances, activeDiseases);
            const r = StructureViewer.showVariants(filtered);
            document.getElementById('ufv-count-text').textContent =
                `${r.varCount} variant${r.varCount !== 1 ? 's' : ''} at ${r.posCount} position${r.posCount !== 1 ? 's' : ''}`;
            // Collect displayed variant positions
            const posSet = new Set();
            filtered.forEach(v => posSet.add(v.position));
            displayedPositions = Array.from(posSet).sort((a, b) => a - b);
        }
    }

    // ================================================================
    // 7) Build filter UI
    // ================================================================
    function buildPTMFilters() {
        const list = document.getElementById('ufv-ptm-list');
        list.innerHTML = '';

        const sorted = Object.entries(ptmGroups).sort((a, b) => b[1].items.length - a[1].items.length);
        sorted.forEach(([cat, group]) => {
            list.appendChild(makeFilterItem(cat, group.color, group.items.length, group.visible, (checked) => {
                group.visible = checked;
                applyMode();
            }));
        });
    }

    function buildVariantFilters() {
        const onVariantPage = isVariantViewerPage();

        // Provenance
        const provList = document.getElementById('ufv-prov-list');
        provList.innerHTML = '';
        const ps = DataProcessor.getProvenanceSummary(variants);
        activeProvenances = new Set(Object.keys(ps));
        Object.entries(ps).forEach(([cat, d]) => {
            provList.appendChild(makeFilterItem(cat, d.color, d.count, true, (checked) => {
                checked ? activeProvenances.add(cat) : activeProvenances.delete(cat);
                applyMode();
            }));
        });

        // Consequence
        const consList = document.getElementById('ufv-cons-list');
        consList.innerHTML = '';
        const cs = DataProcessor.getConsequenceSummary(variants);
        activeConsequences = new Set(Object.keys(cs));
        Object.entries(cs).forEach(([cat, d]) => {
            consList.appendChild(makeFilterItem(cat, d.color, d.count, true, (checked) => {
                checked ? activeConsequences.add(cat) : activeConsequences.delete(cat);
                applyMode();
            }));
        });

        // On variant-viewer page: expand Provenance and Consequence by default
        // On entry page: collapse them
        const provBody = document.getElementById('ufv-prov-body');
        const consBody = document.getElementById('ufv-cons-body');
        const provToggle = document.getElementById('ufv-prov-toggle');
        const consToggle = document.getElementById('ufv-cons-toggle');
        if (onVariantPage) {
            if (provBody) provBody.classList.remove('ufv-collapsed');
            if (consBody) consBody.classList.remove('ufv-collapsed');
            if (provToggle) provToggle.querySelector('.ufv-collapsible-chevron').innerHTML = '&#9660;';
            if (consToggle) consToggle.querySelector('.ufv-collapsible-chevron').innerHTML = '&#9660;';
        } else {
            if (provBody) provBody.classList.add('ufv-collapsed');
            if (consBody) consBody.classList.add('ufv-collapsed');
            if (provToggle) provToggle.querySelector('.ufv-collapsible-chevron').innerHTML = '&#9654;';
            if (consToggle) consToggle.querySelector('.ufv-collapsible-chevron').innerHTML = '&#9654;';
        }

        // Disease section
        const disList = document.getElementById('ufv-dis-list');
        disList.innerHTML = '';
        const disSection = document.getElementById('ufv-dis-section');

        // On variant-viewer page: use ALL diseases from API data
        // On entry page: use scraped h4 disease headings
        const ds = DataProcessor.getDiseaseSummary(variants);
        const useDiseases = onVariantPage
            ? Object.keys(ds).filter(k => k !== 'Unclassified')
            : scrapedDiseases;

        if (useDiseases.length > 0 && disSection) {
            disSection.classList.remove('ufv-hidden');

            const DISEASE_COLORS = [
                '#ef5350', '#42a5f5', '#ab47bc', '#66bb6a', '#ffa726',
                '#26c6da', '#ec407a', '#7e57c2', '#5c6bc0', '#29b6f6',
                '#8d6e63', '#78909c', '#d4e157', '#ff7043',
            ];
            activeDiseases = new Set();
            let colorIdx = 0;

            useDiseases.forEach(diseaseName => {
                let matchedKey = null;
                let count = 0;

                // On variant-viewer, diseaseName IS the key from API
                if (onVariantPage && ds[diseaseName]) {
                    matchedKey = diseaseName;
                    count = ds[diseaseName].count;
                } else {
                    // Scraped name: try partial match against API data
                    for (const [key, val] of Object.entries(ds)) {
                        if (key.toLowerCase().includes(diseaseName.toLowerCase().substring(0, 20)) ||
                            diseaseName.toLowerCase().includes(key.toLowerCase().substring(0, 20))) {
                            matchedKey = key;
                            count = val.count;
                            break;
                        }
                    }
                    if (!matchedKey && ds[diseaseName]) {
                        matchedKey = diseaseName;
                        count = ds[diseaseName].count;
                    }
                    if (!matchedKey) {
                        count = variants.filter(v =>
                            (v.diseases || []).some(d =>
                                d.toLowerCase().includes(diseaseName.toLowerCase().substring(0, 15))
                            )
                        ).length;
                        matchedKey = diseaseName;
                    }
                }

                const color = DISEASE_COLORS[colorIdx % DISEASE_COLORS.length];
                colorIdx++;
                activeDiseases.add(matchedKey);

                disList.appendChild(makeFilterItem(diseaseName, color, count, true, (checked) => {
                    checked ? activeDiseases.add(matchedKey) : activeDiseases.delete(matchedKey);
                    applyMode();
                }));
            });

            // Add "Unclassified" if there are variants with no disease
            const noDisease = variants.filter(v => !v.diseases || v.diseases.length === 0).length;
            if (noDisease > 0) {
                activeDiseases.add('Unclassified');
                disList.appendChild(makeFilterItem('Unclassified', '#9e9e9e', noDisease, true, (checked) => {
                    checked ? activeDiseases.add('Unclassified') : activeDiseases.delete('Unclassified');
                    applyMode();
                }));
            }
        } else {
            if (disSection) disSection.classList.add('ufv-hidden');
            activeDiseases = null;
        }
    }

    function toggleCollapsible(bodyId, toggleId) {
        const body = document.getElementById(bodyId);
        const toggle = document.getElementById(toggleId);
        if (!body || !toggle) return;
        const isCollapsed = body.classList.toggle('ufv-collapsed');
        const chevron = toggle.querySelector('.ufv-collapsible-chevron');
        if (chevron) {
            chevron.innerHTML = isCollapsed ? '&#9654;' : '&#9660;';
        }
    }
    /**
     * Build a detail row using safe DOM construction (no innerHTML).
     * All values are set via textContent to prevent XSS.
     */
    function makeDetailRow(label, value, color) {
        const row = document.createElement('div');
        row.className = 'ufv-detail-row';

        const lbl = document.createElement('span');
        lbl.className = 'ufv-detail-lbl';
        lbl.textContent = label;

        const val = document.createElement('span');
        val.className = 'ufv-detail-val';
        val.textContent = value || '';
        if (color) val.style.color = color;

        row.appendChild(lbl);
        row.appendChild(val);
        return row;
    }

    /**
     * Validate a color string is a safe hex color to prevent style injection.
     * Returns the color if valid, or a default fallback.
     */
    function sanitizeColor(color) {
        if (!color) return '';
        // Allow only hex colors (#abc, #aabbcc, #aabbccdd)
        if (/^#[0-9a-fA-F]{3,8}$/.test(color)) return color;
        // Allow named CSS colors (simple alphanumeric)
        if (/^[a-zA-Z]{3,20}$/.test(color)) return color;
        return '#888';
    }

    function makeFilterItem(label, color, count, checked, onChange) {
        const el = document.createElement('label');
        el.className = 'ufv-filter-item';

        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.checked = checked;
        cb.addEventListener('change', () => onChange(cb.checked));

        const dot = document.createElement('span');
        dot.className = 'ufv-dot';
        dot.style.backgroundColor = color;

        const lbl = document.createElement('span');
        lbl.className = 'ufv-filter-label';
        lbl.textContent = label;

        const cnt = document.createElement('span');
        cnt.className = 'ufv-filter-count';
        cnt.textContent = count;

        el.append(cb, dot, lbl, cnt);
        return el;
    }

    function ptmSetAll(select) {
        Object.values(ptmGroups).forEach(g => g.visible = select);
        document.querySelectorAll('#ufv-ptm-list input[type="checkbox"]').forEach(cb => cb.checked = select);
        applyMode();
    }

    function varSectionSetAll(section, select) {
        if (section === 'consequence') {
            const cs = DataProcessor.getConsequenceSummary(variants);
            if (select) {
                activeConsequences = new Set(Object.keys(cs));
            } else {
                activeConsequences.clear();
            }
            document.querySelectorAll('#ufv-cons-list input[type="checkbox"]').forEach(cb => cb.checked = select);
        } else if (section === 'provenance') {
            const ps = DataProcessor.getProvenanceSummary(variants);
            if (select) {
                activeProvenances = new Set(Object.keys(ps));
            } else {
                activeProvenances.clear();
            }
            document.querySelectorAll('#ufv-prov-list input[type="checkbox"]').forEach(cb => cb.checked = select);
        } else if (section === 'disease') {
            const ds = DataProcessor.getDiseaseSummary(variants);
            if (select) {
                activeDiseases = new Set(Object.keys(ds));
            } else {
                activeDiseases = new Set();
            }
            document.querySelectorAll('#ufv-dis-list input[type="checkbox"]').forEach(cb => cb.checked = select);
        }
        applyMode();
    }

    // ================================================================
    // 7b) Theme cycling
    // ================================================================
    function cycleTheme() {
        const cycle = ['auto', 'light', 'dark'];
        const idx = cycle.indexOf(currentTheme);
        currentTheme = cycle[(idx + 1) % cycle.length];

        const modal = document.querySelector('.ufv-modal');
        if (!modal) return;

        modal.classList.remove('ufv-light', 'ufv-dark');
        if (currentTheme === 'light') {
            modal.classList.add('ufv-light');
        } else if (currentTheme === 'dark') {
            modal.classList.add('ufv-dark');
        }
        // 'auto' = no class, follows prefers-color-scheme

        // Update viewer background to match theme
        const isDark = currentTheme === 'dark' ||
            (currentTheme === 'auto' && window.matchMedia('(prefers-color-scheme: dark)').matches);
        const bgColor = isDark ? '#0c111b' : '#f0f2f5';
        if (StructureViewer.viewer) {
            StructureViewer.viewer.setBackgroundColor(bgColor);
            StructureViewer.viewer.render();
        }

        // Update button title
        const themeBtn = document.getElementById('ufv-btn-theme');
        if (themeBtn) {
            const labels = { auto: 'Theme: Auto', light: 'Theme: Light', dark: 'Theme: Dark' };
            themeBtn.title = labels[currentTheme];
        }
    }

    // ================================================================
    // 8) Hover / Click handlers
    // ================================================================
    function onHover(data, mode, event) {
        const tip = document.getElementById('ufv-tooltip');
        if (!data) {
            tip.classList.remove('show');
            return;
        }

        const hdr = document.getElementById('ufv-tooltip-hdr');
        const body = document.getElementById('ufv-tooltip-body');

        if (mode === 'ptm') {
            hdr.textContent = `Pos ${data.position} — ${data.category}`;
            body.textContent = data.description;
        } else {
            const vv = data.variants || [data];
            hdr.textContent = `Pos ${data.position} — ${data.topConsequence}`;
            const list = vv.slice(0, 3).map(v => `${v.wildType}${v.position}${v.mutant}`).join(', ');
            body.textContent = list + (vv.length > 3 ? ` +${vv.length - 3} more` : '');
        }

        if (event) {
            const wrap = document.querySelector('.ufv-viewer-wrap');
            const rect = wrap.getBoundingClientRect();
            const x = (event.clientX || event.pageX) - rect.left + 14;
            const y = (event.clientY || event.pageY) - rect.top + 14;
            tip.style.left = Math.min(x, rect.width - 280) + 'px';
            tip.style.top  = Math.min(y, rect.height - 60) + 'px';
        }

        tip.classList.add('show');
    }

    function onClick(data, mode) {
        const panel = document.getElementById('ufv-details');
        const title = document.getElementById('ufv-details-title');
        const body  = document.getElementById('ufv-details-body');

        // Clear previous content safely
        body.textContent = '';

        if (mode === 'ptm') {
            title.textContent = `PTM at Position ${Number(data.position) || '?'}`;
            body.appendChild(makeDetailRow('Position', String(data.position)));
            body.appendChild(makeDetailRow('Type', data.type));
            body.appendChild(makeDetailRow('Category', data.category));
            body.appendChild(makeDetailRow('Description', data.description));
        } else {
            const vv = data.variants || [data];
            title.textContent = `${vv.length} Variant${vv.length > 1 ? 's' : ''} at Position ${Number(data.position) || '?'}`;

            vv.slice(0, 12).forEach(v => {
                const item = document.createElement('div');
                item.className = 'ufv-variant-item';

                // Mutation tag
                const mutRow = makeDetailRow('Mutation', '');
                const mutTag = document.createElement('span');
                mutTag.className = 'ufv-mutation-tag';
                mutTag.textContent = `${v.wildType || '?'}${v.position}${v.mutant || '?'}`;
                mutRow.querySelector('.ufv-detail-val').replaceWith(mutTag);
                item.appendChild(mutRow);

                // Consequence (colored)
                item.appendChild(makeDetailRow('Consequence', v.consequence, sanitizeColor(v.consequenceColor)));

                // Provenance (colored)
                item.appendChild(makeDetailRow('Provenance', v.provenance, sanitizeColor(v.provenanceColor)));

                // Description (optional)
                if (v.description) {
                    item.appendChild(makeDetailRow('Details', v.description));
                }

                body.appendChild(item);
            });

            if (vv.length > 12) {
                const more = document.createElement('p');
                more.style.cssText = 'color:#5a6478;padding:6px 0;font-style:italic';
                more.textContent = `+${vv.length - 12} more`;
                body.appendChild(more);
            }
        }

        panel.classList.add('show');
        // Zoom in with stick-and-ball
        StructureViewer.focusResidue(data.position);
    }

    // ================================================================
    // Variant-viewer page button injection
    // ================================================================
    function isVariantViewerPage() {
        return window.location.pathname.includes('/variant-viewer');
    }

    function tryInjectVariantViewerButton() {
        const existing = document.getElementById('ufv-btn-variant-page');
        if (existing && document.body.contains(existing)) return;

        // On variant-viewer page, find the toolbar or any container near the top
        // Try to find the protvista or variant viewer header area
        const targets = [
            document.querySelector('[class*="variant-viewer"] h2'),
            document.querySelector('[class*="variant"] [class*="header"]'),
            document.querySelector('main h2'),
            document.querySelector('main h1'),
            document.querySelector('[class*="toolbar"]'),
        ].filter(Boolean);

        let anchor = targets[0];

        // If no anchor found, create a floating button
        if (!anchor) {
            // Insert button as a fixed floating element
            const btn = createButton('ufv-btn-variant-page', 'View Variants in 3D', () => openModal('variant'));
            btn.style.cssText = 'position:fixed;top:80px;right:20px;z-index:9999;box-shadow:0 4px 16px rgba(0,0,0,0.3);';
            document.body.appendChild(btn);
            // console.log('\[UniProt 3D\] Variant-viewer floating button injected');
            return;
        }

        const btn = createButton('ufv-btn-variant-page', 'View Variants in 3D', () => openModal('variant'));
        anchor.appendChild(btn);
        // console.log('\[UniProt 3D\] Variant-viewer page button injected');
    }

    // ================================================================
    // Boot
    // ================================================================
    console.log('[UniProt 3D Feature Viewer] Content script loaded on:', window.location.href);

    uniprotId = getUniProtId();
    if (uniprotId) {
        console.log(`[UniProt 3D Feature Viewer] Detected protein: ${uniprotId}`);

        // Initial injection based on current page
        handlePageType();

        // --- SPA navigation detection ---
        // UniProt uses React with history.pushState for tab switching.
        // We use three strategies to reliably detect URL changes:

        let lastUrl = window.location.href;

        function onUrlChange() {
            const currentUrl = window.location.href;
            if (currentUrl === lastUrl) return;
            const oldUrl = lastUrl;
            lastUrl = currentUrl;
            // console.log('\[UniProt 3D\] SPA navigation:', oldUrl, '->', currentUrl);

            // If the user navigated to a different protein, the cached
            // structure/annotations belong to the old accession. Drop them so
            // the next openModal() reloads fresh data instead of showing stale.
            const newId = getUniProtId();
            if (newId && newId !== uniprotId) {
                // console.log('\[UniProt 3D\] Protein changed:', uniprotId, '->', newId);
                uniprotId = newId;
                resetProteinState();
            }

            handlePageType();
        }

        function resetProteinState() {
            featuresData = null;
            variationData = null;
            ptms = [];
            ptmGroups = {};
            variants = [];
            activeConsequences = new Set();
            activeProvenances = new Set();
            activeDiseases = null;
            scrapedDiseases = [];
            displayedPositions = [];
            viewerInitialized = false;
            lastPageContext = null;

            // The modal (if built) is showing the previous protein — close it.
            // It will reload from scratch the next time a "View in 3D" button
            // is clicked for the new protein.
            if (overlayEl && overlayEl.style.display !== 'none') {
                closeModal();
            }
        }

        function handlePageType() {
            if (isVariantViewerPage()) {
                // console.log('\[UniProt 3D\] Detected variant-viewer page');
                // Poll to inject — the DOM may not be ready yet
                let attempts = 0;
                const pollId = setInterval(() => {
                    tryInjectVariantViewerButton();
                    attempts++;
                    if (document.getElementById('ufv-btn-variant-page') || attempts > 40) {
                        clearInterval(pollId);
                    }
                }, 1500);
                tryInjectVariantViewerButton();
            } else {
                // Entry page or other — inject PTM + Variant buttons
                injectButtons();
            }
        }

        // Strategy 1: Intercept history.pushState and history.replaceState
        const origPushState = history.pushState;
        const origReplaceState = history.replaceState;
        history.pushState = function() {
            origPushState.apply(this, arguments);
            onUrlChange();
        };
        history.replaceState = function() {
            origReplaceState.apply(this, arguments);
            onUrlChange();
        };

        // Strategy 2: Listen for popstate (back/forward buttons)
        window.addEventListener('popstate', onUrlChange);

        // Strategy 3: Polling fallback (catches any edge cases)
        setInterval(onUrlChange, 1000);

    } else {
        console.log('[UniProt 3D Feature Viewer] No UniProt ID found in URL');
    }

})();
