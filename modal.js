/* global UFVState, UFVApi, UFVExport, UFVAnalysis, DataProcessor, StructureViewer */
const UFVModal = (() => {
    'use strict';

    const ICON_3D = '<svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2"><path d="M12 2L2 7l10 5 10-5-10-5z"/><path d="M2 17l10 5 10-5"/><path d="M2 12l10 5 10-5"/></svg>';
    const ICON_COPY = '<svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg>';
    const ICON_THEME = '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/></svg>';
    const ICON_RESET = '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"><polyline points="23 4 23 10 17 10"/><path d="M20.49 15a9 9 0 1 1-2.12-9.36L23 10"/></svg>';
    const ICON_DOWNLOAD = '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>';
    const ICON_CAMERA = '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M23 19a2 2 0 0 1-2 2H3a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h4l2-3h6l2 3h4a2 2 0 0 1 2 2z"/><circle cx="12" cy="13" r="4"/></svg>';

    let overlayEl = null;
    let modalEl = null;
    let _hostScrollSaved = 0;
    let _hostScrollLock = false;

    // UniProt's entry page scrolls inside a CUSTOM container (an obfuscated-class div, e.g. .vJtX6), NOT the
    // window — window.scrollY stays 0. Find it dynamically (the scrollable ancestor of the page content)
    // rather than by a hard-coded class. (.Gsgt9 is the left-nav active-section INDICATOR, not a scroll
    // container — querying it for scrollTop was a no-op, which is why the background-scroll lock did nothing.)
    function findHostScrollContainer() {
        let el = document.getElementById('function') || document.querySelector('main section[id]') || document.querySelector('section[id]');
        while (el && el !== document.body && el !== document.documentElement) {
            const cs = window.getComputedStyle(el);
            if ((cs.overflowY === 'auto' || cs.overflowY === 'scroll') && el.scrollHeight > el.clientHeight + 4) return el;
            el = el.parentElement;
        }
        return null;
    }
    // Monotonic token: every open() and structure load captures the current value; if a newer
    // call has since bumped it, the older (superseded) flow bails instead of clobbering the UI.
    // This is what prevents the "variant modal shows the PTM header / nothing ever loads" race
    // when the user navigates between proteins or re-opens the modal mid-load.
    let _openSeq = 0;
    let _loadSeq = 0;
    // Constraint-pocket significance threshold (BH-FDR q) controlled by the sensitivity slider.
    // The analysis returns ALL candidates with q-values; this filters what's shown, so moving
    // the slider re-thresholds instantly without recomputing.
    let sensThreshold = 0.10;

    // The displayed constraint-pocket candidate set, filtered to q ≤ sensThreshold.
    function filteredPocketByPos() {
        const res = UFVState.state.analysis.prism;
        if (!res?.byPos) return null;
        const out = new Map();
        res.byPos.forEach((v, pos) => { if (v.q <= sensThreshold) out.set(pos, v); });
        return out;
    }

    function createButton(id, label, onClick) {
        const btn = document.createElement('button');
        btn.id = id;
        btn.className = 'ufv-3d-btn';
        btn.innerHTML = `${ICON_3D} ${label}`;
        btn.addEventListener('click', e => {
            e.preventDefault();
            e.stopPropagation();
            onClick();
        });
        return btn;
    }

    function build() {
        if (overlayEl) return;
        overlayEl = document.createElement('div');
        overlayEl.className = 'ufv-overlay';
        overlayEl.style.display = 'none';
        overlayEl.innerHTML = `
        <div class="ufv-modal">
            <div class="ufv-modal-header">
                <div class="ufv-modal-title">
                    <span class="ufv-badge" id="ufv-id-badge"></span>
                    <h2 id="ufv-modal-heading">3D Feature Viewer</h2>
                </div>
                <div class="ufv-structure-nav">
                    <button class="ufv-nav-btn" id="ufv-structure-prev">&#8592;</button>
                    <div class="ufv-cs" id="ufv-cs"><button class="ufv-cs-btn" id="ufv-cs-btn">Loading&#8230;</button><div class="ufv-cs-drop" id="ufv-cs-drop"></div></div>
                    <button class="ufv-nav-btn" id="ufv-structure-next">&#8594;</button>
                </div>
                <div class="ufv-modal-actions">
                    <button class="ufv-icon-btn" id="ufv-btn-theme" title="Toggle theme">${ICON_THEME}</button>
                    <button class="ufv-icon-btn" id="ufv-btn-reset" title="Reset view">${ICON_RESET}</button>
                    <button class="ufv-icon-btn" id="ufv-btn-protnlm" title="Toggle ProtNLM AI-predicted protein name" style="font-size:11px;font-weight:700;">AI</button>
                    <div class="ufv-dl-wrap" id="ufv-dl-wrap">
                        <button class="ufv-icon-btn" id="ufv-btn-export-pdb" title="Download">${ICON_DOWNLOAD}</button>
                        <div class="ufv-dl-menu" id="ufv-dl-menu">
                            <button class="ufv-dl-opt" id="ufv-dl-pdb">PDB file</button>
                            <button class="ufv-dl-opt" id="ufv-dl-csv">CSV annotation table</button>
                            <button class="ufv-dl-opt" id="ufv-dl-csv-pv">CSV + ProtVar predictions (slower)</button>
                            <button class="ufv-dl-opt" id="ufv-dl-pymol">PyMOL session (.pml)</button>
                            <button class="ufv-dl-opt" id="ufv-dl-vmd">VMD session (.vmd)</button>
                        </div>
                    </div>
                    <button class="ufv-icon-btn" id="ufv-btn-screenshot" title="Screenshot">${ICON_CAMERA}</button>
                    <button class="ufv-close-btn" id="ufv-close" title="Close">&#10005;</button>
                </div>
            </div>
            <div id="ufv-protnlm-banner" class="ufv-hidden" style="padding:6px 14px;font-size:12px;border-bottom:1px solid var(--ufv-border,#ddd);"></div>
            <div class="ufv-body">
                <div class="ufv-left-col">
                <div class="ufv-sequence-wrap" id="ufv-sequence-wrap"></div>
                <div class="ufv-viewer-wrap">
                    <div class="ufv-viewer" id="ufv-mol-viewer"></div>
                    <div class="ufv-legend" id="ufv-legend"></div>
                    <div class="ufv-loading" id="ufv-loading">
                        <div class="ufv-spinner"></div>
                        <div class="ufv-loading-text" id="ufv-loading-text">Loading...</div>
                    </div>
                    <div class="ufv-tooltip" id="ufv-tooltip">
                        <div class="ufv-tooltip-hdr" id="ufv-tooltip-hdr"></div>
                        <div class="ufv-tooltip-body" id="ufv-tooltip-body"></div>
                    </div>
                </div>
                </div>
                <div class="ufv-side">
                    <div class="ufv-view-section">
                        <div class="ufv-cm" id="ufv-cm">
                            <button class="ufv-cm-btn" id="ufv-cm-btn">Default</button>
                            <div class="ufv-cm-drop" id="ufv-cm-drop">
                                <div class="ufv-cm-opt selected" data-value="default">Default</div>
                                <div class="ufv-cm-opt" data-value="plddt">pLDDT confidence</div>
                                <div class="ufv-cm-opt" data-value="bfactor">Experimental B-factor</div>
                                <div class="ufv-cm-opt ufv-hidden" data-value="topology">Membrane topology</div>
                                <div class="ufv-cm-opt" data-value="hotspots">Pathogenic variant hotspots</div>
                                <div class="ufv-cm-opt" data-value="distantContacts">Contact-network centrality</div>
                                <div class="ufv-cm-opt" data-value="alphaMissense">AlphaMissense summary</div>
                                <div class="ufv-cm-opt" data-value="residueBurden">Recurrent phenotype residues</div>
                                <div class="ufv-cm-opt" data-value="prism">Burial-adjusted constraint clusters</div>
                            </div>
                        </div>
                        <div class="ufv-sens-slider ufv-hidden" id="ufv-sens-wrap">
                            <label for="ufv-sens-slider">Sensitivity (FDR q ≤ <span id="ufv-sens-q">0.10</span>)</label>
                            <input type="range" id="ufv-sens-slider" min="1" max="40" value="10">
                        </div>
                    </div>
                    <div id="ufv-ptm-panel" class="ufv-filter-scroll">
                        <div class="ufv-panel-hdr"><h3>PTM Types</h3><div class="ufv-panel-actions"><button class="ufv-sm-btn" id="ufv-ptm-all">All</button><button class="ufv-sm-btn" id="ufv-ptm-none">None</button><button class="ufv-sm-btn ufv-brush-btn" id="ufv-brush-ptm" title="Colour PTM spheres by PTM type">C</button></div></div>
                        <div id="ufv-ptm-list"></div>
                        <div class="ufv-collapsible ufv-hidden" id="ufv-sites-section-ptm"><div class="ufv-collapsible-hdr" id="ufv-sites-ptm-toggle"><span class="ufv-collapsible-chevron">&#9654;</span><span>Sites</span><div class="ufv-section-actions"><button class="ufv-section-btn" id="ufv-sites-ptm-all">All</button><button class="ufv-section-btn" id="ufv-sites-ptm-none">None</button><button class="ufv-section-btn ufv-brush-btn" id="ufv-brush-ptmsite" title="Colour spheres by functional site">C</button></div></div><div class="ufv-collapsible-body ufv-collapsed" id="ufv-sites-ptm-body"><div id="ufv-sites-ptm-list"></div></div></div>
                        <div class="ufv-collapsible ufv-hidden" id="ufv-ligands-section-ptm"><div class="ufv-collapsible-hdr" id="ufv-ligands-ptm-toggle"><span class="ufv-collapsible-chevron">&#9654;</span><span>Ligands</span><div class="ufv-section-actions"><label class="ufv-toggle-switch ufv-ions-toggle" title="Exclude water &amp; ions"><input type="checkbox" id="ufv-ligands-ions-ptm"><span class="ufv-toggle-slider"></span></label><button class="ufv-section-btn" id="ufv-ligands-ptm-all">All</button><button class="ufv-section-btn" id="ufv-ligands-ptm-none">None</button></div></div><div class="ufv-collapsible-body ufv-collapsed" id="ufv-ligands-ptm-body"><div id="ufv-ligands-ptm-list"></div></div></div>
                    </div>
                    <div id="ufv-var-panel" class="ufv-filter-scroll ufv-hidden">
                        <div id="ufv-dis-section" class="ufv-collapsible ufv-hidden"><div class="ufv-collapsible-hdr" id="ufv-dis-toggle"><span class="ufv-collapsible-chevron">&#9660;</span><span>Disease <span class="ufv-section-source">— HumanVar</span></span><div class="ufv-section-actions"><button class="ufv-section-btn" id="ufv-dis-all">All</button><button class="ufv-section-btn" id="ufv-dis-none">None</button><button class="ufv-section-btn ufv-brush-btn" id="ufv-brush-var" title="Colour variant spheres by disease">C</button></div></div><div class="ufv-collapsible-body" id="ufv-dis-body"><div id="ufv-dis-list"></div></div></div>
                        <div class="ufv-collapsible"><div class="ufv-collapsible-hdr" id="ufv-prov-toggle"><span class="ufv-collapsible-chevron">&#9654;</span><span>Provenance</span><div class="ufv-section-actions"><button class="ufv-section-btn" id="ufv-prov-all">All</button><button class="ufv-section-btn" id="ufv-prov-none">None</button><button class="ufv-section-btn ufv-brush-btn" id="ufv-brush-prov" title="Colour by consequence">C</button></div></div><div class="ufv-collapsible-body ufv-collapsed" id="ufv-prov-body"><div id="ufv-prov-list"></div></div></div>
                        <div class="ufv-collapsible"><div class="ufv-collapsible-hdr" id="ufv-cons-toggle"><span class="ufv-collapsible-chevron">&#9654;</span><span>Consequence</span><div class="ufv-section-actions"><button class="ufv-section-btn" id="ufv-cons-all">All</button><button class="ufv-section-btn" id="ufv-cons-none">None</button><button class="ufv-section-btn ufv-brush-btn" id="ufv-brush-dis" title="Colour by consequence">C</button></div></div><div class="ufv-collapsible-body ufv-collapsed" id="ufv-cons-body"><div id="ufv-cons-list"></div></div></div>
                        <div class="ufv-collapsible ufv-hidden" id="ufv-vptm-section"><div class="ufv-collapsible-hdr" id="ufv-vptm-toggle"><span class="ufv-collapsible-chevron">&#9654;</span><span>PTM sites</span><div class="ufv-section-actions"><button class="ufv-section-btn" id="ufv-vptm-all">All</button><button class="ufv-section-btn" id="ufv-vptm-none">None</button><button class="ufv-section-btn ufv-brush-btn" id="ufv-brush-vptm" title="Colour by PTM type">C</button></div></div><div class="ufv-collapsible-body ufv-collapsed" id="ufv-vptm-body"><div id="ufv-vptm-list"></div></div></div>
                        <div class="ufv-collapsible ufv-hidden" id="ufv-sites-section-var"><div class="ufv-collapsible-hdr" id="ufv-sites-var-toggle"><span class="ufv-collapsible-chevron">&#9654;</span><span>Sites</span><div class="ufv-section-actions"><button class="ufv-section-btn" id="ufv-sites-var-all">All</button><button class="ufv-section-btn" id="ufv-sites-var-none">None</button><button class="ufv-section-btn ufv-brush-btn" id="ufv-brush-vsite" title="Colour by consequence">C</button></div></div><div class="ufv-collapsible-body ufv-collapsed" id="ufv-sites-var-body"><div id="ufv-sites-var-list"></div></div></div>
                        <div class="ufv-collapsible ufv-hidden" id="ufv-ligands-section-var"><div class="ufv-collapsible-hdr" id="ufv-ligands-var-toggle"><span class="ufv-collapsible-chevron">&#9654;</span><span>Ligands</span><div class="ufv-section-actions"><label class="ufv-toggle-switch ufv-ions-toggle" title="Exclude water &amp; ions"><input type="checkbox" id="ufv-ligands-ions-var"><span class="ufv-toggle-slider"></span></label><button class="ufv-section-btn" id="ufv-ligands-var-all">All</button><button class="ufv-section-btn" id="ufv-ligands-var-none">None</button></div></div><div class="ufv-collapsible-body ufv-collapsed" id="ufv-ligands-var-body"><div id="ufv-ligands-var-list"></div></div></div>
                    </div>
                    <div id="ufv-feat-panel" class="ufv-filter-scroll ufv-hidden"></div>
                    <div id="ufv-dom-panel" class="ufv-filter-scroll ufv-hidden"></div>
                    <div class="ufv-panel-footer"><span class="ufv-count-text" id="ufv-count-text">-</span><div class="ufv-footer-actions"><div class="ufv-protein-nav" id="ufv-protein-nav"></div><button class="ufv-copy-btn" id="ufv-btn-copy">${ICON_COPY} Copy</button></div></div>
                    <div class="ufv-details" id="ufv-details"><div class="ufv-details-hdr"><h4 id="ufv-details-title">Details</h4><div class="ufv-details-hdr-actions"><label class="ufv-toggle-switch" id="ufv-sphere-toggle" title="Show/hide annotation spheres"><input type="checkbox" id="ufv-sphere-chk" checked><span class="ufv-toggle-slider"></span></label><button class="ufv-details-close" id="ufv-details-close">&#10005;</button></div></div><div class="ufv-details-body" id="ufv-details-body"></div></div>
                </div>
            </div>
        </div>`;
        document.body.appendChild(overlayEl);
        modalEl = overlayEl.querySelector('.ufv-modal');

        // Intercept wheel events at the window capture phase — before any capture-phase
        // listener the host page may have registered on document (e.g. UniProt's custom
        // scroll container).  We stop propagation for non-viewer overlay events so the
        // host page never sees them, then manually forward scroll to the nearest eligible
        // inner element.  For the 3-D viewer we leave propagation intact so 3Dmol's
        // bubble-phase zoom handler can still run.
        window.addEventListener('wheel', function(e) {
            if (!overlayEl || overlayEl.style.display === 'none') return;
            const viewerEl = document.getElementById('ufv-mol-viewer');
            if (viewerEl && viewerEl.contains(e.target)) return; // let 3Dmol handle zoom
            e.stopPropagation();
            e.preventDefault();
            // Manually scroll the nearest scrollable ancestor within the overlay.
            let el = e.target;
            while (el && el !== overlayEl) {
                const cs = window.getComputedStyle(el);
                if ((cs.overflowY === 'auto' || cs.overflowY === 'scroll') &&
                        el.scrollHeight > el.clientHeight) {
                    const down = e.deltaY > 0;
                    if ((down  && el.scrollTop < el.scrollHeight - el.clientHeight - 1) ||
                        (!down && el.scrollTop > 0)) {
                        el.scrollTop += e.deltaY;
                    }
                    return;
                }
                el = el.parentElement;
            }
        }, { capture: true, passive: false });

        // Reactively lock UniProt's custom scroll container while the modal is open. This catches every
        // scroll cause (focus changes, click handlers, wheel events that UniProt's capture handler saw) and
        // resets scrollTop immediately so the background page never jumps.
        (function () {
            const sc = findHostScrollContainer();
            if (!sc) return;
            let _locking = false;
            sc.addEventListener('scroll', function () {
                if (!_hostScrollLock || _locking) return;
                _locking = true;
                sc.scrollTop = _hostScrollSaved;
                _locking = false;
            });
        })();

        bindEvents();
    }

    function bindEvents() {
        byId('ufv-close').addEventListener('click', close);
        // Single overlay-level handler: stops all modal clicks from reaching UniProt's
        // page handlers (which would scroll Gsgt9), and closes dropdowns on outside-click.
        // Only treat it as a backdrop click-to-close when the press STARTED on the backdrop too —
        // otherwise dragging to rotate the structure and releasing outside the viewer would close.
        let _downOnBackdrop = false;
        overlayEl.addEventListener('mousedown', e => { _downOnBackdrop = (e.target === overlayEl); });
        overlayEl.addEventListener('click', e => {
            if (!e.target.closest('#ufv-dl-wrap')) byId('ufv-dl-menu')?.classList.remove('open');
            if (!e.target.closest('#ufv-cs')) byId('ufv-cs')?.classList.remove('open');
            if (!e.target.closest('#ufv-cm')) byId('ufv-cm')?.classList.remove('open');
            // Partners / Other-chains dropdowns: collapse on any click outside an open dropdown (the modal
            // stops propagation here, so the dropdowns' own document listener never sees these clicks).
            if (!e.target.closest('.ufv-dropdown')) overlayEl.querySelectorAll('.ufv-dropdown.open').forEach(d => d.classList.remove('open'));
            if (e.target === overlayEl && _downOnBackdrop) close();
            _downOnBackdrop = false;
            e.stopPropagation();
        });
        document.addEventListener('keydown', e => { if (e.key === 'Escape') close(); });
        byId('ufv-btn-theme').addEventListener('click', cycleTheme);
        byId('ufv-btn-reset').addEventListener('click', () => {
            const s = UFVState.state;
            // Clear the FULL selection (incl. selectedLigand/_inFocusMode) — otherwise applyMode() below
            // re-enters ligand focus with rezoom:false, leaving focus sticks while zoomed out.
            s.selectedResidue = null;
            s.selectedChain = null;
            s.selectedLigand = null;
            s.nearbyResidues = new Set();
            StructureViewer._selectedResi = null;
            StructureViewer._inFocusMode = false;
            _showOtherSpheres = true; byId('ufv-sphere-chk') && (byId('ufv-sphere-chk').checked = true); // back to default (spheres shown)
            byId('ufv-details').classList.remove('show');
            const defaultMode = s.settings.coloringMode || 'default';
            setColorMode(defaultMode);
            // Restore the window's DEFAULT layer/filter state, not just the camera: clear per-variant
            // show/hide overrides and re-run the window setup so spheres the user toggled off (or
            // diseases they narrowed) come back. Without this, reset only un-focused and the window's
            // default spheres stayed hidden.
            _hiddenVariantKeys.clear();
            _forcedVariantKeys.clear();
            setupWindowDefaults(s.currentMode);
            if (s.currentMode === 'variant') buildVariantFilters();
            StructureViewer.resetView(); // unfocus (removes focus rep) + zoom back to the overview
            // Defer the heavy cartoon rebuild so the zoom animation starts immediately
            requestAnimationFrame(() => {
                applyMode();
                renderSequence();
            });
        });
        byId('ufv-btn-protnlm').addEventListener('click', () => {
            const banner = byId('ufv-protnlm-banner');
            const showing = !banner.classList.contains('ufv-hidden');
            if (showing) { banner.classList.add('ufv-hidden'); return; }
            renderProtNLM();
            banner.classList.remove('ufv-hidden');
        });
        byId('ufv-btn-screenshot').addEventListener('click', () => StructureViewer.screenshot());
        byId('ufv-btn-copy').addEventListener('click', copySelection);
        byId('ufv-btn-export-pdb').addEventListener('click', e => {
            e.stopPropagation();
            byId('ufv-dl-menu').classList.toggle('open');
        });
        byId('ufv-dl-pdb').addEventListener('click', () => {
            byId('ufv-dl-menu').classList.remove('open');
            exportPdb();
        });
        byId('ufv-dl-csv').addEventListener('click', () => {
            byId('ufv-dl-menu').classList.remove('open');
            exportCsv();
        });
        byId('ufv-dl-csv-pv').addEventListener('click', () => {
            byId('ufv-dl-menu').classList.remove('open');
            exportCsv(true);
        });
        byId('ufv-dl-pymol').addEventListener('click', () => {
            byId('ufv-dl-menu').classList.remove('open');
            exportSession('pymol');
        });
        byId('ufv-dl-vmd').addEventListener('click', () => {
            byId('ufv-dl-menu').classList.remove('open');
            exportSession('vmd');
        });
        byId('ufv-ptm-all').addEventListener('click', () => ptmSetAll(true));
        byId('ufv-ptm-none').addEventListener('click', () => ptmSetAll(false));
        byId('ufv-cons-all').addEventListener('click', () => varSectionSetAll('consequence', true));
        byId('ufv-cons-none').addEventListener('click', () => varSectionSetAll('consequence', false));
        byId('ufv-brush-ptm')?.addEventListener('click',     e => { e.stopPropagation(); setColorProfile('ptm'); });
        byId('ufv-brush-ptmsite')?.addEventListener('click', e => { e.stopPropagation(); setColorProfile('sites'); }); // PTM-window Sites C → colour by site
        byId('ufv-brush-vptm')?.addEventListener('click',    e => { e.stopPropagation(); setColorProfile('ptm'); });
        byId('ufv-brush-var')?.addEventListener('click',     e => { e.stopPropagation(); setColorProfile('disease'); });
        byId('ufv-brush-dis')?.addEventListener('click',     e => { e.stopPropagation(); setColorProfile('consequence'); });
        byId('ufv-brush-prov')?.addEventListener('click',    e => { e.stopPropagation(); setColorProfile('provenance'); });
        byId('ufv-brush-vsite')?.addEventListener('click',   e => { e.stopPropagation(); setColorProfile('sites'); });
        byId('ufv-prov-all').addEventListener('click', () => varSectionSetAll('provenance', true));
        byId('ufv-prov-none').addEventListener('click', () => varSectionSetAll('provenance', false));
        byId('ufv-dis-all').addEventListener('click', () => varSectionSetAll('disease', true));
        byId('ufv-dis-none').addEventListener('click', () => varSectionSetAll('disease', false));
        byId('ufv-dis-toggle').addEventListener('click', e => { if (!e.target.closest('button')) toggleCollapsible('ufv-dis-body', 'ufv-dis-toggle'); });
        byId('ufv-prov-toggle').addEventListener('click', e => { if (!e.target.closest('button')) toggleCollapsible('ufv-prov-body', 'ufv-prov-toggle'); });
        byId('ufv-cons-toggle').addEventListener('click', e => { if (!e.target.closest('button')) toggleCollapsible('ufv-cons-body', 'ufv-cons-toggle'); });
        byId('ufv-vptm-toggle').addEventListener('click', e => { if (!e.target.closest('button')) toggleCollapsible('ufv-vptm-body', 'ufv-vptm-toggle'); });
        byId('ufv-vptm-all').addEventListener('click', () => variantPtmSetAll(true));
        byId('ufv-vptm-none').addEventListener('click', () => variantPtmSetAll(false));
        byId('ufv-sites-ptm-toggle').addEventListener('click', e => { if (!e.target.closest('button')) toggleCollapsible('ufv-sites-ptm-body', 'ufv-sites-ptm-toggle'); });
        byId('ufv-sites-ptm-all').addEventListener('click', () => sitesSetAll(true));
        byId('ufv-sites-ptm-none').addEventListener('click', () => sitesSetAll(false));
        byId('ufv-sites-var-toggle').addEventListener('click', e => { if (!e.target.closest('button')) toggleCollapsible('ufv-sites-var-body', 'ufv-sites-var-toggle'); });
        byId('ufv-sites-var-all').addEventListener('click', () => sitesSetAll(true));
        byId('ufv-sites-var-none').addEventListener('click', () => sitesSetAll(false));
        const ligHdrClick = (body, tog) => e => { if (!e.target.closest('button') && !e.target.closest('.ufv-toggle-switch')) toggleCollapsible(body, tog); };
        byId('ufv-ligands-ptm-toggle').addEventListener('click', ligHdrClick('ufv-ligands-ptm-body', 'ufv-ligands-ptm-toggle'));
        byId('ufv-ligands-var-toggle').addEventListener('click', ligHdrClick('ufv-ligands-var-body', 'ufv-ligands-var-toggle'));
        byId('ufv-ligands-ptm-all').addEventListener('click', () => ligandsSetAll(true));
        byId('ufv-ligands-ptm-none').addEventListener('click', () => ligandsSetAll(false));
        byId('ufv-ligands-var-all').addEventListener('click', () => ligandsSetAll(true));
        byId('ufv-ligands-var-none').addEventListener('click', () => ligandsSetAll(false));
        // Exclude water & ions — synced between the two panels' toggles.
        // Ions stay in the list (as unchecked rows) so the user can see what was filtered; only the
        // 3D visibility changes via hiddenLigands, matching the per-ligand checkbox behaviour.
        const onIonsToggle = checked => {
            StructureViewer.excludeIons = checked;
            byId('ufv-ligands-ions-ptm').checked = checked;
            byId('ufv-ligands-ions-var').checked = checked;
            const ions = StructureViewer.ION_CODES || new Set();
            UFVState.state.ligands.filter(l => ions.has(l.resn)).forEach(l => {
                checked ? StructureViewer.hiddenLigands.add(ligKey(l)) : StructureViewer.hiddenLigands.delete(ligKey(l));
            });
            StructureViewer._drawLigands();
            buildLigandFilters();
        };
        byId('ufv-ligands-ions-ptm').addEventListener('change', e => onIonsToggle(e.target.checked));
        byId('ufv-ligands-ions-var').addEventListener('change', e => onIonsToggle(e.target.checked));
        byId('ufv-details-close').addEventListener('click', () => {
            const s = UFVState.state;
            byId('ufv-details').classList.remove('show');
            // Closing the panel after a zoom-in returns the structure to its overview.
            const wasFocused = s.selectedResidue != null || s.selectedLigand != null || StructureViewer._inFocusMode;
            s.selectedResidue = null;
            s.selectedChain = null;
            s.selectedLigand = null;
            s.nearbyResidues = new Set();
            StructureViewer._selectedResi = null;
            StructureViewer._inFocusMode = false;
            _showOtherSpheres = true; byId('ufv-sphere-chk') && (byId('ufv-sphere-chk').checked = true); // back to default (spheres shown)
            if (wasFocused && StructureViewer.viewer) {
                StructureViewer.resetView();
                requestAnimationFrame(() => applyMode());
            }
        });
        // Header sphere-visibility toggle: controls whether other annotation spheres stay visible
        // while zoomed into a residue.  Always available (PTM / variant / disease views).
        byId('ufv-sphere-chk').addEventListener('change', e => {
            _showOtherSpheres = e.target.checked;
            // Update ONLY the annotation spheres — NOT a re-focus. Re-focusing (even rezoom:false) rebuilt
            // the ufv-focus sticks and re-ran setCartoon, causing the green flash + camera twitch.
            if (StructureViewer.currentStructure) StructureViewer.setOtherSpheresVisible(_showOtherSpheres);
        });
        byId('ufv-cs-btn').addEventListener('click', e => { e.stopPropagation(); byId('ufv-cs').classList.toggle('open'); });
        byId('ufv-cm-btn').addEventListener('click', e => { e.stopPropagation(); byId('ufv-cm').classList.toggle('open'); });
        byId('ufv-cm-drop').querySelectorAll('.ufv-cm-opt').forEach(opt => {
            opt.addEventListener('click', () => {
                byId('ufv-cm').classList.remove('open');
                const val = opt.dataset.value;
                setColorMode(val);
                deselectDomainsOnColorSwitch(); // clear domain ticks before turning domain-colouring off
                _domainCartoon = false; // explicit colour-mode pick overrides the Domains-C backbone colouring
                syncColorProfileButtons();
                // In-modal coloring is session-only — the persistent startup default lives in the
                // options page (defaults to "Default"), so the viewer always opens on cyan
                // and the expensive constraint-pocket compute is never triggered automatically.
                const isPocket = val === 'prism';
                byId('ufv-sens-wrap').classList.toggle('ufv-hidden', !isPocket);
                if (isPocket) {
                    ensurePocketAnalysis().then(() => applyMode());
                } else {
                    requestAnimationFrame(() => applyMode());
                }
            });
        });
        // Sensitivity slider: update the q label live, but only re-colour on release ('change')
        // — re-colouring on every drag tick caused visible lag.  Re-thresholding the cached
        // candidates is cheap; the cost is the 3Dmol cartoon recolour, so we do it once on release.
        byId('ufv-sens-slider').addEventListener('input', e => {
            sensThreshold = Math.max(0.01, Math.min(0.40, e.target.value / 100));
            byId('ufv-sens-q').textContent = sensThreshold.toFixed(2);
        });
        byId('ufv-sens-slider').addEventListener('change', () => {
            requestAnimationFrame(() => applyMode());
        });
        byId('ufv-structure-prev').addEventListener('click', () => cycleStructure(-1));
        byId('ufv-structure-next').addEventListener('click', () => cycleStructure(1));
    }

    function bindSettings() {
        // Settings moved to extension options page (right-click extension icon → Options)
    }

    async function open(mode, opts = {}) {
        const s = UFVState.state;
        const mySeq = ++_openSeq;
        s.currentMode = mode;
        // Each window's focus drives its default sphere colouring: PTM window → PTM-type colours,
        // Disease & Variants → disease colours, Functional features → site colours, else consequence.
        _colorProfile = { ptm: 'ptm', variant: 'disease', sites: 'sites' }[mode] || 'consequence';
        build();
        syncColorProfileButtons();
        await UFVState.loadSettings();
        if (_openSeq !== mySeq) return; // superseded by a newer open()
        syncSettingsControls();
        byId('ufv-id-badge').textContent = s.uniprotId;
        const HEADINGS = { ptm: 'PTM Viewer', variant: 'Disease & Variants', sites: 'Functional Features', domains: 'Family & Domains', structure: 'Structure', subcellular: 'Subcellular Location' };
        byId('ufv-modal-heading').textContent = HEADINGS[mode] || '3D Feature Viewer';
        // Structure & Subcellular reuse the feature panel — they show the full layer set, all off.
        const featPanelModes = mode === 'sites' || mode === 'structure' || mode === 'subcellular';
        byId('ufv-ptm-panel').classList.toggle('ufv-hidden', mode !== 'ptm');
        byId('ufv-var-panel').classList.toggle('ufv-hidden', mode !== 'variant');
        byId('ufv-feat-panel').classList.toggle('ufv-hidden', !featPanelModes);
        byId('ufv-dom-panel').classList.toggle('ufv-hidden', mode !== 'domains');
        // The Family & Domains window colours the cartoon by domain, so the colour-mode picker
        // isn't meaningful there — hide it. Shown for every other window.
        byId('ufv-cm')?.classList.toggle('ufv-hidden', mode === 'domains');
        _hostScrollLock = false; // reset so focus-scroll and UniProt handlers can settle
        overlayEl.style.display = 'flex';
        document.documentElement.style.overflow = 'hidden';
        document.body.style.overflow = 'hidden';
        // One animation frame is enough for UniProt's focus-scroll and any other deferred scroll handlers
        // to fire and the scroll container to reach its natural resting position before we lock it in place.
        requestAnimationFrame(() => {
            const sc = findHostScrollContainer();
            if (sc) _hostScrollSaved = sc.scrollTop;
            _hostScrollLock = true;
        });
        // Reset any residue/ligand focus + detail state from a previous session so opening a
        // different view (e.g. PTM → Variants) starts from defaults, not a continuation.
        s.selectedResidue = null;
        s.selectedChain = null;
        s.selectedLigand = null;
        s.nearbyResidues = new Set();
        StructureViewer._selectedResi = null;
        StructureViewer._inFocusMode = false;
        StructureViewer.showLigands = true;
        StructureViewer.excludeIons = false;
        byId('ufv-details').classList.remove('show');
        // Immediately reapply the new mode coloring so the viewer doesn't flash stale state
        if (s.loaded && s.annotationsLoaded && StructureViewer.viewer) {
            applyMode();
        }
        // Join any in-progress background prefetch so we don't double-fetch
        if (s.loadingPromise) await s.loadingPromise;
        if (_openSeq !== mySeq) return;
        if (!s.annotationsLoaded) {
            await loadAnnotations();
            if (_openSeq !== mySeq) return;
        }
        buildFilters();
        setupWindowDefaults(s.currentMode); // primary-on / secondary-off defaults + dynamic panel
        if (!s.loaded) {
            loadStructuresAndShow(); // intentionally not awaited — shows viewer async
        } else {
            // Reopening a window starts from a reset view: snap back to the default structure
            // (loadSelectedStructure resets camera/coloring/focus too).
            chooseDefaultStructure();
            renderStructureSelector();
            await loadSelectedStructure();
        }
    }

    async function loadAnnotations() {
        const s = UFVState.state;
        const requestedId = s.uniprotId;
        showLoading('Loading annotations…');
        try {
            const data = await UFVApi.loadFeatureData(requestedId);
            if (UFVState.state.uniprotId !== requestedId) return;
            Object.assign(s, data);
            s.ptmGroups = DataProcessor.groupPTMsByCategory(s.ptms);
            s.activeConsequences = new Set(Object.keys(DataProcessor.getConsequenceSummary(s.variants)));
            s.activeProvenances = new Set(Object.keys(DataProcessor.getProvenanceSummary(s.variants)));
            s.analysis.alphaMissense = UFVAnalysis.aggregateAlphaMissense(s.variants, s.amMap);
            DataProcessor.computeDiseaseColors(s.variants);
            const _ds = DataProcessor.getDiseaseSummary(s.variants);
            _diseaseColorMap = new Map(Object.entries(_ds).map(([n, m]) => [n, m.color]));
            s.annotationsLoaded = true;
            if (!byId('ufv-protnlm-banner')?.classList.contains('ufv-hidden')) renderProtNLM();
            const fc = s.functionContext;
            if (fc) {
                const bits = [];
                if (fc.summary) bits.push('Function: ' + fc.summary);
                if (fc.catalytic?.length) bits.push('Catalytic: ' + fc.catalytic.slice(0, 3).join('; '));
                if (fc.locations?.length) bits.push('Subcellular location: ' + fc.locations.join(', '));
                const h = byId('ufv-modal-heading'); if (h) h.title = bits.join('\n\n');
            }
        } catch (err) {
            showError(err.message || 'Unable to load annotations.');
        }
    }

    async function loadStructuresAndShow() {
        const s = UFVState.state;
        const requestedId = s.uniprotId;
        showLoading('Finding structures…');
        let shown = false;
        // Fast path: paint the canonical AlphaFold model immediately, before the (slower)
        // experimental/isoform/computed discovery finishes. Skip it when the user prefers a
        // different default (experimental / best-coverage) — that needs the full list first.
        const prefersAlphaFold = !s.settings.defaultStructure || s.settings.defaultStructure === 'alphafold';
        try {
            const primary = prefersAlphaFold ? await UFVApi.getPrimaryStructure(requestedId, s.sequence.length) : null;
            if (UFVState.state.uniprotId !== requestedId) return;
            if (primary) {
                s.structures = [primary];
                s.selectedStructureIndex = 0;
                s.loaded = true;
                renderStructureSelector();
                await loadSelectedStructure();
                shown = true;
            }
        } catch (_) { /* fall through to full discovery */ }
        // Full discovery streams in behind the AlphaFold model (or becomes the first display when
        // the protein has no AlphaFold model).
        try {
            const all = await UFVApi.getStructures(requestedId, s.sequence.length);
            if (UFVState.state.uniprotId !== requestedId) return;
            mergeStreamedStructures(all || [], shown);
        } catch (err) {
            if (!shown) showError(err.message || 'Unable to load structures.');
        }
    }

    /**
     * Fold the full structure list in once discovery completes. If the AlphaFold model is already
     * on screen we keep it shown and just grow the selector (no reload, no flash). If nothing was
     * shown yet (no AlphaFold model) we pick the default and load it now.
     */
    function mergeStreamedStructures(all, alreadyShown) {
        const s = UFVState.state;
        if (!all.length) {
            if (!alreadyShown) showError('No 3D structures found for this protein.');
            return;
        }
        // Preserve whatever the user is currently viewing by URL, so its index stays correct even
        // though the list grew (and even if the user switched structures meanwhile).
        const currentUrl = alreadyShown ? UFVState.selectedStructure()?.url : null;
        s.structures = all;
        s.loaded = true;
        if (alreadyShown) {
            const keep = all.findIndex(st => st.url === currentUrl);
            s.selectedStructureIndex = keep >= 0 ? keep : 0;
            renderStructureSelector(); // model already painted — just refresh the list
        } else {
            chooseDefaultStructure();
            renderStructureSelector();
            loadSelectedStructure();
        }
    }

    function chooseDefaultStructure() {
        const s = UFVState.state;
        if (!s.structures.length) return;
        const pref = s.settings.defaultStructure;
        if (pref === 'experimental') {
            const i = s.structures.findIndex(x => x.source === 'PDB');
            s.selectedStructureIndex = i >= 0 ? i : 0;
        } else if (pref === 'coverage') {
            let best = 0;
            s.structures.forEach((x, i) => { if ((x.coverage || 0) > (s.structures[best].coverage || 0)) best = i; });
            s.selectedStructureIndex = best;
        } else {
            s.selectedStructureIndex = 0;
        }
    }

    async function loadSelectedStructure() {
        const s = UFVState.state;
        const requestedId = s.uniprotId;
        const mySeq = ++_loadSeq; // re-entrancy guard: a newer load supersedes this one
        const structure = UFVState.selectedStructure();
        if (!structure) {
            showError('No AlphaFold or mapped PDB structure is available for this protein.');
            return;
        }
        resetViewerTransients(); // clean slate on every structure (re)load / launch
        const alreadyLoaded = StructureViewer.currentStructure?.url === structure.url;
        if (!alreadyLoaded) {
            showLoading(`Loading ${structure.label}...`);
            try {
                if (!StructureViewer.viewer) StructureViewer.init(byId('ufv-mol-viewer'));
                try {
                    await StructureViewer.loadStructure(structure.url, structure);
                } catch (pdbErr) {
                    // PDB format unavailable (e.g. large structures only published as mmCIF)
                    if (structure.cifUrl && !structure.url.toLowerCase().endsWith('.cif')) {
                        console.warn('[UniProt 3D] PDB format failed, retrying with mmCIF:', pdbErr.message);
                        structure.url = structure.cifUrl;
                        await StructureViewer.loadStructure(structure.url, structure);
                    } else {
                        throw pdbErr;
                    }
                }
                // Bail if the protein changed or a newer load started while we were fetching —
                // prevents two interleaved loads from rendering against a half-built model.
                if (UFVState.state.uniprotId !== requestedId || _loadSeq !== mySeq) return;
                // Structure-dependent analyses are (re)computed for THIS structure so the
                // per-chain hotspot / contact-hub results reflect the displayed subunits.
                // Sequence-based residueBurden is identical regardless of structure.
                // (Partner-protein disease residues are folded in afterwards, off the critical
                // path, by augmentHotspotsWithPartners so loading a complex isn't delayed.)
                // The graph analyses (3-D hotspot betweenness, long-range contact hubs) are
                // heavy. Clear any stale results from the previous structure and compute them
                // AFTER this structure paints (scheduleStructureAnalyses), so a large structure
                // can't freeze the UniProt tab on load. The cheap sequence-based burden stays here.
                s.analysis.hotspots = null;
                s.analysis.hotspotsByChain = null;
                s.analysis.hotspotMethod = null;
                s.analysis.distantContacts = null;
                s.analysis.distantContactsByChain = null;
                s.analysis.residueBurden = UFVAnalysis.computeResidueBurden(s.variants);
                // Constraint-pocket analysis is structure-dependent (geometry + PAE) — drop so it
                // recomputes for the newly loaded structure on next selection.
                s.analysis.prism = null;
                byId('ufv-loading').classList.add('hidden');
            } catch (err) {
                if (_loadSeq === mySeq) showError(err.message || 'Unable to load selected structure.');
                return;
            }
        }
        if (_loadSeq !== mySeq) return; // a newer load owns the viewer now
        StructureViewer.hoverCb = onHover;
        StructureViewer.clickCb = onClick;
        StructureViewer.ligandClickCb = onLigandClick;
        // Enumerate ligands present in the loaded model (AlphaFill/SwissModel cofactors etc.).
        s.ligands = StructureViewer.enumerateLigands ? StructureViewer.enumerateLigands() : [];
        // Family & Domains window: a domain the loaded structure doesn't actually resolve can't be
        // drawn, so start it unselected (recomputed whenever the selected structure changes).
        if (s.currentMode === 'domains') applyDomainCoverageDefaults();
        refreshLigandSections();
        annotateAlphaFillLigands(mySeq); // async: adds per-ligand identity + an identity-threshold filter
        // Double-click closes the detail panel and resets the view — same as the × button.
        StructureViewer.dblClickCb = () => { byId('ufv-details-close')?.click(); };
        // The observed-residue cache is built from an async 'atoms' event that arrives AFTER this window
        // is first drawn, so unresolved residues (PTMs/sites not modelled in this structure) aren't greyed
        // yet. Re-grey the current filter window once the cache lands.
        StructureViewer.observedResiCb = () => { if (_loadSeq === mySeq && UFVState.state.currentMode) rebuildCurrentWindow(); };
        // Already-loaded structures skip the loadStructure() that would re-frame the camera, so
        // reset the view here — reopening a window should return to the default framing, not keep
        // wherever the previous session was zoomed/panned.
        if (alreadyLoaded) { StructureViewer.viewer?.zoomTo(); StructureViewer.viewer?.zoom(1.15); }
        updateStructureMeta();
        updateTractabilityNav(); // whole-protein tractability dropdown in the header (cached per protein)
        applyMode();
        // Rebuild the dynamic filter panels so the unresolved-residue greying reflects THIS structure.
        // The builders read visibility from the data model, so toggle state is preserved across the rebuild.
        const cm = s.currentMode;
        if (cm === 'structure' || cm === 'subcellular') buildStructureWindow(cm);
        else if (cm === 'sites' || cm === 'domains') buildFeatureWindow(cm);
        // If the constraint-pocket mode is active, (re)compute it for this structure then recolour.
        if (getColorMode() === 'prism') ensurePocketAnalysis().then(() => applyMode());
        // Heavy graph analyses (hotspots, contact hubs) run off the critical path so the
        // structure shows immediately; they recolour / augment when ready.
        scheduleStructureAnalyses(structure, requestedId, mySeq);
        // Multichain: fetch the OTHER subunits' (partner proteins') annotations and overlay their disease
        // variants, off the critical path. Clicking a partner residue then shows that protein's data.
        loadPartnerAnnotations(structure, requestedId, mySeq);
    }

    // Partner-protein support (other UniProt entries present as separate chains in this structure).
    // No persistent overlay. We fetch the partners' disease variants + PTMs + sites to (a) colour partner
    // neighbours with OUR cutoff/colours when one of OUR residues is focused, (b) list them in Nearby, and
    // (c) drive the opt-in "Partners" accordion (checkboxes for disease variants / PTMs / sites → hoverable
    // spheres). The Partners checkbox does NOT tint chains. Clicking a partner residue gives a header-only
    // minimal panel + a link to that protein's UniProt entry, and zooms it with OUR disease stick colours.
    const _PARTNER_SEV = ['Likely pathogenic or pathogenic', 'Predicted deleterious', 'Uncertain significance', 'Likely benign or benign'];
    let _partnerAnnotations = [];   // [{ chainId, accession, byResi: Map(pdbResi -> {uniPos, variants, ptms, sites}) }]
    let _partnerChainAcc = new Map(); // chain -> accession, for ALL partner chains (even unannotated)
    let _showPartners = false;      // Partners activated? (off by default)
    let _partnerLayers = { variants: false, ptms: false, sites: false }; // which partner sphere layers are on

    async function loadPartnerAnnotations(structure, requestedId, mySeq) {
        _partnerAnnotations = [];
        _partnerChainAcc = new Map();
        StructureViewer._partnerColorMap = new Map();
        (structure?.partners || []).forEach(p => { if (p.chainId && p.accession && !_partnerChainAcc.has(p.chainId)) _partnerChainAcc.set(p.chainId, p.accession); });
        StructureViewer.partnerClickCb = onPartnerClick;
        StructureViewer.setPartnerSpheres?.([]);
        updatePartnerToggle();
        updateProteinNav();
        if (!structure?.partners?.length || !UFVApi.loadPartnerAnnotations) return;
        let data = [];
        try { data = await UFVApi.loadPartnerAnnotations(structure); } catch (_) { data = []; }
        if (UFVState.state.uniprotId !== requestedId || _loadSeq !== mySeq) return;
        _partnerAnnotations = data || [];
        // Partner residue -> our disease colour, so focusing one of OUR interface residues (or a partner
        // residue) colours the partner-chain contacts with OUR consequence colouring.
        const cmap = new Map();
        _partnerAnnotations.forEach(p => p.byResi.forEach((info, pdbResi) => {
            const top = _PARTNER_SEV.map(sev => (info.variants || []).find(v => v.consequence === sev)).find(Boolean) || (info.variants || [])[0];
            if (top && top.consequenceColor) cmap.set(p.chainId + '|' + pdbResi, top.consequenceColor);
        }));
        StructureViewer._partnerColorMap = cmap;
        renderPartnerSpheres();
        if (UFVState.state.selectedResidue != null) applyMode(); // re-colour any active focus's partner contacts
    }

    // Render partner annotation spheres for the enabled layers (disease variants / PTMs / sites). One sphere
    // per partner residue, coloured by the highest-priority enabled annotation, with a hover label.
    function renderPartnerSpheres() {
        if (!_showPartners) { StructureViewer.setPartnerSpheres?.([]); return; }
        const spheres = [];
        _partnerAnnotations.forEach(p => p.byResi.forEach((info, pdbResi) => {
            let color, body;
            if (_partnerLayers.variants) {
                const path = (info.variants || []).filter(v => /pathogenic|deleterious/i.test(v.consequence || ''));
                if (path.length) {
                    // Prefer a variant with a named disease so the sphere is coloured BY DISEASE (not by
                    // pathogenicity consequence); fall back to the most-severe variant's colour.
                    const withDisease = path.find(v => v.diseaseColor && v.diseaseColor !== '#9e9e9e');
                    const top = withDisease || _PARTNER_SEV.map(sev => path.find(v => v.consequence === sev)).find(Boolean) || path[0];
                    color = top.diseaseColor || top.consequenceColor || '#e53935';
                    body = `${path.length} pathogenic variant${path.length === 1 ? '' : 's'}`;
                }
            }
            if (color == null && _partnerLayers.ptms && info.ptms && info.ptms.length) { color = info.ptms[0].color || '#ab47bc'; body = _esc(info.ptms[0].category || info.ptms[0].type || 'PTM'); }
            if (color == null && _partnerLayers.sites && info.sites && info.sites.length) { color = DataProcessor.SITE_COLOR; body = _esc(info.sites[0].description || 'Functional site'); }
            if (color == null) return;
            const label = `<div class="ufv-tip-hdr">${_esc(p.accession)} · ${info.uniPos}</div><div class="ufv-tip-body">${body} (partner protein)</div>`;
            spheres.push({ chain: p.chainId, resi: pdbResi, color, radius: 1.6, label });
        }));
        StructureViewer.setPartnerSpheres?.(spheres);
    }

    // Shared dropdown wiring: toggle .open on the button, close on any outside click (one global listener),
    // and close other open dropdowns. Menu clicks don't bubble so toggling checkboxes keeps it open.
    let _dropdownCloseWired = false;
    function wireDropdown(host, btn) {
        if (!_dropdownCloseWired) { _dropdownCloseWired = true; document.addEventListener('click', () => document.querySelectorAll('.ufv-dropdown.open').forEach(d => d.classList.remove('open'))); }
        btn.addEventListener('click', (e) => {
            e.stopPropagation();
            const open = host.classList.contains('open');
            document.querySelectorAll('.ufv-dropdown.open').forEach(d => { if (d !== host) d.classList.remove('open'); });
            host.classList.toggle('open', !open);
        });
    }

    // Dropdown next to "Copy": the OTHER CHAINS (different proteins) present in this structure. Each entry
    // is a partner chain → opens that protein's UniProt entry. Hidden when there are no partner chains.
    function updateProteinNav() {
        const host = byId('ufv-protein-nav'); if (!host) return;
        const chains = [..._partnerChainAcc.entries()]; // [chain, accession]
        host.textContent = '';
        host.className = 'ufv-protein-nav ufv-dropdown';
        if (!chains.length) return;
        const btn = document.createElement('button');
        btn.type = 'button'; btn.className = 'ufv-dropdown-btn';
        btn.innerHTML = `Other chains (${chains.length}) <span class="ufv-dropdown-caret">▾</span>`;
        const menu = document.createElement('div');
        menu.className = 'ufv-dropdown-menu ufv-dropdown-menu-up';
        chains.forEach(([chain, acc]) => {
            const a = document.createElement('a');
            a.href = `https://www.uniprot.org/uniprotkb/${acc}/entry`; a.target = '_blank'; a.rel = 'noopener noreferrer';
            a.className = 'ufv-dropdown-item';
            a.innerHTML = `<span class="ufv-nav-chain">${_esc(chain)}</span>${_esc(acc)}`;
            menu.appendChild(a);
        });
        wireDropdown(host, btn);
        host.append(btn, menu);
    }

    // "Partners" dropdown in the title row: a button that opens a menu of layer checkboxes (disease variants
    // / PTMs / sites). Checking a layer draws those partner annotations as hoverable spheres (no tint).
    function updatePartnerToggle() {
        const host = byId('ufv-modal-heading')?.parentElement; if (!host) return;
        let wrap = byId('ufv-partner-toggle-wrap');
        if (!_partnerChainAcc.size) { wrap?.remove(); return; }
        if (!wrap) {
            wrap = document.createElement('span');
            wrap.id = 'ufv-partner-toggle-wrap'; wrap.className = 'ufv-dropdown ufv-partner-dropdown';
            const cats = [['variants', 'Disease variants'], ['ptms', 'PTMs'], ['sites', 'Sites']];
            wrap.innerHTML = '<button type="button" class="ufv-dropdown-btn" title="Show annotations of the other proteins in this complex">Partners <span class="ufv-dropdown-caret">▾</span></button>'
                + '<div class="ufv-dropdown-menu">' + cats.map(([k, lbl]) =>
                    `<label class="ufv-dropdown-check"><input type="checkbox" data-layer="${k}"><span>${lbl}</span></label>`).join('') + '</div>';
            host.appendChild(wrap);
            const btn = wrap.querySelector('.ufv-dropdown-btn');
            const menu = wrap.querySelector('.ufv-dropdown-menu');
            wireDropdown(wrap, btn);
            menu.addEventListener('click', e => e.stopPropagation()); // toggling a layer keeps the menu open
            menu.querySelectorAll('input').forEach(inp => inp.addEventListener('change', e => {
                _partnerLayers[e.target.dataset.layer] = e.target.checked;
                _showPartners = _partnerLayers.variants || _partnerLayers.ptms || _partnerLayers.sites;
                wrap.classList.toggle('ufv-active', _showPartners);
                renderPartnerSpheres();
            }));
        }
        wrap.querySelectorAll('input[data-layer]').forEach(i => { i.checked = !!_partnerLayers[i.dataset.layer]; });
        wrap.classList.toggle('ufv-active', _showPartners);
    }

    // Clicking a partner-subunit residue: normal Mol* behaviour (zoom + surroundings, our disease stick
    // colours) and a HEADER-ONLY minimal panel — residue + chain, plus a button to that protein's UniProt
    // entry. No body.
    function onPartnerClick({ chain, resi, resn }) {
        StructureViewer.focusPartnerResidue?.(chain, resi);
        const acc = _partnerChainAcc.get(chain);
        const body = byId('ufv-details-body');
        const titleEl = byId('ufv-details-title');
        const name = (resn || '').trim();
        titleEl.innerHTML = `<span class="ufv-partner-res">${name ? _esc(name) + ' ' : ''}${resi}</span><span class="ufv-title-sub">chain ${_esc(chain)}</span>`;
        if (acc) {
            const a = document.createElement('a');
            a.href = `https://www.uniprot.org/uniprotkb/${acc}/entry`; a.target = '_blank'; a.rel = 'noopener noreferrer';
            a.className = 'ufv-partner-goto'; a.title = `Open ${acc} on UniProt`; a.textContent = `${acc} ↗`;
            titleEl.appendChild(a);
        }
        body.textContent = ''; // header-only: no body for partner residues
        byId('ufv-sphere-toggle')?.classList.add('ufv-hidden'); // sphere show/hide doesn't apply to a partner residue
        byId('ufv-details').classList.add('show');
    }

    // Append partner-protein neighbours to the Nearby list: each shows the partner residue (coloured by its
    // top variant); clicking it focuses + opens that partner residue's panel. partnerList items are
    // { chain, resi, partner, info } as resolved in rebuildList.
    function appendPartnerNeighbours(listEl, partnerList) {
        partnerList.forEach((x, i) => {
            const span = document.createElement('span');
            span.className = 'ufv-nearby-res ufv-nearby-partner';
            span.textContent = `${x.info.variants[0]?.wildType || ''}${x.info.uniPos}`;
            const top = x.info.variants.find(v => /pathogenic|deleterious/i.test(v.consequence || '')) || x.info.variants[0];
            if (top?.consequenceColor) span.style.color = top.consequenceColor;
            span.title = `${x.partner.accession} · chain ${x.chain} (partner protein)`;
            span.addEventListener('click', () => onPartnerClick({ chain: x.chain, resi: x.resi }));
            listEl.appendChild(span);
            if (i < partnerList.length - 1) listEl.appendChild(document.createTextNode(', '));
        });
    }

    /**
     * Computes the structure-dependent graph analyses (3-D variant-enrichment hotspots and
     * long-range contact hubs) after the structure has painted, so a large structure never
     * blocks the main thread during load.  Guarded by the load sequence/accession so a stale
     * deferred run can't overwrite a newer structure's results.  Recolours only if the active
     * mode needs the result, then folds in partner-protein disease residues.
     */
    function scheduleStructureAnalyses(structure, requestedId, mySeq) {
        const run = () => {
            const s = UFVState.state;
            if (s.uniprotId !== requestedId || _loadSeq !== mySeq || !StructureViewer.viewer) return;
            const hotspots = UFVAnalysis.computeHotspots(StructureViewer.viewer, s.variants, structure);
            s.analysis.hotspots = hotspots.merged;
            s.analysis.hotspotsByChain = hotspots.byChain;
            s.analysis.hotspotMethod = hotspots.method;
            const contacts = UFVAnalysis.computeDistantContacts(StructureViewer.viewer, structure, s.variants);
            s.analysis.distantContacts = contacts.merged;
            s.analysis.distantContactsByChain = contacts.byChain;
            if (s.uniprotId !== requestedId || _loadSeq !== mySeq) return;
            const m = getColorMode();
            if (m === 'hotspots' || m === 'distantContacts') applyMode();
            // Fold neighbouring partner-protein disease residues into the hotspot test
            // (network fetch) once the base hotspots exist.
            augmentHotspotsWithPartners(structure, requestedId, mySeq);
        };
        if (typeof requestIdleCallback === 'function') requestIdleCallback(run, { timeout: 1200 });
        else setTimeout(run, 0);
    }

    /**
     * After a hetero-complex structure is shown, fetch the disease residues of its partner
     * proteins and recompute the 3-D hotspots so their proximity is accounted for.  Partner
     * annotations are never displayed; only OUR protein's hotspot tiers may change.
     */
    async function augmentHotspotsWithPartners(structure, requestedId, loadSeq) {
        if (!structure?.partners?.length) return;
        let partnerPoints = [];
        try { partnerPoints = await UFVApi.loadPartnerClassified(structure); } catch (_) { return; }
        if (!partnerPoints.length) return;
        const s = UFVState.state;
        if (s.uniprotId !== requestedId || _loadSeq !== loadSeq || UFVState.selectedStructure() !== structure) return;
        const hotspots = UFVAnalysis.computeHotspots(StructureViewer.viewer, s.variants, structure, 8, partnerPoints);
        s.analysis.hotspots = hotspots.merged;
        s.analysis.hotspotsByChain = hotspots.byChain;
        s.analysis.hotspotMethod = hotspots.method;
        // Re-render only if the user is currently looking at the hotspot view.
        if (getColorMode() === 'hotspots') applyMode();
    }

    function shortMethod(m) {
        if (!m) return '';
        const lc = m.toLowerCase();
        if (lc.includes('electron') || lc.includes('cryo')) return 'EM';
        if (lc.includes('x-ray') || lc.includes('diffraction') || lc.includes('crystallography')) return 'X-ray';
        if (lc.includes('nmr')) return 'NMR';
        if (lc.includes('neutron')) return 'Neutron';
        return m.split(' ').map(w => w[0]?.toUpperCase() || '').join('') || 'Exp';
    }

    function stLabel(st) {
        if (st?.source === 'AlphaFold' && st.isoform) return `AlphaFold ${st.isoform}`;
        if (!st || st.source === 'AlphaFold') return 'AlphaFold';
        if (st.source === 'Computed') {
            return `${st.provider}${st.coverage ? ` (model, ${st.coverage}%)` : ' (model)'}`;
        }
        const chains = st.chainIds?.length > 1 ? st.chainIds.join(',') : (st.chainId || '?');
        const parts = [];
        if (st.method) parts.push(shortMethod(st.method));
        if (st.resolution) parts.push(`${st.resolution}Å`);
        if (st.coverage) parts.push(`${st.coverage}%`);
        const suffix = parts.length ? ` (${parts.join(', ')})` : '';
        // ⚛ flags a cross-species / chimeric construct (proteins from >1 organism) — our mapping covers only
        // part of the chains, so some neighbours fall in unmapped (other-species) regions.
        return `${st.pdbId}-${chains}${suffix}${st.chimeric ? ' ⚛' : ''}`;
    }

    function renderStructureSelector() {
        const drop = byId('ufv-cs-drop');
        drop.textContent = '';
        const idx = UFVState.state.selectedStructureIndex;
        UFVState.state.structures.forEach((st, i) => {
            const opt = document.createElement('div');
            opt.className = 'ufv-cs-opt' + (i === idx ? ' selected' : '');
            opt.textContent = stLabel(st);
            if (st.chimeric) opt.title = chimericTip(st);
            opt.addEventListener('click', async () => {
                byId('ufv-cs').classList.remove('open');
                UFVState.state.selectedStructureIndex = i;
                await loadSelectedStructure();
            });
            drop.appendChild(opt);
        });
        const sel = UFVState.state.structures[idx];
        const btn = byId('ufv-cs-btn');
        if (btn && sel) {
            btn.textContent = stLabel(sel);
            btn.title = chimericTip(sel);
        }
    }

    // Tooltip for the ⚛ chimeric flag — names the organisms when our chain spans >1, else notes the fusion.
    function chimericTip(st) {
        if (!st || !st.chimeric) return '';
        const orgs = st.organisms || [];
        return orgs.length > 1 ? `⚛ Chimeric chain — ${orgs.join(' / ')}` : '⚛ Chimeric chain — our protein is fused with another sequence';
    }

    function getColorMode() {
        return byId('ufv-cm-drop')?.querySelector('.ufv-cm-opt.selected')?.dataset.value || UFVState.state.settings.coloringMode || 'default';
    }

    // The colouring a window should OPEN with (its focus): Subcellular → membrane topology, Structure →
    // pLDDT, everything else → the configured default. Used by both setupWindowDefaults and the per-load
    // transient reset, so re-opening a window restores its intent without needing a manual reset.
    function windowDefaultColorMode() {
        const m = UFVState.state.currentMode;
        if (m === 'subcellular') return 'topology';
        if (m === 'structure') return 'plddt';
        let mode = UFVState.state.settings.coloringMode || 'default';
        if (mode === 'prism' || mode === 'topos') mode = 'default';
        return mode;
    }

    function setColorMode(value) {
        const drop = byId('ufv-cm-drop');
        const btn = byId('ufv-cm-btn');
        // Fall back to 'default' if the requested mode isn't a real option (e.g. a stale
        // saved preference like the long-removed 'ptmVariant') — otherwise nothing would
        // be selected and the UI would silently show no coloring mode.
        if (drop && !drop.querySelector(`.ufv-cm-opt[data-value="${value}"]`)) value = 'default';
        drop?.querySelectorAll('.ufv-cm-opt').forEach(opt => {
            opt.classList.toggle('selected', opt.dataset.value === value);
            if (opt.dataset.value === value && btn) btn.textContent = opt.textContent;
        });
    }

    // Reset transient viewer state to defaults — run on every structure (re)load and on launch
    // so switching structures or reopening the viewer always starts clean.
    function resetViewerTransients() {
        const s = UFVState.state;
        s.selectedResidue = null;
        s.selectedChain = null;
        s.selectedLigand = null;
        s.nearbyResidues = new Set();
        byId('ufv-details')?.classList.remove('show');
        _showOtherSpheres = true;
        const chk = byId('ufv-sphere-chk'); if (chk) chk.checked = true;
        StructureViewer.showLigands = true;    // new structure starts with all ligands shown
        s.analysis.ptmVariantProximity = null; // structure-dependent → recompute on next click
        // Colouring resets to the WINDOW's intended default (Subcellular→topology, Structure→pLDDT,
        // else the configured default) — never the pocket mode. This is what makes re-opening a window
        // restore its focus colouring without a manual reset.
        setColorMode(windowDefaultColorMode());
        byId('ufv-sens-wrap')?.classList.add('ufv-hidden');
    }

    function updateStructureMeta() {
        const st = UFVState.selectedStructure();
        const idx = UFVState.state.selectedStructureIndex;
        // Sync custom structure selector
        const btn = byId('ufv-cs-btn');
        if (btn && st) btn.textContent = stLabel(st);
        document.querySelectorAll('#ufv-cs-drop .ufv-cs-opt').forEach((el, i) => el.classList.toggle('selected', i === idx));
        // Show only relevant color modes for the current structure type
        const isAlphaFold = st?.source === 'AlphaFold';
        const cmDrop = byId('ufv-cm-drop');
        if (cmDrop) {
            cmDrop.querySelector('[data-value="bfactor"]').style.display = isAlphaFold ? 'none' : '';
            cmDrop.querySelector('[data-value="plddt"]').style.display = isAlphaFold ? '' : 'none';
            // Membrane-topology mode only when the entry has topology features.
            const hasTopology = UFVState.state.topology?.length > 0;
            cmDrop.querySelector('[data-value="topology"]').classList.toggle('ufv-hidden', !hasTopology);
            const cur = getColorMode();
            if ((cur === 'bfactor' && isAlphaFold) || (cur === 'plddt' && !isAlphaFold) || (cur === 'topology' && !hasTopology)) {
                setColorMode('default');
            }
        }
        byId('ufv-structure-meta') && (byId('ufv-structure-meta').textContent = '');
        renderStructureSelector();
        renderSequence();
    }

    /**
     * Returns a parenthetical note like " (mapped: 249–444)" when the selected
     * structure only covers a partial stretch of the protein sequence, so users
     * can immediately see why annotations are absent from the uncovered region.
     * Returns '' for AlphaFold, full-coverage structures, or when data is unavailable.
     */
    function getMappedRangeNote() {
        return ''; // mapped-range note removed per request
    }

    // ── Family & Domains window helpers ─────────────────────────────────────────────────────
    function activeDomains() {
        return (UFVState.state.domains || []).filter(d => d.visible !== false);
    }

    // UniProt position → colour for the active *range* domain features. Longer features are laid
    // down first so a shorter, more specific feature (e.g. a region inside a domain) wins on overlap.
    function domainByPos() {
        const map = new Map();
        activeDomains().filter(d => d.isRange)
            .slice().sort((a, b) => (b.endPosition - b.position) - (a.endPosition - a.position))
            .forEach(d => { for (let p = d.position; p <= d.endPosition; p++) map.set(p, d.color); });
        return map;
    }

    // Minimum fraction of a domain's residues that must be MODELLED in the loaded structure for
    // it to be selectable by default. Using mapped (SIFTS) range isn't enough: experimental models
    // often map a wide range but leave big disordered/cytoplasmic stretches unresolved, which would
    // otherwise sit checked-but-invisible. Tunable.
    const DOMAIN_MIN_MODELLED_FRAC = 0.3;

    // Set domain default visibility from what the loaded structure actually resolves. Full-length
    // canonical AlphaFold resolves everything (all selectable); other models only resolve their
    // modelled residues, so a domain lying in an unresolved region starts unselected.
    function applyDomainCoverageDefaults() {
        const s = UFVState.state;
        const st = UFVState.selectedStructure();
        const fullCoverage = !st || (st.source === 'AlphaFold' && !st.isoform);
        const modelled = fullCoverage ? null : new Set(StructureViewer.mappedResidues?.() || []);
        (s.domains || []).forEach(d => { d.visible = domainVisualizable(d, modelled); });
    }

    function domainVisualizable(d, modelled) {
        if (!modelled) return true;            // full-coverage model
        const end = d.endPosition || d.position;
        let n = 0;
        for (let p = d.position; p <= end; p++) if (modelled.has(p)) n++;
        const len = end - d.position + 1;
        return n >= Math.max(1, Math.ceil(len * DOMAIN_MIN_MODELLED_FRAC));
    }

    // Variants belonging to a disease the user selected in the secondary "Disease variants" group.
    // Empty selection ⇒ nothing shown (the group is off/unselected by default).
    function diseaseVariants() {
        const sel = UFVState.state.featDiseases;
        // A variant shows if its disease is selected OR it was individually force-shown, and it isn't hidden.
        const seen = new Set();
        return UFVState.state.variants.filter(v => {
            const k = variantKey(v);
            if (_hiddenVariantKeys.has(k)) return false;
            const inSel = sel && sel.size && (v.diseases || []).some(d => sel.has(d));
            if (!inSel && !_forcedVariantKeys.has(k)) return false;
            if (seen.has(k)) return false; seen.add(k);
            return true;
        });
    }

    // Module-level so dynamically-built sections (e.g. the Disease-variants C button) can call it.
    function setColorProfile(profile) {
        _colorProfile = profile;
        deselectDomainsOnColorSwitch(); // clear domain ticks before turning domain-colouring off
        _domainCartoon = false; // sphere Cs and the domain-backbone C are mutually exclusive
        syncColorProfileButtons();
        applyMode();
    }
    // Rebuild the current window's filter panel from the data model (preserves visibility flags; resets
    // expand/collapse — acceptable for the domain-C mode switch which flips many layers at once).
    function rebuildCurrentWindow() {
        const mode = UFVState.state.currentMode;
        if (mode === 'structure' || mode === 'subcellular') { buildStructureWindow(mode); return; }
        if (mode === 'sites' || mode === 'domains') { buildFeatureWindow(mode); return; }
        if (mode === 'ptm') {
            buildPTMFilters();
            byId('ufv-ptm-panel')?.querySelector('[data-ufv-disease]')?.replaceWith(buildDiseaseVariantSection());
            appendDomainSection('ufv-ptm-panel');
        } else if (mode === 'variant') {
            buildVariantFilters(); // also re-appends the domain section
        }
        buildSiteFilters();
    }
    // C on Family & Domains = a DOMAIN-ONLY view: show all domains, turn every other annotation layer OFF
    // (a "None" for PTMs / sites / disease), paint the off-white backbone by domain. Mutually exclusive
    // with the sphere Cs (those clear it via setColorProfile). Toggling off returns to the window's mode.
    function toggleDomainCartoon() {
        _domainCartoon = !_domainCartoon;
        if (_domainCartoon) {
            const s = UFVState.state;
            (s.domains || []).forEach(d => { d.visible = true; });
            Object.values(s.ptmGroups).forEach(g => setPtmGroupVisible(g, false));
            s.sites.forEach(x => { x.visible = false; });
            s.featDiseases = new Set();
            // rebuildCurrentWindow() re-creates the panel's collapsible sections, which would otherwise
            // snap shut (incl. the Family & Domains section the user just clicked the C in). Snapshot the
            // expand/collapse state by section title and restore it after the rebuild.
            const panelId = currentPanelId();
            const snap = snapshotCollapseStates(panelId);
            rebuildCurrentWindow();
            restoreCollapseStates(panelId, snap);
        }
        syncColorProfileButtons();
        applyMode();
    }

    // Selecting a domain in a Family & Domains checklist turns on backbone-by-domain colouring so the
    // choice is IMMEDIATELY visible (the old behaviour only painted once the C button was pressed). Unlike
    // the C button this is non-destructive — it doesn't switch other layers (PTM/site/disease) off.
    function enableDomainCartoon() {
        _domainCartoon = true;
        syncColorProfileButtons();
        applyMode();
    }
    function disableDomainCartoon() {
        _domainCartoon = false;
        syncColorProfileButtons();
        applyMode();
    }
    // A colour-mode / sphere-C switch turns off domain-backbone colouring. Clear the now-meaningless domain
    // selections (and refresh the checklists) so the boxes don't stay ticked while nothing is painted by
    // domain. No-op unless domain colouring was actually on.
    function deselectDomainsOnColorSwitch() {
        if (!_domainCartoon) return;
        (UFVState.state.domains || []).forEach(d => { d.visible = false; });
        refreshAllDomainSections();
    }

    // The active right-panel id for the current window mode (used to scope collapse-state snapshots).
    function currentPanelId() {
        const m = UFVState.state.currentMode;
        return m === 'ptm' ? 'ufv-ptm-panel'
             : m === 'variant' ? 'ufv-var-panel'
             : m === 'domains' ? 'ufv-dom-panel'
             : 'ufv-feat-panel'; // sites / structure / subcellular
    }
    // Map collapsible-section title → expanded? for every section in a panel (titles are unique per panel).
    const _sectionTitle = box => box.querySelector('.ufv-collapsible-hdr > span:not(.ufv-collapsible-chevron):not(.ufv-section-actions)')?.textContent;
    function snapshotCollapseStates(panelId) {
        const map = new Map();
        byId(panelId)?.querySelectorAll('.ufv-collapsible').forEach(box => {
            const title = _sectionTitle(box);
            if (title) map.set(title, !box.querySelector('.ufv-collapsible-body')?.classList.contains('ufv-collapsed'));
        });
        return map;
    }
    function restoreCollapseStates(panelId, map) {
        byId(panelId)?.querySelectorAll('.ufv-collapsible').forEach(box => {
            const title = _sectionTitle(box);
            if (!title || !map.has(title)) return;
            const body = box.querySelector('.ufv-collapsible-body');
            const chev = box.querySelector('.ufv-collapsible-chevron');
            const expand = map.get(title), collapsed = body?.classList.contains('ufv-collapsed');
            if (expand && collapsed) { body.classList.remove('ufv-collapsed'); if (chev) chev.innerHTML = '&#9660;'; }
            else if (!expand && !collapsed) { body.classList.add('ufv-collapsed'); if (chev) chev.innerHTML = '&#9654;'; }
        });
    }

    // Refresh just the Family & Domains checklist in every panel that hosts one (preserves other
    // sections' expand/collapse state, unlike a full panel rebuild).
    function refreshAllDomainSections() {
        ['ufv-ptm-panel', 'ufv-var-panel', 'ufv-feat-panel', 'ufv-dom-panel'].forEach(pid => {
            const existing = byId(pid)?.querySelector('[data-ufv-domains]');
            if (!existing) return;
            const wasExpanded = !existing.querySelector('.ufv-collapsible-body')?.classList.contains('ufv-collapsed');
            existing.replaceWith(buildDomainSection(wasExpanded));
        });
    }
    function syncColorProfileButtons() {
        // When the domain-backbone C is active no sphere C is lit (they're mutually exclusive).
        const sc = !_domainCartoon;
        const isPtm = sc && _colorProfile === 'ptm';
        const isDis = sc && _colorProfile === 'disease';
        const isCons = sc && _colorProfile === 'consequence';
        const isProv = sc && _colorProfile === 'provenance';
        const isSites = sc && _colorProfile === 'sites';
        byId('ufv-brush-ptm')?.classList.toggle('active', isPtm);
        byId('ufv-brush-ptmsite')?.classList.toggle('active', isSites); // PTM-window Sites C tracks 'sites', not 'ptm'
        byId('ufv-brush-vptm')?.classList.toggle('active', isPtm);
        byId('ufv-brush-var')?.classList.toggle('active', isDis);
        byId('ufv-brush-dis')?.classList.toggle('active', isCons);
        byId('ufv-brush-prov')?.classList.toggle('active', isProv);
        byId('ufv-brush-vsite')?.classList.toggle('active', isSites);
        // Dynamically-built section C buttons (feature/PTM windows) — tagged with data-profile.
        // 'domains' tracks the backbone-colouring toggle, not the sphere _colorProfile.
        document.querySelectorAll('.ufv-brush-dyn').forEach(b => b.classList.toggle('active',
            b.dataset.profile === 'domains' ? _domainCartoon : (sc && b.dataset.profile === _colorProfile)));
    }

    function diseaseSpheres() {
        const sel = UFVState.state.featDiseases;
        return diseaseVariants().map(v => {
            let color;
            if (_colorProfile === 'consequence' || _colorProfile === 'provenance' || _colorProfile === 'sites') {
                color = v.consequenceColor;
            } else {
                // Use color of first disease that's selected; grey for "Other" variants (no named disease match)
                const vDiseases = v.diseases || [];
                const activeDis = vDiseases.find(d => d !== '__other__' && (!sel || !sel.size || sel.has(d)));
                if (activeDis) {
                    color = _diseaseColorMap.get(activeDis) || v.consequenceColor;
                } else {
                    // Variant has no named disease in the active set → "Other" → grey
                    color = '#9e9e9e';
                }
            }
            return { position: v.position, color, hover: v };
        });
    }

    // Secondary "Disease variants" collapsible: a checklist of diseases (off by default) with
    // All/None. Selecting diseases shows their variants as spheres. Shared by the PTM, Functional-
    // features and Family & Domains windows.
    function buildDiseaseVariantSection(expanded = false) {
        const s = UFVState.state;
        const { names, ds } = diseaseNamesToShow();
        const { box, body } = makeCollapsibleSection('Disease variants', {
            onAll: () => { s.featDiseases = new Set(names); buildAllSecondaryDiseaseSections(); setColorProfile('disease'); }, // showing all diseases → colour by disease
            onNone: () => { s.featDiseases = new Set(); buildAllSecondaryDiseaseSections(); reapply(); },
            onBrush: () => setColorProfile('disease'),
            brushProfile: 'disease',
            brushTitle: 'Colour disease-variant spheres by disease',
            expanded,
        });
        box.setAttribute('data-ufv-disease', '1');
        if (!names.length) {
            const empty = document.createElement('div'); empty.className = 'ufv-empty'; empty.textContent = 'No disease-associated variants.';
            body.appendChild(empty);
        } else {
            names.forEach(name => {
                const meta = ds[name] || { color: '#e53935', count: 0 };
                body.appendChild(makeSecondaryDiseaseFilter(name, meta.color, meta.count));
            });
        }
        return box;
    }

    // A disease row in the SECONDARY disease section (featDiseases model): parent toggles the whole
    // disease; expand (chevron) to per-variant rows, each with a checkbox (individually hide/show via
    // _hiddenVariantKeys) and a zoom magnifier. Mirrors the Variant window's disease dropdown so disease
    // variants are individually selectable in every window that hosts this section.
    function makeSecondaryDiseaseFilter(name, color, count) {
        const s = UFVState.state;
        const wrap = document.createElement('div');
        wrap.className = 'ufv-filter-group';
        const children = document.createElement('div');
        children.className = 'ufv-filter-children ufv-collapsed';
        const seen = new Set();
        s.variants.filter(v => (v.diseases || []).includes(name)).sort((a, b) => a.position - b.position).forEach(v => {
            const key = variantKey(v);
            if (seen.has(key)) return; seen.add(key);
            const rowEl = document.createElement('label');
            rowEl.className = 'ufv-filter-item ufv-disease-var';
            const cb = document.createElement('input');
            cb.type = 'checkbox'; cb.className = 'ufv-lig-eye';
            // Checkbox reflects EFFECTIVE visibility (disease selected or force-shown, and not hidden), and
            // toggling it directly controls the sphere — checking force-shows it even when the parent
            // disease is off, unchecking hides it. Matches the Variant window's per-variant model.
            cb.checked = (s.featDiseases.has(name) || _forcedVariantKeys.has(key)) && !_hiddenVariantKeys.has(key);
            cb.addEventListener('change', () => {
                if (cb.checked) { _forcedVariantKeys.add(key); _hiddenVariantKeys.delete(key); setColorProfile('disease'); } // showing a disease variant → colour spheres by disease
                else { _hiddenVariantKeys.add(key); _forcedVariantKeys.delete(key); reapply(); }
            });
            const lbl = document.createElement('span');
            lbl.className = 'ufv-filter-label';
            lbl.textContent = `${v.wildType || ''}${v.position}${v.mutant || ''}`;
            if (v.consequenceColor) lbl.style.color = v.consequenceColor;
            rowEl.append(cb, lbl, makeZoomBtn(`Zoom to ${v.position}`, () => onClick({ position: v.position }, null, s.selectedChain ?? null)));
            children.appendChild(rowEl);
        });
        const top = makeFilterItem(name, color, count, s.featDiseases.has(name), checked => {
            if (checked) { s.featDiseases.add(name); setColorProfile('disease'); } // selecting a disease → colour spheres by disease, not consequence
            else { s.featDiseases.delete(name); reapply(); }
        });
        const chevron = document.createElement('button');
        chevron.type = 'button'; chevron.className = 'ufv-group-chevron'; chevron.innerHTML = '&#9654;';
        chevron.addEventListener('click', e => { e.preventDefault(); e.stopPropagation(); const c = children.classList.toggle('ufv-collapsed'); chevron.innerHTML = c ? '&#9654;' : '&#9660;'; });
        top.appendChild(chevron);
        wrap.append(top, children);
        return wrap;
    }

    // Secondary "Family & Domains" collapsible — a domain checklist with a C that paints the backbone
    // by domain. Shared by the PTM / Variant windows so every window exposes the domain layer. Domain
    // toggles drive domainByPos (which respects d.visible), so they take effect when the backbone-C is on.
    function buildDomainSection(expanded = false) {
        const s = UFVState.state;
        const sec = makeCollapsibleSection('Family & Domains', {
            onAll: () => { (s.domains || []).forEach(d => { d.visible = true; }); refreshAllDomainSections(); enableDomainCartoon(); },
            onNone: () => { (s.domains || []).forEach(d => { d.visible = false; }); refreshAllDomainSections(); disableDomainCartoon(); },
            onBrush: toggleDomainCartoon, brushProfile: 'domains', brushTitle: 'Colour the backbone by domain (and hide other layers)',
            expanded,
        });
        sec.box.setAttribute('data-ufv-domains', '1');
        (s.domains || []).forEach(item => sec.body.appendChild(
            makeFilterItem(featureLabel(item), item.color, '', item.visible !== false,
                checked => { item.visible = checked; if (checked) enableDomainCartoon(); else applyMode(); },
                undefined, { zoomPos: item.position, zoomRange: item.isRange ? [item.position, item.endPosition] : null, chain: s.selectedChain ?? null })));
        return sec.box;
    }

    // Append (or refresh) a Family & Domains section in a static panel that doesn't already have one.
    // Preserves the existing section's expand/collapse state so the domain C (which rebuilds the window)
    // doesn't snap the section shut under the user's cursor.
    function appendDomainSection(panelId) {
        if (!(UFVState.state.domains || []).length) return;
        const panel = byId(panelId);
        const existing = panel?.querySelector('[data-ufv-domains]');
        const wasExpanded = existing && !existing.querySelector('.ufv-collapsible-body')?.classList.contains('ufv-collapsed');
        existing?.remove();
        panel?.appendChild(buildDomainSection(wasExpanded));
    }

    // Re-draw the disease-variant checklist in every window that hosts one (so All/None and
    // cross-window state stay in sync), then re-render the viewer.
    function buildAllSecondaryDiseaseSections() {
        ['ufv-ptm-panel', 'ufv-feat-panel', 'ufv-dom-panel'].forEach(pid => {
            const panel = byId(pid);
            const existing = panel?.querySelector('[data-ufv-disease]');
            if (!existing) return;
            // Preserve the section's expand/collapse state across the rebuild (so All/None doesn't collapse it).
            const wasExpanded = !existing.querySelector('.ufv-collapsible-body')?.classList.contains('ufv-collapsed');
            existing.replaceWith(buildDiseaseVariantSection(wasExpanded));
        });
    }

    // Re-render using the fast PTM path when in the PTM window, else a full applyMode.
    function reapply() {
        if (UFVState.state.currentMode === 'ptm') applyPTMMode();
        else applyMode();
    }

    // A Ligands collapsible built dynamically for the feature windows (the static PTM/Variant
    // panels already have their own). Shows the model's ligands when any are present.
    function buildDynamicLigandSection() {
        const s = UFVState.state;
        if (!s.ligands || !s.ligands.length) return null;
        const { box, body } = makeCollapsibleSection('Ligands', {
            onAll: () => ligandsSetAll(true),
            onNone: () => ligandsSetAll(false),
        });
        const groups = new Map();
        s.ligands.forEach(l => { if (!groups.has(l.resn)) groups.set(l.resn, []); groups.get(l.resn).push(l); });
        if (ligandsHaveIdentity()) body.appendChild(makeLigandThresholdRow());
        groups.forEach((copies, resn) => body.appendChild(makeLigandGroup(resn, copies)));
        return box;
    }

    // Label for a primary feature row: three-letter AA + position for single residues, a span for
    // ranges; then the description.
    function featureLabel(item) {
        const seq = UFVState.state.sequence;
        const single = !item.endPosition || item.endPosition === item.position;
        const aa = single ? (AA1TO3[seq?.[item.position - 1]] || '') : '';
        const loc = single ? `${aa} ${item.position}`.trim() : `${item.position}–${item.endPosition}`;
        return `${loc}: ${item.description}`;
    }

    // Per-window default visibility + dynamic panel build. Feature windows show their own group as
    // the selectable primary list and start PTMs / disease variants collapsed and off.
    // Set a PTM group's master flag AND every item in it (item.visible is what actually drives the
    // sphere — see showPTMs/activePtms), so toggling a layer off truly hides its individual PTMs.
    function setPtmGroupVisible(g, visible) { g.visible = visible; g.items.forEach(it => { it.visible = visible; }); }

    function setupWindowDefaults(mode) {
        const s = UFVState.state;
        s.featDiseases = new Set(); // disease variants always start unselected
        _domainCartoon = false;     // domain-backbone colouring is opt-in per window via the Domains C
        if (mode === 'structure' || mode === 'subcellular') {
            // Structure / Subcellular: the full layer set, ALL off — a clean structure overview the user
            // opts into. Structure colours by pLDDT confidence; Subcellular by membrane topology.
            Object.values(s.ptmGroups).forEach(g => setPtmGroupVisible(g, false));
            s.sites.forEach(x => { x.visible = false; });
            (s.domains || []).forEach(d => { d.visible = false; });
            setColorMode(mode === 'subcellular' ? 'topology' : 'plddt');
            buildStructureWindow(mode);
            return;
        }
        if (mode === 'sites' || mode === 'domains') {
            // Feature windows: own group ON (primary), PTMs + disease variants OFF (secondary).
            Object.values(s.ptmGroups).forEach(g => setPtmGroupVisible(g, false));
            s.sites.forEach(x => { x.visible = (mode === 'sites'); });
            (s.domains || []).forEach(d => { d.visible = (mode === 'domains'); });
            buildFeatureWindow(mode);
        } else {
            // PTM / Variant windows: restore their defaults (a prior feature window may have
            // flipped PTM-group visibility or turned sites on).
            Object.values(s.ptmGroups).forEach(g => setPtmGroupVisible(g, true));
            s.sites.forEach(x => { x.visible = false; });
            if (mode === 'ptm') {
                buildPTMFilters();
                // PTM window also offers disease variants + Family & Domains as collapsed, off-by-default groups.
                const panel = byId('ufv-ptm-panel');
                panel?.querySelector('[data-ufv-disease]')?.remove();
                panel?.appendChild(buildDiseaseVariantSection());
                appendDomainSection('ufv-ptm-panel');
            }
            buildSiteFilters();
        }
    }

    // A chevron section with optional All/None/C buttons. Collapsed by default; pass expanded:true for a
    // window's PRIMARY (focused) category so it opens expanded but is still collapsible — and uses the
    // SAME header chrome/font as every other section (the old <h3> primary header was a different size).
    // onBrush: makes a "C" colour-profile button; brushProfile tags it so syncColorProfileButtons can
    // light it when that profile is active.
    function makeCollapsibleSection(title, { onAll, onNone, onBrush, brushProfile, brushTitle, expanded } = {}) {
        const box = document.createElement('div'); box.className = 'ufv-collapsible';
        const hdr = document.createElement('div'); hdr.className = 'ufv-collapsible-hdr';
        const chev = document.createElement('span'); chev.className = 'ufv-collapsible-chevron'; chev.innerHTML = expanded ? '&#9660;' : '&#9654;';
        const t = document.createElement('span'); t.textContent = title;
        const acts = document.createElement('div'); acts.className = 'ufv-section-actions';
        const mkBtn = (label, fn) => { const b = document.createElement('button'); b.className = 'ufv-section-btn'; b.textContent = label; b.addEventListener('click', e => { e.stopPropagation(); fn(); }); return b; };
        if (onAll) acts.appendChild(mkBtn('All', onAll));
        if (onNone) acts.appendChild(mkBtn('None', onNone));
        if (onBrush) {
            const b = mkBtn('C', onBrush);
            b.classList.add('ufv-brush-btn', 'ufv-brush-dyn');
            b.dataset.profile = brushProfile || '';
            if (brushProfile && _colorProfile === brushProfile) b.classList.add('active');
            if (brushTitle) b.title = brushTitle;
            acts.appendChild(b);
        }
        hdr.append(chev, t, acts);
        const body = document.createElement('div'); body.className = 'ufv-collapsible-body' + (expanded ? '' : ' ufv-collapsed');
        hdr.addEventListener('click', e => { if (e.target.closest('button')) return; const c = body.classList.toggle('ufv-collapsed'); chev.innerHTML = c ? '&#9654;' : '&#9660;'; });
        box.append(hdr, body);
        return { box, body };
    }

    // Build the right-panel content for a feature window: primary selectable list on top, then
    // collapsed (off) PTMs and disease-variant sections, then ligands.
    function buildFeatureWindow(mode) {
        const s = UFVState.state;
        const panel = byId(mode === 'sites' ? 'ufv-feat-panel' : 'ufv-dom-panel');
        if (!panel) return;
        panel.textContent = '';
        const items = mode === 'sites' ? s.sites : (s.domains || []);

        // Primary (focused) category — collapsible, starts EXPANDED, with All/None and a C colour-enforce
        // button (Functional features → colour spheres by site). Uses the SAME header chrome as the
        // secondary sections so the focused-category title is the same size across all windows.
        let primary;
        const setAllPrimary = (vis) => { items.forEach(x => { x.visible = vis; }); primary.body.querySelectorAll('input[type="checkbox"]').forEach(cb => { cb.checked = vis; }); applyMode(); };
        primary = makeCollapsibleSection(mode === 'sites' ? 'Functional features' : 'Family & Domains', {
            onAll: () => setAllPrimary(true),
            onNone: () => setAllPrimary(false),
            ...(mode === 'sites'
                ? { onBrush: () => setColorProfile('sites'), brushProfile: 'sites', brushTitle: 'Colour site spheres by functional site' }
                : { onBrush: toggleDomainCartoon, brushProfile: 'domains', brushTitle: 'Colour the backbone by domain' }),
            expanded: true,
        });
        if (!items.length) {
            const empty = document.createElement('div'); empty.className = 'ufv-empty';
            empty.textContent = mode === 'sites' ? 'No active/binding/metal-site features for this protein.'
                : 'No domain / region / repeat features for this protein.';
            primary.body.appendChild(empty);
        } else {
            items.forEach(item => primary.body.appendChild(
                makeFilterItem(featureLabel(item), item.color, '', item.visible !== false, checked => { item.visible = checked; applyMode(); },
                    undefined, { zoomPos: item.position, zoomRange: item.isRange ? [item.position, item.endPosition] : null, chain: s.selectedChain ?? null })));
        }
        panel.appendChild(primary.box);

        // Secondary: PTMs (collapsed, off by default). All/None sync item.visible (authoritative) in place.
        let ptm;
        const setAllPtm = (vis) => { Object.values(s.ptmGroups).forEach(g => setPtmGroupVisible(g, vis)); ptm.body.querySelectorAll('input[type="checkbox"]').forEach(cb => { cb.checked = vis; }); applyMode(); };
        ptm = makeCollapsibleSection('PTMs', {
            onAll: () => setAllPtm(true),
            onNone: () => setAllPtm(false),
            onBrush: () => setColorProfile('ptm'),
            brushProfile: 'ptm',
            brushTitle: 'Colour PTM spheres by PTM type',
        });
        Object.values(s.ptmGroups).forEach(g => ptm.body.appendChild(
            makeExpandableFilter(g.category, g.color, g.items, !!g.visible, checked => {
                g.visible = checked; g.items.forEach(it => it.visible = checked); applyMode();
            }, item => `Residue ${item.position}: ${item.description}`, applyMode)));
        panel.appendChild(ptm.box);

        // Secondary: disease variants (collapsed, off by default) — a checklist of diseases.
        panel.appendChild(buildDiseaseVariantSection());

        // Secondary: the OTHER feature type, so every feature window exposes BOTH sites and domains.
        // (Sites window → Family & Domains via backbone-C; Domains window → Functional features spheres.)
        const otherItems = mode === 'sites' ? (s.domains || []) : s.sites;
        if (otherItems.length) {
            let otherSec;
            const otherIsDomains = mode === 'sites'; // in the Sites window the "other" layer is domains
            // Selecting domains here should colour the backbone immediately (same fix as the main section).
            const afterOtherOn = () => { if (otherIsDomains) enableDomainCartoon(); else applyMode(); };
            const setAllOther = (vis) => {
                otherItems.forEach(x => { x.visible = vis; });
                otherSec.body.querySelectorAll('input[type="checkbox"]').forEach(cb => { cb.checked = vis; });
                vis ? afterOtherOn() : (otherIsDomains ? disableDomainCartoon() : applyMode());
            };
            otherSec = makeCollapsibleSection(otherIsDomains ? 'Family & Domains' : 'Functional features', {
                onAll: () => setAllOther(true),
                onNone: () => setAllOther(false),
                ...(otherIsDomains
                    ? { onBrush: toggleDomainCartoon, brushProfile: 'domains', brushTitle: 'Colour the backbone by domain (and hide other layers)' }
                    : { onBrush: () => setColorProfile('sites'), brushProfile: 'sites', brushTitle: 'Colour site spheres by functional site' }),
            });
            otherItems.forEach(item => otherSec.body.appendChild(
                makeFilterItem(featureLabel(item), item.color, '', item.visible !== false,
                    checked => { item.visible = checked; checked ? afterOtherOn() : applyMode(); },
                    undefined, { zoomPos: item.position, zoomRange: item.isRange ? [item.position, item.endPosition] : null, chain: s.selectedChain ?? null })));
            panel.appendChild(otherSec.box);
        }

        // Ligands present in the model (e.g. SwissModel/AlphaFill cofactors) — shown whenever any
        // exist, so the window reflects everything in the system, not just the protein.
        const ligSec = buildDynamicLigandSection();
        if (ligSec) panel.appendChild(ligSec);
    }

    // Structure / Subcellular window: the FULL layer set (sites, PTMs, disease variants, domains,
    // ligands) as collapsed sections that all start OFF, so the user begins from a clean structure
    // overview and opts into each layer. This is the unified "every window has every layer" model;
    // the only difference from a feature window is that there is no pre-selected primary group.
    function buildStructureWindow(mode) {
        const s = UFVState.state;
        const panel = byId('ufv-feat-panel');
        if (!panel) return;
        panel.textContent = '';
        const rebuild = () => buildStructureWindow(mode);

        // Generic collapsible layer section over a list of {visible,color,...} items. All/None update the
        // section's checkboxes IN PLACE (not a panel rebuild) so they don't collapse every section.
        const layerSection = (title, items, labelFn, extra) => {
            if (!items.length) return;
            const isDomains = extra && extra.brushProfile === 'domains'; // domains colour the backbone on select
            let sec;
            const setAll = (vis) => { items.forEach(x => { x.visible = vis; }); sec.body.querySelectorAll('input[type="checkbox"]').forEach(cb => { cb.checked = vis; }); isDomains ? (vis ? enableDomainCartoon() : disableDomainCartoon()) : applyMode(); };
            sec = makeCollapsibleSection(title, {
                onAll: () => setAll(true),
                onNone: () => setAll(false),
                ...(extra || {}),
            });
            items.forEach(x => sec.body.appendChild(
                makeFilterItem(labelFn(x), x.color, x.count != null ? x.count : '', x.visible === true,
                    checked => { x.visible = checked; isDomains && checked ? enableDomainCartoon() : applyMode(); },
                    undefined, x.position != null ? { zoomPos: x.position, zoomRange: x.isRange ? [x.position, x.endPosition] : null, chain: s.selectedChain ?? null } : undefined)));
            panel.appendChild(sec.box);
        };

        layerSection('Functional features', s.sites || [], featureLabel,
            { onBrush: () => setColorProfile('sites'), brushProfile: 'sites', brushTitle: 'Colour site spheres by functional site' });
        // PTMs as expandable categories (individual sites + zoom), consistent with the PTM window.
        const ptmGroups = Object.values(s.ptmGroups || {});
        if (ptmGroups.length) {
            let ptmSec;
            const setAllPtm = (vis) => { ptmGroups.forEach(g => { g.visible = vis; g.items.forEach(it => { it.visible = vis; }); }); ptmSec.body.querySelectorAll('input[type="checkbox"]').forEach(cb => { cb.checked = vis; }); applyMode(); };
            ptmSec = makeCollapsibleSection('PTMs', {
                onAll: () => setAllPtm(true),
                onNone: () => setAllPtm(false),
                onBrush: () => setColorProfile('ptm'), brushProfile: 'ptm', brushTitle: 'Colour PTM spheres by PTM type',
            });
            ptmGroups.forEach(g => ptmSec.body.appendChild(
                makeExpandableFilter(g.category, g.color, g.items, g.visible === true, checked => {
                    g.visible = checked; g.items.forEach(it => it.visible = checked); applyMode();
                }, item => `Residue ${item.position}: ${item.description}`, applyMode)));
            panel.appendChild(ptmSec.box);
        }
        panel.appendChild(buildDiseaseVariantSection());       // disease variants (off by default)
        layerSection('Family & Domains', s.domains || [], featureLabel,
            { onBrush: toggleDomainCartoon, brushProfile: 'domains', brushTitle: 'Colour the backbone by domain' });
        const ligSec = buildDynamicLigandSection();
        if (ligSec) panel.appendChild(ligSec);
    }

    // Non-canonical isoform AlphaFold models number residues by the isoform sequence. We carry a
    // real canonical<->isoform mapping (segmented mappedRanges built from the isoform's VSP edits),
    // so annotations now render at the correct residues through the normal mapping path; the only
    // isoform-specific handling left is fetching the isoform's own PAE for pocket analysis.
    function applyMode() {
        const s = UFVState.state;
        if (!StructureViewer.viewer) return;
        const cm = s.currentMode;
        const mode = getColorMode();
        // Additive AND filter across the three axes, PLUS any variants force-shown via an individual
        // disease-dropdown row (those render even when their disease axis is off), MINUS hidden ones.
        const _base = DataProcessor.filterVariants(s.variants, s.activeConsequences, s.activeProvenances, s.activeDiseases);
        const _baseKeys = new Set(_base.map(variantKey));
        const _forced = _forcedVariantKeys.size
            ? s.variants.filter(v => _forcedVariantKeys.has(variantKey(v)) && !_baseKeys.has(variantKey(v)))
            : [];
        const filteredVariants = _base.concat(_forced).filter(v => !_hiddenVariantKeys.has(variantKey(v)));
        // Cartoon colouring. The Family & Domains window colours by domain ranges; every other
        // window uses the chosen colour mode. defer=true skips the intermediate render.
        if (cm === 'domains' || _domainCartoon) {
            StructureViewer.applyCartoonColoring('domains', { domainByPos: domainByPos() }, true);
        } else {
            StructureViewer.applyCartoonColoring(mode, {
                ptms: activePtms(),
                variants: filteredVariants,
                hotspots: s.analysis.hotspots,
                hotspotsByChain: s.analysis.hotspotsByChain,
                distantContacts: s.analysis.distantContacts,
                distantContactsByChain: s.analysis.distantContactsByChain,
                alphaMissense: s.analysis.alphaMissense,
                residueBurden: s.analysis.residueBurden,
                pocketByPos: mode === 'prism' ? filteredPocketByPos() : null,
                topologyByPos: mode === 'topology' ? topologyByPos() : null,
            }, true);
        }
        const rangeNote = getMappedRangeNote();
        const siteList = activeSites();
        const sitePositions = siteList.flatMap(x => x.endPosition && x.endPosition !== x.position ? [x.position, x.endPosition] : [x.position]);
        if (cm === 'ptm') {
            const dsph = diseaseSpheres();
            const n = StructureViewer.showPTMs(s.ptms, s.ptmGroups, siteList, dsph);
            // Copy reflects EVERYTHING shown in the PTM window — PTM sites, functional sites AND any
            // selected disease-variant spheres (previously the diseases were drawn but omitted from copy).
            s.displayedPositions = [...new Set([
                ...activePtms().flatMap(p => p.endPosition && p.endPosition !== p.position ? [p.position, p.endPosition] : [p.position]),
                ...sitePositions,
                ...dsph.map(sp => sp.position),
            ])];
            byId('ufv-count-text').textContent = `${n} PTM site${n === 1 ? '' : 's'}${rangeNote}`;
        } else if (cm === 'variant') {
            const coPtms = activeCoDisplayPtms();
            // Sphere colour follows the active "C" profile the user picked on a section header.
            // showVariants reads v.consequenceColor, so we overwrite that field per profile.
            let variantsColored = filteredVariants;
            if (_colorProfile === 'disease') {
                variantsColored = filteredVariants.map(v => {
                    const ad = s.activeDiseases;
                    const firstActiveDis = (v.diseases || []).find(d => !ad || !ad.size || ad.has(d));
                    const dc = (firstActiveDis && _diseaseColorMap.get(firstActiveDis)) || '#9e9e9e';
                    return { ...v, consequenceColor: dc };
                });
            } else if (_colorProfile === 'provenance') {
                variantsColored = filteredVariants.map(v => ({
                    ...v,
                    consequenceColor: DataProcessor.PROVENANCE_CATEGORIES[v.provenance]?.color || '#9e9e9e',
                }));
            } // 'consequence' / 'sites' / 'ptm' ⇒ keep natural consequenceColor
            const r = StructureViewer.showVariants(variantsColored, coPtms, siteList);
            s.displayedPositions = Array.from(new Set([
                ...filteredVariants.map(v => v.position),
                ...coPtms.flatMap(p => p.endPosition && p.endPosition !== p.position ? [p.position, p.endPosition] : [p.position]),
                ...sitePositions,
            ]));
            const ptmNote = r.ptmCount ? `, ${r.ptmCount} PTM site${r.ptmCount === 1 ? '' : 's'}` : '';
            byId('ufv-count-text').textContent = `${r.varCount} variants at ${r.posCount} positions${ptmNote}${rangeNote}`;
        } else {
            // Feature windows ('sites' / 'domains'): primary feature spheres on top, plus any
            // toggled-on secondary PTMs / disease variants. Domain ranges are already painted on
            // the cartoon above; only single-residue domain features get a sphere.
            const spheres = [];
            if (cm === 'domains') {
                activeDomains().filter(d => !d.isRange).forEach(d => spheres.push({ position: d.position, color: d.color, hover: { position: d.position, description: d.description } }));
            } else {
                siteList.forEach(x => spheres.push({ position: x.position, endPosition: x.endPosition, color: x.color, hover: { position: x.position, description: x.description } }));
            }
            activePtms().forEach(p => spheres.push({ position: p.position, endPosition: p.endPosition, color: p.color, hover: p }));
            diseaseSpheres().forEach(sp => spheres.push(sp));
            StructureViewer.showAnnotationSpheres(spheres);
            s.displayedPositions = [...new Set(spheres.flatMap(sp => sp.endPosition && sp.endPosition !== sp.position ? [sp.position, sp.endPosition] : [sp.position]))];
            if (cm === 'structure' || cm === 'subcellular') {
                byId('ufv-count-text').textContent = `${spheres.length} annotation${spheres.length === 1 ? '' : 's'}${rangeNote}`;
            } else {
                const primaryN = cm === 'domains' ? activeDomains().length : siteList.length;
                const word = cm === 'domains' ? 'domain' : 'site';
                byId('ufv-count-text').textContent = `${primaryN} ${word}${primaryN === 1 ? '' : 's'}${rangeNote}`;
            }
        }
        // Preserve an active focus across coloring changes: re-enter focus on the selected
        // ligand or residue instead of dropping back to the full sphere view.
        if (s.selectedLigand && StructureViewer.currentStructure) {
            const nb = StructureViewer.focusLigand(s.selectedLigand.resn, s.selectedLigand.resi, s.selectedLigand.chain, { showOtherSpheres: _showOtherSpheres, rezoom: false , annotatedResidues: buildAnnotationMap() });
            if (nb) s.nearbyResidues = nb;
        } else if (s.selectedResidue != null && StructureViewer.currentStructure) {
            // rezoom:false — re-applying focus after a filter/colour change must not yank the
            // camera; showOtherSpheres respects the header toggle (so co-displayed PTM spheres
            // appear/disappear correctly while a residue is focused).
            const nearby = StructureViewer.focusResidue(s.selectedResidue, s.selectedChain, { annotatedResidues: buildAnnotationMap() }, { showOtherSpheres: _showOtherSpheres, rezoom: false });
            if (nearby) s.nearbyResidues = nearby;
        }
        renderLegend(mode);
        renderSequence();
    }

    // Fast PTM sphere refresh — skips cartoon rebuild when only sphere visibility changes.
    function applyPTMMode() {
        const s = UFVState.state;
        if (!StructureViewer.viewer) return;
        const n = StructureViewer.refreshPTMDisplay(s.ptms, s.ptmGroups, activeSites(), diseaseSpheres());
        if (n === false) {
            // Focus mode active — need a full rebuild to restore sticks
            applyMode();
            return;
        }
        const sitePositions = activeSites().flatMap(x => x.endPosition && x.endPosition !== x.position ? [x.position, x.endPosition] : [x.position]);
        // Include selected disease-variant spheres too — this fast path is what runs when the user toggles
        // diseases in the PTM window, and copy reads displayedPositions (previously only PTMs + sites).
        s.displayedPositions = [...new Set([
            ...activePtms().flatMap(p => p.endPosition && p.endPosition !== p.position ? [p.position, p.endPosition] : [p.position]),
            ...sitePositions,
            ...diseaseSpheres().map(sp => sp.position),
        ])];
        byId('ufv-count-text').textContent = `${n} PTM site${n === 1 ? '' : 's'}${getMappedRangeNote()}`;
        renderLegend(getColorMode());
        renderSequence();
    }

    function activePtms() {
        const s = UFVState.state;
        // item.visible is authoritative (matches showPTMs) so individually-checked PTMs show even when
        // their group's master toggle is off (Structure/Subcellular windows).
        return s.ptms.filter(p => p.visible !== false && s.ptmGroups[p.category]);
    }

    // UniProt position → topology colour, built from the membrane-topology segments.
    function topologyByPos() {
        const map = new Map();
        (UFVState.state.topology || []).forEach(seg => {
            for (let p = seg.start; p <= seg.end; p++) map.set(p, seg.color);
        });
        return map;
    }

    /**
     * Lazily compute the constraint-pocket analysis for the current structure.
     * Fetches the AlphaFold PAE matrix on demand (only for AlphaFold models — PAE indices
     * align with UniProt positions there) to gate the spatial weights, then runs the heavy
     * permutation test in a Web Worker (off the main thread) with a synchronous fallback.
     */
    async function ensurePocketAnalysis() {
        const s = UFVState.state;
        const st = UFVState.selectedStructure();
        if (!st || !StructureViewer.viewer) return s.analysis.prism;
        if (s.analysis.prism) return s.analysis.prism; // already computed for this structure
        const requestedId = s.uniprotId;
        showLoading('Computing constraint pockets…');
        try {
            let pae = null;
            if (st.source === 'AlphaFold') {
                // Isoform models have their own PAE file (AF-<accession-N>-…); fetch by the
                // isoform accession so PAE indices align with the isoform residue order. Falls
                // back to geometry-only pockets if that PAE isn't available.
                try { pae = await UFVApi.getPaeMatrix(st.isoform || requestedId, st.version); } catch (_) {}
            }
            if (s.uniprotId !== requestedId || UFVState.selectedStructure() !== st) return s.analysis.prism;
            const geometry = StructureViewer.residueGeometry();
            const res = await runPocketAnalysis(geometry, s.amMap, s.sequence, pae);
            if (s.uniprotId !== requestedId || UFVState.selectedStructure() !== st) return s.analysis.prism;
            s.analysis.prism = res;
        } catch (err) {
            s.analysis.prism = { byPos: new Map(), reason: err?.message || 'Computation failed.' };
        } finally {
            byId('ufv-loading').classList.add('hidden');
        }
        return s.analysis.prism;
    }

    // ---- Constraint-pocket compute: short-lived Web Worker with a synchronous fallback ----
    // A fresh worker is created for each compute and ALWAYS terminated when it finishes, so the
    // extension never holds a background worker thread once the analysis is done — no persistent
    // CPU/memory footprint.  A page CSP on uniprot.org could block worker creation; if so (or on
    // any worker error) we mark it broken and compute synchronously instead.
    let _workerBroken = false;

    async function runPocketAnalysis(geometry, amMap, sequence, pae) {
        if (!_workerBroken && typeof chrome !== 'undefined' && chrome.runtime?.getURL) {
            let w = null;
            try {
                w = new Worker(chrome.runtime.getURL('pocket-worker.js'));
                const worker = w;
                return await new Promise((resolve, reject) => {
                    const timer = setTimeout(() => reject(new Error('worker timeout')), 60000);
                    worker.onmessage = (e) => {
                        clearTimeout(timer);
                        const { result, error } = e.data || {};
                        if (error) reject(new Error(error));
                        else resolve({ byPos: new Map(result.byPos), reason: result.reason, hasPae: result.hasPae, n: result.n });
                    };
                    worker.onerror = () => { clearTimeout(timer); reject(new Error('worker error')); };
                    worker.postMessage({ geometry, sequence, amEntries: [...amMap], pae: pae ? { n: pae.n, data: pae.data } : null });
                });
            } catch (_) {
                _workerBroken = true; // don't retry worker creation this session; fall through to sync
            } finally {
                if (w) { try { w.terminate(); } catch (_) {} } // offload: no lingering worker thread
            }
        }
        // Synchronous fallback — paint the spinner before the brief main-thread block.
        await new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)));
        return UFVPocket.computePockets(geometry, amMap, sequence, pae);
    }

    // PTMs the user opted to co-display in the Disease & Variants view (off by default).
    function activeCoDisplayPtms() {
        const s = UFVState.state;
        if (!s.variantPtmCats || s.variantPtmCats.size === 0) return [];
        return s.ptms.filter(p => s.variantPtmCats.has(p.category));
    }

    function buildFilters() {
        buildPTMFilters();
        buildVariantFilters();
        buildSiteFilters();
    }

    // "Site" annotations — rendered as a collapsible group in BOTH the PTM and the variant
    // panels (the same s.sites list, with one shared per-site visibility flag). Off by default.
    function activeSites() {
        return UFVState.state.sites.filter(x => x.visible);
    }

    function buildSiteFilters() {
        const s = UFVState.state;
        const PANEL_PAIRS = [['ufv-sites-section-ptm', 'ufv-sites-ptm-list'], ['ufv-sites-section-var', 'ufv-sites-var-list']];
        PANEL_PAIRS.forEach(([secId, listId]) => {
            const section = byId(secId), list = byId(listId);
            if (!section || !list) return;
            list.textContent = '';
            if (!s.sites.length) { section.classList.add('ufv-hidden'); return; }
            section.classList.remove('ufv-hidden');
            s.sites.forEach((site, idx) => {
                // Label as "VAL 94" (three-letter AA + position), matching the residue panel title.
                const aa = AA1TO3[s.sequence?.[site.position - 1]] || '';
                const loc = site.endPosition && site.endPosition !== site.position
                    ? `${aa} ${site.position}–${site.endPosition}` : `${aa} ${site.position}`.trim();
                list.appendChild(makeFilterItem(`${loc}: ${site.description}`, site.color, '', !!site.visible, checked => {
                    site.visible = checked;
                    // Mirror the change in the other panel's checkbox for this site index.
                    PANEL_PAIRS.forEach(([, otherId]) => {
                        if (otherId === listId) return;
                        const cb = byId(otherId)?.querySelectorAll('input[type="checkbox"]')[idx];
                        if (cb) cb.checked = checked;
                    });
                    applySiteChange();
                }, undefined, { zoomPos: site.position }));
            });
        });
    }

    function sitesSetAll(select) {
        UFVState.state.sites.forEach(site => site.visible = select);
        buildSiteFilters(); // rebuild both lists so their checkboxes reflect the new state
        applySiteChange();
    }

    // Ligands present in the loaded model (AlphaFill cofactors etc.). A collapsible list in both
    // panels; clicking an entry focuses that ligand (zoom, nearby protein residues, others hidden)
    // and opens its chemistry detail panel.
    function ligKey(l) { return (l.chain == null ? '' : l.chain) + '|' + l.resi; }

    // All = show every non-ion ligand (ions respect excludeIons setting); None = hide everything.
    function ligandsSetAll(show) {
        const ions = StructureViewer.ION_CODES || new Set();
        StructureViewer.showLigands = true;
        if (show) {
            const newHidden = new Set();
            if (StructureViewer.excludeIons)
                UFVState.state.ligands.filter(l => ions.has(l.resn)).forEach(l => newHidden.add(ligKey(l)));
            StructureViewer.hiddenLigands = newHidden;
        } else {
            StructureViewer.hiddenLigands = new Set(UFVState.state.ligands.map(ligKey));
        }
        StructureViewer._focusNearbyLigands = null;   // unfocus so the set takes full effect
        UFVState.state.selectedLigand = null;
        byId('ufv-details')?.classList.remove('show');
        StructureViewer._drawLigands();
        refreshLigandSections();
        applyMode();
    }

    // Toggle visibility of a set of ligand copies in one frame update.
    function setLigandGroupVisible(copies, visible) {
        copies.forEach(l => { const k = ligKey(l); visible ? StructureViewer.hiddenLigands.delete(k) : StructureViewer.hiddenLigands.add(k); });
        StructureViewer._drawLigands();
    }

    // ---- AlphaFill identity-threshold filter -----------------------------------------------------
    // AlphaFill packs many overlapping transplanted ligands (the "blob"). Each carries a donor
    // sequence identity; filtering low-identity transplants declutters it, exactly like alphafill.eu.
    let _ligandThreshold = 0.25;            // active identity cutoff (0–1); 0.25 = AlphaFill default (show all)
    const LIGAND_THRESHOLD_STEPS = [25, 30, 40, 50, 60, 70]; // preset buttons (%)

    function ligandsHaveIdentity() { return UFVState.state.ligands.some(l => l.identity != null); }

    // Hide every transplant below the identity cutoff (non-transplant ligands untouched). NO clash
    // filtering — the user wants every AlphaFill transplant loaded; the clash score is shown in the
    // ligand panel purely as information. resetManual=true wipes per-ligand toggles.
    function applyLigandFilters(resetManual) {
        if (resetManual) {
            StructureViewer.hiddenLigands = new Set();
            // Restore ion exclusion that was active before the reset (e.g. AlphaFill async annotation).
            if (StructureViewer.excludeIons) {
                const ions = StructureViewer.ION_CODES || new Set();
                UFVState.state.ligands.filter(l => ions.has(l.resn)).forEach(l => StructureViewer.hiddenLigands.add(ligKey(l)));
            }
        }
        UFVState.state.ligands.forEach(l => {
            if (l.identity == null) return;
            const k = ligKey(l);
            if (l.identity < _ligandThreshold) StructureViewer.hiddenLigands.add(k);
            else StructureViewer.hiddenLigands.delete(k);
        });
        // Do NOT clear _focusNearbyLigands here: the async AlphaFill data fetch calls this AFTER the user
        // may have already focused a ligand, and clearing it would un-hide every transplant (the blob
        // returns, forcing a second click). While focused, _drawLigands keeps only the focused ligand.
        StructureViewer._drawLigands();
    }
    function applyLigandThreshold(threshold, resetManual) { _ligandThreshold = threshold; applyLigandFilters(resetManual); }

    async function annotateAlphaFillLigands(mySeq) {
        const s = UFVState.state;
        const st = UFVState.selectedStructure();
        if (!st || !/alphafill/i.test(st.provider || st.source || st.label || '')) return;
        const map = await UFVApi.getAlphaFillTransplants(s.uniprotId);
        if (_loadSeq !== mySeq || !map) return;            // a newer structure owns the viewer now
        let any = false;
        s.ligands.forEach(l => {
            const meta = map.get(l.chain == null ? '' : String(l.chain));
            if (meta) { l.identity = meta.identity; l.clash = meta.clash; l.donorPdb = meta.donorPdb; any = true; }
        });
        if (!any) return;
        // Identity defaults to 25% = AlphaFill's floor = show every transplant. No clash filtering.
        applyLigandFilters(true);
        refreshLigandSections();
    }

    // "Identity ≥" preset buttons (25/30/40/50/60/70%) prepended to the ligand list for AlphaFill models.
    function makeLigandThresholdRow() {
        const total = UFVState.state.ligands.filter(l => l.identity != null).length;
        const shown = UFVState.state.ligands.filter(l => l.identity != null && l.identity >= _ligandThreshold).length;
        const cur = Math.round(_ligandThreshold * 100);
        const box = document.createElement('div');
        box.className = 'ufv-lig-thresh';
        const lbl = document.createElement('span');
        lbl.className = 'ufv-lig-thresh-lbl';
        lbl.textContent = 'Identity ≥';
        box.appendChild(lbl);
        LIGAND_THRESHOLD_STEPS.forEach(pct => {
            const b = document.createElement('button');
            b.className = 'ufv-thresh-btn' + (pct === cur ? ' active' : '');
            b.textContent = pct;
            b.title = `Show transplants with ≥ ${pct}% donor identity`;
            b.addEventListener('click', () => { applyLigandThreshold(pct / 100, true); buildLigandFilters(); });
            box.appendChild(b);
        });
        const cnt = document.createElement('span');
        cnt.className = 'ufv-lig-thresh-cnt';
        cnt.textContent = `${shown}/${total}`;
        box.appendChild(cnt);
        return box;
    }

    // Rebuild ligand UI everywhere it can appear: the static PTM/Variant panels (buildLigandFilters)
    // and — since ligands are only enumerated after a structure loads — the dynamic feature windows.
    function refreshLigandSections() {
        buildLigandFilters();
        const cm = UFVState.state.currentMode;
        if (cm === 'sites' || cm === 'domains') buildFeatureWindow(cm);
    }

    function buildLigandFilters() {
        const s = UFVState.state;
        // Group copies of the same component (CCD) so multiple chains/copies nest under one entry.
        const groups = new Map();
        s.ligands.forEach(l => { if (!groups.has(l.resn)) groups.set(l.resn, []); groups.get(l.resn).push(l); });
        [['ufv-ligands-section-ptm', 'ufv-ligands-ptm-list'], ['ufv-ligands-section-var', 'ufv-ligands-var-list']].forEach(([secId, listId]) => {
            const section = byId(secId), list = byId(listId);
            if (!section || !list) return;
            list.textContent = '';
            if (!s.ligands.length) { section.classList.add('ufv-hidden'); return; }
            section.classList.remove('ufv-hidden');
            if (ligandsHaveIdentity()) list.appendChild(makeLigandThresholdRow()); // AlphaFill: identity filter
            groups.forEach((copies, resn) => list.appendChild(makeLigandGroup(resn, copies)));
        });
    }

    // Reusable dismissible info popover (M7). Returns an ⓘ icon; clicking it opens a small panel with a
    // title + HTML body (predictor/database explanations, each linking out). The × or a click-away closes
    // it. Stops propagation so it can sit inside a clickable header without toggling it.
    function makeInfoIcon(title, bodyHtml) {
        const icon = document.createElement('span');
        icon.className = 'ufv-info-icon';
        icon.textContent = 'ⓘ';
        icon.setAttribute('aria-label', title);
        icon.addEventListener('click', e => {
            e.stopPropagation(); e.preventDefault();
            const open = document.querySelector('.ufv-info-pop');
            const sameOwner = open && open.dataset.owner === title;
            if (open) open.remove();
            if (sameOwner) return; // clicking the same (i) toggles it off
            const pop = document.createElement('div');
            pop.className = 'ufv-info-pop'; pop.dataset.owner = title;
            pop.innerHTML = `<button class="ufv-info-close" aria-label="Close">×</button><div class="ufv-info-title">${title}</div><div class="ufv-info-body">${bodyHtml}</div>`;
            document.body.appendChild(pop);
            const rect = icon.getBoundingClientRect();
            pop.style.left = Math.max(8, Math.min(rect.left, window.innerWidth - pop.offsetWidth - 8)) + 'px';
            pop.style.top = Math.min(rect.bottom + 6, window.innerHeight - pop.offsetHeight - 8) + 'px';
            // Close on ANY interaction outside the popover or its icon. Capture-phase mousedown fires before
            // (and regardless of) stopPropagation from menus/sections it may be nested in.
            const close = () => { pop.remove(); document.removeEventListener('mousedown', onAway, true); window.removeEventListener('blur', close); };
            const onAway = (ev) => { if (!pop.contains(ev.target) && !icon.contains(ev.target)) close(); };
            pop.querySelector('.ufv-info-close').addEventListener('click', close);
            setTimeout(() => { document.addEventListener('mousedown', onAway, true); window.addEventListener('blur', close); }, 0);
        });
        return icon;
    }

    // Predictor / database explanations shown in the Predictions info popover. Each links to its source.
    const PREDICTOR_INFO = [
        ['AlphaMissense', 'Pathogenicity of a missense substitution (0–1; &gt;0.564 likely pathogenic, &lt;0.34 likely benign).', 'https://alphamissense.hegelab.org'],
        ['EVE', 'Evolutionary model of Variant Effect — unsupervised, from sequence alignments (0–1, higher = more pathogenic).', 'https://evemodel.org'],
        ['CADD', 'Combined Annotation-Dependent Depletion — a phred-scaled deleteriousness score (&gt;20 ≈ top 1% most deleterious genome-wide, &gt;30 ≈ top 0.1%). Genomic, so it covers only substitutions reachable by a single nucleotide change.', 'https://cadd.gs.washington.edu'],
        ['ESM-1b', 'Protein language-model pathogenicity score for the substitution; lower (more negative) = more pathogenic. Scores every missense substitution.', 'https://github.com/facebookresearch/esm'],
        ['FoldX ΔΔG', 'Predicted change in folding free energy (kcal/mol; &gt;2 destabilising, &lt;-1 stabilising, ~0 neutral).', 'https://foldxsuite.crg.eu'],
        ['M3D (Missense3D)', 'Predicted structural consequence of the substitution (damaging vs neutral).', 'http://missense3d.bc.ic.ac.uk'],
        ['Binding pocket', 'Predicted ligand-binding pocket(s) lined by this residue — with buriedness, pocket score and the lining residues (AutoSite, via ProtVar).', 'https://www.ebi.ac.uk/ProtVar'],
        ['Conservation & ProtVar', 'Per-residue evolutionary conservation, and the EVE/ESM1b/FoldX/M3D/pocket values, are served by ProtVar (EMBL-EBI).', 'https://www.ebi.ac.uk/ProtVar'],
    ].map(([name, desc, url]) => `<p><b>${name}</b> — ${desc} <a href="${url}" target="_blank" rel="noopener noreferrer">source ↗</a></p>`).join('');

    // Pocket-metric explanations (info popover on the Pocket section). These are COMPUTATIONAL predictions of
    // cavities on the protein, not the ligands shown and not a specific highlightable region.
    const POCKET_INFO = [
        ['Predicted pocket', 'A predicted cavity that could bind a small molecule, identified with AutoSite on the AlphaFold model and served by ProtVar. Click a pocket’s “res” box to highlight it and list its residues.', 'https://www.ebi.ac.uk/ProtVar'],
        ['Buriedness', 'How enclosed the pocket is (0 = open, 1 = buried).', ''],
        ['Pocket score', 'AutoSite’s ranking of how pronounced and ligand-favourable the cavity is — higher = a larger, better-defined pocket.', ''],
    ].map(([name, desc, url]) => `<p><b>${name}</b> — ${desc}${url ? ` <a href="${url}" target="_blank" rel="noopener noreferrer">source ↗</a>` : ''}</p>`).join('');

    // Small magnifier "zoom to" button used across sidebar lists (ligands, and per-residue entries).
    function makeZoomBtn(title, onClick) {
        const b = document.createElement('button');
        b.className = 'ufv-zoom-btn';
        b.title = title || 'Zoom in';
        b.innerHTML = '<svg viewBox="0 0 16 16" width="12" height="12" aria-hidden="true"><circle cx="6.5" cy="6.5" r="4.2" fill="none" stroke="currentColor" stroke-width="1.6"/><line x1="9.7" y1="9.7" x2="14" y2="14" stroke="currentColor" stroke-width="1.7" stroke-linecap="round"/></svg>';
        b.addEventListener('click', e => { e.stopPropagation(); onClick(); });
        return b;
    }

    // One CCD entry. Single copy → a row with a visibility checkbox, name, location, and a zoom icon.
    // Multiple copies → a header (checkbox toggles all, chevron expands) over per-copy rows, each
    // with its own checkbox + zoom icon.
    function makeLigandGroup(resn, copies) {
        const s = UFVState.state;
        const isSel = l => s.selectedLigand && s.selectedLigand.resn === l.resn && s.selectedLigand.resi === l.resi && s.selectedLigand.chain === l.chain;
        const visible = l => StructureViewer.isLigandVisible(l.chain, l.resi);
        const locText = l => `${l.chain ? l.chain + ' ' : ''}${l.resi}`;
        const focus = l => onLigandClick({ resn: l.resn, resi: l.resi, chain: l.chain });

        const wrap = document.createElement('div');
        wrap.className = 'ufv-lig-group';
        const hdr = document.createElement('div');
        hdr.className = 'ufv-lig-row';
        const cb = document.createElement('input');
        cb.type = 'checkbox'; cb.className = 'ufv-lig-eye';
        cb.checked = copies.some(visible); cb.title = 'Show / hide';
        cb.addEventListener('change', () => { setLigandGroupVisible(copies, cb.checked); buildLigandFilters(); });
        const name = document.createElement('span');
        name.className = 'ufv-lig-ccd';
        name.textContent = resn;
        hdr.append(cb, name);

        // Right-aligned controls: optional AlphaFill identity badge + the zoom icon.
        const idBadge = l => {
            if (l.identity == null) return null;
            const b = document.createElement('span');
            b.className = 'ufv-lig-id';
            b.textContent = `${Math.round(l.identity * 100)}%` + (l.donorPdb ? ` · ${l.donorPdb}` : '');
            b.title = `Donor identity` + (l.donorPdb ? ` from ${l.donorPdb}` : '');
            return b;
        };
        const appendRight = (rowEl, l, zb) => {
            const idb = idBadge(l);
            if (idb) { idb.style.marginLeft = 'auto'; rowEl.append(idb, zb); }
            else { zb.style.marginLeft = 'auto'; rowEl.append(zb); }
        };

        if (copies.length === 1) {
            const l = copies[0];
            // Chain/residue inline right after the name; identity badge + zoom icon stay far right.
            const loc = document.createElement('span');
            loc.className = 'ufv-lig-loc';
            loc.textContent = locText(l);
            if (isSel(l)) hdr.classList.add('selected');
            name.classList.add('ufv-lig-clickable');
            name.addEventListener('click', () => focus(l));
            hdr.append(loc);
            appendRight(hdr, l, makeZoomBtn(`Zoom to ${resn}`, () => focus(l)));
            wrap.appendChild(hdr);
            return wrap;
        }

        // Multiple copies: chevron-expandable nested list.
        const meta = document.createElement('span');
        meta.className = 'ufv-lig-meta';
        meta.textContent = `×${copies.length}`;
        hdr.appendChild(meta);
        const chev = document.createElement('span');
        chev.className = 'ufv-lig-chevron';
        chev.innerHTML = '&#9654;';
        const sub = document.createElement('div');
        sub.className = 'ufv-lig-copies ufv-collapsed';
        const toggleExpand = () => { const c = sub.classList.toggle('ufv-collapsed'); chev.innerHTML = c ? '&#9654;' : '&#9660;'; };
        name.classList.add('ufv-lig-clickable');
        name.addEventListener('click', toggleExpand);
        chev.addEventListener('click', toggleExpand);
        hdr.appendChild(chev);

        copies.forEach(l => {
            const row = document.createElement('div');
            row.className = 'ufv-lig-row ufv-lig-copy' + (isSel(l) ? ' selected' : '');
            const ccb = document.createElement('input');
            ccb.type = 'checkbox'; ccb.className = 'ufv-lig-eye';
            ccb.checked = visible(l); ccb.title = 'Show / hide';
            ccb.addEventListener('change', () => { StructureViewer.setLigandVisible(l.chain, l.resi, ccb.checked); cb.checked = copies.some(visible); });
            const loc = document.createElement('span');
            loc.className = 'ufv-lig-ccd ufv-lig-clickable';
            loc.textContent = locText(l);
            loc.addEventListener('click', () => focus(l));
            row.append(ccb, loc);
            appendRight(row, l, makeZoomBtn(`Zoom to ${resn} ${locText(l)}`, () => focus(l)));
            sub.appendChild(row);
        });
        wrap.append(hdr, sub);
        return wrap;
    }

    function applySiteChange() {
        if (UFVState.state.currentMode === 'ptm') applyPTMMode();
        else applyMode();
    }

    function buildPTMFilters() {
        const list = byId('ufv-ptm-list');
        list.textContent = '';
        Object.entries(UFVState.state.ptmGroups).sort((a, b) => b[1].items.length - a[1].items.length).forEach(([cat, group]) => {
            list.appendChild(makeExpandableFilter(cat, group.color, group.items, group.visible, checked => {
                group.visible = checked;
                group.items.forEach(item => item.visible = checked);
                applyPTMMode();
            }, item => `Residue ${item.position}: ${item.description}`));
        });
        makePanelHeaderCollapsible('ufv-ptm-panel', list);
    }

    // Make a static .ufv-panel-hdr collapsible (chevron + click-to-toggle its list), so a window's
    // primary category collapses like every other section. Idempotent — runs on each rebuild.
    function makePanelHeaderCollapsible(panelId, listEl) {
        const hdr = byId(panelId)?.querySelector('.ufv-panel-hdr');
        if (!hdr || hdr.dataset.ufvCollapsible) return;
        hdr.dataset.ufvCollapsible = '1';
        const chev = document.createElement('span');
        chev.className = 'ufv-collapsible-chevron';
        chev.innerHTML = '&#9660;';
        hdr.insertBefore(chev, hdr.firstChild);
        hdr.addEventListener('click', e => {
            if (e.target.closest('button')) return; // don't toggle when clicking All/None/C
            const collapsed = listEl.classList.toggle('ufv-hidden');
            chev.innerHTML = collapsed ? '&#9654;' : '&#9660;';
        });
    }

    // PTM co-display list inside the Disease & Variants panel — collapsed, all-off by default.
    // Lets the user overlay PTM spheres on top of the variant view without leaving it.
    function buildVariantPtmFilters() {
        const s = UFVState.state;
        const section = byId('ufv-vptm-section');
        const list = byId('ufv-vptm-list');
        if (!section || !list) return;
        list.textContent = '';
        const groups = Object.entries(s.ptmGroups).sort((a, b) => b[1].items.length - a[1].items.length);
        if (groups.length === 0) { section.classList.add('ufv-hidden'); return; }
        section.classList.remove('ufv-hidden');
        groups.forEach(([cat, group]) => {
            list.appendChild(makeFilterItem(cat, group.color, group.items.length, s.variantPtmCats.has(cat), checked => {
                checked ? s.variantPtmCats.add(cat) : s.variantPtmCats.delete(cat);
                applyMode();
            }));
        });
    }

    function variantPtmSetAll(select) {
        const s = UFVState.state;
        s.variantPtmCats = new Set(select ? Object.keys(s.ptmGroups) : []);
        document.querySelectorAll('#ufv-vptm-list input[type="checkbox"]').forEach(cb => { cb.checked = select; });
        applyMode();
    }

    // The disease names to offer (and their summary metadata). On the entry page this is
    // restricted to the diseases that appear as h4 headings in the page's Disease & Variants
    // section; on the variant-viewer page (or when nothing was scraped) it's every disease.
    function diseaseNamesToShow() {
        const s = UFVState.state;
        const ds = DataProcessor.getDiseaseSummary(s.variants);
        // 'Unclassified' = variants with no disease; selecting it never produces spheres (those variants
        // have empty .diseases, so the disease filter never matches them). Drop the dead toggle entirely.
        let names = Object.keys(ds).filter(n => n !== 'Unclassified');
        if (s.pageContext !== 'variant-viewer') {
            const scraped = s.scrapedDiseases || [];
            if (scraped.length > 0) {
                const scrapedIds = new Set(scraped.map(d => d.id).filter(Boolean));
                const scrapedAbbrs = new Set();
                scraped.forEach(d => { const m = (d.label || '').match(/\(([A-Z][A-Z0-9]+)\)\s*$/); if (m) scrapedAbbrs.add(m[1]); });
                // 1-to-1 label → disease ID (avoids a variant's labels all inheriting all its IDs).
                const labelToId = new Map();
                s.variants.forEach(v => (v.diseasePairs || []).forEach(({ id, label }) => { if (id && !labelToId.has(label)) labelToId.set(label, id); }));
                const filtered = names.filter(name => {
                    if (name === 'Unclassified') return false;
                    if (scrapedAbbrs.has(name)) return true;
                    const id = labelToId.get(name);
                    return !!(id && scrapedIds.has(id));
                });
                if (filtered.length > 0) names = filtered;
            }
        }
        return { names, ds };
    }

    function buildVariantFilters() {
        const s = UFVState.state;
        fillFilterList('ufv-prov-list', DataProcessor.getProvenanceSummary(s.variants), s.activeProvenances, 'provenance');
        fillFilterList('ufv-cons-list', DataProcessor.getConsequenceSummary(s.variants), s.activeConsequences, 'consequence');
        buildVariantPtmFilters();
        const { names: diseasesToShow, ds } = diseaseNamesToShow();
        const dis = byId('ufv-dis-section');
        const list = byId('ufv-dis-list');
        list.textContent = '';
        if (diseasesToShow.length) {
            dis.classList.remove('ufv-hidden');
            s.activeDiseases = new Set(diseasesToShow);
            diseasesToShow.forEach(name => {
                const meta = ds[name] || { color: '#9e9e9e', count: 0 };
                list.appendChild(makeDiseaseFilter(name, meta.color, meta.count));
            });
        } else {
            dis.classList.add('ufv-hidden');
            s.activeDiseases = null;
        }
        appendDomainSection('ufv-var-panel'); // Family & Domains layer in the Variant window too
    }

    function fillFilterList(id, summary, activeSet, dim) {
        const list = byId(id);
        list.textContent = '';
        Object.entries(summary).forEach(([label, meta]) => {
            activeSet.add(label);
            list.appendChild(makeFilterItem(label, meta.color, meta.count, true,
                checked => toggleVariantFilter(dim, label, checked), label));
        });
    }

    function renderSequence() {
        const s = UFVState.state;
        const wrap = byId('ufv-sequence-wrap');
        wrap.textContent = '';
        if (!s.sequence) return;
        const st = UFVState.selectedStructure();
        const isAlphaFold = st?.source === 'AlphaFold';
        // Multi-chain experimental structures get one ribbon track per chain: each subunit
        // resolves a (slightly) different set of residues and carries its own structure-
        // dependent hotspot / contact-hub tiers, so they must be shown separately.
        if (st?.chainIds?.length > 1) {
            wrap.classList.add('ufv-multichain');
            st.chainIds.forEach(chain => {
                wrap.appendChild(buildSequenceTrack(chain, StructureViewer.mappedResiduesForChain(chain), isAlphaFold, `Chain ${chain}`));
            });
        } else {
            wrap.classList.remove('ufv-multichain');
            wrap.appendChild(buildSequenceTrack(null, StructureViewer.mappedResidues(), isAlphaFold, null));
        }
    }

    /** Build one sequence-ribbon track (a chain, or the whole structure when chain is null). */
    function buildSequenceTrack(chain, mappedResi, isAlphaFold, label) {
        const s = UFVState.state;
        const hasCoverage = mappedResi !== null;
        const coverage = new Set(mappedResi || s.sequence.split('').map((_, i) => i + 1));
        const ptmPos = new Set(activePtms().map(p => p.position));
        const varPos = new Set(s.variants.map(v => v.position));
        const displayed = new Set(s.displayedPositions);
        // Per-chain structure-dependent tiers reflected directly in the ribbon.
        const mode = getColorMode();
        const tierColors = mode === 'hotspots' ? { strong: '#b71c1c', moderate: '#e64a19', weak: '#ffa726' }
                         : mode === 'distantContacts' ? { strong: '#6a1b9a', moderate: '#ab47bc' } : null;
        let tierMap = null;
        if (mode === 'hotspots') tierMap = chain != null ? s.analysis.hotspotsByChain?.get(chain) : s.analysis.hotspots;
        else if (mode === 'distantContacts') tierMap = chain != null ? s.analysis.distantContactsByChain?.get(chain) : s.analysis.distantContacts;

        const track = document.createElement('div');
        track.className = 'ufv-seq-track';
        if (label) {
            const lab = document.createElement('span');
            lab.className = 'ufv-seq-chain-label';
            lab.textContent = label;
            track.appendChild(lab);
        }
        const row = document.createElement('div');
        row.className = 'ufv-seq-aas';
        s.sequence.split('').forEach((aa, idx) => {
            const pos = idx + 1;
            const span = document.createElement('button');
            span.className = 'ufv-seq-aa';
            span.textContent = aa;
            span.title = `${aa}${pos}`;
            const missing = hasCoverage && !isAlphaFold && !coverage.has(pos);
            if (missing) {
                span.classList.add('missing');
                span.disabled = true;
                span.title = `${aa}${pos} — not resolved in this ${chain != null ? 'chain' : 'structure'}`;
            } else if (coverage.has(pos)) {
                span.classList.add('covered');
            }
            if (ptmPos.has(pos)) span.classList.add('ptm');
            if (varPos.has(pos)) span.classList.add('variant');
            if (displayed.has(pos)) span.classList.add('visible');
            if (s.nearbyResidues.has(pos)) span.classList.add('nearby');
            if (s.selectedResidue === pos) span.classList.add('selected');
            if (!missing && tierMap && tierColors) {
                const tier = tierMap.get(pos);
                if (tier) { span.classList.add('feature'); span.style.boxShadow = `inset 0 -3px 0 ${tierColors[tier]}`; }
            }
            if (!missing) span.addEventListener('click', () => onClick({ position: pos, variants: s.variants.filter(v => v.position === pos) }, 'sequence', chain));
            row.appendChild(span);
        });
        track.appendChild(row);
        return track;
    }

    function renderLegend(mode) {
        const legend = byId('ufv-legend');
        legend.textContent = '';
        if (mode === 'bfactor') {
            const item = document.createElement('span');
            item.className = 'ufv-legend-item';
            const bar = document.createElement('span');
            bar.className = 'ufv-legend-gradient';
            bar.style.background = 'linear-gradient(to right, #313695, #f7f7f7, #d73027)';
            bar.title = 'B-factor: low (rigid, blue) → high (flexible, red)';
            item.append(bar, 'B-factor');
            legend.appendChild(item);
            return;
        }
        const modeItems = {
            plddt: [['#ff7d45', '<50'], ['#ffdb13', '50–70'], ['#65cbf3', '70–90'], ['#0053d6', '>90 pLDDT']],
            hotspots: [['#b9c2cf', 'Not enriched'], ['#ffa726', 'Weak'], ['#e64a19', 'Moderate'], ['#b71c1c', 'Strong hotspot']],
            distantContacts: [['#b9c2cf', 'Low centrality'], ['#ab47bc', 'Moderate'], ['#6a1b9a', 'High centrality']],
            alphaMissense: [['#3d85c8', 'Likely benign (<0.34)'], ['#b9c2cf', 'Ambiguous (0.34–0.564)'], ['#e06666', 'Likely pathogenic (0.564–0.78)'], ['#b71c1c', 'Pathogenic (>0.78)']],
            residueBurden: [['#b9c2cf', 'Not recurrent'], ['#e65100', 'Recurrent phenotype']],
            prism: [['#00897b', 'Buried cluster'], ['#8e24aa', 'Exposed site'], ['#b9c2cf', 'Not significant']],
        };

        // Safe DOM helper — never uses innerHTML so API-derived labels can't inject HTML.
        const appendLegendItem = (color, label) => {
            const row = document.createElement('span');
            row.className = 'ufv-legend-item';
            const swatch = document.createElement('span');
            swatch.style.background = color;
            row.append(swatch, label);
            legend.appendChild(row);
        };

        // Topology legend is built from the entry's actual segments (one swatch per colour).
        if (mode === 'topology') {
            const seen = new Map();
            (UFVState.state.topology || []).forEach(seg => { if (!seen.has(seg.color)) seen.set(seg.color, seg.label); });
            seen.forEach((label, color) => appendLegendItem(color, label));
        }

        // Only show items for explicit coloring modes; default cyan needs no legend unless zoomed in
        (modeItems[mode] || []).forEach(([color, label]) => appendLegendItem(color, label));

        // Be explicit about which hotspot null was used (case-control vs spatial), since the
        // spatial fallback doesn't control for ascertainment bias.
        if (mode === 'hotspots') {
            const m = UFVState.state.analysis.hotspotMethod;
            if (m) {
                const note = document.createElement('span');
                note.className = 'ufv-legend-item ufv-legend-note';
                note.textContent = m === 'spatial' ? 'spatial clustering' : 'case–control enrichment';
                note.title = m === 'spatial'
                    ? 'Pathogenic clustering vs. random placement (no benign controls available); does not adjust for which residues were studied.'
                    : 'Pathogenic vs. benign spatial enrichment (case–control permutation).';
                legend.appendChild(note);
            }
        }

        // Constraint-pocket mode: surface a short status note (candidate count / why empty).
        if (mode === 'prism') {
            const res = UFVState.state.analysis.prism;
            const note = document.createElement('span');
            note.className = 'ufv-legend-item ufv-legend-note';
            if (res?.reason) {
                note.textContent = res.reason;
                note.title = res.reason;
            } else if (res) {
                const n = filteredPocketByPos()?.size || 0;
                note.textContent = `${n} candidate cluster${n === 1 ? '' : 's'} (q≤${sensThreshold.toFixed(2)})`;
                note.title = 'Exploratory heuristic (not a validated predictor): buried, evolutionarily-constrained candidate sites — Getis-Ord Gi* on AlphaMissense residuals vs structural burial.';
            } else {
                note.textContent = 'constraint clusters';
            }
            legend.appendChild(note);
        }

        // In focus/zoom mode, describe what the colours PAINTED ON THE POCKET STICKS mean. The sticks are
        // coloured per-residue by buildAnnotationMap, which (by request) reflects DISEASE-ASSOCIATED VARIANT
        // pathogenicity ONLY — not PTMs or functional sites. Shown for residue AND ligand focus.
        const s = UFVState.state;
        if (s.selectedResidue !== null || s.selectedLigand) {
            // Build the legend ONLY from the variant consequences ACTUALLY painted on the pocket sticks
            // (nearbyResidues) — matching buildAnnotationMap (variants only).
            const pocket = s.nearbyResidues || new Set();
            const varByPos = new Map();
            s.variants.forEach(v => varByPos.set(v.position, v)); // last wins, matching buildAnnotationMap
            const present = []; const seen = new Set();
            const add = (label, color) => { if (label && !seen.has(label)) { seen.add(label); present.push([color, label]); } };
            pocket.forEach(pos => {
                const v = varByPos.get(pos); if (v) add(v.consequence, v.consequenceColor);
            });
            if (present.length) {
                const brk = document.createElement('span'); brk.className = 'ufv-legend-break'; legend.appendChild(brk); // "side chains:" starts its own line below the mode legend
                const note = document.createElement('span');
                note.className = 'ufv-legend-item ufv-legend-note';
                note.textContent = 'side chains:';
                note.title = 'Colours painted on the binding-pocket / neighbour side chains';
                legend.appendChild(note);
                present.forEach(([color, label]) => appendLegendItem(color, label));
            }
        }
    }

    function onHover(data, mode, event) {
        const tip = byId('ufv-tooltip');
        if (!data) {
            tip.classList.remove('show');
            return;
        }
        const pos = data.position || data.resi;
        const aa = UFVState.state.sequence?.[pos - 1] || '';
        byId('ufv-tooltip-hdr').textContent = `${aa}${pos}`;
        byId('ufv-tooltip-body').innerHTML = residueSummary(pos, data);
        if (event) {
            const rect = document.querySelector('.ufv-viewer-wrap').getBoundingClientRect();
            tip.style.left = `${Math.min((event.clientX || 0) - rect.left + 14, rect.width - 280)}px`;
            tip.style.top = `${Math.min((event.clientY || 0) - rect.top + 14, rect.height - 70)}px`;
        }
        tip.classList.add('show');
    }

    // One-letter → three-letter amino acid code (all caps) for the residue title.
    const AA1TO3 = { A:'ALA',C:'CYS',D:'ASP',E:'GLU',F:'PHE',G:'GLY',H:'HIS',I:'ILE',K:'LYS',L:'LEU',M:'MET',N:'ASN',P:'PRO',Q:'GLN',R:'ARG',S:'SER',T:'THR',V:'VAL',W:'TRP',Y:'TYR' };

    // Session-persistent toggle for sphere visibility in focus mode.
    let _showOtherSpheres = true;
    let _showVarEvidence = false;   // persistent: show Review/dbSNP/gnomAD/Genomic in the variant blocks
    let _colorProfile = 'consequence'; // radio-button color prioritization
    let _domainCartoon = false; // C on Family & Domains → colour the backbone by domain in any window
    let _diseaseColorMap = new Map();  // disease name → DISEASE_PALETTE color (from computeDiseaseColors)
    let _hiddenVariantKeys = new Set(); // individually-hidden variants ('position|mutant') from disease dropdown
    let _forcedVariantKeys = new Set(); // individually-SHOWN variants — render even if their disease axis is off
    const variantKey = v => `${v.position}|${v.mutant || ''}`;

    // Helper: make a small CSS toggle switch (<label> wrapping hidden <input> + slider span).
    function makeToggle(checked, title) {
        const lbl = document.createElement('label');
        lbl.className = 'ufv-toggle-switch';
        lbl.title = title || '';
        const chk = document.createElement('input');
        chk.type = 'checkbox';
        chk.checked = checked;
        const sl = document.createElement('span');
        sl.className = 'ufv-toggle-slider';
        lbl.append(chk, sl);
        // Prevent toggle click from bubbling to parent button (collapsible header).
        lbl.addEventListener('click', e => e.stopPropagation());
        return { lbl, chk };
    }

    function onClick(data, _mode, chain = null) {
        const s = UFVState.state;
        const pos = Number(data.position);
        s.selectedResidue = pos;
        s.selectedChain = chain;
        s.selectedLigand = null; // clicking a residue clears any ligand focus
        const annotations = buildAnnotationMap();
        s.nearbyResidues = StructureViewer.focusResidue(pos, chain, { annotatedResidues: annotations }, { showOtherSpheres: _showOtherSpheres }) || new Set([pos]);

        // Lazy-compute PTM–variant proximity (used by CSV export). Per-residue panel uses computeResidueProximity.
        if (!s.analysis.ptmVariantProximity && s.ptms.length && s.variants.length && StructureViewer.viewer) {
            const ptmR = s.settings.proxPtmRadius || 8, varR = s.settings.proxVarRadius || 12;
            s.analysis.ptmVariantProximity = UFVAnalysis.computePtmVariantProximity(s.ptms, s.variants, StructureViewer.residueGeometry(), ptmR, varR);
        }

        const body = byId('ufv-details-body');
        body.textContent = '';
        byId('ufv-sphere-toggle')?.classList.remove('ufv-hidden'); // restore (a partner panel may have hidden it)

        // ── Title: bulbs + "ALA 421" ────────────────────────────────────────────
        const wt = s.sequence?.[pos - 1] || '';
        const titleEl = byId('ufv-details-title');
        titleEl.textContent = '';
        const bulbDefs = [
            { label: 'Pathogenic variant hotspot', color: '#e53935', active: () => s.analysis.hotspots instanceof Map && s.analysis.hotspots.has(pos) },
            { label: 'Recurrent phenotype residue', color: '#e65100', active: () => s.analysis.residueBurden instanceof Set && s.analysis.residueBurden.has(pos) },
            { label: 'Contact-network centrality',  color: '#6a1b9a', active: () => s.analysis.distantContacts instanceof Map && s.analysis.distantContacts.has(pos) },
            { label: 'Burial-adjusted constraint cluster', color: '#00897b', active: () => s.analysis.prism?.byPos instanceof Map && s.analysis.prism.byPos.has(pos) },
        ];
        const bulbRow = document.createElement('span');
        bulbRow.className = 'ufv-bulb-row';
        const bulbInfo = bulbDefs.map(b => `<p><span style="display:inline-block;width:9px;height:9px;border-radius:50%;background:${b.color};margin-right:6px;vertical-align:middle"></span><b>${b.label}</b> — ${b.active() ? 'flagged at this residue' : 'computed; not flagged here'}.</p>`).join('')
            + '<p style="opacity:.75">A filled dot = this residue is flagged by that structural-genomics algorithm; faint = not flagged.</p>';
        bulbDefs.forEach(b => {
            const bulb = document.createElement('span');
            bulb.className = 'ufv-bulb' + (b.active() ? ' ufv-bulb-on' : '');
            bulb.style.setProperty('--bulb-color', b.color);
            bulb.title = b.label;
            bulbRow.appendChild(bulb);
        });
        titleEl.appendChild(document.createTextNode((AA1TO3[wt] || wt) + ' ' + pos));
        titleEl.appendChild(bulbRow);
        bulbRow.appendChild(makeInfoIcon('Structural-genomics flags', bulbInfo)); // explains the coloured dots

        // ── Nearby box (with live distance slider) ──────────────────────────────
        body.appendChild(makeNearbyBox(s.selectedChain));

        // ── PTM annotations ─────────────────────────────────────────────────────
        const ptmsAtPos = s.ptms.filter(p => p.position === pos || p.endPosition === pos);
        ptmsAtPos.forEach(p => body.appendChild(row('PTM', `${p.category}: ${p.description}`, p.color)));

        // ── Site annotations ────────────────────────────────────────────────────
        s.sites.filter(x => x.position === pos || x.endPosition === pos)
            .forEach(x => body.appendChild(row('Site', x.description, x.color)));

        // ── Mutagenesis (experimental) ──────────────────────────────────────────
        (s.mutagenesis || []).filter(mtg => mtg.position <= pos && pos <= mtg.endPosition).forEach(mtg => {
            // No inline colour — #6d4c41 (dark brown) was unreadable in dark mode. Use the theme-adaptive
            // .ufv-detail-val colour; the "Mutagenesis" label already identifies the row.
            const r = row('Mutagenesis', `${mtg.wildType}→${(mtg.mutants || []).join(', ')}: ${mtg.effect}`);
            r.title = 'Experimentally mutated residue (UniProt)';
            body.appendChild(r);
        });

        // ── Variants — collapsible ───────────────────────────────────────────────
        const variants = s.variants.filter(v => v.position === pos).slice(0, 12);
        if (variants.length > 0) {
            if (ptmsAtPos.length > 0) {
                const sep = document.createElement('div'); sep.className = 'ufv-variant-divider'; body.appendChild(sep);
            }
            // Determine top-severity count and color for the collapsed header badge.
            const sevOrder = ['Likely pathogenic or pathogenic', 'Predicted deleterious', 'Uncertain significance', 'Likely benign or benign'];
            const topConsequence = sevOrder.find(sev => variants.some(v => v.consequence === sev));
            const topColor = topConsequence ? (variants.find(v => v.consequence === topConsequence)?.consequenceColor || '#b9c2cf') : '#b9c2cf';
            const topCount = topConsequence ? variants.filter(v => v.consequence === topConsequence).length : 0;
            const countLabel = topCount > 0
                ? `<span class="ufv-am-ratio"><span style="color:${topColor}">${topCount}</span>/${variants.length}</span>`
                : `<span class="ufv-am-ratio">${variants.length}</span>`;
            const varSection = document.createElement('div');
            varSection.className = 'ufv-am-section';
            const varToggle = document.createElement('button');
            varToggle.className = 'ufv-am-toggle';
            varToggle.innerHTML = `<span class="ufv-am-hdr-left">Variants</span>`;
            const varBody = document.createElement('div');
            varBody.className = 'ufv-am-body';
            variants.forEach((v, i) => {
                if (i > 0) { const d = document.createElement('div'); d.className = 'ufv-variant-divider'; varBody.appendChild(d); }
                const vblock = document.createElement('div');
                vblock.className = 'ufv-vblock';
                const vtag = document.createElement('div');
                vtag.className = 'ufv-vtag';
                vtag.textContent = `${v.wildType}${v.position}${v.mutant}`;
                vtag.style.color = v.consequenceColor;
                vblock.appendChild(vtag);
                const sig = v.clinVarSignificance || v.consequence || '—';
                const dis = (v.diseases && v.diseases.length) ? ` (${v.diseases.join('; ')})` : '';
                vblock.appendChild(row('ClinVar', sig + dis));
                // Evidence rows (Review / dbSNP / gnomAD AF / Genomic) — hidden until the evidence
                // toggle is on. Marked with a class so the toggle can flip them without a re-render.
                const evi = (label, value, color) => { const r = row(label, value, color); r.classList.add('ufv-var-evidence'); if (!_showVarEvidence) r.style.display = 'none'; vblock.appendChild(r); };
                if (v.clinVarReviewStatus) evi('Review', v.clinVarReviewStatus);
                if (v.rsIds?.length) evi('dbSNP', v.rsIds.join(', '));
                if (Number.isFinite(v.gnomadAf)) {
                    const af = v.gnomadAf;
                    const txt = af === 0 ? '0 (not observed)' : af >= 0.01 ? `${(af * 100).toFixed(2)}%`
                        : af >= 1e-4 ? `${(af * 100).toFixed(4)}%` : af.toExponential(1);
                    const rarity = af >= 0.01 ? 'common' : af > 0 ? 'rare' : 'absent';
                    evi('gnomAD AF', `${txt} (${rarity})`, af >= 0.01 ? '#43a047' : null);
                } else {
                    evi('gnomAD AF', 'not in gnomAD');
                }
                if (v.genomicLocation) evi('Genomic', v.genomicLocation);
                varBody.appendChild(vblock);
            });
            // Header order: [Variants] [details toggle] [count] [arrow] — the evidence toggle sits BEFORE the
            // count. The label and the count+arrow group both collapse; the toggle is independent.
            const varRight = document.createElement('button');
            varRight.className = 'ufv-am-toggle ufv-am-right-toggle';
            varRight.innerHTML = `<span class="ufv-am-hdr-right">${countLabel}<span class="ufv-am-arrow">▾</span></span>`;
            const collapseVar = () => {
                const open = varBody.classList.toggle('show');
                const arrow = varRight.querySelector('.ufv-am-arrow'); if (arrow) arrow.textContent = open ? '▴' : '▾';
                eviWrap.style.display = open ? '' : 'none'; // the "details" toggle only appears once expanded
            };
            varToggle.addEventListener('click', collapseVar);
            varRight.addEventListener('click', collapseVar);
            // Evidence details: a TOGGLE SWITCH on the variant header. One global switch for
            // Review · dbSNP · gnomAD AF · Genomic across all variant blocks; flips _showVarEvidence and
            // updates every evidence row + every header switch in the panel at once.
            const flipEvidence = (on) => {
                _showVarEvidence = on;
                const body = byId('ufv-details-body');
                body.querySelectorAll('.ufv-var-evidence').forEach(r => { r.style.display = on ? '' : 'none'; });
                body.querySelectorAll('.ufv-evi-chk').forEach(c => { c.checked = on; });
            };
            const eviWrap = document.createElement('div');
            eviWrap.className = 'ufv-evi-switch';
            eviWrap.style.display = 'none'; // variants section starts collapsed → hide the details toggle until expanded
            eviWrap.title = 'Show evidence: review · dbSNP · gnomAD AF · genomic';
            eviWrap.innerHTML = `<span class="ufv-evi-lbl">details</span><label class="ufv-toggle-switch"><input type="checkbox" class="ufv-evi-chk"${_showVarEvidence ? ' checked' : ''}><span class="ufv-toggle-slider"></span></label>`;
            eviWrap.querySelector('input').addEventListener('change', (e) => { e.stopPropagation(); flipEvidence(e.target.checked); });
            const varHdr = document.createElement('div');
            varHdr.className = 'ufv-am-hdr-row';
            varHdr.append(varToggle, eviWrap, varRight);
            varSection.append(varHdr, varBody);
            body.appendChild(varSection);
        }

        // ── PTM–Variant Proximity ────────────────────────────────────────────
        // PTM–Variant Proximity: residue-centric view — all distances are FROM the selected residue.
        // Slider 1 controls which PTMs appear (within X Å); Slider 2 controls pathogenic variants.
        {
            const initPtmR = s.settings.proxPtmRadius || 8;
            const initVarR = s.settings.proxVarRadius || 12;
            const geo = StructureViewer.residueGeometry?.() || [];

            const _compute = (ptmR, varR) =>
                UFVAnalysis.computeResidueProximity(pos, s.ptms, s.variants, geo, ptmR, varR);

            const init = _compute(initPtmR, initVarR);

            if (init.nearbyPtms.length || init.nearbyVariants.length) {
                const proxSection = document.createElement('div');
                proxSection.className = 'ufv-am-section';

                const proxToggle = document.createElement('button');
                proxToggle.className = 'ufv-am-toggle';
                const proxHdrLeft = document.createElement('span');
                proxHdrLeft.className = 'ufv-am-hdr-left';
                proxHdrLeft.appendChild(document.createTextNode('PTM–Variant Proximity'));
                const proxArrow = document.createElement('span');
                proxArrow.className = 'ufv-am-arrow';
                proxArrow.textContent = '▾';
                const proxHdrRight = document.createElement('span');
                proxHdrRight.className = 'ufv-am-hdr-right';
                proxHdrRight.append(proxArrow);
                proxToggle.append(proxHdrLeft, proxHdrRight);

                const proxBody = document.createElement('div');
                proxBody.className = 'ufv-am-body';

                const ptmSlider = document.createElement('input');
                ptmSlider.type = 'range'; ptmSlider.className = 'ufv-nearby-slider';
                ptmSlider.min = 2; ptmSlider.max = 30; ptmSlider.step = 1; ptmSlider.value = initPtmR;
                const ptmValEl = document.createElement('span');
                ptmValEl.className = 'ufv-prox-slider-val'; ptmValEl.textContent = initPtmR + ' Å';

                const varSlider = document.createElement('input');
                varSlider.type = 'range'; varSlider.className = 'ufv-nearby-slider';
                varSlider.min = 2; varSlider.max = 30; varSlider.step = 1; varSlider.value = initVarR;
                const varValEl = document.createElement('span');
                varValEl.className = 'ufv-prox-slider-val'; varValEl.textContent = initVarR + ' Å';

                const makeSliderRow = (lbl, slider, valEl) => {
                    const row = document.createElement('div');
                    row.className = 'ufv-prox-slider-item';
                    const l = document.createElement('span');
                    l.className = 'ufv-detail-lbl'; l.textContent = lbl;
                    row.append(l, slider, valEl);
                    return row;
                };

                // PTM slider sits immediately above the PTMs list
                const ptmContent = document.createElement('div');
                proxBody.append(makeSliderRow('PTMs within', ptmSlider, ptmValEl), ptmContent);

                // Variant slider sits immediately above the variants list
                const varContent = document.createElement('div');
                proxBody.append(makeSliderRow('Variants within', varSlider, varValEl), varContent);

                const renderPtmContent = (nearbyPtms) => {
                    ptmContent.textContent = '';
                    if (!nearbyPtms.length) return;
                    const hdr = document.createElement('div');
                    hdr.className = 'ufv-prox-tier-hdr'; hdr.style.color = '#4db6ac';
                    hdr.textContent = 'Nearby PTMs';
                    ptmContent.appendChild(hdr);
                    const grid = document.createElement('div');
                    grid.className = 'ufv-am-grid';
                    nearbyPtms.forEach(({ ptm, dist }) => {
                        const cell = document.createElement('button');
                        cell.className = 'ufv-am-cell ufv-prox-clickable';
                        cell.title = (ptm.description || ptm.category || '').trim() + ' (click to focus)';
                        const nameSpan = document.createElement('span');
                        nameSpan.className = 'ufv-am-cell-mut';
                        nameSpan.style.color = ptm.color || '#4db6ac';
                        nameSpan.textContent = (s.sequence?.[ptm.position - 1] || '') + ptm.position;
                        const distSpan = document.createElement('span');
                        distSpan.className = 'ufv-am-cell-sc';
                        distSpan.textContent = dist.toFixed(1) + ' Å';
                        cell.append(nameSpan, distSpan);
                        cell.addEventListener('click', () => onClick({ position: ptm.position }, 'focus', s.selectedChain));
                        grid.appendChild(cell);
                    });
                    ptmContent.appendChild(grid);
                };

                const renderVarContent = (nearbyVariants) => {
                    varContent.textContent = '';
                    if (!nearbyVariants.length) return;
                    const hdr = document.createElement('div');
                    hdr.className = 'ufv-prox-tier-hdr'; hdr.style.color = '#ef5350';
                    hdr.textContent = 'Pathogenic variants';
                    varContent.appendChild(hdr);
                    const grid = document.createElement('div');
                    grid.className = 'ufv-am-grid';
                    nearbyVariants.forEach(({ variant: v, dist }) => {
                        const cell = document.createElement('button');
                        cell.className = 'ufv-am-cell ufv-prox-clickable';
                        cell.title = (v.clinVarSignificance || v.consequence || '').trim() + ' (click to focus)';
                        const mutSpan = document.createElement('span');
                        mutSpan.className = 'ufv-am-cell-mut';
                        mutSpan.style.color = v.consequenceColor || '#ef5350';
                        mutSpan.textContent = `${v.wildType || ''}${v.position}${v.mutant || ''}`;
                        const distSpan = document.createElement('span');
                        distSpan.className = 'ufv-am-cell-sc';
                        distSpan.textContent = dist.toFixed(1) + ' Å';
                        cell.append(mutSpan, distSpan);
                        cell.addEventListener('click', () => onClick({ position: v.position }, 'focus', s.selectedChain));
                        grid.appendChild(cell);
                    });
                    varContent.appendChild(grid);
                };

                renderPtmContent(init.nearbyPtms);
                renderVarContent(init.nearbyVariants);

                ptmSlider.addEventListener('input', () => { ptmValEl.textContent = ptmSlider.value + ' Å'; });
                varSlider.addEventListener('input', () => { varValEl.textContent = varSlider.value + ' Å'; });
                ptmSlider.addEventListener('change', () => {
                    const ptmR = Number(ptmSlider.value);
                    s.settings.proxPtmRadius = ptmR; UFVState.saveSettings({ proxPtmRadius: ptmR });
                    renderPtmContent(_compute(ptmR, Number(varSlider.value)).nearbyPtms);
                });
                varSlider.addEventListener('change', () => {
                    const varR = Number(varSlider.value);
                    s.settings.proxVarRadius = varR; UFVState.saveSettings({ proxVarRadius: varR });
                    renderVarContent(_compute(Number(ptmSlider.value), varR).nearbyVariants);
                });

                proxToggle.addEventListener('click', () => {
                    const open = proxBody.classList.toggle('show');
                    proxArrow.textContent = open ? '▴' : '▾';
                });
                proxSection.append(proxToggle, proxBody);
                body.appendChild(proxSection);
            }
        }

        // ── Prediction: AlphaMissense (always) + ProtVar EVE/ESM1b/FoldX (lazy) per substitution ──
        // One boxed table. AM is from the local CSV; the ProtVar columns fetch on first expand.
        {
            const AM_AAS = 'ACDEFGHIKLMNPQRSTVWY';
            const amEntries = [];
            if (s.amMap && wt) {
                for (const mut of AM_AAS) {
                    if (mut === wt) continue;
                    const sc = s.amMap.get(`${wt}${pos}${mut}`);
                    if (Number.isFinite(sc)) amEntries.push({ mut, am: sc });
                }
            }
            const observed = new Set((variants || []).map(v => v.mutant).filter(Boolean));
            // "Binding & pockets" section ABOVE the Predictions dropdown — ONE collapsible holding two
            // sub-blocks: predicted pockets (AutoSite/ProtVar) and experimental binding (PDBe-KB). Both load
            // async; the section is removed if NEITHER source has anything for this residue. (Target
            // tractability is whole-protein, so it lives in the header dropdown, not here.)
            const bindSec = document.createElement('div'); bindSec.className = 'ufv-am-section';
            const bindToggle = document.createElement('button'); bindToggle.className = 'ufv-am-toggle';
            bindToggle.innerHTML = '<span class="ufv-am-hdr-left">Binding &amp; pockets</span><span class="ufv-am-hdr-right"><span class="ufv-am-arrow">▾</span></span>';
            const bindBody = document.createElement('div'); bindBody.className = 'ufv-am-body';
            bindToggle.addEventListener('click', () => { const open = bindBody.classList.toggle('show'); bindToggle.querySelector('.ufv-am-arrow').textContent = open ? '▴' : '▾'; });
            const predSub = document.createElement('div'); const expSub = document.createElement('div');
            bindBody.append(predSub, expSub);
            bindSec.append(bindToggle, bindBody);
            body.appendChild(bindSec);
            let bindPending = 0;
            const finishBind = () => { if (--bindPending <= 0 && !predSub.childElementCount && !expSub.childElementCount) bindSec.remove(); };
            if (UFVApi.fetchProtVarPocket) {
                bindPending++;
                UFVApi.fetchProtVarPocket(s.uniprotId, pos).then(pk => { if (s.selectedResidue === pos) renderBindingPocket(predSub, pk); }).catch(() => {}).finally(finishBind);
            }
            if (UFVApi.fetchPdbeKbContext) {
                bindPending++;
                UFVApi.fetchPdbeKbContext(s.uniprotId).then(ctx => { if (s.selectedResidue === pos) renderStructuralContext(expSub, ctx, pos); }).catch(() => {}).finally(finishBind);
            }
            const predSection = document.createElement('div');
            predSection.className = 'ufv-am-section';
            const predToggle = document.createElement('button');
            predToggle.className = 'ufv-am-toggle';
            const avg = amEntries.length ? amEntries.reduce((a, e) => a + e.am, 0) / amEntries.length : null;
            const avgColor = avg == null ? '' : avg >= 0.564 ? '#ef5350' : avg >= 0.34 ? '#ffa726' : '#66bb6a';
            predToggle.innerHTML = `<span class="ufv-am-hdr-left">Predictions</span><span class="ufv-am-hdr-right">${avg != null ? `<span class="ufv-am-avg" style="color:${avgColor}" title="AlphaMissense mean">${avg.toFixed(3)}</span>` : ''}<span class="ufv-am-arrow">▾</span></span>`;
            predToggle.querySelector('.ufv-am-hdr-left').appendChild(makeInfoIcon('Variant-effect predictors', PREDICTOR_INFO));
            const predBody = document.createElement('div');
            predBody.className = 'ufv-am-body';
            let loaded = false;
            predToggle.addEventListener('click', async () => {
                const open = predBody.classList.toggle('show');
                predToggle.querySelector('.ufv-am-arrow').textContent = open ? '▴' : '▾';
                if (!open || loaded) return;
                loaded = true;
                renderPredictionTable(predBody, amEntries, null, wt, pos, s.uniprotId, observed); // AM-only + loading
                let pv = null;
                try { pv = await UFVApi.fetchProtVar(s.uniprotId, pos, s.sequence?.[pos - 1]); } catch (_e) { pv = null; }
                if (s.selectedResidue !== pos) return;
                renderPredictionTable(predBody, amEntries, pv, wt, pos, s.uniprotId, observed); // full
                // (Pockets are already rendered richly from fetchProtVarPocket into the Binding & pockets
                // section; the full-predictor fetch's legacy pocket field would only downgrade it.)
            });
            predSection.append(predToggle, predBody);
            body.appendChild(predSection);
        }

        byId('ufv-details').classList.add('show');
        renderLegend(getColorMode());
        renderSequence();
    }

    // ProtVar's per-substitution EVE/ESM/AM arrays are the 19 non-wild-type AAs in alphabetical order
    // (verified vs the AlphaMissense CSV). Zip back to {mut: value}; bail on length mismatch.
    const PV_AA = 'ACDEFGHIKLMNPQRSTVWY';
    function pvByMut(raw, wt) {
        if (!Array.isArray(raw) || !wt) return {};
        const muts = [...PV_AA].filter(a => a !== wt);
        if (raw.length !== muts.length) return {};
        const out = {}; muts.forEach((m, i) => { out[m] = raw[i]; });
        return out;
    }

    const amColor = sc => sc >= 0.564 ? '#e53935' : sc >= 0.34 ? '#fb8c00' : '#43a047';

    // Binding-pocket info as a collapsible "Pocket" section (matching the Predictions/Proximity sections),
    // rendered ABOVE the Predictions dropdown — not a dangling labelled row. Only appears when ProtVar flags
    // this residue as lining a predicted pocket (so it's absent on residues that don't line one). Empty clears.
    function renderBindingPocket(container, pk) {
        container.textContent = '';
        if (!pk) return;
        const s = UFVState.state;
        const chain = s.selectedChain ?? null;
        // Pocket-residue chip colours use the SAME annotation map as the Nearby list (site > PTM > variant
        // pathogenicity), so the colours read consistently instead of an ad-hoc accent/disease mix.
        const ann = buildAnnotationMap();
        // Sub-block inside the parent "Binding & pockets" section (the parent owns the collapsible).
        const hdr = document.createElement('div'); hdr.className = 'ufv-sub-hdr';
        hdr.append(document.createTextNode('Predicted pockets'), makeInfoIcon('Predicted pockets (AutoSite)', POCKET_INFO));
        container.appendChild(hdr);
        const body = container; // render the pocket blocks directly into the container

        // Re-focus the residue this panel belongs to (restores the camera + focus sticks we had before a
        // pocket was shown) — called when the last pocket is deselected.
        const restoreResidueFocus = () => {
            if (s.selectedResidue != null && StructureViewer.currentStructure) {
                const nb = StructureViewer.focusResidue(s.selectedResidue, s.selectedChain, { annotatedResidues: buildAnnotationMap() }, { showOtherSpheres: _showOtherSpheres, rezoom: true });
                if (nb) s.nearbyResidues = nb;
            }
        };

        // ONE block per detected pocket (the number of blocks IS the pocket count, so no redundant aggregate
        // header). Clicking a pocket's residues box highlights it in 3-D (translucent surface + annotation-
        // coloured sticks) and lists its residues; opening a different pocket collapses the previous one.
        const pockets = (pk.pockets && pk.pockets.length)
            ? pk.pockets
            : [{ pocketId: null, buriedness: pk.buriedness, score: pk.score, resid: [] }]; // legacy fallback
        const closeFns = []; // mutual-exclusion: collapse every other pocket block when one opens
        pockets.forEach((p, i) => {
            const resids = p.resid || [];
            const block = document.createElement('div'); block.className = 'ufv-pocket-block';
            const lbl = document.createElement('div'); lbl.className = 'ufv-pocket-block-lbl';
            lbl.textContent = p.pocketId != null ? `Pocket ${p.pocketId}` : `Pocket ${i + 1}`;
            block.appendChild(lbl);

            const children = document.createElement('div'); children.className = 'ufv-filter-children ufv-collapsed';
            const bw = document.createElement('div'); bw.className = 'ufv-pocket-stats';
            // Each stat = a small label ABOVE the number, with only the number inside the box (per request).
            const mkBox = (v, l, t, onClickFn) => {
                const wrap = document.createElement('span'); wrap.className = 'ufv-pstat'; wrap.title = t;
                const lbl = document.createElement('span'); lbl.className = 'ufv-pstat-lbl'; lbl.textContent = l;
                const box = document.createElement('span'); box.className = 'ufv-pstat-box'; box.textContent = v;
                wrap.append(lbl, box);
                if (onClickFn) { wrap.classList.add('ufv-stat-box-btn'); wrap.addEventListener('click', onClickFn); }
                return wrap;
            };
            const resBox = mkBox(resids.length, 'res',
                'Residues lining this pocket. Click to highlight it in 3-D (surface + sticks) and list them', null);
            const collapse = () => { children.classList.add('ufv-collapsed'); resBox.classList.remove('ufv-open'); };
            closeFns.push(collapse);
            if (resids.length) {
                resBox.classList.add('ufv-stat-box-btn');
                resBox.addEventListener('click', () => {
                    const willOpen = children.classList.contains('ufv-collapsed');
                    closeFns.forEach(fn => { if (fn !== collapse) fn(); }); // mutual exclusion
                    if (willOpen) {
                        children.classList.remove('ufv-collapsed'); resBox.classList.add('ufv-open');
                        StructureViewer.showPocket?.(resids, chain, ann);
                    } else {
                        collapse();
                        StructureViewer.clearPocket?.();
                        restoreResidueFocus();   // last pocket closed → back to the residue focus we had before
                    }
                });
            }
            bw.appendChild(resBox);
            if (p.buriedness != null) bw.appendChild(mkBox(p.buriedness.toFixed(2), 'buried', 'Buriedness — how enclosed the pocket is (0 = open surface, 1 = fully buried).'));
            if (p.meanPlddt != null) bw.appendChild(mkBox(Math.round(p.meanPlddt), 'pLDDT', 'Mean AlphaFold confidence over the pocket’s residues (0–100).'));
            if (p.radGyration != null) bw.appendChild(mkBox(p.radGyration.toFixed(1), 'Rg', 'Radius of gyration (Å) — the pocket’s spatial extent. Smaller = a tighter, more compact cavity.'));
            if (p.energyPerVol != null) bw.appendChild(mkBox(p.energyPerVol.toFixed(2), 'E/V', 'Interaction-affinity energy per unit pocket volume — higher = a denser, more ligand-favourable cavity.'));
            if (p.score != null) bw.appendChild(mkBox(Math.round(p.score), 'score', POCKET_SCORE_TIP));
            block.appendChild(bw);

            // Pocket character radar (hydrophobic / aromatic / acidic / basic / polar of the lining
            // residues), computed from the sequence — same read as the ligand pocket-evidence panel.
            const axes = pocketCompositionAxes(resids, s.sequence);
            if (axes) {
                const spider = document.createElement('div'); spider.className = 'ufv-pocket-spider ufv-pocket-spider-mini';
                spider.title = axes.map(a => `${Math.round(a.pct * 100)}% ${a.label.toLowerCase()}`).join(' · ');
                spider.appendChild(makeSpiderPlot(axes));
                children.appendChild(spider);
            }
            resids.forEach((rp, j) => {
                const span = document.createElement('span');
                span.className = 'ufv-nearby-res';
                span.textContent = (s.sequence?.[rp - 1] || '') + rp;
                const c = ann.get(rp)?.color; if (c) span.style.color = c; // same colouring as the Nearby list
                span.title = `Focus residue ${rp}`;
                span.addEventListener('click', () => onClick({ position: rp }, null, chain));
                children.appendChild(span);
                if (j < resids.length - 1) children.appendChild(document.createTextNode(', '));
            });
            // Per-pocket motif search — shown only INSIDE the expanded view (i.e. once the residues are
            // selected), at the bottom of the residue list. RCSB Structure Motif Search on this pocket.
            const _st = StructureViewer.currentStructure;
            if (resids.length >= 2 && (_st?.pdbId || isAlphaFoldBased(_st))) {
                const mrow = document.createElement('div'); mrow.className = 'ufv-pocket-motif';
                const ma = document.createElement('a');
                ma.href = 'javascript:void 0'; ma.className = 'ufv-ligand-link';
                ma.textContent = 'Find Similar Motifs ↗';
                ma.title = 'Search the PDB for structures with a similar 3-D arrangement of this pocket’s residues (RCSB Structure Motif Search)';
                ma.addEventListener('click', (e) => {
                    e.preventDefault(); ma.textContent = 'searching…';
                    openPocketSearch(new Set(resids)).then(ok => { ma.textContent = ok ? 'Find Similar Motifs ↗' : 'motif search unavailable'; }).catch(() => { ma.textContent = 'motif search unavailable'; });
                });
                mrow.appendChild(ma); children.appendChild(mrow);
            }
            block.appendChild(children);
            body.appendChild(block);
        });

        // (Attribution to AutoSite / ProtVar lives in the (i) popover.)
    }
    // AutoSite's per-pocket score (ProtVar predicts pockets with AutoSite).
    const POCKET_SCORE_TIP = 'AutoSite pocket score — ranks how pronounced and ligand-favourable the cavity is; higher = a larger, better-defined pocket.';

    // [prototype] PDBe-KB per-residue "Structural context": is this residue a ligand-binding site and/or a
    // protein-protein interface across all experimental PDB structures of this protein. ctx = { ligandByPos,
    // interfaceByPos } from UFVApi.fetchPdbeKbContext; only renders when the residue actually has evidence.
    function renderStructuralContext(container, ctx, pos) {
        container.textContent = '';
        if (!ctx) return;
        const ligs = ctx.ligandByPos.get(pos) || [];
        const ifaces = ctx.interfaceByPos.get(pos) || [];
        if (!ligs.length && !ifaces.length) return;
        // Sub-block inside the parent "Binding & pockets" section (the parent owns the collapsible).
        const hdr = document.createElement('div'); hdr.className = 'ufv-sub-hdr';
        hdr.append(document.createTextNode('Experimental binding'), makeInfoIcon('Experimental binding (PDBe-KB)', STRUCT_CTX_INFO));
        container.appendChild(hdr);
        const body = container;
        const row = (label, valueNode, cls) => {
            const r = document.createElement('div'); r.className = 'ufv-detail-row';
            const l = document.createElement('span'); l.className = 'ufv-detail-lbl'; l.textContent = label;
            const v = document.createElement('span'); v.className = 'ufv-detail-val' + (cls ? ' ' + cls : '');
            if (typeof valueNode === 'string') v.textContent = valueNode; else v.appendChild(valueNode);
            r.append(l, v); body.appendChild(r);
        };
        if (ligs.length) {
            // De-dup ligand codes; most-observed first. Each links to a representative PDB STRUCTURE on RCSB
            // (the actual complex with the ligand bound), not the chemical-dictionary page.
            const seen = new Map();
            ligs.forEach(l => { const e = seen.get(l.code); if (!e || l.pdbCount > e.pdbCount) seen.set(l.code, l); });
            const top = [...seen.values()].sort((a, b) => b.pdbCount - a.pdbCount);
            const wrap = document.createElement('span');
            top.slice(0, 8).forEach((l, i) => {
                if (i) wrap.appendChild(document.createTextNode(', '));
                const a = document.createElement('a');
                a.className = 'ufv-ligand-link'; a.target = '_blank'; a.rel = 'noopener noreferrer';
                a.href = l.pdb ? `https://www.rcsb.org/structure/${String(l.pdb).toUpperCase()}` : `https://www.rcsb.org/ligand/${l.code}`;
                a.textContent = l.code;
                a.title = (l.name || l.code) + ` — bound in ${l.pdbCount} PDB structure${l.pdbCount === 1 ? '' : 's'}` + (l.pdb ? ` (e.g. ${String(l.pdb).toUpperCase()})` : '');
                wrap.appendChild(a);
            });
            if (top.length > 8) wrap.appendChild(document.createTextNode(` +${top.length - 8}`));
            row('Ligands', wrap, 'ufv-ctx-ligand');
        }
        if (ifaces.length) {
            const seen = new Map();
            ifaces.forEach(p => { const e = seen.get(p.accession); if (!e || p.pdbCount > e.pdbCount) seen.set(p.accession, p); });
            const top = [...seen.values()].sort((a, b) => b.pdbCount - a.pdbCount);
            const wrap = document.createElement('span');
            top.slice(0, 4).forEach((p, i) => {
                if (i) wrap.appendChild(document.createTextNode('; '));
                const a = document.createElement('a');
                a.className = 'ufv-ligand-link'; a.target = '_blank'; a.rel = 'noopener noreferrer';
                a.href = `https://www.rcsb.org/structure/${String(p.pdb || '').toUpperCase()}`;
                if (!p.pdb) { a.href = `https://www.uniprot.org/uniprotkb/${p.accession}/entry`; }
                a.textContent = p.name || p.accession;
                a.title = `Interface seen in ${p.pdbCount} PDB complex${p.pdbCount === 1 ? '' : 'es'}` + (p.pdb ? ` (e.g. ${String(p.pdb).toUpperCase()})` : '');
                wrap.appendChild(a);
            });
            if (top.length > 4) wrap.appendChild(document.createTextNode(` +${top.length - 4}`));
            row('PPI interface', wrap, 'ufv-ctx-iface');
        }
    }
    const STRUCT_CTX_INFO = [
        ['Experimental binding (PDBe-KB)', 'Aggregated from every experimental PDB structure of this protein, mapped to UniProt residues — shown whether you’re viewing the AlphaFold model or an experimental structure.', 'https://www.ebi.ac.uk/pdbe/pdbe-kb'],
        ['Ligands', 'Small molecules / ions observed bound at this residue. Each links to a representative PDB structure (RCSB) where it’s bound.', ''],
        ['PPI interface', 'Partner proteins whose interface includes this residue — flagged even when the loaded structure is a monomer.', ''],
    ].map(([name, desc, url]) => `<p><b>${name}</b> — ${desc}${url ? ` <a href="${url}" target="_blank" rel="noopener noreferrer">source ↗</a>` : ''}</p>`).join('');

    // Whole-protein TRACTABILITY (Open Targets) as a HEADER dropdown — sits with the other protein-level
    // controls (Partners / Other chains), not in the per-residue panel. "Tractability" (Open Targets' term),
    // NOT "druggability": it spans small-molecule, antibody AND PROTAC amenability, so "druggability" (a
    // small-molecule, pocket-centric idea) would be a misnomer.
    function updateTractabilityNav() {
        const host = byId('ufv-modal-heading')?.parentElement; if (!host) return;
        byId('ufv-tract-wrap')?.remove(); // reset on (re)load / protein switch
        const acc = UFVState.state.uniprotId;
        if (!acc || !UFVApi.fetchOpenTargets) return;
        UFVApi.fetchOpenTargets(acc).then(ot => {
            if (UFVState.state.uniprotId !== acc) return;          // protein changed while fetching
            if (!ot || !ot.groups || !ot.groups.length) return;     // no tractability → no dropdown
            if (byId('ufv-tract-wrap')) return;
            const wrap = document.createElement('span');
            wrap.id = 'ufv-tract-wrap'; wrap.className = 'ufv-dropdown ufv-tract-dropdown';
            const btn = document.createElement('button');
            btn.type = 'button'; btn.className = 'ufv-dropdown-btn';
            btn.title = 'Open Targets tractability — how amenable this whole protein is to a small molecule, antibody or PROTAC';
            btn.innerHTML = 'Tractability <span class="ufv-dropdown-caret">▾</span>';
            const menu = document.createElement('div'); menu.className = 'ufv-dropdown-menu ufv-tract-menu';
            const hdr = document.createElement('div'); hdr.className = 'ufv-tract-hdr';
            hdr.append(document.createTextNode('Target tractability'), makeInfoIcon('Target tractability (Open Targets)', DRUGGABILITY_INFO));
            menu.appendChild(hdr);
            ot.groups.forEach(g => {
                const block = document.createElement('div'); block.className = 'ufv-drug-block';
                const h = document.createElement('div'); h.className = 'ufv-drug-mod'; h.textContent = g.modality;
                const chips = document.createElement('div'); chips.className = 'ufv-drug-chips';
                g.labels.forEach(lbl => {
                    const chip = document.createElement('span'); chip.className = 'ufv-ctx-chip'; chip.textContent = lbl;
                    if (TRACTABILITY_TIPS[lbl]) chip.title = TRACTABILITY_TIPS[lbl];
                    chips.appendChild(chip);
                });
                block.append(h, chips); menu.appendChild(block);
            });
            const src = document.createElement('div'); src.className = 'ufv-pocket-motif';
            const a = document.createElement('a'); a.className = 'ufv-ligand-link'; a.target = '_blank'; a.rel = 'noopener noreferrer';
            a.href = ot.url; a.textContent = `${ot.symbol} on Open Targets ↗`;
            src.appendChild(a); menu.appendChild(src);
            wrap.append(btn, menu);
            host.appendChild(wrap);
            wireDropdown(wrap, btn);
            menu.addEventListener('click', e => e.stopPropagation()); // clicks inside keep the menu open
        }).catch(() => {});
    }
    // Short hover explanations for Open Targets tractability buckets (the satisfied ones we show as chips).
    const TRACTABILITY_TIPS = {
        'Approved Drug': 'An approved drug of this modality acts on the target.',
        'Advanced Clinical': 'A drug of this modality against the target has reached advanced clinical trials (Phase 2/3).',
        'Phase 1 Clinical': 'A drug of this modality against the target has reached Phase 1.',
        'Structure with Ligand': 'A 3-D structure of the target with a bound small-molecule ligand exists.',
        'High-Quality Ligand': 'A high-quality, drug-like small-molecule ligand is known for the target.',
        'High-Quality Pocket': 'A high-confidence druggable pocket was detected on the structure.',
        'Med-Quality Pocket': 'A medium-confidence druggable pocket was detected.',
        'Druggable Family': 'The target belongs to a protein family historically amenable to small molecules.',
        'Small Molecule Binder': 'A small molecule that binds the target is known (a handle for a PROTAC).',
        'UniProt loc high conf': 'UniProt subcellular location indicates the target is accessible (cell-surface/secreted), high confidence.',
        'UniProt loc med conf': 'UniProt subcellular location indicates accessibility, medium confidence.',
        'GO CC high conf': 'GO cellular-component evidence supports an accessible localisation, high confidence.',
        'GO CC med conf': 'GO cellular-component evidence supports accessibility, medium confidence.',
        'UniProt SigP or TMHMM': 'Signal-peptide / transmembrane prediction suggests the target is accessible to antibodies.',
        'Human Protein Atlas loc': 'Human Protein Atlas localisation supports antibody accessibility.',
        'Literature': 'PROTAC/degrader evidence in the literature for this target.',
        'UniProt Ubiquitination': 'UniProt annotates the target as ubiquitinated (needed for PROTAC degradation).',
        'Database Ubiquitination': 'Databases report the target as ubiquitinated.',
        'Half-life Data': 'Protein half-life data is available (relevant to degrader strategies).',
    };
    // Single-entry popover — no <b>name</b> prefix (it would just repeat the popover's own title).
    const DRUGGABILITY_INFO = '<p>How amenable the whole protein is to a small molecule, antibody or PROTAC. Each chip is a satisfied tractability bucket — hover for detail. <a href="https://platform-docs.opentargets.org/target/tractability" target="_blank" rel="noopener noreferrer">source ↗</a></p>';

    // Amino-acid CLASS composition of a residue set from the canonical sequence (same classification the
    // ligand pocket-evidence panel uses — computed on OUR side, no extra fetch). Drives the pocket radar.
    const _POCKET_AA_CLASS = {};
    'AVLIMFWPGC'.split('').forEach(a => (_POCKET_AA_CLASS[a] = 'hydrophobic'));
    'STNQY'.split('').forEach(a => (_POCKET_AA_CLASS[a] = 'polar'));
    'DE'.split('').forEach(a => (_POCKET_AA_CLASS[a] = 'acidic'));
    'KRH'.split('').forEach(a => (_POCKET_AA_CLASS[a] = 'basic'));
    function pocketCompositionAxes(resids, seq) {
        if (!seq) return null;
        const comp = { hydrophobic: 0, polar: 0, acidic: 0, basic: 0, aromatic: 0, total: 0 };
        resids.forEach(p => { const aa = seq[p - 1]; if (!aa) return; comp.total++; const c = _POCKET_AA_CLASS[aa]; if (c) comp[c]++; if ('FWYH'.includes(aa)) comp.aromatic++; });
        if (!comp.total) return null;
        const frac = n => n / comp.total;
        return [
            { label: 'Hydrophobic', pct: frac(comp.hydrophobic) },
            { label: 'Aromatic', pct: frac(comp.aromatic) },
            { label: 'Acidic', pct: frac(comp.acidic) },
            { label: 'Basic', pct: frac(comp.basic) },
            { label: 'Polar', pct: frac(comp.polar) },
        ];
    }

    // One boxed table: sub | AM (always) | EVE | ESM1b | FoldX (when ProtVar is loaded). pv === null
    // while ProtVar fetches (AM columns render immediately so the table isn't empty).
    function renderPredictionTable(container, amEntries, pv, wt, pos, acc, observed) {
        container.textContent = '';
        const sgn = (x) => { const v = x.toFixed(1); return (x > 0 && v !== '0.0') ? '+' + v : v; }; // no "+0.0"
        const sc = pv?.score || {}, fx = pv?.foldx || {};
        const havePv = !!(pv && (pv.score || pv.foldx));
        // Per-position rows above the table.
        if (havePv && sc.conservation != null) {
            const r = row('Conservation', sc.conservation.toFixed(2));
            r.title = 'Evolutionary conservation across homologues (0 = variable, 1 = invariant)';
            container.appendChild(r);
        }
        if (havePv && sc.m3d) {
            const dmg = String(sc.m3d.prediction || '').toLowerCase().startsWith('damag');
            const feat = (sc.m3d.feature && sc.m3d.feature !== '-') ? ' — ' + sc.m3d.feature : ''; // no " — -" when there's no feature
            const r = row('M3D', `${sc.m3d.prediction || '?'}${feat}`, dmg ? '#e53935' : '#43a047');
            r.title = 'Missense3D: predicted structural consequence';
            container.appendChild(r);
        }
        if (!havePv && pv !== null) { /* fetch finished with no data: AM-only table still shown below */ }

        const amByMut = {}; (amEntries || []).forEach(e => { amByMut[e.mut] = e.am; });
        // CADD is keyed directly by mutant AA (genomic SNV map), unlike EVE/ESM which arrive as ordered arrays.
        const cadd = sc.caddByMut || {};
        const eve = pvByMut(sc.eveRaw, wt), evc = pvByMut(sc.eveClsRaw, wt), esm = pvByMut(sc.esmRaw, wt), ddg = fx.byMut || {};
        const rows = [];
        [...PV_AA].filter(a => a !== wt).forEach(mt => {
            const am = amByMut[mt], e = eve[mt], es = esm[mt], c = cadd[mt], d = ddg[mt];
            if (am == null && e == null && es == null && c == null && d == null) return;
            rows.push({ mt, am, e, ec: evc[mt], es, c, d });
        });
        rows.sort((a, b) => (b.am ?? -2) - (a.am ?? -2) || (b.e ?? -2) - (a.e ?? -2));
        if (!rows.length) { container.appendChild(row('Predictions', 'no per-substitution scores for this position.')); return; }

        // Theme-aware chip (matches the old AlphaMissense cells) so it reads well in dark mode.
        const boxTd = (html, color) => { const td = document.createElement('td'); td.style.cssText = 'background:var(--ufv-bg-hover);color:var(--ufv-text-primary);padding:2px 7px;text-align:center;border-radius:3px;font-variant-numeric:tabular-nums;'; td.innerHTML = html; if (color) td.style.color = color; return td; };
        // Only show predictor columns that have data — EVE covers ~3000 genes, so a dash column would
        // just look like a failure. AM is always present (local CSV).
        const cols = [{ h: 'AM', t: 'AlphaMissense (0–1 pathogenicity)', cell: r => boxTd(r.am != null ? r.am.toFixed(2) : '–', r.am != null ? amColor(r.am) : '') }];
        if (havePv && rows.some(r => r.e != null)) cols.push({ h: 'EVE', t: 'EVE evolutionary model (0–1)', cell: r => boxTd(r.e != null ? r.e.toFixed(2) : '–', r.e != null && String(r.ec || '').toUpperCase() === 'PATHOGENIC' ? '#e53935' : '') });
        if (rows.some(r => r.c != null)) cols.push({ h: 'CADD', t: 'CADD deleteriousness, phred-scaled (>20 ≈ top 1% most deleterious, >30 ≈ top 0.1%). Only for substitutions reachable by a single nucleotide change.', cell: r => boxTd(r.c != null ? r.c.toFixed(1) : '–', r.c == null ? '' : r.c >= 25 ? '#e53935' : r.c >= 20 ? '#fb8c00' : '') });
        if (havePv && rows.some(r => r.es != null)) cols.push({ h: 'ESM1b', t: 'ESM-1b pathogenicity (lower / more negative = more pathogenic)', cell: r => boxTd(r.es != null ? r.es.toFixed(1) : '–', r.es == null ? '' : r.es <= -10 ? '#e53935' : r.es <= -7.5 ? '#fb8c00' : '') });
        // FoldX ΔΔG: destabilising > 0 (red/orange), stabilising < -1 (green), ~neutral (-1…1) uncoloured.
        if (havePv && rows.some(r => r.d != null)) cols.push({ h: 'FoldX', t: 'FoldX ΔΔG kcal/mol; >2 destabilising, <-1 stabilising', cell: r => boxTd(r.d != null ? sgn(r.d) : '–', r.d == null ? '' : r.d >= 2 ? '#e53935' : r.d >= 1 ? '#fb8c00' : r.d <= -1 ? '#43a047' : '') });

        const tbl = document.createElement('table');
        tbl.style.cssText = 'border-collapse:separate;border-spacing:4px;font-size:11px;';
        const head = document.createElement('tr'); head.style.color = 'var(--ufv-text-secondary)';
        head.innerHTML = '<td>sub</td>' + cols.map(c => `<td title="${c.t}">${c.h}</td>`).join('')
            + (pv === null ? '<td colspan="4" style="font-weight:400">EVE / CADD / ESM-1b / FoldX loading…</td>' : '');
        tbl.appendChild(head);
        // Colour the "observed variant" dot by that variant's consequence (mutant AA → consequence colour).
        const varColorByMut = {};
        UFVState.state.variants.filter(v => v.position === pos).forEach(v => { if (v.mutant && v.consequenceColor) varColorByMut[v.mutant] = v.consequenceColor; });
        rows.forEach(r => {
            const tr = document.createElement('tr');
            const subTd = document.createElement('td'); subTd.style.whiteSpace = 'nowrap';
            const a = document.createElement('a');
            a.href = `https://www.ebi.ac.uk/ProtVar/query?search=${acc}+${wt}${pos}${r.mt}`;
            a.target = '_blank'; a.rel = 'noopener'; a.style.cssText = 'color:#7b1fa2;text-decoration:none;';
            a.textContent = `${wt}${pos}${r.mt}`;
            subTd.appendChild(a);
            if (observed && observed.has(r.mt)) {
                subTd.style.fontWeight = '700';
                const dot = document.createElement('span'); dot.textContent = ' •';
                if (varColorByMut[r.mt]) dot.style.color = varColorByMut[r.mt]; // consequence colour
                subTd.appendChild(dot);
            }
            tr.appendChild(subTd);
            cols.forEach(c => tr.appendChild(c.cell(r)));
            tbl.appendChild(tr);
        });
        container.appendChild(tbl);
    }

    // Clicking a ligand/cofactor (e.g. an AlphaFill-transplanted molecule) — NOT a protein
    // residue: no hotspot/variant/AlphaMissense info applies. Focus it (zoom, nearby protein
    // residues as sticks, other ligands hidden) and show a chemistry-only detail panel.
    function onLigandClick(lig) {
        const s = UFVState.state;
        s.selectedResidue = null;
        s.selectedChain = null;
        s.selectedLigand = lig;
        s.nearbyResidues = StructureViewer.focusLigand(lig.resn, lig.resi, lig.chain, { showOtherSpheres: _showOtherSpheres, annotatedResidues: buildAnnotationMap() }) || new Set();
        renderLigandPanel(lig);
        buildLigandFilters(); // refresh the list so the selected ligand is boxed
        renderSequence();
    }

    function copyChip(value) {
        const chip = document.createElement('button');
        chip.className = 'ufv-copy-chip';
        chip.title = 'Copy';
        chip.textContent = value;
        chip.addEventListener('click', async () => {
            const ok = await UFVExport.copyText(value);
            if (ok) { const prev = chip.textContent; chip.textContent = 'Copied'; chip.classList.add('copied'); setTimeout(() => { chip.textContent = prev; chip.classList.remove('copied'); }, 1000); }
        });
        return chip;
    }

    // Exact Tanimoto over two PubChem 881-bit substructure fingerprints.
    function tanimoto881(a, b) {
        if (!a || !b) return 0;
        let inter = 0, union = 0;
        for (let k = 0; k < 881; k++) { inter += a[k] & b[k]; union += a[k] | b[k]; }
        return union ? inter / union : 0;
    }

    // Rank the OTHER ligands loaded in the model by Tanimoto similarity to the focused ligand and
    // show them as a collapsible 3-column grid (like AlphaMissense); clicking one focuses it.
    // Similarity uses the published PubChem 2D substructure fingerprint (CACTVS 881-bit keys).
    async function renderSimilarLigands(lig, focusMeta, sectionEl) {
        const s = UFVState.state;
        const others = new Map(); // CCD → representative instance
        s.ligands.forEach(l => { if (l.resn !== lig.resn && !others.has(l.resn)) others.set(l.resn, l); });
        if (!others.size || !focusMeta?.inchikey) { sectionEl.classList.add('ufv-hidden'); return; }
        const focusFp = await UFVApi.getLigandFingerprint(focusMeta.inchikey);
        if (UFVState.state.selectedLigand !== lig) return;
        if (!focusFp) { sectionEl.classList.add('ufv-hidden'); return; }
        // Cap the candidate set — AlphaFill models can carry dozens of transplanted ligands, and
        // each candidate needs an RCSB + a PubChem lookup (both cached).
        const ranked = (await Promise.all([...others.values()].slice(0, 20).map(async inst => {
            const m = await UFVApi.getLigandInfo(inst.resn);
            const fp = m?.inchikey ? await UFVApi.getLigandFingerprint(m.inchikey) : null;
            return { inst, score: tanimoto881(focusFp, fp) };
        }))).filter(e => e.score > 0).sort((a, b) => b.score - a.score);
        if (UFVState.state.selectedLigand !== lig) return;       // user moved on
        if (!ranked.length) { sectionEl.classList.add('ufv-hidden'); return; }

        sectionEl.className = 'ufv-am-section';
        const toggle = document.createElement('button');
        toggle.className = 'ufv-am-toggle';
        toggle.innerHTML = `<span class="ufv-am-hdr-left">Ligand similarity within structure</span><span class="ufv-am-hdr-right"><span class="ufv-am-ratio">${ranked.length}</span><span class="ufv-am-arrow">▾</span></span>`;
        const bodyEl = document.createElement('div');
        bodyEl.className = 'ufv-am-body';
        const grid = document.createElement('div');
        grid.className = 'ufv-am-grid';
        ranked.forEach(({ inst, score }) => {
            const cell = document.createElement('button');
            cell.className = 'ufv-am-cell ufv-sim-cell';
            cell.title = `CACTVS substructure-fingerprint Tanimoto ${score.toFixed(2)} (click to focus)`;
            cell.innerHTML = `<span class="ufv-am-cell-mut">${_esc(inst.resn)}</span><span class="ufv-am-cell-sc">${score.toFixed(2)}</span>`;
            cell.addEventListener('click', () => onLigandClick({ resn: inst.resn, resi: inst.resi, chain: inst.chain }));
            grid.appendChild(cell);
        });
        bodyEl.appendChild(grid);
        toggle.addEventListener('click', () => {
            const open = bodyEl.classList.toggle('show');
            toggle.querySelector('.ufv-am-arrow').textContent = open ? '▴' : '▾';
        });
        sectionEl.textContent = '';
        toggle.querySelector('.ufv-am-arrow').textContent = '▾';
        sectionEl.append(toggle, bodyEl);
    }

    // Pocket evidence panel (M6): aggregate confidence + known-site context for the residues lining a
    // focused ligand. pLDDT (model B-factor) tells you how trustworthy the pocket geometry is; the
    // functional-site overlap tells you whether the ligand sits at a UniProt-annotated active/binding/
    // metal site. Both come from data already loaded — CA B-factors via residueGeometry, sites via s.sites.
    // Compact SVG radar/spider plot for a few normalised axes (0..1). Used for pocket character, where a
    // string of four percentages reads as clutter — the shape conveys the balance at a glance.
    function makeSpiderPlot(axes) {
        const NS = 'http://www.w3.org/2000/svg', N = axes.length, S = 168, cx = S / 2, cy = S / 2, R = 48;
        const ang = i => -Math.PI / 2 + i * (2 * Math.PI / N);
        const pt = (i, r) => [cx + r * Math.cos(ang(i)), cy + r * Math.sin(ang(i))];
        const mk = (tag, attrs) => { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); return e; };
        const svg = mk('svg', { viewBox: `0 0 ${S} ${S}`, width: S, height: S, class: 'ufv-spider' });
        [0.5, 1].forEach(f => svg.appendChild(mk('polygon', { points: axes.map((_, i) => pt(i, R * f).join(',')).join(' '), fill: 'none', stroke: 'var(--ufv-border)', 'stroke-width': 1 })));
        axes.forEach((a, i) => {
            const [x, y] = pt(i, R);
            svg.appendChild(mk('line', { x1: cx, y1: cy, x2: x, y2: y, stroke: 'var(--ufv-border)', 'stroke-width': 1 }));
            const [lx, ly] = pt(i, R + 15);
            const t = mk('text', { x: lx, y: ly, 'text-anchor': 'middle', 'dominant-baseline': 'middle', fill: 'var(--ufv-text-secondary)', 'font-size': 9, 'font-weight': 600 });
            t.textContent = a.label; svg.appendChild(t);
        });
        svg.appendChild(mk('polygon', { points: axes.map((a, i) => pt(i, R * Math.max(0, Math.min(1, a.pct))).join(',')).join(' '), fill: 'rgba(0,137,123,0.30)', stroke: '#00bfa5', 'stroke-width': 2.2 }));
        axes.forEach((a, i) => { const [x, y] = pt(i, R * Math.max(0, Math.min(1, a.pct))); svg.appendChild(mk('circle', { cx: x, cy: y, r: 2.4, fill: '#00bfa5' })); });
        return svg;
    }

    // Structure-based pocket similarity search: RCSB Structure Motif Search ("strucmotif") for PDB structures
    // with a similar 3-D arrangement of the pocket residues. strucmotif keys on mmCIF label ids, which we
    // derive on demand from PDBe (Mol*'s label_seq_id from the PDBe CIF doesn't match RCSB's). Async because
    // the auth→label map is fetched on click. Needs a real PDB entry + 2–10 mappable residues.
    // AlphaFold-based model (AlphaFold OR AlphaFill — its protein coords ARE the AlphaFold model), usable as
    // an RCSB Computed-Structure-Model (AF_AF<acc>F1) for a strucmotif pocket search. Excludes PDB & other CSMs.
    function isAlphaFoldBased(st) { return !!st && !st.pdbId && (st.source === 'AlphaFold' || /alphafill/i.test(st.provider || st.source || st.label || '')); }

    async function openPocketSearch(pocket) {
        const st = StructureViewer.currentStructure;
        if (!st || !pocket || pocket.size < 2) return false;
        const isPdb = !!st.pdbId, isAF = !isPdb && isAlphaFoldBased(st);
        if (!isPdb && !isAF) return false;
        // Pocket residues with CA geometry, picked COMPACTLY (the 8 closest to the pocket centroid) so the
        // motif stays inside strucmotif's distance threshold rather than spanning the whole site.
        const geoByUni = new Map();
        (StructureViewer.residueGeometry?.() || []).forEach(g => { if (g.uniPos != null && g.ca) geoByUni.set(g.uniPos, g); });
        const pts = Array.from(pocket).map(p => geoByUni.get(p)).filter(Boolean);
        if (pts.length < 2) return false;
        const c = pts.reduce((a, g) => ({ x: a.x + g.ca.x, y: a.y + g.ca.y, z: a.z + g.ca.z }), { x: 0, y: 0, z: 0 });
        c.x /= pts.length; c.y /= pts.length; c.z /= pts.length;
        const d2 = g => (g.ca.x - c.x) ** 2 + (g.ca.y - c.y) ** 2 + (g.ca.z - c.z) ** 2;
        const picked = pts.sort((a, b) => d2(a) - d2(b)).slice(0, 8);
        let entryId, rid;
        if (isPdb) {
            if (!UFVApi.fetchPocketLabelIds) return false;
            const labelMap = await UFVApi.fetchPocketLabelIds(st.pdbId);
            rid = []; for (const g of picked) { const lab = labelMap.get(`${g.chain}|${g.resi}`); if (lab) rid.push({ label_asym_id: lab.labelAsym, label_seq_id: lab.labelSeq }); }
            entryId = st.pdbId.toUpperCase();
        } else {
            // AlphaFold computed model: RCSB id AF_AF<acc>F1; label_seq_id == UniProt position, chain A.
            const acc = (UFVState.state.uniprotId || '').replace(/-\d+$/, '');
            if (!acc) return false;
            entryId = `AF_AF${acc}F1`;
            rid = picked.map(g => ({ label_asym_id: 'A', label_seq_id: g.uniPos }));
        }
        if (rid.length < 2) return false;
        const q = { query: { type: 'terminal', service: 'strucmotif', parameters: { value: { entry_id: entryId, residue_ids: rid } } }, return_type: 'assembly' };
        window.open(`https://www.rcsb.org/search?request=${encodeURIComponent(JSON.stringify(q))}`, '_blank', 'noopener');
        return true;
    }

    function renderPocketEvidence(body, pocket, lig) {
        if (!pocket || !pocket.size || !UFVAnalysis.pocketEvidence) return;
        const s = UFVState.state;
        const geo = StructureViewer.residueGeometry?.() || [];
        const ev = UFVAnalysis.pocketEvidence(pocket, geo, s.sites || [], s.sequence);
        // Closest heavy-atom contact (works for any ligand) — computed here so it lives in the Pocket group.
        let clash = null;
        if (lig && UFVAnalysis.ligandClash) {
            const atoms = StructureViewer.viewer?.getModel?.()?.selectedAtoms?.({}) || [];
            clash = UFVAnalysis.ligandClash(lig.resn, lig.resi, lig.chain, atoms);
        }
        if (ev.meanB == null && !ev.siteHits.length && !ev.composition && !clash) return;

        // Collapsible "Pocket" dropdown, styled like the Predictions / Variants dropdowns (collapsed by default).
        const sec = document.createElement('div'); sec.className = 'ufv-am-section';
        const toggle = document.createElement('button'); toggle.className = 'ufv-am-toggle';
        toggle.innerHTML = '<span class="ufv-am-hdr-left">Pocket</span><span class="ufv-am-hdr-right"><span class="ufv-am-arrow">▾</span></span>';
        const box = document.createElement('div'); box.className = 'ufv-am-body';
        toggle.addEventListener('click', () => { const open = box.classList.toggle('show'); toggle.querySelector('.ufv-am-arrow').textContent = open ? '▴' : '▾'; });
        const row = (label, valueNode, cls, title) => {
            const r = document.createElement('div'); r.className = 'ufv-detail-row';
            const l = document.createElement('span'); l.className = 'ufv-detail-lbl'; l.textContent = label;
            const v = document.createElement('span'); v.className = 'ufv-detail-val' + (cls ? ' ' + cls : '');
            if (typeof valueNode === 'string') v.textContent = valueNode; else v.appendChild(valueNode);
            r.append(l, v); if (title) r.title = title; box.appendChild(r); return r;
        };

        if (clash) {
            const cls = clash.band === 'severe' ? 'high' : clash.band === 'suspicious' ? 'moderate' : 'low';
            row('Closest contact', `${clash.minDist.toFixed(2)} Å (${clash.band})`, 'ufv-clash-' + cls);
        }
        if (ev.meanB != null) {
            // B-factor == pLDDT for AlphaFold/AlphaFill/computed models; crystallographic B for PDB. The
            // residue count moves to the row hover so the value itself reads cleanly (colour = band).
            const isExperimental = /pdb/i.test(StructureViewer.currentStructure?.source || '');
            const title = `mean over ${ev.count} pocket residue${ev.count === 1 ? '' : 's'}`;
            if (isExperimental) {
                row('Mean pocket B-factor', ev.meanB.toFixed(1), null, title);
            } else {
                const band = ev.meanB >= 90 ? 'very-high' : ev.meanB >= 70 ? 'confident' : ev.meanB >= 50 ? 'low' : 'very-low';
                row('Mean pocket pLDDT', ev.meanB.toFixed(1), 'ufv-plddt-' + band, title);
            }
        }
        if (ev.siteHits.length) {
            const txt = ev.siteHits.map(h => `${h.description} (${h.position})`).join('; ');
            row('Functional sites', txt, 'ufv-pocket-site');
        }
        if (ev.composition) {
            // Pocket character as a radar over the five residue classes we score (hydrophobic / aromatic /
            // acidic / basic / polar). Acidic+basic are shown separately (not merged into "charged") so the
            // full set is visible. Classes overlap (e.g. aromatics are also hydrophobic), so axes are
            // independent fractions, not a partition. Exact percentages are in the hover title.
            const c = ev.composition, frac = n => c.total ? n / c.total : 0;
            const axes = [
                { label: 'Hydrophobic', pct: frac(c.hydrophobic) },
                { label: 'Aromatic', pct: frac(c.aromatic) },
                { label: 'Acidic', pct: frac(c.acidic) },
                { label: 'Basic', pct: frac(c.basic) },
                { label: 'Polar', pct: frac(c.polar) },
            ];
            const r = document.createElement('div'); r.className = 'ufv-detail-row ufv-pocket-spider-row';
            const l = document.createElement('span'); l.className = 'ufv-detail-lbl'; l.textContent = 'Pocket character';
            const wrap = document.createElement('div'); wrap.className = 'ufv-pocket-spider';
            wrap.title = axes.map(a => `${Math.round(a.pct * 100)}% ${a.label.toLowerCase()}`).join(' · ');
            wrap.appendChild(makeSpiderPlot(axes));
            r.append(l, wrap); box.appendChild(r);
        }
        // Structure-based pocket similarity search (RCSB strucmotif) — for PDB entries and AlphaFold/AlphaFill.
        const _ps = StructureViewer.currentStructure;
        if ((_ps?.pdbId || isAlphaFoldBased(_ps)) && pocket.size >= 2) {
            const subhdr = document.createElement('div'); subhdr.className = 'ufv-pocket-subhdr'; subhdr.textContent = 'Pocket similarity';
            box.appendChild(subhdr);
            const r = document.createElement('div'); r.className = 'ufv-detail-row';
            const a = document.createElement('a');
            a.href = 'javascript:void 0'; a.className = 'ufv-ligand-link';
            a.textContent = 'Find Similar Motifs ↗';
            a.title = 'Search the PDB for structures with a similar 3-D arrangement of these pocket residues (RCSB Structure Motif Search)';
            a.addEventListener('click', (e) => {
                e.preventDefault(); a.textContent = 'searching…';
                openPocketSearch(pocket).then(ok => { a.textContent = ok ? 'Find Similar Motifs ↗' : 'motif search unavailable'; }).catch(() => { a.textContent = 'motif search unavailable'; });
            });
            r.appendChild(a); box.appendChild(r);
        }
        sec.append(toggle, box);
        body.appendChild(sec);
    }

    function renderLigandPanel(lig) {
        const s = UFVState.state;
        const body = byId('ufv-details-body');
        body.textContent = '';
        byId('ufv-sphere-toggle')?.classList.remove('ufv-hidden'); // restore (a partner panel may have hidden it)
        // Title: the CCD code. When the structure has MULTIPLE copies of this same ligand, add ‹ › arrows to
        // rotate between the copies (each click re-focuses + re-zooms that copy) with an "n/N" counter, so the
        // user can step through e.g. all 4 HEM groups while staying in the zoomed-in pocket view.
        const titleEl = byId('ufv-details-title');
        titleEl.textContent = '';
        const copies = (s.ligands || []).filter(l => l.resn === lig.resn);
        const idx = copies.findIndex(l => l.resi === lig.resi && l.chain === lig.chain);
        if (copies.length > 1 && idx >= 0) {
            const nav = document.createElement('span'); nav.className = 'ufv-lig-nav';
            const mkArrow = (glyph, title, delta) => {
                const b = document.createElement('button');
                b.type = 'button'; b.className = 'ufv-lig-nav-btn'; b.innerHTML = glyph; b.title = title;
                b.addEventListener('click', e => { e.stopPropagation(); onLigandClick(copies[(idx + delta + copies.length) % copies.length]); });
                return b;
            };
            const name = document.createElement('span'); name.className = 'ufv-lig-nav-name'; name.textContent = lig.resn;
            const count = document.createElement('span'); count.className = 'ufv-lig-nav-count'; count.textContent = `${idx + 1}/${copies.length}`;
            nav.append(mkArrow('&#8249;', `Previous ${lig.resn} copy`, -1), name, count, mkArrow('&#8250;', `Next ${lig.resn} copy`, +1));
            titleEl.appendChild(nav);
        } else {
            titleEl.textContent = lig.resn;
        }

        // Ligand label + boxed Nearby slider — same live-distance control as the residue panel,
        // re-running the ligand focus (which also re-hides the surrounding ligands) at the new radius.
        // (The CCD name is already the panel title, so no separate "Ligand" row.)
        const ligRefocus = () => {
            s.nearbyResidues = StructureViewer.focusLigand(lig.resn, lig.resi, lig.chain,
                { showOtherSpheres: _showOtherSpheres, rezoom: false, annotatedResidues: buildAnnotationMap() }) || s.nearbyResidues;
        };
        body.appendChild(makeNearbyBox(lig.chain, ligRefocus));

        // AlphaFill transplant evidence (already on the ligand — shown synchronously): donor identity,
        // clash score with a low/moderate/high band, and the donor PDB. Clash explains poses that pass
        // through the protein (unrefined superposition) — shown for information only, nothing is filtered.
        if (lig.identity != null || lig.clash != null) {
            const aff = document.createElement('div');
            aff.className = 'ufv-ligand-info';
            const affRow = (label, valueNode, cls) => {
                const r = document.createElement('div');
                r.className = 'ufv-detail-row';
                const l = document.createElement('span'); l.className = 'ufv-detail-lbl'; l.textContent = label;
                const v = document.createElement('span'); v.className = 'ufv-detail-val' + (cls ? ' ' + cls : '');
                if (typeof valueNode === 'string') v.textContent = valueNode; else v.appendChild(valueNode);
                r.append(l, v); aff.appendChild(r);
            };
            if (lig.identity != null) {
                affRow('Donor identity', `${Math.round(lig.identity * 100)}%` + (lig.donorPdb ? ` (from ${lig.donorPdb})` : ''));
            }
            if (lig.clash != null) {
                const band = lig.clash <= 0.5 ? 'low' : lig.clash <= 0.9 ? 'moderate' : 'high';
                affRow('Clash', `${lig.clash.toFixed(2)} (${band})`, 'ufv-clash-' + band);
            }
            body.appendChild(aff);
        }

        // Chemistry IDENTITY rows (Name/Formula/MW/…) — populated async from the RCSB CCD.
        const info = document.createElement('div');
        info.className = 'ufv-ligand-info';
        info.innerHTML = `<div class="ufv-detail-row"><span class="ufv-detail-lbl">Loading chemistry…</span></div>`;
        body.appendChild(info);

        // Pocket evidence: BELOW the ligand identity, ABOVE the external links/similarity.
        renderPocketEvidence(body, s.nearbyResidues, lig);

        // Divider + heading for the external compound links (PubChem / DrugBank + 2D/3D similarity).
        // Kept distinct from the pocket-similarity search above (in renderPocketEvidence).
        const simHdr = document.createElement('div');
        simHdr.className = 'ufv-ligand-section-hdr';
        simHdr.textContent = 'Chemical info';
        body.appendChild(simHdr);

        // External links — one UNNAMED row: PubChem · DrugBank · 2D similarity · 3D similarity (filled async).
        const linksRow = document.createElement('div');
        linksRow.className = 'ufv-ligand-info';
        body.appendChild(linksRow);

        const simSection = document.createElement('div'); // populated async by renderSimilarLigands (in-structure)
        simSection.className = 'ufv-am-section ufv-hidden';
        body.appendChild(simSection);

        byId('ufv-details').classList.add('show');
        renderLegend(getColorMode());

        UFVApi.getLigandInfo(lig.resn).then(meta => {
            // Bail if the user moved on to a different ligand/residue meanwhile.
            if (UFVState.state.selectedLigand !== lig) return;
            info.textContent = '';
            renderSimilarLigands(lig, meta, simSection);
            const addRow = (label, valueNode) => {
                const r = document.createElement('div');
                r.className = 'ufv-detail-row';
                const l = document.createElement('span'); l.className = 'ufv-detail-lbl'; l.textContent = label;
                const v = document.createElement('span'); v.className = 'ufv-detail-val';
                if (typeof valueNode === 'string') v.textContent = valueNode; else v.appendChild(valueNode);
                r.append(l, v); info.appendChild(r);
            };
            if (!meta || (!meta.name && !meta.smiles && !meta.inchikey)) {
                addRow('Name', `${lig.resn} (no chemical record found)`);
                return;
            }
            const extLink = (href, text) => {
                const a = document.createElement('a');
                a.href = href; a.target = '_blank'; a.rel = 'noopener noreferrer';
                a.className = 'ufv-ligand-link'; a.textContent = text;
                return a;
            };
            if (meta.name) addRow('Name', meta.name);
            // Chemical descriptors (M6): formula, then molecular weight + H-bond donor/acceptor counts (RCSB).
            if (meta.formula) addRow('Formula', meta.formula);
            if (meta.weight != null) addRow('Mol. weight', `${meta.weight.toFixed(1)} Da`);
            if (meta.hbondDonors != null || meta.hbondAcceptors != null) {
                addRow('H-bond D/A', `${meta.hbondDonors ?? '?'} / ${meta.hbondAcceptors ?? '?'}`);
            }
            if (meta.smiles) addRow('SMILES', copyChip(meta.smiles));
            if (meta.inchikey) addRow('InChIKey', copyChip(meta.inchikey));
            // External links — ONE UNNAMED row of hyperlinks (no "Links"/"Similar" labels):
            // PubChem · DrugBank · 2D similarity · 3D similarity.
            const pubchem = UFVApi.pubchemUrl?.(meta.inchikey);
            const linkEls = [];
            if (pubchem) linkEls.push(extLink(pubchem, 'PubChem'));
            if (meta.drugbank) linkEls.push(extLink(`https://go.drugbank.com/drugs/${meta.drugbank}`, 'DrugBank'));
            // 2D/3D similarity render IMMEDIATELY keyed by InChIKey (PubChem resolves the key); the href is
            // upgraded to the more precise CID-based query once the CID resolves. (Previously these only
            // appeared inside the CID promise, so they were missing whenever that lookup was slow/failed.)
            let sim2d = null, sim3d = null;
            if (meta.inchikey) {
                sim2d = extLink(UFVApi.pubchemSimilarity2dByKey(meta.inchikey), '2D similarity');
                sim3d = extLink(UFVApi.pubchemSimilarity3dByKey(meta.inchikey), '3D similarity');
                linkEls.push(sim2d, sim3d);
            }
            const renderLinks = () => {
                linksRow.textContent = '';
                if (!linkEls.length) return;
                const row = document.createElement('div'); row.className = 'ufv-detail-row ufv-ligand-links';
                linkEls.forEach((a, i) => { if (i) row.appendChild(document.createTextNode(' · ')); row.appendChild(a); });
                linksRow.appendChild(row);
            };
            renderLinks(); // PubChem / DrugBank / 2D / 3D all immediately
            if (meta.inchikey && UFVApi.pubchemCid) {
                UFVApi.pubchemCid(meta.inchikey).then(cid => {
                    if (UFVState.state.selectedLigand !== lig || !cid) return;
                    if (sim2d) sim2d.href = UFVApi.pubchemSimilarity2dUrl(cid); // upgrade to CID-based similarity
                    if (sim3d) sim3d.href = UFVApi.pubchemSimilarity3dUrl(cid);
                });
            }
        });
    }

    function _esc(s) { return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;'); }
    function residueSummary(pos, data = {}) {
        const s = UFVState.state;
        const ptms = s.ptms.filter(p => p.position === pos).map(p => _esc(p.category));
        const vars = s.variants.filter(v => v.position === pos).map(v =>
            `<strong>${_esc(v.wildType)}${v.position}${_esc(v.mutant)}</strong> ${_esc(v.consequence)}`);
        if (data.category) ptms.unshift(_esc(data.category));
        if (data.variants) vars.unshift(...data.variants.map(v =>
            `<strong>${_esc(v.wildType)}${v.position}${_esc(v.mutant)}</strong>`));
        return [...ptms, ...vars].slice(0, 4).join(' | ') || 'No PTM or variant annotation';
    }

    // Per-residue stick/chip colours for the focus view. By request these reflect DISEASE-ASSOCIATED
    // VARIANTS ONLY (pathogenicity colour), not PTMs or functional sites — so the "side chains" legend reads
    // as variant pathogenicity in the pocket, without PTM/site colours leaking in.
    function buildAnnotationMap() {
        const s = UFVState.state;
        const map = new Map();
        s.variants.forEach(v => map.set(v.position, { color: v.consequenceColor }));
        return map;
    }

    // Boxed Nearby section with a live distance slider — replaces the old Position|Nearby two-cell
    // grid in the residue detail panel. Slider adjusts StructureViewer.nearbyDistance on the fly
    // and immediately re-focuses so the 3-D sticks + list update together.
    // refocus(): optional — re-run the focus at the new distance (ligand panels pass a focusLigand
    // version). Defaults to re-focusing the selected residue.
    function makeNearbyBox(chain, refocus) {
        const s = UFVState.state;
        const dist = StructureViewer.nearbyDistance || 5;
        const box = document.createElement('div');
        box.className = 'ufv-nearby-box';
        const hdr = document.createElement('div');
        hdr.className = 'ufv-nearby-hdr';
        const lbl = document.createElement('span');
        lbl.className = 'ufv-detail-lbl';
        lbl.textContent = 'Nearby';
        const slider = document.createElement('input');
        slider.type = 'range'; slider.className = 'ufv-nearby-slider';
        slider.min = 3; slider.max = 15; slider.step = 1; slider.value = dist;
        const distVal = document.createElement('span');
        distVal.className = 'ufv-nearby-dist-val';
        distVal.textContent = dist + ' Å';
        hdr.append(lbl, slider, distVal);
        const listEl = document.createElement('span');
        listEl.className = 'ufv-detail-val ufv-nearby-val';
        const rebuildList = () => {
            listEl.textContent = '';
            // The focused residue itself is in nearbyResidues (so its sticks render) but shouldn't be
            // listed as its own neighbour — exclude it from the "Nearby" list.
            const sorted = Array.from(s.nearbyResidues).filter(p => p !== s.selectedResidue).sort((a, b) => a - b);
            const ann = buildAnnotationMap();
            // multichain: partner-chain neighbours that carry an annotation (disease variant)
            const partnerList = (StructureViewer._nearbyPartners || [])
                .map(pn => { const p = _partnerAnnotations.find(x => x.chainId === pn.chain && x.byResi.has(pn.resi)); return p ? { chain: pn.chain, resi: pn.resi, partner: p, info: p.byResi.get(pn.resi) } : null; })
                .filter(Boolean);
            if (!sorted.length && !partnerList.length) { listEl.textContent = '—'; return; }
            sorted.forEach((p, i) => {
                const span = document.createElement('span');
                span.className = 'ufv-nearby-res';
                span.textContent = (s.sequence?.[p - 1] || '') + p;
                const c = ann.get(p)?.color; if (c) span.style.color = c;
                span.title = `Focus residue ${p}`;
                span.addEventListener('click', () => onClick({ position: p }, null, chain ?? null));
                listEl.appendChild(span);
                if (i < sorted.length - 1) listEl.appendChild(document.createTextNode(', '));
            });
            if (partnerList.length) {
                // Partner-chain neighbours go on their OWN line (they belong to a different protein).
                const line = document.createElement('span');
                line.className = 'ufv-nearby-partner-line';
                const lbl = document.createElement('span'); lbl.className = 'ufv-nearby-chain-lbl'; lbl.textContent = 'other chains: ';
                line.appendChild(lbl);
                appendPartnerNeighbours(line, partnerList);
                listEl.appendChild(line);
            }
        };
        rebuildList();
        slider.addEventListener('input', () => {
            const v = Number(slider.value);
            distVal.textContent = v + ' Å';
            StructureViewer.nearbyDistance = v;
            if (refocus) refocus();
            else if (s.selectedResidue != null) {
                s.nearbyResidues = StructureViewer.focusResidue(s.selectedResidue, s.selectedChain,
                    { annotatedResidues: buildAnnotationMap() },
                    { showOtherSpheres: _showOtherSpheres, rezoom: false }) || s.nearbyResidues;
            }
            rebuildList(); renderSequence();
        });
        slider.addEventListener('change', () => {
            UFVState.saveSettings({ nearbyDistance: Number(slider.value) });
            StructureViewer.rezoomFocus?.(); // re-frame the pocket to the new radius when the user releases the slider
        });
        box.append(hdr, listEl);
        return box;
    }

    // Build a "Nearby" detail-cell whose residue numbers are clickable (clicking re-focuses that
    // residue, in any window). Used by both the residue and ligand detail panels.
    function makeNearbyCell(chain) {
        const s = UFVState.state;
        const annotations = buildAnnotationMap();
        const cell = document.createElement('div');
        cell.className = 'ufv-detail-cell';
        const lbl = document.createElement('span'); lbl.className = 'ufv-detail-lbl'; lbl.textContent = 'Nearby';
        const val = document.createElement('span'); val.className = 'ufv-detail-val ufv-nearby-val';
        const sorted = Array.from(s.nearbyResidues).filter(p => p !== s.selectedResidue).sort((a, b) => a - b);
        if (!sorted.length) val.textContent = '—';
        sorted.forEach((p, i) => {
            const span = document.createElement('span');
            span.className = 'ufv-nearby-res';
            span.textContent = (s.sequence?.[p - 1] || '') + p;
            const c = annotations.get(p)?.color;
            if (c) span.style.color = c;
            span.title = `Focus residue ${p}`;
            span.addEventListener('click', () => onClick({ position: p }, null, chain ?? null));
            val.appendChild(span);
            if (i < sorted.length - 1) val.appendChild(document.createTextNode(', '));
        });
        cell.append(lbl, val);
        return cell;
    }

    // opts.zoomPos (+ opts.chain) — append a magnifier that focuses that residue (like the ligand list),
    // so any single-position annotation row can be zoomed to directly from the sidebar.
    const rangeList = (a, b) => { const out = []; for (let p = a; p <= b; p++) out.push(p); return out; };
    function makeFilterItem(label, color, count, checked, onChange, value, opts) {
        const el = document.createElement('label');
        el.className = 'ufv-filter-item';
        el.innerHTML = `<input type="checkbox"><span class="ufv-dot"></span><span class="ufv-filter-label"></span><span class="ufv-filter-count"></span>`;
        const input = el.querySelector('input');
        input.checked = checked;
        input.dataset.ufvVal = value != null ? value : label; // used by syncVariantFilterChecks
        input.addEventListener('change', e => onChange(e.target.checked));
        el.querySelector('.ufv-dot').style.backgroundColor = color || '#888';
        el.querySelector('.ufv-filter-label').textContent = label;
        el.querySelector('.ufv-filter-count').textContent = count;
        if (opts && opts.zoomPos != null) {
            // Grey out + disable rows whose residue isn't resolved in the loaded structure (can't be drawn).
            if (StructureViewer.isResolved && !StructureViewer.isResolved(opts.zoomPos, opts.chain ?? null)) {
                el.classList.add('ufv-unresolved');
                input.disabled = true;
                el.title = `Residue ${opts.zoomPos} is not resolved in this structure`;
            } else if (opts.zoomRange && opts.zoomRange[1] && opts.zoomRange[1] !== opts.zoomRange[0]) {
                // Range feature (e.g. a domain): the magnifier frames the WHOLE range, not just its start.
                const [a, b] = opts.zoomRange;
                el.appendChild(makeZoomBtn(`Zoom to ${a}–${b}`, () => StructureViewer.frameResidues?.(rangeList(a, b), opts.chain ?? null)));
            } else {
                el.appendChild(makeZoomBtn(`Zoom to ${opts.zoomPos}`, () => onClick({ position: opts.zoomPos }, null, opts.chain ?? null)));
            }
        }
        return el;
    }

    // Disease filter row that expands to list its individual variant residues (each with a zoom icon
    // that focuses that residue). The parent checkbox toggles the disease filter (additive model).
    function makeDiseaseFilter(name, color, count) {
        const s = UFVState.state;
        const box = document.createElement('div');
        box.className = 'ufv-filter-group';
        const children = document.createElement('div');
        children.className = 'ufv-filter-children ufv-collapsed';
        const seen = new Set();
        s.variants
            .filter(v => (v.diseases || []).includes(name))
            .sort((a, b) => a.position - b.position)
            .forEach(v => {
                const key = variantKey(v);
                if (seen.has(key)) return; seen.add(key);
                const rowEl = document.createElement('label');
                rowEl.className = 'ufv-filter-item ufv-disease-var';
                // Individual selection: a checkbox hides/shows just this variant's sphere.
                const cb = document.createElement('input');
                cb.type = 'checkbox'; cb.className = 'ufv-lig-eye';
                cb.checked = !_hiddenVariantKeys.has(key);
                cb.addEventListener('change', () => {
                    // Force-show on check (renders even if the parent disease axis is off); hide on uncheck.
                    if (cb.checked) { _hiddenVariantKeys.delete(key); _forcedVariantKeys.add(key); }
                    else { _forcedVariantKeys.delete(key); _hiddenVariantKeys.add(key); }
                    applyMode();
                });
                const lbl = document.createElement('span');
                lbl.className = 'ufv-filter-label';
                lbl.textContent = `${v.wildType || ''}${v.position}${v.mutant || ''}`;
                if (v.consequenceColor) lbl.style.color = v.consequenceColor;
                // Zoom only via the icon (not the row), so it doesn't fight the checkbox.
                const focus = () => onClick({ position: v.position }, null, s.selectedChain ?? null);
                rowEl.append(cb, lbl, makeZoomBtn(`Zoom to ${v.position}`, focus));
                children.appendChild(rowEl);
            });
        const top = makeFilterItem(name, color, count, true,
            checked => toggleVariantFilter('disease', name, checked), name);
        const chevron = document.createElement('button');
        chevron.type = 'button'; chevron.className = 'ufv-group-chevron'; chevron.innerHTML = '&#9654;';
        chevron.addEventListener('click', e => {
            e.preventDefault(); e.stopPropagation();
            const c = children.classList.toggle('ufv-collapsed');
            chevron.innerHTML = c ? '&#9654;' : '&#9660;';
        });
        top.appendChild(chevron);
        box.append(top, children);
        return box;
    }

    // reRender — how to repaint after a child/parent toggle. Defaults to the fast PTM path (PTM window);
    // other windows pass `reapply`/`applyMode` so expandable PTM categories work there too.
    function makeExpandableFilter(label, color, items, checked, onChange, itemLabel, reRender = applyPTMMode) {
        const box = document.createElement('div');
        box.className = 'ufv-filter-group';

        // Create children first so the parent onChange closure can reference it
        const children = document.createElement('div');
        children.className = 'ufv-filter-children ufv-collapsed';
        items.slice().sort((a, b) => a.position - b.position).forEach(item => {
            const child = makeFilterItem(itemLabel(item), color, item.position, item.visible !== false, val => {
                item.visible = val;
                reRender();
            }, undefined, { zoomPos: item.position });
            children.appendChild(child);
        });

        const top = makeFilterItem(label, color, items.length, checked, (isChecked) => {
            // Sync all child checkbox DOM elements
            children.querySelectorAll('input[type="checkbox"]').forEach(cb => { cb.checked = isChecked; });
            onChange(isChecked);
        });

        // Expand/collapse chevron button
        const chevron = document.createElement('button');
        chevron.type = 'button';
        chevron.className = 'ufv-group-chevron';
        chevron.innerHTML = '&#9654;';
        chevron.addEventListener('click', e => {
            e.preventDefault();
            e.stopPropagation();
            const collapsed = children.classList.toggle('ufv-collapsed');
            chevron.innerHTML = collapsed ? '&#9654;' : '&#9660;';
        });
        top.appendChild(chevron);

        box.append(top, children);
        return box;
    }

    function row(label, value, color) {
        const div = document.createElement('div');
        div.className = 'ufv-detail-row';
        const l = document.createElement('span');
        l.className = 'ufv-detail-lbl';
        l.textContent = label;
        const v = document.createElement('span');
        v.className = 'ufv-detail-val';
        v.textContent = value || '';
        if (color) v.style.color = color;
        div.append(l, v);
        return div;
    }

    function ptmSetAll(select) {
        Object.values(UFVState.state.ptmGroups).forEach(g => {
            g.visible = select;
            g.items.forEach(i => i.visible = select);
        });
        document.querySelectorAll('#ufv-ptm-list input[type="checkbox"]').forEach(cb => cb.checked = select);
        applyPTMMode();
    }

    // ---- Additive variant filter model -----------------------------------------------------------
    // Disease / Consequence / Provenance are three axes over the SAME variants. Checking a value
    // shows that value's variants; to keep the intersection in filterVariants from hiding them, we
    // also auto-select (propagate) the values those variants carry on the OTHER two axes. So checking
    // "DRVT" ticks DRVT's consequence + provenance values too, and the DRVT variants appear. Checking
    // a second value adds its variants on top; unchecking removes only that one value. None ⇒ nothing.
    function variantFilterSet(dim) {
        const s = UFVState.state;
        return dim === 'consequence' ? s.activeConsequences : dim === 'provenance' ? s.activeProvenances : s.activeDiseases;
    }
    function variantMatchesValue(v, dim, value) {
        if (dim === 'disease') return (v.diseases || []).includes(value);
        return v[dim] === value;
    }
    // Propagate a matched variant's values to the OTHER two axes (never the one being toggled — else
    // checking one disease would tick every disease that co-occurs on its multi-disease variants).
    function propagateVariant(v, skipDim) {
        const s = UFVState.state;
        if (skipDim !== 'consequence' && s.activeConsequences && DataProcessor.CONSEQUENCE_CATEGORIES[v.consequence]) s.activeConsequences.add(v.consequence);
        if (skipDim !== 'provenance' && s.activeProvenances && DataProcessor.PROVENANCE_CATEGORIES[v.provenance]) s.activeProvenances.add(v.provenance);
        if (skipDim !== 'disease' && s.activeDiseases) {
            (v.diseases || []).filter(d => d !== '__other__').forEach(d => s.activeDiseases.add(d));
        }
    }
    // Re-tick every variant-filter checkbox from the active sets (no DOM rebuild ⇒ scroll/collapse kept).
    function syncVariantFilterChecks() {
        const s = UFVState.state;
        const apply = (id, set) => {
            if (!set) return;
            document.querySelectorAll(`#${id} input[type="checkbox"]`).forEach(cb => { cb.checked = set.has(cb.dataset.ufvVal); });
        };
        apply('ufv-dis-list', s.activeDiseases);
        apply('ufv-cons-list', s.activeConsequences);
        apply('ufv-prov-list', s.activeProvenances);
    }
    function toggleVariantFilter(dim, value, checked) {
        const set = variantFilterSet(dim);
        if (!set) return;
        if (!checked) { set.delete(value); syncVariantFilterChecks(); applyMode(); return; }
        set.add(value);
        UFVState.state.variants.forEach(v => { if (variantMatchesValue(v, dim, value)) propagateVariant(v, dim); });
        syncVariantFilterChecks();
        applyMode();
    }

    function varSectionSetAll(section, select) {
        const s = UFVState.state;
        const target = variantFilterSet(section);
        if (!target) return;
        const id = section === 'consequence' ? 'ufv-cons-list' : section === 'provenance' ? 'ufv-prov-list' : 'ufv-dis-list';
        if (!select) {
            // None on any axis clears ALL three variant axes: the axes are coupled, so emptying one
            // would hide everything anyway — keeping the other boxes ticked would be misleading.
            s.activeDiseases?.clear();
            s.activeConsequences?.clear();
            s.activeProvenances?.clear();
        } else {
            // Select every value shown in this section, then propagate so the other axes admit them.
            document.querySelectorAll(`#${id} input[type="checkbox"]`).forEach(cb => target.add(cb.dataset.ufvVal));
            s.variants.forEach(v => { if (section === 'disease' || v[section]) propagateVariant(v, section); });
        }
        // All/None override the per-variant force-show/hide overlays: None hides everything, All shows it.
        _forcedVariantKeys.clear();
        _hiddenVariantKeys.clear();
        syncVariantFilterChecks();
        applyMode();
    }

    async function copySelection() {
        const text = UFVExport.formatSelection(UFVState.state.displayedPositions, UFVState.state.settings.copyFormat, UFVState.selectedStructure());
        await UFVExport.copyText(text);
        flashButton('ufv-btn-copy', 'Copied');
    }

    function exportPdb() {
        const st = UFVState.selectedStructure();
        if (!StructureViewer.currentPdbText || !st) return;
        if (StructureViewer.currentFormat === 'mmcif') {
            // B-factor rewriting requires PDB format; CIF structures can still be downloaded as-is
            UFVExport.downloadText(`${UFVState.state.uniprotId}_${st.id || 'structure'}.cif`, StructureViewer.currentPdbText, 'chemical/x-cif');
            return;
        }
        const s = UFVState.state;
        const mode = getColorMode();
        let colorContext;
        if (mode === 'alphaMissense') colorContext = s.analysis.alphaMissense;
        else if (mode === 'hotspots') colorContext = s.analysis.hotspots;
        else if (mode === 'distantContacts') colorContext = s.analysis.distantContacts;
        else if (mode === 'residueBurden') colorContext = s.analysis.residueBurden;
        else if (mode === 'prism') colorContext = s.analysis.prism?.byPos;
        const text = UFVExport.rewritePdbBeta(StructureViewer.currentPdbText, s.displayedPositions, st, mode, colorContext);
        UFVExport.downloadText(`${s.uniprotId}_${st.id || 'structure'}_${mode}.pdb`, text);
    }

    async function exportCsv(withProtVar = false) {
        const s = UFVState.state;
        if (!s.sequence) return;
        // Compute the constraint-pocket values on demand for the loaded structure so the CSV
        // always carries them (the analysis is otherwise only run when its coloring mode is used).
        if (!s.analysis.prism && s.amMap?.size && StructureViewer.viewer) {
            try { await ensurePocketAnalysis(); } catch (_) {}
        }
        // Per-residue ligand contacts (CCD codes within 5 Å) for the loaded structure.
        s.analysis.ligandContacts = (s.ligands?.length && StructureViewer.ligandContactsByResidue)
            ? StructureViewer.ligandContactsByResidue(5) : null;
        // Opt-in: one ProtVar /score request per residue (no bulk endpoint), so it's slow — show progress.
        let protvarByPos = null;
        if (withProtVar) {
            showLoading('ProtVar predictions: 0 / ' + s.sequence.length + ' residues…');
            try {
                protvarByPos = await UFVApi.fetchProtVarAll(s.uniprotId, s.sequence.length,
                    (done, total) => showLoading(`ProtVar predictions: ${done} / ${total} residues…`));
            } catch (_) { protvarByPos = null; }
            byId('ufv-loading').classList.add('hidden');
        }
        const text = UFVExport.buildResidueMatrix(s.sequence, s.ptms, s.ptmGroups || {}, s.variants, s.amMap, s.analysis, UFVState.selectedStructure(), s.sites, s.mutagenesis, protvarByPos);
        UFVExport.downloadText(`${s.uniprotId}_residue_annotations${withProtVar ? '_protvar' : ''}.csv`, text, 'text/csv');
    }

    // Download a self-contained PyMOL/VMD script that recreates the current 3-D view (cartoon
    // colours, annotation spheres, ligands) so the user can keep working from it in that program.
    function exportSession(format) {
        const s = UFVState.state;
        const scene = StructureViewer.getSceneState?.();
        if (!scene) return;
        const st = UFVState.selectedStructure();
        const objName = (st?.id || st?.pdbId || 'structure').toString().replace(/[^A-Za-z0-9_-]/g, '_');
        const stem = `${s.uniprotId}_${objName}`;
        if (format === 'vmd') {
            // VMD scripts ARE Tcl; the .vmd extension makes the file load directly in VMD (double-click /
            // `vmd file.vmd`), which a generic .tcl doesn't signal. Content is unchanged (valid Tcl).
            UFVExport.downloadText(`${stem}.vmd`, UFVExport.buildVmdSession(scene, objName), 'text/plain');
        } else {
            UFVExport.downloadText(`${stem}.pml`, UFVExport.buildPymolSession(scene, objName), 'text/plain');
        }
    }

    async function cycleStructure(delta) {
        const s = UFVState.state;
        if (!s.structures.length) return;
        s.selectedStructureIndex = (s.selectedStructureIndex + delta + s.structures.length) % s.structures.length;
        await loadSelectedStructure();
    }

    function cycleTheme() {
        const s = UFVState.state;
        // Binary light↔dark toggle (was a 3-state auto/light/dark cycle: leaving dark went to 'auto',
        // which on a dark-mode system still LOOKS dark — so it took two clicks to visibly reach light).
        // Resolve the current visible appearance, then flip to the opposite explicit theme.
        const currentlyDark = s.theme === 'dark' || (s.theme !== 'light' && window.matchMedia('(prefers-color-scheme: dark)').matches);
        s.theme = currentlyDark ? 'light' : 'dark';
        modalEl.classList.remove('ufv-light', 'ufv-dark');
        modalEl.classList.add(`ufv-${s.theme}`);
        const dark = s.theme === 'dark';
        StructureViewer.viewer?.setBackgroundColor(dark ? '#0c111b' : '#f0f2f5');
        StructureViewer.viewer?.render();
    }

    function syncSettingsControls() {
        const s = UFVState.state;
        const showExploratory = s.settings.showExploratoryAlgorithms ?? true;
        // Show/hide the 4 exploratory algorithm color modes in the dropdown
        const exploratoryModes = ['hotspots', 'distantContacts', 'residueBurden', 'prism'];
        exploratoryModes.forEach(m => {
            byId('ufv-cm-drop')?.querySelector(`.ufv-cm-opt[data-value="${m}"]`)
                ?.classList.toggle('ufv-hidden', !showExploratory);
        });
        let mode = s.settings.coloringMode;
        // If the saved mode is an exploratory one and exploratory is now off, fall back to default
        if ((!showExploratory && exploratoryModes.includes(mode)) || mode === 'topos') mode = 'default';
        setColorMode(mode);
        byId('ufv-sens-wrap').classList.add('ufv-hidden');
        StructureViewer.nearbyDistance = s.settings.nearbyDistance || 5;
        const fs = s.settings.fontScale || 0;
        modalEl?.classList.toggle('ufv-fs-sm', fs < 0);
        modalEl?.classList.toggle('ufv-fs-lg', fs > 0);
    }

    function switchTab(_name) { /* tabs removed — settings moved to options page */ }

    function toggleCollapsible(bodyId, toggleId) {
        const body = byId(bodyId);
        const collapsed = body.classList.toggle('ufv-collapsed');
        const ch = byId(toggleId).querySelector('.ufv-collapsible-chevron');
        if (ch) ch.innerHTML = collapsed ? '&#9654;' : '&#9660;';
    }

    function showLoading(text) {
        byId('ufv-loading').classList.remove('hidden');
        byId('ufv-loading-text').textContent = text;
    }

    function showError(text) {
        byId('ufv-loading').classList.remove('hidden');
        byId('ufv-loading-text').textContent = text;
    }

    function flashButton(id, text) {
        const btn = byId(id);
        const old = btn.innerHTML;
        btn.textContent = text;
        btn.classList.add('ufv-copied');
        setTimeout(() => {
            btn.innerHTML = old;
            btn.classList.remove('ufv-copied');
        }, 1200);
    }

    function close() {
        _hostScrollLock = false;
        if (overlayEl) overlayEl.style.display = 'none';
        document.documentElement.style.overflow = '';
        document.body.style.overflow = '';
    }

    // Re-render the 3Dmol canvas whenever the user returns to this browser tab.
    // WebGL contexts can be discarded by the browser while the tab is backgrounded;
    // resize() re-calculates canvas dimensions and render() redraws the scene.
    // The camera is deliberately NOT reset here — if the user was zoomed into a residue,
    // switching windows and back should preserve that view (a full WebGL context loss is
    // handled separately by the webglcontextrestored handler, which rebuilds the scene).
    document.addEventListener('visibilitychange', () => {
        if (document.hidden || !overlayEl || overlayEl.style.display === 'none') return;
        const v = StructureViewer.viewer;
        if (!v) return;
        try {
            v.resize();
            v.render();
        } catch (_) { /* WebGL context may be in the process of being restored */ }
    });

    /**
     * Silent background prefetch — called by injector.js immediately on page load
     * so annotation data is ready before the user opens the modal.
     */
    async function prefetchData() {
        const s = UFVState.state;
        if (!s.uniprotId || s.annotationsLoaded || s.loadingPromise) return;
        const requestedId = s.uniprotId;
        const thisPromise = (async () => {
            try {
                await UFVState.loadSettings();
                const data = await UFVApi.loadFeatureData(requestedId);
                if (UFVState.state.uniprotId !== requestedId) return;
                Object.assign(s, data);
                s.ptmGroups = DataProcessor.groupPTMsByCategory(s.ptms);
                s.activeConsequences = new Set(Object.keys(DataProcessor.getConsequenceSummary(s.variants)));
                s.activeProvenances = new Set(Object.keys(DataProcessor.getProvenanceSummary(s.variants)));
                s.analysis.alphaMissense = UFVAnalysis.aggregateAlphaMissense(s.variants, s.amMap);
                DataProcessor.computeDiseaseColors(s.variants);
            const _ds = DataProcessor.getDiseaseSummary(s.variants);
            _diseaseColorMap = new Map(Object.entries(_ds).map(([n, m]) => [n, m.color]));
                s.annotationsLoaded = true;
            } catch (_) {
                // Fail silently — open() will retry via loadAnnotations()
            } finally {
                // Only clear the shared handle if it is still ours: a newer protein's
                // prefetch may have replaced it after a navigation, and nulling that one
                // would make open() skip awaiting the in-flight load.
                if (s.loadingPromise === thisPromise) s.loadingPromise = null;
            }
        })();
        s.loadingPromise = thisPromise;
        return thisPromise;
    }

    function byId(id) {
        return document.getElementById(id);
    }

    // ProtNLM: UniProt's AI model names proteins predicted from sequence. Render the predicted name
    // (flagged as AI) into the header banner, or a note when the entry's name is curated/rule-based.
    function renderProtNLM() {
        const banner = byId('ufv-protnlm-banner');
        if (!banner) return;
        const p = UFVState.state.protnlm;
        if (!p || !p.name) {
            banner.innerHTML = `<span style="opacity:.7;">No protein name in the UniProt entry.</span>`;
            return;
        }
        const esc = (str) => String(str).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
        if (p.isAI) {
            let html = `<b style="color:#7c4dff;">ProtNLM (AI-predicted name):</b> ${esc(p.name)}`;
            if (p.source) html += ` <span style="opacity:.6;">· ${esc(p.source)}</span>`;
            if (p.caution) html += `<br><span style="opacity:.6;font-size:11px;">${esc(p.caution)}</span>`;
            banner.innerHTML = html;
        } else {
            const tag = p.reviewed ? 'curated' : 'not AI-generated';
            banner.innerHTML = `<b>Protein name:</b> ${esc(p.name)} <span style="opacity:.6;">(${tag} — no ProtNLM prediction)</span>`;
        }
    }

    return { createButton, open, close, prefetchData };
})();
