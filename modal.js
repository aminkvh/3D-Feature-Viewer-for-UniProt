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
    let _gsgt9SavedScrollTop = 0;
    let _gsgt9LockActive = false;
    // Monotonic token: every open() and structure load captures the current value; if a newer
    // call has since bumped it, the older (superseded) flow bails instead of clobbering the UI.
    // This is what prevents the "variant modal shows the PTM header / nothing ever loads" race
    // when the user navigates between proteins or re-opens the modal mid-load.
    let _openSeq = 0;
    let _loadSeq = 0;

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
                    <div class="ufv-dl-wrap" id="ufv-dl-wrap">
                        <button class="ufv-icon-btn" id="ufv-btn-export-pdb" title="Download">${ICON_DOWNLOAD}</button>
                        <div class="ufv-dl-menu" id="ufv-dl-menu">
                            <button class="ufv-dl-opt" id="ufv-dl-pdb">PDB file</button>
                            <button class="ufv-dl-opt" id="ufv-dl-csv">CSV annotation table</button>
                        </div>
                    </div>
                    <button class="ufv-icon-btn" id="ufv-btn-screenshot" title="Screenshot">${ICON_CAMERA}</button>
                    <button class="ufv-close-btn" id="ufv-close" title="Close">&#10005;</button>
                </div>
            </div>
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
                            <button class="ufv-cm-btn" id="ufv-cm-btn">Default cyan</button>
                            <div class="ufv-cm-drop" id="ufv-cm-drop">
                                <div class="ufv-cm-opt selected" data-value="default">Default cyan</div>
                                <div class="ufv-cm-opt" data-value="plddt">pLDDT confidence</div>
                                <div class="ufv-cm-opt" data-value="bfactor">Experimental B-factor</div>
                                <div class="ufv-cm-opt" data-value="hotspots">3D enrichment hotspot</div>
                                <div class="ufv-cm-opt" data-value="distantContacts">Long-range contact hub</div>
                                <div class="ufv-cm-opt" data-value="alphaMissense">AlphaMissense score</div>
                                <div class="ufv-cm-opt" data-value="residueBurden">Residue burden hotspot</div>
                            </div>
                        </div>
                    </div>
                    <div id="ufv-ptm-panel" class="ufv-filter-scroll">
                        <div class="ufv-panel-hdr"><h3>PTM Types</h3><div class="ufv-panel-actions"><button class="ufv-sm-btn" id="ufv-ptm-all">All</button><button class="ufv-sm-btn" id="ufv-ptm-none">None</button></div></div>
                        <div id="ufv-ptm-list"></div>
                    </div>
                    <div id="ufv-var-panel" class="ufv-filter-scroll ufv-hidden">
                        <div id="ufv-dis-section" class="ufv-hidden"><div class="ufv-section-title"><span>Disease <span class="ufv-section-source">— HumanVar</span></span><div><button class="ufv-section-btn" id="ufv-dis-all">All</button><button class="ufv-section-btn" id="ufv-dis-none">None</button></div></div><div id="ufv-dis-list"></div></div>
                        <div class="ufv-collapsible"><div class="ufv-collapsible-hdr" id="ufv-prov-toggle"><span class="ufv-collapsible-chevron">&#9654;</span><span>Provenance</span><div class="ufv-section-actions"><button class="ufv-section-btn" id="ufv-prov-all">All</button><button class="ufv-section-btn" id="ufv-prov-none">None</button></div></div><div class="ufv-collapsible-body ufv-collapsed" id="ufv-prov-body"><div id="ufv-prov-list"></div></div></div>
                        <div class="ufv-collapsible"><div class="ufv-collapsible-hdr" id="ufv-cons-toggle"><span class="ufv-collapsible-chevron">&#9654;</span><span>Consequence</span><div class="ufv-section-actions"><button class="ufv-section-btn" id="ufv-cons-all">All</button><button class="ufv-section-btn" id="ufv-cons-none">None</button></div></div><div class="ufv-collapsible-body ufv-collapsed" id="ufv-cons-body"><div id="ufv-cons-list"></div></div></div>
                    </div>
                    <div class="ufv-panel-footer"><span class="ufv-count-text" id="ufv-count-text">-</span><button class="ufv-copy-btn" id="ufv-btn-copy">${ICON_COPY} Copy</button></div>
                    <div class="ufv-details" id="ufv-details"><div class="ufv-details-hdr"><h4 id="ufv-details-title">Details</h4><button class="ufv-details-close" id="ufv-details-close">&#10005;</button></div><div class="ufv-details-body" id="ufv-details-body"></div></div>
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

        // Reactively lock UniProt's Gsgt9 custom-scroll container while the modal is
        // open.  This catches every scroll cause (focus changes, click handlers,
        // wheel events that UniProt's capture handler saw) and resets scrollTop
        // immediately so the background page never jumps.
        (function () {
            const gsgt9 = document.querySelector('.Gsgt9');
            if (!gsgt9) return;
            let _locking = false;
            gsgt9.addEventListener('scroll', function () {
                if (!_gsgt9LockActive || _locking) return;
                _locking = true;
                gsgt9.scrollTop = _gsgt9SavedScrollTop;
                _locking = false;
            });
        })();

        bindEvents();
    }

    function bindEvents() {
        byId('ufv-close').addEventListener('click', close);
        // Single overlay-level handler: stops all modal clicks from reaching UniProt's
        // page handlers (which would scroll Gsgt9), and closes dropdowns on outside-click.
        overlayEl.addEventListener('click', e => {
            if (!e.target.closest('#ufv-dl-wrap')) byId('ufv-dl-menu')?.classList.remove('open');
            if (!e.target.closest('#ufv-cs')) byId('ufv-cs')?.classList.remove('open');
            if (!e.target.closest('#ufv-cm')) byId('ufv-cm')?.classList.remove('open');
            if (e.target === overlayEl) close();
            e.stopPropagation();
        });
        document.addEventListener('keydown', e => { if (e.key === 'Escape') close(); });
        byId('ufv-btn-theme').addEventListener('click', cycleTheme);
        byId('ufv-btn-reset').addEventListener('click', () => {
            const s = UFVState.state;
            s.selectedResidue = null;
            s.nearbyResidues = new Set();
            StructureViewer._selectedResi = null;
            byId('ufv-details').classList.remove('show');
            const defaultMode = s.settings.coloringMode || 'default';
            setColorMode(defaultMode);
            StructureViewer.viewer?.zoomTo({}, 600);
            // Defer the heavy cartoon rebuild so the zoom animation starts immediately
            requestAnimationFrame(() => {
                applyMode();
                renderSequence();
            });
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
        byId('ufv-ptm-all').addEventListener('click', () => ptmSetAll(true));
        byId('ufv-ptm-none').addEventListener('click', () => ptmSetAll(false));
        byId('ufv-cons-all').addEventListener('click', () => varSectionSetAll('consequence', true));
        byId('ufv-cons-none').addEventListener('click', () => varSectionSetAll('consequence', false));
        byId('ufv-prov-all').addEventListener('click', () => varSectionSetAll('provenance', true));
        byId('ufv-prov-none').addEventListener('click', () => varSectionSetAll('provenance', false));
        byId('ufv-dis-all').addEventListener('click', () => varSectionSetAll('disease', true));
        byId('ufv-dis-none').addEventListener('click', () => varSectionSetAll('disease', false));
        byId('ufv-prov-toggle').addEventListener('click', e => { if (!e.target.closest('button')) toggleCollapsible('ufv-prov-body', 'ufv-prov-toggle'); });
        byId('ufv-cons-toggle').addEventListener('click', e => { if (!e.target.closest('button')) toggleCollapsible('ufv-cons-body', 'ufv-cons-toggle'); });
        byId('ufv-details-close').addEventListener('click', () => byId('ufv-details').classList.remove('show'));
        byId('ufv-cs-btn').addEventListener('click', e => { e.stopPropagation(); byId('ufv-cs').classList.toggle('open'); });
        byId('ufv-cm-btn').addEventListener('click', e => { e.stopPropagation(); byId('ufv-cm').classList.toggle('open'); });
        byId('ufv-cm-drop').querySelectorAll('.ufv-cm-opt').forEach(opt => {
            opt.addEventListener('click', () => {
                byId('ufv-cm').classList.remove('open');
                setColorMode(opt.dataset.value);
                UFVState.saveSettings({ coloringMode: opt.dataset.value });
                requestAnimationFrame(() => applyMode());
            });
        });
        byId('ufv-structure-prev').addEventListener('click', () => cycleStructure(-1));
        byId('ufv-structure-next').addEventListener('click', () => cycleStructure(1));
    }

    function bindSettings() {
        // Settings moved to extension options page (right-click extension icon → Options)
    }

    async function open(mode) {
        const s = UFVState.state;
        const mySeq = ++_openSeq;
        s.currentMode = mode;
        build();
        await UFVState.loadSettings();
        if (_openSeq !== mySeq) return; // superseded by a newer open()
        syncSettingsControls();
        byId('ufv-id-badge').textContent = s.uniprotId;
        byId('ufv-modal-heading').textContent = mode === 'ptm' ? 'PTM Viewer' : 'Disease & Variants';
        byId('ufv-ptm-panel').classList.toggle('ufv-hidden', mode !== 'ptm');
        byId('ufv-var-panel').classList.toggle('ufv-hidden', mode !== 'variant');
        _gsgt9LockActive = false; // reset so focus-scroll and UniProt handlers can settle
        overlayEl.style.display = 'flex';
        document.documentElement.style.overflow = 'hidden';
        document.body.style.overflow = 'hidden';
        // One animation frame is enough for UniProt's focus-scroll and any other deferred
        // scroll handlers to fire and Gsgt9 to reach its natural resting position before
        // we lock it in place.
        requestAnimationFrame(() => {
            const g = document.querySelector('.Gsgt9');
            if (g) _gsgt9SavedScrollTop = g.scrollTop;
            _gsgt9LockActive = true;
        });
        // Reset any residue focus / detail state from a previous session
        s.selectedResidue = null;
        s.nearbyResidues = new Set();
        StructureViewer._selectedResi = null;
        StructureViewer._inFocusMode = false;
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
        if (!s.loaded) {
            loadStructuresAndShow(); // intentionally not awaited — shows viewer async
        } else {
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
            s.annotationsLoaded = true;
        } catch (err) {
            showError(err.message || 'Unable to load annotations.');
        }
    }

    async function loadStructuresAndShow() {
        const s = UFVState.state;
        const requestedId = s.uniprotId;
        showLoading('Finding structures…');
        try {
            const structures = await UFVApi.getStructures(requestedId, s.sequence.length);
            if (UFVState.state.uniprotId !== requestedId) return;
            s.structures = structures;
            chooseDefaultStructure();
            s.loaded = true;
            renderStructureSelector();
            await loadSelectedStructure();
        } catch (err) {
            showError(err.message || 'Unable to load structures.');
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
                const hotspots = UFVAnalysis.computeHotspots(StructureViewer.viewer, s.variants, structure);
                s.analysis.hotspots = hotspots.merged;
                s.analysis.hotspotsByChain = hotspots.byChain;
                s.analysis.hotspotMethod = hotspots.method;
                const contacts = UFVAnalysis.computeDistantContacts(StructureViewer.viewer, structure, s.variants);
                s.analysis.distantContacts = contacts.merged;
                s.analysis.distantContactsByChain = contacts.byChain;
                s.analysis.residueBurden = UFVAnalysis.computeResidueBurden(s.variants);
                byId('ufv-loading').classList.add('hidden');
            } catch (err) {
                if (_loadSeq === mySeq) showError(err.message || 'Unable to load selected structure.');
                return;
            }
        }
        if (_loadSeq !== mySeq) return; // a newer load owns the viewer now
        StructureViewer.hoverCb = onHover;
        StructureViewer.clickCb = onClick;
        StructureViewer.dblClickCb = () => {
            const s = UFVState.state;
            s.selectedResidue = null;
            s.nearbyResidues = new Set();
            StructureViewer._selectedResi = null;
            byId('ufv-details').classList.remove('show');
            StructureViewer.viewer?.zoomTo({}, 600);
            // Defer heavy rebuild so the zoom animation starts immediately
            requestAnimationFrame(() => {
                applyMode();
                renderSequence();
            });
        };
        updateStructureMeta();
        applyMode();
        // Fold neighbouring partner-protein disease residues into the hotspot test, off the
        // critical path (network fetch) so opening a complex isn't delayed.
        augmentHotspotsWithPartners(structure, requestedId, mySeq);
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
        if (!st || st.source === 'AlphaFold') return 'AlphaFold';
        const chains = st.chainIds?.length > 1 ? st.chainIds.join(',') : (st.chainId || '?');
        const parts = [];
        if (st.method) parts.push(shortMethod(st.method));
        if (st.resolution) parts.push(`${st.resolution}Å`);
        if (st.coverage) parts.push(`${st.coverage}%`);
        const suffix = parts.length ? ` (${parts.join(', ')})` : '';
        return `${st.pdbId}-${chains}${suffix}`;
    }

    function renderStructureSelector() {
        const drop = byId('ufv-cs-drop');
        drop.textContent = '';
        const idx = UFVState.state.selectedStructureIndex;
        UFVState.state.structures.forEach((st, i) => {
            const opt = document.createElement('div');
            opt.className = 'ufv-cs-opt' + (i === idx ? ' selected' : '');
            opt.textContent = stLabel(st);
            opt.addEventListener('click', async () => {
                byId('ufv-cs').classList.remove('open');
                UFVState.state.selectedStructureIndex = i;
                await loadSelectedStructure();
            });
            drop.appendChild(opt);
        });
        const sel = UFVState.state.structures[idx];
        const btn = byId('ufv-cs-btn');
        if (btn && sel) btn.textContent = stLabel(sel);
    }

    function getColorMode() {
        return byId('ufv-cm-drop')?.querySelector('.ufv-cm-opt.selected')?.dataset.value || UFVState.state.settings.coloringMode || 'default';
    }

    function setColorMode(value) {
        const btn = byId('ufv-cm-btn');
        byId('ufv-cm-drop')?.querySelectorAll('.ufv-cm-opt').forEach(opt => {
            opt.classList.toggle('selected', opt.dataset.value === value);
            if (opt.dataset.value === value && btn) btn.textContent = opt.textContent;
        });
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
            const cur = getColorMode();
            if ((cur === 'bfactor' && isAlphaFold) || (cur === 'plddt' && !isAlphaFold)) {
                setColorMode('default');
                UFVState.saveSettings({ coloringMode: 'default' });
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
        const st = UFVState.selectedStructure();
        const seqLen = UFVState.state.sequence?.length || 0;
        if (!st || st.source === 'AlphaFold' || !st.mappedRanges?.length || seqLen === 0) return '';
        const ranges = st.mappedRanges;
        const totalMapped = ranges.reduce((acc, r) => acc + (r.uniprotEnd - r.uniprotStart + 1), 0);
        if (totalMapped / seqLen > 0.95) return ''; // effectively full coverage — no note needed
        if (ranges.length === 1) return ` (mapped: ${ranges[0].uniprotStart}–${ranges[0].uniprotEnd})`;
        const min = Math.min(...ranges.map(r => r.uniprotStart));
        const max = Math.max(...ranges.map(r => r.uniprotEnd));
        return ` (mapped: ${min}–${max}, ${ranges.length} segments)`;
    }

    function applyMode() {
        const s = UFVState.state;
        if (!StructureViewer.viewer) return;
        const mode = getColorMode();
        const filteredVariants = DataProcessor.filterVariants(s.variants, s.activeConsequences, s.activeProvenances, s.activeDiseases);
        // defer=true skips the intermediate render; showPTMs/showVariants will render once
        StructureViewer.applyCartoonColoring(mode, {
            ptms: activePtms(),
            variants: filteredVariants,
            hotspots: s.analysis.hotspots,
            hotspotsByChain: s.analysis.hotspotsByChain,
            distantContacts: s.analysis.distantContacts,
            distantContactsByChain: s.analysis.distantContactsByChain,
            alphaMissense: s.analysis.alphaMissense,
            residueBurden: s.analysis.residueBurden,
        }, true);
        const rangeNote = getMappedRangeNote();
        if (s.currentMode === 'ptm') {
            const n = StructureViewer.showPTMs(s.ptms, s.ptmGroups);
            s.displayedPositions = activePtms().flatMap(p => p.endPosition && p.endPosition !== p.position ? [p.position, p.endPosition] : [p.position]);
            byId('ufv-count-text').textContent = `${n} PTM site${n === 1 ? '' : 's'}${rangeNote}`;
        } else {
            const r = StructureViewer.showVariants(filteredVariants);
            s.displayedPositions = Array.from(new Set(filteredVariants.map(v => v.position)));
            byId('ufv-count-text').textContent = `${r.varCount} variants at ${r.posCount} positions${rangeNote}`;
        }
        renderLegend(mode);
        renderSequence();
    }

    // Fast PTM sphere refresh — skips cartoon rebuild when only sphere visibility changes.
    function applyPTMMode() {
        const s = UFVState.state;
        if (!StructureViewer.viewer) return;
        const n = StructureViewer.refreshPTMDisplay(s.ptms, s.ptmGroups);
        if (n === false) {
            // Focus mode active — need a full rebuild to restore sticks
            applyMode();
            return;
        }
        s.displayedPositions = activePtms().flatMap(p => p.endPosition && p.endPosition !== p.position ? [p.position, p.endPosition] : [p.position]);
        byId('ufv-count-text').textContent = `${n} PTM site${n === 1 ? '' : 's'}${getMappedRangeNote()}`;
        renderLegend(getColorMode());
        renderSequence();
    }

    function activePtms() {
        const s = UFVState.state;
        return s.ptms.filter(p => s.ptmGroups[p.category]?.visible && p.visible !== false);
    }

    function buildFilters() {
        buildPTMFilters();
        buildVariantFilters();
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
    }

    function buildVariantFilters() {
        const s = UFVState.state;
        fillFilterList('ufv-prov-list', DataProcessor.getProvenanceSummary(s.variants), s.activeProvenances, applyMode);
        fillFilterList('ufv-cons-list', DataProcessor.getConsequenceSummary(s.variants), s.activeConsequences, applyMode);
        const ds = DataProcessor.getDiseaseSummary(s.variants);
        const dis = byId('ufv-dis-section');
        const list = byId('ufv-dis-list');
        list.textContent = '';
        const onVariantPage = s.pageContext === 'variant-viewer';
        let diseasesToShow = Object.keys(ds);
        if (!onVariantPage) {
            // Entry page: restrict disease list to ONLY those that appear as h4 headings
            // in the "Disease & Variants" section of the current UniProt page.
            const scraped = s.scrapedDiseases || [];
            if (scraped.length > 0) {
                // IDs like 'DI-01127' from /diseases/DI-01127 links
                const scrapedIds = new Set(scraped.map(d => d.id).filter(Boolean));
                // Abbreviations like 'DRVT' extracted from labels like "Dravet syndrome (DRVT)"
                const scrapedAbbrs = new Set();
                scraped.forEach(d => {
                    const m = (d.label || '').match(/\(([A-Z][A-Z0-9]+)\)\s*$/);
                    if (m) scrapedAbbrs.add(m[1]);
                });
                // Build label → its own disease ID using the correct 1-to-1 diseasePairs mapping
                // (avoids the bug where all labels on a variant inherit all the variant's IDs)
                const labelToId = new Map();
                s.variants.forEach(v => {
                    (v.diseasePairs || []).forEach(({ id, label }) => {
                        if (id && !labelToId.has(label)) labelToId.set(label, id);
                    });
                });
                const filtered = diseasesToShow.filter(name => {
                    // No 'Unclassified' on entry page — only diseases listed in h4 headings
                    if (name === 'Unclassified') return false;
                    // Exact abbreviation match (e.g. 'DRVT' in scrapedAbbrs)
                    if (scrapedAbbrs.has(name)) return true;
                    // 1-to-1 disease ID match via diseasePairs
                    const id = labelToId.get(name);
                    if (id && scrapedIds.has(id)) return true;
                    return false;
                });
                if (filtered.length > 0) diseasesToShow = filtered;
            }
        }
        if (diseasesToShow.length) {
            dis.classList.remove('ufv-hidden');
            s.activeDiseases = new Set(diseasesToShow);
            diseasesToShow.forEach(name => {
                const meta = ds[name] || { color: '#9e9e9e', count: 0 };
                list.appendChild(makeFilterItem(name, meta.color, meta.count, true, checked => {
                    checked ? s.activeDiseases.add(name) : s.activeDiseases.delete(name);
                    applyMode();
                }));
            });
        } else {
            dis.classList.add('ufv-hidden');
            s.activeDiseases = null;
        }
    }

    function fillFilterList(id, summary, activeSet, onChange) {
        const list = byId(id);
        list.textContent = '';
        Object.entries(summary).forEach(([label, meta]) => {
            activeSet.add(label);
            list.appendChild(makeFilterItem(label, meta.color, meta.count, true, checked => {
                checked ? activeSet.add(label) : activeSet.delete(label);
                onChange();
            }));
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
            hotspots: [['#b9c2cf', 'Not enriched'], ['#ffa726', 'Weak'], ['#e64a19', 'Moderate'], ['#b71c1c', 'Strong']],
            distantContacts: [['#b9c2cf', 'No contact hub'], ['#ab47bc', 'Moderate hub'], ['#6a1b9a', 'Strong hub']],
            alphaMissense: [['#3d85c8', 'Likely benign (<0.34)'], ['#b9c2cf', 'Ambiguous (0.34–0.564)'], ['#e06666', 'Likely pathogenic (0.564–0.78)'], ['#b71c1c', 'Pathogenic (>0.78)']],
            residueBurden: [['#b9c2cf', 'Low burden'], ['#e65100', 'Burden hotspot']],
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

        // In focus/zoom mode always append annotation context on top of coloring legend
        const s = UFVState.state;
        if (s.selectedResidue !== null) {
            if (s.currentMode === 'variant') {
                const summary = DataProcessor.getConsequenceSummary(s.variants);
                Object.entries(summary).forEach(([label, meta]) => {
                    if (!s.activeConsequences || s.activeConsequences.has(label)) {
                        appendLegendItem(meta.color, label);
                    }
                });
            } else if (s.currentMode === 'ptm') {
                Object.entries(s.ptmGroups || {}).forEach(([cat, g]) => {
                    if (g.visible) {
                        const color = DataProcessor.PTM_COLORS[cat] || DataProcessor.PTM_COLORS['default'];
                        appendLegendItem(color, cat);
                    }
                });
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

    function onClick(data, _mode, chain = null) {
        const s = UFVState.state;
        const pos = Number(data.position);
        s.selectedResidue = pos;
        const annotations = buildAnnotationMap();
        s.nearbyResidues = StructureViewer.focusResidue(pos, chain, { annotatedResidues: annotations }) || new Set([pos]);
        const body = byId('ufv-details-body');
        body.textContent = '';
        byId('ufv-details-title').textContent = `Residue ${pos}`;

        // Compact 2-column grid: Position | Nearby
        const topGrid = document.createElement('div');
        topGrid.className = 'ufv-detail-grid';
        const posCell = document.createElement('div');
        posCell.className = 'ufv-detail-cell';
        posCell.innerHTML = `<span class="ufv-detail-lbl">Position</span><span class="ufv-detail-val">${pos}</span>`;
        const nearCell = document.createElement('div');
        nearCell.className = 'ufv-detail-cell';
        nearCell.innerHTML = `<span class="ufv-detail-lbl">Nearby</span><span class="ufv-detail-val ufv-nearby-val">${Array.from(s.nearbyResidues).sort((a, b) => a - b).join(', ')}</span>`;
        topGrid.append(posCell, nearCell);
        body.appendChild(topGrid);

        // PTM annotations
        const ptmsAtPos = s.ptms.filter(p => p.position === pos || p.endPosition === pos);
        ptmsAtPos.forEach(p => body.appendChild(row('PTM', `${p.category}: ${p.description}`, p.color)));

        // Divider between PTM rows and variant rows when both are present
        const variants = s.variants.filter(v => v.position === pos).slice(0, 12);
        if (ptmsAtPos.length > 0 && variants.length > 0) {
            const ptmVarSep = document.createElement('div');
            ptmVarSep.className = 'ufv-variant-divider';
            body.appendChild(ptmVarSep);
        }

        // Variant blocks separated by dividers (AlphaMissense moved to dropdown below)
        variants.forEach((v, i) => {
            if (i > 0) {
                const sep = document.createElement('div');
                sep.className = 'ufv-variant-divider';
                body.appendChild(sep);
            }
            const vblock = document.createElement('div');
            vblock.className = 'ufv-vblock';
            const vtag = document.createElement('div');
            vtag.className = 'ufv-vtag';
            vtag.textContent = `${v.wildType}${v.position}${v.mutant}`;
            vtag.style.color = v.consequenceColor;
            vblock.appendChild(vtag);
            vblock.appendChild(row('ClinVar', v.clinVarSignificance || v.consequence || '—'));
            if (v.clinVarReviewStatus) vblock.appendChild(row('Review', v.clinVarReviewStatus));
            if (v.rsIds?.length) vblock.appendChild(row('dbSNP', v.rsIds.join(', ')));
            body.appendChild(vblock);
        });

        // AlphaMissense full position profile — all 19 possible substitutions from CSV map
        const wt = s.sequence?.[pos - 1] || '';
        const AM_AAS = 'ACDEFGHIKLMNPQRSTVWY';
        const amProfileEntries = [];
        if (s.amMap && wt) {
            for (const mut of AM_AAS) {
                if (mut === wt) continue;
                const key = `${wt}${pos}${mut}`;
                const sc = s.amMap.get(key);
                if (Number.isFinite(sc)) amProfileEntries.push({ mut, score: sc });
            }
        }
        if (amProfileEntries.length > 0) {
            amProfileEntries.sort((a, b) => b.score - a.score);
            const amSection = document.createElement('div');
            amSection.className = 'ufv-am-section';
            const amToggle = document.createElement('button');
            amToggle.className = 'ufv-am-toggle';
            const avgScore = amProfileEntries.reduce((sum, e) => sum + e.score, 0) / amProfileEntries.length;
            const avgColor = avgScore >= 0.564 ? '#ef5350' : avgScore >= 0.34 ? '#ffa726' : '#66bb6a';
            amToggle.innerHTML = `<span class="ufv-am-hdr-left">AlphaMissense<small class="ufv-am-wt">${wt}${pos}</small></span><span class="ufv-am-hdr-right"><span class="ufv-am-avg" style="color:${avgColor}">${avgScore.toFixed(3)}</span><span class="ufv-am-arrow">▾</span></span>`;
            const amBody = document.createElement('div');
            amBody.className = 'ufv-am-body';
            const grid = document.createElement('div');
            grid.className = 'ufv-am-grid';
            amProfileEntries.forEach(({ mut, score }) => {
                const scoreColor = score >= 0.564 ? '#ef5350' : score >= 0.34 ? '#ffa726' : '#66bb6a';
                const cell = document.createElement('div');
                cell.className = 'ufv-am-cell';
                cell.title = score >= 0.564 ? 'Likely pathogenic' : score >= 0.34 ? 'Ambiguous' : 'Likely benign';
                cell.innerHTML = `<span class="ufv-am-cell-mut" style="color:${scoreColor}">${mut}</span><span class="ufv-am-cell-sc" style="color:${scoreColor}">${score.toFixed(3)}</span>`;
                grid.appendChild(cell);
            });
            amBody.appendChild(grid);
            amToggle.addEventListener('click', () => {
                const open = amBody.classList.toggle('show');
                amToggle.querySelector('.ufv-am-arrow').textContent = open ? '▴' : '▾';
            });
            amSection.append(amToggle, amBody);
            body.appendChild(amSection);
        }

        byId('ufv-details').classList.add('show');
        renderLegend(getColorMode());
        renderSequence();
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

    function buildAnnotationMap() {
        const map = new Map();
        activePtms().forEach(p => map.set(p.position, { color: p.color }));
        UFVState.state.variants.forEach(v => map.set(v.position, { color: v.consequenceColor }));
        return map;
    }

    function makeFilterItem(label, color, count, checked, onChange) {
        const el = document.createElement('label');
        el.className = 'ufv-filter-item';
        el.innerHTML = `<input type="checkbox"><span class="ufv-dot"></span><span class="ufv-filter-label"></span><span class="ufv-filter-count"></span>`;
        el.querySelector('input').checked = checked;
        el.querySelector('input').addEventListener('change', e => onChange(e.target.checked));
        el.querySelector('.ufv-dot').style.backgroundColor = color || '#888';
        el.querySelector('.ufv-filter-label').textContent = label;
        el.querySelector('.ufv-filter-count').textContent = count;
        return el;
    }

    function makeExpandableFilter(label, color, items, checked, onChange, itemLabel) {
        const box = document.createElement('div');
        box.className = 'ufv-filter-group';

        // Create children first so the parent onChange closure can reference it
        const children = document.createElement('div');
        children.className = 'ufv-filter-children ufv-collapsed';
        items.slice().sort((a, b) => a.position - b.position).forEach(item => {
            children.appendChild(makeFilterItem(itemLabel(item), color, item.position, item.visible !== false, val => {
                item.visible = val;
                applyPTMMode();
            }));
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

    function varSectionSetAll(section, select) {
        const s = UFVState.state;
        const target = section === 'consequence' ? s.activeConsequences : section === 'provenance' ? s.activeProvenances : s.activeDiseases;
        if (!target) return;
        target.clear();
        const summary = section === 'consequence' ? DataProcessor.getConsequenceSummary(s.variants) : section === 'provenance' ? DataProcessor.getProvenanceSummary(s.variants) : DataProcessor.getDiseaseSummary(s.variants);
        if (select) Object.keys(summary).forEach(k => target.add(k));
        const id = section === 'consequence' ? 'ufv-cons-list' : section === 'provenance' ? 'ufv-prov-list' : 'ufv-dis-list';
        document.querySelectorAll(`#${id} input[type="checkbox"]`).forEach(cb => cb.checked = select);
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
        const text = UFVExport.rewritePdbBeta(StructureViewer.currentPdbText, s.displayedPositions, st, mode, colorContext);
        UFVExport.downloadText(`${s.uniprotId}_${st.id || 'structure'}_${mode}_beta.pdb`, text);
    }

    function exportCsv() {
        const s = UFVState.state;
        if (!s.sequence) return;
        const text = UFVExport.buildResidueMatrix(s.sequence, s.ptms, s.ptmGroups || {}, s.variants, s.amMap, s.analysis, UFVState.selectedStructure());
        UFVExport.downloadText(`${s.uniprotId}_residue_annotations.csv`, text, 'text/csv');
    }

    async function cycleStructure(delta) {
        const s = UFVState.state;
        if (!s.structures.length) return;
        s.selectedStructureIndex = (s.selectedStructureIndex + delta + s.structures.length) % s.structures.length;
        await loadSelectedStructure();
    }

    function cycleTheme() {
        const s = UFVState.state;
        const cycle = ['auto', 'light', 'dark'];
        s.theme = cycle[(cycle.indexOf(s.theme) + 1) % cycle.length];
        modalEl.classList.remove('ufv-light', 'ufv-dark');
        if (s.theme !== 'auto') modalEl.classList.add(`ufv-${s.theme}`);
        const dark = s.theme === 'dark' || (s.theme === 'auto' && window.matchMedia('(prefers-color-scheme: dark)').matches);
        StructureViewer.viewer?.setBackgroundColor(dark ? '#0c111b' : '#f0f2f5');
        StructureViewer.viewer?.render();
    }

    function syncSettingsControls() {
        const set = UFVState.state.settings;
        setColorMode(set.coloringMode);
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
        _gsgt9LockActive = false;
        if (overlayEl) overlayEl.style.display = 'none';
        document.documentElement.style.overflow = '';
        document.body.style.overflow = '';
    }

    // Re-render the 3Dmol canvas whenever the user returns to this browser tab.
    // WebGL contexts can be discarded by the browser while the tab is backgrounded;
    // resize() re-calculates canvas dimensions and render() redraws the scene.
    // Also resets the camera to a full-structure view so the user sees a clean
    // starting position rather than wherever the model was when they switched away.
    document.addEventListener('visibilitychange', () => {
        if (document.hidden || !overlayEl || overlayEl.style.display === 'none') return;
        const v = StructureViewer.viewer;
        if (!v) return;
        try {
            v.resize();
            v.zoomTo({}, 0); // instant reset (0 ms animation) to full-structure view
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

    return { createButton, open, close, prefetchData };
})();
