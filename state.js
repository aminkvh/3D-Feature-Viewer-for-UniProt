const UFVState = (() => {
    'use strict';

    const defaults = {
        copyFormat: 'pymol',
        defaultStructure: 'alphafold',
        coloringMode: 'default',
        coVisualize: true,
        showClinVar: true,
        showAlphaMissense: true,
        showOptionalTracks: false,
    };

    const state = {
        uniprotId: null,
        pageContext: null,
        currentMode: 'ptm',
        theme: 'auto',
        loaded: false,
        annotationsLoaded: false,
        loadingPromise: null,
        featuresData: null,
        variationData: null,
        ptms: [],
        ptmGroups: {},
        variants: [],
        sites: [],
        sitesVisible: false, // "Site" annotation spheres — off by default (opt-in overlay)
        sequence: '',
        structures: [],
        selectedStructureIndex: 0,
        activeConsequences: new Set(),
        activeProvenances: new Set(),
        activeDiseases: null,
        variantPtmCats: new Set(), // PTM categories co-displayed in the Disease & Variants view
        scrapedDiseases: [],
        displayedPositions: [],
        selectedResidue: null,
        selectedChain: null,
        nearbyResidues: new Set(),
        amMap: null,
        settings: { ...defaults },
        analysis: { hotspots: null, hotspotsByChain: null, hotspotMethod: null, distantContacts: null, distantContactsByChain: null, alphaMissense: new Map(), residueBurden: new Set(), prism: null, ptmVariantProximity: null },
    };

    function resetForProtein(id) {
        state.uniprotId = id;
        state.pageContext = null;
        state.loaded = false;
        state.annotationsLoaded = false;
        state.loadingPromise = null;
        state.featuresData = null;
        state.variationData = null;
        state.ptms = [];
        state.ptmGroups = {};
        state.variants = [];
        state.sites = [];
        state.sitesVisible = false;
        state.sequence = '';
        state.structures = [];
        state.selectedStructureIndex = 0;
        state.activeConsequences = new Set();
        state.activeProvenances = new Set();
        state.activeDiseases = null;
        state.variantPtmCats = new Set();
        state.scrapedDiseases = [];
        state.displayedPositions = [];
        state.selectedResidue = null;
        state.selectedChain = null;
        state.nearbyResidues = new Set();
        state.amMap = null;
        state.analysis = { hotspots: null, hotspotsByChain: null, hotspotMethod: null, distantContacts: null, distantContactsByChain: null, alphaMissense: new Map(), residueBurden: new Set(), prism: null, ptmVariantProximity: null };
    }

    async function loadSettings() {
        const loaded = await new Promise(resolve => {
            try {
                if (typeof chrome !== 'undefined' && chrome.storage?.local) {
                    chrome.storage.local.get(defaults, resolve);
                    return;
                }
            } catch (_) {}
            try {
                resolve({ ...defaults, ...JSON.parse(localStorage.getItem('ufv-settings') || '{}') });
            } catch (_) {
                resolve({ ...defaults });
            }
        });
        state.settings = { ...defaults, ...loaded };
        return state.settings;
    }

    async function saveSettings(patch) {
        state.settings = { ...state.settings, ...patch };
        try {
            if (typeof chrome !== 'undefined' && chrome.storage?.local) {
                chrome.storage.local.set(state.settings);
                return state.settings;
            }
        } catch (_) {}
        try {
            localStorage.setItem('ufv-settings', JSON.stringify(state.settings));
        } catch (_) {}
        return state.settings;
    }

    function selectedStructure() {
        return state.structures[state.selectedStructureIndex] || null;
    }

    return { state, resetForProtein, loadSettings, saveSettings, selectedStructure };
})();
