/* UniProt 3D Feature Viewer – Options page */
'use strict';

const DEFAULTS = {
    defaultStructure: 'alphafold',
    coloringMode: 'default',
    copyFormat: 'pymol',
    showOptionalTracks: false,
};

function byId(id) { return document.getElementById(id); }

function loadOptions() {
    const doLoad = (stored) => {
        const s = { ...DEFAULTS, ...stored };
        byId('opt-structure').value = s.defaultStructure;
        byId('opt-color').value = s.coloringMode;
        byId('opt-copy').value = s.copyFormat;
        byId('opt-optional').checked = s.showOptionalTracks;
    };
    if (typeof chrome !== 'undefined' && chrome.storage?.local) {
        chrome.storage.local.get(DEFAULTS, doLoad);
    } else {
        try {
            doLoad(JSON.parse(localStorage.getItem('ufv-settings') || '{}'));
        } catch (_) {
            doLoad({});
        }
    }
}

function saveOptions() {
    const patch = {
        defaultStructure: byId('opt-structure').value,
        coloringMode: byId('opt-color').value,
        copyFormat: byId('opt-copy').value,
        showOptionalTracks: byId('opt-optional').checked,
    };
    const onSaved = () => {
        const el = byId('save-status');
        el.textContent = 'Saved.';
        setTimeout(() => { el.textContent = ''; }, 1800);
    };
    if (typeof chrome !== 'undefined' && chrome.storage?.local) {
        chrome.storage.local.set(patch, onSaved);
    } else {
        try {
            const existing = JSON.parse(localStorage.getItem('ufv-settings') || '{}');
            localStorage.setItem('ufv-settings', JSON.stringify({ ...existing, ...patch }));
        } catch (_) {}
        onSaved();
    }
}

document.addEventListener('DOMContentLoaded', () => {
    loadOptions();
    byId('save-btn').addEventListener('click', saveOptions);
});
