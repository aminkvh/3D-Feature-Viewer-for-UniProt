/* global DataProcessor */
const UFVApi = (() => {
    'use strict';

    const ALPHAFOLD_API = (id) => `https://alphafold.ebi.ac.uk/api/prediction/${id}`;
    const FEATURES_URL = (id) => `https://www.ebi.ac.uk/proteins/api/features/${id}`;
    const VARIATION_URL = (id) => `https://www.ebi.ac.uk/proteins/api/variation/${id}`;
    const PROTEOMICS_PTM_URL = (id) => `https://www.ebi.ac.uk/proteins/api/proteomics-ptm/${id}`;
    const UNIPROT_URL = (id) => `https://rest.uniprot.org/uniprotkb/${id}.json`;
    const AM_CSV_URL = (id) => `https://alphafold.ebi.ac.uk/files/AF-${id}-F1-aa-substitutions.csv`;
    const PAE_URL = (id, ver) => `https://alphafold.ebi.ac.uk/files/AF-${id}-F1-predicted_aligned_error_v${ver}.json`;
    const PDBe_BEST = (id) => `https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/${id}`;
    const PDBe_ENTRY = (pdb) => `https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/${pdb.toLowerCase()}`;
    const PDBe_MOLECULES = (pdb) => `https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/${pdb.toLowerCase()}`;
    const PDBe_SIFTS = (pdb) => `https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/${pdb.toLowerCase()}`;
    const PDBe_RESIDUES = (pdb) => `https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/${pdb.toLowerCase()}`;
    const PDBe_PDB = (pdb) => `https://www.ebi.ac.uk/pdbe/entry-files/download/pdb${pdb.toLowerCase()}.ent`;
    const PDBe_CIF = (pdb) => `https://www.ebi.ac.uk/pdbe/entry-files/download/${pdb.toLowerCase()}.cif`;
    // 3D-Beacons (PDBe-KB) aggregates computed/predicted models (SWISS-MODEL, ModelArchive,
    // AlphaFold, PED, …) for a UniProt accession alongside experimental structures.
    const BEACONS_SUMMARY = (id) => `https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/${id}.json`;
    // RCSB Chemical Component Dictionary entry for a ligand CCD code (e.g. ABU = GABA).
    const LIGAND_CCD = (ccd) => `https://data.rcsb.org/rest/v1/core/chemcomp/${encodeURIComponent(String(ccd).toUpperCase())}`;
    // PubChem PUG-REST: the published 881-bit CACTVS 2D substructure fingerprint, by InChIKey.
    const PUBCHEM_FP = (ikey) => `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/${encodeURIComponent(ikey)}/property/Fingerprint2D/JSON`;
    // Hosts we will actually fetch a computed model file from (reputable model providers only).
    const BEACON_ALLOWED_HOSTS = new Set([
        'alphafold.ebi.ac.uk', 'www.ebi.ac.uk', 'files.rcsb.org', 'swissmodel.expasy.org',
        'www.modelarchive.org', 'modelarchive.org', 'proteinensemble.org', 'pdb-ihm.org',
        'alphafill.eu',
    ]);
    const RCSB_ENTRY_INSTANCES = (pdbId) =>
        'https://data.rcsb.org/graphql?query=' + encodeURIComponent(
            `{entry(entry_id:"${pdbId.toUpperCase()}"){polymer_entities{polymer_entity_instances{rcsb_polymer_entity_instance_container_identifiers{auth_asym_id}rcsb_polymer_instance_info{modeled_residue_count}}}}}`);

    async function fetchJson(url, options = {}) {
        const res = await fetch(url, { headers: { Accept: 'application/json' }, ...options });
        if (!res.ok) throw new Error(`${url} failed: ${res.status}`);
        return res.json();
    }

    async function fetchOptionalJson(url) {
        try {
            return await fetchJson(url);
        } catch (_) {
            return null;
        }
    }

    async function fetchOptionalText(url) {
        try {
            const res = await fetch(url);
            if (!res.ok) return null;
            return res.text();
        } catch (_) {
            return null;
        }
    }

    function parseAmCsv(text) {
        if (!text) return new Map();
        const map = new Map();
        const lines = text.split('\n');
        for (let i = 1; i < lines.length; i++) {
            const line = lines[i].trim();
            if (!line) continue;
            const comma = line.indexOf(',');
            if (comma < 0) continue;
            const variant = line.slice(0, comma);
            const score = Number(line.slice(comma + 1).split(',')[0]);
            if (variant && Number.isFinite(score)) map.set(variant, score);
        }
        return map;
    }

    async function fetchText(url) {
        const res = await fetch(url);
        if (!res.ok) throw new Error(`${url} failed: ${res.status}`);
        const len = parseInt(res.headers.get('content-length') || '0', 10);
        if (len > 50 * 1024 * 1024) throw new Error('Structure file is larger than 50 MB');
        const text = await res.text();
        if (text.length > 50 * 1024 * 1024) throw new Error('Structure file is larger than 50 MB');
        return text;
    }

    async function getAlphaFoldStructure(id, sequenceLength) {
        let entries = [];
        try {
            const data = await fetchJson(ALPHAFOLD_API(id));
            entries = Array.isArray(data) ? data : [];
        } catch (_) {
            entries = [];
        }

        let entry = entries.find(e => e.pdbUrl) || null;
        if (!entry) {
            for (const ver of [6, 5, 4, 3]) {
                const url = `https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v${ver}.pdb`;
                try {
                    const probe = await fetch(url, { method: 'HEAD' });
                    if (probe.ok) {
                        entry = { pdbUrl: url, modelCreatedDate: '', latestVersion: ver };
                        break;
                    }
                } catch (_) {}
            }
        }
        if (!entry?.pdbUrl) return null;
        return {
            id: `AF-${id}`,
            label: `AlphaFold ${id}`,
            source: 'AlphaFold',
            pdbId: null,
            chainId: null,
            url: entry.pdbUrl,
            method: 'Predicted model',
            resolution: null,
            coverage: 100,
            rangeText: sequenceLength ? `1-${sequenceLength}` : 'full sequence',
            mappedRanges: sequenceLength ? [{ uniprotStart: 1, uniprotEnd: sequenceLength, pdbStart: 1, pdbEnd: sequenceLength, chainId: null }] : [],
            otherChains: false,
            mappingStatus: 'AlphaFold numbering follows UniProt sequence positions.',
            version: entry.latestVersion || '',
        };
    }

    function normalizeBestStructure(item, sequenceLength) {
        const start = parseInt(item.uniprot_start || item.unp_start || item.start || 0, 10);
        const end = parseInt(item.uniprot_end || item.unp_end || item.end || 0, 10);
        const pdbStart = parseInt(item.pdb_start || item.start || start || 0, 10);
        const pdbEnd = parseInt(item.pdb_end || item.end || end || 0, 10);
        const covered = start && end ? Math.max(0, end - start + 1) : 0;
        // Prefer PDBe's pre-computed observed-residue fraction (0–1) over the raw range
        // estimate, which counts all residues in the SIFTS mapping range including
        // unresolved ones and therefore overstates coverage.
        const coverage = item.coverage != null
            ? Math.round(item.coverage * 1000) / 10
            : (sequenceLength && covered ? Math.round((covered / sequenceLength) * 1000) / 10 : 0);
        return {
            id: `${String(item.pdb_id || item.pdbId || '').toUpperCase()}_${item.chain_id || item.chain || 'A'}`,
            label: `${String(item.pdb_id || '').toUpperCase()} chain ${item.chain_id || item.chain || '?'}`,
            source: 'PDB',
            pdbId: String(item.pdb_id || '').toUpperCase(),
            chainId: item.chain_id || item.chain || '',
            url: PDBe_PDB(item.pdb_id || ''),
            cifUrl: PDBe_CIF(item.pdb_id || ''),
            method: item.experimental_method || item.method || 'Experimental',
            resolution: item.resolution || null,
            coverage,
            rangeText: start && end ? `${start}-${end}` : 'mapping unavailable',
            mappedRanges: start && end ? [{ uniprotStart: start, uniprotEnd: end, pdbStart, pdbEnd, chainId: item.chain_id || item.chain || '' }] : [],
            otherChains: false,
            mappingStatus: start && end ? 'Mapped with PDBe/SIFTS best-structure data.' : 'Residue mapping unavailable.',
        };
    }

    // Cache residue_listing responses keyed by pdbId (lower-case) so multiple chains
    // of the same structure share a single fetch.
    const _residueListing = new Map(); // pdbId → Promise<json|null>
    function fetchResidueListing(pdbId) {
        const key = pdbId.toLowerCase();
        if (!_residueListing.has(key)) _residueListing.set(key, fetchOptionalJson(PDBe_RESIDUES(pdbId)));
        return _residueListing.get(key);
    }

    // Return the raw residue array for one chain (or [] if not found).
    async function fetchChainResidues(pdbId, chainId) {
        const listing = await fetchResidueListing(pdbId);
        if (!listing) return [];
        const mols = listing[pdbId.toLowerCase()]?.molecules || [];
        for (const mol of mols)
            for (const chain of (mol.chains || []))
                if (chain.chain_id === chainId) return chain.residues || [];
        return [];
    }

    /**
     * Build a Map from SEQRES residue_number → author_residue_number for one chain.
     * Used to correct SIFTS segments where author_residue_number is null.
     */
    async function buildSeqresToAuthorMap(pdbId, chainId) {
        const residues = await fetchChainResidues(pdbId, chainId);
        if (!residues.length) return null;
        const map = new Map();
        for (const r of residues)
            if (r.author_residue_number != null) map.set(r.residue_number, r.author_residue_number);
        return map.size > 0 ? map : null;
    }

    // Cache entry-level RCSB data (one request per PDB entry covers all chains).
    const _rcsbEntryCache = new Map();
    function fetchRcsbEntry(pdbId) {
        const key = pdbId.toLowerCase();
        if (!_rcsbEntryCache.has(key))
            _rcsbEntryCache.set(key, fetchOptionalJson(RCSB_ENTRY_INSTANCES(pdbId)));
        return _rcsbEntryCache.get(key);
    }

    // Cache the per-entry SIFTS mapping (covers our protein AND all partner proteins).
    const _siftsCache = new Map();
    function fetchSifts(pdbId) {
        const key = pdbId.toLowerCase();
        if (!_siftsCache.has(key)) _siftsCache.set(key, fetchOptionalJson(PDBe_SIFTS(pdbId)));
        return _siftsCache.get(key);
    }

    // How many distinct partner accessions we will fetch variant data for (bounds network use
    // on large hetero-complexes).
    const PARTNER_ACCESSION_CAP = 8;

    /**
     * From the per-entry SIFTS mapping, list the OTHER proteins (different UniProt accession)
     * present in this PDB entry and the chains they occupy.  Used so the 3-D hotspot test can
     * account for disease residues on neighbouring partner proteins — without annotating them.
     * Returns [{ accession, chainId, uniprotStart, uniprotEnd, startAuthor, startSeqres }].
     */
    async function extractPartnerMappings(pdbId, uniprotId, ourChainIds) {
        const sifts = await fetchSifts(pdbId);
        const uni = sifts?.[pdbId.toLowerCase()]?.UniProt;
        if (!uni) return [];
        const ourId = (uniprotId || '').toUpperCase();
        const ourBase = ourId.split('-')[0];
        const ours = new Set(ourChainIds || []);
        const out = [];
        for (const [acc, entry] of Object.entries(uni)) {
            const accU = acc.toUpperCase();
            if (accU === ourId || accU === ourBase) continue; // our protein, skip
            (entry.mappings || []).forEach(m => {
                if (ours.has(m.chain_id)) return; // chain already counted as ours
                out.push({
                    accession: acc,
                    chainId: m.chain_id,
                    uniprotStart: m.unp_start,
                    uniprotEnd: m.unp_end,
                    startAuthor: m.start?.author_residue_number ?? null,
                    startSeqres: m.start?.residue_number ?? null,
                });
            });
        }
        return out;
    }

    // Build a UniProt-position → PDB author-residue resolver for one partner chain mapping.
    async function partnerUniToAuthor(pdbId, m) {
        if (m.startAuthor != null) {
            return pos => m.startAuthor + (pos - m.uniprotStart); // linear author numbering
        }
        // Author numbers absent in SIFTS — fall back to the residue_listing seqres→author map.
        const map = await buildSeqresToAuthorMap(pdbId, m.chainId);
        if (map && m.startSeqres != null) {
            return pos => map.get(m.startSeqres + (pos - m.uniprotStart)) ?? null;
        }
        return () => null;
    }

    // Cache classified partner residues per (pdbId, structure identity) so re-opening a
    // structure doesn't refetch partner variant data.
    const _partnerClassifiedCache = new Map();

    /**
     * Fetch variant data for the partner proteins of a loaded structure and map their
     * pathogenic / benign positions to PDB author residue numbers on the partner chains.
     * Returns [{ chainId, pdbResi, path }] — `path` true = pathogenic, false = benign-only.
     * Used by computeHotspots to fold neighbouring-protein disease residues into the spatial
     * enrichment test.  Annotations for these proteins are NOT displayed.
     */
    async function loadPartnerClassified(structure) {
        const partners = structure?.partners || [];
        if (!partners.length) return [];
        const cacheKey = `${structure.pdbId}|${(structure.chainIds || [structure.chainId]).join(',')}`;
        if (_partnerClassifiedCache.has(cacheKey)) return _partnerClassifiedCache.get(cacheKey);

        const promise = (async () => {
            // Group mappings by accession; cap the number of distinct proteins we fetch.
            const byAcc = new Map();
            partners.forEach(p => { if (!byAcc.has(p.accession)) byAcc.set(p.accession, []); byAcc.get(p.accession).push(p); });
            const accs = [...byAcc.keys()].slice(0, PARTNER_ACCESSION_CAP);
            const points = [];
            await Promise.all(accs.map(async acc => {
                const variation = await fetchOptionalJson(VARIATION_URL(acc));
                const variants = DataProcessor.extractVariants(variation);
                if (!variants.length) return;
                const pathSet = new Set();
                variants.forEach(v => { if (/pathogenic|deleterious/i.test(v.consequence || '')) pathSet.add(v.position); });
                const benignSet = new Set();
                variants.forEach(v => { if (!pathSet.has(v.position) && /benign/i.test(v.consequence || '')) benignSet.add(v.position); });
                if (!pathSet.size && !benignSet.size) return;
                for (const m of byAcc.get(acc)) {
                    const resolve = await partnerUniToAuthor(structure.pdbId, m);
                    const add = (set, path) => set.forEach(pos => {
                        if (pos < m.uniprotStart || pos > m.uniprotEnd) return;
                        const author = resolve(pos);
                        if (author != null) points.push({ chainId: m.chainId, pdbResi: author, path });
                    });
                    add(pathSet, true);
                    add(benignSet, false);
                }
            }));
            return points;
        })();
        _partnerClassifiedCache.set(cacheKey, promise);
        return promise;
    }

    // Fetch modelled residue count matching by author chain ID (auth_asym_id),
    // not label_asym_id, so auxiliary/accessory subunits don't get mixed up.
    async function fetchModelledResidueCount(pdbId, chainId) {
        const data = await fetchRcsbEntry(pdbId);
        const entities = data?.data?.entry?.polymer_entities || [];
        for (const entity of entities)
            for (const inst of (entity.polymer_entity_instances || []))
                if (inst.rcsb_polymer_entity_instance_container_identifiers?.auth_asym_id === chainId)
                    return inst.rcsb_polymer_instance_info?.modeled_residue_count ?? null;
        return null;
    }

    async function decoratePdbStructures(structures, uniprotId, sequenceLength) {
        await Promise.all(structures.map(async s => {
            if (!s.pdbId) return;
            const [summary, molecules, sifts] = await Promise.all([
                fetchOptionalJson(PDBe_ENTRY(s.pdbId)),
                fetchOptionalJson(PDBe_MOLECULES(s.pdbId)),
                uniprotId ? fetchSifts(s.pdbId) : Promise.resolve(null),
            ]);
            const sum = summary?.[s.pdbId.toLowerCase()]?.[0];
            if (sum) {
                s.method = sum.experimental_method?.[0] || s.method;
                s.resolution = s.resolution || sum.resolution;
                s.title = sum.title || '';
            }
            const mols = molecules?.[s.pdbId.toLowerCase()] || [];
            const currentChains = new Set([s.chainId]);
            s.otherChains = mols.some(m => (m.in_chains || []).some(c => !currentChains.has(c)));
            const chainMol = mols.find(m => (m.in_chains || []).includes(s.chainId));
            if (chainMol) s.chainName = chainMol.molecule_name?.[0] || chainMol.molecule_type || '';

            // Override mappedRanges with accurate SIFTS per-segment mappings
            if (uniprotId && sifts) {
                const siftsUni = sifts?.[s.pdbId.toLowerCase()]?.UniProt;
                if (siftsUni) {
                    let siftsEntry = siftsUni[uniprotId] || siftsUni[uniprotId.split('-')[0]];
                    if (!siftsEntry) {
                        const key = Object.keys(siftsUni).find(k => k.toUpperCase() === uniprotId.toUpperCase());
                        if (key) siftsEntry = siftsUni[key];
                    }
                    if (siftsEntry?.mappings) {
                        const rawMappings = siftsEntry.mappings.filter(m => m.chain_id === s.chainId);

                        // Detect segments where SIFTS didn't provide author_residue_number.
                        // In that case we look up the residue_listing API to get the exact
                        // author numbering (the number 3Dmol.js reads from ATOM records).
                        const needsCorrection = rawMappings.some(
                            m => m.start?.author_residue_number == null && m.start?.residue_number != null
                        );
                        let seqresToAuthor = null;
                        if (needsCorrection) {
                            seqresToAuthor = await buildSeqresToAuthorMap(s.pdbId, s.chainId);
                        }

                        /** Resolve the correct PDB author residue number for one endpoint. */
                        function resolveAuthorResi(endpoint, unpFallback) {
                            if (endpoint?.author_residue_number != null) return endpoint.author_residue_number;
                            if (seqresToAuthor && endpoint?.residue_number != null) {
                                const author = seqresToAuthor.get(endpoint.residue_number);
                                if (author != null) return author;
                            }
                            return unpFallback; // last resort: identity with UniProt position
                        }

                        const chainMappings = rawMappings
                            .map(m => ({
                                uniprotStart: m.unp_start,
                                uniprotEnd: m.unp_end,
                                pdbStart: resolveAuthorResi(m.start, m.unp_start),
                                pdbEnd: resolveAuthorResi(m.end, m.unp_end),
                                chainId: m.chain_id,
                            }))
                            .sort((a, b) => a.uniprotStart - b.uniprotStart);

                        if (chainMappings.length > 0) {
                            s.mappedRanges = chainMappings;
                            const uniStart = Math.min(...chainMappings.map(r => r.uniprotStart));
                            const uniEnd = Math.max(...chainMappings.map(r => r.uniprotEnd));
                            s.rangeText = `${uniStart}-${uniEnd}`;
                            s.mappingStatus = `Mapped with SIFTS (${chainMappings.length} segment${chainMappings.length > 1 ? 's' : ''})`;
                        }
                    } else {
                        // SIFTS was fetched but has no entry for this protein (e.g. chimeric
                        // or engineered structures like 6CDU).  best_structures returns
                        // sequential SEQRES positions in start/end, which differ from the
                        // author residue numbers actually stored in the PDB ATOM records.
                        // Fetch the residue_listing to build a SEQRES→author map so that
                        // viewer.js can resolve every UniProt position to the correct resi.
                        const seqresToAuthor = await buildSeqresToAuthorMap(s.pdbId, s.chainId);
                        if (seqresToAuthor) {
                            s.seqresToAuthor = seqresToAuthor;
                            s.mappedRanges = s.mappedRanges.map(r => ({
                                ...r,
                                // Save original SEQRES start so per-residue lookup can work.
                                seqresStart: r.pdbStart,
                                // Update range endpoints to author numbers for display / range checks.
                                pdbStart: seqresToAuthor.get(r.pdbStart) ?? r.pdbStart,
                                pdbEnd:   seqresToAuthor.get(r.pdbEnd)   ?? r.pdbEnd,
                            }));
                            s.mappingStatus = 'Mapped via residue_listing (no SIFTS entry for this protein).';
                        }
                    }
                }
            }
            // Use RCSB GraphQL modeled_residue_count for this chain — the most direct
            // and reliable coverage figure available before the structure is loaded.
            if (sequenceLength > 0) {
                try {
                    const modelled = await fetchModelledResidueCount(s.pdbId, s.chainId);
                    if (modelled != null && modelled > 0)
                        s.coverage = Math.round((modelled / sequenceLength) * 1000) / 10;
                } catch (_) {}
            }
        }));
        return structures;
    }

    async function getExperimentalStructures(id, sequenceLength) {
        const data = await fetchOptionalJson(PDBe_BEST(id));
        const list = data?.[id] || data?.[id.toUpperCase()] || [];
        // Deduplicate same (pdbId, chainId) entries from multiple assemblies first
        const dedupMap = new Map();
        for (const item of list.filter(i => i.pdb_id)) {
            const s = normalizeBestStructure(item, sequenceLength);
            if (!s.pdbId || !s.mappedRanges.length) continue;
            const key = s.id;
            const prev = dedupMap.get(key);
            if (!prev || (s.coverage || 0) > (prev.coverage || 0)) dedupMap.set(key, s);
        }
        const decorated = await decoratePdbStructures([...dedupMap.values()], id, sequenceLength);
        // Merge all chains of the same PDB entry that map to this protein into one entry.
        // This handles homodimers / homo-oligomers: show "7QNE-A,D" instead of separate rows.
        const byPdb = new Map();
        for (const s of decorated) {
            if (!byPdb.has(s.pdbId)) {
                byPdb.set(s.pdbId, {
                    ...s,
                    chainIds: [s.chainId],
                    chainMappings: { [s.chainId]: s.mappedRanges },
                    chainSeqresToAuthor: { [s.chainId]: s.seqresToAuthor || null },
                });
            } else {
                const merged = byPdb.get(s.pdbId);
                if (!merged.chainIds.includes(s.chainId)) {
                    merged.chainIds.push(s.chainId);
                    merged.chainMappings[s.chainId] = s.mappedRanges;
                    merged.chainSeqresToAuthor[s.chainId] = s.seqresToAuthor || null;
                    // Promote higher-coverage chain as primary (drives mappedRanges / analysis)
                    if ((s.coverage || 0) > (merged.coverage || 0)) {
                        merged.chainId    = s.chainId;
                        merged.mappedRanges  = s.mappedRanges;
                        merged.coverage   = s.coverage;
                        merged.rangeText  = s.rangeText;
                        merged.mappingStatus = s.mappingStatus;
                    }
                }
            }
        }
        for (const s of byPdb.values()) {
            s.chainIds.sort();
            if (s.chainIds.length > 1) {
                s.id = `${s.pdbId}_${s.chainIds.join(',')}`;
                s.otherChains = false; // all copies of this protein are represented
            }
        }
        // Attach neighbouring (partner) protein chains from the entry SIFTS so the hotspot
        // analysis can account for their disease residues (e.g. other subunits of a GABA-A
        // receptor).  SIFTS is cached, so this adds no extra network beyond what decorate did.
        await Promise.all([...byPdb.values()].map(async s => {
            s.partners = await extractPartnerMappings(s.pdbId, id, s.chainIds);
        }));
        return [...byPdb.values()];
    }

    async function loadFeatureData(id) {
        const [featuresData, variationData, proteomicsPtmData, uniprotData, amCsvText] = await Promise.all([
            fetchOptionalJson(FEATURES_URL(id)),
            fetchOptionalJson(VARIATION_URL(id)),
            fetchOptionalJson(PROTEOMICS_PTM_URL(id)),
            fetchOptionalJson(UNIPROT_URL(id)),
            fetchOptionalText(AM_CSV_URL(id)),
        ]);
        const amMap = parseAmCsv(amCsvText);
        let ptms = DataProcessor.extractPTMs(featuresData);
        if (proteomicsPtmData) {
            const seen = new Set(ptms.map(p => `${p.position}-${p.category}`));
            DataProcessor.extractProteomicsPTMs(proteomicsPtmData).forEach(p => {
                const key = `${p.position}-${p.category}`;
                if (!seen.has(key)) {
                    seen.add(key);
                    ptms.push(p);
                }
            });
        }
        const variants = DataProcessor.extractVariants(variationData, amMap);
        const sites = DataProcessor.extractSites(featuresData);
        const topology = DataProcessor.extractTopology(featuresData);
        const sequence = uniprotData?.sequence?.value || featuresData?.sequence || variationData?.sequence || '';
        return { featuresData, variationData, proteomicsPtmData, uniprotData, ptms, variants, sites, topology, sequence, amMap };
    }

    /**
     * Parse an AlphaFold predicted-aligned-error (PAE) JSON into a flat Float32Array.
     * Handles both the current AFDB format ({ predicted_aligned_error: [[…]] }) and the
     * legacy triplet format ({ residue1, residue2, distance }).  Returns { n, data } where
     * data[(i-1)*n + (j-1)] is the PAE between 1-based residues i and j, or null.
     * Capped at 2500 residues to bound the NxN memory / downstream O(N²) cost in-browser.
     */
    function parsePae(json) {
        if (!json) return null;
        const obj = Array.isArray(json) ? json[0] : json;
        if (!obj) return null;
        const m = obj.predicted_aligned_error;
        if (Array.isArray(m) && Array.isArray(m[0])) {
            const n = m.length;
            if (n === 0 || n > 2500) return null;
            const data = new Float32Array(n * n);
            for (let i = 0; i < n; i++) {
                const row = m[i];
                for (let j = 0; j < n; j++) data[i * n + j] = +row[j] || 0;
            }
            return { n, data };
        }
        if (Array.isArray(obj.distance) && Array.isArray(obj.residue1) && Array.isArray(obj.residue2)) {
            const r1 = obj.residue1, r2 = obj.residue2, d = obj.distance;
            let n = 0;
            for (let k = 0; k < r1.length; k++) { if (r1[k] > n) n = r1[k]; if (r2[k] > n) n = r2[k]; }
            if (n === 0 || n > 2500) return null;
            const data = new Float32Array(n * n);
            for (let k = 0; k < d.length; k++) data[(r1[k] - 1) * n + (r2[k] - 1)] = +d[k] || 0;
            return { n, data };
        }
        return null;
    }

    const _paeCache = new Map(); // id → Promise<{n, data}|null>
    function getPaeMatrix(id, version) {
        const key = id.toUpperCase();
        if (_paeCache.has(key)) return _paeCache.get(key);
        const p = (async () => {
            const versions = [version, 4, 3].filter((v, i, a) => v && a.indexOf(v) === i);
            for (const ver of versions) {
                const parsed = parsePae(await fetchOptionalJson(PAE_URL(id, ver)));
                if (parsed) return parsed;
            }
            return null;
        })();
        _paeCache.set(key, p);
        return p;
    }

    /**
     * Computed / predicted models from the 3D-Beacons network (SWISS-MODEL, ModelArchive, PED,
     * isoform models, …).  Experimental entries (already covered by PDBe best_structures) and
     * AlphaFold DB (already loaded as our primary model) are skipped, as are formats the viewer
     * can't parse (BCIF) and providers we don't allowlist.  Residue numbering for these models is
     * assumed to follow UniProt positions (true for SWISS-MODEL / AlphaFold-derived models).
     */
    async function get3DBeaconsModels(id, sequenceLength) {
        const data = await fetchOptionalJson(BEACONS_SUMMARY(id));
        const list = data?.structures || [];
        const seen = new Set();
        const out = [];
        for (const entry of list) {
            const sm = entry?.summary || entry;
            if (!sm?.model_url) continue;
            const category = String(sm.model_category || '').toUpperCase();
            const provider = sm.provider || 'Computed model';
            if (category.includes('EXPERIMENTAL')) continue;   // covered by best_structures
            // Include AlphaFold DB models EXCEPT the canonical F1 we already load as the primary
            // model — this brings in non-canonical AlphaFold models (e.g. multi-fragment models
            // for very long proteins) that would otherwise be hidden.
            if (/alphafold/i.test(provider) && String(sm.model_identifier || '').toUpperCase() === `AF-${String(id).toUpperCase()}-F1`) continue;
            const fmt = String(sm.model_format || '').toUpperCase();
            if (fmt !== 'PDB' && fmt !== 'MMCIF') continue;     // skip BCIF / unknown
            let host;
            try { host = new URL(sm.model_url).hostname; } catch (_) { continue; }
            if (!BEACON_ALLOWED_HOSTS.has(host)) continue;
            const mid = sm.model_identifier || sm.model_url;
            if (seen.has(mid)) continue;
            seen.add(mid);
            const uStart = parseInt(sm.uniprot_start || 1, 10);
            const uEnd = parseInt(sm.uniprot_end || sequenceLength || 0, 10);
            const cov = sm.coverage != null
                ? Math.round(sm.coverage * 1000) / 10
                : (sequenceLength && uEnd ? Math.round(((uEnd - uStart + 1) / sequenceLength) * 1000) / 10 : null);
            out.push({
                id: `BCN-${mid}`,
                label: provider,
                source: 'Computed',
                provider,
                modelCategory: sm.model_category || '',
                pdbId: null,
                chainId: null,
                url: sm.model_url,
                cifUrl: fmt === 'MMCIF' ? sm.model_url : null,
                format: fmt === 'MMCIF' ? 'mmcif' : 'pdb',
                method: `${provider} model`,
                resolution: sm.resolution || null,
                coverage: cov,
                rangeText: (uStart && uEnd) ? `${uStart}-${uEnd}` : (sequenceLength ? `1-${sequenceLength}` : ''),
                mappedRanges: (uStart && uEnd) ? [{ uniprotStart: uStart, uniprotEnd: uEnd, pdbStart: uStart, pdbEnd: uEnd, chainId: null }] : [],
                otherChains: false,
                mappingStatus: `${provider} computed model — residue numbering assumed to follow UniProt positions.`,
                modelPageUrl: sm.model_page_url || '',
                confidence: sm.confidence_avg_local_score ?? null,
            });
        }
        return out.slice(0, 25); // cap so the structure selector isn't flooded
    }

    /**
     * AlphaFold models for the protein's NON-canonical isoforms (e.g. AF-P35498-2-F1). These are
     * not listed in the 3D-Beacons summary of the canonical accession, so we read the isoform IDs
     * from the UniProt entry and probe AlphaFold DB for each. Residue numbering follows the
     * isoform sequence, so canonical PTM/variant positions may be offset on these models.
     */
    async function getIsoformAlphaFoldStructures(id, sequenceLength) {
        const base = String(id).split('-')[0];
        const u = await fetchOptionalJson(UNIPROT_URL(base));
        const isoIds = new Set();
        (u?.comments || []).forEach(c => {
            if (c.commentType === 'ALTERNATIVE PRODUCTS') (c.isoforms || []).forEach(iso => (iso.isoformIds || []).forEach(iid => isoIds.add(iid)));
        });
        const targets = [...isoIds].filter(iid => iid && iid !== base && iid !== `${base}-1`);
        const out = [];
        // Resolve each isoform's AlphaFold model via its 3D-Beacons summary. We can't HEAD-probe
        // the AlphaFold file directly: alphafold.ebi.ac.uk returns CORS headers on GET but NOT on
        // HEAD, so an in-browser HEAD is blocked. The Beacons summary (www.ebi.ac.uk) is CORS-OK
        // and gives the GET-able model URL.
        await Promise.all(targets.map(async iid => {
            const d = await fetchOptionalJson(BEACONS_SUMMARY(iid));
            const af = (d?.structures || []).map(e => e?.summary || e).find(sm => sm?.model_url && /alphafold/i.test(sm.provider || ''));
            if (!af) return;
            const fmt = String(af.model_format || '').toUpperCase();
            out.push({
                id: `AF-${iid}`, label: `AlphaFold ${iid}`, source: 'AlphaFold', isoform: iid,
                pdbId: null, chainId: null,
                url: af.model_url,
                cifUrl: fmt === 'MMCIF' ? af.model_url : null,
                format: fmt === 'PDB' ? 'pdb' : 'mmcif',
                method: 'Predicted model (isoform)', resolution: null, coverage: null,
                rangeText: `isoform ${iid}`,
                mappedRanges: sequenceLength ? [{ uniprotStart: 1, uniprotEnd: sequenceLength, pdbStart: 1, pdbEnd: sequenceLength, chainId: null }] : [],
                otherChains: false,
                mappingStatus: 'AlphaFold isoform model — residue numbering follows the isoform sequence; canonical annotations may be offset.',
                version: '',
            });
        }));
        return out;
    }

    async function getStructures(id, sequenceLength) {
        const [alphaFold, isoforms, pdbs, beacons] = await Promise.all([
            getAlphaFoldStructure(id, sequenceLength),
            getIsoformAlphaFoldStructures(id, sequenceLength),
            getExperimentalStructures(id, sequenceLength),
            get3DBeaconsModels(id, sequenceLength),
        ]);
        // Order: canonical AlphaFold, then isoform AlphaFold, then experimental, then computed.
        const rank = s => s.source === 'AlphaFold' ? (s.isoform ? 0.5 : 0) : s.source === 'Computed' ? 2 : 1;
        return [alphaFold, ...isoforms, ...pdbs, ...beacons].filter(Boolean).sort((a, b) => {
            if (rank(a) !== rank(b)) return rank(a) - rank(b);
            return (b.coverage || 0) - (a.coverage || 0);
        });
    }

    // Ligand chemistry from the RCSB Chemical Component Dictionary, keyed by 3-letter CCD code.
    // Returns { id, name, formula, smiles, inchikey, drugbank }. Cached (promise) per code.
    const _ligandCache = new Map();
    function getLigandInfo(ccd) {
        if (!ccd) return Promise.resolve(null);
        const key = String(ccd).toUpperCase();
        if (_ligandCache.has(key)) return _ligandCache.get(key);
        const p = (async () => {
            const j = await fetchOptionalJson(LIGAND_CCD(key));
            const d = j?.rcsb_chem_comp_descriptor || {};
            const db = (j?.rcsb_chem_comp_related || []).find(x => x.resource_name === 'DrugBank');
            return {
                id: key,
                name: j?.chem_comp?.name || null,
                formula: j?.chem_comp?.formula || null,
                smiles: d.SMILES_stereo || d.SMILES || null,
                inchikey: d.InChIKey || null,
                drugbank: db?.resource_accession_code || null,
            };
        })();
        _ligandCache.set(key, p);
        return p;
    }

    // PubChem 2D substructure fingerprint (881-bit CACTVS keys) as a Uint8Array(881), by InChIKey.
    // Published/standard fingerprint, decoded client-side; used for rigorous Tanimoto similarity.
    const _fpCache = new Map();
    function b64ToBytes(b64) {
        const s = atob(b64);
        const a = new Uint8Array(s.length);
        for (let i = 0; i < s.length; i++) a[i] = s.charCodeAt(i);
        return a;
    }
    function getLigandFingerprint(inchikey) {
        if (!inchikey) return Promise.resolve(null);
        if (_fpCache.has(inchikey)) return _fpCache.get(inchikey);
        const p = (async () => {
            const j = await fetchOptionalJson(PUBCHEM_FP(inchikey));
            const b64 = j?.PropertyTable?.Properties?.[0]?.Fingerprint2D;
            if (!b64) return null;
            const bin = b64ToBytes(b64); // first 4 bytes are a length prefix, then the 881 bits
            const bits = new Uint8Array(881);
            for (let i = 0; i < 881; i++) { const byte = bin[4 + (i >> 3)] || 0; bits[i] = (byte >> (7 - (i & 7))) & 1; }
            return bits;
        })();
        _fpCache.set(inchikey, p);
        return p;
    }

    return { loadFeatureData, getStructures, fetchText, loadPartnerClassified, getPaeMatrix, getLigandInfo, getLigandFingerprint };
})();
