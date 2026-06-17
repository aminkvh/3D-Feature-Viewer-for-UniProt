/* ============================================
   Data Processing — Chrome Extension
   ============================================ */

// 50-color disease palette — evenly distributed across the visible spectrum, tested for
// deuteranopia/protanopia safety. Expanded well beyond the old 14-color cycling set.
const DISEASE_PALETTE = [
    '#1565c0', '#0288d1', '#0097a7', '#00838f', '#006064',
    '#4527a0', '#7b1fa2', '#9c27b0', '#ba68c8', '#7c4dff',
    '#880e4f', '#c2185b', '#e91e63', '#ec407a', '#f06292',
    '#004d40', '#00695c', '#00897b', '#4db6ac', '#26a69a',
    '#bf360c', '#e65100', '#ef6c00', '#ff8f00', '#ffa000',
    '#311b92', '#3949ab', '#5c6bc0', '#7986cb', '#9fa8da',
    '#01579b', '#0277bd', '#039be5', '#29b6f6', '#80deea',
    '#4e342e', '#6d4c41', '#795548', '#a1887f', '#bcaaa4',
    '#33691e', '#558b2f', '#689f38', '#7cb342', '#9ccc65',
    '#f9a825', '#fbc02d', '#ff6f00', '#ff8f00', '#ffe082',
];

const DataProcessor = {

    // ---- PTM Colors — colorblind-safe palette; no hex shared with CONSEQUENCE_CATEGORIES ----
    PTM_COLORS: {
        'Phosphorylation':              '#8b5cf6',  // violet
        'N6-acetyllysine':              '#06b6d4',  // cyan
        'N6-succinyllysine':            '#0e7490',  // dark cyan
        'N6-methyllysine':              '#16a34a',  // green
        'N6,N6-dimethyllysine':         '#15803d',  // dark green
        'Omega-N-methylarginine':       '#65a30d',  // lime
        'Symmetric dimethylarginine':   '#4ade80',  // light green
        'Asymmetric dimethylarginine':  '#86efac',  // pale green
        'Glycyl lysine isopeptide':     '#c026d3',  // fuchsia
        'Ubiquitination':               '#e879f9',  // magenta
        'Ubiquitinated lysine':         '#e879f9',  // magenta
        'SUMOylation':                  '#9333ea',  // purple
        'ADP-ribosylation':             '#6366f1',  // indigo
        'Glycosylation':                '#be185d',  // rose
        'Disulfide bond':               '#0d9488',  // teal
        'Lipidation':                   '#f472b6',  // pink
        'Cross-link':                   '#a78bfa',  // lavender
        'Acetylation':                  '#22d3ee',  // sky cyan
        'Methylation':                  '#22c55e',  // medium green
        'Citrulline':                   '#2dd4bf',  // turquoise
        'Hydroxyproline':               '#38bdf8',  // light blue
        'Hydroxylation':                '#a3e635',  // yellow-green
        'Nitration':                    '#facc15',  // yellow
        'Pyroglutamic acid':            '#fbbf24',  // amber-yellow
        'S-nitrosocysteine':            '#4f46e5',  // deep indigo
        'Deamidation':                  '#818cf8',  // periwinkle
        'Modified residue':             '#546e7a',  // slate (distinct from consequence grey #9e9e9e)
        'default':                      '#ff8fab',  // light rose
    },

    // ---- Consequence categories — colorblind-safe (vermillion / amber / blue / grey) ----
    // Replaces the old red/green pair that deuteranopes cannot distinguish.
    CONSEQUENCE_CATEGORIES: {
        'Likely pathogenic or pathogenic': { color: '#d55e00' },  // vermillion
        'Predicted deleterious':           { color: '#e6a817' },  // amber
        'Likely benign or benign':         { color: '#0072b2' },  // blue
        'Uncertain significance':          { color: '#9e9e9e' },  // grey
    },

    // ---- Provenance categories ----
    PROVENANCE_CATEGORIES: {
        'UniProt reviewed':    { color: '#42a5f5' },
        'ClinVar':             { color: '#ab47bc' },
        'Large scale studies': { color: '#26c6da' },
    },

    /**
     * Extract PTMs from EBI features API response.
     * Each PTM has a .position (residue number) we render as a Cα sphere.
     */
    extractPTMs(featuresData) {
        if (!featuresData || !featuresData.features) return [];

        const ptmTypes = ['MOD_RES', 'CROSSLNK', 'LIPID', 'CARBOHYD', 'DISULFID'];
        const ptms = [];

        featuresData.features.forEach(f => {
            if (!ptmTypes.includes(f.type)) return;
            const pos = parseInt(f.begin);
            if (isNaN(pos)) return;

            const category = this._categorizePTM(f.description || '', f.type);

            ptms.push({
                position: pos,
                endPosition: parseInt(f.end) || pos,
                type: f.type,
                description: f.description || 'Unknown modification',
                category: category,
                color: this.PTM_COLORS[category] || this.PTM_COLORS['default'],
                evidences: f.evidences || [],
            });
        });

        return ptms;
    },

    /**
     * Extract large-scale PTMs from EBI proteomics-ptm API.
     * These come from PRIDE/PTMeXchange and include phosphorylation etc.
     * The API returns peptide-level data with relative PTM positions.
     */
    extractProteomicsPTMs(proteomicsPtmData) {
        if (!proteomicsPtmData || !proteomicsPtmData.features) return [];
        const ptms = [];
        const seen = new Set(); // deduplicate by position+name

        proteomicsPtmData.features.forEach(f => {
            if (f.type !== 'PROTEOMICS_PTM') return;
            if (!f.ptms || !Array.isArray(f.ptms)) return;

            const peptideStart = parseInt(f.begin);
            if (isNaN(peptideStart)) return;

            f.ptms.forEach(ptm => {
                // ptm.position is relative to the peptide (1-based)
                const absPos = peptideStart + (ptm.position - 1);
                const key = `${absPos}-${ptm.name}`;
                if (seen.has(key)) return;
                seen.add(key);

                const desc = ptm.name || 'Unknown modification';
                const category = this._categorizePTM(desc.toLowerCase(), 'MOD_RES');
                const sources = (ptm.sources || []).join(', ');

                ptms.push({
                    position: absPos,
                    endPosition: absPos,
                    type: 'MOD_RES',
                    description: `${desc} (large scale data, ${sources})`,
                    category: category,
                    color: this.PTM_COLORS[category] || this.PTM_COLORS['default'],
                    evidences: f.evidences || [],
                    isLargeScale: true,
                });
            });
        });

        return ptms;
    },

    /**
     * Extract UniProt "Site" annotations (feature type SITE) — interesting single residues or
     * peptide bonds not covered by other subsections: cleavage sites, protease-inhibitory sites,
     * fusion-protein breakpoints, etc. A site may be one residue (begin === end) or a peptide bond
     * represented by its two flanking residues (begin, end).
     */
    extractSites(featuresData) {
        if (!featuresData || !featuresData.features) return [];
        // Point-like "interesting site" feature types. DNA_BIND is deliberately excluded — it is a
        // region/domain, not a single site, and would render misleading endpoint spheres.
        const SITE_LABELS = { SITE: 'Site', ACT_SITE: 'Active site', BINDING: 'Binding site', METAL: 'Metal binding' };
        const out = [];
        featuresData.features.forEach(f => {
            const typeLabel = SITE_LABELS[f.type];
            if (!typeLabel) return;
            const pos = parseInt(f.begin);
            if (isNaN(pos)) return;
            // SITE descriptions already name the site type (e.g. "Cleavage; by thrombin"). For
            // BINDING/METAL the bound molecule lives in f.ligand (description is often empty), so
            // fold the ligand name in (e.g. "Binding site: ATP" / "Binding site: 4-aminobutanoate").
            const ligandName = f.ligand?.name || '';
            let description;
            if (f.type === 'SITE') {
                description = f.description || 'Site';
            } else {
                const parts = [ligandName, f.description].filter(Boolean);
                description = parts.length ? `${typeLabel}: ${parts.join(' — ')}` : typeLabel;
            }
            out.push({
                position: pos,
                endPosition: parseInt(f.end) || pos,
                description,
                ligand: ligandName || null,
                ligandRef: f.ligand?.dbReference || null,
                category: 'Site',
                color: this.SITE_COLOR,
                evidences: f.evidences || [],
            });
        });
        return out;
    },

    /**
     * Extract "interesting site" features from the UniProt entry JSON (rest.uniprot.org). Richer
     * than the proteins-API features endpoint: binding sites carry the ligand name AND a note
     * (e.g. "ligand shared with the neighboring beta subunit", "agonist"), which the proteins API
     * omits. Used as the preferred site source, merged with the proteins-API sites.
     */
    extractSitesUniProt(uniprotData) {
        if (!uniprotData || !Array.isArray(uniprotData.features)) return [];
        const LABELS = { 'Site': 'Site', 'Active site': 'Active site', 'Binding site': 'Binding site', 'Metal binding': 'Metal binding' };
        const out = [];
        uniprotData.features.forEach(f => {
            const typeLabel = LABELS[f.type];
            if (!typeLabel) return;
            const pos = f.location?.start?.value;
            if (pos == null) return;
            const end = f.location?.end?.value || pos;
            const ligandName = f.ligand?.name || '';
            const ligandNote = f.ligand?.note || '';
            let description;
            if (f.type === 'Site') {
                description = f.description || 'Site';
            } else {
                const parts = [ligandName, ligandNote, f.description].filter(Boolean);
                description = parts.length ? `${typeLabel}: ${parts.join('; ')}` : typeLabel;
            }
            out.push({
                position: pos,
                endPosition: end,
                description,
                ligand: ligandName || null,
                ligandRef: f.ligand?.id || f.featureCrossReferences?.[0]?.id || null,
                category: 'Site',
                color: this.SITE_COLOR,
                evidences: f.evidences || [],
            });
        });
        return out;
    },

    SITE_COLOR: '#fbc02d', // amber — distinct from PTM and variant palettes

    TOPOLOGY_COLORS: {
        transmembrane: '#f9a825', intramembrane: '#ef6c00',
        cytoplasmic: '#1e88e5', extracellular: '#e53935', lumenal: '#43a047',
        periplasmic: '#26a69a', nuclear: '#8e24aa', mitochondrial: '#00897b', other: '#7e57c2',
    },

    /**
     * Extract membrane topology features (TOPO_DOM, TRANSMEM, INTRAMEM) as coloured sequence
     * segments. Returns [{ start, end, type, label, color }] (empty when the entry has none).
     */
    extractTopology(featuresData) {
        if (!featuresData || !featuresData.features) return [];
        const C = this.TOPOLOGY_COLORS;
        const out = [];
        featuresData.features.forEach(f => {
            if (!['TOPO_DOM', 'TRANSMEM', 'INTRAMEM'].includes(f.type)) return;
            const start = parseInt(f.begin), end = parseInt(f.end) || start;
            if (isNaN(start)) return;
            let label, color;
            if (f.type === 'TRANSMEM') { label = 'Transmembrane' + (f.description ? ` (${f.description})` : ''); color = C.transmembrane; }
            else if (f.type === 'INTRAMEM') { label = 'Intramembrane' + (f.description ? ` (${f.description})` : ''); color = C.intramembrane; }
            else {
                const d = (f.description || '').toLowerCase();
                label = f.description || 'Topological domain';
                color = d.includes('cytoplasm') ? C.cytoplasmic
                    : d.includes('extracellular') ? C.extracellular
                    : d.includes('lumen') ? C.lumenal
                    : d.includes('periplasm') ? C.periplasmic
                    : d.includes('nuclear') ? C.nuclear
                    : d.includes('mitochond') ? C.mitochondrial
                    : C.other;
            }
            out.push({ start, end, type: f.type, label, color });
        });
        return out;
    },

    // Feature types that make up the UniProt "Family & Domains" section.
    DOMAIN_TYPES: ['Domain', 'Region', 'Repeat', 'Compositional bias', 'Zinc finger', 'Coiled coil', 'Motif', 'DNA binding'],
    // Qualitatively-distinct palette cycled so every domain feature gets its own colour.
    DOMAIN_PALETTE: ['#1e88e5', '#43a047', '#fb8c00', '#8e24aa', '#00897b', '#e53935', '#3949ab',
        '#c0ca33', '#6d4c41', '#00acc1', '#d81b60', '#5e35b1', '#7cb342', '#f4511e', '#039be5', '#fdd835'],

    /**
     * Extract Family & Domains features (domain / region / repeat / compositional bias / zinc
     * finger / coiled coil / motif / DNA binding) from the UniProt entry JSON. Each gets its own
     * colour. `isRange` is true for multi-residue features (rendered as cartoon colouring) and
     * false for single-residue ones (rendered as a sphere). Returns [] when the entry has none.
     */
    extractDomains(uniprotData) {
        const feats = (uniprotData?.features || []).filter(f => this.DOMAIN_TYPES.includes(f.type));
        const out = [];
        feats.forEach((f, i) => {
            const start = f.location?.start?.value;
            if (start == null) return;
            const end = f.location?.end?.value ?? start;
            const desc = f.description || f.type;
            out.push({
                position: start,
                endPosition: end,
                type: f.type,
                description: f.type === 'Compositional bias' ? `${desc} (compositional bias)`
                    : f.type === desc ? desc : `${desc} (${f.type.toLowerCase()})`,
                color: this.DOMAIN_PALETTE[i % this.DOMAIN_PALETTE.length],
                isRange: end > start,
                visible: true,
            });
        });
        return out.sort((a, b) => a.position - b.position);
    },

    /**
     * UniProt's ProtNLM (and ProtNLM2) AI model names proteins that would otherwise be
     * 'uncharacterized', predicted from sequence alone. Such names appear in the entry's
     * proteinDescription attributed to source 'Google' / id 'ProtNLM' (with a CAUTION note).
     * Returns { name, isAI, source, caution, reviewed } or null. Rule-based automatic names
     * (ARBA/RuleBase/HAMAP) are deliberately NOT flagged as AI.
     */
    extractProtNLM(uniprotData) {
        if (!uniprotData) return null;
        const desc = uniprotData.proteinDescription || {};
        const rec = desc.recommendedName?.fullName;
        const sub = (desc.submissionNames || []).map(s => s.fullName).find(Boolean);
        const fn = (rec && rec.value) ? rec : (sub && sub.value ? sub : null);
        if (!fn || !fn.value) return null;
        let isAI = false, source = '';
        for (const e of (fn.evidences || [])) {
            const tag = `${e.source || ''} ${e.id || ''}`.trim();
            if (tag.toLowerCase().includes('protnlm') || (e.source || '').toLowerCase() === 'google') {
                isAI = true; source = tag; break;
            }
        }
        let caution = '';
        for (const c of (uniprotData.comments || [])) {
            if (c.commentType !== 'CAUTION') continue;
            for (const t of (c.texts || [])) {
                if ((t.value || '').toLowerCase().includes('protnlm')) { caution = t.value; isAI = true; }
            }
        }
        const et = (uniprotData.entryType || '').toLowerCase();
        const reviewed = et.includes('reviewed') && !et.includes('unreviewed');
        return { name: fn.value, isAI, source: source || (isAI ? 'ProtNLM' : ''), caution, reviewed };
    },

    /** UniProt 'Mutagenesis' features — residues experimentally mutated, with the observed effect. */
    extractMutagenesis(uniprotData) {
        const out = [];
        ((uniprotData?.features) || []).filter(f => f.type === 'Mutagenesis').forEach(f => {
            const start = f.location?.start?.value;
            if (start == null) return;
            const end = f.location?.end?.value ?? start;
            const alt = f.alternativeSequence || {};
            const pubmed = (f.evidences || []).filter(e => e.source === 'PubMed' && e.id).map(e => e.id);
            out.push({
                position: start, endPosition: end,
                wildType: alt.originalSequence || '',
                mutants: alt.alternativeSequences || [],
                effect: f.description || '', pubmed, color: '#6d4c41',
            });
        });
        return out.sort((a, b) => a.position - b.position);
    },

    /** Protein-level function context: short summary, subcellular locations, catalytic activity. */
    extractFunction(uniprotData) {
        if (!uniprotData) return null;
        const strip = (t) => (t || '').replace(/\s*\(PubMed:[^)]*\)/g, '').replace(/\s*\(By similarity\)/g, '')
            .replace(/\s*\(Ref\.[^)]*\)/g, '').replace(/\s+/g, ' ').trim();
        let summary = ''; const locations = [], catalytic = [];
        (uniprotData.comments || []).forEach(c => {
            if (c.commentType === 'FUNCTION' && !summary) summary = strip(c.texts?.[0]?.value);
            else if (c.commentType === 'SUBCELLULAR LOCATION') (c.subcellularLocations || []).forEach(sl => {
                const v = sl.location?.value; if (v && !locations.includes(v)) locations.push(v);
            });
            else if (c.commentType === 'CATALYTIC ACTIVITY') { const r = c.reaction?.name; if (r && !catalytic.includes(r)) catalytic.push(r); }
        });
        return (summary || locations.length || catalytic.length) ? { summary, locations, catalytic } : null;
    },

    _categorizePTM(desc, type) {
        const d = desc.toLowerCase();
        if (type === 'DISULFID') return 'Disulfide bond';
        if (type === 'LIPID')    return 'Lipidation';
        if (type === 'CARBOHYD') return 'Glycosylation';
        if (type === 'CROSSLNK') {
            if (d.includes('ubiquitin')) return 'Ubiquitination';
            if (d.includes('sumo'))      return 'SUMOylation';
            return 'Cross-link';
        }
        // Phosphorylation — all variants consolidated into one category
        if (d.includes('phospho'))         return 'Phosphorylation';
        // Acetylation / methylation
        if (d.includes('n6-acetyl'))       return 'N6-acetyllysine';
        if (d.includes('n6-succinyl'))     return 'N6-succinyllysine';
        if (d.includes('n6,n6-dimethyl'))  return 'N6,N6-dimethyllysine';
        if (d.includes('n6-methyl'))       return 'N6-methyllysine';
        if (d.includes('omega-n-methyl'))  return 'Omega-N-methylarginine';
        if (d.includes('symmetric dimethyl')) return 'Symmetric dimethylarginine';
        if (d.includes('asymmetric dimethyl')) return 'Asymmetric dimethylarginine';
        // ADP-ribosylation
        if (d.includes('adp-ribos'))       return 'ADP-ribosylation';
        // Ubiquitinated lysine (MOD_RES type)
        if (d.includes('ubiquitin'))       return 'Ubiquitinated lysine';
        if (d.includes('glycyl lysine isopeptide')) return 'Glycyl lysine isopeptide';
        // Other modified residues
        if (d.includes('citrulline'))      return 'Citrulline';
        if (d.includes('hydroxyproline'))  return 'Hydroxyproline';
        if (d.includes('hydroxylysine'))   return 'Hydroxylation';
        if (d.includes('hydroxylation'))   return 'Hydroxylation';
        if (d.includes('nitrated'))        return 'Nitration';
        if (d.includes('pyroglutam'))      return 'Pyroglutamic acid';
        if (d.includes('s-nitrosocysteine')) return 'S-nitrosocysteine';
        if (d.includes('deamidated'))      return 'Deamidation';
        // Fallback: use the description text before semicolon as category name
        const label = desc.split(';')[0].trim();
        return label || 'Modified residue';
    },

    /**
     * Group PTMs by category, returning { category: { color, items[], visible } }
     */
    groupPTMsByCategory(ptms) {
        const groups = {};
        ptms.forEach(p => {
            if (!groups[p.category]) {
                groups[p.category] = { category: p.category, color: p.color, items: [], visible: true };
            }
            groups[p.category].items.push(p);
        });
        return groups;
    },

    /**
     * Extract variants from EBI variation API response.
     * Each variant has .position (residue number).
     */
    extractVariants(variationData, amMap = null) {
        if (!variationData || !variationData.features) return [];

        const variants = [];
        variationData.features.forEach(f => {
            if (f.type !== 'VARIANT') return;
            const pos = parseInt(f.begin);
            if (isNaN(pos)) return;
            // Keep stop-gain / nonsense variants (alternativeSequence === '*'); only skip
            // entries with no substituted residue at all. AlphaMissense lookups for '*' simply
            // return null (it scores missense only), which the impact extractor handles.
            if (!f.alternativeSequence) return;

            const consequence = this._classifyConsequence(f);
            const provenance  = this._classifyProvenance(f);
            const diseaseData = this._extractDiseaseData(f);
            const impact = this._extractImpactMetadata(f, amMap);

            variants.push({
                position:         pos,
                wildType:         f.wildType || '?',
                mutant:           f.alternativeSequence,
                consequence:      consequence,
                consequenceColor: this.CONSEQUENCE_CATEGORIES[consequence]?.color || '#9e9e9e',
                provenance:       provenance,
                provenanceColor:  this.PROVENANCE_CATEGORIES[provenance]?.color || '#26c6da',
                description:      this._variantDesc(f),
                sourceType:       f.sourceType || 'unknown',
                diseases:         diseaseData.labels,
                diseaseIds:       diseaseData.ids,
                diseasePairs:     diseaseData.pairs,
                clinVarSignificance: impact.clinVarSignificance,
                clinVarReviewStatus: impact.clinVarReviewStatus,
                rsIds:            impact.rsIds,
                alphaMissenseScore: impact.alphaMissenseScore,
                gnomadAf:         impact.gnomadAf,
                genomicLocation:  impact.genomicLocation,
                xrefs:            f.xrefs || [],
            });
        });

        return variants;
    },

    _classifyConsequence(v) {
        const sigs = v.clinicalSignificances || [];
        for (const sig of sigs) {
            const t = (sig.type || '').toLowerCase();
            if (t.includes('pathogenic') && !t.includes('benign') && !t.includes('uncertain'))
                return 'Likely pathogenic or pathogenic';
        }
        for (const sig of sigs) {
            const t = (sig.type || '').toLowerCase();
            if (t.includes('benign') && !t.includes('pathogenic'))
                return 'Likely benign or benign';
        }
        for (const sig of sigs) {
            const t = (sig.type || '').toLowerCase();
            if (t.includes('uncertain'))
                return 'Uncertain significance';
        }
        const predictedDeleterious = (v.predictions || []).some(p => {
            const label = [
                p.predictionValType,
                p.predictionVal,
                p.value,
                p.name,
            ].filter(Boolean).join(' ').toLowerCase();
            return /\b(deleterious|damaging|pathogenic)\b/.test(label)
                && !/\b(benign|tolerated)\b/.test(label);
        });
        if (predictedDeleterious) return 'Predicted deleterious';
        return 'Uncertain significance';
    },

    _classifyProvenance(v) {
        const src = (v.sourceType || '').toLowerCase();
        if (src === 'uniprot' || src === 'mixed') return 'UniProt reviewed';
        const xrefs = v.xrefs || [];
        if (xrefs.some(x => (x.name || '').toLowerCase().includes('clinvar'))) return 'ClinVar';
        return 'Large scale studies';
    },

    _variantDesc(v) {
        if (v.descriptions && v.descriptions.length)
            return v.descriptions.map(d => d.value).join('; ');
        const sigs = v.clinicalSignificances || [];
        if (sigs.length) return sigs.map(s => s.type).join(', ');
        return '';
    },

    /**
     * Extract disease names from variant data.
     * Primary source: association[] array from the EBI variation API.
     * Returns:
     *   labels  – all disease display-labels (abbreviation or full name)
     *   ids     – all disease UniProt IDs (DI-xxxxx) from association entries
     *   pairs   – [{id, label}] with CORRECT 1-to-1 mapping per association entry
     */
    _extractDiseaseData(v) {
        const labels = new Set();
        const ids = new Set();
        const pairs = [];

        // 1) From association array – structured disease data with proper id↔label pairing
        if (v.association && Array.isArray(v.association)) {
            v.association.forEach(a => {
                const id = a.id || a.accession || '';
                if (id) ids.add(id);
                if (a.disease === true && a.name) {
                    const abbrMatch = a.name.match(/\(([A-Z][A-Z0-9]+)\)\s*$/);
                    const label = abbrMatch ? abbrMatch[1] : a.name;
                    labels.add(label);
                    pairs.push({ id: id || '', label }); // correct 1-to-1
                }
            });
        }

        // 2) From descriptions array – fallback labels only, no reliable IDs
        if (v.descriptions && Array.isArray(v.descriptions)) {
            v.descriptions.forEach(d => {
                if (!d.value) return;
                const abbrMatch = d.value.match(/^([A-Z][A-Z0-9]{1,15})(?:;|$)/);
                if (abbrMatch) {
                    const abbr = abbrMatch[1];
                    if (!['ECO', 'MIM', 'NCI', 'VAR', 'PRO'].includes(abbr) && abbr.length > 1) {
                        labels.add(abbr);
                    }
                }
            });
        }

        // 3) Last-resort fallback on flat description text
        if (labels.size === 0) {
            const desc = v.description || '';
            const matches = desc.match(/\bin\s+([A-Z][A-Z0-9]{1,10}(?:\s+and\s+[A-Z][A-Z0-9]{1,10})*)\b/g);
            if (matches) {
                matches.forEach(m => {
                    const parts = m.replace(/^in\s+/, '').split(/\s+and\s+/);
                    parts.forEach(p => {
                        const trimmed = p.trim();
                        if (trimmed && !['dbSNP'].includes(trimmed) && trimmed.length > 1) {
                            labels.add(trimmed);
                        }
                    });
                });
            }
        }

        return { labels: Array.from(labels), ids: Array.from(ids), pairs };
    },

    _extractImpactMetadata(v, amMap = null) {
        const xrefs = v.xrefs || [];
        const sigs = v.clinicalSignificances || [];
        const rsIds = new Set();
        xrefs.forEach(x => {
            const val = x.id || x.name || x.value || '';
            if (/^rs\d+$/i.test(val)) rsIds.add(val);
            if ((x.database || x.name || '').toLowerCase().includes('dbsnp') && val) rsIds.add(val);
        });

        let alphaMissenseScore = null;
        // Primary: look up per-variant score from AlphaMissense CSV (AF-{id}-F1-aa-substitutions.csv)
        if (amMap && v.wildType && v.alternativeSequence && v.alternativeSequence.length === 1) {
            const key = `${v.wildType}${v.begin}${v.alternativeSequence}`;
            const s = amMap.get(key);
            if (Number.isFinite(s)) alphaMissenseScore = s;
        }
        // Fallback: EBI variation API predictions (uses predAlgorithmNameType field)
        if (alphaMissenseScore === null) {
            (v.predictions || []).forEach(p => {
                const label = `${p.predAlgorithmNameType || p.predictionName || ''} ${p.name || ''}`.toLowerCase();
                const score = Number(p.score ?? p.predictionVal ?? p.value);
                if (label.includes('alphamissense') && Number.isFinite(score)) alphaMissenseScore = score;
            });
        }

        // gnomAD allele frequency + genomic HGVS — already present in the EBI variation payload.
        let gnomadAf = null;
        (v.populationFrequencies || []).forEach(pf => {
            if ((pf.source || '').toLowerCase().includes('gnomad') && Number.isFinite(pf.frequency)) {
                gnomadAf = gnomadAf === null ? pf.frequency : Math.max(gnomadAf, pf.frequency);
            }
        });
        const genomicLocation = Array.isArray(v.genomicLocation) ? (v.genomicLocation[0] || null) : (v.genomicLocation || null);

        return {
            clinVarSignificance: sigs.map(s => s.type).filter(Boolean).join(', '),
            clinVarReviewStatus: sigs.map(s => s.reviewStatus || s.review_status).filter(Boolean).join(', '),
            rsIds: Array.from(rsIds),
            alphaMissenseScore,
            gnomadAf,
            genomicLocation,
        };
    },

    /** Summary counts for consequence categories (skip 0-count) */
    getConsequenceSummary(variants) {
        const s = {};
        for (const cat of Object.keys(this.CONSEQUENCE_CATEGORIES))
            s[cat] = { count: 0, color: this.CONSEQUENCE_CATEGORIES[cat].color };
        variants.forEach(v => { if (s[v.consequence]) s[v.consequence].count++; });
        // Remove 0-count entries
        Object.keys(s).forEach(k => { if (s[k].count === 0) delete s[k]; });
        return s;
    },

    /** Summary counts for provenance categories (skip 0-count) */
    getProvenanceSummary(variants) {
        const s = {};
        for (const cat of Object.keys(this.PROVENANCE_CATEGORIES))
            s[cat] = { count: 0, color: this.PROVENANCE_CATEGORIES[cat].color };
        variants.forEach(v => { if (s[v.provenance]) s[v.provenance].count++; });
        Object.keys(s).forEach(k => { if (s[k].count === 0) delete s[k]; });
        return s;
    },

    /** Enrich each variant with a stable diseaseColor (first disease's palette color, or grey). */
    computeDiseaseColors(variants) {
        const colorMap = new Map(); // disease name → color
        let idx = 0;
        variants.forEach(v => {
            (v.diseases || []).forEach(d => {
                if (!colorMap.has(d)) colorMap.set(d, DISEASE_PALETTE[idx++ % DISEASE_PALETTE.length]);
            });
        });
        variants.forEach(v => {
            const firstDisease = (v.diseases || [])[0];
            v.diseaseColor = firstDisease ? colorMap.get(firstDisease) : '#9e9e9e';
        });
        return colorMap;
    },

    /** Get disease summary from variants */
    getDiseaseSummary(variants) {
        const s = {};
        let colorIdx = 0;
        variants.forEach(v => {
            (v.diseases || []).forEach(d => {
                if (!s[d]) {
                    s[d] = { count: 0, color: DISEASE_PALETTE[colorIdx++ % DISEASE_PALETTE.length] };
                }
                s[d].count++;
            });
        });
        // Also count "Unclassified" variants (no disease)
        const noDisease = variants.filter(v => !v.diseases || v.diseases.length === 0).length;
        if (noDisease > 0) {
            s['Unclassified'] = { count: noDisease, color: '#9e9e9e' };
        }
        return s;
    },

    /** Filter variants by active sets (additive model — see modal.js toggleVariantFilter).
     *  Intersection across the three axes, with one key rule: consequence / provenance only block a
     *  variant whose value is a KNOWN category; a missing or uncategorised value is never filtered out
     *  (so disease variants with no consequence/provenance still appear). An empty disease set shows
     *  nothing (additive model). */
    filterVariants(variants, activeConsequences, activeProvenances, activeDiseases) {
        return variants.filter(v => {
            if (this.CONSEQUENCE_CATEGORIES[v.consequence] && !activeConsequences.has(v.consequence)) return false;
            if (this.PROVENANCE_CATEGORIES[v.provenance] && !activeProvenances.has(v.provenance)) return false;
            if (activeDiseases) {
                return (v.diseases || []).some(d => activeDiseases.has(d));
            }
            return true;
        });
    }
};
