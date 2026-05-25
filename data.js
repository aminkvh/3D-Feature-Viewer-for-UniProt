/* ============================================
   Data Processing — Chrome Extension
   ============================================ */

const DataProcessor = {

    // ---- PTM Colors by category ----
    PTM_COLORS: {
        'Phosphorylation':       '#ff6d00',
        'N6-acetyllysine':       '#00e5ff',
        'N6-succinyllysine':     '#00acc1',
        'N6-methyllysine':       '#76ff03',
        'N6,N6-dimethyllysine':  '#64dd17',
        'Omega-N-methylarginine':'#b2ff59',
        'Symmetric dimethylarginine': '#69f0ae',
        'Asymmetric dimethylarginine': '#00e676',
        'Glycyl lysine isopeptide': '#e040fb',
        'Ubiquitination':        '#ea80fc',
        'Ubiquitinated lysine':  '#ea80fc',
        'SUMOylation':           '#d500f9',
        'ADP-ribosylation':      '#ff4081',
        'Glycosylation':         '#b87800',
        'Disulfide bond':        '#26a69a',
        'Lipidation':            '#f48fb1',
        'Cross-link':            '#ce93d8',
        'Acetylation':           '#00e5ff',
        'Methylation':           '#76ff03',
        'Citrulline':            '#7c4dff',
        'Hydroxyproline':        '#00b0ff',
        'Hydroxylation':         '#00b0ff',
        'Nitration':             '#ff6e40',
        'Pyroglutamic acid':     '#ffd740',
        'S-nitrosocysteine':     '#ff5252',
        'Deamidation':           '#448aff',
        'Modified residue':      '#78909c',
        'default':               '#ffab40',
    },

    // ---- Consequence categories ----
    CONSEQUENCE_CATEGORIES: {
        'Likely pathogenic or pathogenic': { color: '#ef5350' },
        'Predicted deleterious':           { color: '#ffa726' },
        'Likely benign or benign':         { color: '#66bb6a' },
        'Uncertain significance':          { color: '#9e9e9e' },
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
            if (!f.alternativeSequence || f.alternativeSequence === '*') return;

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
        if (v.predictions && v.predictions.length > 0) return 'Predicted deleterious';
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

        return {
            clinVarSignificance: sigs.map(s => s.type).filter(Boolean).join(', '),
            clinVarReviewStatus: sigs.map(s => s.reviewStatus || s.review_status).filter(Boolean).join(', '),
            rsIds: Array.from(rsIds),
            alphaMissenseScore,
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

    /** Get disease summary from variants */
    getDiseaseSummary(variants) {
        const DISEASE_COLORS = [
            '#ef5350', '#42a5f5', '#ab47bc', '#66bb6a', '#ffa726',
            '#26c6da', '#ec407a', '#7e57c2', '#5c6bc0', '#29b6f6',
            '#8d6e63', '#78909c', '#d4e157', '#ff7043',
        ];
        const s = {};
        let colorIdx = 0;
        variants.forEach(v => {
            (v.diseases || []).forEach(d => {
                if (!s[d]) {
                    s[d] = { count: 0, color: DISEASE_COLORS[colorIdx % DISEASE_COLORS.length] };
                    colorIdx++;
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

    /** Filter variants by active sets */
    filterVariants(variants, activeConsequences, activeProvenances, activeDiseases) {
        return variants.filter(v => {
            if (!activeConsequences.has(v.consequence)) return false;
            if (!activeProvenances.has(v.provenance)) return false;
            // Disease filter: if activeDiseases is provided, check it
            if (activeDiseases) {
                const vDiseases = v.diseases || [];
                if (vDiseases.length === 0) {
                    if (!activeDiseases.has('Unclassified')) return false;
                } else {
                    if (!vDiseases.some(d => activeDiseases.has(d))) return false;
                }
            }
            return true;
        });
    }
};
