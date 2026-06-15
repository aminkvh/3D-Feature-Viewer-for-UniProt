"""
3D Feature Viewer for UniProt — PyMOL plugin
============================================

Brings UniProt residue-level annotations onto a structure or trajectory already
open in PyMOL, or downloads the AlphaFold model directly, and projects PTMs,
disease variants, ClinVar/AlphaMissense, functional sites, domains/regions, and
membrane topology onto it — the desktop counterpart of the browser extension.

It pulls from the same public resources (UniProtKB, EBI Proteins API, AlphaFold DB,
PDBe/SIFTS) with no API keys and only the Python standard library.

Quick start (PyMOL command line)
--------------------------------
    run ufv_pymol.py                 # load the plugin (or install via Plugin Manager)

    ufv_gui                          # open the control panel (Fetch -> pick structure -> layers)
    ufv_load P35498                  # or just quick-load the AlphaFold model (no auto-projection)

Annotate a structure / trajectory YOU already loaded
----------------------------------------------------
    load mytraj.pdb, traj            # your object (single PDB or a trajectory)
    ufv_fetch P35498                 # fetch the annotations once

    # Tell the plugin how `traj` residue numbers relate to UniProt positions:
    ufv_map traj, P35498, identity              # resi == UniProt position (AlphaFold-style)
    ufv_map traj, P35498, sifts, 7dtd           # map through PDBe/SIFTS for PDB 7DTD
    # ...or define each chain by hand (great for trajectories with custom numbering):
    ufv_chain traj, A, 1                          # chain A resi 1  == UniProt 1
    ufv_chain traj, A, 200, 5, 480                # chain A resi 5  == UniProt 200, valid 5..480

    ufv_ptms traj            # PTM Ca spheres, coloured by category
    ufv_variants traj        # variant spheres, coloured by clinical consequence
    ufv_sites traj           # active / binding / metal sites (amber)
    ufv_domains traj         # colour the cartoon by domain / region / repeat
    ufv_topology traj        # colour by membrane topology
    ufv_alphamissense traj   # colour by mean AlphaMissense pathogenicity
    ufv_clear traj           # remove all UFV colouring/representations

All commands default to the active object when the object name is omitted.
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import re
import sys
import tempfile
import threading
import urllib.request

# This file is dual-purpose: a PyMOL plugin (when PyMOL is present) AND a standalone
# command-line backend used by the VMD plugin (ufv_vmd.tcl) for fetching/mapping. The
# pure data layer below has no PyMOL dependency.
try:
    from pymol import cmd
    _HAS_PYMOL = True
except Exception:  # running outside PyMOL (e.g. as the VMD backend CLI)
    cmd = None
    _HAS_PYMOL = False

# ----------------------------------------------------------------------------------------------
# Endpoints (mirror api.js)
# ----------------------------------------------------------------------------------------------
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{}"
AF_MODEL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v{}.pdb"
FEATURES_URL = "https://www.ebi.ac.uk/proteins/api/features/{}"
VARIATION_URL = "https://www.ebi.ac.uk/proteins/api/variation/{}"
PROTEOMICS_PTM_URL = "https://www.ebi.ac.uk/proteins/api/proteomics-ptm/{}"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.json"
AM_CSV_URL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-aa-substitutions.csv"
PDBe_SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{}"
PDBe_BEST = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{}"
PDBe_PDB = "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{}.ent"
BEACONS_SUMMARY = "https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/{}.json"
LIGAND_CCD = "https://data.rcsb.org/rest/v1/core/chemcomp/{}"
PUBCHEM_FP = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/property/Fingerprint2D/JSON"

_USER_AGENT = "UFV-PyMOL-plugin/1.0 (https://github.com/aminkvh/3D-Feature-Viewer-for-UniProt)"

# ----------------------------------------------------------------------------------------------
# Colour tables (mirror data.js)
# ----------------------------------------------------------------------------------------------
PTM_COLORS = {
    "Phosphorylation": "#ff6d00", "N6-acetyllysine": "#00e5ff", "N6-succinyllysine": "#00acc1",
    "N6-methyllysine": "#76ff03", "N6,N6-dimethyllysine": "#64dd17", "Omega-N-methylarginine": "#b2ff59",
    "Symmetric dimethylarginine": "#69f0ae", "Asymmetric dimethylarginine": "#00e676",
    "Glycyl lysine isopeptide": "#e040fb", "Ubiquitination": "#ea80fc", "Ubiquitinated lysine": "#ea80fc",
    "SUMOylation": "#d500f9", "ADP-ribosylation": "#ff4081", "Glycosylation": "#b87800",
    "Disulfide bond": "#26a69a", "Lipidation": "#f48fb1", "Cross-link": "#ce93d8",
    "Acetylation": "#00e5ff", "Methylation": "#76ff03", "Citrulline": "#7c4dff",
    "Hydroxyproline": "#00b0ff", "Hydroxylation": "#00b0ff", "Nitration": "#ff6e40",
    "Pyroglutamic acid": "#ffd740", "S-nitrosocysteine": "#ff5252", "Deamidation": "#448aff",
    "Modified residue": "#78909c", "default": "#ffab40",
}
CONSEQUENCE_COLORS = {
    "Likely pathogenic or pathogenic": "#ef5350",
    "Predicted deleterious": "#ffa726",
    "Likely benign or benign": "#66bb6a",
    "Uncertain significance": "#9e9e9e",
}
TOPOLOGY_COLORS = {
    "transmembrane": "#f9a825", "intramembrane": "#ef6c00", "cytoplasmic": "#1e88e5",
    "extracellular": "#e53935", "lumenal": "#43a047", "periplasmic": "#26a69a",
    "nuclear": "#8e24aa", "mitochondrial": "#00897b", "other": "#7e57c2",
}
DOMAIN_PALETTE = ["#1e88e5", "#43a047", "#fb8c00", "#8e24aa", "#00897b", "#e53935", "#3949ab",
                  "#c0ca33", "#6d4c41", "#00acc1", "#d81b60", "#5e35b1", "#7cb342", "#f4511e",
                  "#039be5", "#fdd835"]
SITE_COLOR = "#fbc02d"
DOMAIN_TYPES = ["Domain", "Region", "Repeat", "Compositional bias", "Zinc finger",
                "Coiled coil", "Motif", "DNA binding"]
PTM_TYPES = {"MOD_RES", "CROSSLNK", "LIPID", "CARBOHYD", "DISULFID"}
STRUCTURE_LIST_CAP = 60  # max experimental PDB chains shown in the structure selector

# ----------------------------------------------------------------------------------------------
# Module state
# ----------------------------------------------------------------------------------------------
_CACHE = {}          # uniprot_id -> annotations dict
UFV_MAPS = {}        # object name -> mapping descriptor
_STATE = {"uid": None}


# ----------------------------------------------------------------------------------------------
# HTTP helpers
# ----------------------------------------------------------------------------------------------
def _http_get(url):
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT, "Accept": "*/*"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def _get_json(url):
    try:
        return json.loads(_http_get(url).decode("utf-8"))
    except Exception:
        return None


def _get_text(url):
    try:
        return _http_get(url).decode("utf-8", "replace")
    except Exception:
        return None


# ----------------------------------------------------------------------------------------------
# Data extraction (mirror data.js)
# ----------------------------------------------------------------------------------------------
def _categorize_ptm(desc, type_):
    d = (desc or "").lower()
    if type_ == "DISULFID":
        return "Disulfide bond"
    if type_ == "LIPID":
        return "Lipidation"
    if type_ == "CARBOHYD":
        return "Glycosylation"
    if type_ == "CROSSLNK":
        if "ubiquitin" in d:
            return "Ubiquitination"
        if "sumo" in d:
            return "SUMOylation"
        return "Cross-link"
    if "phospho" in d:
        return "Phosphorylation"
    if "n6-acetyl" in d:
        return "N6-acetyllysine"
    if "n6-succinyl" in d:
        return "N6-succinyllysine"
    if "n6,n6-dimethyl" in d:
        return "N6,N6-dimethyllysine"
    if "n6-methyl" in d:
        return "N6-methyllysine"
    if "omega-n-methyl" in d:
        return "Omega-N-methylarginine"
    if "symmetric dimethyl" in d:
        return "Symmetric dimethylarginine"
    if "asymmetric dimethyl" in d:
        return "Asymmetric dimethylarginine"
    if "adp-ribos" in d:
        return "ADP-ribosylation"
    if "ubiquitin" in d:
        return "Ubiquitinated lysine"
    if "glycyl lysine isopeptide" in d:
        return "Glycyl lysine isopeptide"
    if "citrulline" in d:
        return "Citrulline"
    if "hydroxyproline" in d:
        return "Hydroxyproline"
    if "hydroxyl" in d:
        return "Hydroxylation"
    if "nitrated" in d:
        return "Nitration"
    if "pyroglutam" in d:
        return "Pyroglutamic acid"
    if "s-nitrosocysteine" in d:
        return "S-nitrosocysteine"
    if "deamidated" in d:
        return "Deamidation"
    label = (desc or "").split(";")[0].strip()
    return label or "Modified residue"


def _extract_ptms(features):
    out = []
    for f in (features or {}).get("features", []):
        if f.get("type") not in PTM_TYPES:
            continue
        try:
            pos = int(f.get("begin"))
        except (TypeError, ValueError):
            continue
        cat = _categorize_ptm(f.get("description", ""), f.get("type"))
        out.append({
            "position": pos,
            "endPosition": int(f.get("end") or pos),
            "description": f.get("description") or "Unknown modification",
            "category": cat,
            "color": PTM_COLORS.get(cat, PTM_COLORS["default"]),
        })
    return out


def _classify_consequence(v):
    sigs = v.get("clinicalSignificances") or []
    for s in sigs:
        t = (s.get("type") or "").lower()
        if "pathogenic" in t and "benign" not in t and "uncertain" not in t:
            return "Likely pathogenic or pathogenic"
    for s in sigs:
        t = (s.get("type") or "").lower()
        if "benign" in t and "pathogenic" not in t:
            return "Likely benign or benign"
    for s in sigs:
        if "uncertain" in (s.get("type") or "").lower():
            return "Uncertain significance"
    for p in (v.get("predictions") or []):
        label = " ".join(str(p.get(k, "")) for k in
                         ("predictionValType", "predictionVal", "value", "name")).lower()
        if re.search(r"\b(deleterious|damaging|pathogenic)\b", label) and \
           not re.search(r"\b(benign|tolerated)\b", label):
            return "Predicted deleterious"
    return "Uncertain significance"


def _classify_provenance(v):
    src = (v.get("sourceType") or "").lower()
    if src in ("uniprot", "mixed"):
        return "UniProt reviewed"
    for x in (v.get("xrefs") or []):
        if "clinvar" in (x.get("name") or "").lower():
            return "ClinVar"
    return "Large scale studies"


def _disease_labels(v):
    labels = set()
    for a in (v.get("association") or []):
        if a.get("disease") is True and a.get("name"):
            m = re.search(r"\(([A-Z][A-Z0-9]+)\)\s*$", a["name"])
            labels.add(m.group(1) if m else a["name"])
    return sorted(labels)


def _extract_variants(variation, am_map):
    out = []
    for f in (variation or {}).get("features", []):
        if f.get("type") != "VARIANT":
            continue
        try:
            pos = int(f.get("begin"))
        except (TypeError, ValueError):
            continue
        if not f.get("alternativeSequence"):
            continue
        cons = _classify_consequence(f)
        am = None
        if am_map and f.get("wildType") and len(f.get("alternativeSequence", "")) == 1:
            am = am_map.get("{}{}{}".format(f["wildType"], f["begin"], f["alternativeSequence"]))
        sigs = f.get("clinicalSignificances") or []
        rs = set()
        for x in (f.get("xrefs") or []):
            val = x.get("id") or x.get("name") or ""
            if re.match(r"^rs\d+$", str(val), re.I) or "dbsnp" in (x.get("name") or "").lower():
                if val:
                    rs.add(val)
        out.append({
            "position": pos,
            "wildType": f.get("wildType") or "?",
            "mutant": f.get("alternativeSequence"),
            "consequence": cons,
            "consequenceColor": CONSEQUENCE_COLORS.get(cons, "#9e9e9e"),
            "provenance": _classify_provenance(f),
            "reviewed": _classify_provenance(f) == "UniProt reviewed",
            "diseases": _disease_labels(f),
            "clinVar": ", ".join(s.get("type", "") for s in sigs if s.get("type")),
            "clinVarReview": ", ".join(filter(None, (s.get("reviewStatus") or s.get("review_status") or "" for s in sigs))),
            "rsIds": sorted(rs),
            "alphaMissense": am,
        })
    return out


def _extract_sites(features):
    labels = {"SITE": "Site", "ACT_SITE": "Active site", "BINDING": "Binding site", "METAL": "Metal binding"}
    out = []
    for f in (features or {}).get("features", []):
        lab = labels.get(f.get("type"))
        if not lab:
            continue
        try:
            pos = int(f.get("begin"))
        except (TypeError, ValueError):
            continue
        lig = (f.get("ligand") or {}).get("name") or ""
        if f.get("type") == "SITE":
            desc = f.get("description") or "Site"
        else:
            parts = [p for p in (lig, f.get("description")) if p]
            desc = "{}: {}".format(lab, " — ".join(parts)) if parts else lab
        out.append({"position": pos, "endPosition": int(f.get("end") or pos),
                    "description": desc, "color": SITE_COLOR})
    return out


def _extract_sites_uniprot(uniprot):
    labels = {"Site": "Site", "Active site": "Active site", "Binding site": "Binding site",
              "Metal binding": "Metal binding"}
    out = []
    for f in (uniprot or {}).get("features", []):
        lab = labels.get(f.get("type"))
        if not lab:
            continue
        pos = ((f.get("location") or {}).get("start") or {}).get("value")
        if pos is None:
            continue
        end = ((f.get("location") or {}).get("end") or {}).get("value") or pos
        lig = (f.get("ligand") or {}).get("name") or ""
        note = (f.get("ligand") or {}).get("note") or ""
        if f.get("type") == "Site":
            desc = f.get("description") or "Site"
        else:
            parts = [p for p in (lig, note, f.get("description")) if p]
            desc = "{}: {}".format(lab, "; ".join(parts)) if parts else lab
        out.append({"position": pos, "endPosition": end, "description": desc, "color": SITE_COLOR})
    return out


def _extract_topology(features):
    out = []
    for f in (features or {}).get("features", []):
        if f.get("type") not in ("TOPO_DOM", "TRANSMEM", "INTRAMEM"):
            continue
        try:
            start = int(f.get("begin"))
        except (TypeError, ValueError):
            continue
        end = int(f.get("end") or start)
        C = TOPOLOGY_COLORS
        if f.get("type") == "TRANSMEM":
            label, color = "Transmembrane", C["transmembrane"]
        elif f.get("type") == "INTRAMEM":
            label, color = "Intramembrane", C["intramembrane"]
        else:
            d = (f.get("description") or "").lower()
            label = f.get("description") or "Topological domain"
            color = (C["cytoplasmic"] if "cytoplasm" in d else
                     C["extracellular"] if "extracellular" in d else
                     C["lumenal"] if "lumen" in d else
                     C["periplasmic"] if "periplasm" in d else
                     C["nuclear"] if "nuclear" in d else
                     C["mitochondrial"] if "mitochond" in d else C["other"])
        out.append({"start": start, "end": end, "label": label, "color": color})
    return out


def _extract_domains(uniprot):
    out = []
    feats = [f for f in (uniprot or {}).get("features", []) if f.get("type") in DOMAIN_TYPES]
    for i, f in enumerate(feats):
        start = ((f.get("location") or {}).get("start") or {}).get("value")
        if start is None:
            continue
        end = ((f.get("location") or {}).get("end") or {}).get("value") or start
        desc = f.get("description") or f.get("type")
        out.append({"position": start, "endPosition": end, "type": f.get("type"),
                    "description": desc, "color": DOMAIN_PALETTE[i % len(DOMAIN_PALETTE)],
                    "isRange": end > start})
    return sorted(out, key=lambda d: d["position"])


def _parse_am_csv(text):
    """Parse the AlphaFold AlphaMissense substitution CSV into {'<wt><pos><mut>': score}."""
    out = {}
    if not text:
        return out
    for line in text.splitlines()[1:]:
        line = line.strip()
        if not line or "," not in line:
            continue
        var, rest = line.split(",", 1)
        try:
            out[var] = float(rest.split(",")[0])
        except ValueError:
            continue
    return out


def _am_mean_by_position(am_map):
    """Mean AlphaMissense score per residue position (across substitutions)."""
    acc = {}
    for key, score in am_map.items():
        m = re.match(r"^[A-Z](\d+)[A-Z]$", key)
        if not m:
            continue
        pos = int(m.group(1))
        s = acc.setdefault(pos, [0.0, 0])
        s[0] += score
        s[1] += 1
    return {pos: tot / n for pos, (tot, n) in acc.items() if n}


def _compute_burden(variants):
    """Mutation/phenotype burden (mirrors analysis.js): a residue is burden-positive when it has
    >=2 variant records, >=2 distinct disease/phenotype labels, and a composite of the within-
    protein ranks of (record count, distinct-disease count) in the top decile."""
    from collections import defaultdict
    rec = defaultdict(int)
    dis = defaultdict(set)
    for v in variants:
        rec[v["position"]] += 1
        for d in v.get("diseases", []):
            dis[v["position"]].add(d)
    positions = list(rec)
    if not positions:
        return set()

    def rank_map(values):
        order = sorted(set(values))
        return {val: i for i, val in enumerate(order)}
    cr = rank_map([rec[p] for p in positions])
    dr = rank_map([len(dis[p]) for p in positions])
    comp = {p: (cr[rec[p]] + dr[len(dis[p])]) / 2.0 for p in positions}
    cutoff = sorted(comp.values())[int(len(comp) * 0.9)] if comp else 0
    return {p for p in positions
            if rec[p] >= 2 and len(dis[p]) >= 2 and comp[p] >= cutoff}


# ----------------------------------------------------------------------------------------------
# Geometry-based analyses (ported from analysis.js). All operate on a `geom` list of modelled Ca
# atoms — each {uniPos, chain, resi, x, y, z} — so they are unit-testable without PyMOL. The
# permutation RNG is Python's Random(seed); tiers are statistically equivalent to the JS (not
# byte-identical, since the JS uses mulberry32).
# ----------------------------------------------------------------------------------------------
import math as _math
import random as _random

HOTSPOT_PERMUTATIONS = 1000
HOTSPOT_SEED = 0x9e3779b9 & 0xFFFFFFFF
HOTSPOT_MIN_PATH = 3
HOTSPOT_MIN_BENIGN = 2
HUB_Z_STRONG = 3.0
HUB_Z_MODERATE = 2.0
HUB_MAX_CA = 6000
_TIER_RANK = {"strong": 3, "moderate": 2, "weak": 1}


def _d2(a, b):
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return dx * dx + dy * dy + dz * dz


def _adjacency(coords, threshold):
    """Neighbour-index lists for points within `threshold` Å. Vectorised with numpy (block-wise) so
    the O(N²) distance scan that dominated the analyses runs in milliseconds instead of seconds;
    falls back to a pure-Python loop if numpy is unavailable."""
    n = len(coords)
    th2 = threshold * threshold
    adj = [[] for _ in range(n)]
    try:
        import numpy as np
        pts = np.asarray(coords, dtype=np.float64)
        block = 512
        for i0 in range(0, n, block):
            i1 = min(i0 + block, n)
            d2 = ((pts[i0:i1, None, :] - pts[None, :, :]) ** 2).sum(-1)
            for ii in range(i1 - i0):
                i = i0 + ii
                for j in np.nonzero(d2[ii] <= th2)[0].tolist():
                    if j > i:
                        adj[i].append(j)
                        adj[j].append(i)
        return adj
    except Exception:
        for i in range(n):
            ax, ay, az = coords[i]
            for j in range(i + 1, n):
                bx, by, bz = coords[j]
                if (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2 <= th2:
                    adj[i].append(j)
                    adj[j].append(i)
        return adj


def _ca_by_uni(geom):
    m = {}
    for g in geom:
        if g["uniPos"] is not None and g["uniPos"] not in m:
            m[g["uniPos"]] = g
    return m


def _is_pathogenic(v):
    s = (v.get("consequence") or v.get("clinVar") or "").lower()
    return "pathogenic" in s or "deleterious" in s


def _residue_neighborhood(geom, center_uni, threshold=8.0):
    """[(uniPos, distance)] for residues within `threshold` Å of center_uni's Cα, nearest first."""
    cab = _ca_by_uni(geom)
    c = cab.get(center_uni)
    if not c:
        return []
    th2 = threshold * threshold
    out = [(uni, _math.sqrt(_d2(c, g))) for uni, g in cab.items()
           if uni != center_uni and _d2(c, g) <= th2]
    out.sort(key=lambda t: t[1])
    return out


def _ptm_variant_proximity(ptms, variants, geom):
    """Port of analysis.js computePtmVariantProximity: per PTM, variants in tiers (same residue /
    <=8 Å pathogenic / <=12 Å), nearest variant, counts."""
    cab = _ca_by_uni(geom)
    vbypos = {}
    for v in variants:
        vbypos.setdefault(v["position"], []).append(v)
    result, seen = {}, set()
    for ptm in ptms:
        pp = ptm["position"]
        if pp in seen:
            continue
        seen.add(pp)
        pc = cab.get(pp)
        if not pc:
            continue
        tier1 = list(vbypos.get(pp, []))
        tier2, tier3 = [], []
        nearest, nearest_var = float("inf"), None
        for uni, g in cab.items():
            if uni == pp:
                continue
            d2 = _d2(pc, g)
            if d2 > 144.0:  # 12 Å
                continue
            vh = vbypos.get(uni)
            if not vh:
                continue
            dist = _math.sqrt(d2)
            for v in vh:
                (tier2 if (dist <= 8 and _is_pathogenic(v)) else tier3).append((v, dist))
                if dist < nearest:
                    nearest = dist
                    nearest_var = "%s%s%s" % (v.get("wildType", ""), v["position"], v.get("mutant", ""))
        if tier1 and nearest > 0:
            v0 = tier1[0]
            nearest, nearest_var = 0.0, "%s%s%s" % (v0.get("wildType", ""), pp, v0.get("mutant", ""))
        tier = 1 if tier1 else 2 if tier2 else 3 if tier3 else None
        if tier is not None:
            result[pp] = {
                "tier": tier, "tier1": tier1, "tier2": tier2, "tier3": tier3,
                "pathCount8A": len(tier2),
                "nearbyCount8A": len(tier1) + len(tier2) + sum(1 for _, d in tier3 if d <= 8),
                "nearestDist": None if nearest == float("inf") else nearest,
                "nearestVariant": nearest_var,
            }
    return result


HUB_SAMPLE_SOURCES = 400  # Brandes is O(V·E); sampling sources makes it O(k·E) with ~same z-scores


def _betweenness_hubs(geom, threshold=8.0):
    """Port of analysis.js betweennessHubs: Brandes' betweenness on the Cα contact graph, tiers by
    absolute z-score. For large chains betweenness is estimated from a random sample of source
    nodes (k=HUB_SAMPLE_SOURCES) — the tiering uses z-scores, which the sampling preserves — so a
    2000-residue chain drops from ~5 s to ~1 s. Returns {uniPos: 'strong'|'moderate'}."""
    nodes = [g for g in geom if g["uniPos"] is not None]
    n = len(nodes)
    if n < 8 or n > HUB_MAX_CA:
        return {}
    adj = _adjacency([(g["x"], g["y"], g["z"]) for g in nodes], threshold)
    bc = [0.0] * n
    if n > HUB_SAMPLE_SOURCES:
        sources = _random.Random(HOTSPOT_SEED).sample(range(n), HUB_SAMPLE_SOURCES)
    else:
        sources = range(n)
    for s in sources:
        stack, pred = [], [[] for _ in range(n)]
        sigma = [0.0] * n
        sigma[s] = 1.0
        dist = [-1] * n
        dist[s] = 0
        queue, qi = [s], 0
        while qi < len(queue):
            v = queue[qi]
            qi += 1
            stack.append(v)
            for w in adj[v]:
                if dist[w] < 0:
                    dist[w] = dist[v] + 1
                    queue.append(w)
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    pred[w].append(v)
        delta = [0.0] * n
        for idx in range(len(stack) - 1, -1, -1):
            w = stack[idx]
            for v in pred[w]:
                delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w])
            if w != s:
                bc[w] += delta[w]
    scale = 1.0 / ((n - 1) * (n - 2)) if n > 2 else 0.0
    bc = [b * scale for b in bc]
    mean = sum(bc) / n
    std = _math.sqrt(sum((b - mean) ** 2 for b in bc) / n)
    if std <= 0:
        return {}
    out = {}
    for i in range(n):
        if bc[i] <= 0:
            continue
        z = (bc[i] - mean) / std
        if z >= HUB_Z_STRONG:
            out[nodes[i]["uniPos"]] = "strong"
        elif z >= HUB_Z_MODERATE:
            out[nodes[i]["uniPos"]] = "moderate"
    return out


def _bh_fdr(centres):
    """In-place Benjamini-Hochberg q-values; `centres` sorted ascending by 'p'."""
    n = len(centres)
    for i, c in enumerate(centres):
        c["q"] = min(c["p"] * n / (i + 1), 1.0)
    for i in range(n - 2, -1, -1):
        centres[i]["q"] = min(centres[i]["q"], centres[i + 1]["q"])


def _compute_hotspots(geom, variants, threshold=8.0):
    """Port of analysis.js computeHotspots: case-control label-permutation when benign controls are
    available, else a spatial-placement null. Returns {uniPos: 'strong'|'moderate'|'weak'}."""
    path = {v["position"] for v in variants if _is_pathogenic(v)}
    benign = {v["position"] for v in variants
              if v["position"] not in path and "benign" in (v.get("consequence") or "").lower()}
    if len(path) < HOTSPOT_MIN_PATH:
        return {}
    cab = _ca_by_uni(geom)
    th2 = threshold * threshold
    rng = _random.Random(HOTSPOT_SEED)

    if len(benign) >= HOTSPOT_MIN_BENIGN:  # case-control
        items = [(p, True, cab[p]) for p in path if p in cab] + \
                [(p, False, cab[p]) for p in benign if p in cab]
        k = len(items)
        if k < HOTSPOT_MIN_PATH + 1:
            return {}
        adj = _adjacency([(it[2]["x"], it[2]["y"], it[2]["z"]) for it in items], threshold)
        neigh = [[i] + adj[i] for i in range(k)]  # include self, as in the JS
        labels = [it[1] for it in items]
        local = [sum(1 for j in neigh[i] if labels[j]) for i in range(k)]
        base = sum(1 for l in labels if l) / k
        ge = [0] * k
        perm = labels[:]
        for _ in range(HOTSPOT_PERMUTATIONS):
            rng.shuffle(perm)
            for i in range(k):
                if sum(1 for j in neigh[i] if perm[j]) >= local[i]:
                    ge[i] += 1
        centres = []
        for i in range(k):
            if not items[i][1]:
                continue
            tot = len(neigh[i])
            frac = local[i] / tot if tot else 0
            centres.append({"pos": items[i][0], "p": (ge[i] + 1) / (HOTSPOT_PERMUTATIONS + 1),
                            "local": local[i], "frac": frac, "er": frac / max(base, 1e-9)})
        centres.sort(key=lambda c: c["p"])
        _bh_fdr(centres)
        out = {}
        for c in centres:
            if c["q"] <= 0.10 and c["local"] >= 3 and c["er"] >= 2.0:
                out[c["pos"]] = "strong"
            elif c["q"] <= 0.25 and c["local"] >= 2 and c["er"] >= 1.5:
                out[c["pos"]] = "moderate"
            elif c["q"] <= 0.25 and c["local"] >= 2 and c["frac"] >= 0.40:
                out[c["pos"]] = "weak"
        return out

    # spatial-placement null (no benign controls needed)
    universe = [g for g in geom if g["uniPos"] is not None]
    n = len(universe)
    if n < 10 or n > 6000:
        return {}
    idx_by_uni = {u["uniPos"]: i for i, u in enumerate(universe)}
    path_idx = [idx_by_uni[p] for p in path if p in idx_by_uni]
    m = len(path_idx)
    if m < HOTSPOT_MIN_PATH:
        return {}
    neigh = _adjacency([(u["x"], u["y"], u["z"]) for u in universe], threshold)
    is_path = [0] * n
    for i in path_idx:
        is_path[i] = 1
    obs = [sum(1 for j in neigh[i] if is_path[j]) for i in range(n)]
    ge = {ci: 0 for ci in path_idx}
    pool = list(range(n))
    for _ in range(HOTSPOT_PERMUTATIONS):
        sim = [0] * n
        for i in range(m):
            r = i + rng.randrange(n - i)
            pool[i], pool[r] = pool[r], pool[i]
            sim[pool[i]] = 1
        for ci in path_idx:
            if sum(1 for j in neigh[ci] if sim[j]) >= obs[ci]:
                ge[ci] += 1
    centres = [{"ci": ci, "p": (ge[ci] + 1) / (HOTSPOT_PERMUTATIONS + 1), "den": obs[ci]} for ci in path_idx]
    centres.sort(key=lambda c: c["p"])
    _bh_fdr(centres)
    out = {}
    for c in centres:
        uni = universe[c["ci"]]["uniPos"]
        if c["q"] <= 0.10 and c["den"] >= 3:
            out[uni] = "strong"
        elif c["q"] <= 0.25 and c["den"] >= 2:
            out[uni] = "moderate"
    return out


def _per_chain_merge(geom, fn):
    """Run a per-chain analysis `fn(geom_of_chain)` and merge tiers across chains (strongest wins)."""
    chains = {}
    for g in geom:
        chains.setdefault(g["chain"], []).append(g)
    merged = {}
    for gs in chains.values():
        for uni, tier in fn(gs).items():
            if uni is None:
                continue
            if uni not in merged or _TIER_RANK.get(tier, 0) > _TIER_RANK.get(merged[uni], 0):
                merged[uni] = tier
    return merged


def _am_by_pos(am_map):
    """Index AlphaMissense by position once: {pos: [(wt, mut, score)] highest-first}. Avoids
    rescanning the whole substitution map (tens of thousands of keys) on every residue report."""
    out = {}
    for key, score in am_map.items():
        mt = re.match(r"^([A-Z])(\d+)([A-Z])$", key)
        if mt:
            out.setdefault(int(mt.group(2)), []).append((mt.group(1), mt.group(3), score))
    for p in out:
        out[p].sort(key=lambda t: -t[2])
    return out


def _am_profile(am_map, pos):
    """Per-substitution AlphaMissense scores at a position: [(wt, mut, score)], highest first."""
    return _am_by_pos(am_map).get(pos, [])


# Constraint-pocket (port of algorithms.js UFVPocket.computePockets). Getis-Ord Gi* on the
# residual of per-residue AlphaMissense mean after regressing out structural burial (coordination
# number) — residues more evolutionarily constrained than their burial explains, clustered in 3D.
PRISM_HSE = 13.0     # Å coordination-sphere radius (burial)
PRISM_CUT = 13.0     # Å spatial-weight neighbourhood
PRISM_SIGMA = 5.0    # Å Gaussian bandwidth
PRISM_MAXN = 3000


def _perm_count(n):
    return 999 if n <= 600 else 599 if n <= 1200 else 399


def _loess(X, Y, np, f=0.3):
    n = len(X)
    k = max(2, min(n, int(round(f * n))))
    pred = np.empty(n)
    for i in range(n):
        d = np.abs(X - X[i])
        h = np.partition(d, min(k - 1, n - 1))[min(k - 1, n - 1)] or 1e-9
        u = d / h
        w = np.where(u < 1, (1 - u ** 3) ** 3, 0.0)
        sw = w.sum(); swx = (w * X).sum(); swy = (w * Y).sum()
        swxx = (w * X * X).sum(); swxy = (w * X * Y).sum()
        denom = sw * swxx - swx * swx
        if abs(denom) < 1e-12:
            pred[i] = swy / sw if sw else Y[i]
        else:
            b = (sw * swxy - swx * swy) / denom
            pred[i] = (swy - b * swx) / sw + b * X[i]
    return pred


def _compute_pockets(geom, am_mean, q_threshold=0.10):
    """{uniPos: 'pocket'|'exposed'} for constraint-pocket candidates with BH-FDR q <= threshold.
    Burial context is the whole model (all chains); scored residues are those with an AlphaMissense
    mean. PAE is not used (pure Euclidean weights), unlike the web app's optional PAE gating."""
    try:
        import numpy as np
    except Exception:
        return {}
    C = np.array([(g["x"], g["y"], g["z"]) for g in geom], float)
    recs = [(i, g) for i, g in enumerate(geom) if g["uniPos"] is not None and g["uniPos"] in am_mean]
    n = len(recs)
    if n < 12 or n > PRISM_MAXN:
        return {}
    ctx_idx = np.array([i for i, _ in recs])
    P = C[ctx_idx]
    am = np.array([am_mean[g["uniPos"]] for _, g in recs], float)

    by_cr = {(g["chain"], g["resi"]): k for k, (_, g) in enumerate(recs)}
    dirs = np.zeros((n, 3))
    for k, (_, g) in enumerate(recs):
        pv = by_cr.get((g["chain"], g["resi"] - 1))
        nx = by_cr.get((g["chain"], g["resi"] + 1))
        v = None
        if pv is not None and nx is not None:
            v = 2 * P[k] - P[pv] - P[nx]
        elif pv is not None:
            v = P[k] - P[pv]
        elif nx is not None:
            v = P[k] - P[nx]
        if v is not None:
            L = float(np.linalg.norm(v)) or 1.0
            dirs[k] = v / L

    hse2 = PRISM_HSE ** 2
    cn = np.zeros(n); hse_up = np.zeros(n); hse_down = np.zeros(n)
    blk = 256
    for i0 in range(0, n, blk):
        i1 = min(i0 + blk, n)
        d = C[None, :, :] - P[i0:i1, None, :]           # (b, nctx, 3)
        d2 = (d ** 2).sum(-1)
        within = d2 <= hse2
        for ii in range(i1 - i0):
            k = i0 + ii
            within[ii, ctx_idx[k]] = False
            cn[k] = within[ii].sum()
            if dirs[k].any():
                proj = (d[ii] * dirs[k]).sum(-1)
                hse_up[k] = (within[ii] & (proj > 0)).sum()
                hse_down[k] = (within[ii] & (proj <= 0)).sum()

    cut2 = PRISM_CUT ** 2
    two_sig2 = 2 * PRISM_SIGMA ** 2
    W = np.zeros((n, n))
    for i0 in range(0, n, blk):
        i1 = min(i0 + blk, n)
        d2 = ((P[i0:i1, None, :] - P[None, :, :]) ** 2).sum(-1)
        W[i0:i1] = np.where(d2 <= cut2, np.exp(-d2 / two_sig2), 0.0)
    np.fill_diagonal(W, 0.0)
    W[W < 1e-6] = 0.0

    rp = am - _loess(cn, am, np, 0.3)
    M = W + np.eye(n)                 # Gi*_i = rp[i] + Σ w_ij rp[j]
    gi_obs = M.dot(rp)

    perm = _perm_count(n)
    rng = np.random.RandomState(HOTSPOT_SEED)
    Rp = np.empty((n, perm))
    for p in range(perm):
        Rp[:, p] = rng.permutation(rp)
    ge = (M.dot(Rp) >= gi_obs[:, None]).sum(1)

    med_cn = float(np.median(cn))
    best = {}  # uniPos -> (idx, p)
    for k, (_, g) in enumerate(recs):
        if rp[k] <= 0:
            continue
        p = (ge[k] + 1) / (perm + 1)
        u = g["uniPos"]
        if u not in best or p < best[u][1]:
            best[u] = (k, p)
    cands = sorted(best.items(), key=lambda kv: kv[1][1])  # by p
    nc = len(cands)
    qs = {}
    for rank, (u, (k, p)) in enumerate(cands):
        qs[u] = min(p * nc / (rank + 1), 1.0)
    keys = [u for u, _ in cands]
    for r in range(nc - 2, -1, -1):
        qs[keys[r]] = min(qs[keys[r]], qs[keys[r + 1]])

    out = {}
    for u, (k, p) in cands:
        if qs[u] <= q_threshold:
            concavity = hse_up[k] / (hse_down[k] + 1)
            out[u] = "pocket" if (cn[k] >= med_cn or concavity >= 1.3) else "exposed"
    return out


# ----------------------------------------------------------------------------------------------
# Ligand chemistry + Tanimoto similarity (ported from api.js getLigandInfo / getLigandFingerprint)
# ----------------------------------------------------------------------------------------------
_LIG_CACHE = {}
_FP_CACHE = {}


def ligand_info(ccd):
    """RCSB Chemical Component Dictionary entry for a CCD code: name/formula/SMILES/InChIKey/DrugBank."""
    key = str(ccd).upper()
    if key in _LIG_CACHE:
        return _LIG_CACHE[key]
    import urllib.parse
    j = _get_json(LIGAND_CCD.format(urllib.parse.quote(key))) or {}
    d = j.get("rcsb_chem_comp_descriptor") or {}
    db = next((x for x in (j.get("rcsb_chem_comp_related") or []) if x.get("resource_name") == "DrugBank"), None)
    info = {"id": key, "name": (j.get("chem_comp") or {}).get("name"),
            "formula": (j.get("chem_comp") or {}).get("formula"),
            "smiles": d.get("SMILES_stereo") or d.get("SMILES"), "inchikey": d.get("InChIKey"),
            "drugbank": (db or {}).get("resource_accession_code")}
    _LIG_CACHE[key] = info
    return info


def ligand_fingerprint(inchikey):
    """PubChem 881-bit CACTVS 2D substructure fingerprint (decoded) for an InChIKey, or None."""
    if not inchikey:
        return None
    if inchikey in _FP_CACHE:
        return _FP_CACHE[inchikey]
    import base64
    import urllib.parse
    j = _get_json(PUBCHEM_FP.format(urllib.parse.quote(inchikey)))
    b64 = (((j or {}).get("PropertyTable") or {}).get("Properties") or [{}])[0].get("Fingerprint2D")
    fp = None
    if b64:
        raw = base64.b64decode(b64)  # first 4 bytes are a length prefix, then 881 bits
        fp = bytes(((raw[4 + (i >> 3)] if 4 + (i >> 3) < len(raw) else 0) >> (7 - (i & 7))) & 1 for i in range(881))
    _FP_CACHE[inchikey] = fp
    return fp


def tanimoto(a, b):
    """Tanimoto coefficient between two bit vectors (bytes of 0/1), or None if either is missing."""
    if not a or not b:
        return None
    inter = union = 0
    for x, y in zip(a, b):
        if x and y:
            inter += 1
        if x or y:
            union += 1
    return inter / union if union else 0.0


# ----------------------------------------------------------------------------------------------
# Fetch orchestration
# ----------------------------------------------------------------------------------------------
def fetch_annotations(uid, force=False):
    """Fetch and cache every annotation layer for a UniProt accession."""
    uid = uid.strip().upper()
    if not force and uid in _CACHE:
        return _CACHE[uid]
    print("[UFV] fetching annotations for %s ..." % uid)
    # Fetch the independent endpoints in parallel — the AlphaMissense CSV is large, so doing them
    # concurrently cuts wall time from the sum to the slowest single request.
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as ex:
        f_features = ex.submit(_get_json, FEATURES_URL.format(uid))
        f_variation = ex.submit(_get_json, VARIATION_URL.format(uid))
        f_uniprot = ex.submit(_get_json, UNIPROT_URL.format(uid))
        f_am = ex.submit(_get_text, AM_CSV_URL.format(uid))
        features = f_features.result()
        variation = f_variation.result()
        uniprot = f_uniprot.result()
        am_map = _parse_am_csv(f_am.result())

    sites_uni = _extract_sites_uniprot(uniprot)
    sites_api = _extract_sites(features)
    seen = {}
    for s in sites_uni + sites_api:
        seen.setdefault((s["position"], s["endPosition"]), s)
    sites = sorted(seen.values(), key=lambda s: s["position"])

    ann = {
        "uid": uid,
        "sequence": ((uniprot or {}).get("sequence") or {}).get("value", ""),
        "ptms": _extract_ptms(features),
        "variants": _extract_variants(variation, am_map),
        "sites": sites,
        "topology": _extract_topology(features),
        "domains": _extract_domains(uniprot),
        "amMean": _am_mean_by_position(am_map),
        "amMap": am_map,
        "amByPos": _am_by_pos(am_map),
    }
    ann["burden"] = _compute_burden(ann["variants"])
    _CACHE[uid] = ann
    _STATE["uid"] = uid
    print("[UFV] %s: %d PTMs, %d variants, %d sites, %d domains, %d topology segments, AM=%s"
          % (uid, len(ann["ptms"]), len(ann["variants"]), len(ann["sites"]),
             len(ann["domains"]), len(ann["topology"]), "yes" if ann["amMean"] else "no"))
    return ann


def _alphafold_url(uid):
    data = _get_json(ALPHAFOLD_API.format(uid))
    if isinstance(data, list):
        for e in data:
            if e.get("pdbUrl"):
                return e["pdbUrl"]
    for ver in (6, 5, 4, 3):
        url = AF_MODEL.format(uid, ver)
        try:
            req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": _USER_AGENT})
            with urllib.request.urlopen(req, timeout=30) as r:
                if r.status == 200:
                    return url
        except Exception:
            continue
    return None


def _sifts_segments(pdb_id, uid):
    """{chain: [{u_start, u_end, pdb_start, pdb_end}]} from PDBe/SIFTS for this accession."""
    data = _get_json(PDBe_SIFTS.format(pdb_id.lower()))
    uni = (((data or {}).get(pdb_id.lower()) or {}).get("UniProt") or {})
    entry = uni.get(uid) or uni.get(uid.split("-")[0])
    if not entry:
        for k in uni:
            if k.upper() == uid.upper():
                entry = uni[k]
                break
    if not entry:
        return {}
    segs = {}
    for m in entry.get("mappings", []):
        ch = m.get("chain_id")
        start = m.get("start") or {}
        end = m.get("end") or {}
        pdb_start = start.get("author_residue_number")
        if pdb_start is None:
            pdb_start = start.get("residue_number", m.get("unp_start"))
        pdb_end = end.get("author_residue_number")
        if pdb_end is None:
            pdb_end = end.get("residue_number", m.get("unp_end"))
        segs.setdefault(ch, []).append({
            "u_start": m.get("unp_start"), "u_end": m.get("unp_end"),
            "pdb_start": pdb_start, "pdb_end": pdb_end,
        })
    return segs


def _isoform_structures(uid):
    """Non-canonical isoform AlphaFold models (e.g. AF-P35498-2). Their AlphaFold model URL comes
    from each isoform's 3D-Beacons summary. Numbering follows the isoform sequence; canonical
    annotations are assumed identical except in spliced regions (no VSP remapping here)."""
    base = uid.split("-")[0]
    u = _get_json(UNIPROT_URL.format(base)) or {}
    iso_ids = []
    for c in (u.get("comments") or []):
        if c.get("commentType") == "ALTERNATIVE PRODUCTS":
            for iso in (c.get("isoforms") or []):
                for iid in (iso.get("isoformIds") or []):
                    if iid and iid not in (base, base + "-1"):
                        iso_ids.append(iid)
    out = []
    for iid in iso_ids:
        d = _get_json(BEACONS_SUMMARY.format(iid)) or {}
        af = next((sm for e in (d.get("structures") or []) for sm in [e.get("summary") or e]
                   if sm.get("model_url") and re.search("alphafold", sm.get("provider", ""), re.I)), None)
        if not af:
            continue
        fmt = str(af.get("model_format", "")).upper()
        out.append({"key": "AF-%s" % iid, "label": "AlphaFold isoform %s" % iid, "source": "AlphaFold",
                    "url": af["model_url"], "fmt": "pdb" if fmt == "PDB" else "mmcif",
                    "pdbId": None, "chainId": None, "numbering": "identity", "cov": None})
    return out


def list_structures(uid, seqlen=0):
    """Available structures for an accession (for the structure selector): the canonical AlphaFold
    model, non-canonical isoform AlphaFold models, experimental PDB chains (PDBe best_structures,
    numbered via SIFTS on load), and computed models (3D-Beacons). Mirrors api.js getStructures."""
    uid = uid.upper()
    out = []
    # AlphaFold URL, isoform models, PDBe best_structures and 3D-Beacons summary fetched concurrently.
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as ex:
        f_af = ex.submit(_alphafold_url, uid)
        f_iso = ex.submit(_isoform_structures, uid)
        f_best = ex.submit(_get_json, PDBe_BEST.format(uid))
        f_summ = ex.submit(_get_json, BEACONS_SUMMARY.format(uid))
        af = f_af.result()
        isoforms = f_iso.result()
        best = f_best.result() or {}
        summ_data = f_summ.result() or {}
    if af:
        out.append({"key": "AF-%s" % uid, "label": "AlphaFold — predicted, full length",
                    "source": "AlphaFold", "url": af, "fmt": "pdb",
                    "pdbId": None, "chainId": None, "numbering": "identity", "cov": 100.0})
    out.extend(isoforms)
    items = best.get(uid) or best.get(uid.upper()) or []
    seen = set()
    exp = []
    for it in items:
        pid = (it.get("pdb_id") or "").upper()
        ch = it.get("chain_id") or "A"
        if not pid or (pid, ch) in seen:
            continue
        seen.add((pid, ch))
        cov = it.get("coverage")
        covpct = round(cov * 100, 1) if cov is not None else None
        res = it.get("resolution")
        # Don't show the SIFTS-range coverage here (it overstates); the real modelled coverage is
        # computed and shown after the structure loads.
        exp.append({"key": "%s_%s" % (pid, ch),
                    "label": "%s chain %s — %s%s" % (pid, ch, it.get("experimental_method", "experimental"),
                                                     (" %.1f Å" % res) if res else ""),
                    "source": "PDB", "url": PDBe_PDB.format(pid.lower()), "fmt": "pdb",
                    "pdbId": pid, "chainId": ch, "numbering": "sifts", "cov": covpct,
                    "_res": res or 99.0})
    # Proteins like insulin map to hundreds of PDB chains; sort best-coverage/-resolution first
    # and cap the list so the selector stays usable.
    exp.sort(key=lambda s: (-(s["cov"] or 0), s["_res"]))
    out.extend(exp[:STRUCTURE_LIST_CAP])
    for e in (summ_data.get("structures") or []):
        sm = e.get("summary") or e
        murl = sm.get("model_url")
        if not murl:
            continue
        if "EXPERIMENTAL" in str(sm.get("model_category", "")).upper():
            continue
        prov = sm.get("provider", "Computed")
        if re.search("alphafold", prov, re.I) and \
           str(sm.get("model_identifier", "")).upper() == "AF-%s-F1" % uid:
            continue
        fmt = str(sm.get("model_format", "")).upper()
        if fmt not in ("PDB", "MMCIF"):
            continue
        cov = sm.get("coverage")
        out.append({"key": "BCN-%s" % (sm.get("model_identifier") or murl),
                    "label": "%s — computed model" % prov, "source": "Computed", "url": murl,
                    "fmt": "pdb" if fmt == "PDB" else "mmcif", "pdbId": None, "chainId": None,
                    "numbering": "identity", "cov": round(cov * 100, 1) if cov is not None else None})
    return out


# ----------------------------------------------------------------------------------------------
# Residue-number mapping (UniProt position -> object chain+resi)
# ----------------------------------------------------------------------------------------------
def _set_identity_map(obj):
    UFV_MAPS[obj] = {"mode": "identity", "chains": {}}


def _set_sifts_map(obj, pdb_id, uid):
    segs = _sifts_segments(pdb_id, uid)
    if not segs:
        print("[UFV] no SIFTS mapping found for %s / %s — falling back to identity." % (pdb_id, uid))
        _set_identity_map(obj)
        return False
    UFV_MAPS[obj] = {"mode": "sifts", "segments": segs}
    print("[UFV] %s mapped via SIFTS for %s: chains %s" % (obj, pdb_id.upper(), ", ".join(sorted(segs))))
    return True


def _ensure_manual_map(obj):
    mp = UFV_MAPS.get(obj)
    if not mp or mp.get("mode") != "manual":
        mp = {"mode": "manual", "chains": {}}
        UFV_MAPS[obj] = mp
    return mp


def _uni_to_resi(obj, chain, pos):
    """Map a UniProt position to the object's author residue number on `chain` (or None)."""
    mp = UFV_MAPS.get(obj)
    if not mp:
        return pos  # default identity
    mode = mp["mode"]
    if mode == "identity":
        return pos
    if mode == "manual":
        c = mp["chains"].get(chain)
        if not c:
            return None
        resi = c["resi_start"] + (pos - c["u_start"])
        if resi < c["resi_start"]:
            return None
        if c["resi_end"] is not None and resi > c["resi_end"]:
            return None
        return resi
    if mode == "sifts":
        for seg in mp["segments"].get(chain, []):
            if seg["u_start"] <= pos <= seg["u_end"]:
                return seg["pdb_start"] + (pos - seg["u_start"])
        return None
    return pos


def _mapped_chains(obj):
    mp = UFV_MAPS.get(obj)
    if not mp:
        return cmd.get_chains(obj)
    if mp["mode"] == "sifts":
        return list(mp["segments"].keys())
    if mp["mode"] == "manual":
        return list(mp["chains"].keys())
    return cmd.get_chains(obj) or [""]


def _resi_to_uni(obj, chain, resi):
    """Inverse of _uni_to_resi: object author residue -> UniProt position (or None)."""
    try:
        resi = int(resi)
    except (TypeError, ValueError):
        return None
    mp = UFV_MAPS.get(obj)
    if not mp or mp["mode"] == "identity":
        return resi
    if mp["mode"] == "manual":
        c = mp["chains"].get(chain)
        if not c:
            return None
        if resi < c["resi_start"]:
            return None
        if c["resi_end"] is not None and resi > c["resi_end"]:
            return None
        return c["u_start"] + (resi - c["resi_start"])
    if mp["mode"] == "sifts":
        for seg in mp["segments"].get(chain, []):
            lo, hi = sorted((seg["pdb_start"], seg["pdb_end"]))
            if lo <= resi <= hi:
                return seg["u_start"] + (resi - seg["pdb_start"])
        return None
    return resi


# ----------------------------------------------------------------------------------------------
# Projection helpers
# ----------------------------------------------------------------------------------------------
_COLOR_REGISTERED = set()


def _color_name(hex_str):
    h = hex_str.lstrip("#")
    name = "ufv_%s" % h
    if name not in _COLOR_REGISTERED:   # set_color is cheap but this is called per colour per redraw
        cmd.set_color(name, [int(h[0:2], 16) / 255.0, int(h[2:4], 16) / 255.0, int(h[4:6], 16) / 255.0])
        _COLOR_REGISTERED.add(name)
    return name


def _resolve_object(obj):
    if obj:
        return obj
    objs = cmd.get_object_list()
    for o in objs:
        if o in UFV_MAPS:
            return o
    return objs[-1] if objs else None


def _resolve_uid(uid):
    return (uid or _STATE.get("uid"))


def _collapse_resi(resis):
    """Collapse a sorted int list into PyMOL resi tokens with ranges: [1,2,3,5] -> '1-3+5'. Keeps
    selection strings short so PyMOL doesn't choke parsing thousands of '+'-separated residues
    (a 260-residue domain becomes one '100-360' token instead of 260)."""
    out, i, n = [], 0, len(resis)
    while i < n:
        j = i
        while j + 1 < n and resis[j + 1] == resis[j] + 1:
            j += 1
        out.append(str(resis[i]) if i == j else "%d-%d" % (resis[i], resis[j]))
        i = j + 1
    return "+".join(out)


def _sel_for_positions(obj, positions, ca_only=False, target=None):
    """Build a compact PyMOL selection covering the given UniProt positions. Numbering is resolved
    with `obj`'s map; the selection is rooted at `target` (defaults to `obj`) so it can address a
    copied sphere layer object that shares the same chain/resi numbering."""
    per_chain = {}
    for ch in _mapped_chains(obj):
        resis = []
        for p in positions:
            r = _uni_to_resi(obj, ch, p)
            if r is not None:
                resis.append(int(r))
        if resis:
            per_chain[ch] = sorted(set(resis))
    if not per_chain:
        return None
    clauses = []
    for ch, resis in per_chain.items():
        rsel = "resi " + _collapse_resi(resis)
        clauses.append("(chain {} and {})".format(ch, rsel) if ch else "({})".format(rsel))
    sel = "({}) and ({})".format(target or obj, " or ".join(clauses))
    if ca_only:
        sel += " and name CA"
    return sel


# Each point (sphere) layer is drawn on its OWN object (a copy of just the annotated Cα atoms, no
# cartoon) named ufv_<obj>_<tag>. This keeps the main structure's cartoon out of every sphere
# rebuild — large variant sets no longer crawl — and keeps the object panel to a few clean names.
_SPHERE_STATE = {}    # obj -> { tag -> {color_hex: [uniprot_positions]} }  (active sphere layers)
_CARTOON_LAYER = {}   # obj -> active cartoon-colouring mode (informational)
POINT_TAGS = {"ptms": "ptm", "variants": "var", "sites": "site"}
CARTOON_TAGS = ("domains", "topology", "alphamissense", "burden", "plddt", "bfactor")


def _cons_token(cons):
    c = (cons or "").lower()
    if "pathogenic" in c:
        return "pathogenic"
    if "benign" in c:
        return "benign"
    if "deleterious" in c:
        return "deleterious"
    return "uncertain"


@contextlib.contextmanager
def _batch():
    """Suspend PyMOL scene updates while running many cmd ops, then rebuild once. Also disables
    auto_zoom for the duration: PyMOL zooms to any newly created object/selection by default, so
    creating a marker layer (cmd.create) used to yank the camera ("it zooms when I enable PTMs").
    The user's auto_zoom setting is restored afterwards."""
    az = None
    try:
        az = cmd.get("auto_zoom")
        cmd.set("auto_zoom", 0)
        cmd.set("suspend_updates", "on")
    except Exception:
        pass
    try:
        yield
    finally:
        try:
            cmd.set("suspend_updates", "off")
            if az is not None:
                cmd.set("auto_zoom", az)
            # NB: no cmd.refresh() here. refresh() forces a synchronous re-render of the WHOLE
            # scene (every loaded structure) on every action — the main reason actions felt slow
            # with several structures loaded. Turning suspend_updates off already marks the scene
            # dirty; PyMOL repaints on its normal cycle.
        except Exception:
            pass


def _ufv_layer_obj(obj, tag):
    return "ufv_%s_%s" % (re.sub(r"[^A-Za-z0-9]", "_", obj or "x"), tag)


def _ca_scaffold(obj):
    """A persistent, hidden CA-only copy of `obj` (ufv_<obj>_ca), built once. Marker layers are then
    carved out of THIS small object instead of re-scanning the full structure (which is slow for big
    models, e.g. AlphaFill) on every toggle. Returns the scaffold name, or `obj` if it can't be built."""
    name = _ufv_layer_obj(obj, "ca")
    if name not in (cmd.get_object_list() or []):
        try:
            cmd.create(name, "(%s) and name CA and polymer" % obj)
            cmd.hide("everything", name)
            cmd.disable(name)
        except Exception:
            return obj
    return name


def _tune_sphere_obj(lobj):
    # GPU impostor spheres render thousands of points cheaply.
    for k, v in (("sphere_mode", 9), ("sphere_quality", 1), ("sphere_scale", 0.8)):
        try:
            cmd.set(k, v, lobj)
        except Exception:
            pass


def _set_sphere_layer(obj, tag, groups):
    """Render a sphere layer on its OWN lightweight object — a copy of just the annotated Cα atoms,
    with no cartoon. Showing / colouring / re-rendering thousands of points then never rebuilds the
    main structure's cartoon, which is what made large variant sets crawl. groups = {color:
    [positions]} to show, or None to remove the layer."""
    lobj = _ufv_layer_obj(obj, tag)
    try:
        cmd.delete(lobj)
    except Exception:
        pass
    st = _SPHERE_STATE.setdefault(obj, {})
    if not groups:
        st.pop(tag, None)
        return 0
    allpos = sorted({p for positions in groups.values() for p in positions})
    scaffold = _ca_scaffold(obj)                 # carve the layer out of the small CA-only copy
    sel = _sel_for_positions(obj, allpos, ca_only=True, target=scaffold)
    if not sel:
        st.pop(tag, None)
        return 0
    st[tag] = groups
    # Map every UniProt position to its colour index once, then colour the (small) layer object in a
    # SINGLE alter pass. Re-deriving a big resi-list selection per consequence colour and letting
    # PyMOL re-parse it was the remaining cost for large variant sets — alter walks only the copied
    # Cα atoms (a few hundred/thousand) with no selection-string parsing.
    pos_color = {}
    for color, positions in groups.items():
        idx = cmd.get_color_index(_color_name(color))
        for p in positions:
            pos_color[int(p)] = idx
    base_idx = cmd.get_color_index("gray80")

    def _idx(chain, resi):
        return pos_color.get(_resi_to_uni(obj, chain, resi), base_idx)

    with _batch():
        cmd.create(lobj, sel)            # copy only the Cα atoms into a separate, cartoon-free object
        cmd.hide("everything", lobj)
        cmd.show("spheres", lobj)
        _tune_sphere_obj(lobj)
        cmd.alter(lobj, "color = _idx(chain, resi)", space={"_idx": _idx})
        cmd.recolor(lobj)
    return len(allpos)


def _hide_layer(obj, tag):
    _set_sphere_layer(obj, tag, None)


def _delete_layers(obj):
    """Remove every sphere-layer object belonging to `obj`."""
    prefix = _ufv_layer_obj(obj, "")
    for name in list(cmd.get_object_list() or []):
        if name.startswith(prefix):
            try:
                cmd.delete(name)
            except Exception:
                pass
    _SPHERE_STATE.pop(obj, None)


def _reset_cartoon(obj):
    """Drop any cartoon-colour layer: recolour neutral grey and remove domain single-residue spheres."""
    cmd.color("gray80", obj)
    _set_sphere_layer(obj, "dom", None)
    _CARTOON_LAYER.pop(obj, None)


def _set_cartoon_layer(obj, name):
    _reset_cartoon(obj)
    cmd.show("cartoon", obj)
    _CARTOON_LAYER[obj] = name


# ----------------------------------------------------------------------------------------------
# Projection commands  (each is independently toggleable; pair with ufv_hide to remove)
# ----------------------------------------------------------------------------------------------
def ufv_ptms(obj=None, uid=None):
    """ufv_ptms [object [, uniprot_id]] — PTM Ca spheres coloured by category."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    groups = {}
    for p in ann["ptms"]:
        groups.setdefault(p["color"], []).append(p["position"])
        if p["endPosition"] != p["position"]:
            groups[p["color"]].append(p["endPosition"])
    _set_sphere_layer(obj, "ptm", groups)
    print("[UFV] %s: PTM layer on (%d categories)." % (obj, len(groups)))


def _variant_groups(ann, tokens, reviewed_only=False):
    groups = {}
    for v in ann["variants"]:
        if reviewed_only and not v.get("reviewed"):
            continue
        if tokens and _cons_token(v["consequence"]) not in tokens:
            continue
        groups.setdefault(v["consequenceColor"], []).append(v["position"])
    return groups


def ufv_variants(obj=None, uid=None, only="pathogenic"):
    """ufv_variants [object [, uniprot_id [, only]]]
    Variant Ca spheres coloured by clinical consequence. `only` filters by consequence and may
    list several (space/comma separated): pathogenic | benign | uncertain | deleterious | all.
    Defaults to pathogenic to avoid flooding the view; use 'all' for every variant."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    only = (only or "").lower()
    tokens = None if only in ("", "all") else set(re.split(r"[ ,]+", only.strip()))
    _set_sphere_layer(obj, "var", _variant_groups(ann, tokens))
    print("[UFV] %s: variant layer on (%s)." % (obj, only or "all"))


def ufv_sites(obj=None, uid=None):
    """ufv_sites [object [, uniprot_id]] — active / binding / metal sites as amber spheres."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    positions = []
    for s in ann["sites"]:
        positions.append(s["position"])
        if s["endPosition"] != s["position"]:
            positions.append(s["endPosition"])
    _set_sphere_layer(obj, "site", {SITE_COLOR: positions} if positions else None)
    print("[UFV] %s: site layer on (%d sites)." % (obj, len(ann["sites"])))


def _ligand_sel(obj):
    return "(%s) and not polymer and not solvent" % obj


def ufv_ligands(obj=None, uid=None):
    """ufv_ligands [object] — show bound ligands/cofactors (non-polymer, non-water) as sticks,
    coloured by element, with ions as small spheres."""
    obj = _resolve_object(obj)
    sel = _ligand_sel(obj)
    n = cmd.count_atoms(sel)
    if not n:
        print("[UFV] %s: no ligands/cofactors present." % obj)
        return 0
    with _batch():
        cmd.show("sticks", sel)
        ions = "(%s) and inorganic" % obj
        if cmd.count_atoms(ions):
            cmd.show("spheres", ions)
            cmd.set("sphere_scale", 0.3, ions)
        try:
            cmd.color("green", "(%s) and elem C" % sel)        # ligand carbons green (stand out from grey protein)
            cmd.color("atomic", "(%s) and (not elem C)" % sel)  # CPK colour the rest
        except Exception:
            pass
    print("[UFV] %s: ligands shown (%d atoms)." % (obj, n))
    return n


def ufv_ligands_hide(obj=None):
    """ufv_ligands_hide [object] — hide ligand/cofactor sticks and ion spheres."""
    obj = _resolve_object(obj)
    with _batch():
        cmd.hide("sticks", _ligand_sel(obj))
        cmd.hide("spheres", "(%s) and inorganic" % obj)


def enumerate_ligands(obj):
    """Ligand/cofactor components in the structure -> {resn: instance_count} (distinct chain+resi
    copies, no water) so multiple copies of the same ligand are counted."""
    rows = []
    try:
        cmd.iterate(_ligand_sel(obj), "rows.append((resn, chain, resi))", space={"rows": rows})
    except Exception:
        pass
    counts = {}
    for resn in set(rows):
        counts[resn[0]] = counts.get(resn[0], 0) + 1
    return counts


def ligand_chemistry(obj):
    """For each distinct ligand in the structure, fetch CCD chemistry + PubChem fingerprint and rank
    the others by Tanimoto similarity. Returns {resn: {info..., 'similar': [(resn, score)]}}.
    Network-heavy — call off the UI thread."""
    codes = sorted(enumerate_ligands(obj))
    info = {c: ligand_info(c) for c in codes}
    fps = {c: ligand_fingerprint(info[c].get("inchikey")) for c in codes}
    out = {}
    for c in codes:
        sims = []
        for d in codes:
            if d == c:
                continue
            t = tanimoto(fps.get(c), fps.get(d))
            if t is not None:
                sims.append((d, round(t, 3)))
        sims.sort(key=lambda x: -x[1])
        out[c] = dict(info[c], similar=sims)
    return out


def ligand_instances(obj, resn):
    """Distinct copies of a ligand by name -> sorted [(chain, resi), ...] so the GUI can step through
    each occurrence of the same component (e.g. several HEM/ATP in one structure)."""
    rows = []
    try:
        cmd.iterate("(%s) and resn %s and not solvent" % (obj, resn),
                    "rows.append((chain, resi))", space={"rows": rows})
    except Exception:
        pass
    seen, out = set(), []
    for ch, resi in rows:
        key = (ch, resi)
        if key not in seen:
            seen.add(key); out.append(key)
    out.sort(key=lambda x: (x[0], int(re.sub(r"[^0-9-]", "", x[1]) or 0)))
    return out


def ufv_ligand_focus(obj=None, resn=None, chain=None, resi=None):
    """ufv_ligand_focus [object,] resn [, chain, resi] — zoom to a ligand and show it + its 5 Å pocket
    as sticks. With chain+resi, focus that single copy; otherwise all copies of the named component."""
    obj = _resolve_object(obj)
    if chain is not None and resi is not None:
        lig = "(%s) and resn %s and chain %s and resi %s" % (obj, resn, chain, resi)
    else:
        lig = "(%s) and resn %s" % (obj, resn)
    if not cmd.count_atoms(lig):
        return
    pocket = "byres ((%s) around 5) and polymer" % lig
    with _batch():
        _clear_focus(obj)                       # drop any previous residue/ligand focus sticks (+ label)
        cmd.show("sticks", lig)
        cmd.color("green", "%s and elem C" % lig)
        cmd.set("stick_radius", 0.18, lig)      # fatten the focused copy so it reads as 'the selected one'
        cmd.show("sticks", pocket)
        cmd.set("stick_color", "grey70", pocket)
        if chain is not None and resi is not None:   # 3D label so the user can see which copy is selected
            try:
                cmd.delete(_UFV_LIGLABEL)
                cmd.pseudoatom(_UFV_LIGLABEL, selection=lig,
                               label="%s %s/%s" % (resn, chain, resi))
                cmd.set("label_size", 16, _UFV_LIGLABEL)
                cmd.set("label_color", "yellow", _UFV_LIGLABEL)
            except Exception:
                pass
    _FOCUS_SHOWN[obj] = {"sticks": "(%s) or (%s)" % (lig, pocket), "sel": None,
                         "disabled": [], "label": _UFV_LIGLABEL, "lig": lig}
    cmd.zoom(lig, 6, animate=0)


def ufv_domains(obj=None, uid=None):
    """ufv_domains [object [, uniprot_id]] — colour the cartoon by domain / region / repeat."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    _set_cartoon_layer(obj, "domains")
    dom_spheres = {}
    pos_color = {}
    n = 0
    for d in ann["domains"]:
        if d["isRange"]:
            for p in range(d["position"], d["endPosition"] + 1):
                pos_color[p] = d["color"]
        else:
            dom_spheres.setdefault(d["color"], []).append(d["position"])
        n += 1
    if pos_color:
        _recolor_positions(obj, pos_color)
    _set_sphere_layer(obj, "dom", dom_spheres or None)
    print("[UFV] %s: coloured %d domain/region features." % (obj, n))


def ufv_topology(obj=None, uid=None):
    """ufv_topology [object [, uniprot_id]] — colour the cartoon by membrane topology."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    _set_cartoon_layer(obj, "topology")
    pos_color = {}
    for t in ann["topology"]:
        for p in range(t["start"], t["end"] + 1):
            pos_color[p] = t["color"]
    _recolor_positions(obj, pos_color)
    print("[UFV] %s: coloured %d topology segments." % (obj, len(ann["topology"])))


def ufv_alphamissense(obj=None, uid=None):
    """ufv_alphamissense [object [, uniprot_id]] — colour by mean AlphaMissense pathogenicity."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    am = ann["amMean"]
    if not am:
        print("[UFV] no AlphaMissense scores available for %s." % ann["uid"])
        return
    _set_cartoon_layer(obj, "alphamissense")
    pos_color = {}
    for pos, avg in am.items():
        pos_color[pos] = ("#b71c1c" if avg >= 0.78 else "#e06666" if avg >= 0.564
                          else "#b9c2cf" if avg >= 0.34 else "#3d85c8")
    _recolor_positions(obj, pos_color, base_hex="#b9c2cf")
    print("[UFV] %s: coloured by AlphaMissense (%d scored positions)." % (obj, len(am)))


def ufv_burden(obj=None, uid=None):
    """ufv_burden [object [, uniprot_id]] — colour mutation/phenotype burden-positive residues."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    _set_cartoon_layer(obj, "burden")
    _recolor_positions(obj, {p: "#e65100" for p in ann["burden"]})
    print("[UFV] %s: %d burden-positive residues." % (obj, len(ann["burden"])))


def ufv_plddt(obj=None, uid=None):
    """ufv_plddt [object] — colour by AlphaFold pLDDT (the B-factor column), exact AFDB bins.
    Only meaningful for predicted models (the B column is pLDDT there)."""
    obj = _resolve_object(obj)
    if not _structure_is_predicted(obj):
        print("[UFV] pLDDT applies to predicted models only (%s is experimental)." % obj)
        return
    _set_cartoon_layer(obj, "plddt")
    with _batch():
        # AlphaFold DB bins: >=90 dark blue, 70-90 light blue, 50-70 yellow, <50 orange.
        cmd.color(_color_name("#ff7d45"), "(%s) and b<50" % obj)
        cmd.color(_color_name("#ffdb13"), "(%s) and (b>50 or b=50) and b<70" % obj)
        cmd.color(_color_name("#65cbf3"), "(%s) and (b>70 or b=70) and b<90" % obj)
        cmd.color(_color_name("#0053d6"), "(%s) and (b>90 or b=90)" % obj)
    print("[UFV] %s: coloured by pLDDT." % obj)


def ufv_bfactor(obj=None, uid=None):
    """ufv_bfactor [object] — colour by crystallographic/EM B-factor with the extension's exact
    blue-white-red gradient (#313695 -> #f7f7f7 -> #d73027), clamped 0-100. Experimental only
    (for predicted models the B column is pLDDT — use that mode instead)."""
    obj = _resolve_object(obj)
    if _structure_is_predicted(obj):
        print("[UFV] B-factor applies to experimental structures only (%s is a predicted model)." % obj)
        return
    _set_cartoon_layer(obj, "bfactor")
    try:
        cmd.spectrum("b", "0x313695 0xf7f7f7 0xd73027", obj, minimum=0, maximum=100)
    except Exception:
        cmd.spectrum("b", "blue_white_red", obj)
    print("[UFV] %s: coloured by B-factor." % obj)


# ----------------------------------------------------------------------------------------------
# Geometry + structure-dependent analyses (hotspots, contact hubs) + residue report / focus
# ----------------------------------------------------------------------------------------------
_GEOM_CACHE = {}    # obj -> geom list
_FOCUS_SHOWN = {}   # obj -> stick selection currently shown by a focus


def _ca_geometry(obj, refresh=False):
    """Modelled Cα atoms of `obj` as [{uniPos, chain, resi, x, y, z}] (cached per object)."""
    if not refresh and obj in _GEOM_CACHE:
        return _GEOM_CACHE[obj]
    rows = []
    try:
        cmd.iterate_state(1, "(%s) and name CA and polymer" % obj,
                          "rows.append((chain, resi, x, y, z))", space={"rows": rows})
    except Exception:
        rows = []
    geom = [{"uniPos": _resi_to_uni(obj, ch, resi), "chain": ch, "resi": resi, "x": x, "y": y, "z": z}
            for (ch, resi, x, y, z) in rows]
    _GEOM_CACHE[obj] = geom
    return geom


def _recolor_positions(obj, pos_color, base_hex=None):
    """Colour `obj`'s polymer atoms from a {uniProtPos: hex} map in ONE `alter` pass. This replaces
    the old 'one cmd.color per colour group with a big resi-list selection' pattern, whose repeated
    selection parsing over the whole structure was the slow part of every cartoon-colouring mode.
    `alter` walks the atoms once in C; `_resi_to_uni` is memoised per residue. Atoms not in the map
    keep their colour, or take `base_hex` if given. Returns how many positions were coloured."""
    idx_map = {}
    for pos, hexc in pos_color.items():
        idx_map[int(pos)] = cmd.get_color_index(_color_name(hexc))
    if not idx_map and base_hex is None:
        return 0
    base_idx = cmd.get_color_index(_color_name(base_hex)) if base_hex else None
    memo = {}

    def _idx(chain, resi, color):
        key = (chain, resi)
        if key not in memo:
            w = idx_map.get(_resi_to_uni(obj, chain, resi))
            if w is None:
                w = base_idx
            memo[key] = w
        v = memo[key]
        return color if v is None else v

    with _batch():
        cmd.alter("(%s) and polymer" % obj, "color = _idx(chain, resi, color)", space={"_idx": _idx})
        cmd.recolor(obj)
    return len(idx_map)


def _color_tiers(obj, tiers, colors):
    pos_color = {}
    for uni, tier in tiers.items():
        c = colors.get(tier)
        if c:
            pos_color[uni] = c
    _recolor_positions(obj, pos_color)
    return len(tiers)


def ufv_hotspots(obj=None, uid=None):
    """ufv_hotspots [object [, uniprot_id]] — colour pathogenic-enrichment hotspot tiers."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    _set_cartoon_layer(obj, "hotspots")
    tiers = _per_chain_merge(_ca_geometry(obj), lambda gs: _compute_hotspots(gs, ann["variants"]))
    n = _color_tiers(obj, tiers, {"strong": "#b71c1c", "moderate": "#e64a19", "weak": "#ffa726"})
    print("[UFV] %s: %d hotspot residues." % (obj, n))


def ufv_contacthubs(obj=None, uid=None):
    """ufv_contacthubs [object] — colour long-range contact-hub tiers (Cα betweenness)."""
    obj = _resolve_object(obj)
    _set_cartoon_layer(obj, "contacthubs")
    tiers = _per_chain_merge(_ca_geometry(obj), _betweenness_hubs)
    n = _color_tiers(obj, tiers, {"strong": "#6a1b9a", "moderate": "#ab47bc"})
    print("[UFV] %s: %d contact-hub residues." % (obj, n))


def ufv_pockets(obj=None, uid=None):
    """ufv_pockets [object [, uniprot_id]] — constraint-pocket prioritisation: residues more
    evolutionarily constrained (AlphaMissense) than their structural burial explains, clustered in
    3D. Buried/concave candidates are teal ('pocket'), exposed ones purple."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    if not ann.get("amMean"):
        print("[UFV] no AlphaMissense data — constraint-pocket needs it.")
        return
    _set_cartoon_layer(obj, "pockets")
    cats = _compute_pockets(_ca_geometry(obj), ann["amMean"])
    n = _color_tiers(obj, cats, {"pocket": "#00897b", "exposed": "#8e24aa"})
    print("[UFV] %s: %d constraint-pocket residues." % (obj, n))


def residue_report(obj, uid, uni_pos):
    """Everything the extension's detail panel shows for one residue, as a plain dict."""
    ann = fetch_annotations(uid)
    geom = _ca_geometry(obj)
    seq = ann["sequence"]
    uni_pos = int(uni_pos)
    aa = seq[uni_pos - 1] if 0 < uni_pos <= len(seq) else "?"
    ptm_pos = {p["position"] for p in ann["ptms"]}
    var_pos = {v["position"]: [] for v in ann["variants"]}
    for v in ann["variants"]:
        var_pos[v["position"]].append(v)
    near = []
    for nuni, dist in _residue_neighborhood(geom, uni_pos, 12.0):
        tags = []
        if nuni in ptm_pos:
            tags.append("PTM")
        if nuni in var_pos:
            tags.append(_cons_token(var_pos[nuni][0]["consequence"]))
        near.append({"pos": nuni, "dist": round(dist, 1), "tags": tags})
    # ligands/cofactors within 5 Å of this residue (read from the structure)
    ligs = []
    rsel = _sel_for_positions(obj, [uni_pos])
    if rsel:
        try:
            cmd.iterate("(%s) and not polymer and not solvent and ((%s) around 5)" % (obj, rsel),
                        "ligs.append(resn)", space={"ligs": ligs})
        except Exception:
            pass
    return {
        "uid": ann["uid"], "pos": uni_pos, "aa": aa,
        "ptms": [p for p in ann["ptms"] if p["position"] <= uni_pos <= p["endPosition"]],
        "variants": var_pos.get(uni_pos, []),
        "sites": [s for s in ann["sites"] if s["position"] <= uni_pos <= s["endPosition"]],
        "domains": [d for d in ann["domains"] if d["position"] <= uni_pos <= d["endPosition"]],
        "amMean": ann["amMean"].get(uni_pos),
        "amProfile": ann.get("amByPos", {}).get(uni_pos, []),  # full 19-substitution profile
        "nearby": near,
        "nearbyLigands": sorted(set(ligs)),
    }


def format_report(rep):
    """Render a residue_report dict as a text block for the GUI / console."""
    if not rep:
        return ""
    L = ["%s  position %d  (%s)" % (rep["uid"], rep["pos"], rep["aa"])]
    for p in rep["ptms"]:
        L.append("  PTM: %s" % p["description"])
    for s in rep["sites"]:
        L.append("  Site: %s" % s["description"])
    for d in rep["domains"]:
        L.append("  Domain: %s" % d["description"])
    for v in rep["variants"]:
        extra = []
        if v.get("clinVar"):
            extra.append("ClinVar: %s" % v["clinVar"])
        if v.get("alphaMissense") is not None:
            extra.append("AM %.2f" % v["alphaMissense"])
        if v.get("diseases"):
            extra.append("; ".join(v["diseases"]))
        L.append("  Variant %s%d%s — %s%s" % (v.get("wildType", ""), rep["pos"], v.get("mutant", ""),
                                              v["consequence"], (" [" + " | ".join(extra) + "]") if extra else ""))
    if rep["amMean"] is not None:
        L.append("  AlphaMissense (mean): %.3f" % rep["amMean"])
    if rep["amProfile"]:
        L.append("  AM top: " + ", ".join("%s%s %.2f" % (wt, mt, sc) for wt, mt, sc in rep["amProfile"]))
    if rep.get("nearbyLigands"):
        L.append("  Nearby ligands (≤5Å): " + ", ".join(rep["nearbyLigands"]))
    if rep["nearby"]:
        near = ", ".join("%d (%.1fÅ%s)" % (n["pos"], n["dist"], " " + "/".join(n["tags"]) if n["tags"] else "")
                         for n in rep["nearby"][:12])
        L.append("  Nearby ≤12Å: " + near)
    return "\n".join(L)


def _am_color(score):
    """AlphaMissense colour: red likely-pathogenic (>=0.564), green likely-benign (<=0.34), amber
    ambiguous between — matching the extension's substitution grid."""
    if score is None:
        return "#888"
    return "#d32f2f" if score >= 0.564 else "#4caf50" if score <= 0.34 else "#f5a623"


def format_report_html(rep, expanded=False):
    """Clean label/value layout for the GUI detail panel, mirroring the extension's residue card:
    a header, then colour-coded section headers with aligned 'label : value' rows (no indentation).
    `expanded` controls whether each variant's evidence (ClinVar / review / dbSNP / disease) is shown
    or folded behind a toggle link. Nearby residues are clickable (ufv:res: anchors)."""
    if not rep:
        return ""

    def esc(s):
        return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    T = ['<table cellspacing="0" cellpadding="3" '
         'style="font-family:sans-serif; font-size:12px; border-collapse:collapse;">']
    T.append('<tr><td colspan="2" style="font-size:15px; padding-bottom:2px;">'
             '<b>%s%d</b></td></tr>' % (esc(rep["aa"]), rep["pos"]))

    def hdr(title, color):
        T.append('<tr><td colspan="2" style="border-top:1px solid #ccc; padding-top:5px;">'
                 '<b style="color:%s;">%s</b></td></tr>' % (color, esc(title)))

    def row(label, value, vcolor=None, bold=False):
        # label is a controlled string / pre-built HTML (variant mutations); value is caller-escaped
        v = ('<span style="color:%s;">%s</span>' % (vcolor, value)) if vcolor else value
        lab = ('<b>%s</b>' % label) if bold else ('<span style="color:#888;">%s</span>' % label)
        T.append('<tr><td style="vertical-align:top; white-space:nowrap;">%s</td>'
                 '<td style="vertical-align:top;">%s</td></tr>' % (lab, v))

    if rep["variants"]:
        hdr("Variants", "#c62828")
        has_evidence = any(v.get("clinVar") or v.get("clinVarReview") or v.get("rsIds") or v.get("diseases")
                           for v in rep["variants"])
        for v in rep["variants"]:
            mut = "%s%d%s" % (esc(v.get("wildType", "")), rep["pos"], esc(v.get("mutant", "")))
            c = v.get("consequenceColor", "#9e9e9e")
            row('<span style="color:%s;">%s</span>' % (c, mut), esc(v["consequence"]), c, bold=True)
            if v.get("alphaMissense") is not None:
                row("AlphaMissense", "%.2f" % v["alphaMissense"], _am_color(v["alphaMissense"]))
            if expanded:  # evidence rows, folded by default
                if v.get("clinVar"):
                    row("ClinVar", esc(v["clinVar"]))
                if v.get("clinVarReview"):
                    row("Review", esc(v["clinVarReview"]))
                if v.get("rsIds"):
                    row("dbSNP", esc(", ".join(v["rsIds"])))
                if v.get("diseases"):
                    row("Disease", esc("; ".join(v["diseases"])))
        if has_evidence:
            link = "▾ hide evidence" if expanded else "▸ show evidence"
            T.append('<tr><td colspan="2" style="padding-top:2px;">'
                     '<a href="ufv:evi" style="color:#1565c0; text-decoration:none; font-size:11px;">%s</a>'
                     '</td></tr>' % link)

    if rep["ptms"] or rep["sites"] or rep["domains"]:
        hdr("Features", "#1565c0")
        for p in rep["ptms"]:
            row("PTM", esc(p["description"]), p.get("color"))
        for s in rep["sites"]:
            row("Site", esc(s["description"]), SITE_COLOR)
        for d in rep["domains"]:
            row("Domain", esc(d["description"]), d.get("color"))

    if rep["amMean"] is not None:
        hdr("AlphaMissense", "#6a1b9a")
        row("Mean", "%.3f" % rep["amMean"], _am_color(rep["amMean"]))
        if rep.get("amProfile"):
            # full substitution grid, 10 per row, coloured by pathogenicity
            cells = ['<tr><td colspan="2"><table cellspacing="2"><tr>']
            for i, (w, m, sc) in enumerate(rep["amProfile"]):
                if i and i % 10 == 0:
                    cells.append('</tr><tr>')
                cells.append('<td style="background:#f3f3f3; padding:2px 4px; text-align:center;">'
                             '<b style="color:%s;">%s</b><br><span style="color:%s; font-size:10px;">%.2f</span></td>'
                             % (_am_color(sc), esc(m), _am_color(sc), sc))
            cells.append('</tr></table></td></tr>')
            T.append("".join(cells))

    if rep.get("nearbyLigands"):
        hdr("Nearby ligands (≤5 Å)", "#6d4c41")
        row("Ligands", esc(", ".join(rep["nearbyLigands"])))

    if rep["nearby"]:
        tagcol = {"pathogenic": "#ef5350", "deleterious": "#ffa726", "benign": "#66bb6a"}

        def chip(n):
            col = "#37474f"
            for t in n["tags"]:
                col = "#b87800" if t == "PTM" else tagcol.get(t, col)
            tip = ("  " + "/".join(n["tags"])) if n["tags"] else ""
            # clickable: jumps the report to that residue (ufv:res anchor)
            return ('<a href="ufv:res:%d" title="residue %d%s" style="color:%s; text-decoration:none;">'
                    '%d&nbsp;(%.1f)</a>' % (n["pos"], n["pos"], esc(tip), col, n["pos"], n["dist"]))

        inner = [n for n in rep["nearby"] if n["dist"] <= 8.0]
        outer = [n for n in rep["nearby"] if 8.0 < n["dist"] <= 12.0]
        hdr("Nearby residues", "#37474f")
        if inner:
            row("≤ 8 Å", ", ".join(chip(n) for n in inner[:24]))
        if outer:
            row("8–12 Å", ", ".join(chip(n) for n in outer[:24]))

    T.append("</table>")
    return "".join(T)


def format_ligand_html(resn, info, instances=1, copies=None, cur=0):
    """Clean label/value table for a ligand's chemistry + Tanimoto similars (extension-like).
    When the component occurs multiple times, `copies` = [(chain, resi), ...] and `cur` is the focused
    copy index; a ◂ i/N ▸ stepper (ufv:lig anchors) lets the user circulate the occurrences."""
    def esc(s):
        return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    n = instances if not copies else len(copies)
    T = ['<table cellspacing="0" cellpadding="3" style="font-family:sans-serif; font-size:12px;">']
    title = "%s%s" % (esc(resn), (" ×%d" % n) if n > 1 else "")
    T.append('<tr><td colspan="2" style="font-size:14px;"><b style="color:#2e7d32;">%s</b> '
             '<span style="color:#888;">(CCD)</span></td></tr>' % title)
    if copies and len(copies) > 1:
        ch, resi = copies[cur % len(copies)]
        T.append('<tr><td colspan="2" style="padding:3px 0;">'
                 '<span style="background:#e8f0fe; padding:2px 6px; border-radius:3px;">'
                 '<a href="ufv:lig:prev" style="color:#1565c0; text-decoration:none; font-size:14px;">◂</a>'
                 '&nbsp;<b>showing copy %d of %d</b> &mdash; <b style="color:#2e7d32;">chain %s · %s</b>&nbsp;'
                 '<a href="ufv:lig:next" style="color:#1565c0; text-decoration:none; font-size:14px;">▸</a>'
                 '</span> <span style="color:#888; font-size:10px;">(highlighted in 3D)</span>'
                 '</td></tr>' % (cur % len(copies) + 1, len(copies), esc(ch or "?"), esc(resi)))

    def row(label, value):
        T.append('<tr><td style="color:#888; vertical-align:top; white-space:nowrap;">%s</td>'
                 '<td style="vertical-align:top;">%s</td></tr>' % (esc(label), value))
    if info.get("name"):
        row("Name", esc(info["name"]))
    if info.get("formula"):
        row("Formula", esc(info["formula"]))
    if info.get("smiles"):
        row("SMILES", '<span style="font-family:monospace; font-size:11px;">%s</span>' % esc(info["smiles"]))
    if info.get("inchikey"):
        row("InChIKey", '<span style="font-family:monospace; font-size:11px;">%s</span>' % esc(info["inchikey"]))
    if info.get("drugbank"):
        db = esc(info["drugbank"])
        row("DrugBank", '<a href="https://go.drugbank.com/drugs/%s" '
                        'style="color:#1565c0; text-decoration:none;">%s ↗</a>' % (db, db))
    sims = info.get("similar") or []
    if sims:
        T.append('<tr><td colspan="2" style="border-top:1px solid #ccc; padding-top:5px;">'
                 '<b style="color:#6d4c41;">Similar in structure (Tanimoto)</b></td></tr>')
        cells = ['<tr><td colspan="2"><table cellspacing="2"><tr>']
        for i, (d, s) in enumerate(sims[:12]):
            if i and i % 3 == 0:
                cells.append('</tr><tr>')
            cells.append('<td style="background:#f3f3f3; padding:2px 6px;">'
                         '<b style="color:#2e7d32;">%s</b>&nbsp;%.2f</td>' % (esc(d), s))
        cells.append('</tr></table></td></tr>')
        T.append("".join(cells))
    T.append("</table>")
    return "".join(T)


def _annotation_color_map(ann):
    """uniprot pos -> annotation colour (variant > PTM > site), matching the extension's focus."""
    m = {}
    for s in ann.get("sites", []):
        m[s["position"]] = s["color"]
    for p in ann.get("ptms", []):
        m[p["position"]] = p["color"]
    for v in ann.get("variants", []):
        m[v["position"]] = v["consequenceColor"]
    return m


_DIMMED = set()  # objects whose cartoon is currently dimmed for focus
_UFV_LIGLABEL = "ufv_liglabel"  # transient pseudoatom labelling the focused ligand copy


def _clear_focus(obj, restore_spheres=True):
    """Remove the current focus sticks and re-enable the sphere layers it hid (keeps the dim)."""
    try:
        cmd.delete(_UFV_LIGLABEL)
    except Exception:
        pass
    prev = _FOCUS_SHOWN.pop(obj, None)
    if not prev:
        return
    try:
        cmd.hide("sticks", prev["sticks"])
        cmd.unset("stick_color", prev["sticks"])
        if prev.get("lig"):
            cmd.unset("stick_radius", prev["lig"])  # un-fatten the previously focused ligand copy
        if prev.get("sel"):
            cmd.hide("spheres", prev["sel"])
            cmd.color("gray80", prev["sel"])  # undo the element colouring of the selected residue
    except Exception:
        pass
    if restore_spheres:
        for lobj in prev.get("disabled", []):
            try:
                cmd.enable(lobj)
            except Exception:
                pass


def ufv_focus(obj=None, uni_pos=None):
    """ufv_focus [object,] uniprot_pos — zoom into a residue: show it and its 5 Å neighbourhood as
    sticks coloured by annotation, dim the cartoon, and hide the sphere markers there. The cartoon
    is dimmed only ONCE per focus session (changing cartoon_transparency rebuilds the whole cartoon
    mesh, which is what made each zoom slow) — subsequent residue clicks just move the sticks."""
    obj = _resolve_object(obj)
    ann = _CACHE.get(_resolve_uid(None)) or {}
    sel = _sel_for_positions(obj, [int(uni_pos)])
    if not sel:
        print("[UFV] position %s not modelled in %s." % (uni_pos, obj))
        return
    nb = "byres ((%s) around 5) and polymer" % sel
    focus_sel = "(%s) or (%s)" % (sel, nb)
    cmap = _annotation_color_map(ann)
    with _batch():
        _clear_focus(obj)                       # clear previous sticks (keep dim — no rebuild)
        if obj not in _DIMMED:                   # dim once: the only cartoon rebuild
            cmd.set("cartoon_transparency", 0.55, obj)
            _DIMMED.add(obj)
        cmd.show("sticks", focus_sel)
        cmd.set("stick_color", "grey70", focus_sel)
        nb_atoms = []
        cmd.iterate(nb + " and name CA", "nb_atoms.append((chain, resi))", space={"nb_atoms": nb_atoms})
        by_color, per_chain = {}, {}
        for ch, resi in set(nb_atoms):
            per_chain.setdefault(ch, []).append(resi)
            c = cmap.get(_resi_to_uni(obj, ch, resi))
            if c:
                by_color.setdefault(c, {}).setdefault(ch, []).append(resi)
        for c, chres in by_color.items():
            clauses = " or ".join("(chain %s and resi %s)" % (ch, "+".join(rs)) for ch, rs in chres.items())
            cmd.set("stick_color", _color_name(c), "(%s) and (%s)" % (obj, clauses))
        # The selected residue: ball-and-stick with atom-type colours. Carbons take the residue's
        # DISEASE annotation colour (variant/PTM/site) or neutral grey — never the active cartoon
        # colouring (hub/burden/pocket), which we explicitly override here so it can't bleed onto the
        # sidechain. Heteroatoms keep element colours.
        dcol = cmap.get(int(uni_pos))
        ccol = _color_name(dcol) if dcol else "gray85"
        cmd.show("spheres", sel)
        cmd.set("sphere_scale", 0.28, sel)
        cmd.color("atomic", "(%s) and (not elem C)" % sel)
        cmd.color(ccol, "(%s) and elem C" % sel)
        cmd.unset("stick_color", sel)  # let element/disease colours show on the selected residue's sticks
        # Hide ALL annotation sphere layers while zoomed in (declutters; cmd.disable is instant —
        # no rebuild). Restored on reset / next clear.
        disabled = []
        objlist = cmd.get_object_list() or []
        for tag in ("ptm", "var", "site", "dom", "filt"):
            lobj = _ufv_layer_obj(obj, tag)
            if lobj in objlist:
                try:
                    cmd.disable(lobj)
                    disabled.append(lobj)
                except Exception:
                    pass
    _FOCUS_SHOWN[obj] = {"sticks": focus_sel, "sel": sel, "disabled": disabled}
    cmd.zoom(focus_sel, 3, animate=0)


def ufv_resetview(obj=None):
    """ufv_resetview [object] — clear focus sticks, restore spheres + cartoon, and zoom back out."""
    obj = _resolve_object(obj)
    with _batch():
        _clear_focus(obj)
        if obj in _DIMMED:
            try:
                cmd.set("cartoon_transparency", 0.0, obj)
            except Exception:
                pass
            _DIMMED.discard(obj)
    cmd.zoom(obj)


def ufv_report(uni_pos, obj=None, uid=None):
    """ufv_report uniprot_pos [, object [, uid]] — print the residue report + focus the residue."""
    obj = _resolve_object(obj)
    uid = _resolve_uid(uid)
    rep = residue_report(obj, uid, uni_pos)
    print(format_report(rep))
    ufv_focus(obj, uni_pos)
    return rep


def _loaded_structures():
    return [o for o in (cmd.get_object_list() or []) if not o.startswith("ufv_")]


def ufv_align(reference=None):
    """ufv_align [reference] — superpose all loaded protein structures onto a reference (default the
    first loaded), so several structures of the protein overlay. Uses cmd.align."""
    objs = _loaded_structures()
    if len(objs) < 2:
        print("[UFV] need at least 2 loaded structures to align.")
        return
    ref = reference or objs[0]
    for o in objs:
        if o == ref:
            continue
        try:
            r = cmd.align(o, ref)
            print("[UFV] aligned %s onto %s: RMSD %.2f Å over %d atoms." % (o, ref, r[0], r[1]))
        except Exception as e:
            print("[UFV] align %s failed: %s" % (o, e))


def ufv_structures(uid=None):
    """ufv_structures [uniprot_id] — list the structures available for an accession."""
    uid = _resolve_uid(uid)
    fetch_annotations(uid)
    structs = list_structures(uid)
    print("[UFV] %d structures for %s:" % (len(structs), uid))
    for s in structs:
        print("    %-22s %s" % (s["key"], s["label"]))
    _STATE["structures"] = structs
    return structs


def ufv_use(key=None, uid=None):
    """ufv_use key [, uniprot_id] — download and load a structure listed by ufv_structures,
    setting its residue numbering automatically (AlphaFold/computed = identity, PDB = SIFTS)."""
    uid = _resolve_uid(uid)
    structs = _STATE.get("structures") or ufv_structures(uid)
    match = next((s for s in structs if s["key"].lower() == str(key).lower()), None)
    if not match:
        print("[UFV] no structure '%s' — run ufv_structures to list them." % key)
        return
    return _load_structure(match, uid)


def _load_structure(struct, uid, name=None):
    name = name or re.sub(r"[^A-Za-z0-9_]", "_", struct["key"])
    print("[UFV] downloading %s ..." % struct["url"])
    data = _get_text(struct["url"])
    if not data:
        print("[UFV] download failed.")
        return None
    ext = "cif" if struct.get("fmt") == "mmcif" else "pdb"
    tmp = os.path.join(tempfile.gettempdir(), "ufv_%s.%s" % (name, ext))
    with open(tmp, "w") as fh:
        fh.write(data)
    cmd.load(tmp, name)
    try:
        os.remove(tmp)
    except OSError:
        pass
    cmd.hide("everything", name)
    cmd.show("cartoon", name)
    cmd.color("gray80", name)
    if struct.get("numbering") == "sifts" and struct.get("pdbId"):
        _set_sifts_map(name, struct["pdbId"], uid)
    else:
        _set_identity_map(name)
    _GEOM_CACHE.pop(name, None)
    _STATE["obj"] = name
    _STATE.setdefault("sources", {})[name] = struct.get("source", "PDB")
    cmd.orient(name)
    print("[UFV] loaded %s (%s). Numbering: %s." % (name, struct["label"], struct.get("numbering")))
    return name


def ufv_hide(obj=None, layer=None):
    """ufv_hide [object [, layer]] — hide one layer: ptms | variants | sites | domains |
    topology | alphamissense | cartoon. With no layer, hides every UFV overlay (like ufv_clear)."""
    obj = _resolve_object(obj)
    layer = (layer or "").lower()
    if not layer:
        ufv_clear(obj)
        return
    if layer in CARTOON_TAGS or layer == "cartoon":
        _reset_cartoon(obj)
    else:
        _hide_layer(obj, POINT_TAGS.get(layer, layer))
    print("[UFV] %s: hid %s layer." % (obj, layer))


def ufv_clear(obj=None):
    """ufv_clear [object] — remove all UFV overlays and reset colour."""
    obj = _resolve_object(obj)
    _delete_layers(obj)
    _CARTOON_LAYER.pop(obj, None)
    _FOCUS_SHOWN.pop(obj, None)
    if obj in _DIMMED:
        try:
            cmd.set("cartoon_transparency", 0.0, obj)
        except Exception:
            pass
        _DIMMED.discard(obj)
    if obj:
        cmd.color("gray80", obj)
    print("[UFV] cleared UFV overlays on %s." % obj)


def _structure_is_predicted(obj):
    """True for AlphaFold/computed models (B-factor column is pLDDT), False for experimental."""
    return (_STATE.get("sources") or {}).get(obj, "PDB") != "PDB"


def actual_coverage(obj, seqlen):
    """Fraction of the UniProt sequence actually modelled in `obj` (CA atoms with a mapped UniProt
    position) — the real coverage, like the extension's modelled-residue count (not the SIFTS range)."""
    if not seqlen:
        return None
    covered = {g["uniPos"] for g in _ca_geometry(obj) if g["uniPos"] is not None}
    return round(len(covered) / seqlen * 100, 1)


# ----------------------------------------------------------------------------------------------
# Mapping / loading commands
# ----------------------------------------------------------------------------------------------
def ufv_fetch(uid):
    """ufv_fetch uniprot_id — download and cache all annotation layers for an accession."""
    fetch_annotations(uid, force=True)


def ufv_map(obj, uid, mode="identity", pdb_id=None):
    """ufv_map object, uniprot_id [, mode [, pdb_id]]
    Establish how `object` residue numbers relate to UniProt positions.
      mode = identity        resi == UniProt position (AlphaFold / UniProt-numbered models)
      mode = sifts, <pdb>    map through PDBe/SIFTS for the given PDB id"""
    obj = _resolve_object(obj)
    fetch_annotations(uid)
    if mode.lower() == "sifts":
        pid = pdb_id or uid  # allow `ufv_map obj, P35498, sifts, 7dtd`
        _set_sifts_map(obj, pid, uid.upper())
    else:
        _set_identity_map(obj)
        print("[UFV] %s mapped with identity numbering (resi == UniProt position)." % obj)


def ufv_chain(obj, chain, uniprot_start, resi_start=None, resi_end=None):
    """ufv_chain object, chain, uniprot_start [, resi_start [, resi_end]]
    Manually anchor one chain's numbering: residue `resi_start` (default = uniprot_start)
    corresponds to UniProt position `uniprot_start`, linearly. `resi_end` bounds the chain.
    Ideal for trajectories with custom numbering or to define each chain's start/finish."""
    obj = _resolve_object(obj)
    mp = _ensure_manual_map(obj)
    u_start = int(uniprot_start)
    r_start = int(resi_start) if resi_start not in (None, "", "None") else u_start
    r_end = int(resi_end) if resi_end not in (None, "", "None") else None
    mp["chains"][chain] = {"u_start": u_start, "resi_start": r_start, "resi_end": r_end}
    print("[UFV] %s chain %s: resi %d == UniProt %d%s"
          % (obj, chain, r_start, u_start, (" (..%d)" % r_end) if r_end else ""))


def ufv_load(uid, name=None):
    """ufv_load uniprot_id [, name] — download the AlphaFold model, load it as a plain grey
    cartoon (identity numbering), and fetch the annotations. Nothing is projected automatically:
    choose layers from ufv_gui or the ufv_* commands so the view stays legible."""
    uid = uid.strip().upper()
    name = name or uid
    url = _alphafold_url(uid)
    if not url:
        print("[UFV] no AlphaFold model found for %s." % uid)
        return
    print("[UFV] downloading %s ..." % url)
    data = _get_text(url)
    if not data:
        print("[UFV] failed to download model.")
        return
    tmp = os.path.join(tempfile.gettempdir(), "ufv_%s.pdb" % uid)
    with open(tmp, "w") as fh:
        fh.write(data)
    cmd.load(tmp, name)
    try:
        os.remove(tmp)
    except OSError:
        pass
    cmd.hide("everything", name)
    cmd.show("cartoon", name)
    cmd.color("gray80", name)
    _set_identity_map(name)
    _GEOM_CACHE.pop(name, None)
    _STATE["obj"] = name
    _STATE.setdefault("sources", {})[name] = "AlphaFold"
    fetch_annotations(uid)
    cmd.orient(name)
    print("[UFV] loaded %s. Open the panel (ufv_gui) or add layers: ufv_ptms / "
          "ufv_variants / ufv_sites / ufv_domains / ufv_topology / ufv_alphamissense %s" % (name, name))


def ufv_info(uid=None):
    """ufv_info [uniprot_id] — print a summary of the cached annotations."""
    ann = fetch_annotations(_resolve_uid(uid))
    print("[UFV] %s  seq=%daa  PTMs=%d  variants=%d  sites=%d  domains=%d  topology=%d"
          % (ann["uid"], len(ann["sequence"]), len(ann["ptms"]), len(ann["variants"]),
             len(ann["sites"]), len(ann["domains"]), len(ann["topology"])))
    npath = sum(1 for v in ann["variants"] if "pathogenic" in v["consequence"].lower())
    print("[UFV] pathogenic/likely-pathogenic variants: %d" % npath)


# ----------------------------------------------------------------------------------------------
# Register commands
# ----------------------------------------------------------------------------------------------
if _HAS_PYMOL:
    for _name, _fn in [
        ("ufv_load", ufv_load), ("ufv_fetch", ufv_fetch), ("ufv_map", ufv_map),
        ("ufv_chain", ufv_chain), ("ufv_ptms", ufv_ptms), ("ufv_variants", ufv_variants),
        ("ufv_sites", ufv_sites), ("ufv_domains", ufv_domains), ("ufv_topology", ufv_topology),
        ("ufv_alphamissense", ufv_alphamissense), ("ufv_burden", ufv_burden),
        ("ufv_plddt", ufv_plddt), ("ufv_bfactor", ufv_bfactor),
        ("ufv_hotspots", ufv_hotspots), ("ufv_contacthubs", ufv_contacthubs), ("ufv_pockets", ufv_pockets),
        ("ufv_structures", ufv_structures), ("ufv_use", ufv_use),
        ("ufv_report", ufv_report), ("ufv_focus", ufv_focus), ("ufv_resetview", ufv_resetview),
        ("ufv_align", ufv_align),
        ("ufv_ligands", ufv_ligands), ("ufv_ligands_hide", ufv_ligands_hide),
        ("ufv_ligand_focus", ufv_ligand_focus),
        ("ufv_hide", ufv_hide), ("ufv_clear", ufv_clear), ("ufv_info", ufv_info),
    ]:
        cmd.extend(_name, _fn)


# ----------------------------------------------------------------------------------------------
# Qt control panel (PyMOL is Qt-based) + Plugin Manager hook
# ----------------------------------------------------------------------------------------------
_gui_ref = None  # keep a reference so the window isn't garbage-collected


def ufv_gui():
    """ufv_gui - open the Qt control panel (visualization driver + report generator)."""
    global _gui_ref
    try:
        from pymol.Qt import QtWidgets, QtCore, QtGui
    except Exception as e:
        print("[UFV] Qt GUI unavailable (%s). Use the ufv_* commands instead." % e)
        return
    cls = _build_panel_class(QtWidgets, QtCore, QtGui)
    _gui_ref = cls()
    _gui_ref.show()
    return _gui_ref


def _build_panel_class(QtWidgets, QtCore, QtGui):
    from collections import Counter

    class _Worker(QtCore.QThread):
        done = QtCore.Signal(object)

        def __init__(self, fn):
            super(_Worker, self).__init__()
            self._fn = fn

        def run(self):
            try:
                self.done.emit(self._fn())
            except Exception as exc:
                self.done.emit(exc)

    class _UFVPanel(QtWidgets.QWidget):
        def __init__(self):
            super(_UFVPanel, self).__init__()
            self._workers = set()
            self._pick_timer = None
            self._picker = None
            self.setWindowTitle("3D Feature Viewer for UniProt")
            self.resize(540, 790)
            self.setMinimumWidth(500)

            outer = QtWidgets.QVBoxLayout(self); outer.setContentsMargins(0, 0, 0, 0)
            scroll = QtWidgets.QScrollArea(); scroll.setWidgetResizable(True); outer.addWidget(scroll)
            bodyw = QtWidgets.QWidget(); scroll.setWidget(bodyw)
            lay = QtWidgets.QVBoxLayout(bodyw); lay.setContentsMargins(8, 8, 8, 8); lay.setSpacing(8)

            # ── Structure ───────────────────────────────────────────────
            sg = QtWidgets.QGroupBox("Structure"); sgl = QtWidgets.QVBoxLayout(sg); lay.addWidget(sg)
            r1 = QtWidgets.QHBoxLayout(); sgl.addLayout(r1)
            r1.addWidget(QtWidgets.QLabel("UniProt"))
            self.uid_edit = QtWidgets.QLineEdit(_STATE.get("uid") or ""); r1.addWidget(self.uid_edit, 1)
            self.fetch_btn = QtWidgets.QPushButton("Fetch"); r1.addWidget(self.fetch_btn)
            self.fetch_btn.clicked.connect(self.on_fetch)
            self.num_btn = QtWidgets.QToolButton(); self.num_btn.setText("№")
            self.num_btn.setToolTip("Residue numbering (for structures/trajectories you loaded yourself)")
            self.num_btn.clicked.connect(self.open_numbering_dialog); r1.addWidget(self.num_btn)
            r2 = QtWidgets.QHBoxLayout(); sgl.addLayout(r2)
            self.struct_combo = QtWidgets.QComboBox(); r2.addWidget(self.struct_combo, 1)
            self.load_btn = QtWidgets.QPushButton("Load"); r2.addWidget(self.load_btn)
            self.load_btn.clicked.connect(self.on_load_selected)
            self.loadall_btn = QtWidgets.QPushButton("Load all")
            self.loadall_btn.setToolTip("Download and load every listed structure")
            self.loadall_btn.clicked.connect(self.on_load_all); r2.addWidget(self.loadall_btn)
            self.align_btn = QtWidgets.QPushButton("Align"); self.align_btn.setToolTip("Superpose all loaded structures")
            self.align_btn.clicked.connect(lambda: ufv_align()); r2.addWidget(self.align_btn)
            hr = QtWidgets.QHBoxLayout(); sgl.addLayout(hr)
            hr.addWidget(QtWidgets.QLabel("Loaded structures"))
            self.obj_label = QtWidgets.QLabel(""); self.obj_label.setStyleSheet("color:#357; font-size:11px;")
            hr.addWidget(self.obj_label, 1)
            self.cb_apply_all = QtWidgets.QCheckBox("apply to all")
            self.cb_apply_all.setToolTip("Apply layers / colouring / zoom to every loaded structure (e.g. aligned copies)")
            hr.addWidget(self.cb_apply_all)
            self.obj_list = QtWidgets.QListWidget(); self.obj_list.setMaximumHeight(96)
            sgl.addWidget(self.obj_list)
            self.status = QtWidgets.QLabel(""); self.status.setStyleSheet("color:#a33; font-size:11px;")
            sgl.addWidget(self.status)

            # ── Layers ──────────────────────────────────────────────────
            lg = QtWidgets.QGroupBox("Layers (markers)"); lgl = QtWidgets.QVBoxLayout(lg); lay.addWidget(lg)
            lr = QtWidgets.QHBoxLayout(); lgl.addLayout(lr)
            self.cb_ptm = QtWidgets.QCheckBox("PTMs"); self.cb_site = QtWidgets.QCheckBox("Sites")
            self.cb_lig = QtWidgets.QCheckBox("Ligands")
            for cb in (self.cb_ptm, self.cb_site, self.cb_lig):
                lr.addWidget(cb)
            lr.addStretch(1)
            self.cb_ptm.clicked.connect(lambda: self.toggle_points("ptm", self.cb_ptm, ufv_ptms))
            self.cb_site.clicked.connect(lambda: self.toggle_points("site", self.cb_site, ufv_sites))
            self.cb_lig.clicked.connect(self.toggle_ligands)
            vr = QtWidgets.QHBoxLayout(); lgl.addLayout(vr)
            self.cb_var = QtWidgets.QCheckBox("Variants"); vr.addWidget(self.cb_var)
            self.cb_reviewed = QtWidgets.QCheckBox("reviewed only"); vr.addWidget(self.cb_reviewed)
            vr.addStretch(1)
            cr = QtWidgets.QHBoxLayout(); lgl.addLayout(cr)
            self.cons_cb = {}
            for tok, lbl, col in [("pathogenic", "Path", "#ef5350"), ("deleterious", "Delet", "#ffa726"),
                                  ("benign", "Benign", "#66bb6a"), ("uncertain", "Uncert", "#9e9e9e")]:
                c = QtWidgets.QCheckBox(lbl); c.setChecked(tok in ("pathogenic", "deleterious"))
                c.setStyleSheet("color:%s;" % col); cr.addWidget(c); self.cons_cb[tok] = c
                c.clicked.connect(self.on_var_filter)
            cr.addStretch(1)
            self.cb_var.clicked.connect(self.refresh_variants)
            self.cb_reviewed.clicked.connect(self.on_var_filter)

            # ── Cartoon colouring ───────────────────────────────────────
            cg = QtWidgets.QGroupBox("Cartoon colouring"); cgl = QtWidgets.QHBoxLayout(cg); lay.addWidget(cg)
            self.cart = QtWidgets.QComboBox()
            self.cart.addItems(["None", "Domains", "Topology", "pLDDT", "B-factor", "AlphaMissense",
                                "Burden", "Hotspots", "Contact hubs", "Constraint pocket"])
            cgl.addWidget(self.cart)
            self.cart.currentIndexChanged.connect(lambda _i: self.apply_cartoon())

            # ── Annotations & residue report ────────────────────────────
            ag = QtWidgets.QGroupBox("Annotations & residue report"); agl = QtWidgets.QVBoxLayout(ag)
            lay.addWidget(ag, 1)
            tr = QtWidgets.QHBoxLayout(); agl.addLayout(tr)
            self.list_kind = QtWidgets.QComboBox()
            self.list_kind.addItems(["PTMs", "Variants", "Sites", "Domains", "Ligands"])
            self.list_kind.currentIndexChanged.connect(lambda _i: self.refresh_table())
            tr.addWidget(self.list_kind)
            self.filter_edit = QtWidgets.QLineEdit(); self.filter_edit.setPlaceholderText("filter…")
            # Debounce: rebuilding a multi-thousand-row table on every keystroke is what made the
            # filter feel laggy. Coalesce keystrokes and refresh once the user pauses (250 ms).
            self._filter_timer = QtCore.QTimer(self); self._filter_timer.setSingleShot(True)
            self._filter_timer.setInterval(250)
            self._filter_timer.timeout.connect(self.refresh_table)
            self.filter_edit.textChanged.connect(lambda _t: self._filter_timer.start())
            tr.addWidget(self.filter_edit, 1)

            self.table = QtWidgets.QTableWidget(0, 2)
            self.table.setHorizontalHeaderLabels(["Pos", "Annotation"])
            self.table.horizontalHeader().setStretchLastSection(True)
            self.table.verticalHeader().setVisible(False)
            self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.table.setMinimumHeight(150)
            self.table.cellClicked.connect(self.on_row)            # click a row = zoom in
            self.table.cellDoubleClicked.connect(self.on_row)
            agl.addWidget(self.table, 1)

            hint = QtWidgets.QLabel("Click a row to zoom in · turn on Pick 3D and click an atom in the viewer")
            hint.setStyleSheet("color:#888; font-size:10px;"); agl.addWidget(hint)
            br = QtWidgets.QHBoxLayout(); agl.addLayout(br)
            self.show_filt_btn = QtWidgets.QPushButton("Show"); self.show_filt_btn.setToolTip("Show all filtered rows as markers")
            self.show_filt_btn.clicked.connect(self.show_filtered); br.addWidget(self.show_filt_btn)
            self.hide_filt_btn = QtWidgets.QPushButton("Hide"); self.hide_filt_btn.setToolTip("Hide the filtered-row markers")
            self.hide_filt_btn.clicked.connect(self.hide_filtered); br.addWidget(self.hide_filt_btn)
            self.cb_pick = QtWidgets.QCheckBox("Pick 3D"); br.addWidget(self.cb_pick)
            self.cb_pick.setToolTip("Click any atom in the PyMOL viewer to zoom to its residue/ligand")
            self.cb_pick.clicked.connect(self.toggle_pick)
            br.addStretch(1)
            self.reset_btn = QtWidgets.QPushButton("Reset view"); br.addWidget(self.reset_btn)
            self.reset_btn.clicked.connect(lambda: ufv_resetview(self.obj()))

            self.detail = QtWidgets.QTextBrowser(); self.detail.setMinimumHeight(150)
            self.detail.setOpenLinks(False); self.detail.setOpenExternalLinks(False)
            self.detail.anchorClicked.connect(self.on_anchor)
            agl.addWidget(self.detail)

            self.info = QtWidgets.QLabel(""); self.info.setWordWrap(True); self.info.setStyleSheet("font-size:11px; color:#555;")
            lay.addWidget(self.info)
            self.clear_btn = QtWidgets.QPushButton("Clear all overlays"); lay.addWidget(self.clear_btn)
            self.clear_btn.clicked.connect(self.on_clear)

            self._QtGui = QtGui
            self._last_pos = None
            self._sel = None               # unified current selection: {"kind":"res","pos":n} or
                                           # {"kind":"lig","resn":r} — drives 'Zoom selected'
            self._last_rep = None          # cached residue_report dict for re-rendering (evidence toggle)
            self._show_evidence = False    # variant-evidence fold state in the detail panel
            self._lig_resn = None          # ligand component currently shown
            self._lig_copies = []          # [(chain, resi), ...] occurrences of that component
            self._lig_idx = 0              # focused occurrence
            self._lig_chem = {}            # cached chemistry for the current ligand
            if self.cur_uid() and self.cur_uid() in _CACHE:
                self.after_fetch()

        # ---- helpers ----
        def cur_uid(self):
            return self.uid_edit.text().strip().upper() or _STATE.get("uid")

        def obj(self):
            return _STATE.get("obj") or _resolve_object(None)

        def _async(self, work, apply, status):
            self.status.setText(status); self.setEnabled(False)
            wk = _Worker(work)

            def on_done(res):
                self.setEnabled(True); self.status.setText(""); self._workers.discard(wk)
                if isinstance(res, Exception):
                    self.status.setText("Error: %s" % res); return
                if apply:
                    apply(res)
            wk.done.connect(on_done); self._workers.add(wk); wk.start()

        def update_info(self):
            ann = _CACHE.get(self.cur_uid() or "")
            if not ann:
                self.info.setText("Not fetched."); return
            cons = Counter(_cons_token(v["consequence"]) for v in ann["variants"])
            self.info.setText("%s · %d aa · PTM %d · Var %d (path %d) · Site %d · Dom %d · AM %s"
                              % (ann["uid"], len(ann["sequence"]), len(ann["ptms"]), len(ann["variants"]),
                                 cons.get("pathogenic", 0), len(ann["sites"]), len(ann["domains"]),
                                 "y" if ann["amMean"] else "n"))

        def loaded_objs(self):
            return [o for o in (cmd.get_object_list() or []) if not o.startswith("ufv_")]

        def target_objs(self):
            """Objects layer/colour/zoom operations act on: all loaded if 'apply to all', else active."""
            if self.cb_apply_all.isChecked():
                return self.loaded_objs()
            o = self.obj()
            return [o] if o else []

        def refresh_objs(self):
            active = self.obj()
            self.obj_list.clear()
            for o in self.loaded_objs():
                item = QtWidgets.QListWidgetItem()
                w = self._obj_row(o, active)
                item.setSizeHint(w.sizeHint())
                self.obj_list.addItem(item)
                self.obj_list.setItemWidget(item, w)
        obj_label_refresh = refresh_objs

        def _obj_row(self, name, active):
            w = QtWidgets.QWidget(); h = QtWidgets.QHBoxLayout(w); h.setContentsMargins(2, 0, 2, 0); h.setSpacing(4)
            vis = QtWidgets.QCheckBox(); vis.setChecked(True); vis.setToolTip("Show / hide")
            vis.clicked.connect(lambda checked, n=name: self.toggle_visible(n, checked))
            h.addWidget(vis)
            btn = QtWidgets.QPushButton(name); btn.setFlat(True)
            btn.setStyleSheet("text-align:left; border:none;" + ("font-weight:bold; color:#1565c0;" if name == active else ""))
            btn.setToolTip("Click to make active")
            btn.clicked.connect(lambda _=False, n=name: self.set_active(n))
            h.addWidget(btn, 1)
            x = QtWidgets.QToolButton(); x.setText("✕"); x.setToolTip("Remove structure")
            x.clicked.connect(lambda _=False, n=name: self.remove_structure(n))
            h.addWidget(x)
            return w

        def set_active(self, name):
            _STATE["obj"] = name
            self.update_colour_modes()
            ann = _CACHE.get(self.cur_uid()) or {}
            cov = actual_coverage(name, len(ann.get("sequence", "")))
            self.obj_label.setText(("%s · %.0f%% modelled" % (name, cov)) if cov is not None else name)
            self.refresh_objs(); self.refresh_table()

        def toggle_visible(self, name, checked):
            try:
                (cmd.enable if checked else cmd.disable)(name)
                for tag in ("ptm", "var", "site", "dom", "filt"):
                    lobj = _ufv_layer_obj(name, tag)
                    if lobj in (cmd.get_object_list() or []):
                        (cmd.enable if checked else cmd.disable)(lobj)
            except Exception:
                pass

        def remove_structure(self, name):
            if not name:
                return
            ufv_clear(name)
            try:
                cmd.delete(name)
            except Exception:
                pass
            (_STATE.get("sources") or {}).pop(name, None)
            _GEOM_CACHE.pop(name, None)
            if _STATE.get("obj") == name:
                _STATE["obj"] = (self.loaded_objs() or [None])[0]
            self.refresh_objs()

        # ---- fetch / structures ----
        def on_fetch(self):
            uid = self.cur_uid()
            if not uid:
                self.status.setText("Enter an accession."); return
            self._async(lambda: (fetch_annotations(uid), list_structures(uid)),
                        lambda res: self.after_fetch(res[1]), "Fetching %s ..." % uid)

        def after_fetch(self, structs=None):
            structs = structs if structs is not None else (_STATE.get("structures") or [])
            _STATE["structures"] = structs
            self.struct_combo.clear()
            for s in structs:
                self.struct_combo.addItem(s["label"], s)
            self.update_info(); self.refresh_table(); self.obj_label_refresh()

        def on_load_selected(self):
            s = self.struct_combo.currentData()
            if not s:
                self.status.setText("Fetch, then choose a structure."); return
            uid = self.cur_uid()
            self._async(lambda: _prepare_structure(s, uid),
                        lambda prep: self._finish_load(s, prep, uid), "Downloading %s ..." % s["label"])

        def _finish_load(self, struct, prep, uid):
            if not prep or not prep.get("path"):
                self.status.setText("Download failed."); return
            name = _load_from_path(struct, prep["path"], uid, prep.get("segments"))
            for cb in (self.cb_ptm, self.cb_site, self.cb_lig, self.cb_var):
                cb.setChecked(False)
            self.cart.setCurrentIndex(0)
            self.set_active(name)

        def on_load_all(self):
            structs = _STATE.get("structures") or []
            if not structs:
                self.status.setText("Fetch, then load."); return
            uid = self.cur_uid()

            def work():
                out = []
                for s in structs:
                    try:
                        out.append((s, _prepare_structure(s, uid)))
                    except Exception:
                        out.append((s, None))
                return out

            def apply(results):
                first = None
                for s, prep in results:
                    if prep and prep.get("path"):
                        try:
                            nm = _load_from_path(s, prep["path"], uid, prep.get("segments"))
                            first = first or nm
                        except Exception:
                            pass
                for cb in (self.cb_ptm, self.cb_site, self.cb_lig, self.cb_var):
                    cb.setChecked(False)
                self.cart.setCurrentIndex(0)
                if first:
                    self.set_active(first)
                self.refresh_objs()
                self.status.setText("Loaded %d of %d structures." % (len(self.loaded_objs()), len(results)))
            self._async(work, apply, "Downloading %d structures ..." % len(structs))

        def update_colour_modes(self):
            predicted = _structure_is_predicted(self.obj())
            model = self.cart.model()
            for i in range(self.cart.count()):
                txt = self.cart.itemText(i)
                ok = True
                if txt == "pLDDT":
                    ok = predicted
                elif txt == "B-factor":
                    ok = not predicted
                model.item(i).setEnabled(ok)

        # ---- layers ----
        def toggle_points(self, tag, cb, fn):
            objs = self.target_objs()
            if not objs or not _CACHE.get(self.cur_uid()):
                cb.setChecked(False); self.status.setText("Fetch + load a structure first."); return
            for o in objs:
                fn(o, self.cur_uid()) if cb.isChecked() else _hide_layer(o, tag)

        def toggle_ligands(self):
            objs = self.target_objs()
            if not objs:
                self.cb_lig.setChecked(False); return
            for o in objs:
                ufv_ligands(o) if self.cb_lig.isChecked() else ufv_ligands_hide(o)

        def on_var_filter(self):
            if not self.cb_var.isChecked():
                self.cb_var.setChecked(True)
            self.refresh_variants()

        def refresh_variants(self):
            objs = self.target_objs()
            if not objs:
                self.cb_var.setChecked(False); return
            if not self.cb_var.isChecked():
                for o in objs:
                    _hide_layer(o, "var")
                return
            ann = _CACHE.get(self.cur_uid())
            if not ann:
                self.cb_var.setChecked(False); self.status.setText("Fetch first."); return
            toks = {t for t, c in self.cons_cb.items() if c.isChecked()}
            groups = _variant_groups(ann, toks or None, self.cb_reviewed.isChecked())
            for o in objs:
                _set_sphere_layer(o, "var", groups)

        def apply_cartoon(self):
            objs = self.target_objs(); choice = self.cart.currentText()
            if not objs:
                return
            if choice == "None":
                for o in objs:
                    _reset_cartoon(o)
                return
            uid = self.cur_uid()
            if choice not in ("pLDDT", "B-factor") and not _CACHE.get(uid):
                self.status.setText("Fetch first."); self.cart.setCurrentIndex(0); return
            cheap = {"Domains": ufv_domains, "Topology": ufv_topology, "pLDDT": ufv_plddt,
                     "B-factor": ufv_bfactor, "AlphaMissense": ufv_alphamissense, "Burden": ufv_burden}
            if choice in cheap:
                for o in objs:
                    cheap[choice](o, uid)
                return
            ann = _CACHE.get(uid) or {}
            if choice == "Constraint pocket" and not ann.get("amMean"):
                self.status.setText("Constraint pocket needs AlphaMissense data."); self.cart.setCurrentIndex(0); return
            geoms = {o: _ca_geometry(o) for o in objs}    # read coords on the UI thread
            if choice == "Hotspots":
                colors = {"strong": "#b71c1c", "moderate": "#e64a19", "weak": "#ffa726"}
                work = lambda: {o: _per_chain_merge(g, lambda gs: _compute_hotspots(gs, ann.get("variants", []))) for o, g in geoms.items()}
            elif choice == "Contact hubs":
                colors = {"strong": "#6a1b9a", "moderate": "#ab47bc"}
                work = lambda: {o: _per_chain_merge(g, _betweenness_hubs) for o, g in geoms.items()}
            else:
                colors = {"pocket": "#00897b", "exposed": "#8e24aa"}
                work = lambda: {o: _compute_pockets(g, ann["amMean"]) for o, g in geoms.items()}

            def apply(res):
                for o, tiers in res.items():
                    _set_cartoon_layer(o, choice.lower()); _color_tiers(o, tiers, colors)
                self.status.setText("%s applied to %d structure(s)." % (choice, len(res)))
            self._async(work, apply, "Computing %s ..." % choice)

        # ---- list table + report ----
        def refresh_table(self):
            ann = _CACHE.get(self.cur_uid())
            self.table.setRowCount(0)
            kind = self.list_kind.currentText()
            flt = self.filter_edit.text().strip().lower()
            rows = []
            if not ann and kind != "Ligands":
                self._filtered_rows = []
                return
            if kind == "PTMs":
                for p in ann["ptms"]:
                    rows.append((p["position"], "%s — %s" % (p["category"], p["description"]), p["color"]))
            elif kind == "Variants":
                for v in ann["variants"]:
                    rows.append((v["position"], "%s%d%s · %s%s" % (v.get("wildType", ""), v["position"],
                                 v.get("mutant", ""), v["consequence"],
                                 (" · " + "; ".join(v["diseases"])) if v.get("diseases") else ""),
                                 v["consequenceColor"]))
            elif kind == "Sites":
                for s in ann["sites"]:
                    rows.append((s["position"], s["description"], SITE_COLOR))
            elif kind == "Domains":
                for d in ann["domains"]:
                    rng = "%d-%d" % (d["position"], d["endPosition"]) if d["isRange"] else str(d["position"])
                    rows.append((d["position"], "%s · %s" % (rng, d["description"]), d["color"]))
            else:
                for resn, cnt in sorted(enumerate_ligands(self.obj() or "").items()):
                    nm = (_LIG_CACHE.get(resn) or {}).get("name") or "(click for chemistry)"
                    tag = " ×%d" % cnt if cnt > 1 else ""
                    rows.append((resn, "%s%s — %s" % (resn, tag, nm), "#2e7d32"))
            if flt:
                rows = [r for r in rows if flt in r[1].lower() or flt in str(r[0])]
            # Which UniProt positions are actually modelled in the active structure (grey out the rest).
            modelled = None
            o = self.obj()
            if o and kind != "Ligands":
                modelled = {g["uniPos"] for g in _ca_geometry(o) if g["uniPos"] is not None}
            self._filtered_rows = [(r[0], r[2]) for r in rows
                                   if modelled is None or not isinstance(r[0], int) or r[0] in modelled]
            # Cap how many rows we actually build as widgets. Creating thousands of QTableWidgetItems
            # (e.g. the full variant set) is slow and freezes the panel; the marker overlay and
            # "Show all filtered" still use the complete set above. Narrow with the filter box to see more.
            MAX_ROWS = 600
            shown = rows[:MAX_ROWS]
            self.table.setRowCount(len(shown))
            QtGui = self._QtGui
            for r, (pos, text, color) in enumerate(shown):
                it0 = QtWidgets.QTableWidgetItem(str(pos)); it1 = QtWidgets.QTableWidgetItem(text)
                unmod = modelled is not None and isinstance(pos, int) and pos not in modelled
                if unmod:
                    grey = QtGui.QBrush(QtGui.QColor("#bbbbbb"))
                    it0.setForeground(grey); it1.setForeground(grey)
                    it0.setFlags(QtCore.Qt.NoItemFlags); it1.setFlags(QtCore.Qt.NoItemFlags)
                    tip = text + "  (not modelled in this structure)"
                else:
                    brush = QtGui.QBrush(QtGui.QColor(color))
                    it0.setForeground(brush); it1.setForeground(brush)
                    tip = text
                it0.setToolTip(tip); it1.setToolTip(tip)
                self.table.setItem(r, 0, it0); self.table.setItem(r, 1, it1)
            if len(rows) > len(shown):
                self.status.setText("Showing %d of %d — use the filter to narrow." % (len(shown), len(rows)))

        def _row_pos(self, r):
            it = self.table.item(r, 0)
            return it.text() if it else None

        def on_row(self, r, _c):
            """Clicking a row zooms straight in (residue or ligand)."""
            val = self._row_pos(r)
            if val is None:
                return
            if self.list_kind.currentText() == "Ligands":
                self.on_ligand(val); return
            try:
                self.report_residue(int(val), focus=True)
            except ValueError:
                pass

        def report_residue(self, uni_pos, focus=False):
            objs, uid = self.target_objs(), self.cur_uid()
            if not objs or not uid:
                return
            self._last_pos = uni_pos
            self._sel = {"kind": "res", "pos": uni_pos}
            if focus:  # zoom FIRST so a slow/failing report compute can never swallow the double-click
                for o in objs:  # aligned copies overlap, so the last zoom frames them all
                    try:
                        ufv_focus(o, uni_pos)
                    except Exception as exc:
                        self.status.setText("Focus error: %s" % exc)
            try:
                self._last_rep = residue_report(objs[0], uid, uni_pos)
                self._lig_resn = None      # leaving ligand context
                self._show_evidence = False
                self.detail.setHtml(format_report_html(self._last_rep, self._show_evidence))
            except Exception as exc:
                self._last_rep = None
                self.detail.setHtml("<i>Report error: %s</i>" % exc)

        def on_anchor(self, url):
            """Handle in-report links: ufv:res:<n> (jump to residue), ufv:evi (toggle evidence),
            ufv:lig:prev / ufv:lig:next (circulate same-named ligand copies)."""
            s = url.toString() if hasattr(url, "toString") else str(url)
            if s.startswith("http://") or s.startswith("https://"):
                try:
                    import webbrowser
                    webbrowser.open(s)
                except Exception:
                    pass
                return
            if s.startswith("ufv:res:"):
                try:
                    self.report_residue(int(s.rsplit(":", 1)[1]), focus=True)
                except ValueError:
                    pass
            elif s == "ufv:evi":
                self._show_evidence = not self._show_evidence
                if self._last_rep:
                    self.detail.setHtml(format_report_html(self._last_rep, self._show_evidence))
            elif s in ("ufv:lig:prev", "ufv:lig:next") and self._lig_copies:
                self._lig_idx = (self._lig_idx + (1 if s.endswith("next") else -1)) % len(self._lig_copies)
                ch, resi = self._lig_copies[self._lig_idx]
                o = self.obj()
                if o:
                    ufv_ligand_focus(o, self._lig_resn, ch, resi)
                self._render_ligand()

        def on_ligand(self, resn, pick=None):
            o = self.obj()
            if not o:
                return
            self._lig_resn = resn
            self._sel = {"kind": "lig", "resn": resn}
            self._lig_copies = ligand_instances(o, resn)
            self._lig_idx = 0
            if pick is not None and pick in self._lig_copies:
                self._lig_idx = self._lig_copies.index(pick)
            if self._lig_copies:
                ch, resi = self._lig_copies[self._lig_idx]
                ufv_ligand_focus(o, resn, ch, resi)
            else:
                ufv_ligand_focus(o, resn)
            self.detail.setHtml("<b>%s</b> — fetching chemistry…" % resn)
            self._async(lambda: ligand_chemistry(o),
                        lambda chem: self.show_ligand(resn, chem), "Fetching ligand chemistry ...")

        def show_ligand(self, resn, chem):
            self._lig_chem = chem
            self._render_ligand()
            self.refresh_table()

        def _render_ligand(self):
            if not self._lig_resn:
                return
            info = (self._lig_chem or {}).get(self._lig_resn, {})
            self.detail.setHtml(format_ligand_html(self._lig_resn, info,
                                                   instances=len(self._lig_copies) or 1,
                                                   copies=self._lig_copies, cur=self._lig_idx))

        # ---- 3D picking ----
        # Two capture paths feed one handler so picking works regardless of PyMOL build:
        #   1. a native Wizard (do_select / do_pick) — fires the instant you click an atom;
        #   2. a QTimer that watches the "sele" PyMOL creates on a click — a fallback if (1) is quiet.
        # _handle_pick is debounced so the two paths never double-act, ALWAYS zooms (it does not depend
        # on an accession being fetched), and writes a status line so it's obvious the click registered.
        def toggle_pick(self):
            if self.cb_pick.isChecked():
                self._install_picker()
                t = getattr(self, "_pick_timer", None)
                if t is None:
                    t = QtCore.QTimer(self); t.setInterval(200)
                    t.timeout.connect(self._poll_pick); self._pick_timer = t
                try:
                    cmd.deselect()
                except Exception:
                    pass
                t.start()
                self.status.setText("Pick 3D on — click any atom (protein or ligand) to zoom in.")
            else:
                self._remove_picker()
                t = getattr(self, "_pick_timer", None)
                if t is not None:
                    t.stop()
                self.status.setText("")

        def _install_picker(self):
            try:
                from pymol.wizard import Wizard
            except Exception:
                self._picker = None
                return
            panel = self

            class _Picker(Wizard):
                def get_prompt(self):
                    return ["Click any atom to zoom to its residue / ligand."]

                def get_panel(self):
                    return [[1, "UFV pick", ""], [2, "Stop picking", "cmd.set_wizard()"]]

                def _grab(self, selection):
                    atoms = []
                    try:
                        m = self.cmd.get_model(selection)
                        atoms = list(getattr(m, "atom", []))
                    except Exception:
                        pass
                    try:
                        self.cmd.deselect()
                    except Exception:
                        pass
                    if atoms:
                        a = atoms[0]
                        panel._handle_pick(a.chain, a.resi, a.resn, int(getattr(a, "hetatm", 0) or 0))

                def do_select(self, selection):
                    self._grab(selection)
                    try:
                        self.cmd.refresh_wizard()
                    except Exception:
                        pass
                    return None

                def do_pick(self, bondFlag):
                    self._grab("pk1")
                    try:
                        self.cmd.unpick()
                    except Exception:
                        pass
                    return None

            try:
                self._picker = _Picker()
                cmd.set_wizard(self._picker)
            except Exception:
                self._picker = None

        def _remove_picker(self):
            if getattr(self, "_picker", None) is not None:
                try:
                    cmd.set_wizard()
                except Exception:
                    pass
                self._picker = None

        def _poll_pick(self):
            try:
                if cmd.count_atoms("sele") == 0:
                    return
            except Exception:
                return
            info = []
            try:
                cmd.iterate("sele", "info.append((chain, resi, resn, hetatm, model))", space={"info": info})
                cmd.deselect()      # clear the click so it isn't re-read and the pink dots disappear
            except Exception:
                return
            info = [t for t in info if t[4] != _UFV_LIGLABEL]   # ignore the focus label pseudoatom
            if info:
                ch, resi, resn, het, _model = info[0]
                self._handle_pick(ch, resi, resn, het)

        def _handle_pick(self, ch, resi, resn, het):
            import time
            key = (ch, resi, resn, het)
            now = time.time()
            if key == getattr(self, "_last_pick", None) and now - getattr(self, "_last_pick_t", 0) < 0.6:
                return                                  # debounce: the two capture paths firing together
            self._last_pick, self._last_pick_t = key, now
            o = self.obj()
            if not o:
                self.status.setText("Pick: load a structure first."); return
            if het:
                self.on_ligand(resn, pick=(ch, resi))
                self.status.setText("Picked ligand %s (%s/%s)" % (resn, ch, resi)); return
            uni = _resi_to_uni(o, ch, resi)
            if uni is None:
                self.status.setText("Picked %s/%s — not mapped to a UniProt position." % (ch, resi)); return
            self._sel = {"kind": "res", "pos": uni}; self._last_pos = uni
            try:                                        # ALWAYS zoom — independent of fetched annotations
                for ob in self.target_objs() or [o]:
                    ufv_focus(ob, uni)
            except Exception as exc:
                self.status.setText("Focus error: %s" % exc); return
            uid = self.cur_uid()
            if uid:                                     # detail panel is best-effort
                try:
                    self._last_rep = residue_report(o, uid, uni)
                    self._show_evidence = False
                    self.detail.setHtml(format_report_html(self._last_rep, self._show_evidence))
                except Exception:
                    pass
            self.status.setText("Picked residue %s%d" % (resn[:1] if resn else "", uni))

        # ---- filtered show/hide ----
        def show_filtered(self):
            o = self.obj()
            rows = [(p, c) for p, c in getattr(self, "_filtered_rows", []) if isinstance(p, int)]
            if not o or not rows:
                return
            groups = {}
            for pos, color in rows:
                groups.setdefault(color, []).append(pos)
            _set_sphere_layer(o, "filt", groups)

        def hide_filtered(self):
            o = self.obj()
            if o:
                _hide_layer(o, "filt")

        # ---- numbering dialog ----
        def open_numbering_dialog(self):
            o, uid = self.obj(), self.cur_uid()
            if not o or not uid:
                self.status.setText("Load/select an object first."); return
            dlg = QtWidgets.QDialog(self); dlg.setWindowTitle("Residue numbering — %s" % o)
            g = QtWidgets.QGridLayout(dlg)
            rb_id = QtWidgets.QRadioButton("Identity (resi == UniProt)"); rb_id.setChecked(True)
            rb_si = QtWidgets.QRadioButton("SIFTS, PDB:"); pdb = QtWidgets.QLineEdit(); pdb.setMaximumWidth(80)
            rb_man = QtWidgets.QRadioButton("Manual chain:")
            ch = QtWidgets.QLineEdit("A"); ch.setMaximumWidth(34)
            us = QtWidgets.QLineEdit(); us.setPlaceholderText("UniProt@"); us.setMaximumWidth(72)
            rs = QtWidgets.QLineEdit(); rs.setPlaceholderText("resi@"); rs.setMaximumWidth(58)
            re = QtWidgets.QLineEdit(); re.setPlaceholderText("end"); re.setMaximumWidth(52)
            g.addWidget(rb_id, 0, 0, 1, 5)
            g.addWidget(rb_si, 1, 0); g.addWidget(pdb, 1, 1)
            g.addWidget(rb_man, 2, 0); g.addWidget(ch, 2, 1); g.addWidget(us, 2, 2); g.addWidget(rs, 2, 3); g.addWidget(re, 2, 4)
            bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
            g.addWidget(bb, 3, 0, 1, 5)

            def apply_and_close():
                if rb_si.isChecked():
                    ufv_map(o, uid, "sifts", pdb.text().strip())
                elif rb_man.isChecked():
                    ufv_chain(o, ch.text().strip(), us.text().strip(), rs.text().strip() or None, re.text().strip() or None)
                else:
                    ufv_map(o, uid, "identity")
                _GEOM_CACHE.pop(o, None); dlg.accept()
            bb.accepted.connect(apply_and_close); bb.rejected.connect(dlg.reject)
            dlg.exec_()

        def on_clear(self):
            ufv_clear(self.obj())
            for c in (self.cb_ptm, self.cb_site, self.cb_lig, self.cb_var):
                c.setChecked(False)
            self.cart.setCurrentIndex(0); self.detail.setHtml("")

        def closeEvent(self, ev):
            t = getattr(self, "_pick_timer", None)
            if t is not None:
                t.stop()
            self._remove_picker()
            super(_UFVPanel, self).closeEvent(ev)

    return _UFVPanel


def _prepare_structure(struct, uid):
    data = _get_text(struct["url"])
    if not data:
        return None
    ext = "cif" if struct.get("fmt") == "mmcif" else "pdb"
    name = re.sub(r"[^A-Za-z0-9_]", "_", struct["key"])
    path = os.path.join(tempfile.gettempdir(), "ufv_%s.%s" % (name, ext))
    with open(path, "w") as fh:
        fh.write(data)
    segments = None
    if struct.get("numbering") == "sifts" and struct.get("pdbId"):
        segments = _sifts_segments(struct["pdbId"], uid.upper()) or None
    return {"path": path, "segments": segments}


def _load_from_path(struct, path, uid, segments=None):
    name = re.sub(r"[^A-Za-z0-9_]", "_", struct["key"])
    cmd.load(path, name)
    try:
        os.remove(path)
    except OSError:
        pass
    cmd.hide("everything", name)
    cmd.show("cartoon", name)
    cmd.color("gray80", name)
    if segments:
        UFV_MAPS[name] = {"mode": "sifts", "segments": segments}
    elif struct.get("numbering") == "sifts" and struct.get("pdbId"):
        _set_sifts_map(name, struct["pdbId"], uid)
    else:
        _set_identity_map(name)
    _GEOM_CACHE.pop(name, None)
    _STATE["obj"] = name
    _STATE.setdefault("sources", {})[name] = struct.get("source", "PDB")
    cmd.orient(name)
    return name


def __init_plugin__(app=None):
    """PyMOL Plugin Manager entry point."""
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt("3D Feature Viewer for UniProt", ufv_gui)
    except Exception:
        pass


if _HAS_PYMOL:
    cmd.extend("ufv_gui", ufv_gui)
    print("[UFV] 3D Feature Viewer for UniProt loaded. Open the panel with: ufv_gui  |  commands: "
          "ufv_fetch, ufv_structures, ufv_use, ufv_load, ufv_map, ufv_chain, ufv_ptms, ufv_variants, "
          "ufv_sites, ufv_domains, ufv_topology, ufv_alphamissense, ufv_burden, ufv_plddt, ufv_bfactor, "
          "ufv_hotspots, ufv_contacthubs, ufv_pockets, ufv_ligands, ufv_report, ufv_focus, "
          "ufv_resetview, ufv_align, ufv_hide, ufv_clear, ufv_info")


# ----------------------------------------------------------------------------------------------
# Command-line backend for the VMD plugin: emit annotations / SIFTS maps as Tcl, or download a
# model.  All progress output is suppressed so stdout carries only the requested payload, which
# the Tcl side captures with `exec`.
# ----------------------------------------------------------------------------------------------
def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _emit_tcl_annotations(uid):
    ann = _quiet(fetch_annotations, uid)
    out = ["set ::ufv::a_seqlen %d" % len(ann["sequence"])]

    ptm_items = []
    for p in ann["ptms"]:
        ptm_items.append("{%d %s}" % (p["position"], p["color"]))
        if p["endPosition"] != p["position"]:
            ptm_items.append("{%d %s}" % (p["endPosition"], p["color"]))
    out.append("set ::ufv::a_ptms { %s }" % " ".join(ptm_items))

    def _cons_token(c):
        c = c.lower()
        if "pathogenic" in c:
            return "pathogenic"
        if "benign" in c:
            return "benign"
        if "deleterious" in c:
            return "deleterious"
        return "uncertain"
    var_items = ["{%d %s %s}" % (v["position"], v["consequenceColor"], _cons_token(v["consequence"]))
                 for v in ann["variants"]]
    out.append("set ::ufv::a_variants { %s }" % " ".join(var_items))

    site_items = []
    for s in ann["sites"]:
        site_items.append("{%d %s}" % (s["position"], s["color"]))
        if s["endPosition"] != s["position"]:
            site_items.append("{%d %s}" % (s["endPosition"], s["color"]))
    out.append("set ::ufv::a_sites { %s }" % " ".join(site_items))

    dom_items = ["{%d %d %s %d}" % (d["position"], d["endPosition"], d["color"], 1 if d["isRange"] else 0)
                 for d in ann["domains"]]
    out.append("set ::ufv::a_domains { %s }" % " ".join(dom_items))

    topo_items = ["{%d %d %s}" % (t["start"], t["end"], t["color"]) for t in ann["topology"]]
    out.append("set ::ufv::a_topology { %s }" % " ".join(topo_items))

    am_items = ["{%d %.4f}" % (pos, avg) for pos, avg in sorted(ann["amMean"].items())]
    out.append("set ::ufv::a_am { %s }" % " ".join(am_items))

    sys.stdout.write("\n".join(out) + "\n")


def _emit_tcl_sifts(pdb_id, uid):
    segs = _quiet(_sifts_segments, pdb_id, uid.upper())
    out = ["set ::ufv::s_chains { %s }" % " ".join(segs.keys())]
    for ch, segl in segs.items():
        items = ["{%s %s %s %s}" % (s["u_start"], s["u_end"], s["pdb_start"], s["pdb_end"]) for s in segl]
        out.append("set ::ufv::s_%s { %s }" % (ch, " ".join(items)))
    sys.stdout.write("\n".join(out) + "\n")


def _download_model(uid):
    url = _quiet(_alphafold_url, uid)
    if not url:
        sys.exit("no AlphaFold model for %s" % uid)
    data = _quiet(_get_text, url)
    if not data:
        sys.exit("download failed for %s" % uid)
    path = os.path.join(tempfile.gettempdir(), "ufv_%s.pdb" % uid.upper())
    with open(path, "w") as fh:
        fh.write(data)
    sys.stdout.write(path + "\n")


if __name__ == "__main__":
    _args = sys.argv[1:]
    if len(_args) >= 2 and _args[0] == "--emit-tcl":
        _emit_tcl_annotations(_args[1])
    elif len(_args) >= 3 and _args[0] == "--emit-sifts":
        _emit_tcl_sifts(_args[1], _args[2])
    elif len(_args) >= 2 and _args[0] == "--download":
        _download_model(_args[1])
    else:
        sys.stderr.write(
            "3D Feature Viewer for UniProt — backend for the VMD plugin\n"
            "usage:\n"
            "  python ufv_pymol.py --emit-tcl   <uniprot_id>\n"
            "  python ufv_pymol.py --emit-sifts <pdb_id> <uniprot_id>\n"
            "  python ufv_pymol.py --download   <uniprot_id>\n")
        sys.exit(2)
