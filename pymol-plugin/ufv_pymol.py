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


def _betweenness_hubs(geom, threshold=8.0):
    """Port of analysis.js betweennessHubs: Brandes' betweenness on the Cα contact graph, tiers by
    absolute z-score. Returns {uniPos: 'strong'|'moderate'}."""
    nodes = [g for g in geom if g["uniPos"] is not None]
    n = len(nodes)
    if n < 8 or n > HUB_MAX_CA:
        return {}
    th2 = threshold * threshold
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if _d2(nodes[i], nodes[j]) <= th2:
                adj[i].append(j)
                adj[j].append(i)
    bc = [0.0] * n
    for s in range(n):
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
        neigh = [[i] for i in range(k)]
        for i in range(k):
            for j in range(i + 1, k):
                if _d2(items[i][2], items[j][2]) <= th2:
                    neigh[i].append(j)
                    neigh[j].append(i)
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
    neigh = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if _d2(universe[i], universe[j]) <= th2:
                neigh[i].append(j)
                neigh[j].append(i)
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


def _am_profile(am_map, pos):
    """Per-substitution AlphaMissense scores at a position: [(wt, mut, score)], highest first."""
    out = []
    for key, score in am_map.items():
        mt = re.match(r"^([A-Z])(\d+)([A-Z])$", key)
        if mt and int(mt.group(2)) == pos:
            out.append((mt.group(1), mt.group(3), score))
    out.sort(key=lambda t: -t[2])
    return out


# ----------------------------------------------------------------------------------------------
# Fetch orchestration
# ----------------------------------------------------------------------------------------------
def fetch_annotations(uid, force=False):
    """Fetch and cache every annotation layer for a UniProt accession."""
    uid = uid.strip().upper()
    if not force and uid in _CACHE:
        return _CACHE[uid]
    print("[UFV] fetching annotations for %s ..." % uid)
    features = _get_json(FEATURES_URL.format(uid))
    variation = _get_json(VARIATION_URL.format(uid))
    uniprot = _get_json(UNIPROT_URL.format(uid))
    am_map = _parse_am_csv(_get_text(AM_CSV_URL.format(uid)))

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


def list_structures(uid, seqlen=0):
    """Available structures for an accession (for the structure selector): the canonical AlphaFold
    model, experimental PDB chains (PDBe best_structures, numbered via SIFTS on load), and computed
    models (3D-Beacons). Each item carries how to number it. Mirrors api.js getStructures (the
    common cases)."""
    uid = uid.upper()
    out = []
    af = _alphafold_url(uid)
    if af:
        out.append({"key": "AF-%s" % uid, "label": "AlphaFold — predicted, full length",
                    "source": "AlphaFold", "url": af, "fmt": "pdb",
                    "pdbId": None, "chainId": None, "numbering": "identity", "cov": 100.0})
    best = _get_json(PDBe_BEST.format(uid)) or {}
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
        exp.append({"key": "%s_%s" % (pid, ch),
                    "label": "%s chain %s%s — %s" % (pid, ch, " (%.0f%%)" % covpct if covpct else "",
                                                     it.get("experimental_method", "experimental")),
                    "source": "PDB", "url": PDBe_PDB.format(pid.lower()), "fmt": "pdb",
                    "pdbId": pid, "chainId": ch, "numbering": "sifts", "cov": covpct,
                    "_res": it.get("resolution") or 99.0})
    # Proteins like insulin map to hundreds of PDB chains; sort best-coverage/-resolution first
    # and cap the list so the selector stays usable.
    exp.sort(key=lambda s: (-(s["cov"] or 0), s["_res"]))
    out.extend(exp[:STRUCTURE_LIST_CAP])
    summ = _get_json(BEACONS_SUMMARY.format(uid)) or {}
    for e in (summ.get("structures") or []):
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
def _color_name(hex_str):
    h = hex_str.lstrip("#")
    name = "ufv_%s" % h
    cmd.set_color(name, [int(h[0:2], 16) / 255.0, int(h[2:4], 16) / 255.0, int(h[4:6], 16) / 255.0])
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


def _sel_for_positions(obj, positions, ca_only=False):
    """Build a compact PyMOL selection covering the given UniProt positions on the mapped object."""
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
    sel = "({}) and ({})".format(obj, " or ".join(clauses))
    if ca_only:
        sel += " and name CA"
    return sel


# Selection-free overlay model. Point (sphere) layers are kept as state only; redrawing hides the
# previously-shown UFV spheres and re-shows the active layers. We colour the *representation*
# (sphere_color), never the atom, so cartoon colouring (which follows atom colour) is untouched —
# this is what fixes "showing variants recolours the ribbon". No named selections are created, so
# the PyMOL object panel stays clean.
_SPHERE_STATE = {}    # obj -> { tag -> {color_hex: [uniprot_positions]} }  (active sphere layers)
_SPHERE_SHOWN = {}    # obj -> last inline selection of shown UFV spheres (so we can hide them)
_CARTOON_LAYER = {}   # obj -> active cartoon-colouring mode (informational)
_SPHERE_ORDER = ["site", "dom", "ptm", "var"]  # later tags win colour on a shared Ca atom
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
    """Suspend PyMOL scene updates while running many cmd ops, then rebuild once. Without this,
    every show/set/color triggers a full rebuild — painting thousands of variant spheres one
    selection at a time froze the viewer."""
    try:
        cmd.set("suspend_updates", "on")
    except Exception:
        pass
    try:
        yield
    finally:
        try:
            cmd.set("suspend_updates", "off")
            cmd.refresh()
        except Exception:
            pass


_SPHERE_TUNED = set()  # objects we've already set fast sphere settings on


def _tune_spheres(obj):
    if obj in _SPHERE_TUNED:
        return
    # GPU impostor spheres render thousands of points cheaply; avoids the slow CPU sphere geometry.
    for k, v in (("sphere_mode", 9), ("sphere_quality", 1), ("sphere_scale", 0.5)):
        try:
            cmd.set(k, v, obj)
        except Exception:
            pass
    _SPHERE_TUNED.add(obj)


def _redraw_spheres(obj):
    """Hide the previously-shown UFV spheres, then re-show every active sphere layer (batched)."""
    _tune_spheres(obj)
    with _batch():
        prev = _SPHERE_SHOWN.get(obj)
        if prev:
            try:
                cmd.hide("spheres", prev)
            except Exception:
                pass
        state = _SPHERE_STATE.get(obj, {})
        shown = set()
        for tag in _SPHERE_ORDER:
            for color, positions in (state.get(tag) or {}).items():
                sel = _sel_for_positions(obj, positions, ca_only=True)
                if not sel:
                    continue
                cmd.show("spheres", sel)
                cmd.set("sphere_color", _color_name(color), sel)
                shown.update(positions)
        _SPHERE_SHOWN[obj] = _sel_for_positions(obj, sorted(shown), ca_only=True) if shown else None
        return len(shown)


def _set_sphere_layer(obj, tag, groups):
    """Turn a sphere layer on (groups = {color: [positions]}) or off (groups = None)."""
    st = _SPHERE_STATE.setdefault(obj, {})
    if groups:
        st[tag] = groups
    else:
        st.pop(tag, None)
    return _redraw_spheres(obj)


def _hide_layer(obj, tag):
    _set_sphere_layer(obj, tag, None)


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


def ufv_domains(obj=None, uid=None):
    """ufv_domains [object [, uniprot_id]] — colour the cartoon by domain / region / repeat."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    _set_cartoon_layer(obj, "domains")
    dom_spheres = {}
    n = 0
    with _batch():
        for d in ann["domains"]:
            positions = list(range(d["position"], d["endPosition"] + 1))
            if d["isRange"]:
                sel = _sel_for_positions(obj, positions)
                if not sel:
                    continue
                cmd.color(_color_name(d["color"]), sel)  # cartoon = atom colour (intended here)
            else:
                dom_spheres.setdefault(d["color"], []).append(d["position"])
            n += 1
    _set_sphere_layer(obj, "dom", dom_spheres or None)
    print("[UFV] %s: coloured %d domain/region features." % (obj, n))


def ufv_topology(obj=None, uid=None):
    """ufv_topology [object [, uniprot_id]] — colour the cartoon by membrane topology."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    _set_cartoon_layer(obj, "topology")
    n = 0
    with _batch():
        for t in ann["topology"]:
            sel = _sel_for_positions(obj, list(range(t["start"], t["end"] + 1)))
            if not sel:
                continue
            cmd.color(_color_name(t["color"]), sel)
            n += 1
    print("[UFV] %s: coloured %d topology segments." % (obj, n))


def ufv_alphamissense(obj=None, uid=None):
    """ufv_alphamissense [object [, uniprot_id]] — colour by mean AlphaMissense pathogenicity."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    am = ann["amMean"]
    if not am:
        print("[UFV] no AlphaMissense scores available for %s." % ann["uid"])
        return
    _set_cartoon_layer(obj, "alphamissense")
    buckets = {"#3d85c8": [], "#b9c2cf": [], "#e06666": [], "#b71c1c": []}
    for pos, avg in am.items():
        c = "#b71c1c" if avg >= 0.78 else "#e06666" if avg >= 0.564 else "#b9c2cf" if avg >= 0.34 else "#3d85c8"
        buckets[c].append(pos)
    with _batch():
        cmd.color(_color_name("#b9c2cf"), obj)
        for color, positions in buckets.items():
            sel = _sel_for_positions(obj, positions)
            if sel:
                cmd.color(_color_name(color), sel)
    print("[UFV] %s: coloured by AlphaMissense (%d scored positions)." % (obj, len(am)))


def ufv_burden(obj=None, uid=None):
    """ufv_burden [object [, uniprot_id]] — colour mutation/phenotype burden-positive residues."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    _set_cartoon_layer(obj, "burden")
    sel = _sel_for_positions(obj, sorted(ann["burden"]))
    if sel:
        cmd.color(_color_name("#e65100"), sel)
    print("[UFV] %s: %d burden-positive residues." % (obj, len(ann["burden"])))


def ufv_plddt(obj=None, uid=None):
    """ufv_plddt [object] — colour by AlphaFold pLDDT (the B-factor column of an AF model)."""
    obj = _resolve_object(obj)
    _set_cartoon_layer(obj, "plddt")
    cmd.color(_color_name("#ff7d45"), "(%s) and b<50" % obj)
    cmd.color(_color_name("#ffdb13"), "(%s) and b>=50 and b<70" % obj)
    cmd.color(_color_name("#65cbf3"), "(%s) and b>=70 and b<90" % obj)
    cmd.color(_color_name("#0053d6"), "(%s) and b>=90" % obj)
    print("[UFV] %s: coloured by pLDDT." % obj)


def ufv_bfactor(obj=None, uid=None):
    """ufv_bfactor [object] — colour by B-factor (blue low → red high)."""
    obj = _resolve_object(obj)
    _set_cartoon_layer(obj, "bfactor")
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


def _color_tiers(obj, tiers, colors):
    groups = {}
    for uni, tier in tiers.items():
        c = colors.get(tier)
        if c:
            groups.setdefault(c, []).append(uni)
    with _batch():
        for c, positions in groups.items():
            sel = _sel_for_positions(obj, positions)
            if sel:
                cmd.color(_color_name(c), sel)
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
    return {
        "uid": ann["uid"], "pos": uni_pos, "aa": aa,
        "ptms": [p for p in ann["ptms"] if p["position"] <= uni_pos <= p["endPosition"]],
        "variants": var_pos.get(uni_pos, []),
        "sites": [s for s in ann["sites"] if s["position"] <= uni_pos <= s["endPosition"]],
        "domains": [d for d in ann["domains"] if d["position"] <= uni_pos <= d["endPosition"]],
        "amMean": ann["amMean"].get(uni_pos),
        "amProfile": _am_profile(ann["amMap"], uni_pos)[:6],
        "nearby": near,
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
    if rep["nearby"]:
        near = ", ".join("%d (%.1fÅ%s)" % (n["pos"], n["dist"], " " + "/".join(n["tags"]) if n["tags"] else "")
                         for n in rep["nearby"][:12])
        L.append("  Nearby ≤12Å: " + near)
    return "\n".join(L)


def ufv_focus(obj=None, uni_pos=None):
    """ufv_focus [object,] uniprot_pos — zoom to a residue and show its 5 Å neighbourhood as sticks."""
    obj = _resolve_object(obj)
    sel = _sel_for_positions(obj, [int(uni_pos)])
    if not sel:
        print("[UFV] position %s not modelled in %s." % (uni_pos, obj))
        return
    prev = _FOCUS_SHOWN.get(obj)
    if prev:
        try:
            cmd.hide("sticks", prev)
        except Exception:
            pass
    show = "(%s) or (byres ((%s) around 5) and polymer)" % (sel, sel)
    cmd.show("sticks", show)
    _FOCUS_SHOWN[obj] = show
    cmd.zoom(show, 3)


def ufv_resetview(obj=None):
    """ufv_resetview [object] — clear focus sticks and zoom back out to the whole structure."""
    obj = _resolve_object(obj)
    prev = _FOCUS_SHOWN.pop(obj, None)
    if prev:
        try:
            cmd.hide("sticks", prev)
        except Exception:
            pass
    cmd.zoom(obj)


def ufv_report(uni_pos, obj=None, uid=None):
    """ufv_report uniprot_pos [, object [, uid]] — print the residue report + focus the residue."""
    obj = _resolve_object(obj)
    uid = _resolve_uid(uid)
    rep = residue_report(obj, uid, uni_pos)
    print(format_report(rep))
    ufv_focus(obj, uni_pos)
    return rep


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
    _SPHERE_STATE.pop(obj, None)
    _SPHERE_SHOWN.pop(obj, None)
    _CARTOON_LAYER.pop(obj, None)
    _FOCUS_SHOWN.pop(obj, None)
    if obj:
        cmd.hide("spheres", obj)
        cmd.color("gray80", obj)
    print("[UFV] cleared UFV overlays on %s." % obj)


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
        ("ufv_hotspots", ufv_hotspots), ("ufv_contacthubs", ufv_contacthubs),
        ("ufv_structures", ufv_structures), ("ufv_use", ufv_use),
        ("ufv_report", ufv_report), ("ufv_focus", ufv_focus), ("ufv_resetview", ufv_resetview),
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
        from pymol.Qt import QtWidgets, QtCore
    except Exception as e:
        print("[UFV] Qt GUI unavailable (%s). Use the ufv_* commands instead." % e)
        return
    cls = _build_panel_class(QtWidgets, QtCore)
    _gui_ref = cls(QtWidgets, QtCore)
    _gui_ref.show()
    return _gui_ref


def _build_panel_class(QtWidgets, QtCore):
    from collections import Counter

    class _Worker(QtCore.QThread):
        # Runs a pure-Python callable off the UI thread (no cmd.* calls inside!).
        done = QtCore.Signal(object)

        def __init__(self, fn):
            super(_Worker, self).__init__()
            self._fn = fn

        def run(self):
            try:
                self.done.emit(self._fn())
            except Exception as exc:  # surface errors instead of crashing the thread
                self.done.emit(exc)

    class _UFVPanel(QtWidgets.QWidget):
        def __init__(self, _qtw, _qtc):
            super(_UFVPanel, self).__init__()
            self._workers = set()
            self._wizard = None
            self.setWindowTitle("3D Feature Viewer for UniProt")
            self.setMinimumWidth(400)
            lay = QtWidgets.QVBoxLayout(self)

            # accession + fetch
            top = QtWidgets.QGridLayout(); lay.addLayout(top)
            top.addWidget(QtWidgets.QLabel("UniProt"), 0, 0)
            self.uid_edit = QtWidgets.QLineEdit(_STATE.get("uid") or ""); top.addWidget(self.uid_edit, 0, 1)
            self.fetch_btn = QtWidgets.QPushButton("Fetch"); top.addWidget(self.fetch_btn, 0, 2)
            self.fetch_btn.clicked.connect(self.on_fetch)

            # structure
            sbox = QtWidgets.QGroupBox("Structure"); sl = QtWidgets.QGridLayout(sbox); lay.addWidget(sbox)
            self.struct_combo = QtWidgets.QComboBox(); sl.addWidget(self.struct_combo, 0, 0, 1, 2)
            self.load_btn = QtWidgets.QPushButton("Load selected"); sl.addWidget(self.load_btn, 1, 0)
            self.load_btn.clicked.connect(self.on_load_selected)
            self.obj_label = QtWidgets.QLabel("Object: -"); sl.addWidget(self.obj_label, 1, 1)

            # layers
            lbox = QtWidgets.QGroupBox("Layers"); gl = QtWidgets.QGridLayout(lbox); lay.addWidget(lbox)
            self.cb_ptm = QtWidgets.QCheckBox("PTMs"); self.cb_site = QtWidgets.QCheckBox("Functional sites")
            gl.addWidget(self.cb_ptm, 0, 0); gl.addWidget(self.cb_site, 0, 1)
            self.cb_ptm.clicked.connect(lambda: self.toggle_points("ptm", self.cb_ptm, ufv_ptms))
            self.cb_site.clicked.connect(lambda: self.toggle_points("site", self.cb_site, ufv_sites))
            self.cb_var = QtWidgets.QCheckBox("Disease variants"); gl.addWidget(self.cb_var, 1, 0)
            self.cb_reviewed = QtWidgets.QCheckBox("reviewed only"); gl.addWidget(self.cb_reviewed, 1, 1)
            self.cons_cb = {}
            for i, (tok, lbl) in enumerate([("pathogenic", "Pathogenic"), ("deleterious", "Pred. deleterious"),
                                            ("benign", "Benign"), ("uncertain", "Uncertain")]):
                c = QtWidgets.QCheckBox(lbl); c.setChecked(tok in ("pathogenic", "deleterious"))
                gl.addWidget(c, 2 + i // 2, i % 2); self.cons_cb[tok] = c
                c.clicked.connect(self.refresh_variants)
            self.cb_var.clicked.connect(self.refresh_variants)
            self.cb_reviewed.clicked.connect(self.refresh_variants)

            # cartoon colouring
            cbox = QtWidgets.QGroupBox("Cartoon colouring"); cl = QtWidgets.QHBoxLayout(cbox); lay.addWidget(cbox)
            self.cart = QtWidgets.QComboBox()
            self.cart.addItems(["None", "Domains", "Topology", "pLDDT", "B-factor", "AlphaMissense",
                                "Burden", "Hotspots", "Contact hubs"])
            cl.addWidget(self.cart)
            self.cart.currentIndexChanged.connect(lambda _i: self.apply_cartoon())

            # annotation list + detail (the report generator)
            rbox = QtWidgets.QGroupBox("Annotations / residue report"); rl = QtWidgets.QVBoxLayout(rbox); lay.addWidget(rbox)
            row = QtWidgets.QHBoxLayout(); rl.addLayout(row)
            self.list_kind = QtWidgets.QComboBox(); self.list_kind.addItems(["PTMs", "Variants", "Sites", "Domains"])
            self.list_kind.currentIndexChanged.connect(lambda _i: self.refresh_table())
            row.addWidget(QtWidgets.QLabel("List:")); row.addWidget(self.list_kind)
            self.cb_pick = QtWidgets.QCheckBox("Pick in 3D"); row.addWidget(self.cb_pick)
            self.cb_pick.clicked.connect(self.toggle_pick)
            row.addStretch(1)
            self.reset_btn = QtWidgets.QPushButton("Reset view"); row.addWidget(self.reset_btn)
            self.reset_btn.clicked.connect(lambda: (ufv_resetview(self.obj())))
            self.table = QtWidgets.QTableWidget(0, 3)
            self.table.setHorizontalHeaderLabels(["Pos", "What", "Detail"])
            self.table.horizontalHeader().setStretchLastSection(True)
            self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.table.setMaximumHeight(150)
            self.table.cellClicked.connect(self.on_row)
            rl.addWidget(self.table)
            self.detail = QtWidgets.QPlainTextEdit(); self.detail.setReadOnly(True); self.detail.setMaximumHeight(150)
            rl.addWidget(self.detail)

            # advanced numbering (collapsible)
            self.adv_btn = QtWidgets.QToolButton(); self.adv_btn.setText("Advanced numbering (loaded / trajectory)")
            self.adv_btn.setCheckable(True); self.adv_btn.setStyleSheet("QToolButton{border:none;}")
            self.adv_btn.setArrowType(QtCore.Qt.RightArrow)
            lay.addWidget(self.adv_btn)
            self.adv = QtWidgets.QWidget(); al = QtWidgets.QGridLayout(self.adv); self.adv.setVisible(False)
            lay.addWidget(self.adv)
            self.rb_id = QtWidgets.QRadioButton("Identity"); self.rb_id.setChecked(True)
            self.rb_si = QtWidgets.QRadioButton("SIFTS, PDB:")
            self.rb_man = QtWidgets.QRadioButton("Manual chain:")
            self.pdb_edit = QtWidgets.QLineEdit(); self.pdb_edit.setMaximumWidth(70)
            al.addWidget(self.rb_id, 0, 0, 1, 4)
            al.addWidget(self.rb_si, 1, 0); al.addWidget(self.pdb_edit, 1, 1)
            al.addWidget(self.rb_man, 2, 0)
            self.ch_e = QtWidgets.QLineEdit("A"); self.ch_e.setMaximumWidth(34)
            self.us_e = QtWidgets.QLineEdit(); self.us_e.setPlaceholderText("UniProt@"); self.us_e.setMaximumWidth(72)
            self.rs_e = QtWidgets.QLineEdit(); self.rs_e.setPlaceholderText("resi@"); self.rs_e.setMaximumWidth(58)
            self.re_e = QtWidgets.QLineEdit(); self.re_e.setPlaceholderText("end"); self.re_e.setMaximumWidth(52)
            mh = QtWidgets.QHBoxLayout()
            for x in (self.ch_e, self.us_e, self.rs_e, self.re_e):
                mh.addWidget(x)
            al.addLayout(mh, 2, 1, 1, 3)
            apply_btn = QtWidgets.QPushButton("Apply numbering"); al.addWidget(apply_btn, 3, 0, 1, 4)
            apply_btn.clicked.connect(self.apply_map)

            def _toggle_adv():
                vis = self.adv_btn.isChecked()
                self.adv.setVisible(vis)
                self.adv_btn.setArrowType(QtCore.Qt.DownArrow if vis else QtCore.Qt.RightArrow)
            self.adv_btn.clicked.connect(_toggle_adv)

            self.status = QtWidgets.QLabel(""); self.status.setStyleSheet("color:#a33; font-size:11px;")
            lay.addWidget(self.status)
            self.info = QtWidgets.QLabel(""); self.info.setWordWrap(True); self.info.setStyleSheet("font-size:11px;")
            lay.addWidget(self.info)
            self.clear_btn = QtWidgets.QPushButton("Clear all overlays"); lay.addWidget(self.clear_btn)
            self.clear_btn.clicked.connect(self.on_clear)

            if self.cur_uid() and self.cur_uid() in _CACHE:
                self.after_fetch()

        # ---- helpers ----
        def cur_uid(self):
            return self.uid_edit.text().strip().upper() or _STATE.get("uid")

        def obj(self):
            return _STATE.get("obj") or _resolve_object(None)

        def _async(self, work, apply, status):
            self.status.setText(status)
            self.setEnabled(False)
            wk = _Worker(work)

            def on_done(res):
                self.setEnabled(True)
                self.status.setText("")
                self._workers.discard(wk)
                if isinstance(res, Exception):
                    self.status.setText("Error: %s" % res)
                    return
                if apply:
                    apply(res)
            wk.done.connect(on_done)
            self._workers.add(wk)
            wk.start()

        def update_info(self):
            ann = _CACHE.get(self.cur_uid() or "")
            if not ann:
                self.info.setText("Not fetched."); return
            cons = Counter(_cons_token(v["consequence"]) for v in ann["variants"])
            reviewed = sum(1 for v in ann["variants"] if v.get("reviewed"))
            self.info.setText(
                "%s - %d aa | PTMs %d  Sites %d  Domains %d  Topology %d  Burden %d\n"
                "Variants %d (reviewed %d; path %d, delet %d, benign %d, uncertain %d)  AlphaMissense %s"
                % (ann["uid"], len(ann["sequence"]), len(ann["ptms"]), len(ann["sites"]),
                   len(ann["domains"]), len(ann["topology"]), len(ann["burden"]), len(ann["variants"]),
                   reviewed, cons.get("pathogenic", 0), cons.get("deleterious", 0),
                   cons.get("benign", 0), cons.get("uncertain", 0), "yes" if ann["amMean"] else "no"))

        def obj_label_refresh(self):
            self.obj_label.setText("Object: %s" % (self.obj() or "-"))

        # ---- fetch / structures ----
        def on_fetch(self):
            uid = self.cur_uid()
            if not uid:
                self.status.setText("Enter an accession."); return
            self._async(lambda: (fetch_annotations(uid), list_structures(uid)),
                        lambda res: self.after_fetch(res[1]), "Fetching %s ..." % uid)

        def after_fetch(self, structs=None):
            if structs is None:
                structs = _STATE.get("structures") or []
            else:
                _STATE["structures"] = structs
            self.struct_combo.clear()
            for s in structs:
                self.struct_combo.addItem(s["label"], s)
            self.update_info()
            self.refresh_table()
            self.obj_label_refresh()

        def on_load_selected(self):
            s = self.struct_combo.currentData()
            if not s:
                self.status.setText("Press Fetch, then choose a structure."); return
            uid = self.cur_uid()
            # Download AND fetch SIFTS off the UI thread (both are network) so load can't hang.
            self._async(lambda: _prepare_structure(s, uid),
                        lambda prep: self._finish_load(s, prep, uid),
                        "Downloading %s ..." % s["label"])

        def _finish_load(self, struct, prep, uid):
            if not prep or not prep.get("path"):
                self.status.setText("Download failed."); return
            _load_from_path(struct, prep["path"], uid, prep.get("segments"))
            self.obj_label_refresh()

        # ---- layers ----
        def toggle_points(self, tag, cb, fn):
            o = self.obj()
            if not o or not _CACHE.get(self.cur_uid()):
                cb.setChecked(False)
                self.status.setText("Fetch an accession and load a structure first.")
                return
            if cb.isChecked():
                fn(o, self.cur_uid())
            else:
                _hide_layer(o, tag)

        def refresh_variants(self):
            o = self.obj()
            if not o:
                self.cb_var.setChecked(False); return
            if not self.cb_var.isChecked():
                _hide_layer(o, "var"); return
            ann = _CACHE.get(self.cur_uid())
            if not ann:
                self.cb_var.setChecked(False); return
            toks = {t for t, c in self.cons_cb.items() if c.isChecked()}
            _set_sphere_layer(o, "var", _variant_groups(ann, toks or None, self.cb_reviewed.isChecked()))

        def apply_cartoon(self):
            o = self.obj(); choice = self.cart.currentText()
            if not o:
                return
            if choice == "None":
                _reset_cartoon(o); return
            uid = self.cur_uid()
            # pLDDT / B-factor read the structure only; the rest need cached annotations.
            if choice not in ("pLDDT", "B-factor") and not _CACHE.get(uid):
                self.status.setText("Fetch an accession first.")
                self.cart.setCurrentIndex(0); return
            cheap = {"Domains": ufv_domains, "Topology": ufv_topology, "pLDDT": ufv_plddt,
                     "B-factor": ufv_bfactor, "AlphaMissense": ufv_alphamissense, "Burden": ufv_burden}
            if choice in cheap:
                cheap[choice](o, uid); return
            ann = _CACHE.get(uid) or {}
            geom = _ca_geometry(o)  # main thread (reads structure)
            if choice == "Hotspots":
                work = lambda: _per_chain_merge(geom, lambda gs: _compute_hotspots(gs, ann.get("variants", [])))
                colors = {"strong": "#b71c1c", "moderate": "#e64a19", "weak": "#ffa726"}
            else:  # Contact hubs
                work = lambda: _per_chain_merge(geom, _betweenness_hubs)
                colors = {"strong": "#6a1b9a", "moderate": "#ab47bc"}

            def apply(tiers):
                _set_cartoon_layer(o, choice.lower())
                _color_tiers(o, tiers, colors)
                self.status.setText("%s: %d residues." % (choice, len(tiers)))
            self._async(work, apply, "Computing %s ..." % choice)

        # ---- list table + report ----
        def refresh_table(self):
            ann = _CACHE.get(self.cur_uid())
            self.table.setRowCount(0)
            if not ann:
                return
            kind = self.list_kind.currentText()
            rows = []
            if kind == "PTMs":
                for p in ann["ptms"]:
                    rows.append((p["position"], p["category"], p["description"]))
            elif kind == "Variants":
                for v in ann["variants"]:
                    rows.append((v["position"], "%s%d%s" % (v.get("wildType", ""), v["position"], v.get("mutant", "")),
                                 v["consequence"] + (" | " + "; ".join(v["diseases"]) if v.get("diseases") else "")))
            elif kind == "Sites":
                for s in ann["sites"]:
                    rows.append((s["position"], "site", s["description"]))
            else:
                for d in ann["domains"]:
                    rng = "%d-%d" % (d["position"], d["endPosition"]) if d["isRange"] else str(d["position"])
                    rows.append((d["position"], rng, d["description"]))
            self.table.setRowCount(len(rows))
            for r, (pos, what, detail) in enumerate(rows):
                self.table.setItem(r, 0, QtWidgets.QTableWidgetItem(str(pos)))
                self.table.setItem(r, 1, QtWidgets.QTableWidgetItem(str(what)))
                self.table.setItem(r, 2, QtWidgets.QTableWidgetItem(str(detail)))

        def on_row(self, r, _c):
            item = self.table.item(r, 0)
            if not item:
                return
            try:
                pos = int(item.text())
            except ValueError:
                return
            self.report_residue(pos)

        def report_residue(self, uni_pos):
            o = self.obj(); uid = self.cur_uid()
            if not o or not uid:
                return
            rep = residue_report(o, uid, uni_pos)
            self.detail.setPlainText(format_report(rep))
            ufv_focus(o, uni_pos)

        # ---- 3D picking ----
        def toggle_pick(self):
            if self.cb_pick.isChecked():
                self._install_wizard()
            else:
                self._remove_wizard()

        def _install_wizard(self):
            try:
                from pymol.wizard import Wizard
            except Exception as e:
                self.status.setText("3D pick unavailable (%s)" % e); self.cb_pick.setChecked(False); return
            panel = self

            class _PickW(Wizard):
                def do_select(self, name):
                    info = []
                    try:
                        cmd.iterate(name, "info.append((chain, resi))", space={"info": info})
                    except Exception:
                        pass
                    cmd.delete(name)
                    self._report(info)
                    return None

                def do_pick(self, bondFlag):
                    info = []
                    try:
                        cmd.iterate("pk1", "info.append((chain, resi))", space={"info": info})
                    except Exception:
                        pass
                    cmd.unpick()
                    self._report(info)
                    return None

                def _report(self, info):
                    if not info:
                        return
                    ch, resi = info[0]
                    uni = _resi_to_uni(panel.obj(), ch, resi)
                    if uni is not None:
                        panel.report_residue(uni)

                def get_prompt(self):
                    return ["Click an atom to report its residue (UFV)."]
            self._wizard = _PickW()
            cmd.set_wizard(self._wizard)

        def _remove_wizard(self):
            try:
                cmd.unset_wizard()
            except Exception:
                try:
                    cmd.set_wizard()
                except Exception:
                    pass
            self._wizard = None

        # ---- numbering / clear ----
        def apply_map(self):
            o, uid = self.obj(), self.cur_uid()
            if not o or not uid:
                self.status.setText("Set accession and load/select an object first."); return
            if self.rb_si.isChecked():
                ufv_map(o, uid, "sifts", self.pdb_edit.text().strip())
            elif self.rb_man.isChecked():
                ufv_chain(o, self.ch_e.text().strip(), self.us_e.text().strip(),
                          self.rs_e.text().strip() or None, self.re_e.text().strip() or None)
            else:
                ufv_map(o, uid, "identity")
            _GEOM_CACHE.pop(o, None)

        def on_clear(self):
            ufv_clear(self.obj())
            for c in (self.cb_ptm, self.cb_site, self.cb_var):
                c.setChecked(False)
            self.cart.setCurrentIndex(0)
            self.detail.setPlainText("")

        def closeEvent(self, ev):
            self._remove_wizard()
            super(_UFVPanel, self).closeEvent(ev)

    return _UFVPanel


# Network prep runs off the UI thread (no cmd.* here): download the coordinates AND, for a PDB,
# fetch the SIFTS segments — so neither stalls the panel. _load_from_path then only does cmd ops.
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
        _set_sifts_map(name, struct["pdbId"], uid)  # fallback (will fetch synchronously)
    else:
        _set_identity_map(name)
    _GEOM_CACHE.pop(name, None)
    _STATE["obj"] = name
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
          "ufv_hotspots, ufv_contacthubs, ufv_report, ufv_focus, ufv_resetview, ufv_hide, ufv_clear, ufv_info")


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
