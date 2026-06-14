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

    ufv_load P35498                  # download AlphaFold model + annotate it
    ufv_gui                          # open the graphical panel

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

import json
import os
import re
import sys
import tempfile
import threading
import urllib.request

from pymol import cmd

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
    }
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


def _sel_for_positions(obj, positions, ca_only=False):
    """Build a PyMOL selection covering the given UniProt positions on the mapped object."""
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
        rsel = "resi " + "+".join(str(r) for r in resis)
        clauses.append("(chain {} and {})".format(ch, rsel) if ch else "({})".format(rsel))
    sel = "({}) and ({})".format(obj, " or ".join(clauses))
    if ca_only:
        sel += " and name CA"
    return sel


def _show_point_groups(obj, groups, tag):
    """groups: dict color_hex -> [positions]; draws coloured Ca spheres."""
    cmd.set("sphere_scale", 0.5, obj)
    total = 0
    for color, positions in groups.items():
        sel = _sel_for_positions(obj, positions, ca_only=True)
        if not sel:
            continue
        name = "ufv_%s_%s" % (tag, color.lstrip("#"))
        cmd.select(name, sel)
        cmd.show("spheres", name)
        cmd.color(_color_name(color), name)
        total += cmd.count_atoms(name)
        cmd.deselect()
    return total


# ----------------------------------------------------------------------------------------------
# Projection commands
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
    n = _show_point_groups(obj, groups, "ptm")
    print("[UFV] %s: drew %d PTM spheres." % (obj, n))


def ufv_variants(obj=None, uid=None, only=None):
    """ufv_variants [object [, uniprot_id [, only]]]
    Variant Ca spheres coloured by clinical consequence.
    only: 'pathogenic' | 'benign' | 'uncertain' | 'deleterious' to filter (optional)."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    flt = (only or "").lower()
    groups = {}
    for v in ann["variants"]:
        cons = v["consequence"].lower()
        if flt:
            if flt in ("pathogenic",) and "pathogenic" not in cons:
                continue
            if flt in ("benign",) and "benign" not in cons:
                continue
            if flt in ("uncertain",) and "uncertain" not in cons:
                continue
            if flt in ("deleterious",) and "deleterious" not in cons:
                continue
        groups.setdefault(v["consequenceColor"], []).append(v["position"])
    n = _show_point_groups(obj, groups, "var")
    print("[UFV] %s: drew %d variant spheres%s." % (obj, n, " (%s)" % only if only else ""))


def ufv_sites(obj=None, uid=None):
    """ufv_sites [object [, uniprot_id]] — active / binding / metal sites as amber spheres."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    positions = []
    for s in ann["sites"]:
        positions.append(s["position"])
        if s["endPosition"] != s["position"]:
            positions.append(s["endPosition"])
    n = _show_point_groups(obj, {SITE_COLOR: positions}, "site")
    print("[UFV] %s: drew %d functional-site spheres." % (obj, n))


def ufv_domains(obj=None, uid=None):
    """ufv_domains [object [, uniprot_id]] — colour the cartoon by domain / region / repeat."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    cmd.show("cartoon", obj)
    n = 0
    for i, d in enumerate(ann["domains"]):
        positions = list(range(d["position"], d["endPosition"] + 1))
        sel = _sel_for_positions(obj, positions)
        if not sel:
            continue
        name = "ufv_dom_%d" % i
        cmd.select(name, sel)
        if d["isRange"]:
            cmd.color(_color_name(d["color"]), name)
        else:
            cmd.show("spheres", name + " and name CA")
            cmd.color(_color_name(d["color"]), name + " and name CA")
        n += 1
        cmd.deselect()
    print("[UFV] %s: coloured %d domain/region features." % (obj, n))


def ufv_topology(obj=None, uid=None):
    """ufv_topology [object [, uniprot_id]] — colour the cartoon by membrane topology."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    cmd.show("cartoon", obj)
    n = 0
    for i, t in enumerate(ann["topology"]):
        sel = _sel_for_positions(obj, list(range(t["start"], t["end"] + 1)))
        if not sel:
            continue
        name = "ufv_topo_%d" % i
        cmd.select(name, sel)
        cmd.color(_color_name(t["color"]), name)
        n += 1
        cmd.deselect()
    print("[UFV] %s: coloured %d topology segments." % (obj, n))


def ufv_alphamissense(obj=None, uid=None):
    """ufv_alphamissense [object [, uniprot_id]] — colour by mean AlphaMissense pathogenicity."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    am = ann["amMean"]
    if not am:
        print("[UFV] no AlphaMissense scores available for %s." % ann["uid"])
        return
    cmd.show("cartoon", obj)
    cmd.color(_color_name("#b9c2cf"), obj)  # neutral background
    buckets = {"#3d85c8": [], "#b9c2cf": [], "#e06666": [], "#b71c1c": []}
    for pos, avg in am.items():
        c = "#b71c1c" if avg >= 0.78 else "#e06666" if avg >= 0.564 else "#b9c2cf" if avg >= 0.34 else "#3d85c8"
        buckets[c].append(pos)
    for color, positions in buckets.items():
        sel = _sel_for_positions(obj, positions)
        if sel:
            cmd.color(_color_name(color), sel)
    print("[UFV] %s: coloured by AlphaMissense (%d scored positions)." % (obj, len(am)))


def ufv_clear(obj=None):
    """ufv_clear [object] — remove UFV selections / representations and reset colour."""
    obj = _resolve_object(obj)
    for name in list(cmd.get_names("selections")):
        if name.startswith("ufv_"):
            cmd.delete(name)
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
    """ufv_load uniprot_id [, name] — download the AlphaFold model, load it, and project
    PTMs + variants + sites (identity numbering). The full annotation set is cached for the
    other ufv_* layer commands."""
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
    fetch_annotations(uid)
    ufv_ptms(name)
    ufv_variants(name)
    ufv_sites(name)
    cmd.orient(name)
    print("[UFV] loaded and annotated %s. Try: ufv_domains %s | ufv_topology %s | ufv_alphamissense %s"
          % (name, name, name, name))


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
for _name, _fn in [
    ("ufv_load", ufv_load), ("ufv_fetch", ufv_fetch), ("ufv_map", ufv_map),
    ("ufv_chain", ufv_chain), ("ufv_ptms", ufv_ptms), ("ufv_variants", ufv_variants),
    ("ufv_sites", ufv_sites), ("ufv_domains", ufv_domains), ("ufv_topology", ufv_topology),
    ("ufv_alphamissense", ufv_alphamissense), ("ufv_clear", ufv_clear), ("ufv_info", ufv_info),
]:
    cmd.extend(_name, _fn)


# ----------------------------------------------------------------------------------------------
# Optional Tk GUI + Plugin Manager hook
# ----------------------------------------------------------------------------------------------
def ufv_gui():
    """ufv_gui — open the graphical control panel."""
    try:
        import tkinter as tk
        from tkinter import ttk, messagebox
    except Exception as e:
        print("[UFV] Tk GUI unavailable (%s). Use the ufv_* commands instead." % e)
        return

    win = tk.Toplevel()
    win.title("3D Feature Viewer for UniProt")
    win.geometry("420x520")

    frm = ttk.Frame(win, padding=10)
    frm.pack(fill="both", expand=True)

    ttk.Label(frm, text="UniProt accession").grid(row=0, column=0, sticky="w")
    uid_var = tk.StringVar(value=_STATE.get("uid") or "")
    ttk.Entry(frm, textvariable=uid_var, width=18).grid(row=0, column=1, sticky="w", pady=2)

    def _obj():
        objs = cmd.get_object_list()
        return obj_var.get() or (objs[-1] if objs else None)

    ttk.Label(frm, text="Target object").grid(row=1, column=0, sticky="w")
    obj_var = tk.StringVar()
    obj_box = ttk.Combobox(frm, textvariable=obj_var, width=16,
                           values=cmd.get_object_list())
    obj_box.grid(row=1, column=1, sticky="w", pady=2)

    def refresh_objs():
        obj_box["values"] = cmd.get_object_list()
    ttk.Button(frm, text="↻", width=3, command=refresh_objs).grid(row=1, column=2, sticky="w")

    # Numbering
    ttk.Separator(frm).grid(row=2, column=0, columnspan=3, sticky="ew", pady=8)
    ttk.Label(frm, text="Residue numbering", font=("", 10, "bold")).grid(row=3, column=0, columnspan=3, sticky="w")
    mode_var = tk.StringVar(value="identity")
    ttk.Radiobutton(frm, text="Identity (resi == UniProt)", variable=mode_var, value="identity").grid(row=4, column=0, columnspan=3, sticky="w")
    sifts_row = ttk.Frame(frm); sifts_row.grid(row=5, column=0, columnspan=3, sticky="w")
    ttk.Radiobutton(sifts_row, text="SIFTS by PDB id:", variable=mode_var, value="sifts").pack(side="left")
    pdb_var = tk.StringVar()
    ttk.Entry(sifts_row, textvariable=pdb_var, width=8).pack(side="left", padx=4)

    man_row = ttk.Frame(frm); man_row.grid(row=6, column=0, columnspan=3, sticky="w")
    ttk.Radiobutton(man_row, text="Manual chain:", variable=mode_var, value="manual").pack(side="left")
    ch_var = tk.StringVar(value="A"); us_var = tk.StringVar(); rs_var = tk.StringVar(); re_var = tk.StringVar()
    for lbl, v, w in [("chain", ch_var, 3), ("UniProt@", us_var, 6), ("resi@", rs_var, 5), ("end", re_var, 5)]:
        ttk.Label(man_row, text=lbl).pack(side="left")
        ttk.Entry(man_row, textvariable=v, width=w).pack(side="left", padx=1)

    def apply_map():
        obj = _obj()
        uid = uid_var.get().strip()
        if not obj or not uid:
            messagebox.showwarning("UFV", "Set a UniProt accession and a target object.")
            return
        m = mode_var.get()
        if m == "sifts":
            ufv_map(obj, uid, "sifts", pdb_var.get().strip())
        elif m == "manual":
            ufv_chain(obj, ch_var.get().strip(), us_var.get().strip(),
                      rs_var.get().strip() or None, re_var.get().strip() or None)
        else:
            ufv_map(obj, uid, "identity")
    ttk.Button(frm, text="Apply numbering", command=apply_map).grid(row=7, column=0, columnspan=3, sticky="w", pady=4)

    # Layers
    ttk.Separator(frm).grid(row=8, column=0, columnspan=3, sticky="ew", pady=8)
    ttk.Label(frm, text="Annotation layers", font=("", 10, "bold")).grid(row=9, column=0, columnspan=3, sticky="w")

    def run(fn):
        obj = _obj(); uid = uid_var.get().strip()
        if not uid:
            messagebox.showwarning("UFV", "Enter a UniProt accession.")
            return
        threading.Thread(target=lambda: fn(obj, uid), daemon=True).start()

    layer_btns = [
        ("PTMs", ufv_ptms), ("Disease variants", ufv_variants), ("Functional sites", ufv_sites),
        ("Domains / regions", ufv_domains), ("Membrane topology", ufv_topology),
        ("AlphaMissense", ufv_alphamissense),
    ]
    for i, (label, fn) in enumerate(layer_btns):
        ttk.Button(frm, text=label, width=20, command=lambda f=fn: run(f)).grid(
            row=10 + i // 2, column=i % 2, sticky="w", pady=2, padx=2)

    ttk.Separator(frm).grid(row=14, column=0, columnspan=3, sticky="ew", pady=8)
    bottom = ttk.Frame(frm); bottom.grid(row=15, column=0, columnspan=3, sticky="w")
    ttk.Button(bottom, text="Download AlphaFold + annotate",
               command=lambda: threading.Thread(target=lambda: ufv_load(uid_var.get().strip()), daemon=True).start()).pack(side="left")
    ttk.Button(bottom, text="Clear", command=lambda: ufv_clear(_obj())).pack(side="left", padx=6)


cmd.extend("ufv_gui", ufv_gui)


def __init_plugin__(app=None):
    """PyMOL Plugin Manager entry point."""
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt("3D Feature Viewer for UniProt", ufv_gui)
    except Exception:
        pass


print("[UFV] 3D Feature Viewer for UniProt loaded. Commands: ufv_load, ufv_fetch, ufv_map, "
      "ufv_chain, ufv_ptms, ufv_variants, ufv_sites, ufv_domains, ufv_topology, "
      "ufv_alphamissense, ufv_clear, ufv_info, ufv_gui")
