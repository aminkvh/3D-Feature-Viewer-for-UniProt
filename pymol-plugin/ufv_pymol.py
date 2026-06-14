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


# Per-object, per-layer registry of the named selections a layer created, so each layer can be
# toggled off independently without disturbing the others.
_LAYER_SELS = {}      # obj -> { tag -> [selection names] }
_CARTOON_LAYER = {}   # obj -> name of the active cartoon-colouring layer (domains/topology/am)

# Point (sphere) layers are independent and can coexist; cartoon-colour layers are mutually
# exclusive (only one cartoon colouring at a time), matching the extension's single colour mode.
POINT_TAGS = {"ptms": "ptm", "variants": "var", "sites": "site"}
CARTOON_TAGS = ("domains", "topology", "alphamissense")


def _cons_token(cons):
    c = (cons or "").lower()
    if "pathogenic" in c:
        return "pathogenic"
    if "benign" in c:
        return "benign"
    if "deleterious" in c:
        return "deleterious"
    return "uncertain"


def _hide_layer(obj, tag):
    """Remove the named selections of one point layer and hide their spheres."""
    for name in _LAYER_SELS.get(obj, {}).pop(tag, []):
        try:
            cmd.hide("spheres", name)
            cmd.delete(name)
        except Exception:
            pass


def _show_point_groups(obj, groups, tag):
    """groups: dict color_hex -> [positions]; (re)draws this layer as coloured Ca spheres."""
    _hide_layer(obj, tag)
    cmd.set("sphere_scale", 0.5, obj)
    names, total = [], 0
    for color, positions in groups.items():
        sel = _sel_for_positions(obj, positions, ca_only=True)
        if not sel:
            continue
        name = "ufv_%s_%s_%s" % (obj, tag, color.lstrip("#"))
        cmd.select(name, sel)
        cmd.show("spheres", name)
        cmd.color(_color_name(color), name)
        total += cmd.count_atoms(name)
        names.append(name)
        cmd.deselect()
    _LAYER_SELS.setdefault(obj, {})[tag] = names
    return total


def _reset_cartoon(obj):
    """Drop any cartoon-colour layer: recolour neutral grey and remove its single-residue spheres."""
    cmd.color("gray80", obj)
    _hide_layer(obj, "dom_pts")
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
    n = _show_point_groups(obj, groups, "ptm")
    print("[UFV] %s: %d PTM spheres." % (obj, n))


def _variant_groups(ann, tokens):
    groups = {}
    for v in ann["variants"]:
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
    n = _show_point_groups(obj, _variant_groups(ann, tokens), "var")
    print("[UFV] %s: %d variant spheres (%s)." % (obj, n, only or "all"))


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
    print("[UFV] %s: %d functional-site spheres." % (obj, n))


def ufv_domains(obj=None, uid=None):
    """ufv_domains [object [, uniprot_id]] — colour the cartoon by domain / region / repeat."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    _set_cartoon_layer(obj, "domains")
    names, n = [], 0
    for i, d in enumerate(ann["domains"]):
        positions = list(range(d["position"], d["endPosition"] + 1))
        if d["isRange"]:
            sel = _sel_for_positions(obj, positions)
            if not sel:
                continue
            cmd.color(_color_name(d["color"]), sel)
        else:
            sel = _sel_for_positions(obj, positions, ca_only=True)
            if not sel:
                continue
            name = "ufv_%s_dompt_%d" % (obj, i)
            cmd.select(name, sel)
            cmd.show("spheres", name)
            cmd.color(_color_name(d["color"]), name)
            names.append(name)
            cmd.deselect()
        n += 1
    _LAYER_SELS.setdefault(obj, {})["dom_pts"] = names
    print("[UFV] %s: coloured %d domain/region features." % (obj, n))


def ufv_topology(obj=None, uid=None):
    """ufv_topology [object [, uniprot_id]] — colour the cartoon by membrane topology."""
    obj = _resolve_object(obj)
    ann = fetch_annotations(_resolve_uid(uid))
    _set_cartoon_layer(obj, "topology")
    n = 0
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
    cmd.color(_color_name("#b9c2cf"), obj)
    buckets = {"#3d85c8": [], "#b9c2cf": [], "#e06666": [], "#b71c1c": []}
    for pos, avg in am.items():
        c = "#b71c1c" if avg >= 0.78 else "#e06666" if avg >= 0.564 else "#b9c2cf" if avg >= 0.34 else "#3d85c8"
        buckets[c].append(pos)
    for color, positions in buckets.items():
        sel = _sel_for_positions(obj, positions)
        if sel:
            cmd.color(_color_name(color), sel)
    print("[UFV] %s: coloured by AlphaMissense (%d scored positions)." % (obj, len(am)))


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
    for name in list(cmd.get_names("selections")):
        if name.startswith("ufv_"):
            cmd.delete(name)
    _LAYER_SELS.pop(obj, None)
    _CARTOON_LAYER.pop(obj, None)
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
        ("ufv_alphamissense", ufv_alphamissense), ("ufv_hide", ufv_hide),
        ("ufv_clear", ufv_clear), ("ufv_info", ufv_info),
    ]:
        cmd.extend(_name, _fn)


# ----------------------------------------------------------------------------------------------
# Qt control panel (PyMOL is Qt-based; this replaces the old Tk dialog) + Plugin Manager hook
# ----------------------------------------------------------------------------------------------
_gui_ref = None  # keep a reference so the window isn't garbage-collected


def ufv_gui():
    """ufv_gui - open the Qt control panel: pick layers, filter variants, read annotation counts."""
    global _gui_ref
    try:
        from pymol.Qt import QtWidgets
    except Exception as e:
        print("[UFV] Qt GUI unavailable (%s). Use the ufv_* commands instead." % e)
        return
    from collections import Counter

    w = QtWidgets.QWidget()
    w.setWindowTitle("3D Feature Viewer for UniProt")
    w.setMinimumWidth(370)
    lay = QtWidgets.QVBoxLayout(w)

    top = QtWidgets.QGridLayout(); lay.addLayout(top)
    top.addWidget(QtWidgets.QLabel("UniProt"), 0, 0)
    uid_edit = QtWidgets.QLineEdit(_STATE.get("uid") or ""); top.addWidget(uid_edit, 0, 1)
    fetch_btn = QtWidgets.QPushButton("Fetch"); top.addWidget(fetch_btn, 0, 2)
    load_btn = QtWidgets.QPushButton("Load AF"); top.addWidget(load_btn, 0, 3)
    top.addWidget(QtWidgets.QLabel("Object"), 1, 0)
    obj_combo = QtWidgets.QComboBox(); obj_combo.setEditable(True); top.addWidget(obj_combo, 1, 1, 1, 2)
    refresh_btn = QtWidgets.QPushButton("R"); refresh_btn.setMaximumWidth(28); top.addWidget(refresh_btn, 1, 3)

    def refresh_objs():
        cur = obj_combo.currentText()
        obj_combo.clear()
        obj_combo.addItems(cmd.get_object_list() if _HAS_PYMOL else [])
        if cur:
            obj_combo.setEditText(cur)
    refresh_objs()

    def cur_obj():
        return obj_combo.currentText() or _resolve_object(None)

    def cur_uid():
        return uid_edit.text().strip().upper() or _STATE.get("uid")

    info = QtWidgets.QLabel("Enter an accession and press Fetch.")
    info.setWordWrap(True)
    info.setStyleSheet("font-size:11px;")

    def update_info():
        ann = _CACHE.get(cur_uid() or "")
        if not ann:
            info.setText("Not fetched.")
            return
        cons = Counter(_cons_token(v["consequence"]) for v in ann["variants"])
        info.setText(
            "%s - %d aa\nPTMs %d   Sites %d   Domains %d   Topology %d\n"
            "Variants %d  (pathogenic %d, deleterious %d, benign %d, uncertain %d)\n"
            "AlphaMissense: %s"
            % (ann["uid"], len(ann["sequence"]), len(ann["ptms"]), len(ann["sites"]),
               len(ann["domains"]), len(ann["topology"]), len(ann["variants"]),
               cons.get("pathogenic", 0), cons.get("deleterious", 0), cons.get("benign", 0),
               cons.get("uncertain", 0), "yes" if ann["amMean"] else "no"))

    def need_uid():
        uid = cur_uid()
        if not uid:
            info.setText("Enter an accession and press Fetch first.")
            return False
        fetch_annotations(uid)
        update_info()
        return True

    def do_fetch():
        if need_uid():
            refresh_objs()
    fetch_btn.clicked.connect(do_fetch)
    refresh_btn.clicked.connect(refresh_objs)

    def do_load():
        uid = uid_edit.text().strip().upper()
        if not uid:
            return
        ufv_load(uid)
        refresh_objs()
        obj_combo.setEditText(uid)
        update_info()
    load_btn.clicked.connect(do_load)

    numbox = QtWidgets.QGroupBox("Residue numbering"); nl = QtWidgets.QGridLayout(numbox); lay.addWidget(numbox)
    rb_id = QtWidgets.QRadioButton("Identity (resi == UniProt)"); rb_id.setChecked(True)
    rb_si = QtWidgets.QRadioButton("SIFTS, PDB:")
    rb_man = QtWidgets.QRadioButton("Manual chain:")
    pdb_edit = QtWidgets.QLineEdit(); pdb_edit.setMaximumWidth(70)
    nl.addWidget(rb_id, 0, 0, 1, 4)
    nl.addWidget(rb_si, 1, 0); nl.addWidget(pdb_edit, 1, 1)
    nl.addWidget(rb_man, 2, 0)
    ch_e = QtWidgets.QLineEdit("A"); ch_e.setMaximumWidth(34)
    us_e = QtWidgets.QLineEdit(); us_e.setPlaceholderText("UniProt@"); us_e.setMaximumWidth(72)
    rs_e = QtWidgets.QLineEdit(); rs_e.setPlaceholderText("resi@"); rs_e.setMaximumWidth(58)
    re_e = QtWidgets.QLineEdit(); re_e.setPlaceholderText("end"); re_e.setMaximumWidth(52)
    mh = QtWidgets.QHBoxLayout()
    for x in (ch_e, us_e, rs_e, re_e):
        mh.addWidget(x)
    nl.addLayout(mh, 2, 1, 1, 3)
    apply_btn = QtWidgets.QPushButton("Apply numbering"); nl.addWidget(apply_btn, 3, 0, 1, 4)

    def apply_map():
        obj, uid = cur_obj(), cur_uid()
        if not obj or not uid:
            info.setText("Set an accession and object first.")
            return
        if rb_si.isChecked():
            ufv_map(obj, uid, "sifts", pdb_edit.text().strip())
        elif rb_man.isChecked():
            ufv_chain(obj, ch_e.text().strip(), us_e.text().strip(),
                      rs_e.text().strip() or None, re_e.text().strip() or None)
        else:
            ufv_map(obj, uid, "identity")
    apply_btn.clicked.connect(apply_map)

    lbox = QtWidgets.QGroupBox("Annotation layers"); ll = QtWidgets.QGridLayout(lbox); lay.addWidget(lbox)
    cb_ptm = QtWidgets.QCheckBox("PTMs"); cb_site = QtWidgets.QCheckBox("Functional sites")
    ll.addWidget(cb_ptm, 0, 0); ll.addWidget(cb_site, 0, 1)

    def tog_ptm():
        obj = cur_obj()
        if cb_ptm.isChecked():
            if not need_uid():
                cb_ptm.setChecked(False); return
            ufv_ptms(obj, cur_uid())
        else:
            _hide_layer(obj, "ptm")
    cb_ptm.clicked.connect(tog_ptm)

    def tog_site():
        obj = cur_obj()
        if cb_site.isChecked():
            if not need_uid():
                cb_site.setChecked(False); return
            ufv_sites(obj, cur_uid())
        else:
            _hide_layer(obj, "site")
    cb_site.clicked.connect(tog_site)

    vbox = QtWidgets.QGroupBox("Disease variants"); vl = QtWidgets.QGridLayout(vbox); lay.addWidget(vbox)
    cb_var = QtWidgets.QCheckBox("Show variants"); vl.addWidget(cb_var, 0, 0, 1, 2)
    cons_cb = {}
    for i, (tok, lbl) in enumerate([("pathogenic", "Pathogenic"), ("deleterious", "Pred. deleterious"),
                                    ("benign", "Benign"), ("uncertain", "Uncertain")]):
        c = QtWidgets.QCheckBox(lbl); c.setChecked(tok in ("pathogenic", "deleterious"))
        vl.addWidget(c, 1 + i // 2, i % 2); cons_cb[tok] = c

    def refresh_vars():
        obj = cur_obj()
        if not cb_var.isChecked():
            _hide_layer(obj, "var")
            return
        if not need_uid():
            cb_var.setChecked(False); return
        toks = {t for t, c in cons_cb.items() if c.isChecked()}
        ann = fetch_annotations(cur_uid())
        _show_point_groups(obj, _variant_groups(ann, toks or None), "var")
    cb_var.clicked.connect(refresh_vars)
    for c in cons_cb.values():
        c.clicked.connect(refresh_vars)

    cbox = QtWidgets.QGroupBox("Cartoon colouring"); cl = QtWidgets.QHBoxLayout(cbox); lay.addWidget(cbox)
    cart = QtWidgets.QComboBox(); cart.addItems(["None", "Domains", "Topology", "AlphaMissense"]); cl.addWidget(cart)

    def apply_cart():
        obj = cur_obj(); choice = cart.currentText()
        if choice == "None":
            _reset_cartoon(obj)
            return
        if not need_uid():
            cart.setCurrentIndex(0); return
        {"Domains": ufv_domains, "Topology": ufv_topology, "AlphaMissense": ufv_alphamissense}[choice](obj, cur_uid())
    cart.currentIndexChanged.connect(lambda _i: apply_cart())

    lay.addWidget(info)
    clr = QtWidgets.QPushButton("Clear all overlays"); lay.addWidget(clr)

    def do_clear():
        ufv_clear(cur_obj())
        for c in (cb_ptm, cb_site, cb_var):
            c.setChecked(False)
        cart.setCurrentIndex(0)
    clr.clicked.connect(do_clear)

    if cur_uid():
        update_info()
    w.show()
    _gui_ref = w
    return w



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
          "ufv_load, ufv_fetch, ufv_map, ufv_chain, ufv_ptms, ufv_variants, ufv_sites, ufv_domains, "
          "ufv_topology, ufv_alphamissense, ufv_hide, ufv_clear, ufv_info")


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
