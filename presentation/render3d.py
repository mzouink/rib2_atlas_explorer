"""Recolor the example 3D thumbnails like the website: cartoon by 2A3 reactivity
(blue protected -> red reactive) + tertiary-motif residues as colored sticks.
Run with: pymol -cq render3d.py"""
import json
import os
import re
from pymol import cmd

ROOT = "/groups/das/home/zouinkhim/atlas_explorer"
PRES = f"{ROOT}/presentation"
CFG = json.load(open(f"{ROOT}/config.json")); BASES = CFG["struct_bases"]
MOT = json.load(open(f"{ROOT}/data/motifs.json"))
A = json.load(open(f"{PRES}/analysis.json"))

MOTIF_COLORS = {
    "A_MINOR": "#d1495b", "TL_RECEPTOR": "#e8862e", "TETRALOOP_TL_RECEPTOR": "#e8862e",
    "UA_HANDLE": "#8c2f39", "T_LOOP": "#c879c8", "INTERCALATED_T_LOOP": "#9b5fb0",
    "GA_MINOR": "#e8a598", "PLATFORM": "#edae49", "TANDEM_GA_SHEARED": "#c9a96b",
    "TANDEM_GA_WATSON_CRICK": "#c9a96b", "U_TURN": "#2e6f95", "Z_TURN": "#5b7c99",
    "GNRA_TETRALOOP": "#16a0a0", "LOOP_E_SUBMOTIF": "#6a4c93", "BULGED_G": "#3a7d44"}
DEFAULT = "#9aa7b3"


def rgb(h):
    h = h.lstrip("#"); return [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]


def cif(sid):
    lib = sid.split("-")[1].replace("ribonanza2", "").upper()
    return os.path.join(BASES["AE" if lib in "ABCDE" else "FGH"], sid + ".cif")


def res_ints(s):
    out = []
    for tok in re.split(r"[,\s]+", s.strip()):
        m = re.match(r"[A-Za-z]*:?(\d+)(?:-(\d+))?$", tok)
        if m:
            a = int(m.group(1)); b = int(m.group(2)) if m.group(2) else a
            out += list(range(a, b + 1))
    return out


ids = []
for m in A["metrics"]:
    for k in ("good", "bad"):
        if m[k]:
            ids.append(m[k]["id"])
hp = f"{PRES}/hvt.json"
if os.path.exists(hp):
    H = json.load(open(hp))
    for r in ("route1", "route2"):
        ids += [e["id"] for e in H.get(r, [])]
ids = list(dict.fromkeys(ids))

for sid in ids:
    cmd.reinitialize()
    cmd.bg_color("white"); cmd.set("ray_opaque_background", 1); cmd.set("ray_shadows", 0)
    cmd.set("cartoon_ring_mode", 3); cmd.set("stick_radius", 0.3)
    cmd.load(cif(sid), "m")
    cmd.hide("everything"); cmd.show("cartoon")
    rj = json.load(open(f"{ROOT}/dist/react/{sid}.json"))
    a23 = rj.get("a23")
    if a23:
        bd = {str(i + 1): float(v) for i, v in enumerate(a23) if v is not None}
        cmd.alter("m", "b=bd.get(resi,0.0)", space={"bd": bd})
        cmd.spectrum("b", "blue_white_red", "m", -0.3, 1.0)
    else:
        cmd.color("grey70", "m")
    for mtype, resstr in MOT.get(sid, []):
        rs = res_ints(resstr)
        if not rs:
            continue
        cmd.set_color("mc", rgb(MOTIF_COLORS.get(mtype, DEFAULT)))
        sel = "m and resi " + "+".join(map(str, rs))
        cmd.show("sticks", sel); cmd.color("mc", sel)
    cmd.orient("m")
    cmd.ray(720, 560)
    cmd.png(f"{PRES}/thumbs/{sid}.png", dpi=120)
    print("3d", sid)
print("done")
