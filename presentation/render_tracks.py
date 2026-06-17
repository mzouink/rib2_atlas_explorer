#!/usr/bin/env python
"""Render the sequence · chemmap · motif tracks panel for each deck example,
reusing align_figure.make_tracks (same palette as the website/decks). Run in rna env."""
import json
import os
import re
import sys
import numpy as np

ROOT = "/groups/das/home/zouinkhim/atlas_explorer"
PRES = f"{ROOT}/presentation"
AF_DIR = "/groups/das/home/zouinkhim/ribonanza_inf_aws/RNAnix/lsf/20260611_rna_motif_atlas"
sys.path.insert(0, AF_DIR)
import align_figure as AF      # noqa: E402
import viz_style as V          # noqa: E402

# match the website's default base palette (A gold · C green · G red · U blue)
V.NUC_COLORS = {"A": "#E8A317", "C": "#2E7D32", "G": "#D32F2F", "U": "#1F6FB2", "T": "#1F6FB2"}

os.makedirs(f"{PRES}/tracks", exist_ok=True)
SPANS = json.load(open(f"{ROOT}/data/motifs.json"))
A = json.load(open(f"{PRES}/analysis.json"))


def res_ints(s):
    out = []
    for tok in re.split(r"[,\s]+", s.strip()):
        m = re.match(r"[A-Za-z]*:?(\d+)(?:-(\d+))?$", tok)
        if not m:
            continue
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
    rj = json.load(open(f"{ROOT}/dist/react/{sid}.json"))
    seq = rj["seq"] or ""
    dlen = len(seq)
    if not dlen:
        print("skip (no seq):", sid); continue

    def arr(key):
        v = rj.get(key)
        a = np.full(dlen, np.nan) if not v else np.array([x if x is not None else np.nan for x in v], float)
        if len(a) < dlen:
            a = np.concatenate([a, np.full(dlen - len(a), np.nan)])
        return a[:dlen]
    sn = [(x if x is not None else float("nan")) for x in (rj.get("sn") or [float("nan")] * 2)]
    motifs = []
    for mtype, resstr in SPANS.get(sid, []):
        rs = [d for d in res_ints(resstr) if 1 <= d <= dlen]
        if rs:
            motifs.append((mtype, rs))
    d = dict(seq=seq, dms=arr("dms"), a23=arr("a23"), sn=sn, motifs=motifs, dlen=dlen)
    AF.make_tracks(sid, d, f"{PRES}/tracks/{sid}.png")
    print("tracks", sid)
print("done:", len(ids), "examples")
