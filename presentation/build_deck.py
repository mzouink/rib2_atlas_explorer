#!/usr/bin/env python
"""Analyze the filters/metrics: per-metric drop %, a good + a dropped example each,
and a metric-combination drop matrix. Writes analysis.json, matrix.png, render.pml.
Run in rna/base env."""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRES = os.path.join(ROOT, "presentation")
CFG = json.load(open(os.path.join(ROOT, "config.json")))
BASES = CFG["struct_bases"]
D = json.load(open(os.path.join(ROOT, "data", "folds.json")))
N = len(D)


def cif(sid):
    lib = sid.split("-")[1].replace("ribonanza2", "").upper()
    return os.path.join(BASES["AE" if lib in "ABCDE" else "FGH"], sid + ".cif")


# metric: (name, short, keep, val, higher_better, blurb)
METRICS = [
    ("Length ≥ 40 nt", "len", lambda f: (f["length"] or 0) >= 40, lambda f: f["length"], True,
     "Short constructs are mostly trivial hairpins; we keep design regions ≥ 40 nt."),
    ("pLDDT ≥ 90", "pLDDT", lambda f: (f["plddt"] or 0) >= 90, lambda f: f["plddt"], True,
     "Model confidence. The pool bar is 85; tightening to 90 keeps the most confident folds."),
    ("Novelty (best_tm1 ≤ 0.40)", "novel", lambda f: f["best_tm1"] is not None and f["best_tm1"] <= 0.40,
     lambda f: f["best_tm1"], False,
     "TM1 to the closest v341 PDB-RNA chain; lower = more novel. ≤ 0.40 = clearly unlike known folds."),
    ("SHAPE-supported", "SHAPE", lambda f: f["shape_ok"] == 1, lambda f: f["shape_agr"], True,
     "Chemical mapping backs the fold: tertiary residues protected, or reactivity agrees with pairing."),
    ("SHAPE agreement ≥ 0.2", "agr", lambda f: f["shape_agr"] is not None and f["shape_agr"] >= 0.2,
     lambda f: f["shape_agr"], True,
     "Correlation of 2A3 reactivity with unpaired positions; + = mapping supports the predicted fold."),
    ("Compactness ≥ 0.4", "compact", lambda f: f["contact_ratio"] is not None and f["contact_ratio"] >= 0.4,
     lambda f: f["contact_ratio"], True,
     "C1′ contact ratio (globularity proxy); higher = more compact, extended chains score near 0."),
    ("Base-paired ≥ 0.5", "paired", lambda f: f["bp_fraction"] is not None and f["bp_fraction"] >= 0.5,
     lambda f: f["bp_fraction"], True,
     "Fraction of paired positions in the predicted secondary structure; higher = more structured."),
    ("≥ 1 tertiary motif", "motif", lambda f: f["n_tert"] >= 1, lambda f: f["n_tert"], True,
     "Carries a Rosetta tertiary motif (A-minor, TL-receptor, T-loop, …) — a real 3D interaction."),
]


def pick(cands, val, want_max, named=True):
    cands = [f for f in cands if val(f) is not None and (f["name"] if named else True)]
    if not cands:
        return None
    return (max if want_max else min)(cands, key=val)


def ex(f, val):
    if not f:
        return None
    return {"id": f["id"], "name": f["name"], "letter": f["letter"], "length": f["length"],
            "plddt": f["plddt"], "value": round(val(f), 3) if val(f) is not None else None,
            "best_tm1": f["best_tm1"], "near": f["near"], "near_title": f["near_title"],
            "shape_agr": f["shape_agr"], "contact_ratio": f["contact_ratio"],
            "bp_fraction": f["bp_fraction"], "n_tert": f["n_tert"], "cif": cif(f["id"])}


out = []
for name, short, keep, val, hib, blurb in METRICS:
    keepers = [f for f in D if keep(f)]
    droppers = [f for f in D if not keep(f)]
    good = pick(keepers, val, want_max=hib) or pick(keepers, val, want_max=hib, named=False)
    bad = pick(droppers, val, want_max=not hib) or pick(droppers, val, want_max=not hib, named=False)
    out.append({"name": name, "short": short, "blurb": blurb,
                "drop": len(droppers), "drop_pct": round(100 * len(droppers) / N, 1),
                "keep": len(keepers), "good": ex(good, val), "bad": ex(bad, val)})

# combination matrix: % dropped applying both filters
labels = [m["short"] for m in out]
keeps = [m[2] for m in METRICS]
mat = [[round(100 * (1 - sum(1 for f in D if keeps[i](f) and keeps[j](f)) / N), 1)
        for j in range(len(METRICS))] for i in range(len(METRICS))]

json.dump({"n": N, "metrics": out, "matrix": mat, "labels": labels},
          open(f"{PRES}/analysis.json", "w"), indent=1)

import itertools
import csv as _csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- combination truth-table over 6 core filters (all 2^6 on/off combinations) ----
CORE = [0, 1, 2, 3, 5, 7]   # len, pLDDT, novel, SHAPE, compact, motif
core_names = [METRICS[i][1] for i in CORE]
core_keep = [METRICS[i][2] for i in CORE]
combos = []
for bits in itertools.product([0, 1], repeat=len(CORE)):
    cnt = 0
    for f in D:
        if all((not b) or core_keep[k](f) for k, b in enumerate(bits)):
            cnt += 1
    combos.append((bits, cnt, round(100 * cnt / N, 1)))
combos.sort(key=lambda r: -r[1])
with open(f"{PRES}/combinations.csv", "w", newline="") as fo:
    w = _csv.writer(fo)
    w.writerow(core_names + ["pct_left", "nb_left"])
    for bits, cnt, pct in combos:
        w.writerow(list(bits) + [pct, cnt])

# render the truth table as one figure
fig, ax = plt.subplots(figsize=(7.6, 8.2)); ax.axis("off")
col_labels = core_names + ["% left", "# left"]
cell_text, cell_colours = [], []
for bits, cnt, pct in combos:
    cell_text.append(["✓" if b else "·" for b in bits] + [f"{pct:.1f}", f"{cnt:,}"])
    g = plt.cm.YlGnBu(0.15 + 0.7 * cnt / N)
    cell_colours.append(["#dcEFdc" if b else "#ffffff" for b in bits] + [g, g])
t = ax.table(cellText=cell_text, colLabels=col_labels, cellColours=cell_colours,
             cellLoc="center", bbox=[0, 0, 1, 1])
t.auto_set_font_size(False); t.set_fontsize(6.0)
for (r, c), cell in t.get_celld().items():
    cell.set_edgecolor("#e3e8ec")
    if r == 0:
        cell.set_text_props(weight="bold", color="white"); cell.set_facecolor("#2e6f95")
fig.suptitle("All 64 on/off combinations of 6 core filters · sorted by # kept", fontsize=10, y=0.985)
fig.subplots_adjust(top=0.93, bottom=0.01, left=0.02, right=0.98)
fig.savefig(f"{PRES}/combos.png", dpi=170); plt.close(fig)

# ---- scaling figure: estimated wall-clock per pipeline step ----
steps = ["rna_motif", "SHAPE gate", "compact / paired", "USalign novelty"]
h414 = [4, 3, 3, 36]                 # hours over the 414k high-confidence set (~200 slots)
hfull = [220, 170, 170, 2000]        # hours over the full ~23M atlas (~55x, same 200 slots)
lab414 = ["~4 h", "~3 h", "~3 h", "~1–2 d"]
labfull = ["~1 wk", "~5 d", "~1 wk", "~2–3 mo"]
fig, ax = plt.subplots(figsize=(8.8, 4.7))
x = np.arange(len(steps)); w = 0.38
b1 = ax.bar(x - w / 2, h414, w, label="414k high-confidence", color="#7bb6d6")
b2 = ax.bar(x + w / 2, hfull, w, label="full atlas ~23M", color="#b8c6d0")
b1[3].set_color("#c0402c"); b2[3].set_color("#7a160c")
ax.set_yscale("log"); ax.set_ylabel("estimated wall-clock (hours, log scale)")
ax.set_xticks(x); ax.set_xticklabels(steps, fontsize=11)
for bars, labs in [(b1, lab414), (b2, labfull)]:
    for rect, l in zip(bars, labs):
        ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() * 1.06, l,
                ha="center", va="bottom", fontsize=8)
ax.legend(loc="upper left", fontsize=9)
ax.set_title("Cost per pipeline step (~200-slot LSF array) — USalign novelty is the bottleneck", fontsize=11.5)
ax.set_ylim(1, 5000)
fig.tight_layout(); fig.savefig(f"{PRES}/scaling.png", dpi=150); plt.close(fig)
print(f"combos: {len(combos)} rows -> combos.png + combinations.csv; scaling.png written")

# ---- high-value-target funnel + two routes ----
def _novel(f): return f["best_tm1"] is not None and f["best_tm1"] <= 0.40
def _shp(f): return f["shape_ok"] == 1
def _full(f): return ((f["length"] or 0) >= 100 and (f["plddt"] or 0) >= 88 and _shp(f)
                      and (f["contact_ratio"] or 0) >= 0.5 and _novel(f) and f["n_tert"] >= 1)
funnel = [("all novel folds", N),
          ("≥ 80 nt", sum((f["length"] or 0) >= 80 for f in D)),
          ("≥ 100 nt", sum((f["length"] or 0) >= 100 for f in D)),
          ("≥100 nt + full gate", sum(_full(f) for f in D))]
r1 = [f for f in D if (f["length"] or 0) >= 80 and (f["plddt"] or 0) >= 88 and _shp(f)
      and (f["contact_ratio"] or 0) >= 0.5 and _novel(f) and f["n_tert"] >= 1]
r1.sort(key=lambda f: -(f["length"] or 0))
r2 = [f for f in D if _novel(f) and _shp(f) and f["pseudoknot"] == 1 and f["n_rare"] >= 1]
r2.sort(key=lambda f: (-f["n_rare"], -(f["shape_agr"] if f["shape_agr"] is not None else -9)))


def _ej(f):
    return {k: f[k] for k in ("id", "name", "length", "plddt", "best_tm1", "shape_agr",
                              "contact_ratio", "bp_fraction", "n_tert", "n_rare", "pseudoknot", "near", "near_title")}


json.dump({"funnel": funnel, "route1_n": len(r1), "route2_n": len(r2),
           "route1": [_ej(f) for f in r1[:4]], "route2": [_ej(f) for f in r2[:4]]},
          open(f"{PRES}/hvt.json", "w"), indent=1)

fig, ax = plt.subplots(figsize=(7.4, 3.7))
flabels = [x[0] for x in funnel]; fvals = [x[1] for x in funnel]
ys = list(range(len(fvals)))[::-1]
ax.barh(ys, fvals, color=["#7bb6d6", "#7bb6d6", "#e0a44a", "#c0402c"])
ax.set_yticks(ys); ax.set_yticklabels(flabels, fontsize=11)
ax.set_xscale("log"); ax.set_xlabel("# folds (log scale)"); ax.set_xlim(0.6, 1.5e4)
for yi, v in zip(ys, fvals):
    ax.text(v * 1.25, yi, f"{v:,}", va="center", fontsize=12, weight="bold")
ax.set_title("The funnel collapses on size — 1 fold clears the full gate at ≥100 nt", fontsize=11)
fig.tight_layout(); fig.savefig(f"{PRES}/funnel.png", dpi=150); plt.close(fig)
print(f"hvt: route1={len(r1)} route2={len(r2)}; funnel.png written")

# pymol render script for the example structures
ids = []
for m in out:
    for k in ("good", "bad"):
        if m[k]:
            ids.append((m[k]["id"], m[k]["cif"]))
seen = set()
with open(f"{PRES}/render.pml", "w") as o:
    o.write("bg_color white\nset ray_opaque_background, 1\nset cartoon_ring_mode, 3\n")
    for sid, path in ids:
        if sid in seen:
            continue
        seen.add(sid)
        o.write(f"load {path}, m\nhide everything\nshow cartoon\nspectrum count, rainbow, m\n")
        o.write("orient\nset ray_shadows, 0\n")
        o.write(f"ray 700,560\npng {PRES}/thumbs/{sid}.png, dpi=120\ndelete m\n")

print(f"analysis.json: {len(out)} metrics; matrix {len(labels)}x{len(labels)}; {len(seen)} thumbs to render")
print("drops:", {m["short"]: m["drop_pct"] for m in out})
