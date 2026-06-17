#!/usr/bin/env python
"""Merge the LSF analysis outputs (motifs, novelty, secondary structure) into each
add-on dataset's explorer files, matching the ribo2 schema exactly.

Updates dist/datasets/<ds>/data/:
  folds.json   <- bp_fraction, pseudoknot, motifs, n_tert, n_rare, best_tm1, near,
                  near_title, is_novel_v341 (recomputed from best_tm1 < 0.45)
  motifs.json  <- {sid: [[motif_type, residues], ...]}   (deep-view spans / 3D sticks)
  pairing.json <- {sid: dbn}                              (deep-view pairing track)

Inputs:
  ss.json per dataset (already written by analyze_ss.py)
  rna_motif chunks    (letter column carries the dataset name; demux not needed -- ids are unique)
  novelty chunks      (best_tm1_v341, nearest_known)
"""
import csv
import glob
import json
import os
import re
import urllib.request
import urllib.parse

ROOT = os.path.dirname(os.path.abspath(__file__))
LSF = "/groups/das/home/zouinkhim/ribonanza_inf_aws/RNAnix/lsf/20260617_atlas_addon_analysis/results"
MAN = "/groups/das/home/joshic/RNAnix/release_data/distillation_atlases"
DATASETS = ["pseudolabels", "openknot", "rfam_pdb130", "rfam_pdb240"]
TERT = {"A_MINOR", "TL_RECEPTOR", "UA_HANDLE", "T_LOOP", "GA_MINOR", "PLATFORM",
        "TANDEM_GA_SHEARED", "TANDEM_GA_WATSON_CRICK", "TETRALOOP_TL_RECEPTOR"}
RARE_TERT = {"TL_RECEPTOR", "GA_MINOR", "T_LOOP", "TETRALOOP_TL_RECEPTOR", "UA_HANDLE"}
CACHE = os.path.join(ROOT, ".rcsb_titles.json")


def derive_name(ds, sid):
    """Human-readable name from the seq_id, where the id carries one.
    rfam_pdb: 'RF01807_GIR1_branching_ribozyme_URS..._5793_1-188' -> 'GIR1 branching ribozyme'.
    openknot: '13608507_Eli_Fisker_..._975_0pad25_libraryready'   -> the design description.
    pseudolabels: ids are bare hashes -> no recoverable name."""
    if ds.startswith("rfam"):
        toks = sid.split("_")
        if toks and re.fullmatch(r"RF\d+", toks[0]):
            toks = toks[1:]
        while toks and (re.fullmatch(r"URS[0-9A-Za-z]+", toks[-1]) or re.fullmatch(r"\d+", toks[-1])
                        or re.fullmatch(r"\d+-\d+", toks[-1]) or re.fullmatch(r"[A-Z0-9]{4}", toks[-1])):
            toks = toks[:-1]
        return " ".join(toks).strip()
    if ds == "openknot":
        s = re.sub(r"^\d+_", "", sid)
        s = re.sub(r"_\d+_0pad\d+_libraryready$", "", s)
        s = s.strip('"').replace("_", " ").strip()
        return s[:90]
    return ""


def load_motifs():
    spans, agg = {}, {}
    for fn in sorted(glob.glob(f"{LSF}/motifs/chunk_*.tsv")):
        with open(fn) as f:
            for m in csv.DictReader(f, delimiter="\t"):
                mt = m.get("motif_type"); sid = m.get("seq_id")
                if not mt or not sid or sid.startswith("#"):
                    continue
                spans.setdefault(sid, []).append([mt, m.get("residues", "")])
                a = agg.setdefault(sid, {"motifs": set(), "tert": set(), "rare": set()})
                a["motifs"].add(mt)
                if mt in TERT:
                    a["tert"].add(mt)
                if mt in RARE_TERT:
                    a["rare"].add(mt)
    return spans, agg


def load_novelty():
    nov = {}
    for fn in sorted(glob.glob(f"{LSF}/novelty/chunk_*.tsv")):
        with open(fn) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                try:
                    tm = float(r["best_tm1_v341"])
                except (KeyError, ValueError, TypeError):
                    continue
                near = r.get("nearest_known", "")
                if tm <= 0 or not near:   # -1 error / no USalign hit (LALI<20) -> unscored, not "novel"
                    continue
                nov[r["seq_id"]] = (round(tm, 4), near)
    return nov


def load_novel_set(ds):
    s = set()
    p = f"{MAN}/{ds}_manifest_novel.tsv"
    if os.path.exists(p):
        with open(p) as f:
            next(f, None)
            for ln in f:
                if ln.strip():
                    s.add(ln.split("\t", 1)[0])
    return s


def pdb_titles(near_ids):
    titles = {}
    if os.path.exists(CACHE):
        titles = json.load(open(CACHE))
    want = {n.split("_")[0].upper() for n in near_ids if n}
    missing = sorted(want - set(titles))
    for i in range(0, len(missing), 50):
        batch = missing[i:i + 50]
        q = '{entries(entry_ids:[%s]){rcsb_id struct{title}}}' % ",".join(f'"{b}"' for b in batch)
        try:
            url = "https://data.rcsb.org/graphql?query=" + urllib.parse.quote(q)
            with urllib.request.urlopen(url, timeout=30) as r:
                data = json.load(r)
            for e in data.get("data", {}).get("entries", []) or []:
                titles[e["rcsb_id"].upper()] = (e.get("struct") or {}).get("title", "")
        except Exception as ex:
            print(f"  title fetch failed for batch {i//50}: {ex}")
            break
    json.dump(titles, open(CACHE, "w"))
    return titles


def main():
    spans, agg = load_motifs()
    nov = load_novelty()
    print(f"motifs: {len(spans)} folds w/ motifs; novelty: {len(nov)} scored")
    titles = pdb_titles([v[1] for v in nov.values()])

    for ds in DATASETS:
        dd = f"{ROOT}/dist/datasets/{ds}"
        folds = json.load(open(f"{dd}/data/folds.json"))
        ss = json.load(open(f"{dd}/data/ss.json")) if os.path.exists(f"{dd}/data/ss.json") else {}
        novel_set = load_novel_set(ds)
        ms, pr = {}, {}
        n_ss = n_mot = n_nov = 0
        n_name = 0
        for fd in folds:
            sid = fd["id"]
            if not fd.get("name"):
                nm = derive_name(ds, sid)
                if nm:
                    fd["name"] = nm; n_name += 1
            if sid in ss and ss[sid][0] is not None:
                bpf, pk, dbn = ss[sid]
                fd["bp_fraction"] = bpf; fd["pseudoknot"] = pk
                pr[sid] = dbn; n_ss += 1
            a = agg.get(sid)
            if a:
                fd["motifs"] = sorted(a["motifs"]); fd["n_tert"] = len(a["tert"]); fd["n_rare"] = len(a["rare"])
                if sid in spans:
                    ms[sid] = spans[sid]
                n_mot += 1
            if sid in nov:
                tm, near = nov[sid]
                fd["best_tm1"] = tm; fd["near"] = near
                fd["near_title"] = titles.get(near.split("_")[0].upper(), "")
                fd["is_novel_v341"] = 1 if tm < 0.45 else 0
                n_nov += 1
            else:   # no usable USalign hit -> leave best_tm1 unscored; novelty from the manifest flag
                fd["best_tm1"] = None; fd["near"] = ""; fd["near_title"] = ""
                fd["is_novel_v341"] = 1 if sid in novel_set else 0
        json.dump(folds, open(f"{dd}/data/folds.json", "w"), separators=(",", ":"))
        json.dump(ms, open(f"{dd}/data/motifs.json", "w"), separators=(",", ":"))
        json.dump(pr, open(f"{dd}/data/pairing.json", "w"), separators=(",", ":"))
        print(f"{ds:14} n={len(folds)}  ss={n_ss} motifs={n_mot} novelty={n_nov} names={n_name}")


if __name__ == "__main__":
    main()
