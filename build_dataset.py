#!/usr/bin/env python
"""Build one manifest-based atlas (Ribo-1 pseudolabel / OpenKnot / RFAM-PDB) into the
explorer's dataset layout: dist/datasets/<name>/data/folds.json + structs/<key>.pdb (gz).

These manifests carry confidence + clustering only (no motifs/SHAPE/novelty-TM/bp), so the
record fills what's available and computes C1' compactness from the structure. Run in rna env.

  python build_dataset.py --name pseudolabels --manifest <.tsv> --novel <_novel.tsv> --label "Ribo-1 pseudolabel"
"""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil

ROOT = os.path.dirname(os.path.abspath(__file__))


def key_of(sid):
    k = re.sub(r"[^A-Za-z0-9_-]+", "_", sid).strip("_")[:70]
    return f"{k}_{hashlib.md5(sid.encode()).hexdigest()[:6]}"


def contact_ratio(path):
    import gemmi
    import numpy as np
    try:
        st = gemmi.read_structure(path)
    except Exception:
        return None
    pts = []
    for r in st[0][0]:
        for a in r:
            if a.name in ("C1'", "C1*"):
                pts.append([a.pos.x, a.pos.y, a.pos.z]); break
    n = len(pts)
    if n < 2:
        return None
    P = np.asarray(pts)
    D = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
    idx = np.arange(n)
    sep = np.abs(idx[:, None] - idx[None, :])
    return round(int(np.triu((D <= 8.0) & (sep >= 6), 1).sum()) / n, 4)


def fl(x):
    try:
        v = float(x); return v if v == v else None
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--novel", default="")
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", default=os.path.join(ROOT, "dist", "datasets"))
    args = ap.parse_args()
    od = os.path.join(args.out, args.name)
    os.makedirs(f"{od}/data", exist_ok=True)
    os.makedirs(f"{od}/structs", exist_ok=True)

    novel = set()
    if args.novel and os.path.exists(args.novel):
        with open(args.novel) as f:
            next(f)
            novel = {ln.split("\t", 1)[0] for ln in f if ln.strip()}

    folds, n_struct, n_cr = [], 0, 0
    with open(args.manifest) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sid = r["seq_id"]; k = key_of(sid)
            p = r.get("pdb_off_relaxed", "")
            cr = None
            if p and os.path.exists(p):
                shutil.copyfile(p, f"{od}/structs/{k}.pdb")   # gz bytes; upload with content-encoding gzip
                n_struct += 1
                cr = contact_ratio(p)
                if cr is not None:
                    n_cr += 1
            folds.append({
                "id": sid, "key": k, "name": "", "letter": "", "source": args.label,
                "sublibrary": "", "length": int(r["length"]) if r.get("length", "").isdigit() else None,
                "plddt": fl(r.get("mean_plddt_on")) or fl(r.get("mean_plddt_off")),
                "ptm": fl(r.get("mean_ptm_on")) or fl(r.get("mean_ptm_off")),
                "gpde": fl(r.get("mean_gpde_on")) or fl(r.get("mean_gpde_off")),
                "clashscore": None, "n_tert": 0, "n_rare": 0, "motifs": [], "pseudoknot": 0,
                "ss_class": "", "r2a3": None, "shape_agr": None, "mean_prot_2a3": None, "shape_ok": 0,
                "openknot": None, "overlap_ae": None,
                "is_novel_v341": 1 if sid in novel else 0,
                "best_tm1": None, "near": "", "near_title": "", "score": None,
                "contact_ratio": cr, "bp_fraction": None, "in_shortlist": 0,
                "seq_cluster_size": int(r["seq_cluster_size"]) if r.get("seq_cluster_size", "").isdigit() else None,
                "struct_rep": 1 if r.get("struct_is_representative") == "1" else 0,
            })
    folds.sort(key=lambda x: -(x["plddt"] or 0))
    json.dump(folds, open(f"{od}/data/folds.json", "w"), separators=(",", ":"))
    nov = sum(f["is_novel_v341"] for f in folds)
    print(f"{args.name}: {len(folds)} folds, {n_struct} structs, {n_cr} compactness, {nov} novel -> {od}")


if __name__ == "__main__":
    main()
