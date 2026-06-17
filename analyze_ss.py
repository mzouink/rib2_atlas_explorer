#!/usr/bin/env python
"""Derive secondary structure from a predicted 3D RNA structure (no external SS tool).

For each staged dataset fold we detect Watson-Crick / wobble base pairs geometrically
from the predicted coordinates, then emit a pseudoknot-aware dot-bracket plus scalar
features the explorer can filter/sort on:

    bp_fraction  = 2 * n_pairs / length      (fraction of bases that are paired)
    pseudoknot   = 1 if any base pairs cross  (assigned to a non-primary bracket level)
    dbn          = dot-bracket string (() primary, [] {} <> for crossing pairs)

Base-pair rule: >=2 hydrogen-bond-distance (<3.5 A) contacts between the two bases'
Watson-Crick edge atoms (catches WC A-U / G-C and G-U wobble). Greedy one-pair-per-base
assignment by contact count then distance. Reads the gz PDBs already staged under
dist/datasets/<name>/structs (keyed by the folds.json `key`); never touches /nrs.

    python analyze_ss.py --name openknot --workers 48
"""
import argparse
import gzip
import json
import os
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.abspath(__file__))

# Watson-Crick edge atoms per base (donors/acceptors that H-bond in a pair).
EDGE = {
    "A": ["N1", "N6"],          "G": ["N1", "N2", "O6"],
    "C": ["N3", "N4", "O2"],    "U": ["N3", "O4", "O2"],
}
MIN_LOOP = 3        # |i-j| >= 4 (at least 3 unpaired between paired bases)
HB = 3.5            # H-bond contact cutoff (A)
HB2 = HB * HB


def _edge_coords(res):
    out = {}
    nm = res.name.strip().upper()
    nm = {"RA": "A", "RG": "G", "RC": "C", "RU": "U", "ADE": "A", "GUA": "G", "CYT": "C", "URA": "U"}.get(nm, nm)
    if nm not in EDGE:
        return nm, out, None
    c1 = None
    for a in res:
        an = a.name.strip()
        if an in ("C1'", "C1*"):
            c1 = (a.pos.x, a.pos.y, a.pos.z)
        if an in EDGE[nm]:
            out[an] = (a.pos.x, a.pos.y, a.pos.z)
    return nm, out, c1


def _bp(ei, ej):
    """>=2 edge-atom contacts within H-bond distance -> base-paired."""
    n = 0
    for a in ei.values():
        for b in ej.values():
            dx = a[0] - b[0]; dy = a[1] - b[1]; dz = a[2] - b[2]
            if dx * dx + dy * dy + dz * dz <= HB2:
                n += 1
                if n >= 2:
                    return True
    return False


def _crosses(p, q):
    a, b = p; c, d = q
    return (a < c < b < d) or (c < a < d < b)


def _dbn(n, pairs):
    """Assign pairs to bracket levels (pseudoknots -> higher levels)."""
    OPEN = "([{<"; CLOSE = ")]}>"
    dbn = ["."] * n
    pk = 0
    levels = []  # list of lists of placed pairs per level
    for (i, j) in sorted(pairs):
        lvl = 0
        while True:
            if lvl >= len(levels):
                levels.append([]);
            if not any(_crosses((i, j), pq) for pq in levels[lvl]):
                levels[lvl].append((i, j))
                o = OPEN[lvl] if lvl < len(OPEN) else OPEN[-1]
                c = CLOSE[lvl] if lvl < len(CLOSE) else CLOSE[-1]
                dbn[i] = o; dbn[j] = c
                if lvl > 0:
                    pk = 1
                break
            lvl += 1
    return "".join(dbn), pk


def _process(args):
    sid, path = args
    try:
        import gemmi
        if path.endswith(".gz") or _is_gz(path):
            import tempfile
            with gzip.open(path, "rt") as fh:
                txt = fh.read()
            st = gemmi.read_pdb_string(txt)
        else:
            st = gemmi.read_structure(path)
        chain = st[0][0]
        res = [r for r in chain]
        edges, c1s = [], []
        for r in res:
            nm, e, c1 = _edge_coords(r)
            edges.append(e); c1s.append(c1)
        n = len(res)
        cand = []
        for i in range(n):
            if not edges[i] or c1s[i] is None:
                continue
            for j in range(i + MIN_LOOP + 1, n):
                if not edges[j] or c1s[j] is None:
                    continue
                # cheap C1'-C1' prefilter (paired bases sit ~8-12 A apart)
                dx = c1s[i][0] - c1s[j][0]; dy = c1s[i][1] - c1s[j][1]; dz = c1s[i][2] - c1s[j][2]
                d2 = dx * dx + dy * dy + dz * dz
                if d2 < 49 or d2 > 196:
                    continue
                if _bp(edges[i], edges[j]):
                    cand.append((d2, i, j))
        # greedy one pair per base, closest first
        cand.sort()
        used = [False] * n
        pairs = []
        for _, i, j in cand:
            if not used[i] and not used[j]:
                used[i] = used[j] = True
                pairs.append((i, j))
        dbn, pk = _dbn(n, pairs)
        bpf = round(2 * len(pairs) / n, 4) if n else 0.0
        return (sid, bpf, pk, dbn)
    except Exception as e:
        return (sid, None, None, f"ERR:{type(e).__name__}")


def _is_gz(path):
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--workers", type=int, default=48)
    args = ap.parse_args()
    dd = os.path.join(ROOT, "dist", "datasets", args.name)
    folds = json.load(open(f"{dd}/data/folds.json"))
    work = []
    for f in folds:
        p = os.path.join(dd, "structs", f["key"] + ".pdb")
        if os.path.exists(p):
            work.append((f["id"], p))
    print(f"{args.name}: {len(work)}/{len(folds)} structures")
    res = {}
    with Pool(args.workers) as pool:
        for sid, bpf, pk, dbn in pool.imap_unordered(_process, work, chunksize=16):
            res[sid] = (bpf, pk, dbn)
    nerr = sum(1 for v in res.values() if v[0] is None)
    out = f"{dd}/data/ss.json"
    json.dump(res, open(out, "w"), separators=(",", ":"))
    pkn = sum(1 for v in res.values() if v[1] == 1)
    print(f"{args.name}: wrote {out}  (errors={nerr}, pseudoknots={pkn})")


if __name__ == "__main__":
    main()
