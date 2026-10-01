"""ph_compare.py — clique complex vs Dowker (actual co-occurrence) complex:
persistent homology (H0, H1) on the *same* 1-skeleton filtration.

Both complexes share vertices and edges (the one-mode projection with pair
co-occurrence counts).  They differ only in which triangles exist and when:

  clique complex  : triangle {a,b,c} appears as soon as its three edges exist
                    (this is what cograph::build_simplicial / any weighted
                    clique-complex filtration does)
  Dowker complex  : triangle {a,b,c} appears only when some hyperedge
                    (simplex in the ScHoLP data) actually contains all three

Filtration modes
  count : sublevel of f = -count (strong co-occurrences enter first)
          edge f = -#hyperedges containing the pair
          Dowker triangle f = -#hyperedges containing the triple
          clique triangle f = max over its three edges
  time  : sublevel of first-appearance time
          edge f = first time the pair co-occur
          Dowker triangle f = first time the triple co-occur
          clique triangle f = max over its three edges  (cumulative network)

Since #triple <= min #pair and t_first(triple) >= max t_first(pair), the
Dowker triangle never appears before the clique triangle, and K_D(t) ⊆ K_C(t).
Hence beta_1^C(t) <= beta_1^D(t): the difference is the number of H1 classes
that the clique complex kills with triangles that never co-occurred (phantom
fillings).

Pure Python/numpy implementation (GF(2) column reduction with bigint bitsets),
because gudhi/ripser are unavailable in this environment.  Dimension cap 2, so
H0 and H1 are exact; H2 is not computed.
"""
import argparse, gzip, os, sys, time, json
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(root, name):
    d = os.path.join(root, name, name)
    rd = lambda s: np.loadtxt(gzip.open(f"{d}-{s}.txt.gz"), dtype=np.int64)
    nv, sx, tm = rd("nverts"), rd("simplices"), rd("times")
    starts = np.concatenate([[0], np.cumsum(nv)[:-1]])
    return nv, sx, tm, starts


def build(nv, sx, tm, starts, maxk=25, mode="count", split=1.0):
    """Return edge dict {(u,v): f}, dowker triangle dict {(a,b,c): f}."""
    order = np.argsort(tm, kind="stable")
    cut = int(split * len(order))
    idx = order[:cut]
    ew, tw = {}, {}
    for i in idx:
        k = nv[i]
        if k < 2 or k > maxk:
            continue
        s = sorted(set(sx[starts[i]:starts[i] + k].tolist()))
        if len(s) < 2:
            continue
        t = float(tm[i])
        for e in combinations(s, 2):
            if mode == "count":
                ew[e] = ew.get(e, 0) - 1
            else:
                if e not in ew or t < ew[e]:
                    ew[e] = t
        if len(s) >= 3:
            for tr in combinations(s, 3):
                if mode == "count":
                    tw[tr] = tw.get(tr, 0) - 1
                else:
                    if tr not in tw or t < tw[tr]:
                        tw[tr] = t
    return ew, tw


def clique_triangles(ew):
    adj = defaultdict(set)
    for u, v in ew:
        adj[u].add(v); adj[v].add(u)
    out = {}
    for (u, v), f in ew.items():
        for w in adj[u] & adj[v]:
            if w > v:
                out[(u, v, w)] = max(f, ew[(u, w)], ew[(v, w)])
    return out


class UF:
    def __init__(self):
        self.p = {}
    def find(self, x):
        p = self.p
        p.setdefault(x, x)
        while p[x] != x:
            p[x] = p[p[x]]; x = p[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        self.p[ra] = rb
        return True


def persistence(ew, tri, vertices):
    """H0 and H1 bars.  ew: {edge: f}, tri: {triangle: f}.
    Returns dict with h0 (finite deaths list), h1 bars [(birth, death|inf)],
    plus the sorted edge list and positive-edge set for cross-complex matching."""
    edges = sorted(ew, key=lambda e: (ew[e], e))
    eid = {e: i for i, e in enumerate(edges)}
    # H0 via union-find; positive edges create H1 classes
    uf = UF()
    h0_deaths, positive = [], []
    for e in edges:
        if uf.union(*e):
            h0_deaths.append(ew[e])
        else:
            positive.append(e)
    # H1: reduce boundary_2 columns in filtration order
    tris = sorted(tri, key=lambda t: (tri[t], t))
    pivots = {}          # pivot edge id -> reduced column bitset
    death = {}           # birth edge -> death value
    for t in tris:
        a, b, c = t
        col = (1 << eid[(a, b)]) | (1 << eid[(a, c)]) | (1 << eid[(b, c)])
        while col:
            p = col.bit_length() - 1
            q = pivots.get(p)
            if q is None:
                pivots[p] = col
                death[edges[p]] = tri[t]
                break
            col ^= q
    h1 = [(ew[e], death.get(e, np.inf)) for e in positive]
    ncomp = len(vertices) - len(h0_deaths)
    return dict(h0_deaths=h0_deaths, h1=h1, n_components=ncomp,
                death=death, positive=positive)


def betti1_curve(h1, grid):
    b = np.zeros(len(grid), dtype=np.int64)
    for birth, d in h1:
        b += (grid >= birth) & (grid < d)
    return b


def run(root, name, mode, maxk, split, outdir):
    t0 = time.time()
    nv, sx, tm, starts = load(root, name)
    ew, tw = build(nv, sx, tm, starts, maxk=maxk, mode=mode, split=split)
    verts = set()
    for u, v in ew:
        verts.add(u); verts.add(v)
    tc = clique_triangles(ew)
    print(f"[{name}/{mode}] V={len(verts)} E={len(ew)} tri_clique={len(tc)} "
          f"tri_dowker={len(tw)} phantom_tri={len(tc)-len(tw)} "
          f"({(len(tc)-len(tw))/max(len(tc),1):.1%})  {time.time()-t0:.0f}s", flush=True)
    PC = persistence(ew, tc, verts)
    PD = persistence(ew, tw, verts)
    print(f"  reduced  {time.time()-t0:.0f}s", flush=True)

    # per-birth-edge comparison (same positive edges in both complexes)
    rows = []
    for e in PC["positive"]:
        dc = PC["death"].get(e, np.inf)
        dd = PD["death"].get(e, np.inf)
        rows.append((ew[e], dc, dd))
    rows = np.array(rows, dtype=float)
    births, dC, dD = rows[:, 0], rows[:, 1], rows[:, 2]
    finC, finD = np.isfinite(dC), np.isfinite(dD)
    # classes: (i) killed in both, (ii) killed only by clique (phantom-filled hole),
    # (iii) essential in both
    both = finC & finD
    phantom_only = finC & ~finD
    essential = ~finC & ~finD
    # among (i): clique death earlier than dowker death (premature filling)
    premature = both & (dC < dD)

    grid = np.unique(np.concatenate([births, dC[finC], dD[finD]]))
    bC = betti1_curve(PC["h1"], grid)
    bD = betti1_curve(PD["h1"], grid)

    summary = dict(
        dataset=name, mode=mode, maxk=maxk, split=split,
        V=len(verts), E=len(ew), tri_clique=len(tc), tri_dowker=len(tw),
        phantom_triangles=len(tc) - len(tw),
        phantom_triangle_frac=(len(tc) - len(tw)) / max(len(tc), 1),
        h0_components=PC["n_components"],
        h1_births=len(births),
        h1_finite_clique=int(finC.sum()), h1_finite_dowker=int(finD.sum()),
        h1_essential_clique=int((~finC).sum()), h1_essential_dowker=int((~finD).sum()),
        h1_killed_both=int(both.sum()),
        h1_phantom_only=int(phantom_only.sum()),        # true hole, filled by phantom
        h1_premature=int(premature.sum()),              # filled, but earlier in clique
        h1_essential_both=int(essential.sum()),
        betti1_max_clique=int(bC.max()), betti1_max_dowker=int(bD.max()),
        betti1_gap_max=int((bD - bC).max()),
        betti1_gap_L1=float(np.abs(bD - bC).sum()),
        mean_life_clique=float(np.mean(dC[finC] - births[finC])) if finC.any() else None,
        mean_life_dowker=float(np.mean(dD[finD] - births[finD])) if finD.any() else None,
        seconds=time.time() - t0,
    )
    os.makedirs(outdir, exist_ok=True)
    stem = os.path.join(outdir, f"{name}_{mode}")
    np.savetxt(stem + "_bars.csv", rows, delimiter=",",
               header="birth,death_clique,death_dowker", comments="")
    np.savetxt(stem + "_betti1.csv", np.column_stack([grid, bC, bD]), delimiter=",",
               header="f,betti1_clique,betti1_dowker", comments="")
    with open(stem + "_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    print("  " + json.dumps({k: summary[k] for k in
          ["phantom_triangle_frac", "h1_births", "h1_finite_clique", "h1_finite_dowker",
           "h1_phantom_only", "h1_premature", "h1_essential_both",
           "betti1_max_clique", "betti1_max_dowker", "betti1_gap_max"]}), flush=True)
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("datasets", nargs="+")
    ap.add_argument("--root", default=os.path.join(HERE, "data"))
    ap.add_argument("--mode", default="count", choices=["count", "time"])
    ap.add_argument("--maxk", type=int, default=25)
    ap.add_argument("--split", type=float, default=1.0)
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    a = ap.parse_args()
    for d in a.datasets:
        run(a.root, d, a.mode, a.maxk, a.split, a.out)
