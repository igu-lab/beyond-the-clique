"""null_model.py — ハイパーグラフ configuration model（次数・辺サイズ保存）による
Dowker vs クリーク比較指標のヌル分布。

ハイパーエッジのリスト（頂点集合のリスト）を入力に取り、
  観測: phantom 率、max β1^D、max β1^C、核の次元（クリークだけが潰す H1 クラス数）、β1 曲線
  ヌル: 上の指標を N サンプルについて計算 → z 値・経験的 p 値・β1 曲線の 2.5/97.5% 帯
を出す。スワップ法（2 つのハイパーエッジから 1 頂点ずつ交換、重複が出れば却下）。

使い方（モジュールとして）: from null_model import run_null; run_null(hyperedges, N=200, out=...)
CLI: python3 src/null_model.py --scholp data/processed/human_long-w3 --N 200 --out results/null_<stamp>
"""
import argparse, gzip, json, os, random, sys, time
from collections import Counter
from itertools import combinations
import numpy as np
FRAC_GRID = 101

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ph_compare import clique_triangles, persistence, betti1_curve


def edges_from_hyperedges(H):
    ew, tw = {}, {}
    for S in H:
        s = sorted(set(S))
        if len(s) < 2:
            continue
        for e in combinations(s, 2):
            ew[e] = ew.get(e, 0) - 1
        for t in combinations(s, 3):
            tw[t] = tw.get(t, 0) - 1
    return ew, tw


def analyze(H, grid=None):
    ew, tw = edges_from_hyperedges(H)
    verts = {v for e in ew for v in e}
    tc = clique_triangles(ew)
    PC = persistence(ew, tc, verts)
    PD = persistence(ew, tw, verts)
    dC = np.array([PC["death"].get(e, np.inf) for e in PC["positive"]])
    dD = np.array([PD["death"].get(e, np.inf) for e in PC["positive"]])
    births = np.array([ew[e] for e in PC["positive"]])
    if grid is None:
        grid = np.unique(np.concatenate([births, dC[np.isfinite(dC)], dD[np.isfinite(dD)]]))
    bC, bD = betti1_curve(PC["h1"], grid), betti1_curve(PD["h1"], grid)
    # 辺の割合で正規化した曲線：x = 1-骨格に入っている辺の割合（0..1）。ヌルと観測で count の分布が違っても比較できる
    fvals = np.sort(np.array(list(ew.values()), dtype=float))
    xs = np.linspace(0.0, 1.0, FRAC_GRID)
    fx = np.array([fvals[min(int(np.ceil(x * len(fvals))) - 1, len(fvals) - 1)] if x > 0 else fvals[0] - 1 for x in xs])
    bC_frac, bD_frac = betti1_curve(PC["h1"], fx), betti1_curve(PD["h1"], fx)
    return dict(xs=xs, bC_frac=bC_frac, bD_frac=bD_frac,
        phantom_frac=(len(tc) - len(tw)) / max(len(tc), 1),
        tri_clique=len(tc), tri_dowker=len(tw),
        max_b1_C=int(bC.max()) if len(bC) else 0, max_b1_D=int(bD.max()) if len(bD) else 0,
        kernel_dim=int((np.isfinite(dC) & ~np.isfinite(dD)).sum()),   # 本当の穴をクリークが埋めた数
        premature=int((np.isfinite(dC) & np.isfinite(dD) & (dC < dD)).sum()),
        gap_L1=float(np.abs(bD - bC).sum()),
        grid=grid, bC=bC, bD=bD)


def shuffle(H, n_swaps=None, rng=random):
    """次数と辺サイズを保存するスワップ。H は list of list。"""
    H = [list(S) for S in H]
    m = len(H)
    n_swaps = n_swaps or 10 * sum(len(S) for S in H)
    done = 0
    while done < n_swaps:
        i, j = rng.randrange(m), rng.randrange(m)
        if i == j:
            continue
        a, b = rng.randrange(len(H[i])), rng.randrange(len(H[j]))
        u, v = H[i][a], H[j][b]
        if u == v or v in H[i] or u in H[j]:
            continue
        H[i][a], H[j][b] = v, u
        done += 1
    return H


def run_null(H, N=200, seed=0, out=None, label="data", verbose=True, chain=False):
    """chain=False（既定, 2026-09-25 変更）: 各サンプルを元データから独立に並べ替える（サンプル間の自己相関なし）。
    chain=True: 旧来の方式（前のサンプルから続けて並べ替える）。"""
    rng = random.Random(seed)
    t0 = time.time()
    obs = analyze(H)
    grid = obs["grid"]
    keys = ["phantom_frac", "max_b1_C", "max_b1_D", "kernel_dim", "premature", "gap_L1"]
    null = {k: [] for k in keys}
    curves, curves_frac, curvesC_frac = [], [], []
    Hs = [list(S) for S in H]
    for r in range(N):
        Hs = shuffle(Hs if chain else H, rng=rng)   # chain=False: 毎回元データから独立に
        a = analyze(Hs, grid=grid)
        for k in keys:
            null[k].append(a[k])
        curves.append(a["bD"])
        curves_frac.append(a["bD_frac"]); curvesC_frac.append(a["bC_frac"])
        if verbose and (r + 1) % max(N // 10, 1) == 0:
            print(f"  [{label}] null {r+1}/{N}  {time.time()-t0:.1f}s", flush=True)
    curves = np.array(curves)
    lo, hi, med = np.percentile(curves, 2.5, axis=0), np.percentile(curves, 97.5, axis=0), np.median(curves, axis=0)
    summ = dict(label=label, N=N, n_hyperedges=len(H), n_vertices=len({v for S in H for v in S}), seconds=time.time() - t0,
                sampling="chain" if chain else "independent", p_min_attainable=1.0 / (N + 1))
    for k in keys:
        arr = np.array(null[k], dtype=float)
        mu, sd = arr.mean(), arr.std(ddof=1) if N > 1 else 0.0
        z = (obs[k] - mu) / sd if sd > 0 else float("nan")
        p_hi = (np.sum(arr >= obs[k]) + 1) / (N + 1)
        p_lo = (np.sum(arr <= obs[k]) + 1) / (N + 1)
        summ[k] = dict(obs=float(obs[k]), null_mean=float(mu), null_sd=float(sd), z=float(z), p_ge=float(p_hi), p_le=float(p_lo))
    outside = np.mean((obs["bD"] < lo) | (obs["bD"] > hi))
    summ["b1D_curve_frac_outside_95band"] = float(outside)
    summ["null_samples"] = {k: [float(v) for v in null[k]] for k in keys}
    CF = np.array(curves_frac); CCF = np.array(curvesC_frac)
    flo, fhi, fmed = np.percentile(CF, 2.5, axis=0), np.percentile(CF, 97.5, axis=0), np.median(CF, axis=0)
    clo, chi, cmed = np.percentile(CCF, 2.5, axis=0), np.percentile(CCF, 97.5, axis=0), np.median(CCF, axis=0)
    summ["b1D_fraccurve_frac_outside_95band"] = float(np.mean((obs["bD_frac"] < flo) | (obs["bD_frac"] > fhi)))
    if out:
        os.makedirs(out, exist_ok=True)
        json.dump(summ, open(os.path.join(out, f"{label}_null_summary.json"), "w"), indent=1)
        np.savetxt(os.path.join(out, f"{label}_null_b1curve.csv"),
                   np.column_stack([grid, obs["bC"], obs["bD"], lo, med, hi]), delimiter=",",
                   header="f,b1_clique_obs,b1_dowker_obs,null_lo2.5,null_median,null_hi97.5", comments="")
        np.savetxt(os.path.join(out, f"{label}_null_b1curve_edgefrac.csv"),
                   np.column_stack([obs["xs"], obs["bC_frac"], obs["bD_frac"], flo, fmed, fhi, clo, cmed, chi]), delimiter=",",
                   header="edge_frac,b1_clique_obs,b1_dowker_obs,nullD_lo2.5,nullD_median,nullD_hi97.5,nullC_lo2.5,nullC_median,nullC_hi97.5", comments="")
    return summ, obs, (grid, lo, med, hi)


def load_scholp(d):
    name = os.path.basename(d.rstrip("/"))
    base = os.path.join(d, name)
    rd = lambda s: np.loadtxt(gzip.open(f"{base}-{s}.txt.gz"), dtype=np.int64)
    nv, sx = rd("nverts"), rd("simplices")
    H, p = [], 0
    for k in nv:
        H.append(sx[p:p + k].tolist()); p += k
    return name, H


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scholp", nargs="*", default=[])
    ap.add_argument("--N", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chain", action="store_true", help="旧来の連鎖サンプリング")
    a = ap.parse_args()
    for d in a.scholp:
        name, H = load_scholp(d)
        s, _, _ = run_null(H, N=a.N, seed=a.seed, out=a.out, label=name, chain=a.chain)
        print(json.dumps({k: (v if not isinstance(v, dict) else {kk: round(vv, 3) for kk, vv in v.items()}) for k, v in s.items()}, ensure_ascii=False))
