"""check_nestimate_clique.py — 付録 A の検証：クリーク複体の再実装が Nestimate::build_simplicial と一致するか．

R 側（scripts/r/run_saqr_pipeline.R）が書き出した <name>_sc_<mat>_t<thr>_simplices.csv（build_simplicial の全単体）と，
同じ重み行列 W から Python で作ったクリーク複体（pmax(W, t(W)) >= thr を辺とし，max_dim = 3 までの全クリーク）を
単体の集合として突き合わせる．

使い方: python3 scripts/check_nestimate_clique.py --r-dir results/r_saqr
"""
import argparse, glob, itertools as it, os, re, sys
import numpy as np, pandas as pd


def clique_complex(W, thr, max_dim=3):
    n = W.shape[0]; S = np.maximum(W, W.T)
    adj = {i: {j for j in range(n) if j != i and S[i, j] >= thr} for i in range(n)}
    out = {(i,) for i in range(n)}
    for k in range(2, max_dim + 2):
        for c in it.combinations(range(n), k):
            if all(b in adj[a] for a, b in it.combinations(c, 2)):
                out.add(c)
    return {tuple(i + 1 for i in c) for c in out}  # 1-based (R の列番号)


def r_simplices(path):
    df = pd.read_csv(path, dtype={"nodes": str})
    return {tuple(sorted(int(x) for x in s.split("|"))) for s in df["nodes"]}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--r-dir", required=True); a = ap.parse_args()
    pat = re.compile(r"(?P<name>.+)_sc_(?P<mat>.+)_t(?P<thr>[0-9.]+)_simplices\.csv$")
    rows = []
    for p in sorted(glob.glob(os.path.join(a.r_dir, "*_sc_*_simplices.csv"))):
        m = pat.match(os.path.basename(p)); name, mat, thr = m["name"], m["mat"], float(m["thr"])
        wf = os.path.join(a.r_dir, f"{name}_W_{mat}.csv")
        W = pd.read_csv(wf, index_col=0).values
        P, R = clique_complex(W, thr), r_simplices(p)
        rows.append(dict(dataset=name, matrix=mat, threshold=thr, n_python=len(P), n_R=len(R),
                         only_python=len(P - R), only_R=len(R - P), identical=P == R))
    df = pd.DataFrame(rows); print(df.to_string(index=False))
    print(f"\n{int(df.identical.sum())} / {len(df)} conditions identical")
    sys.exit(0 if df.identical.all() else 1)
