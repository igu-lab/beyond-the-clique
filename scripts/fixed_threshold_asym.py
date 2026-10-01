"""fixed_threshold_asym.py — arXiv v1 表 2（閾値固定版）と有向性の定量。
(1) Saqr らの build_simplicial(tna, type="clique", threshold=t) と同じ 1-骨格（pmax 対称化した遷移確率 >= t）上で、
    クリーク複体と Dowker 複体（窓 w3/w4 の中で実際に共起した組だけを単体にする）の三角形の数と β0, β1 を比べる。
    クリーク側の三角形数が R の f_vector と一致することを確認する。
(2) 遷移行列の非対称性：両方向とも正のペアで max/min 比が k 倍以上の割合、片方向だけのペアの割合、
    反対称成分の割合 ||W - W^T||_F / ||W + W^T||_F。
使い方: python3 scripts/fixed_threshold_asym.py --out results/fixed_threshold [--r-dir results/r_saqr]
"""
import argparse, gzip, itertools as it, json, os
import numpy as np, pandas as pd
R = "results/r_saqr"

def rank2(rows):
    piv = {}; r = 0
    for v in rows:
        while v:
            h = v.bit_length() - 1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; r += 1; break
    return r

def betti01(V, E, T):
    vi = {v: i for i, v in enumerate(sorted(V))}; ei = {e: i for i, e in enumerate(sorted(E))}
    r1 = rank2([(1 << vi[a]) | (1 << vi[b]) for a, b in ei])
    r2 = rank2([sum(1 << ei[f] for f in it.combinations(t, 2)) for t in sorted(T)])
    return len(V) - r1, len(E) - r1 - r2

def windows(name, w):
    d = f"data/processed/{name}-w{w}/{name}-w{w}"
    nv = np.loadtxt(gzip.open(d + "-nverts.txt.gz"), dtype=int); sx = np.loadtxt(gzip.open(d + "-simplices.txt.gz"), dtype=int)
    out, p = set(), 0
    for k in nv:
        s = tuple(sorted(set(sx[p:p + k].tolist()))); p += k
        for t in it.combinations(s, 3): out.add(t)
    return out  # 1-based code ids (codes.json 順 = 列順)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--r-dir", default=R)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); R = a.r_dir
    rows, asym = [], []
    for name in ["group_regulation", "human_long"]:
        W = pd.read_csv(f"{R}/{name}_W_relative.csv", index_col=0).values
        n = W.shape[0]; S = np.maximum(W, W.T)
        # 非対称性
        pos = [(i, j) for i in range(n) for j in range(i + 1, n) if W[i, j] > 0 or W[j, i] > 0]
        both = [(i, j) for i, j in pos if W[i, j] > 0 and W[j, i] > 0]
        ratio = np.array([max(W[i, j], W[j, i]) / min(W[i, j], W[j, i]) for i, j in both])
        off = ~np.eye(n, dtype=bool)
        anti = np.linalg.norm((W - W.T)[off]) / np.linalg.norm((W + W.T)[off])
        asym.append(dict(dataset=name, pairs=len(pos), one_way=len(pos) - len(both),
                         ratio_ge2=float(np.mean(ratio >= 2)), ratio_ge3=float(np.mean(ratio >= 3)), ratio_ge5=float(np.mean(ratio >= 5)),
                         median_ratio=float(np.median(ratio)), antisym_share=float(anti)))
        for w in [3, 4]:
            tri_obs = windows(name, w)
            for t in [0.02, 0.05, 0.10, 0.20]:
                V = set(range(1, n + 1))
                E = {(i + 1, j + 1) for i in range(n) for j in range(i + 1, n) if S[i, j] >= t}
                TC = {c for c in it.combinations(sorted(V), 3) if all(e in E for e in it.combinations(c, 2))}
                TD = TC & tri_obs
                bC, bD = betti01(V, E, TC), betti01(V, E, TD)
                fv = pd.read_csv(f"{R}/{name}_sc_relative_t{t:.3f}_summary.csv", index_col=0).loc["f_vector", "value"].split("|")
                rows.append(dict(dataset=name, window=w, threshold=t, edges=len(E), tri_clique=len(TC), tri_clique_R=int(fv[2]) if len(fv) > 2 else 0,
                                 tri_dowker=len(TD), phantom=len(TC) - len(TD), b0_C=bC[0], b1_C=bC[1], b0_D=bD[0], b1_D=bD[1]))
    df = pd.DataFrame(rows); df.to_csv(os.path.join(a.out, "table2_fixed_threshold.csv"), index=False)
    da = pd.DataFrame(asym).round(6)   # 2026-10-01: 6 桁に丸める（BLAS の違いによる最後の桁の差をなくし，どの計算機でも同じ CSV にする）
    da.to_csv(os.path.join(a.out, "asymmetry.csv"), index=False)
    print(df.to_string(index=False)); print(da.round(3).to_string(index=False))
