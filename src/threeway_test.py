"""threeway_test.py — ペアの共起を保ったうえでの既約な三体相互作用の検定（2026-09-29 作成）。

各三つ組 {a,b,c}（三辺とも共起している＝クリーク複体の三角形）について、ハイパーエッジ（観測された集まり）を
単位とする 2×2×2 分割表（a, b, c の在・不在）を作り、「三次の交互作用なし」の対数線形モデル
（すべての二次周辺＝ペアの共起を保存）を IPF で当てはめる。観測された三つ組の共起数 n_abc を期待値 m_abc と比べ、
尤度比統計量 G²（自由度 1）と符号（正＝ペアから期待されるより三者が一緒に起きる）を返す。

  n_abc > m_abc かつ有意 → 三者の共起はペアの寄せ集めでは説明できない（既約な正の三体相互作用）
  Dowker の三角形（n_abc ≥ 1）でも n_abc < m_abc なら、ペアから期待されるより「三者一緒」が少ない
  phantom 三角形（n_abc = 0）で m_abc が大きい → ペアからは三者の共起が期待されるのに一度も起きない

使い方: from threeway_test import threeway_table; df = threeway_table(H)
"""
import itertools as it
from collections import Counter
import numpy as np
import pandas as pd
from math import erfc, sqrt


def counts(H):
    n1, n2, n3 = Counter(), Counter(), Counter()
    for S in H:
        s = sorted(set(S))
        n1.update(s); n2.update(it.combinations(s, 2)); n3.update(it.combinations(s, 3))
    return n1, n2, n3


def ipf_no3(tab, iters=200, tol=1e-10):
    """tab: (T,2,2,2) 観測度数。二次周辺を保つ三次交互作用なしモデルの期待度数（IPF）。"""
    m = np.ones_like(tab, dtype=float)
    for _ in range(iters):
        old = m.copy()
        for ax in (1, 2, 3):   # ax を潰した二次周辺（残り 2 変数）を合わせる
            obs = tab.sum(axis=ax, keepdims=True); fit = m.sum(axis=ax, keepdims=True)
            m = m * np.divide(obs, fit, out=np.zeros_like(fit), where=fit > 0)
        if np.max(np.abs(m - old)) < tol:
            break
    return m


def threeway_table(H, triples=None):
    N = len(H)
    n1, n2, n3 = counts(H)
    if triples is None:   # クリーク複体の三角形（三辺とも共起）
        adj = {}
        for (a, b) in n2:
            adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
        triples = [(a, b, c) for (a, b) in n2 for c in adj[a] & adj[b] if c > b]
    T = np.array(triples)
    a, b, c = T[:, 0], T[:, 1], T[:, 2]
    g1 = lambda x: np.array([n1[v] for v in x], float)
    g2 = lambda x, y: np.array([n2[(min(u, v), max(u, v))] for u, v in zip(x, y)], float)
    na, nb, nc = g1(a), g1(b), g1(c)
    nab, nac, nbc = g2(a, b), g2(a, c), g2(b, c)
    nabc = np.array([n3[tuple(t)] for t in map(tuple, T)], float)
    tab = np.zeros((len(T), 2, 2, 2))
    tab[:, 1, 1, 1] = nabc
    tab[:, 1, 1, 0] = nab - nabc; tab[:, 1, 0, 1] = nac - nabc; tab[:, 0, 1, 1] = nbc - nabc
    tab[:, 1, 0, 0] = na - nab - nac + nabc; tab[:, 0, 1, 0] = nb - nab - nbc + nabc
    tab[:, 0, 0, 1] = nc - nac - nbc + nabc
    tab[:, 0, 0, 0] = N - (na + nb + nc - nab - nac - nbc + nabc)
    m = ipf_no3(tab)
    with np.errstate(divide="ignore", invalid="ignore"):
        G2 = 2 * np.nansum(np.where(tab > 0, tab * np.log(tab / m), 0.0), axis=(1, 2, 3))
    df = pd.DataFrame(dict(a=a, b=b, c=c, n_abc=nabc, m_abc=m[:, 1, 1, 1], G2=G2))
    df["p"] = [erfc(sqrt(max(g, 0.0) / 2.0)) for g in df.G2]   # 自由度 1 の χ² の上側確率
    df["sign"] = np.sign(df.n_abc - df.m_abc)
    df["dowker"] = df.n_abc > 0
    return df
