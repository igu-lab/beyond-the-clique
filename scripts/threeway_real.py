"""threeway_real.py — 実データに三次交互作用の検定（src/threeway_test.py）を当て、要約と三つ組ごとの表を保存する。
使い方: python3 scripts/threeway_real.py --scholp-root <ScHoLP-Data> --out results/threeway [--only ...]
TNA の窓データは窓が重なるので観測が独立でない → p 値は記述的な目安として扱う。
"""
import argparse, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "src"))
from threeway_test import threeway_table
from null_model import load_scholp
import numpy as np
SCHOLP = ["tags-matheducators", "threads-matheducators", "tags-academia", "tags-physics"]
TNA = ["group_regulation-w3", "human_long-w3", "group_regulation-w4", "human_long-w4"]
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--scholp-root", required=True)
    ap.add_argument("--processed-root", default=os.path.join(HERE, "..", "data", "processed"))
    ap.add_argument("--out", required=True); ap.add_argument("--only", nargs="*"); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for n, d in [(n, os.path.join(a.scholp_root, n)) for n in SCHOLP] + [(n, os.path.join(a.processed_root, n)) for n in TNA]:
        if a.only and n not in a.only: continue
        _, H = load_scholp(d); t = time.time(); df = threeway_table(H)
        bonf = 0.05 / len(df); sig = df.p < bonf; D = df.dowker
        s = dict(dataset=n, hyperedges=len(H), clique_triangles=len(df), dowker_triangles=int(D.sum()), phantom=int((~D).sum()),
                 bonferroni_alpha=bonf,
                 dowker_pos_sig=int((D & sig & (df.sign > 0)).sum()), dowker_neg_sig=int((D & sig & (df.sign < 0)).sum()),
                 phantom_neg_sig=int((~D & sig).sum()),
                 phantom_expected_ge1=int((~D & (df.m_abc >= 1)).sum()),
                 phantom_expected_sum=float(df.m_abc[~D].sum()),
                 dowker_obs_sum=float(df.n_abc[D].sum()), dowker_expected_sum=float(df.m_abc[D].sum()),
                 seconds=round(time.time() - t, 1))
        json.dump(s, open(os.path.join(a.out, f"{n}_threeway_summary.json"), "w"), indent=1)
        df.to_csv(os.path.join(a.out, f"{n}_threeway_triples.csv.gz"), index=False, compression="gzip")
        print(json.dumps(s), flush=True)
