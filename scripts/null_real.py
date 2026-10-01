"""null_real.py — 実データ 6 本（ScHoLP 4 本＋TNA 同梱 2 本）に次数・サイズ保存ヌルモデルを当てる（arXiv v1 用）。
独立サンプリング（src/null_model.py の chain=False）。
使い方: python3 scripts/null_real.py --scholp-root <ScHoLP-Data のパス> --out results/null_real [--N 200] [--only name ...]
"""
import argparse, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from null_model import run_null, load_scholp

SCHOLP = ["tags-matheducators", "threads-matheducators", "tags-academia", "tags-physics"]
TNA = ["group_regulation-w3", "human_long-w3", "group_regulation-w4", "human_long-w4"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scholp-root", required=True)
    ap.add_argument("--processed-root", default=os.path.join(HERE, "..", "data", "processed"))
    ap.add_argument("--out", required=True); ap.add_argument("--N", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0); ap.add_argument("--only", nargs="*")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    jobs = [(n, os.path.join(a.scholp_root, n)) for n in SCHOLP] + [(n, os.path.join(a.processed_root, n)) for n in TNA]
    for n, d in jobs:
        if a.only and n not in a.only: continue
        _, H = load_scholp(d)
        t = time.time()
        s, _, _ = run_null(H, N=a.N, seed=a.seed, out=a.out, label=n)
        print(n, "done", round(time.time() - t, 1), "s", json.dumps({k: {kk: round(vv, 3) for kk, vv in s[k].items()}
              for k in ["phantom_frac", "max_b1_C", "max_b1_D", "kernel_dim", "premature", "gap_L1"]}), flush=True)
