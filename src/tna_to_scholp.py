"""tna_to_scholp.py — R 側（scripts/r/run_saqr_pipeline.R）が書き出した長形式 CSV
(actor, time, action) を ScHoLP 形式（<name>-nverts/-simplices/-times .txt.gz）に変換する。

単位（ハイパーエッジ）の切り方 --unit:
  seq   : 行為者（系列）ごとに出現したコードの集合を 1 つのハイパーエッジにする
          （Nestimate::build_network(method="co_occurrence") が対で数えている単位と同じ）
  w<k>  : 長さ k のスライディング窓内のコード集合（k=2,3,4,...）

出力: <out>/<name>-<unit>/<name>-<unit>-{nverts,simplices,times}.txt.gz と codes.json（頂点番号→コード名）
使い方: python3 src/tna_to_scholp.py results/r_saqr/group_regulation_long.csv --unit seq w3 w4 --out data/processed
"""
import argparse, csv, gzip, json, os
from collections import defaultdict


def read_long(path):
    seqs = defaultdict(list)
    for row in csv.DictReader(open(path)):
        seqs[row["actor"]].append((float(row["time"]), row["action"]))
    return {a: [x for _, x in sorted(v)] for a, v in seqs.items()}


def units(seqs, unit):
    """yield (set_of_codes, time_stamp)"""
    for ai, (a, s) in enumerate(sorted(seqs.items())):
        if unit == "seq":
            yield set(s), ai
        else:
            k = int(unit[1:])
            for i in range(0, max(len(s) - k + 1, 1)):
                yield set(s[i:i + k]), ai * 10000 + i


def write_scholp(seqs, unit, name, out):
    codes = sorted({c for s in seqs.values() for c in s})
    idx = {c: i + 1 for i, c in enumerate(codes)}
    d = os.path.join(out, f"{name}-{unit}")
    os.makedirs(d, exist_ok=True)
    base = os.path.join(d, f"{name}-{unit}")
    n = 0
    with gzip.open(base + "-nverts.txt.gz", "wt") as fv, gzip.open(base + "-simplices.txt.gz", "wt") as fs, gzip.open(base + "-times.txt.gz", "wt") as ft:
        for S, t in units(seqs, unit):
            if len(S) < 1:
                continue
            fv.write(f"{len(S)}\n")
            for c in sorted(S, key=lambda c: idx[c]):
                fs.write(f"{idx[c]}\n")
            ft.write(f"{t}\n")
            n += 1
    json.dump({str(v): k for k, v in idx.items()}, open(os.path.join(d, "codes.json"), "w"), ensure_ascii=False, indent=1)
    print(f"{name}-{unit}: {n} hyperedges, {len(codes)} codes -> {d}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--unit", nargs="+", default=["seq", "w3", "w4"])
    ap.add_argument("--name", default=None)
    ap.add_argument("--out", default="data/processed")
    a = ap.parse_args()
    name = a.name or os.path.basename(a.csv).replace("_long.csv", "")
    seqs = read_long(a.csv)
    for u in a.unit:
        write_scholp(seqs, u, name, a.out)
