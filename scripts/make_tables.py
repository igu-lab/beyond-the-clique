"""make_tables.py（threads-matheducators はハイパーエッジのサイズ上限なし＝null_model・threeway と同じ条件の結果を使う）
 — arXiv v1 の表 1（構成の比較とヌルモデル）・表 2（閾値固定版）・表 3（三次交互作用）を CSV と LaTeX で出す。
使い方（プロジェクトルートで）: python3 scripts/make_tables.py --out results/tables 
（2026-10-01: --stamp を省くと固定名 table1_construction_null.csv などで保存）
"""
import argparse, json, os
import pandas as pd
DS = [("tags-matheducators", "tags-matheducators"), ("threads-matheducators", "threads-matheducators"),
      ("tags-academia", "tags-academia"), ("tags-physics", "tags-physics"),
      ("group_regulation-w3", "group\\_regulation ($w{=}3$)"), ("human_long-w3", "human\\_long ($w{=}3$)")]
PH = {"tags-matheducators": "results/ph_compare_se", "threads-matheducators": "results/ph_compare_se_nocap",
      "tags-academia": "results/ph_compare_se", "tags-physics": "results/ph_compare_se",
      "group_regulation-w3": "results/ph_compare_tna", "human_long-w3": "results/ph_compare_tna"}
NULL = ["results/null_real", "results/null_real/physics_N200"]
TW = "results/threeway"

def nullsum(k):
    for d in NULL:
        f = os.path.join(d, f"{k}_null_summary.json")
        if os.path.exists(f): return json.load(open(f))
    return None

def zfmt(s, key):
    if s is None: return "--"
    z = s[key]["z"]
    return "n/a" if z != z else f"{z:.1f}"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--stamp", default=""); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    r1 = []
    for k, lab in DS:
        p = json.load(open(os.path.join(PH[k], f"{k}_count_summary.json"))); s = nullsum(k)
        r1.append(dict(dataset=lab, V=p["V"], E=p["E"], tri_clique=p["tri_clique"], tri_dowker=p["tri_dowker"],
                       phantom_rate=round(p["phantom_triangle_frac"], 3),
                       phantom_rate_null=("--" if s is None else round(s["phantom_frac"]["null_mean"], 3)),
                       max_b1_clique=p["betti1_max_clique"], max_b1_dowker=p["betti1_max_dowker"],
                       max_b1_dowker_null=("--" if s is None else round(s["max_b1_D"]["null_mean"], 1)),
                       z_dowker=zfmt(s, "max_b1_D"), z_clique=zfmt(s, "max_b1_C"),
                       kernel_dim=p["h1_phantom_only"], clique_earlier=p["h1_premature"], h1_classes=p["h1_births"],
                       null_N=("--" if s is None else s["N"])))
    t1 = pd.DataFrame(r1)
    r3 = []
    for k, lab in DS:
        f = os.path.join(TW, f"{k}_threeway_summary.json")
        if not os.path.exists(f): r3.append(dict(dataset=lab)); continue
        s = json.load(open(f))
        r3.append(dict(dataset=lab, clique_triangles=s["clique_triangles"], dowker_triangles=s["dowker_triangles"],
                       OE_all=round(s["dowker_obs_sum"] / (s["dowker_expected_sum"] + s["phantom_expected_sum"]), 2),  # 全クリーク三角形（phantom は観測 0）
                       OE_dowker=round(s["dowker_obs_sum"] / s["dowker_expected_sum"], 2),  # n_abc>=1 で選んだ三角形（選択の偏りで上に振れる）
                       pos_sig=s["dowker_pos_sig"], neg_sig=s["dowker_neg_sig"] + s["phantom_neg_sig"],
                       phantom_expected_triple_events=round(s["phantom_expected_sum"], 1)))
    t3 = pd.DataFrame(r3)
    t2 = pd.read_csv("results/fixed_threshold/table2_fixed_threshold.csv")
    t2 = t2[t2.window == 3].drop(columns=["tri_clique_R"])
    sfx = f"_{a.stamp}" if a.stamp else ""
    for name, t in (("table1_construction_null", t1), ("table2_fixed_threshold", t2), ("table3_threeway", t3)):
        t.to_csv(os.path.join(a.out, f"{name}{sfx}.csv"), index=False)
        with open(os.path.join(a.out, f"{name}{sfx}.tex"), "w") as fh:   # booktabs 形式（jinja2 に依存しない）
            cols = [c.replace("_", "\\_") for c in t.columns]
            fh.write("\\begin{tabular}{l" + "r" * (len(cols) - 1) + "}\n\\toprule\n" + " & ".join(cols) + " \\\\\n\\midrule\n")
            for _, r in t.iterrows():
                fh.write(" & ".join("--" if (isinstance(v, float) and v != v) else (f"{v:,}" if isinstance(v, int) else str(v)) for v in r.values) + " \\\\\n")
            fh.write("\\bottomrule\n\\end{tabular}\n")
        print(name); print(t.to_string(index=False))
