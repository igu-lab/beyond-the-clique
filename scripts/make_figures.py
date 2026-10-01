"""make_figures.py — 論文の図 2〜4（ファイル名は fig1〜fig3）．2026-09-29 作成，2026-09-30 にプロジェクトの図の指示に合わせて改訂．
出力は固定名（--stamp を付けたときだけ日時付き）
図 1: β₁ 曲線（観測とヌルの 95% 帯、Dowker とクリーク）6 パネル
図 2: H₁ クラスの内訳（両方で同時に消える／クリークが早く埋める／クリークだけが埋める＝核／両方で残る）
図 3: max β₁ のヌルに対する z 値（Dowker とクリーク）
図 fig_phantom_maxb1（2026-10-01 追加，論文の図 2）: a 見かけの三角形の割合，b 最大 β₁
スタイル: 背景色・グリッドなし，上と右の枠線なし，seaborn 'colorblind'，セリフ体（Times 系），
          軸ラベル 12pt・タイトル 14pt・凡例 10pt，重なる要素は alpha 0.7〜0.8，300 dpi で .pdf と .svg を保存
使い方（プロジェクトルートで）:
  python3 scripts/make_figures.py
  （既定: --null results/null_real
          --bars results/ph_compare_se_nocap results/ph_compare_se results/ph_compare_tna
          --out results/figures）
"""
import argparse, json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

CB = sns.color_palette("colorblind").as_hex()
C_D, C_C = CB[0], "#c0392b"              # Dowker（青）, clique（赤）．2026-09-30_1840: 橙→赤
INK, INK2, BAND = "#000000", "#404040", "#e6e6e6"
CAT4 = [CB[0], C_C, CB[2], "#8a8a8a"]   # 2026-09-30_1840: 「クリークが先」をクリークの赤に，「どちらでも埋まらない」を灰に（橙・朱をなくす）
SERIF = ["Times New Roman", "Liberation Serif", "TeX Gyre Termes", "Nimbus Roman", "DejaVu Serif"]
plt.rcParams.update({
    "font.family": "serif", "font.serif": SERIF, "mathtext.fontset": "stix",
    "axes.labelsize": 12, "axes.titlesize": 14, "legend.fontsize": 10,
    "xtick.labelsize": 10, "ytick.labelsize": 10,
    "axes.facecolor": "white", "figure.facecolor": "white", "axes.grid": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.linewidth": 0.8, "pdf.fonttype": 42, "svg.fonttype": "path", "savefig.dpi": 300})
ALPHA = 0.8

def save(fig, out):
    for ext in ("pdf", "svg"): fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")

DS = [("tags-matheducators", "tags-matheducators"), ("threads-matheducators", "threads-matheducators"),
      ("tags-academia", "tags-academia"), ("tags-physics", "tags-physics"),
      ("group_regulation-w3", "group_regulation (w=3)"), ("human_long-w3", "human_long (w=3)")]

def find(nulldir, name):
    """ヌルの結果を nulldir 直下か、その下の physics_N200/ から探す"""
    for d in (nulldir, os.path.join(nulldir, "physics_N200")):
        f = os.path.join(d, name)
        if os.path.exists(f): return f
    return os.path.join(nulldir, name)

def style(ax):
    ax.grid(False)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

def fig1(nulldir, out):
    fig, axes = plt.subplots(2, 3, figsize=(9.0, 6.2), constrained_layout=True)
    for ax, (key, title) in zip(axes.flat, DS):
        f = find(nulldir, f"{key}_null_b1curve_edgefrac.csv")
        pending = not os.path.exists(f)
        d = pd.read_csv(f if not pending else os.path.join(nulldir, f"{key}_obs_edgefrac.csv"))
        x = d.edge_frac
        if not pending:
            ax.fill_between(x, d["nullD_lo2.5"], d["nullD_hi97.5"], color=C_D, alpha=0.2, lw=0)
            ax.fill_between(x, d["nullC_lo2.5"], d["nullC_hi97.5"], color=C_C, alpha=0.2, lw=0)
            ax.plot(x, d.nullD_median, color=C_D, lw=1.0, ls=(0, (3, 2)), alpha=ALPHA)
        ax.plot(x, d.b1_dowker_obs, color=C_D, lw=1.8, alpha=ALPHA, label="Dowker (observed)")
        ax.plot(x, d.b1_clique_obs, color=C_C, lw=1.8, alpha=ALPHA, label="clique (observed)")
        ax.set_title(title + ("  [null pending]" if pending else ""), color=INK, loc="left", fontsize=14)
        ax.set_xlim(0, 1); ax.set_ylim(bottom=0); style(ax)
    for ax in axes[1]: ax.set_xlabel("fraction of edges in the 1-skeleton\n(strongest co-occurrences first)")
    for ax in axes[:, 0]: ax.set_ylabel(r"$\beta_1$")
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=C_D, lw=1.8), Line2D([], [], color=C_C, lw=1.8),
         Patch(color=C_D, alpha=0.2), Line2D([], [], color=C_D, lw=1.0, ls=(0, (3, 2))), Patch(color=C_C, alpha=0.2)]
    fig.legend(h, ["Dowker, observed", "clique, observed", "Dowker, null 95%", "Dowker, null median", "clique, null 95%"],
               loc="outside upper center", ncol=3, frameon=False, fontsize=10)
    save(fig, out)
    plt.close(fig)

def classes(bars):
    b = pd.read_csv(bars); dC, dD = b.death_clique.values, b.death_dowker.values
    fC, fD = np.isfinite(dC), np.isfinite(dD)
    both = fC & fD
    return {"filled at the same step": int((both & (dC == dD)).sum()),
            "clique fills earlier": int((both & (dC < dD)).sum()),
            "only clique fills (kernel)": int((fC & ~fD).sum()),
            "never filled in either": int((~fC & ~fD).sum())}

def txt(hexcol):
    """背景色に対して読みやすい文字色（相対輝度で白か黒）"""
    r, g, b = [int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in (r, g, b)]
    return "black" if 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2] > 0.18 else "white"

def fig2(barsdir, out):
    rows = []
    for key, title in DS:
        for bd in barsdir:
            f = os.path.join(bd, f"{key}_count_bars.csv")
            if os.path.exists(f): rows.append((title, classes(f))); break
    cats = list(rows[0][1].keys())
    fig, ax = plt.subplots(figsize=(9.0, 3.8), constrained_layout=True)
    for i, (title, c) in enumerate(rows[::-1]):
        tot = sum(c.values()); left = 0
        for j, k in enumerate(cats):
            w = c[k] / tot
            ax.barh(i, w, left=left, color=CAT4[j], height=0.62, edgecolor="white", linewidth=1.0)
            if w > 0.07: ax.text(left + w / 2, i, f"{c[k]:,}", ha="center", va="center", fontsize=10, color=txt(CAT4[j]))
            left += w
        ax.text(1.01, i, f"$n$ = {tot:,}", va="center", fontsize=10, color=INK2)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([t for t, _ in rows[::-1]])
    ax.set_xlim(0, 1); ax.set_xlabel(r"share of $H_1$ classes born in the common 1-skeleton")
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=12)
    from matplotlib.patches import Patch
    fig.legend([Patch(color=CAT4[j]) for j in range(4)], cats, loc="outside upper center", ncol=2, frameon=False, fontsize=10)
    save(fig, out)
    plt.close(fig)
    return rows

def fig3(nulldir, out):
    fig, ax = plt.subplots(figsize=(6.4, 3.9), constrained_layout=True)
    ax.axvspan(-1.96, 1.96, color=BAND, lw=0, zorder=0)   # 有意でない範囲（データ上の帯）
    labels = []
    for i, (key, title) in enumerate(DS[::-1]):
        f = find(nulldir, f"{key}_null_summary.json"); labels.append(title)
        if not os.path.exists(f):
            ax.text(0, i, "null pending", ha="center", va="center", fontsize=10, color=INK2); continue
        s = json.load(open(f))
        ax.axhline(i, color="#808080", lw=0.8, ls=(0, (1, 2)), zorder=1.5)   # 2026-09-30_1900: 行の補助線
        # 2026-10-01: 印を上下にずらさず補助線の上に載せ，Dowker とクリークを細い線で結ぶ（ダンベル型）
        zs = [s[k]["z"] if (s[k]["z"] is not None and np.isfinite(s[k]["z"])) else 0.0 for k in ("max_b1_D", "max_b1_C")]
        ax.plot(zs, [i, i], color="#9a9a9a", lw=1.0, solid_capstyle="butt", zorder=2)
        for k, col, mk, zo in (("max_b1_C", C_C, "D", 3), ("max_b1_D", C_D, "o", 4)):
            z = s[k]["z"]
            if z is None or not np.isfinite(z):
                ax.plot(0, i, marker=mk, ms=7, mfc="white", mec=col, mew=1.4, ls="", zorder=zo)
            else:
                ax.plot(z, i, marker=mk, ms=7, color=col, mec="white", mew=0.8, ls="", zorder=zo)
    ax.set_xscale("symlog", linthresh=2)
    tk = [-200, -50, -10, -2, 0, 2, 10, 100]; ax.set_xticks(tk); ax.set_xticklabels([str(t).replace("-", "\u2212") for t in tk]); ax.minorticks_off()
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel(r"$z$ of max $\beta_1$ against the degree- and size-preserving null")
    style(ax); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", color=C_D, ls="", ms=7), Line2D([], [], marker="D", color=C_C, ls="", ms=7),
         Line2D([], [], marker="o", mfc="white", mec=INK2, ls="", ms=7)]
    fig.legend(h, ["Dowker", "clique", r"$\beta_1=0$ in data and null"], loc="outside upper center", ncol=3, frameon=False, fontsize=10)
    save(fig, out)
    plt.close(fig)

def fig_phantom(table, out):
    """2026-10-01 追加: 結果の最初の図．a 見かけの三角形の割合（観測とヌル平均），b 最大 β₁（クリークと Dowker）．
    数値は表1・表2と同じ CSV から読む"""
    t = pd.read_csv(table)
    t["key"] = [r.replace("\\_", "_").replace("($w{=}3$)", "(w=3)") for r in t.dataset]
    t = t.set_index("key").loc[[title for _, title in DS]]
    n = len(t); y = np.arange(n)[::-1]          # 上から DS の順
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.9), sharey=True, constrained_layout=True,
                               gridspec_kw={"width_ratios": [1, 1], "wspace": 0.08})
    for ax in (a, b):
        style(ax); ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
        ax.axhline(1.5, color="#808080", lw=0.8, ls=(0, (1, 2)))   # Stack Exchange と TNA の区切り
    # a: 見かけの三角形の割合
    a.barh(y, t.phantom_rate, height=0.6, color=C_C, alpha=ALPHA, lw=0, label="observed")
    a.errorbar(t.phantom_rate_null, y, yerr=0.36, fmt="none", ecolor=INK, elinewidth=2.0, capsize=0, label="null mean")
    for yi, v, nm in zip(y, t.phantom_rate, t.phantom_rate_null):
        a.text(max(v, nm) + 0.03, yi, f"{v:.2f}", va="center", ha="left", fontsize=10, color=INK2)
    a.set_xlim(0, 1.0); a.set_xticks([0, 0.25, 0.5, 0.75, 1.0]); a.set_xticklabels(["0", "0.25", "0.50", "0.75", "1"])
    a.set_xlabel("phantom rate\n(share of clique triangles never observed together)")
    a.set_yticks(y); a.set_yticklabels(t.index, fontsize=12)
    a.set_title("a  Phantom triangles of the clique complex", loc="left", fontsize=12, fontweight="bold")
    a.legend(loc="lower right", frameon=False, fontsize=10, handlelength=1.2)
    # b: 最大 β₁（ダンベル型，symlog で 0 も表示）
    for yi, (c, d) in zip(y, zip(t.max_b1_clique, t.max_b1_dowker)):
        b.plot([c, d], [yi, yi], color="#9a9a9a", lw=1.0, zorder=2)
        if c == 0: b.plot(c, yi, marker="D", ms=7, mfc="white", mec=C_C, mew=1.4, ls="", zorder=3)
        else:      b.plot(c, yi, marker="D", ms=7, color=C_C, mec="white", mew=0.8, ls="", zorder=3)
        b.plot(d, yi, marker="o", ms=7, color=C_D, mec="white", mew=0.8, ls="", zorder=4)
        b.text(c * 0.78 if c > 1 else -0.25, yi, f"{c:,}", va="center", ha="right", fontsize=10, color=C_C)
        b.text(d * 1.3, yi, f"{d:,}", va="center", ha="left", fontsize=10, color=C_D)
    b.set_xscale("symlog", linthresh=1, linscale=0.6)
    b.set_xlim(-0.6, 1e5); tk = [0, 1, 10, 100, 1000, 10000, 100000]
    b.set_xticks(tk); b.set_xticklabels(["0", "1", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"]); b.minorticks_off()
    b.set_xlabel(r"max $\beta_1$ along the filtration (log scale)")
    b.set_title(r"b  Holes found in the same data", loc="left", fontsize=12, fontweight="bold")
    from matplotlib.lines import Line2D
    b.legend([Line2D([], [], marker="D", color=C_C, ls="", ms=7), Line2D([], [], marker="o", color=C_D, ls="", ms=7)],
             ["clique complex", "Dowker complex"], loc="lower right", frameon=False, fontsize=10, handletextpad=0.3)
    save(fig, out)
    plt.close(fig)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", default="results/null_real")
    ap.add_argument("--bars", nargs="+", default=["results/ph_compare_se_nocap",
                    "results/ph_compare_se", "results/ph_compare_tna"])
    ap.add_argument("--table", default="results/tables/table1_construction_null.csv")
    ap.add_argument("--out", default="results/figures")
    ap.add_argument("--stamp", default="", help="付けると fig*_<stamp> という日時付きの名前で保存（通常は不要）")
    a = ap.parse_args()
    sfx = f"_{a.stamp}" if a.stamp else ""
    os.makedirs(a.out, exist_ok=True)
    fig_phantom(a.table, os.path.join(a.out, f"fig_phantom_maxb1{sfx}"))
    fig1(a.null, os.path.join(a.out, f"fig1_betti1_null{sfx}"))
    rows = fig2(a.bars, os.path.join(a.out, f"fig2_h1_classes{sfx}"))
    fig3(a.null, os.path.join(a.out, f"fig3_zscores{sfx}"))
    for t, c in rows: print(t, c)
