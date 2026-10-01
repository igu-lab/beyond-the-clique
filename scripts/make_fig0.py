"""make_fig0.py — 説明図（図 0）：「3 人で話した」と「2 人ずつ 3 回話した」。英語版と日本語版。
2026-09-30_1148 改訂: プロジェクトの図の指示に合わせる（seaborn 'colorblind'，セリフ体（英: Times 系，和: Noto Serif CJK JP），
列見出し 14pt・その他のラベル 12pt・判定の枠 10pt，塗りの alpha 0.7，300 dpi で .pdf と .svg を保存）。色は図 1〜3 と同じ。
使い方（プロジェクトルートで）: python3 scripts/make_fig0.py   （既定 --out results/figures，固定名で保存．--stamp は任意）
"""
import argparse, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Polygon, Ellipse, Circle, FancyBboxPatch
CB = sns.color_palette("colorblind").as_hex()
C_D, C_C, INK, INK2, EDGE = CB[0], "#c0392b", "#000000", "#404040", "#7f7f7f"   # 2026-09-30_1840: クリークを橙→赤に
GFILL, GEDGE = "#e8e8e8", "#9a9a9a"      # 観測された集まりの楕円（無彩色）
SERIF_EN = ["Times New Roman", "Liberation Serif", "TeX Gyre Termes", "Nimbus Roman"]
plt.rcParams.update({"mathtext.fontset": "stix", "pdf.fonttype": 42, "svg.fonttype": "path", "savefig.dpi": 300})
T = {
 "en": dict(font=SERIF_EN, cols=["Observed", "Pairwise network", "Clique complex", "Dowker complex"],
            sub=["", "(pairs only)", "fill if all 3 edges exist", "fill only if all 3\nactually co-occurred"],
            rows=["Case 1", "Case 2"], obs=["one group: {A,B,C}", "three pairs: {A,B}, {B,C}, {A,C}"],
            hole="hole", verdict=["same graph", "both filled", "filled vs. hole"],
            verdict2=["cannot tell apart", "cannot tell apart", "can tell apart"]),
 "ja": dict(font=["Noto Serif CJK JP"], cols=["観測されたこと", "ペアのネットワーク", "クリーク複体", "Dowker 複体"],
            sub=["", "（ペアだけを見る）", "3 辺がそろえば塗る", "3 人が実際に一緒の\nときだけ塗る"],
            rows=["現実 1", "現実 2"], obs=["集まり 1 回：{A,B,C}", "集まり 3 回：{A,B}, {B,C}, {A,C}"],
            hole="穴", verdict=["同じ図になる", "両方とも塗られる", "上は塗られ，下は穴"],
            verdict2=["区別できない", "区別できない", "区別できる"]),
}
P = {"A": (0.0, 0.0), "B": (1.0, 0.0), "C": (0.5, 0.87)}

def tri(ax, x0, y0, fill=None, hole=None):
    pts = [(x0 + P[k][0], y0 + P[k][1]) for k in "ABC"]
    if fill: ax.add_patch(Polygon(pts, closed=True, fc=fill, ec="none", alpha=0.7))
    for i in range(3):
        a, b = pts[i], pts[(i + 1) % 3]; ax.plot([a[0], b[0]], [a[1], b[1]], color=EDGE, lw=2.2, zorder=2)
    for k, (x, y) in zip("ABC", pts):
        ax.add_patch(Circle((x, y), 0.19, fc="white", ec=INK, lw=1.4, zorder=3))
        ax.text(x, y, k, ha="center", va="center", fontsize=12, fontweight="bold", color=INK, zorder=4)
    if hole: ax.text(x0 + 0.5, y0 + 0.3, hole, ha="center", va="center", fontsize=12, color=C_D, fontweight="bold")

def main(lang, out):
    t = T[lang]; plt.rcParams["font.family"] = "serif"; plt.rcParams["font.serif"] = t["font"]
    # Noto Serif CJK は CFF 形式の OpenType のため，Type 42 だと PDF で不整合警告が出る．和文版だけ Type 3 で埋め込む
    plt.rcParams["pdf.fonttype"] = 3 if lang == "ja" else 42
    fig, ax = plt.subplots(figsize=(10.0, 5.6)); ax.set_xlim(-0.7, 9.1); ax.set_ylim(-1.7, 3.85); ax.axis("off"); ax.set_aspect("equal")
    xs = [0.6, 2.9, 5.2, 7.5]
    for x, c, s, col in zip(xs, t["cols"], t["sub"], [INK, INK, C_C, C_D]):
        ax.text(x + 0.5, 3.55, c, ha="center", fontsize=14, color=col, fontweight="bold")
        if s: ax.text(x + 0.5, 3.3, s, ha="center", va="top", fontsize=12, color=INK2, linespacing=1.2)
    ys = [1.55, -0.35]
    for r, y in enumerate(ys):
        ax.text(-0.5, y + 0.43, t["rows"][r], rotation=90, ha="center", va="center", fontsize=14, fontweight="bold", color=INK)
        if r == 0:
            ax.add_patch(Ellipse((xs[0] + 0.5, y + 0.43), 1.8, 0.75, fc=GFILL, ec=GEDGE))
            for k, dx in zip("ABC", (-0.45, 0, 0.45)):
                ax.add_patch(Circle((xs[0] + 0.5 + dx, y + 0.43), 0.17, fc="white", ec=INK, lw=1.2))
                ax.text(xs[0] + 0.5 + dx, y + 0.43, k, ha="center", va="center", fontsize=12, fontweight="bold")
        else:
            for (cx, cy), pair in zip([(0.0, 0.75), (1.0, 0.75), (0.5, 0.05)], ["AB", "BC", "AC"]):
                ax.add_patch(Ellipse((xs[0] + cx, y + cy), 1.0, 0.55, fc=GFILL, ec=GEDGE))
                for k, dx in zip(pair, (-0.2, 0.2)):
                    ax.add_patch(Circle((xs[0] + cx + dx, y + cy), 0.15, fc="white", ec=INK, lw=1.1))
                    ax.text(xs[0] + cx + dx, y + cy, k, ha="center", va="center", fontsize=10, fontweight="bold")
        ax.text(xs[0] + 0.5, y - (0.05 if r == 0 else 0.25), t["obs"][r], ha="center", va="top", fontsize=12, color=INK2)
        tri(ax, xs[1], y)
        tri(ax, xs[2], y, fill=C_C)
        tri(ax, xs[3], y, fill=C_D if r == 0 else None, hole=t["hole"] if r == 1 else None)
    for x, v, v2, col in zip(xs[1:], t["verdict"], t["verdict2"], [INK2, C_C, C_D]):
        ax.add_patch(FancyBboxPatch((x - 0.5, -1.65), 2.0, 0.72, boxstyle="round,pad=0.02", fc="white", ec=col, lw=1.3))
        ax.text(x + 0.5, -1.12, v, ha="center", va="center", fontsize=10, color=col)
        ax.text(x + 0.5, -1.42, v2, ha="center", va="center", fontsize=10, color=col, fontweight="bold")
    for ext in ("pdf", "svg"): fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="results/figures"); ap.add_argument("--stamp", default=""); a = ap.parse_args()
    sfx = f"_{a.stamp}" if a.stamp else ""
    os.makedirs(a.out, exist_ok=True)
    for lang in ("en", "ja"): main(lang, os.path.join(a.out, f"fig0_core_idea_{lang}{sfx}"))
