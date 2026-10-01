# Beyond the Clique（クリーク複体を超えて）

[English README](README.md)

次の論文のコードです．

> 安武 公一，辻 若菜，井上 仁（2026）．
> **Beyond the Clique: Comparing Clique and Dowker Complexes for Co-occurrence Data in Learning Analytics.**
> arXiv：*公開後に追記*

学習分析では，共起データ（談話の窓に現れるコード，スレッドの参加者，投稿のタグなど）を単体複体で表し，パーシステントホモロジーで分析することが増えています．標準的な**クリーク複体**はペアのネットワークだけで決まるため，「3要素が同時に共起した」ことと「3つのペアがそれぞれ別々に共起した」ことを区別できません．**Dowker 複体**は，実際に一緒に観測された組を単体とします．このリポジトリには，Stack Exchange の 4 つのデータセットと R パッケージ `tna`・`Nestimate` の例示データで二つの複体を比較し，論文のすべての図表を再生成するコードが入っています．

## 構成

| パス | 内容 |
|---|---|
| `src/ph_compare.py` | 同じ 1-骨格の上でのクリーク複体と Dowker 複体のパーシステントホモロジー（H0, H1）．Python/NumPy のみ |
| `src/null_model.py` | 次数とサイズを保存するハイパーグラフのヌルモデル |
| `src/threeway_test.py` | 三次の交互作用の検定（すべてのペアの周辺を保つ対数線形モデル） |
| `src/stackexchange_to_scholp.py` | Stack Exchange の `Posts.xml` → ScHoLP 形式（`nverts`／`simplices`／`times`） |
| `src/tna_to_scholp.py` | TNA の系列（長形式）→ 窓に切った ScHoLP 形式 |
| `scripts/r/run_saqr_pipeline.R` | `tna` → `Nestimate::build_simplicial`／`persistent_homology` を実行し，結果を書き出す |
| `scripts/check_nestimate_clique.py` | ここでのクリーク複体の構成が `Nestimate::build_simplicial` と一致することの確認 |
| `scripts/null_real.py`，`scripts/threeway_real.py`，`scripts/fixed_threshold_asym.py` | 実データの分析 |
| `scripts/make_tables.py`，`scripts/make_figures.py`，`scripts/make_fig0.py` | 論文の表と図 |
| `scripts/run_all.py` | すべてを順に実行 |
| `results/` | 論文で使った結果（下記） |

## 必要なもの

- Python 3.10 以上と `numpy`，`pandas`，`matplotlib`，`seaborn`（`pip install -r requirements.txt` または `uv sync`）
- TNA データを使う場合のみ：R と `tna` 1.2.3，`Nestimate` 0.8.5
- 図 1 の日本語版はフォント Noto Serif CJK JP を使います（英語版は Times 系のセリフ体）

## データ

データはこのリポジトリでは**再配布していません**．

- **Stack Exchange**：[2025 年 12 月 31 日付の Stack Exchange データダンプ](https://archive.org/details/stackexchange_20251231)（CC BY-SA）のうち，Mathematics Educators，Academia，Physics の各サイトの `Posts.xml`．次のように展開してください．

  ```
  data/stackexchange/matheducators/Posts.xml
  data/stackexchange/academia/Posts.xml
  data/stackexchange/physics/Posts.xml
  ```

  各サイトから `tags-<site>`（質問ごとのタグの集合）と `threads-<site>`（質問者と回答者の集合）の 2 つを作ります．変換後の集まり（ハイパーエッジ）の数は次のとおりです．

  | データセット | 集まりの数 |
  |---|---:|
  | tags-matheducators | 3,157 |
  | threads-matheducators | 3,572 |
  | tags-academia | 39,019 |
  | tags-physics | 223,639 |

- **TNA の例示データ**：`tna` 1.2.3 の `group_regulation` と `Nestimate` 0.8.5 の `human_long`．`scripts/r/run_saqr_pipeline.R` が書き出し，`src/tna_to_scholp.py` が各系列を長さ *w* の窓に切ります（*w* より短い系列は一つの集まりにします）．

## 論文の再現

リポジトリのルートで：

```bash
pip install -r requirements.txt

# 同梱の結果から図表だけを作る（数秒）
python3 scripts/run_all.py --steps tables figures

# 生データからすべてを作る
python3 scripts/run_all.py --se-dir data/stackexchange
```

`--steps` で段階を選べます：`se tna check ph null threeway fixed tables figures`．時間がかかるのはヌルモデルで，N = 200 のとき `tags-physics` はノート PC で約 6 時間かかりました（ほかのデータセットは数分）．`--N` でサンプル数を変えられます．乱数の種はすべて 0 です．

| 段階 | 出力 | 用途 |
|---|---|---|
| `se` | `data/ScHoLP-Data/` | Stack Exchange のデータセット |
| `tna` | `results/r_saqr/`，`data/processed/` | TNA のデータセットと `Nestimate` の出力 |
| `check` | （画面に表示） | 付録 A：クリーク複体の構成が `Nestimate::build_simplicial` と一致（32 条件すべてで同一） |
| `ph` | `results/ph_compare_se/`，`results/ph_compare_se_nocap/`，`results/ph_compare_tna/` | 見かけの三角形，最大 β₁，H₁ クラスの行方 |
| `null` | `results/null_real/` | ヌルモデル（N = 200） |
| `threeway` | `results/threeway/` | 三次の交互作用の検定 |
| `fixed` | `results/fixed_threshold/` | 閾値を固定した分析と遷移行列の非対称性 |
| `tables` | `results/tables/` | `table1_construction_null`（構成の比較とヌルモデル），`table2_fixed_threshold`，`table3_threeway` |
| `figures` | `results/figures/` | `fig0_core_idea_{en,ja}`（図 1），`fig_phantom_maxb1`（図 2），`fig1_betti1_null`（図 3），`fig2_h1_classes`（図 4），`fig3_zscores`（図 5） |

付録の `w = 4` の結果は，`results/ph_compare_tna/` と `results/null_real/` の `*-w4_*` のファイルです．`results/null_real/physics_N200/` は別に実行した `tags-physics` のヌルモデルです．三次の交互作用の三つ組ごとのファイル（`*_triples.csv.gz`，最大 100 MB）は含めていません．`run_all.py --steps threeway` で再生成できます．

## ライセンス

コードは [MIT ライセンス](LICENSE) で公開します．Stack Exchange から作ったデータは投稿者による CC BY-SA のもとにあります．変換したデータを共有する場合は，そのライセンスと帰属表示を保ってください．

## 引用

[CITATION.cff](CITATION.cff) を参照してください．arXiv の番号はプレプリントの公開後に追記します．

## 謝辞

本研究は JSPS 科研費 JP23K22315，JP23K17619，JP25K00845，JP25K21951 の助成を受けたものです．
