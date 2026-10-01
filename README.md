# Beyond the Clique

[日本語版 README](README.ja.md)

Code for the paper

> Koichi Yasutake, Wakana Tsuji, and Hitoshi Inoue (2026).
> **Beyond the Clique: Comparing Clique and Dowker Complexes for Co-occurrence Data in Learning Analytics.**
> arXiv: *to be added*

Learning analytics often turns co-occurrence data (codes in a window of discourse, participants in a thread, tags on a post) into a simplicial complex and analyses it with persistent homology. The usual choice is the **clique complex**, which is determined by the pairwise network alone and therefore cannot tell "three elements co-occurred together" from "each of the three pairs co-occurred separately". The **Dowker complex** instead takes as simplices the sets that were actually observed together. This repository contains the code that compares the two on four Stack Exchange data sets and on the example data of the R packages `tna` and `Nestimate`, and regenerates every figure and table of the paper.

## Repository layout

| Path | Contents |
|---|---|
| `src/ph_compare.py` | Clique vs Dowker persistent homology (H0, H1) on the same 1-skeleton; pure Python/NumPy |
| `src/null_model.py` | Degree- and size-preserving hypergraph null model |
| `src/threeway_test.py` | Three-way interaction test (log-linear model that keeps all pairwise margins) |
| `src/stackexchange_to_scholp.py` | Stack Exchange `Posts.xml` → ScHoLP format (`nverts` / `simplices` / `times`) |
| `src/tna_to_scholp.py` | TNA sequences (long format) → windows in ScHoLP format |
| `scripts/r/run_saqr_pipeline.R` | Runs `tna` → `Nestimate::build_simplicial` / `persistent_homology` and exports the results |
| `scripts/check_nestimate_clique.py` | Checks that the clique construction here equals `Nestimate::build_simplicial` |
| `scripts/null_real.py`, `scripts/threeway_real.py`, `scripts/fixed_threshold_asym.py` | Analyses on the real data |
| `scripts/make_tables.py`, `scripts/make_figures.py`, `scripts/make_fig0.py` | Tables and figures of the paper |
| `scripts/run_all.py` | Runs everything in order |
| `results/` | Results used in the paper (see below) |

Some code comments are in Japanese.

## Requirements

- Python ≥ 3.10 with `numpy`, `pandas`, `matplotlib`, `seaborn` (`pip install -r requirements.txt`, or `uv sync`)
- For the TNA data only: R with `tna` 1.2.3 and `Nestimate` 0.8.5
- Figure 1 (Japanese version) uses the font Noto Serif CJK JP; the English figures use a Times-like serif font

## Data

The data are **not redistributed** in this repository.

- **Stack Exchange.** `Posts.xml` of the Mathematics Educators, Academia and Physics sites in the [Stack Exchange data dump of 31 December 2025](https://archive.org/details/stackexchange_20251231) (CC BY-SA). Extract them as

  ```
  data/stackexchange/matheducators/Posts.xml
  data/stackexchange/academia/Posts.xml
  data/stackexchange/physics/Posts.xml
  ```

  Each site gives two data sets: `tags-<site>` (one question = its set of tags) and `threads-<site>` (one question = asker ∪ answerers). After conversion the number of groups (hyperedges) should be

  | data set | groups |
  |---|---:|
  | tags-matheducators | 3,157 |
  | threads-matheducators | 3,572 |
  | tags-academia | 39,019 |
  | tags-physics | 223,639 |

- **TNA example data.** `group_regulation` from `tna` 1.2.3 and `human_long` from `Nestimate` 0.8.5. `scripts/r/run_saqr_pipeline.R` exports them; `src/tna_to_scholp.py` cuts each sequence into windows of length *w* (a sequence shorter than *w* forms a single group).

## Reproducing the paper

From the repository root:

```bash
pip install -r requirements.txt

# figures and tables from the included results (seconds)
python3 scripts/run_all.py --steps tables figures

# everything from the raw data
python3 scripts/run_all.py --se-dir data/stackexchange
```

`--steps` selects steps: `se tna check ph null threeway fixed tables figures`. The null model is the slow step: with N = 200 samples, `tags-physics` took about 6 hours on a laptop; the other data sets take minutes. `--N` changes the number of samples. All random steps use seed 0.

| Step | Output | Used for |
|---|---|---|
| `se` | `data/ScHoLP-Data/` | Stack Exchange data sets |
| `tna` | `results/r_saqr/`, `data/processed/` | TNA data sets and the `Nestimate` output |
| `check` | (printed) | Appendix A: the clique construction equals `Nestimate::build_simplicial` (32 / 32 conditions identical) |
| `ph` | `results/ph_compare_se/`, `results/ph_compare_se_nocap/`, `results/ph_compare_tna/` | Phantom triangles, max β₁, fate of H₁ classes |
| `null` | `results/null_real/` | Null model (N = 200) |
| `threeway` | `results/threeway/` | Three-way interaction test |
| `fixed` | `results/fixed_threshold/` | Fixed-threshold analysis and asymmetry of the transition matrices |
| `tables` | `results/tables/` | `table1_construction_null` (construction comparison and null model), `table2_fixed_threshold`, `table3_threeway` |
| `figures` | `results/figures/` | `fig0_core_idea_{en,ja}` (Fig. 1), `fig_phantom_maxb1` (Fig. 2), `fig1_betti1_null` (Fig. 3), `fig2_h1_classes` (Fig. 4), `fig3_zscores` (Fig. 5) |

The `w = 4` results of the appendix are in the `*-w4_*` files of `results/ph_compare_tna/` and `results/null_real/`. `results/null_real/physics_N200/` holds the `tags-physics` null model, which was run separately. The per-triple files of the three-way test (`*_triples.csv.gz`, up to 100 MB) are not included; `run_all.py --steps threeway` regenerates them.

## License

The code is released under the [MIT License](LICENSE). Data derived from Stack Exchange are licensed under CC BY-SA by their contributors; if you share converted data, keep that license and the attribution.

## Citation

See [CITATION.cff](CITATION.cff). The arXiv identifier will be added when the preprint is posted.

## Acknowledgments

This work was supported by JSPS KAKENHI Grant Numbers JP23K22315, JP23K17619, JP25K00845 and JP25K21951.
