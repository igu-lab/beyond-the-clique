"""run_all.py — regenerate every figure and table of the paper from the raw data.

Run from the repository root:

    python3 scripts/run_all.py --se-dir data/stackexchange            # everything
    python3 scripts/run_all.py --se-dir data/stackexchange --steps se ph tables figures
    python3 scripts/run_all.py --steps tables figures                  # from the included results only

Steps (in order):
  se        Stack Exchange Posts.xml -> ScHoLP format              data/ScHoLP-Data/
  tna       R: tna / Nestimate example data, build_simplicial        results/r_saqr/
            Python: sequences -> windows (seq, w2, w3, w4)           data/processed/
  check     clique re-implementation vs Nestimate::build_simplicial  (prints, exits 1 on mismatch)
  ph        clique vs Dowker persistent homology                    results/ph_compare_{se,se_nocap,tna}/
  null      degree- and size-preserving null model (slow)           results/null_real/
  threeway  three-way interaction test                              results/threeway/
  fixed     fixed-threshold analysis and asymmetry                  results/fixed_threshold/
  tables    Tables 1-3                                              results/tables/
  figures   Figures 1-5                                             results/figures/
"""
import argparse, os, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
SITES = ["matheducators", "academia", "physics"]
SE_COUNT = ["tags-matheducators", "threads-matheducators", "tags-academia", "tags-physics"]
TNA = [f"{n}-{u}" for n in ("group_regulation", "human_long") for u in ("seq", "w2", "w3", "w4")]
STEPS = ["se", "tna", "check", "ph", "null", "threeway", "fixed", "tables", "figures"]


def sh(*cmd):
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--se-dir", help="directory with <site>/Posts.xml for " + ", ".join(SITES))
    ap.add_argument("--steps", nargs="+", default=STEPS, choices=STEPS)
    ap.add_argument("--N", type=int, default=200, help="null-model samples (paper: 200)")
    a = ap.parse_args()
    scholp, processed, rdir = "data/ScHoLP-Data", "data/processed", "results/r_saqr"

    if "se" in a.steps:
        if not a.se_dir:
            ap.error("--se-dir is required for the 'se' step")
        for s in SITES:
            sh(PY, "src/stackexchange_to_scholp.py", os.path.join(a.se_dir, s, "Posts.xml"), s, "--out", scholp)
    if "tna" in a.steps:
        sh("Rscript", "scripts/r/run_saqr_pipeline.R", rdir)
        for n in ("group_regulation", "human_long"):
            sh(PY, "src/tna_to_scholp.py", f"{rdir}/{n}_long.csv", "--unit", "seq", "w2", "w3", "w4", "--out", processed)
    if "check" in a.steps:
        sh(PY, "scripts/check_nestimate_clique.py", "--r-dir", rdir)
    if "ph" in a.steps:
        sh(PY, "src/ph_compare.py", *SE_COUNT, "--root", scholp, "--maxk", "25", "--out", "results/ph_compare_se")
        sh(PY, "src/ph_compare.py", "threads-matheducators", "--root", scholp, "--maxk", "1000", "--out", "results/ph_compare_se_nocap")
        sh(PY, "src/ph_compare.py", *TNA, "--root", processed, "--maxk", "25", "--out", "results/ph_compare_tna")
    if "null" in a.steps:
        sh(PY, "scripts/null_real.py", "--scholp-root", scholp, "--processed-root", processed, "--out", "results/null_real", "--N", str(a.N))
    if "threeway" in a.steps:
        sh(PY, "scripts/threeway_real.py", "--scholp-root", scholp, "--processed-root", processed, "--out", "results/threeway")
    if "fixed" in a.steps:
        sh(PY, "scripts/fixed_threshold_asym.py", "--r-dir", rdir, "--out", "results/fixed_threshold")
    if "tables" in a.steps:
        sh(PY, "scripts/make_tables.py", "--out", "results/tables")
    if "figures" in a.steps:
        sh(PY, "scripts/make_fig0.py")
        sh(PY, "scripts/make_figures.py")


if __name__ == "__main__":
    main()
