#!/usr/bin/env Rscript
# run_saqr_pipeline.R — Saqr らの手順（tna → Nestimate::build_simplicial / persistent_homology）を
# そのまま走らせ、Python 側（Dowker 複体）との突き合わせに必要な出力を CSV で書き出す。
#
# 実行（プロジェクトルートで）:
#   Rscript scripts/r/run_saqr_pipeline.R
# または RStudio で「Source」。出力先は results/r_saqr_YYYY-MM-DD_HHMM/（既存は上書きしない）。
#   出力先を指定する場合: Rscript scripts/r/run_saqr_pipeline.R results/r_saqr
#
# 書き出すもの（データセットごとに接頭辞 <name>_）:
#   <name>_long.csv              actor, time, action の長形式（Python 側で窓を切るための原データ）
#   <name>_W_relative.csv        tna の遷移行列（行=from, 列=to, 相対頻度）
#   <name>_W_relative_sym.csv    pmax(W, t(W))（Nestimate が内部で使う対称化）
#   <name>_W_cooc_w<k>.csv       Nestimate::build_network(method="co_occurrence", window_size=k) の重み
#   <name>_sc_<mat>_t<thr>_simplices.csv  build_simplicial の全単体（dim, nodes を "|" 区切り）
#   <name>_sc_<mat>_t<thr>_summary.csv    f-vector, Betti, Euler
#   <name>_ph_<mat>_betti_curve.csv       persistent_homology の Betti 曲線
#   <name>_ph_<mat>_persistence.csv       birth/death 対
#   sessionInfo.txt, log.txt

suppressPackageStartupMessages({
  library(tna)
  library(Nestimate)
})

stamp   <- format(Sys.time(), "%Y-%m-%d_%H%M")
if (!file.exists("pyproject.toml")) stop("プロジェクトルート（pyproject.toml のある場所）で実行してください: setwd(\"<beyond_the_clique>\")")
args    <- commandArgs(trailingOnly = TRUE)   # 2026-10-01: 出力先を引数で指定できるようにした（省略時は日時付き）
outdir  <- if (length(args) >= 1) args[1] else file.path("results", paste0("r_saqr_", stamp))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
logf    <- file.path(outdir, "log.txt")
logmsg  <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S "), sprintf(...)); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
wcsv    <- function(x, name) { p <- file.path(outdir, name); write.csv(x, p, row.names = FALSE); logmsg("wrote %s (%d rows)", name, NROW(x)) }
wmat    <- function(M, name) { p <- file.path(outdir, name); write.csv(as.data.frame(as.matrix(M)), p, row.names = TRUE); logmsg("wrote %s (%dx%d)", name, nrow(M), ncol(M)) }

## ---- 単体複体オブジェクトの汎用書き出し -------------------------------------
dump_sc <- function(sc, prefix) {
  saveRDS(sc, file.path(outdir, paste0(prefix, ".rds")))
  writeLines(capture.output(print(sc)), file.path(outdir, paste0(prefix, "_print.txt")))
  writeLines(capture.output(str(sc, max.level = 2)), file.path(outdir, paste0(prefix, "_str.txt")))
  # 単体リスト：sc$simplices が list（各要素がノード名/添字のベクトル）である想定。違えば str を見て直す。
  simp <- NULL
  for (nm in c("simplices", "simplex_list", "cliques")) if (!is.null(sc[[nm]])) { simp <- sc[[nm]]; break }
  if (is.list(simp)) {
    if (!is.null(names(simp)) && all(grepl("^(dim|d)?[0-9]+$", names(simp)))) simp <- unlist(simp, recursive = FALSE)
    df <- data.frame(
      dim   = vapply(simp, function(s) length(s) - 1L, integer(1)),
      nodes = vapply(simp, function(s) paste(s, collapse = "|"), character(1)),
      stringsAsFactors = FALSE)
    wcsv(df, paste0(prefix, "_simplices.csv"))
  } else {
    logmsg("  [%s] simplices not found as list; see _str.txt", prefix)
  }
  summ <- list()
  for (nm in c("f_vector", "fvector", "f", "betti", "betti_numbers", "euler", "euler_characteristic", "dimension", "n_simplices", "threshold"))
    if (!is.null(sc[[nm]]) && is.atomic(sc[[nm]])) summ[[nm]] <- paste(sc[[nm]], collapse = "|")
  if (length(summ)) wcsv(data.frame(key = names(summ), value = unlist(summ)), paste0(prefix, "_summary.csv"))
}

dump_ph <- function(ph, prefix) {
  saveRDS(ph, file.path(outdir, paste0(prefix, ".rds")))
  writeLines(capture.output(print(ph)), file.path(outdir, paste0(prefix, "_print.txt")))
  if (!is.null(ph$betti_curve)) wcsv(ph$betti_curve, paste0(prefix, "_betti_curve.csv"))
  if (!is.null(ph$persistence)) wcsv(ph$persistence, paste0(prefix, "_persistence.csv"))
  if (!is.null(ph$thresholds))  wcsv(data.frame(threshold = ph$thresholds), paste0(prefix, "_thresholds.csv"))
}

## ---- 1 データセット分の処理 -------------------------------------------------
run_one <- function(name, long, thresholds = c(0.02, 0.05, 0.10, 0.20), windows = c(2L, 3L, 4L)) {
  logmsg("==== %s: %d rows, %d actors, %d codes", name, nrow(long), length(unique(long$actor)), length(unique(long$action)))
  wcsv(long, paste0(name, "_long.csv"))

  ## (a) tna の遷移行列（Saqr らの標準の入口）
  wide <- tryCatch({
    sp <- split(long$action, long$actor)
    L  <- max(lengths(sp))
    as.data.frame(do.call(rbind, lapply(sp, function(s) c(s, rep(NA, L - length(s))))), stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  model <- tna::tna(wide)
  W <- model$weights
  wmat(W, paste0(name, "_W_relative.csv"))
  wmat(pmax(W, t(W)), paste0(name, "_W_relative_sym.csv"))

  mats <- list(relative = W)

  ## (b) Nestimate の共起ネットワーク（窓幅ごと）— Dowker の窓と同じ対象
  for (k in windows) {
    Wc <- tryCatch({
      n <- Nestimate::build_network(long, method = "co_occurrence", actor = "actor", action = "action",
                                    time = "time", window_size = k, mode = "non-overlapping")
      n$weights
    }, error = function(e) { logmsg("  co_occurrence w=%d failed: %s", k, conditionMessage(e)); NULL })
    if (!is.null(Wc)) { wmat(Wc, sprintf("%s_W_cooc_w%d.csv", name, k)); mats[[sprintf("cooc_w%d", k)]] <- Wc }
  }

  ## (c) 各行列について build_simplicial（固定閾値）と persistent_homology（閾値掃引）
  for (mn in names(mats)) {
    M <- mats[[mn]]
    for (thr in thresholds) {
      sc <- tryCatch(Nestimate::build_simplicial(M, type = "clique", threshold = thr, max_dim = 3L),
                     error = function(e) { logmsg("  build_simplicial %s t=%.3f failed: %s", mn, thr, conditionMessage(e)); NULL })
      if (!is.null(sc)) dump_sc(sc, sprintf("%s_sc_%s_t%.3f", name, mn, thr))
    }
    ph <- tryCatch(Nestimate::persistent_homology(M, n_steps = 40L, max_dim = 2L, type = "clique"),
                   error = function(e) { logmsg("  persistent_homology %s failed: %s", mn, conditionMessage(e)); NULL })
    if (!is.null(ph)) dump_ph(ph, sprintf("%s_ph_%s", name, mn))
  }
}

## ---- データセット ------------------------------------------------------------

# 1. tna::group_regulation（wide → long）
gr <- tna::group_regulation
gr_long <- do.call(rbind, lapply(seq_len(nrow(gr)), function(i) {
  s <- as.character(unlist(gr[i, ])); s <- s[!is.na(s) & s != ""]
  if (!length(s)) return(NULL)
  data.frame(actor = paste0("g", i), time = seq_along(s), action = s, stringsAsFactors = FALSE)
}))
tryCatch(run_one("group_regulation", gr_long), error = function(e) logmsg("group_regulation failed: %s", conditionMessage(e)))

# 2. Nestimate::human_long（チュートリアルで使われている長形式データ）
hl <- tryCatch({ data("human_long", package = "Nestimate"); get("human_long") }, error = function(e) NULL)
if (!is.null(hl)) {
  logmsg("human_long columns: %s", paste(names(hl), collapse = ", "))
  ac <- intersect(c("session_id", "actor", "session", "id"), names(hl))[1]
  cd <- intersect(c("code", "action", "state"), names(hl))[1]
  tm <- intersect(c("timestamp", "time", "t"), names(hl))[1]
  hl_long <- data.frame(actor = as.character(hl[[ac]]), time = hl[[tm]], action = as.character(hl[[cd]]), stringsAsFactors = FALSE)
  hl_long <- hl_long[order(hl_long$actor, hl_long$time), ]
  tryCatch(run_one("human_long", hl_long), error = function(e) logmsg("human_long failed: %s", conditionMessage(e)))
} else logmsg("human_long not found in Nestimate; skipped")

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
logmsg("done -> %s", outdir)
