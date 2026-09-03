## LENS B - Distributional / higher-moment differences between conditions
## audiobook vs focus, paired within-subject (n=6 control subjects; AD pooled)
## For each subject x metric: log variance ratio log(var_focus/var_audiobook),
##   two-sample KS statistic, IQR difference. Aggregate across subjects:
##   Wilcoxon signed-rank of per-subject logVarRatio vs 0, paired t (more powerful
##   with n=6), sign test, and a KS permutation test combined via Stouffer.
## BH-FDR across gamma metrics (airflow_morphology = the breathing manipulation,
##   reported as positive controls, excluded from the gamma FDR family).

suppressWarnings(suppressMessages({
  library(readr); library(dplyr); library(stringr)
}))
set.seed(12345)

base   <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
tabdir <- file.path(base, "out/tables")
figdir <- file.path(base, "out/figs/taskcmp")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

d    <- read.csv(file.path(tabdir, "taskcmp_perbreath.csv"), stringsAsFactors = FALSE)
mets <- read.csv(file.path(tabdir, "taskcmp_metric_list.csv"), stringsAsFactors = FALSE)

idcols  <- c("subject","sessID","sessNum","condition","breathIdx")
metcols <- setdiff(names(d), idcols)
stopifnot(all(metcols %in% mets$metric))
famOf   <- setNames(mets$family, mets$metric)

subs   <- sort(unique(d$subject))
minN   <- 8          # min non-NA breaths per condition per subject
NPERM  <- 2000       # KS permutation reps

## ---- fast two-sample KS statistic (no ties correction; consistent under perm) ----
ksstat <- function(x, y){
  n <- length(x); m <- length(y)
  v <- c(x, y); o <- order(v)
  d <- cumsum(c(rep(1/n, n), rep(-1/m, m))[o])
  max(abs(d))
}

## ---- per subject x metric distributional contrasts ----
rows <- list()
for (mc in metcols){
  for (s in subs){
    sub <- d[d$subject == s, ]
    ab  <- sub[[mc]][sub$condition == "audiobook"]
    fo  <- sub[[mc]][sub$condition == "focus"]
    ab  <- ab[is.finite(ab)]; fo <- fo[is.finite(fo)]
    if (length(ab) < minN || length(fo) < minN) next
    va <- var(ab); vf <- var(fo)
    lvr <- if (is.finite(va) && is.finite(vf) && va > 0 && vf > 0) log(vf/va) else NA_real_
    ks  <- ksstat(ab, fo)                               # raw KS (location+scale+shape)
    abc <- ab - median(ab); foc <- fo - median(fo)      # median-centered -> spread/shape only
    ksc <- ksstat(abc, foc)
    iqrd <- IQR(fo, na.rm = TRUE) - IQR(ab, na.rm = TRUE)
    ## KS permutation p (label shuffle within subject), one-sided (KS larger = more diff)
    pool <- c(ab, fo); poolc <- c(abc, foc)
    N <- length(pool); n1 <- length(ab)
    perm <- numeric(NPERM); permc <- numeric(NPERM)
    for (p in seq_len(NPERM)){
      idx <- sample.int(N, n1)
      perm[p]  <- ksstat(pool[idx],  pool[-idx])
      permc[p] <- ksstat(poolc[idx], poolc[-idx])
    }
    ks_p  <- (sum(perm  >= ks ) + 1) / (NPERM + 1)
    ksc_p <- (sum(permc >= ksc) + 1) / (NPERM + 1)
    rows[[length(rows)+1]] <- data.frame(
      metric = mc, subject = s,
      nAB = length(ab), nFO = length(fo),
      logVarRatio = lvr, KS = ks, ks_p = ks_p,
      KS_centered = ksc, ksc_p = ksc_p, IQRdiff = iqrd,
      stringsAsFactors = FALSE)
  }
}
persub <- do.call(rbind, rows)
write_csv(persub, file.path(tabdir, "taskcmp_B_distribution_persubject.csv"))

## ---- Stouffer combine of one-sided KS perm p-values across subjects ----
stouffer <- function(pv){
  pv <- pv[is.finite(pv)]
  k  <- length(pv); if (k < 3) return(NA_real_)
  pv <- pmin(pmax(pv, 1e-6), 1 - 1e-6)
  z  <- qnorm(1 - pv)
  1 - pnorm(sum(z) / sqrt(k))
}

## ---- aggregate across subjects per metric ----
agg <- lapply(metcols, function(mc){
  sub <- persub[persub$metric == mc, ]
  lvr <- sub$logVarRatio[is.finite(sub$logVarRatio)]
  nsu <- length(lvr)
  wp <- tp <- sp <- NA_real_
  if (nsu >= 4){
    wp <- suppressWarnings(wilcox.test(lvr, mu = 0, exact = FALSE)$p.value)
    tp <- tryCatch(t.test(lvr, mu = 0)$p.value, error = function(e) NA_real_)
    npos <- sum(lvr > 0); nneg <- sum(lvr < 0); nnz <- npos + nneg
    sp <- if (nnz > 0) binom.test(max(npos, nneg), nnz, 0.5)$p.value else NA_real_
  }
  data.frame(
    metric = mc,
    family = famOf[[mc]],
    n_subj = nsu,
    mean_logVarRatio = if (nsu > 0) mean(lvr) else NA_real_,
    sd_logVarRatio   = if (nsu > 1) sd(lvr) else NA_real_,
    n_pos = sum(lvr > 0), n_neg = sum(lvr < 0),
    varRatio_p       = wp,           # Wilcoxon signed-rank (per spec)
    varRatio_t_p     = tp,           # paired t on logVarRatio (more powerful at n=6)
    varRatio_sign_p  = sp,           # sign test
    mean_KS          = mean(sub$KS),
    mean_KS_centered = mean(sub$KS_centered),
    mean_IQRdiff     = mean(sub$IQRdiff),
    KS_combined_p          = stouffer(sub$ks_p),
    KS_centered_combined_p = stouffer(sub$ksc_p),
    stringsAsFactors = FALSE)
})
agg <- do.call(rbind, agg)

agg$is_manipulation <- agg$family == "airflow_morphology"
agg$direction <- ifelse(is.na(agg$mean_logVarRatio), "n/a",
                 ifelse(agg$mean_logVarRatio > 0, "higher_focus", "higher_audiobook"))

## ---- FDR within the gamma family (exclude airflow manipulation metrics) ----
gmask <- !agg$is_manipulation
agg$varRatio_p_fdr             <- NA_real_
agg$varRatio_t_p_fdr           <- NA_real_
agg$KS_combined_p_fdr          <- NA_real_
agg$KS_centered_combined_p_fdr <- NA_real_
agg$varRatio_p_fdr[gmask]             <- p.adjust(agg$varRatio_p[gmask],             "BH")
agg$varRatio_t_p_fdr[gmask]           <- p.adjust(agg$varRatio_t_p[gmask],           "BH")
agg$KS_combined_p_fdr[gmask]          <- p.adjust(agg$KS_combined_p[gmask],          "BH")
agg$KS_centered_combined_p_fdr[gmask] <- p.adjust(agg$KS_centered_combined_p[gmask], "BH")

## order: gamma metrics by paired-t p first, then airflow controls
agg <- agg[order(agg$is_manipulation, agg$varRatio_t_p), ]

outcols <- c("metric","family","is_manipulation","n_subj","n_pos","n_neg",
             "mean_logVarRatio","sd_logVarRatio","direction",
             "varRatio_p","varRatio_p_fdr",
             "varRatio_t_p","varRatio_t_p_fdr",
             "varRatio_sign_p",
             "mean_KS","mean_KS_centered","mean_IQRdiff",
             "KS_combined_p","KS_combined_p_fdr",
             "KS_centered_combined_p","KS_centered_combined_p_fdr")
write_csv(agg[, outcols], file.path(tabdir, "taskcmp_B_distribution.csv"))

## ---- console summary ----
cat("=== LENS B summary ===\n")
cat(sprintf("metrics analysed: %d (gamma %d, airflow-manip %d)\n",
            nrow(agg), sum(gmask), sum(!gmask)))
cat(sprintf("min Wilcoxon signed-rank p possible at n=6: %.5f (2/2^6)\n", 2/2^6))
gam <- agg[!agg$is_manipulation, ]
cat("\n-- Wilcoxon signed-rank (varRatio_p) FDR<0.10 in gamma set:\n")
print(gam[which(gam$varRatio_p_fdr < 0.10),
          c("metric","family","mean_logVarRatio","varRatio_p","varRatio_p_fdr","direction")])
cat("\n-- Paired-t on logVarRatio, FDR<0.10 in gamma set:\n")
print(gam[which(gam$varRatio_t_p_fdr < 0.10),
          c("metric","family","mean_logVarRatio","varRatio_t_p","varRatio_t_p_fdr","direction")])
cat("\n-- KS permutation (raw, Stouffer) FDR<0.10 in gamma set:\n")
kr <- gam[which(gam$KS_combined_p_fdr < 0.10), ]
kr <- kr[order(kr$KS_combined_p_fdr), ]
print(kr[, c("metric","family","mean_KS","KS_combined_p","KS_combined_p_fdr","direction")])
cat(sprintf("  (%d gamma metrics survive raw-KS FDR<0.10)\n", nrow(kr)))
cat("\n-- CENTERED KS permutation (spread/shape only, Stouffer) FDR<0.10 in gamma set:\n")
kc <- gam[which(gam$KS_centered_combined_p_fdr < 0.10), ]
kc <- kc[order(kc$KS_centered_combined_p_fdr), ]
print(kc[, c("metric","family","mean_KS_centered","KS_centered_combined_p","KS_centered_combined_p_fdr","direction")])
cat(sprintf("  (%d gamma metrics survive CENTERED-KS FDR<0.10 = genuine higher-moment diffs)\n", nrow(kc)))
cat("\n-- Airflow manipulation positive controls (raw p, not in gamma FDR):\n")
print(agg[agg$is_manipulation,
          c("metric","mean_logVarRatio","varRatio_p","varRatio_t_p","mean_KS","direction")])
cat("\n-- Top 15 gamma metrics by |mean_logVarRatio|:\n")
gtop <- gam[order(-abs(gam$mean_logVarRatio)), ]
print(head(gtop[,c("metric","family","mean_logVarRatio","n_pos","n_neg",
                    "varRatio_p","varRatio_t_p","mean_KS","direction")], 15))

## ---- PLOT: top spread differences (per-subject logVarRatio dots + mean) ----
plot_top <- function(rank_by, fname, ttl, nshow = 16){
  g <- gam[order(rank_by[match(gam$metric, agg$metric[!agg$is_manipulation])]), ]
  ## simpler: order gam by supplied vector aligned to gam rows
  invisible(NULL)
}
## rank gamma metrics by paired-t p (fallback mean|lvr|); show per-subject points
gr <- gam
ordv <- gr$varRatio_t_p; ordv[is.na(ordv)] <- 1
gr <- gr[order(ordv), ]
nshow <- min(16, nrow(gr))
sel <- gr$metric[seq_len(nshow)]

png(file.path(figdir, "B_top_logVarRatio.png"), width = 1100, height = 900, res = 120)
op <- par(mar = c(5, 9, 4, 2))
yl <- length(sel)
plot(NA, xlim = range(persub$logVarRatio[is.finite(persub$logVarRatio)], na.rm = TRUE),
     ylim = c(0.5, yl + 0.5), yaxt = "n", xlab = "log(var_focus / var_audiobook)",
     ylab = "", main = ttl <- "Lens B: per-subject spread ratio (top gamma metrics)")
abline(v = 0, col = "gray50", lwd = 2)
for (i in seq_len(nshow)){
  m  <- sel[i]; yy <- yl - i + 1
  pv <- persub$logVarRatio[persub$metric == m]
  pv <- pv[is.finite(pv)]
  points(pv, rep(yy, length(pv)), pch = 21, bg = "#4C78A8", col = "white", cex = 1.3)
  mm <- mean(pv)
  points(mm, yy, pch = 18, col = "#E4572E", cex = 1.9)
  tp <- gr$varRatio_t_p[i]
  axis(2, at = yy, labels = sprintf("%s", m), las = 1, cex.axis = 0.8)
  text(par("usr")[2], yy, sprintf("t p=%.3f", tp), pos = 2, cex = 0.7, col = "gray30")
}
legend("topleft", legend = c("subject", "mean"), pch = c(21, 18),
       pt.bg = c("#4C78A8", NA), col = c("white", "#E4572E"), pt.cex = c(1.3, 1.9), bg = "white")
par(op); dev.off()

## PLOT 2: mean KS by metric (all metrics, colored by family), top annotated
png(file.path(figdir, "B_meanKS_by_family.png"), width = 1200, height = 700, res = 120)
op <- par(mar = c(5, 5, 4, 2))
ga <- agg[!agg$is_manipulation, ]
ga <- ga[order(-ga$mean_KS), ]
cols <- as.integer(factor(ga$family))
plot(ga$mean_KS, type = "h", col = "#B0B0B0",
     xlab = "gamma metric (ranked by mean KS)", ylab = "mean two-sample KS statistic",
     main = "Lens B: mean per-subject KS distance (audiobook vs focus)")
points(seq_len(nrow(ga)), ga$mean_KS, pch = 19, col = rainbow(max(cols))[cols], cex = 0.8)
topk <- head(seq_len(nrow(ga)), 8)
text(topk, ga$mean_KS[topk], ga$metric[topk], pos = 4, cex = 0.65, srt = 0)
legend("topright", legend = levels(factor(ga$family)),
       col = rainbow(max(cols)), pch = 19, cex = 0.8, bg = "white")
par(op); dev.off()

cat("\nWrote:\n",
    file.path(tabdir, "taskcmp_B_distribution.csv"), "\n",
    file.path(tabdir, "taskcmp_B_distribution_persubject.csv"), "\n",
    file.path(figdir, "B_top_logVarRatio.png"), "\n",
    file.path(figdir, "B_meanKS_by_family.png"), "\n")
