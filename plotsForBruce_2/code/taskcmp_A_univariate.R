## LENS A — PAIRED UNIVARIATE across all 93 metrics
## audiobook vs focus, paired within 6 control subjects
suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
})

tabDir <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/tables"
figDir <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/figs/taskcmp"
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)

sm  <- read_csv(file.path(tabDir, "taskcmp_subject_means.csv"), show_col_types = FALSE)
ml  <- read_csv(file.path(tabDir, "taskcmp_metric_list.csv"),   show_col_types = FALSE)

metrics <- ml$metric
famOf   <- setNames(ml$family, ml$metric)
stopifnot(length(metrics) == 93)

subjects <- sort(unique(sm$subject))
cat("subjects:", paste(subjects, collapse=", "), " n=", length(subjects), "\n")

## generic paired test for one aggregate (mean or median column suffix)
run_lens <- function(agg) {
  suffix <- paste0("__", agg)
  res <- map_dfr(metrics, function(m) {
    col <- paste0(m, suffix)
    if (!col %in% names(sm)) return(NULL)
    w <- sm %>% select(subject, condition, val = all_of(col)) %>%
      pivot_wider(names_from = condition, values_from = val) %>%
      arrange(subject)
    # ensure both conditions present for all subjects
    ab <- w$audiobook; fo <- w$focus
    ok <- is.finite(ab) & is.finite(fo)
    ab <- ab[ok]; fo <- fo[ok]
    n  <- length(ab)
    d  <- ab - fo
    md <- mean(d)
    sdd <- sd(d)
    dz <- if (sdd > 0) md / sdd else NA_real_
    # paired t
    tp <- NA_real_; tstat <- NA_real_
    if (n >= 2 && sdd > 0) {
      tt <- tryCatch(t.test(ab, fo, paired = TRUE), error = function(e) NULL)
      if (!is.null(tt)) { tp <- tt$p.value; tstat <- unname(tt$statistic) }
    }
    # wilcoxon signed rank (exact for n=6)
    wp <- NA_real_
    if (n >= 2 && any(d != 0)) {
      wp <- tryCatch(suppressWarnings(wilcox.test(ab, fo, paired = TRUE))$p.value,
                     error = function(e) NA_real_)
    }
    tibble(
      metric        = m,
      family        = famOf[[m]],
      n             = n,
      mean_audiobook = mean(ab),
      mean_focus     = mean(fo),
      mean_diff      = md,
      dz            = dz,
      t_stat        = tstat,
      t_p           = tp,
      wilcox_p      = wp,
      direction     = ifelse(is.na(md) | md == 0, "n/a",
                             ifelse(md > 0, "higher_audiobook", "higher_focus"))
    )
  })
  res$t_p_fdr <- p.adjust(res$t_p, method = "BH")
  res %>% arrange(desc(abs(dz)))
}

meanRes <- run_lens("mean")
medRes  <- run_lens("median")

## primary output = mean-based, with median robustness columns joined
out <- meanRes %>%
  select(metric, family, mean_audiobook, mean_focus, mean_diff, dz,
         t_stat, t_p, t_p_fdr, wilcox_p, direction) %>%
  left_join(
    medRes %>% select(metric,
                      median_dz = dz, median_t_p = t_p,
                      median_t_p_fdr = t_p_fdr, median_wilcox_p = wilcox_p,
                      median_direction = direction),
    by = "metric"
  ) %>%
  arrange(desc(abs(dz)))

write_csv(out, file.path(tabDir, "taskcmp_A_univariate.csv"))

## ---- console summary ----
cat("\n==== LENS A summary (mean columns) ====\n")
nFDR05 <- sum(out$t_p_fdr < 0.05, na.rm = TRUE)
nFDR10 <- sum(out$t_p_fdr < 0.10, na.rm = TRUE)
nRawP05 <- sum(out$t_p < 0.05, na.rm = TRUE)
nBigDz  <- sum(abs(out$dz) > 0.8, na.rm = TRUE)
cat(sprintf("metrics tested: %d\n", nrow(out)))
cat(sprintf("raw t p<0.05: %d  |  FDR<0.10: %d  |  FDR<0.05: %d  |  |dz|>0.8: %d\n",
            nRawP05, nFDR10, nFDR05, nBigDz))

cat("\nTop 15 by |dz|:\n")
print(out %>% mutate(across(where(is.numeric), ~round(.,3))) %>%
        select(metric, family, mean_diff, dz, t_p, t_p_fdr, wilcox_p, direction, median_dz) %>%
        head(15), n = 15)

## agreement between t and wilcoxon/median on the flagged set
flagged <- out %>% filter(t_p_fdr < 0.10 | abs(dz) > 0.8)
cat(sprintf("\nFlagged (FDR<0.10 OR |dz|>0.8): %d metrics\n", nrow(flagged)))
if (nrow(flagged) > 0) {
  print(flagged %>% mutate(across(where(is.numeric), ~round(.,3))) %>%
          select(metric, family, dz, t_p, t_p_fdr, wilcox_p, median_dz, median_wilcox_p, direction))
}

## ---- plot: top 12 by |dz|, paired dots per subject ----
top12 <- head(out$metric, 12)
plotdat <- sm %>%
  select(subject, condition, all_of(paste0(top12, "__mean"))) %>%
  pivot_longer(cols = ends_with("__mean"), names_to = "metric", values_to = "val") %>%
  mutate(metric = str_remove(metric, "__mean$"),
         metric = factor(metric, levels = top12))

png(file.path(figDir, "A_univariate_top12.png"), width = 1500, height = 1100, res = 130)
op <- par(mfrow = c(3,4), mar = c(3,3.2,2.4,0.6), mgp = c(1.9,0.6,0))
for (m in top12) {
  d <- plotdat %>% filter(metric == m)
  wd <- d %>% pivot_wider(names_from = condition, values_from = val)
  ylim <- range(c(wd$audiobook, wd$focus), na.rm = TRUE)
  if (diff(ylim) == 0) ylim <- ylim + c(-1,1)
  plot(NA, xlim = c(0.7, 2.3), ylim = ylim, xaxt = "n",
       xlab = "", ylab = "", main = sprintf("%s\n(dz=%.2f, p=%.3f, q=%.3f)",
       m, out$dz[out$metric==m], out$t_p[out$metric==m], out$t_p_fdr[out$metric==m]),
       cex.main = 0.85)
  axis(1, at = c(1,2), labels = c("audio", "focus"))
  cols <- rainbow(length(subjects))
  for (i in seq_along(subjects)) {
    s <- subjects[i]
    row <- wd %>% filter(subject == s)
    if (nrow(row) == 1) {
      segments(1, row$audiobook, 2, row$focus, col = cols[i], lwd = 1.5)
      points(c(1,2), c(row$audiobook, row$focus), pch = 19, col = cols[i], cex = 0.9)
    }
  }
}
par(op)
dev.off()

cat("\nWrote:", file.path(tabDir, "taskcmp_A_univariate.csv"), "\n")
cat("Wrote:", file.path(figDir, "A_univariate_top12.png"), "\n")
