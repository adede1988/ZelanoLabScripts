## LENS C -- Multivariate within-subject decoding of condition (audiobook vs focus)
## Per subject: stratified 5-fold CV logistic (elastic-net glmnet), robust to p>n.
## All imputation/scaling/lambda-selection INSIDE the training fold (no leakage).
## Within-subject label-permutation null for empirical per-subject p on balanced accuracy.
## Two feature sets: gamma_only, gamma+airflow.
suppressMessages({
  library(dplyr); library(readr); library(stringr); library(glmnet); library(pROC)
})

set.seed(42)
NPERM <- as.integer(Sys.getenv("NPERM", "200"))
KFOLD <- 5
ALPHA <- 0.5   # elastic net

tblDir <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/tables/"
figDir <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/figs/taskcmp/"
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)

pb <- suppressMessages(read_csv(paste0(tblDir, "taskcmp_perbreath.csv")))
ml <- suppressMessages(read_csv(paste0(tblDir, "taskcmp_metric_list.csv")))

airflow      <- ml$metric[ml$family == "airflow_morphology"]
gammaMetrics <- ml$metric[ml$family != "airflow_morphology"]
allMetrics   <- ml$metric
stopifnot(all(allMetrics %in% names(pb)))

fsets <- list(gamma_only = gammaMetrics, gamma_airflow = allMetrics)
subjects <- sort(unique(pb$subject))

## ---- helpers -----------------------------------------------------------------
# stratified k-fold assignment given binary label vector y (0/1)
make_folds <- function(y, k = KFOLD) {
  fold <- integer(length(y))
  for (lev in unique(y)) {
    idx <- which(y == lev)
    idx <- sample(idx)
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

# fixed-direction AUC: predictor = P(focus); focus is the case (y==1)
auc_fixed <- function(prob, y) {
  as.numeric(pROC::auc(response = y, predictor = prob,
                       levels = c(0, 1), direction = "<", quiet = TRUE))
}
# balanced accuracy at 0.5 threshold (mean of sensitivity & specificity)
bacc_fixed <- function(prob, y, thr = 0.5) {
  pred <- as.integer(prob >= thr)
  sens <- mean(pred[y == 1] == 1)   # focus recalled
  spec <- mean(pred[y == 0] == 0)   # audiobook recalled
  mean(c(sens, spec))
}

# Run full 5-fold CV, return pooled out-of-fold P(focus). folds precomputed.
run_cv_oof <- function(X, y, folds, k = KFOLD) {
  n <- length(y)
  oof <- rep(NA_real_, n)
  for (f in seq_len(k)) {
    te <- which(folds == f); tr <- which(folds != f)
    if (length(unique(y[tr])) < 2) { oof[te] <- mean(y[tr]); next }
    Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
    ytr <- y[tr]
    # train-fold medians for imputation
    med <- apply(Xtr, 2, median, na.rm = TRUE)
    for (j in seq_len(ncol(Xtr))) {
      na <- is.na(Xtr[, j]); if (any(na)) Xtr[na, j] <- med[j]
    }
    mu  <- colMeans(Xtr)
    sdv <- apply(Xtr, 2, sd)
    keep <- which(is.finite(mu) & is.finite(sdv) & sdv > 0 & is.finite(med))
    if (length(keep) < 2) { oof[te] <- mean(ytr); next }
    Xtr <- scale(Xtr[, keep, drop = FALSE], center = mu[keep], scale = sdv[keep])
    Xte <- Xte[, keep, drop = FALSE]
    medk <- med[keep]
    for (j in seq_len(ncol(Xte))) {
      na <- is.na(Xte[, j]); if (any(na)) Xte[na, j] <- medk[j]
    }
    Xte <- scale(Xte, center = mu[keep], scale = sdv[keep])
    fit <- tryCatch(
      cv.glmnet(Xtr, ytr, family = "binomial", alpha = ALPHA, nfolds = 5),
      error = function(e) NULL)
    if (is.null(fit)) { oof[te] <- mean(ytr); next }
    oof[te] <- as.numeric(predict(fit, newx = Xte, s = "lambda.min", type = "response"))
  }
  oof
}

## ---- main loop ---------------------------------------------------------------
res <- list()
auc_perm_p_store <- list()  # keep both perm p flavors
t0 <- Sys.time()
for (fs in names(fsets)) {
  feats <- fsets[[fs]]
  for (s in subjects) {
    d <- pb[pb$subject == s, ]
    y <- as.integer(d$condition == "focus")
    X <- as.matrix(d[, feats])
    storage.mode(X) <- "double"

    folds <- make_folds(y)
    oof <- run_cv_oof(X, y, folds)
    obs_auc  <- auc_fixed(oof, y)
    obs_bacc <- bacc_fixed(oof, y)

    # permutation null (shuffle labels within subject, re-stratify, re-CV)
    perm_bacc <- numeric(NPERM); perm_auc <- numeric(NPERM)
    for (p in seq_len(NPERM)) {
      yp <- sample(y)
      fp <- make_folds(yp)
      op <- run_cv_oof(X, yp, fp)
      perm_bacc[p] <- bacc_fixed(op, yp)
      perm_auc[p]  <- auc_fixed(op, yp)
    }
    perm_p_bacc <- (1 + sum(perm_bacc >= obs_bacc)) / (1 + NPERM)
    perm_p_auc  <- (1 + sum(perm_auc  >= obs_auc )) / (1 + NPERM)

    res[[length(res) + 1]] <- data.frame(
      subject = s, feature_set = fs,
      auc = obs_auc, balanced_acc = obs_bacc,
      perm_p = perm_p_bacc, perm_p_auc = perm_p_auc,
      perm_bacc_mean = mean(perm_bacc), perm_auc_mean = mean(perm_auc),
      nBreath = nrow(d), nFocus = sum(y), nAudiobook = sum(y == 0),
      stringsAsFactors = FALSE)
    cat(sprintf("[%s | %-13s] AUC=%.3f bacc=%.3f perm_p(bacc)=%.4f perm_p(auc)=%.4f  (%s)\n",
                s, fs, obs_auc, obs_bacc, perm_p_bacc, perm_p_auc,
                format(Sys.time() - t0, digits = 3)))
    # incremental checkpoint so partial progress survives
    write_csv(bind_rows(res), paste0(tblDir, "taskcmp_C_decoding.csv"))
  }
}
R <- bind_rows(res)

## ---- group-level combination -------------------------------------------------
fisher_p <- function(ps) {
  ps <- pmin(pmax(ps, .Machine$double.eps), 1)
  stat <- -2 * sum(log(ps)); df <- 2 * length(ps)
  pchisq(stat, df, lower.tail = FALSE)
}

grp <- lapply(names(fsets), function(fs) {
  sub <- R[R$feature_set == fs, ]
  # one-sample Wilcoxon signed-rank of AUC vs 0.5 (one-sided greater)
  w_auc  <- suppressWarnings(wilcox.test(sub$auc,          mu = 0.5, alternative = "greater"))
  w_bacc <- suppressWarnings(wilcox.test(sub$balanced_acc, mu = 0.5, alternative = "greater"))
  data.frame(
    feature_set = fs,
    mean_auc = mean(sub$auc), median_auc = median(sub$auc),
    mean_bacc = mean(sub$balanced_acc),
    n_subj_auc_gt_0.5 = sum(sub$auc > 0.5),
    wilcox_auc_p  = w_auc$p.value,
    wilcox_bacc_p = w_bacc$p.value,
    fisher_perm_p_bacc = fisher_p(sub$perm_p),
    fisher_perm_p_auc  = fisher_p(sub$perm_p_auc),
    stringsAsFactors = FALSE)
})
G <- bind_rows(grp)

## ---- write outputs -----------------------------------------------------------
out_csv <- paste0(tblDir, "taskcmp_C_decoding.csv")
write_csv(R, out_csv)
write_csv(G, paste0(tblDir, "taskcmp_C_decoding_group.csv"))

cat("\n===== PER-SUBJECT =====\n"); print(R, row.names = FALSE)
cat("\n===== GROUP =====\n"); print(G, row.names = FALSE)

## ---- plot: per-subject AUC by feature set ------------------------------------
png(paste0(figDir, "C_decoding_auc_by_subject.png"), width = 1100, height = 700, res = 130)
op <- par(mar = c(5, 5, 4, 8), xpd = FALSE)
subs <- subjects
xpos <- seq_along(subs)
cols <- c(gamma_only = "#2166ac", gamma_airflow = "#b2182b")
plot(NA, xlim = c(0.5, length(subs) + 0.5), ylim = c(0.3, 1),
     xaxt = "n", xlab = "subject", ylab = "cross-validated AUC",
     main = "Lens C: within-subject decoding of condition (audiobook vs focus)")
abline(h = 0.5, lty = 2, col = "grey40")
for (fs in names(fsets)) {
  sub <- R[R$feature_set == fs, ]
  sub <- sub[match(subs, sub$subject), ]
  off <- if (fs == "gamma_only") -0.12 else 0.12
  points(xpos + off, sub$auc, pch = 19, col = cols[fs], cex = 1.5)
  # mark permutation-significant (perm_p < .05) with a ring
  sig <- which(sub$perm_p_auc < 0.05)
  if (length(sig)) points((xpos + off)[sig], sub$auc[sig], pch = 1, cex = 2.6, col = cols[fs], lwd = 2)
  segments(xpos + off, R[R$feature_set==fs,][match(subs,R[R$feature_set==fs,]$subject),]$perm_auc_mean,
           xpos + off, sub$auc, col = cols[fs], lwd = 1)
}
# group means
for (fs in names(fsets)) {
  gm <- G$mean_auc[G$feature_set == fs]
  off <- if (fs == "gamma_only") -0.12 else 0.12
  segments(0.5, gm, length(subs) + 0.5, gm, col = cols[fs], lty = 3, lwd = 1.5)
}
axis(1, at = xpos, labels = subs)
par(xpd = TRUE)
legend(length(subs) + 0.7, 0.95, legend = c("gamma_only", "gamma+airflow",
       "chance (0.5)", "perm p<.05 (ring)"),
       col = c(cols["gamma_only"], cols["gamma_airflow"], "grey40", "black"),
       pch = c(19, 19, NA, 1), lty = c(NA, NA, 2, NA), bty = "n", cex = 0.9)
par(op); dev.off()

cat("\nWrote:", out_csv, "\n")
cat("Elapsed:", format(Sys.time() - t0, digits = 4), "\n")
