# =====================================================================
# LENS E — AIRFLOW-CONTROLLED GAMMA
# Within-control (n=6) AUDIOBOOK vs FOCUSED breathing, paired.
# (1) Which airflow_morphology metrics differ (paired t, subject unit)
# (2) For each GAMMA metric: does the condition effect SURVIVE airflow
#     control? per-breath LMM metric ~ condition (+/- airflow) + (1|subject)
# =====================================================================
suppressWarnings(suppressMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
  library(lme4);  library(lmerTest); library(ggplot2)
}))

tblDir <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/tables"
figDir <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/figs/taskcmp"
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)

pb   <- read_csv(file.path(tblDir, "taskcmp_perbreath.csv"), show_col_types = FALSE)
sm   <- read_csv(file.path(tblDir, "taskcmp_subject_means.csv"), show_col_types = FALSE)
mlst <- read_csv(file.path(tblDir, "taskcmp_metric_list.csv"), show_col_types = FALSE)

airflow_metrics <- mlst$metric[mlst$family == "airflow_morphology"]
gamma_metrics   <- mlst$metric[mlst$family != "airflow_morphology"]
covars          <- c("inhaleVolume","inhaleDuration","peakInspFlow","breathLength")
fam_of <- setNames(mlst$family, mlst$metric)

cond_lv <- c("audiobook","focus")               # audiobook = reference
pb$condition <- factor(pb$condition, levels = cond_lv)

# ---------------------------------------------------------------------
# PART 1 — BREATHING MANIPULATION: paired t on airflow metrics (subject unit, n=6)
# uses subject_means (AD already pooled to one subject/condition)
# ---------------------------------------------------------------------
manip <- map_dfr(airflow_metrics, function(m){
  col <- paste0(m, "__mean")
  w <- sm %>% select(subject, condition, val = all_of(col)) %>%
    pivot_wider(names_from = condition, values_from = val) %>%
    filter(!is.na(audiobook) & !is.na(focus))
  d  <- w$focus - w$audiobook                    # focus - audiobook
  n  <- nrow(w)
  tt <- tryCatch(t.test(w$focus, w$audiobook, paired = TRUE), error = function(e) NULL)
  dz <- mean(d) / sd(d)                          # Cohen's dz (paired)
  data.frame(metric = m, family = "airflow_morphology", n = n,
             mean_audiobook = mean(w$audiobook), mean_focus = mean(w$focus),
             mean_diff = mean(d),
             t = if(is.null(tt)) NA else unname(tt$statistic),
             df = if(is.null(tt)) NA else unname(tt$parameter),
             p_raw = if(is.null(tt)) NA else tt$p.value,
             dz = dz,
             direction = ifelse(mean(d) > 0, "higher_focus", "higher_audiobook"),
             stringsAsFactors = FALSE)
})
manip$p_fdr <- p.adjust(manip$p_raw, method = "BH")
manip <- manip %>% arrange(p_raw)
write_csv(manip, file.path(tblDir, "taskcmp_E_airflow_manipulation.csv"))
cat("=== PART 1: BREATHING MANIPULATION (airflow metrics, paired t, n=6) ===\n")
print(manip %>% mutate(across(where(is.numeric), ~round(.,4))), row.names = FALSE)

# ---------------------------------------------------------------------
# PART 2 — GAMMA metrics: condition effect BEFORE vs AFTER airflow control
# per-breath LMM.  raw : metric ~ condition + (1|subject)
#                  adj : metric ~ condition + z(airflow covars) + (1|subject)
# both fit on the SAME complete-case rows for a fair before/after compare.
# ---------------------------------------------------------------------
# global z-score of covariates (helps convergence; does not change cond effect)
pb <- pb %>% mutate(across(all_of(covars), ~ as.numeric(scale(.)), .names = "{.col}_z"))
covz <- paste0(covars, "_z")

fit_one <- function(m){
  keep <- c("subject","condition", m, covz)
  d <- pb[, keep]
  names(d)[3] <- "y"
  d <- d[stats::complete.cases(d), , drop = FALSE]
  out <- data.frame(metric = m, family = unname(fam_of[m]),
                    n_breath = nrow(d), n_subj = NA_integer_,
                    binary = FALSE, method = NA_character_,
                    cond_beta_raw = NA, cond_se_raw = NA, cond_p_raw = NA,
                    cond_beta_adj = NA, cond_se_adj = NA, cond_p_adj = NA,
                    cond_df_adj = NA,                    # Satterthwaite df (exposes pseudorep)
                    cond_p_adj_rs = NA, rs_singular = NA, # random-SLOPE adj model (honest LMM)
                    cond_beta_adj_subj = NA, cond_p_adj_subj = NA, # subject-level paired t (n=6)
                    sd_y = NA, stringsAsFactors = FALSE)
  if (nrow(d) < 10 || length(unique(d$condition)) < 2 ||
      length(unique(na.omit(d$y))) < 2) return(out)
  out$n_subj <- length(unique(d$subject))
  out$sd_y   <- sd(d$y, na.rm = TRUE)
  out$binary <- length(unique(na.omit(d$y))) <= 2

  # both conditions must be present within >=2 subjects for a paired-ish test
  perSubj <- d %>% group_by(subject) %>% summarise(k = n_distinct(condition), .groups="drop")
  bothN   <- sum(perSubj$k == 2)

  get_cond <- function(fit){
    co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(co)) return(c(NA,NA,NA,NA))
    rn <- grep("^condition", rownames(co), value = TRUE)
    if (length(rn) == 0) return(c(NA,NA,NA,NA))
    est <- co[rn[1], "Estimate"]
    se  <- co[rn[1], "Std. Error"]
    pc  <- if ("Pr(>|t|)" %in% colnames(co)) co[rn[1], "Pr(>|t|)"] else NA
    df  <- if ("df"       %in% colnames(co)) co[rn[1], "df"]       else NA
    c(est, se, pc, df)
  }

  # subject-level paired t on airflow-RESIDUALIZED means (the honest n=6 test)
  subj_paired_adj <- function(dd){
    ols <- tryCatch(lm(y ~ inhaleVolume_z + inhaleDuration_z + peakInspFlow_z +
                         breathLength_z, data = dd), error = function(e) NULL)
    if (is.null(ols)) return(c(NA,NA))
    dd$res <- resid(ols)
    w <- dd %>% group_by(subject, condition) %>%
      summarise(mv = mean(res), .groups="drop") %>%
      pivot_wider(names_from = condition, values_from = mv) %>%
      filter(!is.na(audiobook) & !is.na(focus))
    if (nrow(w) < 3) return(c(NA,NA))
    da <- w$focus - w$audiobook
    tt <- tryCatch(t.test(w$focus, w$audiobook, paired=TRUE), error=function(e) NULL)
    c(mean(da), if(is.null(tt)) NA else tt$p.value)
  }

  lmer_ok <- (bothN >= 2 && out$n_subj >= 3)
  if (lmer_ok) {
    f0 <- tryCatch(suppressWarnings(suppressMessages(
            lmerTest::lmer(y ~ condition + (1|subject), data = d, REML = TRUE))),
            error = function(e) NULL)
    f1 <- tryCatch(suppressWarnings(suppressMessages(
            lmerTest::lmer(y ~ condition + inhaleVolume_z + inhaleDuration_z +
                             peakInspFlow_z + breathLength_z + (1|subject),
                           data = d, REML = TRUE))),
            error = function(e) NULL)
    if (!is.null(f0) && !is.null(f1)) {
      out$method <- "lmer"
      r0 <- get_cond(f0); r1 <- get_cond(f1)
      out$cond_beta_raw <- r0[1]; out$cond_se_raw <- r0[2]; out$cond_p_raw <- r0[3]
      out$cond_beta_adj <- r1[1]; out$cond_se_adj <- r1[2]; out$cond_p_adj <- r1[3]
      out$cond_df_adj   <- r1[4]
      # honest robustness: random-SLOPE adj model (condition effect varies by subject)
      frs <- tryCatch(suppressWarnings(suppressMessages(
               lmerTest::lmer(y ~ condition + inhaleVolume_z + inhaleDuration_z +
                                peakInspFlow_z + breathLength_z + (1 + condition|subject),
                              data = d, REML = TRUE))),
               error = function(e) NULL)
      if (!is.null(frs)) {
        rr <- get_cond(frs); out$cond_p_adj_rs <- rr[3]
        out$rs_singular <- isSingular(frs, tol = 1e-4)
      }
      # honest robustness: subject-level paired t on residualized means (n=6)
      sp <- subj_paired_adj(d)
      out$cond_beta_adj_subj <- sp[1]; out$cond_p_adj_subj <- sp[2]
      return(out)
    }
  }

  # ---- FALLBACK: per-subject residualize-then-paired-t ----
  # subject-level condition means; adj = residualize y on covars (pooled OLS) first
  out$method <- "resid_pairedt"
  # RAW subject means paired t
  sm_raw <- d %>% group_by(subject, condition) %>%
    summarise(mv = mean(y), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = mv) %>%
    filter(!is.na(audiobook) & !is.na(focus))
  if (nrow(sm_raw) >= 2) {
    dr <- sm_raw$focus - sm_raw$audiobook
    tt <- tryCatch(t.test(sm_raw$focus, sm_raw$audiobook, paired = TRUE),
                   error = function(e) NULL)
    out$cond_beta_raw <- mean(dr)
    out$cond_se_raw   <- sd(dr)/sqrt(length(dr))
    out$cond_p_raw    <- if (is.null(tt)) NA else tt$p.value
  }
  # ADJ: residualize y on covariates, then subject means paired t
  dd <- d
  ols <- tryCatch(lm(y ~ inhaleVolume_z + inhaleDuration_z + peakInspFlow_z +
                       breathLength_z, data = dd), error = function(e) NULL)
  if (!is.null(ols)) {
    dd$res <- resid(ols)
    sm_adj <- dd %>% group_by(subject, condition) %>%
      summarise(mv = mean(res), .groups = "drop") %>%
      pivot_wider(names_from = condition, values_from = mv) %>%
      filter(!is.na(audiobook) & !is.na(focus))
    if (nrow(sm_adj) >= 2) {
      da <- sm_adj$focus - sm_adj$audiobook
      tt <- tryCatch(t.test(sm_adj$focus, sm_adj$audiobook, paired = TRUE),
                     error = function(e) NULL)
      out$cond_beta_adj <- mean(da)
      out$cond_se_adj   <- sd(da)/sqrt(length(da))
      out$cond_p_adj    <- if (is.null(tt)) NA else tt$p.value
    }
  }
  out
}

res <- map_dfr(gamma_metrics, fit_one)

# FDR (BH) across gamma metrics
res$cond_p_raw_fdr      <- p.adjust(res$cond_p_raw,      method = "BH")
res$cond_p_adj_fdr      <- p.adjust(res$cond_p_adj,      method = "BH")
res$cond_p_adj_rs_fdr   <- p.adjust(res$cond_p_adj_rs,   method = "BH")
res$cond_p_adj_subj_fdr <- p.adjust(res$cond_p_adj_subj, method = "BH")

# attenuation: condition effect present raw but gone after airflow control
res$airflow_attenuated <- with(res, !is.na(cond_p_raw) & !is.na(cond_p_adj) &
                                 cond_p_raw < 0.05 & cond_p_adj >= 0.05)
# survives: still significant after airflow control
res$survives_adj_raw05 <- with(res, !is.na(cond_p_adj) & cond_p_adj < 0.05)
res$survives_adj_fdr   <- with(res, !is.na(cond_p_adj_fdr) & cond_p_adj_fdr < 0.05)
res$direction_adj <- ifelse(is.na(res$cond_beta_adj), "n/a",
                       ifelse(res$cond_beta_adj > 0, "higher_focus", "higher_audiobook"))
# standardized condition effect (adj) for cross-metric comparison
res$cond_d_adj <- res$cond_beta_adj / res$sd_y

res <- res %>% arrange(cond_p_adj)

# honest survivor: significant under airflow-adjusted random-intercept LMM (FDR),
# random-SLOPE LMM (uncorrected .05), AND subject-level paired-t (uncorrected .05)
res$survives_honest <- with(res,
  survives_adj_fdr &
  !is.na(cond_p_adj_rs)   & cond_p_adj_rs   < 0.05 &
  !is.na(cond_p_adj_subj) & cond_p_adj_subj < 0.05)

# ---- write the spec'd table (plus extras) ----
outCols <- c("metric","family","method","n_breath","n_subj","binary",
             "cond_beta_raw","cond_p_raw","cond_p_raw_fdr",
             "cond_beta_adj","cond_p_adj","cond_p_adj_fdr","cond_df_adj",
             "cond_p_adj_rs","cond_p_adj_rs_fdr","rs_singular",
             "cond_beta_adj_subj","cond_p_adj_subj","cond_p_adj_subj_fdr",
             "cond_d_adj","direction_adj",
             "airflow_attenuated","survives_adj_raw05","survives_adj_fdr",
             "survives_honest")
write_csv(res[, outCols], file.path(tblDir, "taskcmp_E_airflow.csv"))

# ---------------------------------------------------------------------
# SUMMARY to stdout
# ---------------------------------------------------------------------
cat("\n=== PART 2: GAMMA condition effect, before vs after airflow control ===\n")
cat(sprintf("gamma metrics tested: %d  (method lmer: %d, fallback: %d)\n",
            nrow(res), sum(res$method=="lmer", na.rm=TRUE),
            sum(res$method=="resid_pairedt", na.rm=TRUE)))
cat(sprintf("raw   p<.05 (uncorr): %d ;  FDR<.05: %d\n",
            sum(res$cond_p_raw<0.05, na.rm=TRUE), sum(res$cond_p_raw_fdr<0.05, na.rm=TRUE)))
cat(sprintf("ADJ   p<.05 (uncorr): %d ;  FDR<.05: %d   <-- survive airflow control\n",
            sum(res$survives_adj_raw05, na.rm=TRUE), sum(res$survives_adj_fdr, na.rm=TRUE)))
cat(sprintf("attenuated by airflow (raw sig -> adj n.s.): %d\n", sum(res$airflow_attenuated, na.rm=TRUE)))
cat("\n-- HONEST triangulation of the airflow-ADJUSTED condition effect --\n")
cat(sprintf("  random-intercept LMM (df~breaths, OPTIMISTIC): FDR<.05 = %d\n", sum(res$survives_adj_fdr, na.rm=TRUE)))
cat(sprintf("  random-SLOPE LMM (condition|subject, honest) : uncorr<.05 = %d ; FDR<.05 = %d ; singular fits = %d\n",
            sum(res$cond_p_adj_rs<0.05, na.rm=TRUE), sum(res$cond_p_adj_rs_fdr<0.05, na.rm=TRUE),
            sum(res$rs_singular, na.rm=TRUE)))
cat(sprintf("  subject-level paired-t (n=6, most honest) : uncorr<.05 = %d ; FDR<.05 = %d\n",
            sum(res$cond_p_adj_subj<0.05, na.rm=TRUE), sum(res$cond_p_adj_subj_fdr<0.05, na.rm=TRUE)))
cat(sprintf("  survives ALL THREE (FDR-LMM & slope<.05 & subj-t<.05): %d\n", sum(res$survives_honest, na.rm=TRUE)))
cat(sprintf("  median Satterthwaite df of condition (adj): %.0f  (n_subj=6 -> anything >>5 = pseudoreplication)\n",
            median(res$cond_df_adj, na.rm=TRUE)))

cat("\n-- top 20 gamma metrics by ADJ p (after airflow control) --\n")
show <- res %>% transmute(metric, family, method,
                          beta_raw=round(cond_beta_raw,3), p_raw=round(cond_p_raw,4),
                          beta_adj=round(cond_beta_adj,3), p_adj=round(cond_p_adj,4),
                          p_adj_fdr=round(cond_p_adj_fdr,4), d_adj=round(cond_d_adj,3),
                          dir=direction_adj, atten=airflow_attenuated,
                          surv_fdr=survives_adj_fdr) %>% head(20)
print(as.data.frame(show), row.names = FALSE)

cat("\n-- metrics surviving ALL THREE honest lenses (the defensible gamma differences) --\n")
hs <- res %>% filter(survives_honest) %>%
  transmute(metric, family, beta_adj=round(cond_beta_adj,3), d_adj=round(cond_d_adj,3),
            p_adj_fdr=round(cond_p_adj_fdr,4), p_slope=round(cond_p_adj_rs,4),
            p_subjt=round(cond_p_adj_subj,4), dir=direction_adj)
if (nrow(hs)==0) cat("  (none)\n") else print(as.data.frame(hs), row.names = FALSE)

# ---------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------
# Plot A: breathing manipulation (paired subject dots per airflow metric, standardized)
smA <- sm %>% select(subject, condition, ends_with("__mean"))
mlong <- map_dfr(airflow_metrics, function(m){
  col <- paste0(m,"__mean")
  x <- smA[[col]]
  z <- (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)   # z within metric for comparable axis
  data.frame(subject = smA$subject, condition = smA$condition, metric = m, z = z)
})
mlong$metric <- factor(mlong$metric, levels = manip$metric)  # order by significance
mlong$sig <- manip$p_fdr[match(mlong$metric, manip$metric)] < 0.05
pA <- ggplot(mlong, aes(condition, z, group = subject)) +
  geom_line(alpha=.35, colour="grey55") +
  geom_point(aes(colour=condition), size=2) +
  facet_wrap(~metric, nrow=2) +
  scale_colour_manual(values=c(audiobook="#4477AA", focus="#EE6677")) +
  labs(title="Lens E Part 1: Breathing manipulation (airflow metrics, n=6 paired)",
       subtitle="z-scored within metric; lines = subjects. Bold facet strips would be FDR-sig.",
       y="z (within metric)", x=NULL) +
  theme_bw(base_size=10) + theme(legend.position="none")
ggsave(file.path(figDir,"E_airflow_manipulation.png"), pA, width=11, height=5, dpi=130)

# Plot B: gamma condition effect before vs after airflow control
pdat <- res %>% filter(!is.na(cond_p_raw) & !is.na(cond_p_adj)) %>%
  mutate(nlp_raw = -log10(cond_p_raw), nlp_adj = -log10(cond_p_adj),
         status = case_when(survives_adj_fdr ~ "survives (FDR)",
                            survives_adj_raw05 ~ "survives (raw p<.05)",
                            airflow_attenuated ~ "attenuated by airflow",
                            TRUE ~ "n.s. either way"))
lim <- max(c(pdat$nlp_raw, pdat$nlp_adj, -log10(0.05)), na.rm=TRUE)*1.05
pB <- ggplot(pdat, aes(nlp_raw, nlp_adj, colour=status)) +
  geom_abline(slope=1, intercept=0, colour="grey70", linetype=2) +
  geom_hline(yintercept=-log10(0.05), colour="grey60", linetype=3) +
  geom_vline(xintercept=-log10(0.05), colour="grey60", linetype=3) +
  geom_point(aes(shape=family), size=2.4, alpha=.85) +
  scale_shape_manual(values=c(16,17,15,3,7,8,4,5)) +
  scale_colour_manual(values=c("survives (FDR)"="#228833","survives (raw p<.05)"="#66CCEE",
                               "attenuated by airflow"="#EE6677","n.s. either way"="grey65")) +
  coord_equal(xlim=c(0,lim), ylim=c(0,lim)) +
  labs(title="Lens E Part 2: gamma condition effect before vs after airflow control",
       subtitle="each point = one gamma metric; per-breath LMM (1|subject). Dotted = p=.05.",
       x="-log10 p  (condition, NO airflow control)",
       y="-log10 p  (condition, WITH airflow control)") +
  theme_bw(base_size=10)
ggsave(file.path(figDir,"E_gamma_beforeafter.png"), pB, width=8.5, height=6.5, dpi=130)

cat("\nWROTE:\n",
    file.path(tblDir,"taskcmp_E_airflow.csv"), "\n",
    file.path(tblDir,"taskcmp_E_airflow_manipulation.csv"), "\n",
    file.path(figDir,"E_airflow_manipulation.png"), "\n",
    file.path(figDir,"E_gamma_beforeafter.png"), "\n")
