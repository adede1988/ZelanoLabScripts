# Lens D — respiration-gamma coupling, paired within-control (audiobook vs focusedBreathing)
# n=6 control subjects; AD's two sessions pooled per subject.
suppressWarnings(suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
}))

inCSV  <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/tables/coupling_session.csv"
outCSV <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/tables/taskcmp_D_coupling.csv"
figDir <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/figs/taskcmp"
dir.create(figDir, showWarnings = FALSE, recursive = TRUE)

d <- read_csv(inCSV, show_col_types = FALSE)

# --- filter: control cohort, the two breathing-task conditions
d <- d %>%
  filter(cohort %in% c("OBE","Control"),
         taskRow %in% c("audiobook","focusedBreathing")) %>%
  mutate(condition = ifelse(taskRow == "audiobook", "audiobook", "focus"))

# ---- circular helpers ----
wrap_pi <- function(a) ((a + pi) %% (2*pi)) - pi          # -> [-pi, pi)
circ_mean <- function(a) atan2(mean(sin(a)), mean(cos(a)))
circ_R    <- function(a) sqrt(mean(cos(a))^2 + mean(sin(a))^2)  # mean resultant length

rayleigh_test <- function(a){
  n <- length(a); Rbar <- circ_R(a); Z <- n * Rbar^2
  p <- exp(-Z) * (1 + (2*Z - Z^2)/(4*n) -
                    (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4)/(288*n^2))
  p <- min(max(p, 0), 1)
  list(n=n, Rbar=Rbar, meanDir=circ_mean(a), Z=Z, p=p)
}
# V-test toward mu0 (concentration at a specified direction)
v_test <- function(a, mu0=0){
  n <- length(a); Rbar <- circ_R(a); md <- circ_mean(a)
  Vbar <- Rbar * cos(md - mu0)
  u <- Vbar * sqrt(2*n)
  p <- 1 - pnorm(u)              # one-sided
  list(Vbar=Vbar, u=u, p=p)
}

# ---- pool per subject x condition ----
# scalar measures: arithmetic mean over sessions; phase: circular mean over sessions
scal_measures <- c("coup_MI","coup_resultantLen","coup_inhExhRatio")

subj_scalar <- d %>%
  group_by(participant, condition) %>%
  summarise(across(all_of(scal_measures), ~mean(.x, na.rm=TRUE)),
            nSess = n(), .groups="drop")

subj_phase <- d %>%
  group_by(participant, condition) %>%
  summarise(prefPhase = circ_mean(coup_prefPhaseRad), .groups="drop")

subjects <- sort(unique(d$participant))
cat("Control subjects (n=", length(subjects), "): ", paste(subjects, collapse=", "), "\n", sep="")
# session counts per subject
print(d %>% count(participant, condition) %>% pivot_wider(names_from=condition, values_from=n))

# ---- paired tests on scalar measures ----
mk_wide <- function(df, val){
  df %>% select(participant, condition, !!val) %>%
    pivot_wider(names_from=condition, values_from=!!val) %>%
    arrange(participant)
}

# degeneracy guard: source file assigns one session-level coupling value to BOTH
# audiobook & focusedBreathing rows, so per-subject diffs may be identically 0.
TOL <- 1e-12
res_rows <- list()
for (m in scal_measures){
  w <- mk_wide(subj_scalar, m)
  w <- w %>% filter(!is.na(audiobook) & !is.na(focus))
  a <- w$audiobook; f <- w$focus; diff <- a - f
  degenerate <- max(abs(diff)) < TOL
  if (degenerate){
    dz <- NA_real_; tp <- NA_real_; wp <- NA_real_
    note <- "identical audiobook/focus values in source (no condition contrast)"
    dir <- "n/a"
  } else {
    dz <- mean(diff) / sd(diff)
    tp  <- tryCatch(t.test(a, f, paired=TRUE)$p.value, error=function(e) NA_real_)
    wp  <- tryCatch(suppressWarnings(wilcox.test(a, f, paired=TRUE)$p.value), error=function(e) NA_real_)
    note <- ""
    dir <- ifelse(mean(a) > mean(f), "higher_audiobook", "higher_focus")
  }
  res_rows[[m]] <- tibble(
    measure       = m,
    n             = nrow(w),
    mean_audiobook= mean(a),
    mean_focus    = mean(f),
    mean_diff     = mean(diff),
    dz            = dz,
    t_p           = tp,
    wilcox_p      = wp,
    rayleigh_p    = NA_real_,
    vtest_p       = NA_real_,
    mean_shift_rad= NA_real_,
    resultant_len = NA_real_,
    direction     = dir,
    note          = note
  )
}

# ---- preferred-phase shift (audiobook - focus) ----
wp_ph <- subj_phase %>%
  pivot_wider(names_from=condition, values_from=prefPhase) %>%
  arrange(participant) %>%
  filter(!is.na(audiobook) & !is.na(focus)) %>%
  mutate(shift = wrap_pi(audiobook - focus))
shifts <- wp_ph$shift
phase_degenerate <- max(abs(shifts)) < TOL   # all shifts exactly 0 -> circular tests meaningless
ray <- rayleigh_test(shifts)
vt  <- v_test(shifts, mu0=0)
meanShift <- circ_mean(shifts)
if (phase_degenerate){
  ray_p_rep <- NA_real_; vt_p_rep <- NA_real_
  ph_note <- "identical audiobook/focus prefPhase in source (all shifts=0; Rayleigh/V undefined)"
} else {
  ray_p_rep <- ray$p; vt_p_rep <- vt$p; ph_note <- ""
}

res_rows[["coup_prefPhaseRad"]] <- tibble(
  measure="coup_prefPhaseRad", n=length(shifts),
  mean_audiobook = circ_mean(wp_ph$audiobook),
  mean_focus     = circ_mean(wp_ph$focus),
  mean_diff      = meanShift,
  dz=NA_real_, t_p=NA_real_, wilcox_p=NA_real_,
  rayleigh_p = ray_p_rep, vtest_p = vt_p_rep,
  mean_shift_rad = meanShift, resultant_len = ray$Rbar,
  direction = "n/a", note = ph_note
)

res <- bind_rows(res_rows)

# ---- FDR across primary tests: 3 paired-t on scalar + Rayleigh on phase ----
# (NA where degenerate; p.adjust ignores NA, so only real tests are corrected)
res$p_fdr      <- p.adjust(res$t_p, method="BH")
res$p_fdr[res$measure=="coup_prefPhaseRad"] <- p.adjust(res$rayleigh_p, method="BH")[res$measure=="coup_prefPhaseRad"]
res$wilcox_fdr <- p.adjust(res$wilcox_p, method="BH")

res <- res %>% select(measure, n, mean_audiobook, mean_focus, mean_diff, dz,
                      t_p, p_fdr, wilcox_p, wilcox_fdr, rayleigh_p, vtest_p,
                      mean_shift_rad, resultant_len, direction, note)
write_csv(res, outCSV)
cat("\n==== RESULTS ====\n"); print(as.data.frame(res), digits=4)
cat("\nPer-subject phase shifts (rad):\n")
print(as.data.frame(wp_ph), digits=4)

# ---- plots ----
png(file.path(figDir, "D_coupling_paired.png"), width=1400, height=460, res=120)
op <- par(mfrow=c(1,4), mar=c(4,4,3,1))
pair_plot <- function(m, ttl){
  w <- mk_wide(subj_scalar, m) %>% filter(!is.na(audiobook)&!is.na(focus))
  ylim <- range(c(w$audiobook, w$focus))
  plot(NA, xlim=c(0.7,2.3), ylim=ylim, xaxt="n", xlab="", ylab=m, main=ttl)
  axis(1, at=c(1,2), labels=c("audio","focus"))
  for(i in 1:nrow(w)) lines(c(1,2), c(w$audiobook[i], w$focus[i]), col="grey60")
  points(rep(1,nrow(w)), w$audiobook, pch=19, col="steelblue")
  points(rep(2,nrow(w)), w$focus,     pch=19, col="firebrick")
  segments(0.85,mean(w$audiobook),1.15,mean(w$audiobook),lwd=3,col="steelblue")
  segments(1.85,mean(w$focus),2.15,mean(w$focus),lwd=3,col="firebrick")
}
pair_plot("coup_MI","MI")
pair_plot("coup_resultantLen","resultantLen")
pair_plot("coup_inhExhRatio","inhExhRatio")
# circular plot of phase shifts
ph_main <- if (phase_degenerate) "prefPhase shift\n(all shifts=0: no contrast)" else sprintf("prefPhase shift\nRayleigh p=%.3f", ray$p)
plot(NA, xlim=c(-1.2,1.2), ylim=c(-1.2,1.2), asp=1, xaxt="n", yaxt="n",
     xlab="", ylab="", main=ph_main)
th <- seq(0,2*pi,length.out=200); lines(cos(th), sin(th), col="grey70")
segments(0,0, cos(0),0, col="grey85"); segments(0,0,0,sin(0),col="grey85")
if (phase_degenerate){
  points(0,0,pch=19,col="darkorange",cex=1.4)
  text(0,-0.5,"6 subjects\nall at origin\n(shift=0)",col="darkorange",cex=0.9)
} else {
  for(s in shifts) arrows(0,0, cos(s), sin(s), length=0.08, col="darkorange")
  arrows(0,0, ray$Rbar*cos(meanShift), ray$Rbar*sin(meanShift), length=0.14, lwd=3, col="black")
}
par(op); dev.off()
cat("\nWrote:", outCSV, "\nWrote fig:", file.path(figDir,"D_coupling_paired.png"), "\n")
