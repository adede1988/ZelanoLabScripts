# behavioral_composite.R — per-session olfactory scores, composite index, responder split, plots.
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
tdir <- file.path(proj,"out","tables"); fdir <- file.path(proj,"out","figs")
dir.create(fdir, showWarnings=FALSE, recursive=TRUE)

beh <- read_csv(file.path(tdir,"behavioral_long.csv"), show_col_types=FALSE)
idx <- read_csv(file.path(tdir,"session_index.csv"), show_col_types=FALSE) %>%
  filter(cohort %in% c("Dupi","OBE"), onDisk==1) %>%
  distinct(sessID, cohort, group, participant, sessNum)

# wide per session: one column per metric (metrics are unique per task, so no collision)
w <- beh %>% select(sessID, metric, value) %>%
  pivot_wider(names_from=metric, values_from=value, values_fn=~mean(.x, na.rm=TRUE))
sess <- idx %>% left_join(w, by="sessID")
# ensure expected metric columns exist (robust to partial runs / absent tasks)
needcols <- c("cue_d","cue_HR","cue_FA","cue_hits","cue_misses","cue_fas","cue_crs",
              "thresh_none","thresh_low","thresh_high","thresh_low_cal","thresh_high_cal",
              "O15_acc","O15_score")
for(nm in needcols) if(!nm %in% names(sess)) sess[[nm]] <- NA_real_

# ---- composite (Olfaction Index): mean of present & non-zero of the 3 normalized subscores
sess <- sess %>% mutate(
  cs_O15    = O15_acc,
  cs_thresh = thresh_low_cal/400,
  cs_cue    = cue_d/3.5,
  nAll3     = (!is.na(cue_d)) + (!is.na(thresh_low_cal)) + (!is.na(O15_acc)),
  cs_n      = (!is.na(cs_O15)    & cs_O15    !=0) +
              (!is.na(cs_thresh) & cs_thresh !=0) +
              (!is.na(cs_cue)    & cs_cue    !=0),
  composite = ifelse(cs_n>0,
                     (coalesce(cs_O15,0)+coalesce(cs_thresh,0)+coalesce(cs_cue,0))/cs_n, NA_real_),
  composite_all3 = ifelse(nAll3==3, composite, NA_real_)  # composite requiring all 3 tasks present
)
write_csv(sess %>% select(-starts_with("cs_")), file.path(tdir,"session_scores.csv"))

# ---- xpos factor: session1/2/3/control ----
sess <- sess %>% mutate(
  xpos = ifelse(cohort=="OBE","Control", paste0("S",sessNum)),
  xpos = factor(xpos, levels=c("S1","S2","S3","Control"))
)

# ================= responder / non-responder split (Dupi, composite S2 vs S1) =================
dupiComp <- sess %>% filter(cohort=="Dupi") %>%
  select(participant, sessNum, composite, composite_all3)
compWide <- dupiComp %>%
  pivot_wider(names_from=sessNum, values_from=c(composite, composite_all3), names_sep="_S")
# use all-3 composite when available at both S1 & S2; else fall back to partial composite
getpair <- function(df){
  s1a <- df[["composite_all3_S1"]]; s2a <- df[["composite_all3_S2"]]
  s1p <- df[["composite_S1"]];      s2p <- df[["composite_S2"]]
  s1 <- ifelse(!is.na(s1a) & !is.na(s2a), s1a, s1p)
  s2 <- ifelse(!is.na(s1a) & !is.na(s2a), s2a, s2p)
  list(s1=s1, s2=s2, usedAll3 = (!is.na(s1a) & !is.na(s2a)))
}
pr <- getpair(compWide)
resp <- compWide %>% transmute(
  participant,
  composite_S1 = pr$s1, composite_S2 = pr$s2, usedAll3 = pr$usedAll3,
  delta = composite_S2 - composite_S1,
  class = case_when(
    is.na(composite_S1) | is.na(composite_S2) ~ "unclassified (missing S1 or S2 composite)",
    delta > 0 ~ "responder",
    TRUE ~ "non-responder")
) %>% arrange(desc(delta))
write_csv(resp, file.path(tdir,"responder_table.csv"))

# map class back onto sessions (for coloring). Controls -> "control"; unclassified Dupi -> "unclassified"
classMap <- resp %>% select(participant, class)
sess <- sess %>% left_join(classMap, by="participant") %>%
  mutate(grpColor = case_when(
    cohort=="OBE" ~ "control",
    class=="responder" ~ "responder",
    class=="non-responder" ~ "non-responder",
    TRUE ~ "unclassified"))
write_csv(sess %>% select(sessID,cohort,group,participant,sessNum,xpos,grpColor,
                          cue_d,cue_HR,cue_FA,thresh_none,thresh_low,thresh_high,
                          thresh_low_cal,thresh_high_cal,O15_score,O15_acc,composite,composite_all3),
          file.path(tdir,"session_scores_labeled.csv"))

# ================= plots: performance by session/control, per task =================
pal <- c(responder="#1b9e77", `non-responder`="#d95f02", control="#7570b3", unclassified="#999999")
plotMetric <- function(df, ycol, ylab, fname){
  d <- df %>% filter(!is.na(.data[[ycol]]))
  if(nrow(d)==0) return(invisible())
  ms <- d %>% group_by(xpos) %>%
    summarise(m=mean(.data[[ycol]],na.rm=TRUE),
              se=sd(.data[[ycol]],na.rm=TRUE)/sqrt(sum(!is.na(.data[[ycol]]))),
              n=sum(!is.na(.data[[ycol]])), .groups="drop")
  p <- ggplot(d, aes(x=xpos, y=.data[[ycol]])) +
    geom_hline(yintercept=0, color="grey85", linewidth=.3) +
    geom_jitter(aes(color=grpColor), width=.12, height=0, size=2, alpha=.75) +
    geom_errorbar(data=ms, aes(x=xpos, ymin=m-se, ymax=m+se), inherit.aes=FALSE, width=.16, linewidth=.5) +
    geom_point(data=ms, aes(x=xpos, y=m), inherit.aes=FALSE, size=3.2, shape=95, color="black") +
    geom_text(data=ms, aes(x=xpos, y=m, label=paste0("n=",n)), inherit.aes=FALSE,
              vjust=-1.1, size=3, color="grey30") +
    scale_color_manual(values=pal, name="group") +
    labs(x=NULL, y=ylab, title=ylab) + theme_bw(base_size=12) +
    theme(legend.position="bottom")
  ggsave(file.path(fdir,fname), p, width=6, height=4.2, dpi=130)
}
plotMetric(sess, "cue_d",           "cueTask d'",                        "beh_cue_d.png")
plotMetric(sess, "thresh_low_cal",  "threshTask low-PEA (air-calibrated)","beh_thresh_lowcal.png")
plotMetric(sess, "thresh_high_cal", "threshTask high-PEA (air-calibrated)","beh_thresh_highcal.png")
plotMetric(sess, "O15_acc",         "O15 identification accuracy",        "beh_o15_acc.png")
plotMetric(sess, "composite",       "Olfaction composite index",          "beh_composite.png")
plotMetric(sess, "composite_all3",  "Olfaction composite (all-3 tasks)",  "beh_composite_all3.png")

# ---- console summary ----
cat("=== Responder / non-responder split (Dupi, composite S2 vs S1) ===\n")
print(as.data.frame(resp), row.names=FALSE)
cat("\n=== per-xpos composite means ===\n")
print(sess %>% group_by(xpos) %>%
        summarise(n=sum(!is.na(composite)), mean_composite=mean(composite,na.rm=TRUE)) %>% as.data.frame(),
      row.names=FALSE)
cat("\nWrote session_scores.csv, session_scores_labeled.csv, responder_table.csv and beh_*.png\n")
