# make_review_packet.R — assemble a concise markdown packet of methods+results for reviewers.
suppressWarnings(suppressMessages({library(dplyr); library(readr); library(stringr); library(knitr)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
T <- file.path(proj,"out","tables")
rd <- function(f){ p<-file.path(T,f); if(file.exists(p)) suppressMessages(read_csv(p,show_col_types=FALSE)) else NULL }
o <- c()
add <- function(...) o <<- c(o, sprintf(...))
tbl <- function(df,n=15){ if(is.null(df)) return("_(missing)_"); paste(capture.output(print(kable(head(as.data.frame(df),n)))), collapse="\n") }

add("# Dupi Olfactory Gamma — Review Packet\n")
add("## Design\n")
add("Cohorts: Dupi intervention (S1/2/3) vs OBE controls; EEG excluded (no nasal electrode). Best macBP channel per recording = max peak flattened power (FOOOF) in 25-58 Hz. Superlet TFR (fractional adaptive, c1=3, order [3 30]), validated to 4e-14 vs reference faslt. Ridge via Frequency_ridge_tracking forward-backward tracker (ported to MATLAB) on per-breath within-frequency z (myChanZscore, no baseline). Per-breath FOOOF-lite (aperiodic exponent/offset + gamma peak presence). Session-level respiration-gamma coupling (Tort MI, preferred phase, resultant length, inhale/exhale ratio). Group stats aggregate to session level first.\n")

add("## Inventory\n")
idx <- rd("session_index.csv")
if(!is.null(idx)){ s <- idx |> filter(cohort %in% c("Dupi","OBE"), onDisk==1) |> count(task,cohort)
  add("On-disk finals:\n\n%s\n", tbl(s,20)) }
add("Missing (sheet->disk): %s\n", { m<-rd("INVENTORY_missing_finals.csv"); if(is.null(m)||nrow(m)==0) "none" else paste(m$sessID, collapse=", ") })

add("\n## Behavioral & responder split\n")
add("Per-session composite means by X-position:\n\n")
ss <- rd("session_scores_labeled.csv")
if(!is.null(ss)){
  bx <- ss |> group_by(xpos) |> summarise(n=sum(is.finite(composite)), mean_comp=round(mean(composite,na.rm=TRUE),3))
  add("%s\n", tbl(bx,10))
}
add("\nResponder table:\n\n%s\n", tbl(rd("responder_table.csv"),20))

add("\n## macBP gamma selection\n")
b <- rd("macbp_best.csv")
if(!is.null(b)){
  add("Selected channel had FOOOF-detectable gamma peak in %d/%d recordings (%.0f%%). Mean non-selected channels with a detectable peak: %.2f.\n",
      sum(b$bestGammaDetected), nrow(b), 100*mean(b$bestGammaDetected), mean(b$nNonMaxDetected))
  add("Selected-channel spikeFrac summary (noise proxy): median %.4f, >0.02 in %d recordings.\n",
      median(b$bestSpikeFrac,na.rm=TRUE), sum(b$bestSpikeFrac>0.02,na.rm=TRUE))
}

add("\n## Group spectrograms\n")
add("5 task-rows (cueTask, threshTask, O15, audiobook, focusedBreathing) x 5 columns (control, S2/3 responder, S1 responder, S2/3 non-responder, S1 non-responder). See figs/spectrograms_5x5.png. Group map = mean over breaths of per-breath single-trial within-frequency z (no baseline).\n")

add("\n## Gamma measures — goodness ranking (KEY RESULT)\n")
add("Goodness = |max separation effect| / (1 + |control CV|). Top measures:\n\n")
G <- rd("gamma_top_measures.csv")
if(!is.null(G)){
  disp <- G |> transmute(task,metric,control_CV=round(control_CV,3),
    sep_control_dupi=round(sep_control_dupi,2), sep_resp_nonresp=round(sep_resp_nonresp,2),
    recovery_rho_resp=round(recovery_rho_resp,2), goodness=round(goodness,2)) |> arrange(desc(goodness))
  add("%s\n", tbl(disp,30))
}
add("\n## Coupling summary\n")
add("Session-level respiration-gamma coupling per recording in out/gamma/coupling/*.csv (coup_MI, coup_prefPhaseRad, coup_resultantLen, coup_rayleighP, coup_inhExhRatio).\n")

writeLines(o, file.path(proj,"out","reviews","review_packet.md"))
cat("wrote out/reviews/review_packet.md\n")
