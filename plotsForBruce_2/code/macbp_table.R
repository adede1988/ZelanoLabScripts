# macbp_table.R — per-recording macBP 25-58 Hz peak flattened power table + FOOOF detectability.
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(stringr)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
tdir <- file.path(proj,"out","tables")

g <- read_csv(file.path(tdir,"macbp_gamma_long.csv"), show_col_types=FALSE)
b <- read_csv(file.path(tdir,"macbp_best.csv"), show_col_types=FALSE)

# wide: peak flattened power (dB over aperiodic) per macBP channel, one row per recording
wide <- g %>% select(sessID, task, cohort, group, participant, sessNum, chanLabel, peakFlatDb) %>%
  pivot_wider(names_from=chanLabel, values_from=peakFlatDb) %>%
  arrange(cohort, participant, sessNum, task)

# attach best-channel + detectability info
wide2 <- wide %>% left_join(
  b %>% select(sessID, task, bestLabel, bestPeakFlatDb, bestGammaDetected,
               bestSpikeFrac, nNonMaxDetected, nMacBPDetected, nMacBP),
  by=c("sessID","task"))
write_csv(wide2, file.path(tdir,"macbp_table_wide.csv"))

# detectability among non-max channels: also make a per-recording detail
det <- g %>% group_by(sessID, task) %>%
  summarise(nMacBP=n(),
            nDetected=sum(gammaDetected),
            .groups="drop")

# ================= summaries the prompt asks for =================
nrec <- nrow(b)
cat("=== macBP gamma selection: ", nrec, " recordings (session x task) ===\n\n", sep="")

cat("--- Was the SELECTED (max peak-flattened) channel a FOOOF-detectable gamma peak? ---\n")
print(b %>% count(bestGammaDetected) %>%
        mutate(pct=round(100*n/sum(n),1)) %>% as.data.frame(), row.names=FALSE)

cat("\n--- Non-max channels with a detectable FOOOF gamma peak (per recording) ---\n")
print(b %>% count(nNonMaxDetected) %>% as.data.frame(), row.names=FALSE)
cat(sprintf("  mean non-max detected per recording: %.2f (of mean %.1f non-max channels)\n",
            mean(b$nNonMaxDetected), mean(b$nMacBP-1)))

cat("\n--- by task: fraction of recordings where selected channel had detectable gamma ---\n")
print(b %>% group_by(task) %>%
        summarise(n=n(), selDetected=sum(bestGammaDetected),
                  pctSelDetected=round(100*mean(bestGammaDetected),1),
                  meanNonMaxDet=round(mean(nNonMaxDetected),2)) %>% as.data.frame(),
      row.names=FALSE)

cat("\n--- selected-channel noisiness (spikeFrac) distribution ---\n")
print(summary(b$bestSpikeFrac))
cat("  recordings where selected channel spikeFrac>0.02 (possibly noise-driven):",
    sum(b$bestSpikeFrac>0.02, na.rm=TRUE), "\n")

cat("\nWrote macbp_table_wide.csv\n")
