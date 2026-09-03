# inventory_report.R — participant x session x task grid + missing flags.
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(stringr)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
tdir <- file.path(proj,"out","tables")
idx <- read_csv(file.path(tdir,"session_index.csv"), show_col_types=FALSE)

taskAbbr <- c(breathingTask="B", cueTask="C", threshTask="T", O15="O")
idx <- idx %>% mutate(ab = taskAbbr[task])

# ---- 1. Missing finals (sheet-listed but no preproc on disk) ----
missing <- idx %>% filter(missing==1) %>% select(task, sessID, cohort, group, participant, sessNum)
write_csv(missing, file.path(tdir,"INVENTORY_missing_finals.csv"))

# ---- 2. Dupi participant x session x task grid ----
dupi <- idx %>% filter(cohort=="Dupi", onDisk==1)
# collapse tasks present per (participant, sessNum)
grid <- dupi %>% group_by(participant, sessNum) %>%
  summarise(tasks = paste(sort(ab), collapse=""), nTask=n(), .groups="drop") %>%
  arrange(participant, sessNum)
# wide: one row per participant, columns S1/S2/S3 = task string
wide <- grid %>% mutate(scol=paste0("S",sessNum)) %>%
  select(participant, scol, tasks) %>%
  pivot_wider(names_from=scol, values_from=tasks, values_fill="")
# ensure S1..S3 cols
for (s in c("S1","S2","S3")) if (!s %in% names(wide)) wide[[s]] <- ""
wide <- wide %>% select(participant, S1, S2, S3)
nsess <- dupi %>% distinct(participant, sessNum) %>% count(participant, name="nSessions")
wide <- wide %>% left_join(nsess, by="participant") %>% arrange(participant)
write_csv(wide, file.path(tdir,"INVENTORY_dupi_grid.csv"))

# per-task Dupi session counts by session number
dupiCounts <- dupi %>% count(task, sessNum) %>% arrange(task, sessNum)
write_csv(dupiCounts, file.path(tdir,"INVENTORY_dupi_taskcounts.csv"))

# ---- 3. Control (OBE) inventory ----
obe <- idx %>% filter(cohort=="OBE", onDisk==1)
obeGrid <- obe %>% group_by(sessID, participant) %>%
  summarise(tasks=paste(sort(ab), collapse=""), nTask=n(), .groups="drop") %>% arrange(sessID)
write_csv(obeGrid, file.path(tdir,"INVENTORY_control_grid.csv"))

obeTaskCounts <- obe %>% count(task, name="nControl")

# ---- 4. EEG (excluded) note ----
eeg <- idx %>% filter(cohort=="EEG", onDisk==1)
eegN <- eeg %>% count(task, name="nEEG_excluded")

# ---- print summary ----
cat("=========== INVENTORY SUMMARY ===========\n\n")
cat("Analysis cohorts: Dupi (intervention, S1/2/3) + OBE (control). EEG EXCLUDED (no nasal electrode).\n\n")

cat("--- Missing finals (in sheet, not on disk) ---\n")
if (nrow(missing)==0) cat("  none\n") else print(as.data.frame(missing))
cat("\n")

cat("--- Dupi participants:", n_distinct(dupi$participant), "---\n")
print(as.data.frame(wide), row.names=FALSE)
cat("\n(letters: B=breathing C=cue T=thresh O=O15)\n\n")

cat("--- Dupi finals per task x session ---\n")
print(dupiCounts %>% pivot_wider(names_from=sessNum, values_from=n, values_fill=0) %>% as.data.frame(), row.names=FALSE)
cat("\n")

cat("--- OBE control finals per task ---\n")
print(as.data.frame(obeTaskCounts), row.names=FALSE)
cat("  distinct OBE control sessions:", n_distinct(obe$sessID), "\n\n")

cat("--- EEG (excluded) per task ---\n")
print(as.data.frame(eegN), row.names=FALSE)

cat("\nWrote INVENTORY_*.csv to", tdir, "\n")
