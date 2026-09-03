# inventory_reverse.R — session FOLDERS on disk that are NOT in the sheet-driven index.
# Catches raw folders present on disk but not marked raw-extracted (so excluded by applyParams).
suppressWarnings(suppressMessages({library(dplyr); library(readr); library(stringr)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
tdir <- file.path(proj,"out","tables")
idx <- read_csv(file.path(tdir,"session_index.csv"), show_col_types=FALSE)
sheetIDs <- unique(idx$sessID)

roots <- c(Dupi="R:/Neurology/Zelano_Lab/Lab_Common/Dupi",
           OBE ="R:/Neurology/Zelano_Lab/Lab_Common/OBEControl")
sufMap <- c(breathing="_breathingPreProc.mat", cue="_cueTaskPreproc.mat",
            thresh="_PEA_threshold_preproc.mat", O15="_O15preproc.mat")

rows <- list()
for (cn in names(roots)) {
  root <- roots[[cn]]
  dirs <- list.dirs(root, recursive=FALSE, full.names=FALSE)
  # session-like: start with 6 digits + _ + cohort token
  sess <- dirs[str_detect(dirs, "^\\d{6}_(Dupi|OBE)_")]
  for (id in sess) {
    inSheet <- id %in% sheetIDs
    pdir <- file.path(root, id, "preProc")
    finals <- c()
    if (dir.exists(pdir)) {
      for (tk in names(sufMap)) {
        f <- file.path(pdir, paste0(id, sufMap[[tk]]))
        if (file.exists(f)) finals <- c(finals, tk)
      }
    }
    rows[[length(rows)+1]] <- tibble(cohort=cn, sessID=id, inSheet=inSheet,
                                     hasPreProcDir=dir.exists(pdir),
                                     finalsOnDisk=paste(finals, collapse=","),
                                     nFinals=length(finals))
  }
}
allDisk <- bind_rows(rows)
notInSheet <- allDisk %>% filter(!inSheet) %>% arrange(cohort, sessID)
write_csv(notInSheet, file.path(tdir,"INVENTORY_ondisk_not_in_sheet.csv"))
write_csv(allDisk,    file.path(tdir,"INVENTORY_all_disk_folders.csv"))

cat("=== Session folders ON DISK but NOT in sheet-driven analysis set ===\n")
cat("(present as raw folders; excluded because Raw Data Extracted blank/INCOMPLETE, or no finals)\n\n")
if (nrow(notInSheet)==0) cat("  none\n") else print(as.data.frame(notInSheet), row.names=FALSE)
cat("\nTotal disk session folders scanned:", nrow(allDisk),
    "| in sheet:", sum(allDisk$inSheet), "| not in sheet:", nrow(notInSheet), "\n")
