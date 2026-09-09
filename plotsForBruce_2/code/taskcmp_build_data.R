# taskcmp_build_data.R — assemble the control-only audiobook-vs-focusedBreathing
# per-breath dataset used by the task-comparison sweep. Writes:
#   out/tables/taskcmp_perbreath.csv       (long: one row per breath, both conditions)
#   out/tables/taskcmp_subject_means.csv   (subject x condition mean/median/sd/n per metric)
#   out/tables/taskcmp_metric_list.csv     (metric name + family)
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
pbdir <- file.path(proj,"out","gamma","perbreath"); tdir <- file.path(proj,"out","tables")

read_cond <- function(cond){
  fs <- list.files(pbdir, pattern=paste0("OBE.*__",cond,"\\.csv$"), full.names=TRUE)
  fs |> map(function(f){
    d <- suppressMessages(read_csv(f, show_col_types=FALSE))
    d$condition <- ifelse(cond=="audiobook","audiobook","focus")
    d
  }) |> list_rbind()
}
pb <- bind_rows(read_cond("audiobook"), read_cond("focusedBreathing"))
# QC: keep good breaths (breathing has goodBreath; treat NA as good)
pb <- pb |> filter(is.na(goodBreath) | goodBreath==1)
# subject id already in `participant`; AD has two sessions (AD_1, AD_2) -> same participant
pb <- pb |> mutate(subject=participant)

id_cols <- c("sessID","task","cohort","group","participant","sessNum","breathIdx",
             "onsetSample","goodBreath","nEv","condition","subject")
metric_cols <- setdiff(names(pb), id_cols)
metric_cols <- metric_cols[map_lgl(metric_cols, ~is.numeric(pb[[.x]]))]

# metric family tagging
fam <- function(m){
  if (m %in% c("inhalePeakMs","returnCrossMs","exhaleStartMs","exhaleTroughMs","endSymMs",
               "breathLenMs","inhaleVolume","inhaleDuration","peakInspFlow","breathLength")) return("airflow_morphology")
  if (str_detect(m,"^w[1-5]_")) return("window_gamma")
  if (str_detect(m,"^p[1-4]_")) return("segment_gamma")
  if (str_detect(m,"^ap")) return("aperiodic")
  if (str_detect(m,"chirp")) return("chirp")
  if (str_detect(m,"[Bb]urst|dutyCycle|timeAbove|anyBurst|nBursts")) return("burst")
  if (str_detect(m,"freq")) return("ridge_frequency")
  if (str_detect(m,"peak|gammaBump|gammaPeak|rpow|aucZ|bandDb|rpowZ")) return("gamma_power")
  "other"
}
ml <- tibble(metric=metric_cols, family=map_chr(metric_cols, fam))
write_csv(ml, file.path(tdir,"taskcmp_metric_list.csv"))

# long per-breath
write_csv(pb |> select(subject, sessID, sessNum, condition, breathIdx, all_of(metric_cols)),
          file.path(tdir,"taskcmp_perbreath.csv"))

# subject x condition summaries (pool sessions within subject; AD_1+AD_2 -> AD)
sm <- pb |> group_by(subject, condition) |>
  summarise(across(all_of(metric_cols),
                   list(mean=~mean(.x,na.rm=TRUE), median=~median(.x,na.rm=TRUE),
                        sd=~sd(.x,na.rm=TRUE)), .names="{.col}__{.fn}"),
            nBreath=n(), .groups="drop")
write_csv(sm, file.path(tdir,"taskcmp_subject_means.csv"))

cat("subjects x condition:\n"); print(pb |> count(subject, condition) |> pivot_wider(names_from=condition, values_from=n))
cat(sprintf("\n%d metrics across %d families; %d breaths total (%d audiobook, %d focus)\n",
            length(metric_cols), n_distinct(ml$family), nrow(pb),
            sum(pb$condition=="audiobook"), sum(pb$condition=="focus")))
cat("families:\n"); print(ml |> count(family))
