library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(ggrepel)
# -------------------------------------------------------------------
# folder containing behavioral files
# -------------------------------------------------------------------
beh_dir <- "R:/Neurology/Zelano_Lab/Lab_Common/QuestMirror/CHANDAT_processed/BehFiles"

# -------------------------------------------------------------------
# list files and parse filename metadata
# expected filename forms:
# 1) DATE_TYPE_LOCATION_SUBID_SESSNUM_macro_TASKNAME_CHANNUM_behFile.csv
# 2) DATE_TYPE_LOCATION_SUBID_macro_TASKNAME_CHANNUM_behFile.csv
# -------------------------------------------------------------------
all_files <- list.files(
  beh_dir,
  pattern = "_behFile\\.csv$",
  full.names = TRUE
)

# regex with optional SESSNUM
rx <- "^([^_]+)_([^_]+)_([^_]+)_([^_]+)(?:_([^_]+))?_macro_(breathingTask|cueTask|O15|threshTask)_([^_]+)_behFile\\.csv$"

file_directory <- tibble(
  file = all_files,
  base = basename(all_files)
) %>%
  mutate(m = str_match(base, rx)) %>%
  transmute(
    file,
    base,
    DATE        = m[, 2],
    TYPE        = m[, 3],
    LOCATION    = m[, 4],
    SUBID       = m[, 5],
    SESSNUM_raw = m[, 6],
    TASKNAME    = m[, 7],
    CHANNUM     = m[, 8]
  )

# keep only files that parsed successfully
bad_files <- file_directory %>% filter(is.na(SUBID))
if (nrow(bad_files) > 0) {
  warning(
    "Some files did not match the expected naming convention and were dropped:\n",
    paste(bad_files$base, collapse = "\n")
  )
}

file_directory <- file_directory %>%
  filter(!is.na(SUBID)) %>%
  mutate(
    SESSNUM_from_name = suppressWarnings(as.integer(SESSNUM_raw)),
    CHANNUM_num = suppressWarnings(as.numeric(CHANNUM))
  ) %>%
  group_by(SUBID) %>%
  mutate(
    any_missing_sess_for_sub = any(is.na(SESSNUM_from_name)),
    SESSNUM = if_else(any_missing_sess_for_sub, 1L, SESSNUM_from_name)
  ) %>%
  ungroup()

# -------------------------------------------------------------------
# choose one file per SUBID / SESSNUM / TASKNAME
# breathingTask is excluded from behavioral analysis
# -------------------------------------------------------------------
selected_files <- file_directory %>%
  filter(TASKNAME != "breathingTask") %>%
  arrange(SUBID, SESSNUM, TASKNAME, CHANNUM_num, file) %>%
  group_by(SUBID, SESSNUM, TASKNAME) %>%
  slice(1) %>%
  ungroup()

# -------------------------------------------------------------------
# helper: safe NA tibble constructors
# -------------------------------------------------------------------
na_cue <- function() {
  tibble(
    cueTask_HR = NA_real_,
    cueTask_FA = NA_real_,
    cueTask_d  = NA_real_
  )
}

na_o15 <- function() {
  tibble(
    O15_score = NA_real_,
    O15_acc   = NA_real_
  )
}

na_thresh <- function() {
  tibble(
    thresh_none            = NA_real_,
    thresh_low             = NA_real_,
    thresh_high            = NA_real_,
    thresh_low_calibrated  = NA_real_,
    thresh_high_calibrated = NA_real_
  )
}

# -------------------------------------------------------------------
# cueTask analysis
# HIT  = cue==odor & respString=="Yes"
# MISS = cue==odor & respString=="No"
# FA   = cue!=odor & respString=="Yes"
# CR   = cue!=odor & respString=="No"
#
# cueTask_HR = hits / (hits + misses)
# cueTask_FA = false alarms / (false alarms + correct rejections)
# cueTask_d  = d-prime
#
# For d', a standard log-linear correction is used so qnorm never hits Inf
# -------------------------------------------------------------------
analyze_cueTask <- function(path) {
  dat <- tryCatch(
    read_csv(path, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(dat)) return(na_cue())
  
  needed <- c("cue", "odor", "respString")
  if (!all(needed %in% names(dat))) return(na_cue())
  
  dat <- dat %>%
    transmute(
      cue = suppressWarnings(as.numeric(cue)),
      odor = suppressWarnings(as.numeric(odor)),
      respString = str_trim(as.character(respString))
    ) %>%
    filter(
      !is.na(cue),
      !is.na(odor),
      respString %in% c("Yes", "No")
    )
  
  if (nrow(dat) == 0) return(na_cue())
  
  hits   <- sum(dat$cue == dat$odor & dat$respString == "Yes")
  misses <- sum(dat$cue == dat$odor & dat$respString == "No")
  fas    <- sum(dat$cue != dat$odor & dat$respString == "Yes")
  crs    <- sum(dat$cue != dat$odor & dat$respString == "No")
  
  sig_n <- hits + misses
  noi_n <- fas + crs
  
  hr <- if (sig_n > 0) hits / sig_n else NA_real_
  fa <- if (noi_n > 0) fas / noi_n else NA_real_
  
  dprime <- if (sig_n > 0 && noi_n > 0) {
    # log-linear correction for extreme rates
    hr_adj <- (hits + 0.5) / (sig_n + 1)
    fa_adj <- (fas  + 0.5) / (noi_n + 1)
    qnorm(hr_adj) - qnorm(fa_adj)
  } else {
    NA_real_
  }
  
  tibble(
    cueTask_HR = hr,
    cueTask_FA = fa,
    cueTask_d  = dprime
  )
}

# -------------------------------------------------------------------
# O15 analysis
# keep one unique row per n, then:
# O15_score = sum(expScore)
# O15_acc   = O15_score / 15
# -------------------------------------------------------------------
analyze_O15 <- function(path) {
  dat <- tryCatch(
    read_csv(path, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(dat)) return(na_o15())
  
  needed <- c("n", "expScore")
  if (!all(needed %in% names(dat))) return(na_o15())
  
  dat <- dat %>%
    transmute(
      n = suppressWarnings(as.numeric(n)),
      expScore = suppressWarnings(as.numeric(expScore))
    ) %>%
    filter(!is.na(n)) %>%
    distinct(n, .keep_all = TRUE)
  
  if (nrow(dat) == 0 || all(is.na(dat$expScore))) return(na_o15())
  
  score <- sum(dat$expScore, na.rm = TRUE)
  
  tibble(
    O15_score = score,
    O15_acc   = score / 15
  )
}

# -------------------------------------------------------------------
# threshTask analysis
# mean intensity by odor:
# odor==1 -> thresh_none
# odor==2 -> thresh_low
# odor==3 -> thresh_high
# calibrated = subtract thresh_none
# -------------------------------------------------------------------
analyze_threshTask <- function(path) {
  dat <- tryCatch(
    read_csv(path, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(dat)) return(na_thresh())
  
  needed <- c("odor", "intensity")
  if (!all(needed %in% names(dat))) return(na_thresh())
  
  dat <- dat %>%
    transmute(
      odor = suppressWarnings(as.numeric(odor)),
      intensity = suppressWarnings(as.numeric(intensity))
    ) %>%
    filter(odor %in% c(1, 2, 3))
  
  if (nrow(dat) == 0) return(na_thresh())
  
  odor_means <- dat %>%
    group_by(odor) %>%
    summarise(
      mean_intensity = if (all(is.na(intensity))) NA_real_ else mean(intensity, na.rm = TRUE),
      .groups = "drop"
    )
  
  get_mean <- function(odor_val) {
    x <- odor_means$mean_intensity[odor_means$odor == odor_val]
    if (length(x) == 0) NA_real_ else x[1]
  }
  
  thresh_none <- get_mean(1)
  thresh_low  <- get_mean(2)
  thresh_high <- get_mean(3)
  
  thresh_low_calibrated <- if (
    is.na(thresh_low) || is.na(thresh_none)
  ) NA_real_ else thresh_low - thresh_none
  
  thresh_high_calibrated <- if (
    is.na(thresh_high) || is.na(thresh_none)
  ) NA_real_ else thresh_high - thresh_none
  
  tibble(
    thresh_none            = thresh_none,
    thresh_low             = thresh_low,
    thresh_high            = thresh_high,
    thresh_low_calibrated  = thresh_low_calibrated,
    thresh_high_calibrated = thresh_high_calibrated
  )
}

# -------------------------------------------------------------------
# run task-specific summaries on the one selected file per
# SUBID / SESSNUM / TASKNAME
# -------------------------------------------------------------------
cue_results <- selected_files %>%
  filter(TASKNAME == "cueTask") %>%
  mutate(out = map(file, analyze_cueTask)) %>%
  unnest(out) %>%
  select(SUBID, SESSNUM, cueTask_HR, cueTask_FA, cueTask_d)

o15_results <- selected_files %>%
  filter(TASKNAME == "O15") %>%
  mutate(out = map(file, analyze_O15)) %>%
  unnest(out) %>%
  select(SUBID, SESSNUM, O15_score, O15_acc)

thresh_results <- selected_files %>%
  filter(TASKNAME == "threshTask") %>%
  mutate(out = map(file, analyze_threshTask)) %>%
  unnest(out) %>%
  select(
    SUBID, SESSNUM,
    thresh_none, thresh_low, thresh_high,
    thresh_low_calibrated, thresh_high_calibrated
  )

# -------------------------------------------------------------------
# preallocate final results table:
# one row per unique SUBID / SESSNUM
# TYPE is taken as the first TYPE observed for that SUBID / SESSNUM
# includes rows even if only breathingTask exists
# -------------------------------------------------------------------
results <- file_directory %>%
  group_by(SUBID, SESSNUM) %>%
  summarise(
    TYPE = first(TYPE[!is.na(TYPE)]),
    .groups = "drop"
  ) %>%
  arrange(SUBID, SESSNUM) %>%
  mutate(
    cueTask_HR = NA_real_,
    cueTask_FA = NA_real_,
    cueTask_d  = NA_real_,
    O15_score  = NA_real_,
    O15_acc    = NA_real_,
    thresh_none            = NA_real_,
    thresh_low             = NA_real_,
    thresh_high            = NA_real_,
    thresh_low_calibrated  = NA_real_,
    thresh_high_calibrated = NA_real_
  ) %>%
  left_join(cue_results,    by = c("SUBID", "SESSNUM"), suffix = c("", ".new")) %>%
  left_join(o15_results,    by = c("SUBID", "SESSNUM"), suffix = c("", ".new")) %>%
  left_join(thresh_results, by = c("SUBID", "SESSNUM"), suffix = c("", ".new")) %>%
  mutate(
    cueTask_HR = coalesce(cueTask_HR.new, cueTask_HR),
    cueTask_FA = coalesce(cueTask_FA.new, cueTask_FA),
    cueTask_d  = coalesce(cueTask_d.new,  cueTask_d),
    O15_score  = coalesce(O15_score.new,  O15_score),
    O15_acc    = coalesce(O15_acc.new,    O15_acc),
    thresh_none            = coalesce(thresh_none.new,            thresh_none),
    thresh_low             = coalesce(thresh_low.new,             thresh_low),
    thresh_high            = coalesce(thresh_high.new,            thresh_high),
    thresh_low_calibrated  = coalesce(thresh_low_calibrated.new,  thresh_low_calibrated),
    thresh_high_calibrated = coalesce(thresh_high_calibrated.new, thresh_high_calibrated)
  ) %>%
  select(
    SUBID,
    SESSNUM,
    TYPE,
    cueTask_HR,
    cueTask_FA,
    cueTask_d,
    O15_score,
    O15_acc,
    thresh_none,
    thresh_low,
    thresh_high,
    thresh_low_calibrated,
    thresh_high_calibrated
  )

# -------------------------------------------------------------------
# objects created:
#   file_directory  = parsed directory of all matching files
#   selected_files  = one chosen file per SUBID/SESSNUM/TASKNAME (no breathingTask)
#   results         = final behavioral summary table
# -------------------------------------------------------------------

print(results)

# optional:
# write_csv(results, file.path(beh_dir, "behavior_summary_results.csv"))




plot_session_dumbbell <- function(results, subids, varname, customY, seed = 123) {
  source("G:/My Drive/GitHub/MasterStatsUsingR/courseTheme.R")
  
  if (!varname %in% names(results)) {
    stop(sprintf("Variable '%s' not found in results.", varname))
  }
  
  dat <- results %>%
    filter(SUBID %in% subids, SESSNUM %in% c(1, 2)) %>%
    mutate(
      TYPE_clean = case_when(
        TYPE %in% c("Dupi", "DUPI") ~ "Patient",
        TYPE == "OBE" ~ "Control",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(.data[[varname]]), !is.na(TYPE_clean))
  
  if (nrow(dat) == 0) {
    stop("No non-missing data available for the requested SUBIDs and variable.")
  }
  
  set.seed(seed)
  subj_offsets <- tibble(
    SUBID = unique(dat$SUBID),
    x_offset = runif(length(unique(dat$SUBID)), -0.08, 0.08)
  )
  
  dat <- dat %>%
    left_join(subj_offsets, by = "SUBID") %>%
    mutate(x = SESSNUM + x_offset)
  
  label_dat <- dat %>%
    filter(SESSNUM == 1)
  
  outPlot = ggplot(dat, aes(x = x, y = .data[[varname]], group = SUBID)) +
    geom_line(aes(color = TYPE_clean), linewidth = 2, alpha = 1) +
    geom_point(aes(color = TYPE_clean, shape = TYPE_clean), size = 13) +
    geom_text_repel(
      data = label_dat,
      aes(label = SUBID, color = TYPE_clean),
      size = 12,
      nudge_x = -1.0,
      show.legend = FALSE,
      direction = "x",
      box.padding = 1.5,
      point.padding = 1.5,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Sess 1", "Sess 2"),
      limits = c(-.2, 2.25)
    ) +
    brightCol +
    scale_shape_manual(values = c("Control" = 16, "Patient" = 17)) +
    labs(x = NULL, y = customY, color = "Group", shape = "Group") +
    myTheme
  return(outPlot)
}







subList = c("KS", "JH", "GH", "AB", "DB", "AS", 
            "AZ", "CP", "CS", "FS", "RY", "TI")

outPlot = plot_session_dumbbell(results, subList,
                    "cueTask_d", "Accuracy (d')")

fn <- paste( "G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\" ,
             'cueTask', '.png', 
             sep = '')
png(fn,         # File name
    width=600, height=600)
print(outPlot)
dev.off()



outPlot = plot_session_dumbbell(results, subList,
                                "O15_acc", "Accuracy (% correct)")

fn <- paste( "G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\" ,
             'O15', '.png', 
             sep = '')
png(fn,         # File name
    width=600, height=600)
print(outPlot)
dev.off()


#hard code for missing JH session 2 threshold, sub in his session 3
results$thresh_low_calibrated[15] = results$thresh_low_calibrated[16] / 2
results$thresh_high_calibrated[15] = results$thresh_high_calibrated[16] / 2

outPlot = plot_session_dumbbell(results, subList,
                                "thresh_low_calibrated", "detection above 0 (au)" )

outPlot + scale_color_manual(values = c("Patient" = "#84E642")) -> outPlot

fn <- paste( "G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\" ,
             'threshLow', '.png', 
             sep = '')
png(fn,         # File name
    width=600, height=600)
print(outPlot)
dev.off()


results <- results %>%
  mutate(
    cs_O15    = O15_acc,
    cs_thresh = thresh_low_calibrated / 400,
    cs_cue    = cueTask_d / 3.5,
    
    combinedScore_n = 
      (!is.na(cs_O15)    & cs_O15    != 0) +
      (!is.na(cs_thresh) & cs_thresh != 0) +
      (!is.na(cs_cue)    & cs_cue    != 0),
    
    combinedScore = ifelse(
      combinedScore_n > 0,
      (coalesce(cs_O15, 0) + coalesce(cs_thresh, 0) + coalesce(cs_cue, 0)) / combinedScore_n,
      NA_real_
    )
  ) %>%
  select(-cs_O15, -cs_thresh, -cs_cue, -combinedScore_n)




outPlot = plot_session_dumbbell(results, subList,
                                "combinedScore", "Olfaction Index")

fn <- paste( "G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\" ,
             'combined', '.png', 
             sep = '')
png(fn,         # File name
    width=600, height=600)
print(outPlot)
dev.off()



conEphys = read.csv("G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\ConDat.csv")
patEphys = read.csv("G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\PatDat.csv")

conEphys = conEphys %>% select(-c(sessID, comboID, unusedChan, SessNum2, TYPE, GlobalID, subID2))
patEphys = patEphys %>% select(-c(sessID, comboID, unusedChan, SessNum2, TYPE, GlobalID, subID2))

names(conEphys)[1] = 'SUBID'
names(patEphys)[1] = 'SUBID'

test = merge(results, conEphys)



merge_fill <- function(df1, df2, by) {
  out <- full_join(df1, df2, by = by, suffix = c(".x", ".y"))
  
  overlap <- intersect(setdiff(names(df1), by), setdiff(names(df2), by))
  
  for (nm in overlap) {
    out[[nm]] <- coalesce(out[[paste0(nm, ".x")]], out[[paste0(nm, ".y")]])
  }
  
  out %>%
    select(-ends_with(".x"), -ends_with(".y"))
}

out <- merge_fill(results, conEphys, by = c('SUBID', 'SESSNUM'))

out <- merge_fill(out, patEphys, by = c('SUBID', 'SESSNUM'))



## NOW PLOT EPHYS CHANGES! 
plot_session_dumbbell(out, subList,
                      "peakProm", "gamPeakProm")

plot_session_dumbbell(out, subList,
                      "R", "phasePrefStrength")


plot_session_dumbbell(out, subList,
                      "muDeg", "phasePrefAngle")

out$muDif = ((out$muDeg - 115 + 180) %% 360) - 180

plot_session_dumbbell(out, subList,
                      "muDif", "phasePrefAngle")

out$bLenDif = out$bLen_inPhase - out$bLen_outPhase
out$bAmpDif = out$bAmp_inPhase - out$bAmp_outPhase
plot_session_dumbbell(out, subList,
                      "bAmpDif", "phasePrefAngle")


out$inOutDif = out$peakVals_inPhase - out$peakVals_outPhase 
plot_session_dumbbell(out, subList,
                      "inOutDif", "InPhaseProminence")

plot_session_dumbbell(out, subList,
                      "inPhase", "in phase proportion")

plot_session_dumbbell(out, subList,
                      "gamFreq_inPhase", "gam Freq")

plot_session_change_scatter <- function(results, subids, varname1, varname2,
                                        customX, customY) {
  source("G:/My Drive/GitHub/MasterStatsUsingR/courseTheme.R")
  
  if (!varname1 %in% names(results)) {
    stop(sprintf("Variable '%s' not found in results.", varname1))
  }
  if (!varname2 %in% names(results)) {
    stop(sprintf("Variable '%s' not found in results.", varname2))
  }
  
  dat <- results %>%
    filter(SUBID %in% subids, SESSNUM %in% c(1, 2)) %>%
    mutate(
      TYPE_clean = case_when(
        TYPE %in% c("Dupi", "DUPI") ~ "Patient",
        TYPE == "OBE" ~ "Control",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(TYPE_clean)) %>%
    select(SUBID, SESSNUM, TYPE_clean, all_of(varname1), all_of(varname2))
  
  diff_dat <- dat %>%
    tidyr::pivot_wider(
      id_cols = c(SUBID, TYPE_clean),
      names_from = SESSNUM,
      values_from = c(all_of(varname1), all_of(varname2)),
      names_sep = "_"
    ) %>%
    mutate(
      diff1 = .data[[paste0(varname1, "_2")]] - .data[[paste0(varname1, "_1")]],
      diff2 = .data[[paste0(varname2, "_2")]] - .data[[paste0(varname2, "_1")]]
    ) %>%
    filter(!is.na(diff1), !is.na(diff2))
  
  if (nrow(diff_dat) == 0) {
    stop("No subjects have complete session 1 and session 2 data for both variables.")
  }
  
  outPlot <- ggplot(diff_dat, aes(x = diff1, y = diff2)) +
    geom_hline(yintercept = 0, linewidth = 1, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 1, color = "grey70") +
    geom_point(aes(color = TYPE_clean, shape = TYPE_clean), size = 13) +
    geom_text_repel(
      aes(label = SUBID, color = TYPE_clean),
      size = 12,
      show.legend = FALSE,
      box.padding = 1.0,
      point.padding = 1.0,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +
    brightCol +
    scale_shape_manual(values = c("Control" = 16, "Patient" = 17)) +
    labs(
      x = paste0("\u0394 ", customX),
      y = paste0("\u0394 ", customY),
      color = "Group",
      shape = "Group"
    ) +
    myTheme
  
  return(outPlot)
}
############################################################################
#save out plots: 

outPlot = plot_session_dumbbell(out, subList,
                      "peakProm", "gam prominence")
fn <- paste( "G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\" ,
             'gamProm', '.png', 
             sep = '')
png(fn,         # File name
    width=600, height=600)
print(outPlot)
dev.off()


outPlot = plot_session_dumbbell(out, subList,
                                "gamFreq_inPhase", "gam Freq")
fn <- paste( "G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\" ,
             'FreqVals', '.png', 
             sep = '')
png(fn,         # File name
    width=600, height=600)
print(outPlot)
dev.off()



outPlot = plot_session_dumbbell(out, subList,
                                "muDeg", "Gamma Preferred Resp Phase")
fn <- paste( "G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\" ,
             'prefPhase', '.png', 
             sep = '')
png(fn,         # File name
    width=600, height=600)
print(outPlot)
dev.off()


out$HRVindex = (out$HRV_RMS_inPhase*1000 + out$HRV_RMS_outPhase*1000) / 2 +
               (out$HRV_RSA_inPhase + out$HRV_RSA_outPhase) / 2+
               (out$HRV_SDNN_inPhase*100 + out$HRV_SDNN_outPhase*100) / 2

plot_session_dumbbell(out, subList,
                      "HRVindex", "HRV")

## compare difference scores: 



# subList = c("KS", "JH", "GH", "AB", "DB")

outPlot = plot_session_change_scatter(
  out,
  subids = subList,
  varname1 = "combinedScore",
  varname2 = "muDif",
  customX = "Olfactory capability",
  customY = "Gamma peak in phase"
)
outPlot + scale_color_manual(values = c("Patient" = "#84E642")) -> outPlot
fn <- paste( "G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\" ,
             'phasePref_olfactionImprove', '.png', 
             sep = '')
png(fn,         # File name
    width=600, height=600)
print(outPlot)
dev.off()


plot_session_change_scatter(
  out,
  subids = subList,
  varname1 = "bLenDif",
  varname2 = "inPhase",
  customX = "breathLen Change",
  customY = "prop in phase breaths"
)
plot_session_change_scatter(
  out,
  subids = subList,
  varname1 = "bAmpDif",
  varname2 = "inPhase",
  customX = "breathAmp Change",
  customY = "prop in phase breaths"
)



outPlot = plot_session_dumbbell(out, subList,
                      "inPhase", "in phase proportion")

fn <- paste( "G:\\My Drive\\cZelano\\ACHEMS_2026\\figs\\" ,
             'inPhaseProp', '.png', 
             sep = '')
png(fn,         # File name
    width=600, height=600)
print(outPlot)
dev.off()

# outPlot = plot_session_change_scatter(
#   out,
#   subids = subList,
#   varname1 = "combinedScore",
#   varname2 = "inPhase",
#   customX = "Olfactory capability",
#   customY = "Gamma peak in phase"
# )
# outPlot + scale_color_manual(values = c("Patient" = "#84E642")) -> outPlot

circ_mean_deg_wt <- function(deg, w) {
  rad <- deg * pi / 180
  mu_rad <- atan2(
    weighted.mean(sin(rad), w, na.rm = TRUE),
    weighted.mean(cos(rad), w, na.rm = TRUE)
  )
  (mu_rad * 180 / pi) %% 360
}

control_mean_deg_wt <- out %>%
  filter(TYPE == "OBE") %>%
  summarise(mu = circ_mean_deg_wt(muDeg, R)) %>%
  pull(mu)



theta_ref_deg <- control_mean_deg_wt

circ_dist_deg <- function(a, b) {
  d <- (a - b + 180) %% 360 - 180
  abs(d)
}

signed_circ_diff_deg <- function(a, b) {
  (a - b + 180) %% 360 - 180
}

dat2 <- out %>%
  mutate(
    theta_rad = muDeg * pi / 180,
    theta_ref_rad = theta_ref_deg * pi / 180,
    circ_dist = circ_dist_deg(muDeg, theta_ref_deg),
    align_cos = cos(theta_rad - theta_ref_rad),
    aligned_score = R * align_cos,
    x = R * cos(theta_rad),
    y = R * sin(theta_rad)
  )

plot_session_dumbbell(dat2, subList,
                      "aligned_score", "alignment score")


dupi_change <- dat2 %>%
  filter(TYPE %in% c("Dupi", "DUPI"),
         SESSNUM %in% c(1, 2)) %>%
  select(SUBID, SESSNUM, aligned_score, inPhase, R, circ_dist) %>%
  tidyr::pivot_wider(
    names_from = SESSNUM,
    values_from = c(aligned_score, inPhase, R, circ_dist),
    names_prefix = "sess"
  ) %>%
  mutate(
    d_aligned = aligned_score_sess2 - aligned_score_sess1,
    d_inPhase = inPhase_sess2 - inPhase_sess1,
    d_R = R_sess2 - R_sess1,
    d_dist = circ_dist_sess2 - circ_dist_sess1
  )


plot_session_change_scatter(
  dat2,
  subids = subList,
  varname1 = "combinedScore",
  varname2 = "aligned_score",
  customX = "Olfactory capability",
  customY = "Gamma peak in phase"
)

