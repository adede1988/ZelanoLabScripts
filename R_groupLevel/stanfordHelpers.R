
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




merge_fill <- function(df1, df2, by) {
  out <- full_join(df1, df2, by = by, suffix = c(".x", ".y"))
  
  overlap <- intersect(setdiff(names(df1), by), setdiff(names(df2), by))
  
  for (nm in overlap) {
    out[[nm]] <- coalesce(out[[paste0(nm, ".x")]], out[[paste0(nm, ".y")]])
  }
  
  out %>%
    select(-ends_with(".x"), -ends_with(".y"))
}




plot_session_dumbbell <- function(results, subids, varname, customY,
                                  test_type = c("t", "wilcox"),
                                  seed = 123) {
  source("G:/My Drive/GitHub/MasterStatsUsingR/courseTheme.R")
  
  test_type <- match.arg(test_type)
  
  if (!varname %in% names(results)) {
    stop(sprintf("Variable '%s' not found in results.", varname))
  }
  
  dat <- results %>%
    dplyr::filter(SUBID %in% subids, SESSNUM %in% c(1, 2)) %>%
    dplyr::mutate(
      TYPE_clean = dplyr::case_when(
        TYPE %in% c("Dupi", "DUPI") ~ "Patient",
        TYPE == "OBE" ~ "Control",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(.data[[varname]]), !is.na(TYPE_clean)) %>%
    dplyr::group_by(SUBID, SESSNUM, TYPE_clean) %>%
    dplyr::summarise(value = mean(.data[[varname]], na.rm = TRUE), .groups = "drop")
  
  if (nrow(dat) == 0) {
    stop("No non-missing data available for the requested SUBIDs and variable.")
  }
  
  # --- paired test within patients only ---
  patient_wide <- dat %>%
    dplyr::filter(TYPE_clean == "Patient") %>%
    dplyr::select(SUBID, SESSNUM, value) %>%
    tidyr::pivot_wider(
      names_from = SESSNUM,
      values_from = value,
      names_prefix = "sess"
    ) %>%
    tidyr::drop_na(sess1, sess2)
  
  if (nrow(patient_wide) < 2) {
    test_res <- NULL
    test_label <- "Patient paired test not run: fewer than 2 complete pairs"
  } else {
    if (test_type == "t") {
      test_res <- t.test(patient_wide$sess2, patient_wide$sess1, paired = TRUE)
      test_label <- sprintf(
        "Patients: paired t-test, t(%0.1f) = %.3f, p = %.4g, n = %d",
        unname(test_res$parameter),
        unname(test_res$statistic),
        test_res$p.value,
        nrow(patient_wide)
      )
    } else {
      test_res <- wilcox.test(patient_wide$sess2, patient_wide$sess1, paired = TRUE, exact = FALSE)
      test_label <- sprintf(
        "Patients: paired Wilcoxon, V = %.3f, p = %.4g, n = %d",
        unname(test_res$statistic),
        test_res$p.value,
        nrow(patient_wide)
      )
    }
  }
  
  set.seed(seed)
  subj_offsets <- tibble::tibble(
    SUBID = unique(dat$SUBID),
    x_offset = runif(length(unique(dat$SUBID)), -0.08, 0.08)
  )
  
  dat <- dat %>%
    dplyr::left_join(subj_offsets, by = "SUBID") %>%
    dplyr::mutate(x = SESSNUM + x_offset)
  
  label_dat <- dat %>%
    dplyr::filter(SESSNUM == 1)
  
  outPlot <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = value, group = SUBID)) +
    ggplot2::geom_line(ggplot2::aes(color = TYPE_clean), linewidth = 2, alpha = 1) +
    ggplot2::geom_point(ggplot2::aes(color = TYPE_clean, shape = TYPE_clean), size = 13) +
    ggrepel::geom_text_repel(
      data = label_dat,
      ggplot2::aes(label = SUBID, color = TYPE_clean),
      size = 12,
      nudge_x = -1.0,
      show.legend = FALSE,
      direction = "x",
      box.padding = 1.5,
      point.padding = 1.5,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Sess 1", "Sess 2"),
      limits = c(-.2, 2.25)
    ) +
    scale_color_manual(values = c("Patient" = "#84E642", 'Control' = '#20CDF2'))+
    ggplot2::scale_shape_manual(values = c("Control" = 16, "Patient" = 17)) +
    ggplot2::labs(
      x = NULL,
      y = customY,
      color = "Group",
      shape = "Group",
      subtitle = test_label
    ) +
    myTheme
  
  return(list(
    plot = outPlot,
    test = test_res,
    paired_data = patient_wide
  ))
}


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




circ_mean_deg_wt <- function(deg, w) {
  rad <- deg * pi / 180
  mu_rad <- atan2(
    weighted.mean(sin(rad), w, na.rm = TRUE),
    weighted.mean(cos(rad), w, na.rm = TRUE)
  )
  (mu_rad * 180 / pi) %% 360
}



circ_dist_deg <- function(a, b) {
  d <- (a - b + 180) %% 360 - 180
  abs(d)
}

signed_circ_diff_deg <- function(a, b) {
  (a - b + 180) %% 360 - 180
}