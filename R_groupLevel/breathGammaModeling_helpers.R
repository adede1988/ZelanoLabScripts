
# --------------------------
# Helpers
# --------------------------
parse_file_meta <- function(path) {
  fname <- basename(path)
  base  <- str_remove(fname, "_behFile\\.csv$")
  chan  <- as.integer(str_extract(base, "(?<=_)[0-9]{2}$"))
  stem  <- str_remove(base, "_[0-9]{2}$")
  tibble(file = path, file_name = fname, base = base, stem = stem, channel = chan)
}

mad_robust <- function(x) {
  # comparable to MATLAB robust-ish scaling if you want it; adjust if desired
  # here: raw MAD (constant=1)
  mad(x, constant = 1, na.rm = TRUE)
}

load_all_beh <- function(in_dir) {
  files <- list.files(in_dir, pattern = "_behFile\\.csv$", full.names = TRUE)
  meta  <- map_dfr(files, parse_file_meta)
  
  dat <- meta %>%
    mutate(dat = map(file, ~ readr::read_csv(.x, show_col_types = FALSE)))
  
  allDat = dat[1,] %>% unnest(cols = c(dat))
  
  for(ii in 2:dim(dat)[1]){
    newDat = dat[ii,] %>% unnest(cols = c(dat))
    if(mode(newDat$response) == 'character'){
      newDat$respString = newDat$response
      newDat$response = NA
      newDat$response = as.numeric(newDat$response)
    }
    allDat = bind_rows(list(allDat, newDat))
    
  }

  
  # enforce expected columns (light touch; no heavy checking)
  allDat %>%
    mutate(
      stem    = factor(stem),
      channel = factor(channel),
      use     = as.integer(use)
    ) -> allDat
  return(allDat)
}

to_long_epoch_metrics <- function(dat_qc) {
  # Long-format for epoch-wise metrics
  dat_qc %>%
    pivot_longer(
      cols = matches("^(inRise|inFall|exRise|exFall|pause)(Freq|FitR2|PromLocal|PromWide)$"),
      names_to = c("epoch_code", "metric"),
      names_pattern = "^(inRise|inFall|exRise|exFall|pause)(.*)$",
      values_to = "value"
    ) %>%
    mutate(
      epoch = recode(epoch_code, !!!epoch_key),
      epoch = factor(epoch, levels = unname(epoch_key))
    )
}

score_channels <- function(dat_qc, dat_long) {
  # Breath-level -> channel summaries for quick ranking
  # (You can change the score definition to match your intuition.)
  fit_summ <- dat_long %>%
    filter(metric == "FitR2") %>%
    group_by(stem, channel) %>%
    summarise(
      n_fit   = sum(is.finite(value)),
      prop_good_fit = mean(value >= r2_thresh, na.rm = TRUE),
      med_fit = median(value, na.rm = TRUE),
      mad_fit = mad_robust(value),
      .groups = "drop"
    )

  prom_local_summ <- dat_long %>%
    filter(metric == "PromLocal") %>%
    group_by(stem, channel) %>%
    summarise(
      n_prom_local = sum(is.finite(value)),
      med_prom_local = median(value, na.rm = TRUE),
      mad_prom_local = mad_robust(value),
      .groups = "drop"
    )

  prom_wide_summ <- dat_long %>%
    filter(metric == "PromWide") %>%
    group_by(stem, channel) %>%
    summarise(
      n_prom_wide = sum(is.finite(value)),
      med_prom_wide = median(value, na.rm = TRUE),
      mad_prom_wide = mad_robust(value),
      .groups = "drop"
    )

  burst_summ <- dat_qc %>%
    group_by(stem, channel) %>%
    summarise(
      med_prom_burst_f = median(prom_burst_f, na.rm = TRUE),
      mad_prom_burst_f = mad_robust(prom_burst_f),
      med_prom_burst_t = median(prom_burst_t, na.rm = TRUE),
      mad_prom_burst_t = mad_robust(prom_burst_t),
      .groups = "drop"
    )

  summ <- fit_summ %>%
    full_join(prom_local_summ, by = c("stem","channel")) %>%
    full_join(prom_wide_summ,  by = c("stem","channel")) %>%
    full_join(burst_summ,      by = c("stem","channel")) %>%
    group_by(stem) %>%
    mutate(
      # z-scores within stem so sessions are comparable across channels
      z_fit   = scale(med_fit)[,1],
      z_prop  = scale(prop_good_fit)[,1],
      z_ploc  = scale(med_prom_local)[,1],
      z_pwid  = scale(med_prom_wide)[,1],
      z_burst = scale(med_prom_burst_t)[,1],

      # "consistency" bonus: lower MAD => better (negate z of MAD)
      z_cons  = -scale(mad_prom_local)[,1],

      # composite score (tweak weights!)
      score = 0.30*z_prop + 0.25*z_fit + 0.25*z_ploc + 0.10*z_burst + 0.10*z_cons
    ) %>%
    ungroup()

  summ
}

make_plots_for_stem <- function(dat_stem, out_dir, stem_name) {
  dat_qc   <- dat_stem %>% filter(use == 1)
  dat_long <- to_long_epoch_metrics(dat_qc)

  # ---- Plot 1: Proportion of "good" FitR2 per epoch (counts above bars)
  p_r2_prop <- dat_long %>%
    filter(metric == "FitR2") %>%
    group_by(channel, epoch) %>%
    summarise(
      n = sum(is.finite(value)),
      n_good = sum(value >= r2_thresh, na.rm = TRUE),
      prop_good = ifelse(n > 0, n_good / n, NA_real_),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = channel, y = prop_good, fill = epoch)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.85) +
    geom_text(
      aes(label = n_good),
      position = position_dodge(width = 0.9),
      vjust = -0.2, size = 3
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0(stem_name, " | FitR2: proportion of good fits (>= ", r2_thresh, ")"),
      x = "Channel", y = "Proportion (QC breaths)"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

  # ---- Plot 2: FitR2 distributions (epoch on x, facet by channel)
  p_r2_dist <- dat_long %>%
    filter(metric == "FitR2") %>%
    ggplot(aes(x = epoch, y = value)) +
    geom_violin(trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.6) +
    facet_wrap(~ channel, nrow = 1) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "FitR2 distributions by epoch (QC breaths)", x = NULL, y = "FitR2") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  # ---- Plot 3: Prominence distributions by epoch (Local vs Wide), facet by channel
  p_prom <- dat_long %>%
    filter(metric %in% c("PromLocal","PromWide")) %>%
    mutate(metric = factor(metric, levels = c("PromLocal","PromWide"))) %>%
    ggplot(aes(x = epoch, y = value)) +
    geom_violin(trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.6) +
    facet_grid(metric ~ channel) +
    labs(title = "Epoch-wise prominence (QC breaths)", x = NULL, y = "Prominence (log10 units)") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  # ---- Plot 4: Burst prominence summaries (prom_burst_f / prom_burst_t)
  p_burst <- dat_qc %>%
    select(channel, prom_burst_f, prom_burst_t) %>%
    pivot_longer(cols = c(prom_burst_f, prom_burst_t),
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
                           prom_burst_f = "prom_burst_f",
                           prom_burst_t = "prom_burst_t")) %>%
    ggplot(aes(x = channel, y = value)) +
    geom_violin(trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.6) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    labs(title = "Burst prominence summaries (QC breaths)", x = "Channel", y = "Value") +
    theme_classic(base_size = 12)

  # ---- Plot 5: Channel-level summary scatter (prominence vs detectability)
  summ <- score_channels(dat_qc, dat_long)

  p_scatter <- summ %>%
    ggplot(aes(x = med_prom_local, y = prop_good_fit, label = channel)) +
    geom_point(aes(size = med_prom_burst_t, color = score), alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
    labs(
      title = "Channel summary: median PromLocal vs gamma detectability (FitR2 pass-rate)",
      x = "Median PromLocal (QC breaths, pooled epochs)",
      y = paste0("Proportion FitR2 >= ", r2_thresh),
      size = "Median prom_burst_t",
      color = "Composite score"
    ) +
    theme_classic(base_size = 12)

  # Save a multi-page PDF per stem
  pdf_file <- file.path(out_dir, paste0(stem_name, "_channelQC.pdf"))
  grDevices::pdf(pdf_file, width = 13, height = 8.5)
  print(p_r2_prop)
  print(p_r2_dist)
  print(p_prom)
  print(p_burst)
  print(p_scatter)
  grDevices::dev.off()

  # Save the summary table (so you can sort in R quickly)
  summ_file <- file.path(out_dir, paste0(stem_name, "_channelQC_summary.csv"))
  readr::write_csv(summ, summ_file)

  invisible(list(pdf = pdf_file, summary = summ_file, summary_tbl = summ))
}

