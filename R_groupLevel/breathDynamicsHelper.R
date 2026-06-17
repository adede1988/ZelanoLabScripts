

# -------------------------
# helper functions
# -------------------------

first_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  x[1]
}

safe_cv <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(m) || m == 0) return(NA_real_)
  s / m
}

safe_rmssd <- function(x) {
  d <- diff(x)
  if (sum(!is.na(d)) == 0) return(NA_real_)
  sqrt(mean(d^2, na.rm = TRUE))
}

safe_cor <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

safe_div <- function(num, den) {
  if (!is.finite(den) || is.na(den) || den <= 0) return(NA_real_)
  num / den
}

robust_high <- function(x, k = 3) {
  med_x <- median(x, na.rm = TRUE)
  mad_x <- mad(x, na.rm = TRUE)
  
  out <- rep(FALSE, length(x))
  
  if (!is.finite(med_x) || !is.finite(mad_x) || mad_x == 0) {
    return(out)
  }
  
  out <- x > med_x + k * mad_x
  out[is.na(out)] <- FALSE
  out
}

run_lengths <- function(x, value = TRUE) {
  x <- as.logical(x)
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return(integer(0))
  
  r <- rle(x)
  r$lengths[r$values == value]
}

longest_run <- function(x, value = TRUE) {
  runs <- run_lengths(x, value)
  if (length(runs) == 0) return(0L)
  max(runs)
}

mean_run <- function(x, value = TRUE) {
  runs <- run_lengths(x, value)
  if (length(runs) == 0) return(0)
  mean(runs)
}

transition_rate <- function(x) {
  x <- as.logical(x)
  x <- x[!is.na(x)]
  
  if (length(x) < 2) return(NA_real_)
  
  mean(diff(as.integer(x)) != 0)
}

clock_time_sec <- function(start_samp, end_samp, fs = 500) {
  ok <- is.finite(start_samp) & is.finite(end_samp)
  if (!any(ok)) return(NA_real_)
  
  total_samp <- max(end_samp[ok], na.rm = TRUE) - min(start_samp[ok], na.rm = TRUE)
  total_samp / fs
}


# -----------------------------
# helpers for stats text
# -----------------------------
p_fmt <- function(p) {
  if (is.na(p)) {
    "NA"
  } else if (p < .001) {
    "< .001 \n"
  } else {
    paste0(sub("^0", "", sprintf("%.3f", p)), '\n')
  }
}

short_contrast <- function(x) {
  x %>%
    gsub("No Wave", "No", .) %>%
    gsub("Fast Wave", "Fast", .) %>%
    gsub("Slow Wave", "Slow", .)
}

make_stats_text <- function(dat, yvar) {
  form <- as.formula(paste0(yvar, " ~ combinedCondition + (1 | subjectID)"))
  fit  <- lmer(form, data = dat, REML = FALSE)
  
  aov_tab <- as.data.frame(anova(fit))
  overall_row <- aov_tab["combinedCondition", ]
  
  overall_txt <- paste0(
    "Overall: F(",
    overall_row$NumDF, ", ",
    sprintf("%.1f", overall_row$DenDF), ") = ",
    sprintf("%.2f", overall_row[["F value"]]), ", p ",
    p_fmt(overall_row[["Pr(>F)"]]), '\n'
  )
  
  pw <- emmeans(fit, pairwise ~ combinedCondition, adjust = "holm")$contrasts %>%
    as.data.frame()
  
  pair_txt <- paste0(
    short_contrast(pw$contrast),
    ": p ",
    sapply(pw$p.value, p_fmt),
    collapse = "; "
  )
  
  paste0(overall_txt, "\nPairs (Holm): ", pair_txt)
}

# -----------------------------
# generic panel maker
# -----------------------------
make_panel <- function(plot_data, line_data, yvar, title, subtitle, ylab,
                       percent_axis = FALSE, y_limits = NULL) {
  
  stat_label <- make_stats_text(plot_data, yvar)
  
  p <- ggplot(plot_data, aes(x = combinedCondition, y = .data[[yvar]])) +
    geom_boxplot(width = 0.65, outlier.shape = NA) +
    geom_line(
      data = line_data,
      aes(group = subjectID),
      alpha = 0.35,
      color = "gray40"
    ) +
    geom_point(
      aes(color = subjectID),
      alpha = 0.8,
      size = 2
    ) +
    annotate(
      "text",
      x = 1.02,
      y = Inf,
      label = stat_label,
      hjust = 0,
      vjust = 1.1,
      size = 3.2,
      lineheight = 1.1
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      y = ylab
    ) +
    coord_cartesian(clip = "off") +
    fig_theme
  
  if (percent_axis) {
    if (is.null(y_limits)) y_limits <- c(0, 1)
    p <- p +
      scale_y_continuous(
        labels = percent_format(accuracy = 1),
        limits = y_limits,
        expand = expansion(mult = c(0.02, 0.18))
      )
  } else {
    p <- p +
      scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.22))
      )
  }
  
  p
}