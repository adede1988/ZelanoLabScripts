library(tidyverse)

df = read.csv("R:/Neurology/Zelano_Lab/Lab_Common/Adam/Dupi_processing/QC_summary_metrics.csv")

outPath = 'R:\\Neurology\\Zelano_Lab\\Lab_Common\\Adam\\Dupi_processing\\groupStatFigs\\'

# ---- clean / normalize IDs ----
df <- df %>%
  mutate(
    participant_ID = as.character(participant_ID),
    participant_ID = str_replace_all(participant_ID, "DUPI|dupi", "Dupi"),
    # hard-coded patch: TB was misnamed as TPB (with NMH) in session 1
    participant_ID = str_replace(participant_ID, "^250811_Dupi_NMH_TPB$", "250811_Dupi_NMH_TB"),
    session_num = as.integer(session_num),
    performance = as.numeric(performance),
    performanceType = as.character(performanceType)
  )

# Keep only the rows that actually have performance
perf <- df %>%
  filter(!is.na(performanceType), session_num %in% c(1, 2)) %>%
  select(participant_ID, session_num, performanceType, performance) %>%
  group_by(participant_ID, session_num, performanceType) %>%
  summarise(performance = mean(performance, na.rm = TRUE), .groups = "drop")

# helper: identify who has both sessions (for drawing lines)
perf <- perf %>%
  group_by(participant_ID, performanceType) %>%
  mutate(has_both = n_distinct(session_num) == 2) %>%
  ungroup()

library(tidyverse)

# ---- build a consistent participant color palette (reuse across figures) ----
id_levels <- perf %>%
  distinct(participant_ID) %>%
  arrange(participant_ID) %>%
  pull(participant_ID)

# 7 maximally-different hard-coded colors (Okabe–Ito + black)
hard7 <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#000a73", # purple
  "#F0E442", # yellow
  "#009200", # green
  "#994E00", # vermillion
  "#cC00A7"  # pink
)

# Map participants to these 7 colors (recycles if >7 IDs)
pal <- setNames(rep(hard7, length.out = length(id_levels)), id_levels)

# ---- add per-participant-session jitter so points AND line endpoints jitter consistently ----
set.seed(1)  # for reproducible jitter
perf_j <- perf %>%
  mutate(session_num = factor(session_num, levels = c(1, 2))) %>%
  mutate(x_num = as.numeric(session_num)) 
perf_j = perf_j %>%
  group_by(participant_ID, performanceType, session_num) %>%
  summarise(performance = mean(performance, na.rm = TRUE),
            x_num = x_num, .groups = "drop") %>%
  mutate(x_j = x_num + runif(n(), -0.08, 0.08))  # jitter width

# recompute has_both after summarising
perf_j <- perf_j %>%
  group_by(participant_ID, performanceType) %>%
  mutate(has_both = n_distinct(session_num) == 2) %>%
  ungroup()

# -----------------------------
# d' plot with jittered endpoints + consistent palette
# -----------------------------
p_dprime <- perf_j %>%
  filter(performanceType == "d_prime") %>%
  ggplot(aes(x = x_j, y = performance, group = participant_ID, color = participant_ID)) +
  geom_line(data = ~ filter(.x, has_both), alpha = 0.6, linewidth = 0.7) +
  geom_point(size = 5) +
  scale_x_continuous(breaks = c(1, 2), labels = c("1", "2")) +
  scale_color_manual(values = pal, drop = FALSE) +
  labs(x = "Session", y = "d'", title = "CueTask performance (d') from Session 1 to 2") +
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  )

p_dprime


print(p_dprime)

# -----------------------------
# Plot 2: accuracy
# -----------------------------
p_acc <- perf_j %>%
  filter(performanceType == "accuracy") %>%
  ggplot(aes(x = x_j, y = performance, group = participant_ID, color = participant_ID)) +
  geom_line(data = ~ filter(.x, has_both), alpha = 0.6, linewidth = 0.7) +
  geom_point(size = 5) +
  scale_x_continuous(breaks = c(1, 2), labels = c("1", "2")) +
  scale_color_manual(values = pal, drop = FALSE) +
  labs(x = "Session", y = "Accuracy (mean expScore across targets)",
       title = "O15 performance (accuracy) from Session 1 to 2") +
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  )

print(p_acc)

ggsave(paste0(outPath, "performance_dprime_session1to2.png"), p_dprime, width = 7, height = 5, dpi = 300)
ggsave(paste0(outPath, "performance_accuracy_session1to2.png"), p_acc, width = 7, height = 5, dpi = 300)




# ---- compute session2 - session1 changes for each measure ----
chg <- df %>%
  filter(session_num %in% c(1, 2),
         performanceType %in% c("d_prime", "accuracy"),
         !is.na(performance)) %>%
  group_by(participant_ID, performanceType, session_num) %>%
  summarise(performance = mean(performance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = session_num, values_from = performance,
              names_prefix = "sess") %>%
  mutate(change = sess2 - sess1) %>%
  select(participant_ID, performanceType, change) %>%
  pivot_wider(names_from = performanceType, values_from = change,
              names_prefix = "change_")

# keep only participants with both changes available
chg_complete <- chg %>%
  filter(!is.na(change_d_prime) & !is.na(change_accuracy))

# ---- scatter: d' change vs accuracy change ----
p <- ggplot(chg_complete, aes(x = change_d_prime, y = change_accuracy)) +
  geom_hline(yintercept = 0, linewidth = 0.6, alpha = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.6, alpha = 0.6) +
  geom_point(size = 2) +
  # geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Δ d' (Session 2 − Session 1)",
    y = "Δ Accuracy (Session 2 − Session 1)",
    title = "Change–change relationship across participants"
  ) +
  theme_bw(base_size = 12)

print(p)




# ---- keep + reshape peakFreq1 across time windows ----
peak1_long <- df %>%
  select(participant_ID, session_num, task,
         base_peakFreq1, early_peakFreq1, late_peakFreq1) %>%
  pivot_longer(
    cols = c(base_peakFreq1, early_peakFreq1, late_peakFreq1),
    names_to = "period",
    values_to = "peakFreq1"
  ) %>%
  mutate(
    period = recode(period,
                    base_peakFreq1  = "base",
                    early_peakFreq1 = "early",
                    late_peakFreq1  = "late"
    ),
    period = factor(period, levels = c("base", "early", "late")),
    peakFreq1 = as.numeric(peakFreq1),
    peakFreq1 = if_else(peakFreq1 >= 20 & peakFreq1 <= 50, peakFreq1, NA_real_)
  )

# ---- plot: line per participant (separate lines per session if you have 1+2) ----
p <- ggplot(
  peak1_long,
  aes(x = period, y = peakFreq1,
      group = interaction(participant_ID, session_num),
      linestyle = session_num,
      color = participant_ID, 
      shape = as.factor(session_num))
) +
  scale_color_manual(values = pal, drop = FALSE) +
  geom_line(na.rm = TRUE, alpha = 0.8) +
  geom_point(na.rm = TRUE, size = 6) +
  facet_wrap(~ task) +
  labs(
    x = NULL,
    y = "Peak frequency 1 (Hz)",
    title = "FOOOF peakFreq1 across base/early/late (values outside 20–50 Hz omitted)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  )

p

ggsave(paste0(outPath, "lowGammaPeakAcrossResp.png"), p, width = 10, height = 5, dpi = 300)


# ---- compute breathing-task mean per participant x session, then subtract from cue/O15 ----
breath_mean <- peak1_long %>%
  filter(task == "breathingPreProc") %>%
  group_by(participant_ID, session_num) %>%
  summarise(breath_mean = mean(peakFreq1, na.rm = TRUE), .groups = "drop") %>%
  mutate(breath_mean = if_else(is.nan(breath_mean), NA_real_, breath_mean))

peak1_adj <- peak1_long %>%
  left_join(breath_mean, by = c("participant_ID", "session_num")) %>%
  mutate(peakFreq1_adj = peakFreq1 - breath_mean) %>%
  filter(task %in% c("cueTaskPreProc", "O15preproc"))

# ---- plot adjusted values (cue + O15 only) ----
p <- ggplot(
  peak1_adj,
  aes(x = period, y = peakFreq1_adj,
      group = interaction(participant_ID, session_num),
      linetype = as.factor(session_num),
      color = participant_ID,
      shape = as.factor(session_num))
) +
  scale_color_manual(values = pal, drop = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.6, alpha = 0.6) +
  geom_line(na.rm = TRUE, alpha = 0.8) +
  geom_point(na.rm = TRUE, size = 2) +
  facet_wrap(~ task) +
  labs(
    x = NULL,
    y = "PeakFreq1 (Hz) relative to breathing mean",
    title = "FOOOF peakFreq1 across base/early/late (cue & O15; breathing mean subtracted)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  )

p



# get a score for how much early is rising relative to sides: 
peak1_adj %>% group_by(participant_ID, session_num, task) %>% 
  summarize(freqJump = peakFreq1[2] - mean(peakFreq1[c(1,3)], na.rm = T)) -> peakJump


# -----------------------------
# 1) Compute freqJump per participant x session x task
# -----------------------------
peakJump <- peak1_adj %>%
  mutate(period = factor(period, levels = c("base","early","late"))) %>%
  arrange(participant_ID, session_num, task, period) %>%
  group_by(participant_ID, session_num, task) %>%
  summarise(
    freqJump = {
      v <- peakFreq1_adj
      v[2] - mean(v[c(1,3)], na.rm = TRUE)
    },
    .groups = "drop"
  )

# -----------------------------
# 2) Pull performance per participant x session x task
# -----------------------------
perf_task <- df %>%
  mutate(
    participant_ID = as.character(participant_ID),
    participant_ID = str_replace_all(participant_ID, "DUPI|dupi", "Dupi"),
    participant_ID = str_replace(participant_ID, "^250811_Dupi_NMH_TPB$", "250811_Dupi_NMH_TB"),
    session_num = as.integer(session_num),
    performanceType = as.character(performanceType),
    performance = as.numeric(performance)
  ) %>%
  filter(task %in% c("O15preproc", "cueTaskPreProc")) %>%
  filter((task == "O15preproc" & performanceType == "accuracy") |
           (task == "cueTaskPreProc" & performanceType == "d_prime")) %>%
  group_by(participant_ID, session_num, task) %>%
  summarise(performance = mean(performance, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 3) Join, compute Session2-Session1 deltas, plot delta~delta
# -----------------------------
delta_dat <- peakJump %>%
  inner_join(perf_task, by = c("participant_ID","session_num","task")) %>%
  filter(session_num %in% c(1,2)) %>%
  pivot_wider(
    names_from = session_num,
    values_from = c(freqJump, performance),
    names_prefix = "sess"
  ) 

delta_dat %>%
  mutate(
    d_freqJump = freqJump_sess2 - freqJump_sess1,
    d_perf     = performance_sess2 - performance_sess1
  ) -> delta_dat

# -----------------------------
# 4) Plot for O15 (accuracy deltas)
# -----------------------------
p_o15 <- delta_dat %>%
  filter(task == "O15preproc", !is.na(freqJump_sess1)) %>%
  ggplot(aes(x = freqJump_sess1, y = d_perf, color = participant_ID)) +
  geom_hline(yintercept = 0, linewidth = 0.6, alpha = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.6, alpha = 0.6) +
  geom_point(size = 6, alpha = 0.9) +
  # geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "black") +
  scale_color_manual(values = pal, drop = FALSE) +
  labs(
    x = "Δ freqJump (Session 2 − Session 1)",
    y = "Δ accuracy (Session 2 − Session 1)",
    title = "O15: change in performance vs change in gamma peak jump",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")

p_o15

# -----------------------------
# 5) Plot for cueTask (d' deltas)
# -----------------------------
p_cue <- delta_dat %>%
  filter(task == "cueTaskPreProc", !is.na(freqJump_sess1)) %>%
  ggplot(aes(x = freqJump_sess1, y = d_perf, color = participant_ID)) +
  geom_hline(yintercept = 0, linewidth = 0.6, alpha = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.6, alpha = 0.6) +
  geom_point(size = 6, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "black") +
  scale_color_manual(values = pal, drop = FALSE) +
  labs(
    x = "Δ freqJump (Session 2 − Session 1)",
    y = "Δ d' (Session 2 − Session 1)",
    title = "cueTask: change in performance vs change in gamma peak jump",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")

p_cue


