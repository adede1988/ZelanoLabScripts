# taskcmp_D_analysis.R — group boxplots + paired stats for the timing/late-window
# metrics, the within-participant Dupi S1->S2 change, and its correlation with the
# olfactory-composite change. Bands: theta (4-8) & gamma (25-58). Metrics:
#   peakLatMs (time-to-peak of band-z), earlyBandMax (100-500 ms), lateBandMax (3000-4000 ms).
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2); library(purrr)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
T <- file.path(proj,"out","tables"); Fg <- file.path(proj,"out","figs","taskcmp")
bm <- read_csv(file.path(T,"taskcmp_bandmetrics.csv"), show_col_types=FALSE)
sc <- read_csv(file.path(T,"session_scores.csv"), show_col_types=FALSE) |> select(participant, sessNum, composite)
mets <- c("peakLatMs","earlyBandMax","lateBandMax"); bands <- c("theta","gamma")
grp3 <- c("Control","DupiS1","DupiS2")
bm <- bm |> mutate(group=factor(group, levels=c("Control","DupiS1","DupiS2","DupiS3")),
                   condition=factor(condition, levels=c("audiobook","focus")))
pair_stat <- function(a,f){ ok<-is.finite(a)&is.finite(f); a<-a[ok]; f<-f[ok]; d<-a-f; n<-length(d)
  if(n<3) return(c(n=n,dz=NA,t_p=NA,w_p=NA))
  c(n=n, dz=mean(d)/sd(d), t_p=tryCatch(t.test(a,f,paired=TRUE)$p.value,error=function(e)NA),
    w_p=tryCatch(wilcox.test(a,f,paired=TRUE)$p.value,error=function(e)NA)) }

# ---------- (1) group boxplots + paired audiobook-vs-focus stats ----------
d3 <- bm |> filter(group %in% grp3)
bl <- d3 |> pivot_longer(all_of(mets), names_to="metric", values_to="val") |>
  mutate(metric=factor(metric, levels=mets))
pb <- ggplot(bl, aes(group, val, fill=condition)) +
  geom_boxplot(outlier.size=0.6, width=0.65, position=position_dodge(0.7)) +
  facet_grid(metric ~ band, scales="free_y", switch="y") +
  scale_fill_manual(values=c(audiobook="#3366cc", focus="#d9480f")) +
  labs(x=NULL, y=NULL, title="Timing / late-window gamma & theta metrics: audiobook vs focus by group") +
  theme_bw(base_size=13) + theme(legend.position="top", strip.placement="outside",
    strip.background=element_rect(fill="grey92"), plot.title=element_text(face="bold", size=14))
ggsave(file.path(Fg,"boxes_group_condition.png"), pb, width=11, height=10, dpi=120)

BX <- list()
for(b in bands) for(mt in mets) for(g in grp3){
  s <- d3 |> filter(band==b, group==g)
  wide <- s |> select(participant, condition, val=all_of(mt)) |>
    group_by(participant, condition) |> summarise(val=mean(val,na.rm=TRUE), .groups="drop") |>
    pivot_wider(names_from=condition, values_from=val)
  st <- pair_stat(wide$audiobook, wide$focus)
  BX[[length(BX)+1]] <- tibble(band=b, metric=mt, group=g, n=st["n"],
     mean_audiobook=mean(wide$audiobook,na.rm=TRUE), mean_focus=mean(wide$focus,na.rm=TRUE),
     dz=st["dz"], t_p=st["t_p"], wilcox_p=st["w_p"]) }
BXt <- bind_rows(BX); write_csv(BXt, file.path(T,"taskcmp_D_boxstats.csv"))
cat("=== paired audiobook-vs-focus by group (focus earlier => peakLat +; focus higher => bandMax -) ===\n")
print(as.data.frame(BXt |> mutate(across(where(is.numeric),~round(.x,3)))), row.names=FALSE)

# ---------- (2) contrasts per session, then Dupi S1->S2 within-participant change ----------
# contrast sign: peakLat -> audiobook - focus (+ = focus earlier); bandMax -> focus - audiobook (+ = focus higher)
contr <- bm |> select(participant, group, sessNum, cohort, condition, band, all_of(mets)) |>
  pivot_longer(all_of(mets), names_to="metric", values_to="val") |>
  pivot_wider(names_from=condition, values_from=val) |>
  mutate(contrast = ifelse(metric=="peakLatMs", audiobook - focus, focus - audiobook))
# Dupi sessions with both S1 and S2 (per participant, per band, per metric)
dupi <- contr |> filter(cohort=="Dupi", sessNum %in% c(1,2))
s1s2 <- dupi |> select(participant, band, metric, sessNum, contrast) |>
  pivot_wider(names_from=sessNum, values_from=contrast, names_prefix="S") |>
  filter(is.finite(S1) & is.finite(S2)) |> mutate(delta = S2 - S1)
write_csv(s1s2, file.path(T,"taskcmp_D_s1s2.csv"))
S1S2 <- s1s2 |> group_by(band, metric) |>
  summarise(n=n(), mean_S1=mean(S1), mean_S2=mean(S2), mean_delta=mean(delta),
            dz=mean(delta)/sd(delta),
            t_p=tryCatch(t.test(S2,S1,paired=TRUE)$p.value,error=function(e)NA),
            wilcox_p=tryCatch(wilcox.test(S2,S1,paired=TRUE)$p.value,error=function(e)NA), .groups="drop")
write_csv(S1S2, file.path(T,"taskcmp_D_s1s2_stats.csv"))
cat("\n=== Dupi S1->S2 change in the audiobook-vs-focus contrast (paired within participant) ===\n")
print(as.data.frame(S1S2 |> mutate(across(where(is.numeric),~round(.x,3)))), row.names=FALSE)

# S1->S2 paired plot (contrast at S1 vs S2, participant lines)
s1s2L <- s1s2 |> pivot_longer(c(S1,S2), names_to="session", values_to="contrast")
ps <- ggplot(s1s2L, aes(session, contrast, group=participant)) +
  geom_line(color="grey60") + geom_point(color="grey40", size=1.8) +
  stat_summary(aes(group=1), fun=mean, geom="line", color="black", linewidth=1.4) +
  stat_summary(aes(group=1), fun=mean, geom="point", color="black", size=3) +
  facet_grid(metric ~ band, scales="free_y", switch="y") + geom_hline(yintercept=0, linetype=3) +
  labs(x=NULL, y="audiobook-vs-focus contrast (+ = focus earlier / higher)",
       title="Dupi within-participant S1\u2192S2 change in the focus-vs-audiobook effect") +
  theme_bw(base_size=13) + theme(strip.placement="outside", strip.background=element_rect(fill="grey92"),
    plot.title=element_text(face="bold", size=14))
ggsave(file.path(Fg,"s1s2_change.png"), ps, width=10, height=9, dpi=120)

# ---------- (3) correlate S1->S2 neural change with olfactory change ----------
dOlf <- sc |> filter(sessNum %in% c(1,2)) |>
  group_by(participant, sessNum) |> summarise(composite=mean(composite, na.rm=TRUE), .groups="drop") |>
  pivot_wider(names_from=sessNum, values_from=composite, names_prefix="olf_S") |>
  filter(is.finite(olf_S1) & is.finite(olf_S2)) |> mutate(dOlf = olf_S2 - olf_S1)
cc <- s1s2 |> left_join(dOlf |> select(participant, dOlf), by="participant") |> filter(is.finite(dOlf))
COR <- cc |> group_by(band, metric) |>
  summarise(n=n(),
            r_pearson=tryCatch(cor(delta, dOlf), error=function(e)NA),
            p_pearson=tryCatch(cor.test(delta, dOlf)$p.value, error=function(e)NA),
            r_spear=tryCatch(cor(delta, dOlf, method="spearman"), error=function(e)NA),
            p_spear=tryCatch(suppressWarnings(cor.test(delta, dOlf, method="spearman")$p.value), error=function(e)NA),
            .groups="drop")
write_csv(COR, file.path(T,"taskcmp_D_corr.csv"))
cat("\n=== corr: S1->S2 neural contrast change vs S1->S2 olfaction (composite) change ===\n")
print(as.data.frame(COR |> mutate(across(where(is.numeric),~round(.x,3)))), row.names=FALSE)

lab_layer <- if(requireNamespace("ggrepel", quietly=TRUE)) ggrepel::geom_text_repel(aes(label=participant), size=3, max.overlaps=20) else geom_text(aes(label=participant), size=3, vjust=-0.6)
pc <- ggplot(cc, aes(dOlf, delta)) + geom_point(size=2.4, color="#8856a7") +
  geom_smooth(method="lm", se=FALSE, color="grey40", linewidth=0.8) + lab_layer +
  facet_grid(metric ~ band, scales="free", switch="y") +
  labs(x="\u0394 olfactory composite (S2 \u2212 S1)", y="\u0394 focus-vs-audiobook contrast (S2 \u2212 S1)",
       title="Does the S1\u2192S2 change in the focus effect track olfactory change?") +
  theme_bw(base_size=13) + theme(strip.placement="outside", strip.background=element_rect(fill="grey92"),
    plot.title=element_text(face="bold", size=13))
ggsave(file.path(Fg,"corr_olfaction.png"), pc, width=10, height=9, dpi=120)
cat("\nwrote taskcmp_D_boxstats.csv, taskcmp_D_s1s2.csv, taskcmp_D_s1s2_stats.csv, taskcmp_D_corr.csv + 3 figs\n")
