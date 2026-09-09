# breathing_deficit.R — audiobook & focusedBreathing: which per-breath gamma measures
# separate Dupi from control AND track across sessions. Writes a stats table + a
# 4-metric x 2-task figure (S1/S2/S3/Control). Self-contained; sheet/data driven.
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr); library(ggplot2)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if(length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
T <- file.path(proj,"out","tables"); Fg <- file.path(proj,"out","figs")
dir.create(Fg, showWarnings=FALSE, recursive=TRUE)
sess <- read_csv(file.path(T,"gamma_session_level.csv"), show_col_types=FALSE)
resp <- read_csv(file.path(T,"responder_table.csv"), show_col_types=FALSE) |> transmute(participant, cls=class)

metric_cols <- names(sess)[str_detect(names(sess), "__(mean|var)$")]
metric_cols <- metric_cols[map_lgl(metric_cols, ~is.numeric(sess[[.x]]))]
pooled_d <- function(a,b){ a<-a[is.finite(a)]; b<-b[is.finite(b)]; if(length(a)<2||length(b)<2) return(NA)
  sp<-sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2)); if(!is.finite(sp)||sp==0) return(NA); (mean(a)-mean(b))/sp }

# ---- full dual-criterion sweep (separation + tracking) per task ----
sweep_task <- function(tk){
  d <- sess |> filter(task==tk); ctrlD <- d |> filter(cohort=="Control"); dupiD <- d |> filter(cohort=="Dupi", sessNum %in% 1:3)
  map_dfr(metric_cols, function(mc){
    cv <- ctrlD[[mc]]; cv<-cv[is.finite(cv)]
    dd <- dupiD |> select(participant, sessNum, val=all_of(mc)) |> filter(is.finite(val))
    if(length(cv)<3 || nrow(dd)<6) return(NULL)
    pm <- dd |> group_by(participant) |> summarise(v=mean(val), .groups="drop") |> pull(v)
    sm <- dd |> group_by(sessNum) |> summarise(m=mean(val), .groups="drop"); g <- function(s){x<-sm$m[sm$sessNum==s]; if(length(x))x else NA}
    d_s1 <- pooled_d(dd$val[dd$sessNum==1], cv); d_s3 <- pooled_d(dd$val[dd$sessNum==3], cv)
    tibble(task=tk, metric=mc,
           ctrl_mean=mean(cv), ctrl_sd=sd(cv), dupi_sd=sd(dd$val),
           d_sep=pooled_d(pm, cv), p_sep=tryCatch(t.test(pm,cv)$p.value,error=function(e)NA),
           rho=tryCatch(suppressWarnings(cor(dd$val,dd$sessNum,method="spearman")),error=function(e)NA),
           p_rho=tryCatch(suppressWarnings(cor.test(dd$val,dd$sessNum,method="spearman")$p.value),error=function(e)NA),
           S1=g(1), S2=g(2), S3=g(3), d_S1_ctrl=d_s1, d_S3_ctrl=d_s3, recover_gain=abs(d_s1)-abs(d_s3))
  }) |> mutate(p_sep_fdr=p.adjust(p_sep,"BH"), p_rho_fdr=p.adjust(p_rho,"BH"))
}
allstats <- bind_rows(lapply(c("audiobook","focusedBreathing"), sweep_task))
write_csv(allstats, file.path(T,"breathing_deficit_stats.csv"))

# ---- highlighted-measure summary (for the report table) ----
hi <- c("dutyCycle__mean","timeAboveMs__mean","gammaPeakPresent__mean","chirpSlope__mean")
summ <- allstats |> filter(metric %in% hi) |>
  transmute(task, metric,
            ctrl=sprintf("%.3g\u00b1%.2g", ctrl_mean, ctrl_sd),
            S1=signif(S1,3), S2=signif(S2,3), S3=signif(S3,3),
            d_sep=round(d_sep,2), p_sep=signif(p_sep,2), p_sep_fdr=signif(p_sep_fdr,2),
            rho=round(rho,2), p_rho=signif(p_rho,2),
            SDratio_ctrl_dupi=round(ctrl_sd/dupi_sd,2),
            chirp_S1toS2=ifelse(metric=="chirpSlope__mean", round(S2-S1,2), NA),
            chirp_S2toS3=ifelse(metric=="chirpSlope__mean", round(S3-S2,2), NA))
write_csv(summ, file.path(T,"breathing_deficit_summary.csv"))
cat("=== highlighted-measure summary ===\n"); print(as.data.frame(summ), row.names=FALSE)

# ---- figure ----
labs <- c(dutyCycle__mean="Bursting: duty cycle\n(fraction of breath, ridge z>3)",
          timeAboveMs__mean="Bursting: time above thr.\n(ms of gamma burst / breath)",
          gammaPeakPresent__mean="Periodic peak rate\n(frac. breaths w/ FOOOF gamma bump)",
          chirpSlope__mean="Chirp slope\n(ridge Hz/s, 0-500 ms)")
d <- sess |> filter(task %in% c("audiobook","focusedBreathing")) |> left_join(resp, by="participant") |>
  mutate(grp=case_when(cohort=="Control"~"control", cls=="responder"~"responder", cls=="non-responder"~"non-responder", TRUE~"unclassified"),
         xpos=factor(ifelse(cohort=="Control","Control",paste0("S",sessNum)), levels=c("S1","S2","S3","Control")),
         task=factor(task, levels=c("audiobook","focusedBreathing"), labels=c("Audiobook","Focused breathing")))
long <- d |> select(task, xpos, grp, all_of(hi)) |> pivot_longer(all_of(hi), names_to="metric", values_to="val") |>
  filter(is.finite(val)) |> mutate(metric=factor(labs[metric], levels=labs))
ms <- long |> group_by(task, metric, xpos) |> summarise(m=mean(val), se=sd(val)/sqrt(sum(is.finite(val))), .groups="drop")
pal <- c(responder="#1b9e77", `non-responder`="#d95f02", control="#7570b3", unclassified="#999999")
p <- ggplot(long, aes(xpos, val)) +
  geom_vline(xintercept=3.5, linetype="dashed", color="grey60") +
  geom_jitter(aes(color=grp), width=.14, height=0, size=2.1, alpha=.75) +
  geom_errorbar(data=ms, aes(xpos, ymin=m-se, ymax=m+se), inherit.aes=FALSE, width=.28, linewidth=.7) +
  geom_point(data=ms, aes(xpos, m), inherit.aes=FALSE, shape=95, size=9) +
  geom_line(data=ms |> filter(xpos!="Control"), aes(as.numeric(xpos), m), inherit.aes=FALSE, linewidth=.6, color="grey20") +
  facet_grid(metric ~ task, scales="free_y", switch="y") + scale_color_manual(values=pal, name="") +
  labs(x=NULL, y=NULL, title="Gamma bursting / periodic deficit & chirp slope \u2014 breathing tasks",
       subtitle="Point = one session; bar = mean \u00b1 SE; black line = Dupi S1\u2192S2\u2192S3; dashed = Dupi vs OBE control") +
  theme_bw(base_size=15) +
  theme(legend.position="bottom", legend.text=element_text(size=13), strip.text.y.left=element_text(angle=0, size=12.5, face="bold"),
        strip.text.x=element_text(size=15, face="bold"), strip.placement="outside", strip.background.y=element_rect(fill="grey95"),
        axis.text=element_text(size=13), plot.title=element_text(size=18, face="bold"),
        plot.subtitle=element_text(size=12.5, color="grey30"), panel.spacing=unit(0.9,"lines"))
ggsave(file.path(Fg,"breathing_deficit_tracking.png"), p, width=12.5, height=13.5, dpi=135)
cat("\nwrote out/figs/breathing_deficit_tracking.png, out/tables/breathing_deficit_stats.csv, breathing_deficit_summary.csv\n")
