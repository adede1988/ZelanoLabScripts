# grant_freqJitter_interaction.R — within-trial gamma-frequency modulation, audiobook vs
# focused breathing, Control vs Dupi session-1. freqJitter = SD of the tracked 25-58 Hz
# ridge frequency over 0-1000 ms post-inhale (per breath; session mean). Paired condition
# effect per group + group x condition interaction. Writes out/grant_freqJitter_interaction.png.
# AD_2 (72 h fasted) excluded from Control. Self-contained; data driven.
suppressWarnings(suppressMessages({library(dplyr);library(readr);library(tidyr);library(purrr);library(ggplot2)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if(length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
pbDir <- file.path(proj,"out","gamma","perbreath")
pb <- list.files(pbDir, pattern="(audiobook|focusedBreathing)[.]csv$", full.names=TRUE)
d <- pb |> map(~suppressMessages(read_csv(.x, show_col_types=FALSE))) |> list_rbind() |>
  filter(is.na(goodBreath)|goodBreath==1, sessID != "260326_OBE_NWU_AD_2") |>
  mutate(grp=case_when(cohort=="OBE" ~ "Control",
                       cohort=="Dupi" & sessNum==1 ~ "Dupi S1", TRUE ~ NA_character_)) |>
  filter(!is.na(grp))
sess <- d |> group_by(grp, sessID, task) |>
  summarise(freqJitter=mean(freqJitter,na.rm=TRUE), .groups="drop") |>
  mutate(cond=recode(task, audiobook="audiobook", focusedBreathing="focused"),
         cond=factor(cond, levels=c("audiobook","focused")),
         grp=factor(grp, levels=c("Control","Dupi S1")))
both <- sess |> group_by(sessID) |> filter(n()==2) |> ungroup()   # sessions with both conditions
gm <- both |> group_by(grp,cond) |> summarise(m=mean(freqJitter), se=sd(freqJitter)/sqrt(n()), .groups="drop")
lab <- both |> select(grp,sessID,cond,freqJitter) |> pivot_wider(names_from=cond,values_from=freqJitter) |>
  group_by(grp) |> summarise(n=n(), dz=mean(focused-audiobook)/sd(focused-audiobook),
                             p=t.test(focused,audiobook,paired=TRUE)$p.value,
                             y=max(both$freqJitter), .groups="drop") |>
  mutate(txt=sprintf("n=%d, dz=%.2f, p=%.3f", n, dz, p))
# group x condition interaction. Condition is WITHIN-session, so the correct test is a
# two-sample t on the per-session condition effect (focused - audiobook) across groups; the
# ordinary lm(~grp*cond) is invalid here (it treats a session's two conditions as independent).
delt <- both |> select(grp,sessID,cond,freqJitter) |> pivot_wider(names_from=cond,values_from=freqJitter) |>
  mutate(delta=focused-audiobook)
tt <- t.test(delta ~ grp, data=delt)                    # Welch, Control vs Dupi S1
pint <- tt$p.value
cat(sprintf("interaction (Welch t on per-session delta): t=%.2f df=%.1f p=%.4f\n", tt$statistic, tt$parameter, pint))
print(as.data.frame(lab[,c("grp","n","dz","p")]))
p <- ggplot(both, aes(cond, freqJitter, group=sessID)) +
  geom_line(color="grey70", linewidth=.7) +
  geom_point(color="grey45", size=2.2) +
  geom_line(data=gm, aes(cond,m,group=1), inherit.aes=FALSE, color="firebrick", linewidth=2) +
  geom_errorbar(data=gm, aes(cond,ymin=m-se,ymax=m+se), inherit.aes=FALSE, color="firebrick", width=.15, linewidth=1) +
  geom_point(data=gm, aes(cond,m), inherit.aes=FALSE, color="firebrick", size=4) +
  geom_text(data=lab, aes(x=1.5,y=y*1.005,label=txt), inherit.aes=FALSE, size=6, fontface="bold") +
  facet_wrap(~grp) +
  labs(x=NULL, y="Ridge-frequency SD (Hz)  [freqJitter]",
       title="Within-trial gamma frequency modulation: audiobook vs focused breathing",
       subtitle=sprintf("Each line = one session (mean per breath); red = group mean \u00b1 SE. Interaction p = %.4f.", pint)) +
  theme_bw(base_size=20) +
  theme(strip.text=element_text(size=22,face="bold"), plot.title=element_text(size=22,face="bold"),
        plot.subtitle=element_text(size=17,color="grey30"), axis.text=element_text(size=18))
ggsave(file.path(proj,"out","grant_freqJitter_interaction.png"), p, width=12, height=6.5, dpi=150)
cat("wrote out/grant_freqJitter_interaction.png\n")
