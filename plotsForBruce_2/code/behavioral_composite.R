# behavioral_composite.R — per-session olfactory scores, composite index,
# SUSTAINED-RECOVERY responder split, and line-based olfaction plots.
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
tdir <- file.path(proj,"out","tables"); fdir <- file.path(proj,"out","figs")
dir.create(fdir, showWarnings=FALSE, recursive=TRUE)
T_RECOVER <- 0.40   # sustained-recovery threshold (~2/3 of control-mean composite)

beh <- read_csv(file.path(tdir,"behavioral_long.csv"), show_col_types=FALSE)
idx <- read_csv(file.path(tdir,"session_index.csv"), show_col_types=FALSE) %>%
  filter(cohort %in% c("Dupi","OBE"), onDisk==1) %>%
  distinct(sessID, cohort, group, participant, sessNum)

w <- beh %>% select(sessID, metric, value) %>%
  pivot_wider(names_from=metric, values_from=value, values_fn=~mean(.x, na.rm=TRUE))
sess <- idx %>% left_join(w, by="sessID")
needcols <- c("cue_d","cue_HR","cue_FA","cue_hits","cue_misses","cue_fas","cue_crs",
              "thresh_none","thresh_low","thresh_high","thresh_low_cal","thresh_high_cal","O15_acc","O15_score")
for(nm in needcols) if(!nm %in% names(sess)) sess[[nm]] <- NA_real_

sess <- sess %>% mutate(
  cs_O15=O15_acc, cs_thresh=thresh_low_cal/400, cs_cue=cue_d/3.5,
  nAll3=(!is.na(cue_d))+(!is.na(thresh_low_cal))+(!is.na(O15_acc)),
  cs_n=(!is.na(cs_O15)&cs_O15!=0)+(!is.na(cs_thresh)&cs_thresh!=0)+(!is.na(cs_cue)&cs_cue!=0),
  composite=ifelse(cs_n>0,(coalesce(cs_O15,0)+coalesce(cs_thresh,0)+coalesce(cs_cue,0))/cs_n,NA_real_),
  composite_all3=ifelse(nAll3==3,composite,NA_real_))
write_csv(sess %>% select(-starts_with("cs_")), file.path(tdir,"session_scores.csv"))

# ================= SUSTAINED-RECOVERY responder split =================
# Responder = final-session composite >= T_RECOVER AND > session-1 composite.
# Non-responder = ended below threshold or declined. Unclassified = single session.
traj <- sess %>% filter(cohort=="Dupi", is.finite(composite)) %>%
  select(participant, sessNum, composite) %>% arrange(participant, sessNum)
resp <- traj %>% group_by(participant) %>%
  summarise(nSess=n_distinct(sessNum),
            firstSess=min(sessNum), lastSess=max(sessNum),
            composite_S1=composite[sessNum==min(sessNum)][1],
            composite_final=composite[sessNum==max(sessNum)][1], .groups="drop") %>%
  mutate(delta=composite_final-composite_S1,
         class=case_when(nSess<2 ~ "unclassified (single session)",
                         composite_final>=T_RECOVER & composite_final>composite_S1 ~ "responder",
                         TRUE ~ "non-responder")) %>%
  arrange(class, desc(composite_final))
write_csv(resp, file.path(tdir,"responder_table.csv"))

classMap <- resp %>% select(participant, class)
sess <- sess %>% left_join(classMap, by="participant") %>%
  mutate(grpColor=case_when(cohort=="OBE" ~ "control",
                            class=="responder" ~ "responder",
                            class=="non-responder" ~ "non-responder", TRUE ~ "unclassified"),
         xpos=ifelse(cohort=="OBE","Control", paste0("S",sessNum)),
         xpos=factor(xpos, levels=c("S1","S2","S3","Control")),
         xnum=ifelse(cohort=="OBE", 4, sessNum))
write_csv(sess %>% select(sessID,cohort,group,participant,sessNum,xpos,xnum,grpColor,class,
                          cue_d,cue_HR,cue_FA,thresh_none,thresh_low,thresh_high,
                          thresh_low_cal,thresh_high_cal,O15_score,O15_acc,composite,composite_all3),
          file.path(tdir,"session_scores_labeled.csv"))

# ================= line-based olfaction plots =================
# distinct color per Dupi participant; solid=responder / dashed=non-responder / dotted=unclassified;
# labeled at last point; controls = grey points (single session). Session means (black) overlaid.
partpal <- setNames(scales::hue_pal()(length(unique(resp$participant))), sort(unique(resp$participant)))
ltymap <- c(responder="solid", `non-responder`="dashed", `unclassified (single session)`="dotted", control="blank")

plotMetric <- function(df, ycol, ylab, fname){
  dd <- df %>% filter(is.finite(.data[[ycol]]))
  if(nrow(dd)==0) return(invisible())
  dup <- dd %>% filter(cohort=="Dupi")
  ctl <- dd %>% filter(cohort=="OBE")
  ms  <- dd %>% group_by(xnum) %>% summarise(m=mean(.data[[ycol]],na.rm=TRUE),
              se=sd(.data[[ycol]],na.rm=TRUE)/sqrt(sum(is.finite(.data[[ycol]]))),
              n=sum(is.finite(.data[[ycol]])), .groups="drop")
  labs <- dup %>% group_by(participant) %>% filter(xnum==max(xnum)) %>% ungroup()
  p <- ggplot() +
    geom_line(data=dup, aes(xnum, .data[[ycol]], group=participant, color=participant, linetype=class), linewidth=.75, alpha=.9) +
    geom_point(data=dup, aes(xnum, .data[[ycol]], color=participant), size=1.9) +
    geom_point(data=ctl, aes(xnum, .data[[ycol]]), color="grey45", size=1.9, alpha=.6, position=position_jitter(width=.06,height=0)) +
    geom_errorbar(data=ms, aes(xnum, ymin=m-se, ymax=m+se), width=.14, linewidth=.5, color="black") +
    geom_point(data=ms, aes(xnum, m), shape=95, size=6, color="black") +
    geom_text(data=ms, aes(xnum, m, label=paste0("n=",n)), vjust=2.1, size=2.6, color="grey25") +
    geom_text(data=labs, aes(xnum, .data[[ycol]], label=participant, color=participant), size=2.6, nudge_x=.12, hjust=0, show.legend=FALSE) +
    scale_color_manual(values=partpal, guide="none") +
    scale_linetype_manual(values=ltymap, name="Dupi class",
                          breaks=c("responder","non-responder","unclassified (single session)")) +
    scale_x_continuous(breaks=1:4, labels=c("S1","S2","S3","Control"), limits=c(.8,4.4)) +
    labs(x=NULL, y=ylab, title=ylab) + theme_bw(base_size=12) +
    theme(legend.position="bottom")
  ggsave(file.path(fdir,fname), p, width=6.5, height=4.4, dpi=130)
}
plotMetric(sess,"cue_d","cueTask d'","beh_cue_d.png")
plotMetric(sess,"thresh_low_cal","threshTask low-PEA (air-calibrated)","beh_thresh_lowcal.png")
plotMetric(sess,"thresh_high_cal","threshTask high-PEA (air-calibrated)","beh_thresh_highcal.png")
plotMetric(sess,"O15_acc","O15 identification accuracy","beh_o15_acc.png")
plotMetric(sess,"composite","Olfaction composite index","beh_composite.png")
plotMetric(sess,"composite_all3","Olfaction composite (all-3 tasks)","beh_composite_all3.png")

cat("=== Responder / non-responder split (sustained recovery, threshold",T_RECOVER,") ===\n")
print(as.data.frame(resp %>% mutate(across(where(is.numeric),~round(.x,3)))), row.names=FALSE)
cat("\n=== per-xpos composite means ===\n")
print(sess %>% group_by(xpos) %>% summarise(n=sum(is.finite(composite)),mean_composite=round(mean(composite,na.rm=TRUE),3)) %>% as.data.frame(), row.names=FALSE)
cat("\nwrote session_scores.csv, session_scores_labeled.csv, responder_table.csv, beh_*.png\n")
