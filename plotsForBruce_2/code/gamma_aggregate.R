# gamma_aggregate.R — per-breath gamma metrics -> session-level -> group plots + goodness ranking.
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2); library(purrr)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
tdir <- file.path(proj,"out","tables"); gdir <- file.path(proj,"out","gamma")
pbdir <- file.path(gdir,"perbreath"); fdir <- file.path(proj,"out","figs","gamma")
dir.create(fdir, showWarnings=FALSE, recursive=TRUE)

# ---- load all per-breath csvs ----
files <- list.files(pbdir, pattern="\\.csv$", full.names=TRUE)
stopifnot(length(files)>0)
pb <- files |> map(~suppressMessages(read_csv(.x, show_col_types=FALSE))) |> list_rbind()
cat("loaded", nrow(pb), "breaths across", length(files), "recordings\n")

id_cols <- c("sessID","task","cohort","group","participant","sessNum","breathIdx","onsetSample","goodBreath","nEv")
metric_cols <- setdiff(names(pb), id_cols)
metric_cols <- metric_cols[map_lgl(metric_cols, ~is.numeric(pb[[.x]]))]

resp <- read_csv(file.path(tdir,"responder_table.csv"), show_col_types=FALSE) |>
  transmute(participant, class)

# ---- QC: use goodBreath where available (breathing); sniff tasks goodBreath=1 ----
pbq <- pb |> filter(is.na(goodBreath) | goodBreath==1)

# ---- session-level aggregation: mean AND variance of each metric across breaths ----
sess_mean <- pbq |> group_by(sessID, task, cohort, participant, sessNum) |>
  summarise(nBreath=n(), across(all_of(metric_cols), ~mean(.x, na.rm=TRUE), .names="{.col}__mean"), .groups="drop")
sess_var  <- pbq |> group_by(sessID, task) |>
  summarise(across(all_of(metric_cols), ~var(.x, na.rm=TRUE), .names="{.col}__var"), .groups="drop")
sess <- sess_mean |> left_join(sess_var, by=c("sessID","task")) |>
  left_join(resp, by="participant") |>
  mutate(cohort=ifelse(cohort=="OBE","Control",cohort),
         grpColor = case_when(cohort=="Control"~"control",
                              class=="responder"~"responder",
                              class=="non-responder"~"non-responder",
                              TRUE~"unclassified"),
         xpos = ifelse(cohort=="Control","Control", paste0("S",sessNum)),
         xpos = factor(xpos, levels=c("S1","S2","S3","Control")))
write_csv(sess, file.path(tdir,"gamma_session_level.csv"))

# ================= goodness ranking =================
agg_metric_cols <- names(sess)[str_detect(names(sess), "__(mean|var)$")]
d_between <- function(a,b){ a<-a[is.finite(a)]; b<-b[is.finite(b)]
  if(length(a)<2||length(b)<2) return(NA_real_)
  sp <- sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2))
  if(!is.finite(sp)||sp==0) return(NA_real_); (mean(a)-mean(b))/sp }

goodness <- list()
for (tk in unique(sess$task)) {
  st <- sess |> filter(task==tk)
  ctrl <- st |> filter(cohort=="Control")
  dupi <- st |> filter(cohort=="Dupi")
  rr   <- dupi |> filter(grpColor=="responder")
  nn   <- dupi |> filter(grpColor=="non-responder")
  for (mc in agg_metric_cols) {
    cv_ctrl <- { v<-ctrl[[mc]]; v<-v[is.finite(v)]
      if(length(v)<2) NA_real_ else sd(v)/(abs(mean(v))+1e-9) }
    sd_ctrl <- { v<-ctrl[[mc]]; v<-v[is.finite(v)]; if(length(v)<2) NA_real_ else sd(v) }
    sep_cd <- d_between(dupi[[mc]], ctrl[[mc]])
    sep_rn <- d_between(rr[[mc]], nn[[mc]])
    # recovery: Spearman corr of metric with sessNum within Dupi (and within responders)
    rec_all <- { x<-dupi$sessNum; y<-dupi[[mc]]; ok<-is.finite(x)&is.finite(y)
      if(sum(ok)<4) NA_real_ else suppressWarnings(cor(x[ok],y[ok],method="spearman")) }
    rec_resp <- { x<-rr$sessNum; y<-rr[[mc]]; ok<-is.finite(x)&is.finite(y)
      if(sum(ok)<4) NA_real_ else suppressWarnings(cor(x[ok],y[ok],method="spearman")) }
    goodness[[length(goodness)+1]] <- tibble(
      task=tk, metric=mc, n_ctrl=sum(is.finite(ctrl[[mc]])), n_dupi=sum(is.finite(dupi[[mc]])),
      control_mean=mean(ctrl[[mc]],na.rm=TRUE), control_SD=sd_ctrl, control_CV=cv_ctrl,
      sep_control_dupi=sep_cd, sep_resp_nonresp=sep_rn,
      recovery_rho_dupi=rec_all, recovery_rho_resp=rec_resp,
      max_separation=max(abs(c(sep_cd,sep_rn,rec_all,rec_resp)),na.rm=TRUE))
  }
}
G <- list_rbind(goodness)
# rank: low control CV AND high max separation
G <- G |> mutate(
  control_CV_rank = rank(abs(control_CV), na.last="keep"),
  # goodness score: reward large separation, penalize control variability
  goodness = abs(max_separation) / (1 + abs(control_CV))) |>
  arrange(desc(goodness))
write_csv(G, file.path(tdir,"gamma_goodness_ranking.csv"))

# top measures overall + per criterion
top_overall <- G |> filter(is.finite(goodness)) |> slice_head(n=40)
write_csv(top_overall, file.path(tdir,"gamma_top_measures.csv"))

# ================= plots: key + top metrics =================
pal <- c(responder="#1b9e77", `non-responder`="#d95f02", control="#7570b3", unclassified="#999999")
plot_metric_faceted <- function(mc, fname){
  d <- sess |> select(task, xpos, grpColor, val=all_of(mc)) |> filter(is.finite(val))
  if(nrow(d)==0) return(invisible())
  ms <- d |> group_by(task,xpos) |>
    summarise(m=mean(val,na.rm=TRUE), se=sd(val,na.rm=TRUE)/sqrt(sum(is.finite(val))),
              n=sum(is.finite(val)), .groups="drop")
  p <- ggplot(d, aes(xpos, val)) +
    geom_jitter(aes(color=grpColor), width=.13, height=0, size=1.4, alpha=.7) +
    geom_errorbar(data=ms, aes(xpos, ymin=m-se, ymax=m+se), inherit.aes=FALSE, width=.2) +
    geom_point(data=ms, aes(xpos, m), inherit.aes=FALSE, shape=95, size=4) +
    facet_wrap(~task, scales="free_y", nrow=1) +
    scale_color_manual(values=pal, name="") +
    labs(x=NULL, y=mc, title=mc) + theme_bw(base_size=10) +
    theme(legend.position="bottom", axis.text.x=element_text(size=7))
  ggsave(file.path(fdir,fname), p, width=13, height=3.4, dpi=120)
}
key_metrics <- c("peakZ__mean","peakLatMs__mean","peakFreq__mean","timeAboveMs__mean",
  "dutyCycle__mean","nBursts__mean","anyBurst__mean","gammaBumpDb__mean","gammaPeakPresent__mean",
  "apExp__mean","chirpSlope__mean","freqSpan__mean","freqJitter__mean",
  "w2_rpowZ__mean","w2_rfreqPW__mean","w3_rpowZ__mean","p1_rfreqPW__mean","p1_rpowZ__mean",
  "peakZ__var","w2_rfreqPW__var","peakFreq__var")
key_metrics <- intersect(key_metrics, names(sess))
for (mc in key_metrics) plot_metric_faceted(mc, paste0("gm_", str_replace_all(mc,"[^A-Za-z0-9]","_"), ".png"))
# top-ranked measures too
for (mc in unique(top_overall$metric)[1:min(20,nrow(top_overall))])
  if(mc %in% names(sess)) plot_metric_faceted(mc, paste0("top_", str_replace_all(mc,"[^A-Za-z0-9]","_"), ".png"))

cat("\n=== TOP 20 gamma measures (low control CV x high separation) ===\n")
print(as.data.frame(top_overall |> select(task,metric,control_CV,sep_control_dupi,sep_resp_nonresp,recovery_rho_resp,goodness) |> slice_head(n=20)), row.names=FALSE)
cat("\nwrote gamma_session_level.csv, gamma_goodness_ranking.csv, gamma_top_measures.csv, and figs/gamma/*.png\n")
