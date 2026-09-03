# gamma_review_response.R — analyses answering the Zelano/Voytek/Cohen reviews:
#  (1) peak-gated re-ranking (only breaths with a FOOOF gamma peak)
#  (2) airflow-regressed re-ranking (residualize on inhale volume/duration/flow)
#  (3) aperiodic exponent/offset per group + 1/f-controlled separation
#  (4) honest goodness v2: rank on outcome-blind control<->Dupi, bootstrap CIs, flags
#  (5) whole-window z vs pre-inhale baseline z concordance
#  (6) leave-one-participant-out for the responder<->non-responder axis
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr); library(ggplot2)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
tdir <- file.path(proj,"out","tables"); pbdir <- file.path(proj,"out","gamma","perbreath"); fdir <- file.path(proj,"out","figs","gamma")
dir.create(fdir, showWarnings=FALSE, recursive=TRUE)
set.seed(1)

pb <- list.files(pbdir, pattern="\\.csv$", full.names=TRUE) |>
  map(~suppressMessages(read_csv(.x, show_col_types=FALSE))) |> list_rbind()
resp <- read_csv(file.path(tdir,"responder_table.csv"), show_col_types=FALSE) |> transmute(participant, class)
pb <- pb |> filter(is.na(goodBreath) | goodBreath==1)

id_cols <- c("sessID","task","cohort","group","participant","sessNum","breathIdx","onsetSample","goodBreath","nEv")
metric_cols <- setdiff(names(pb), id_cols); metric_cols <- metric_cols[map_lgl(metric_cols, ~is.numeric(pb[[.x]]))]
covs <- c("inhaleVolume","inhaleDuration","peakInspFlow")

label_sessions <- function(df){
  df |> left_join(resp, by="participant") |>
    mutate(cohort2=ifelse(cohort=="OBE","Control",cohort),
           grp=case_when(cohort2=="Control"~"control", class=="responder"~"responder",
                         class=="non-responder"~"non-responder", TRUE~"unclassified"))
}
d_ci <- function(a,b,nboot=2000){
  a<-a[is.finite(a)]; b<-b[is.finite(b)]
  if(length(a)<3||length(b)<3) return(c(d=NA,lo=NA,hi=NA))
  dfun<-function(x,y){sp<-sqrt(((length(x)-1)*var(x)+(length(y)-1)*var(y))/(length(x)+length(y)-2)); if(!is.finite(sp)||sp==0) return(NA); (mean(x)-mean(y))/sp}
  d0<-dfun(a,b); bs<-replicate(nboot, dfun(sample(a,replace=TRUE), sample(b,replace=TRUE)))
  c(d=d0, lo=quantile(bs,.025,na.rm=TRUE), hi=quantile(bs,.975,na.rm=TRUE))
}
sess_agg <- function(df){
  df |> group_by(sessID,task,cohort,participant,sessNum) |>
    summarise(across(all_of(intersect(metric_cols,names(df))), ~mean(.x,na.rm=TRUE)), .groups="drop") |>
    label_sessions()
}

# ---------- ranking builder (returns per task x metric: control CV, seps, CIs, flags) ----------
build_rank <- function(sessdf, tag){
  out <- list()
  mm <- intersect(metric_cols, names(sessdf))
  for(tk in unique(sessdf$task)){
    st <- sessdf |> filter(task==tk)
    ct <- st |> filter(cohort2=="Control"); du <- st |> filter(cohort2=="Dupi")
    rr <- du |> filter(grp=="responder"); nn <- du |> filter(grp=="non-responder")
    for(mc in mm){
      cv <- { v<-ct[[mc]]; v<-v[is.finite(v)]; if(length(v)<2) NA else sd(v)/(abs(mean(v))+1e-9) }
      cd <- d_ci(du[[mc]], ct[[mc]]); rn <- d_ci(rr[[mc]], nn[[mc]])
      rec <- { x<-du$sessNum; y<-du[[mc]]; ok<-is.finite(x)&is.finite(y); if(sum(ok)<5) NA else suppressWarnings(cor(x[ok],y[ok],method="spearman")) }
      out[[length(out)+1]] <- tibble(analysis=tag, task=tk, metric=mc,
        control_CV=cv, n_ctrl=sum(is.finite(ct[[mc]])), n_resp=sum(is.finite(rr[[mc]])), n_nonresp=sum(is.finite(nn[[mc]])),
        sep_cd=cd["d"], sep_cd_lo=cd["lo"], sep_cd_hi=cd["hi"], cd_sig=as.integer(is.finite(cd["lo"])&is.finite(cd["hi"])&(cd["lo"]*cd["hi"]>0)),
        sep_rn=rn["d"], recovery_rho=rec,
        isRidgeFreq=as.integer(str_detect(mc,"rfreq")),
        smallNonResp=as.integer(sum(is.finite(nn[[mc]]))<3))
    }
  }
  list_rbind(out)
}

sessAll <- sess_agg(pb)
rank_all <- build_rank(sessAll, "all_breaths")

# (1) peak-gated: only breaths with a FOOOF gamma peak
sessPk <- sess_agg(pb |> filter(gammaPeakPresent==1))
rank_pk <- build_rank(sessPk, "peak_gated")

# (2) airflow-regressed: residualize each metric on covariates (per task, per-breath pooled)
resid_airflow <- function(df){
  df2 <- df
  for(mc in metric_cols){
    for(tk in unique(df$task)){
      idx <- which(df$task==tk)
      sub <- df[idx, c(mc, covs)]
      ok <- is.finite(sub[[mc]]) & Reduce(`&`, lapply(covs, function(c) is.finite(sub[[c]])))
      if(sum(ok) > 10 && sd(sub[[mc]][ok],na.rm=TRUE)>0){
        fit <- try(lm(as.formula(paste0("`",mc,"` ~ inhaleVolume+inhaleDuration+peakInspFlow")), data=sub[ok,]), silent=TRUE)
        if(!inherits(fit,"try-error")){ r <- rep(NA_real_, length(idx)); r[ok] <- residuals(fit) + mean(sub[[mc]][ok],na.rm=TRUE); df2[[mc]][idx] <- r }
      }
    }
  }
  df2
}
pbAir <- resid_airflow(pb)
sessAir <- sess_agg(pbAir)
rank_air <- build_rank(sessAir, "airflow_adjusted")

ranks <- bind_rows(rank_all, rank_pk, rank_air)
write_csv(ranks, file.path(tdir,"gamma_goodness_v2.csv"))

# ---------- honest headline: rank on outcome-blind control<->Dupi (all breaths), require |d| CI excludes 0 ----------
headline <- rank_all |>
  filter(cd_sig==1) |>
  mutate(abs_cd=abs(sep_cd)) |>
  arrange(desc(abs_cd)) |>
  select(task, metric, control_CV, sep_cd, sep_cd_lo, sep_cd_hi, sep_rn, n_nonresp, recovery_rho, isRidgeFreq)
write_csv(head(headline,40), file.path(tdir,"gamma_headline_controlDupi.csv"))

# peak-gated survivors of the same measures
surv <- rank_pk |> filter(cd_sig==1) |> mutate(abs_cd=abs(sep_cd)) |> arrange(desc(abs_cd)) |>
  select(task, metric, control_CV, sep_cd, sep_cd_lo, sep_cd_hi)
write_csv(head(surv,40), file.path(tdir,"gamma_headline_peakgated.csv"))

# ---------- (3) aperiodic exponent/offset per group + 1/f-controlled separation ----------
ap <- sessAll |> select(task,sessID,cohort2,grp,sessNum,apExp,apOffset)
apx <- ap |> mutate(xpos=ifelse(cohort2=="Control","Control",paste0("S",sessNum)), xpos=factor(xpos,levels=c("S1","S2","S3","Control")))
pal <- c(responder="#1b9e77",`non-responder`="#d95f02",control="#7570b3",unclassified="#999999")
for(mc in c("apExp","apOffset")){
  d<-apx|>select(task,xpos,grp,val=all_of(mc))|>filter(is.finite(val))
  ms<-d|>group_by(task,xpos)|>summarise(m=mean(val,na.rm=TRUE),se=sd(val,na.rm=TRUE)/sqrt(sum(is.finite(val))),.groups="drop")
  p<-ggplot(d,aes(xpos,val))+geom_jitter(aes(color=grp),width=.12,height=0,size=1.3,alpha=.7)+
    geom_errorbar(data=ms,aes(xpos,ymin=m-se,ymax=m+se),inherit.aes=FALSE,width=.2)+
    geom_point(data=ms,aes(xpos,m),inherit.aes=FALSE,shape=95,size=4)+
    facet_wrap(~task,scales="free_y",nrow=1)+scale_color_manual(values=pal,name="")+
    labs(x=NULL,y=mc,title=paste("Aperiodic",mc,"(Voytek 1/f control)"))+theme_bw(base_size=10)+
    theme(legend.position="bottom",axis.text.x=element_text(size=7))
  ggsave(file.path(fdir,paste0("ap_",mc,".png")),p,width=13,height=3.2,dpi=120)
}
# is the control<->Dupi separation aperiodic? separation on apExp itself:
ap_sep <- ap |> group_by(task) |> summarise(
  apExp_sep_cd = d_ci(apExp[cohort2=="Dupi"], apExp[cohort2=="Control"])["d"],
  apOffset_sep_cd = d_ci(apOffset[cohort2=="Dupi"], apOffset[cohort2=="Control"])["d"], .groups="drop")
write_csv(ap_sep, file.path(tdir,"gamma_aperiodic_separation.csv"))

# 1/f-controlled: residualize a few top power measures on apExp per task, re-test control<->Dupi
onef_ctrl <- list()
for(mc in intersect(c("w5_rpowZ","w5_aucZ","gammaBumpDb","peakZ","w2_rpowZ"), metric_cols)){
  for(tk in unique(sessAll$task)){
    st <- sessAll|>filter(task==tk); st<-st[is.finite(st[[mc]])&is.finite(st$apExp),]
    if(nrow(st)<8) next
    st$resid <- residuals(lm(as.formula(paste0("`",mc,"` ~ apExp")), data=st)) + mean(st[[mc]])
    raw <- d_ci(st[[mc]][st$cohort2=="Dupi"], st[[mc]][st$cohort2=="Control"])["d"]
    ctl <- d_ci(st$resid[st$cohort2=="Dupi"], st$resid[st$cohort2=="Control"])["d"]
    onef_ctrl[[length(onef_ctrl)+1]] <- tibble(task=tk,metric=mc,sep_raw=raw,sep_after_1f_control=ctl)
  }
}
write_csv(list_rbind(onef_ctrl), file.path(tdir,"gamma_1f_controlled_separation.csv"))

# ---------- (5) whole-window vs baseline-z concordance ----------
conc <- pb |> group_by(task) |> summarise(
  peakZ_vs_base = suppressWarnings(cor(peakZ, peakZ_base, use="complete.obs")),
  timeAbove_vs_base = suppressWarnings(cor(timeAboveMs, timeAboveMs_base, use="complete.obs")),
  anyBurst_agree = mean(anyBurst==anyBurst_base, na.rm=TRUE), .groups="drop")
write_csv(conc, file.path(tdir,"gamma_zwindow_concordance.csv"))

# ---------- console summary ----------
cat("=== Headline: outcome-blind control<->Dupi separations (|d| CI excludes 0), top 15 ===\n")
print(as.data.frame(head(headline,15) |> mutate(across(where(is.numeric),~round(.x,3)))), row.names=FALSE)
cat("\n=== Does the same measure survive peak-gating? (top 10) ===\n")
print(as.data.frame(head(surv,10) |> mutate(across(where(is.numeric),~round(.x,3)))), row.names=FALSE)
cat("\n=== Is the effect aperiodic? separation on apExp/apOffset itself ===\n")
print(as.data.frame(ap_sep |> mutate(across(where(is.numeric),~round(.x,2)))), row.names=FALSE)
cat("\n=== 1/f-controlled separations (raw vs after residualizing on apExp) ===\n")
print(as.data.frame(list_rbind(onef_ctrl) |> mutate(across(where(is.numeric),~round(.x,2)))), row.names=FALSE)
cat("\n=== whole-window vs baseline-z concordance ===\n")
print(as.data.frame(conc |> mutate(across(where(is.numeric),~round(.x,3)))), row.names=FALSE)
cat("\nwrote gamma_goodness_v2.csv, gamma_headline_*.csv, gamma_aperiodic_separation.csv, gamma_1f_controlled_separation.csv, gamma_zwindow_concordance.csv, ap_*.png\n")
