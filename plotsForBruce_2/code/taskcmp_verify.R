## Adversarial verification of taskcmp lenses A-E
suppressMessages({library(readr);library(dplyr);library(tidyr);library(stringr)})
tab <- "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/tables/"
pb <- suppressMessages(read_csv(paste0(tab,"taskcmp_perbreath.csv")))
sm <- suppressMessages(read_csv(paste0(tab,"taskcmp_subject_means.csv")))
ml <- suppressMessages(read_csv(paste0(tab,"taskcmp_metric_list.csv")))

cat("=== 0. sanity: subjects, sessions, condition counts ===\n")
print(pb %>% group_by(subject,sessID,condition) %>% summarise(n=n(),.groups="drop") %>%
      arrange(subject,condition), n=50)

cat("\n=== 1. TEMPORAL BLOCKING within session (Lens C confound) ===\n")
## Is condition confounded with breathIdx (time) within a session?
## If breathIdx ranges of audiobook vs focus are disjoint -> blocked -> random CV inflates AUC.
blk <- pb %>% group_by(subject,sessID) %>%
  summarise(
    ab_min=min(breathIdx[condition=="audiobook"]), ab_max=max(breathIdx[condition=="audiobook"]),
    fo_min=min(breathIdx[condition=="focus"]),     fo_max=max(breathIdx[condition=="focus"]),
    # overlap of the two index ranges
    overlap = max(0, min(ab_max,fo_max) - max(ab_min,fo_min)),
    span    = max(ab_max,fo_max)-min(ab_min,fo_min),
    # rank-biserial: how separable is condition by breathIdx (0.5=interleaved,1=fully blocked)
    auc_time = {
      y <- as.integer(condition=="focus"); x <- breathIdx
      r <- rank(x); n1<-sum(y==1); n0<-sum(y==0)
      (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0)
    },
    .groups="drop") %>%
  mutate(time_sep = pmax(auc_time,1-auc_time))
print(as.data.frame(blk), row.names=FALSE)
cat(sprintf("\nMean |time-separability AUC| = %.3f (0.5=interleaved, 1.0=fully time-blocked)\n",
            mean(blk$time_sep)))

cat("\n=== 2. LOSO robustness for top paired hits (subject means, audiobook-focus) ===\n")
loso <- function(metric){
  col <- paste0(metric,"__mean")
  w <- sm %>% select(subject,condition,val=all_of(col)) %>%
    pivot_wider(names_from=condition,values_from=val) %>% arrange(subject)
  d <- w$audiobook - w$focus
  subs <- w$subject
  full_p <- t.test(w$audiobook,w$focus,paired=TRUE)$p.value
  full_dz <- mean(d)/sd(d)
  cat(sprintf("\n%s : full dz=%.2f  t_p=%.4f  wilcox_p=%.4f  (n=%d)\n",
              metric, full_dz, full_p,
              suppressWarnings(wilcox.test(w$audiobook,w$focus,paired=TRUE)$p.value), length(d)))
  cat("  per-subject diff (audiobook-focus): ",
      paste(sprintf("%s=%.2f",subs,d),collapse="  "),"\n")
  for(i in seq_along(subs)){
    di <- d[-i]
    p <- t.test(w$audiobook[-i],w$focus[-i],paired=TRUE)$p.value
    cat(sprintf("  drop %-3s: dz=%+.2f  t_p=%.4f%s\n", subs[i], mean(di)/sd(di), p,
                ifelse(p<0.05,"","  <-- NS")))
  }
}
for(m in c("p3_aucZ","w5_rpowDb","freqJitter")) loso(m)

cat("\n=== 3. Lens A FDR recomputation across all 93 (independent check) ===\n")
metrics <- ml$metric
tp <- sapply(metrics,function(m){
  col<-paste0(m,"__mean")
  w<-sm%>%select(subject,condition,val=all_of(col))%>%
    pivot_wider(names_from=condition,values_from=val)
  ab<-w$audiobook;fo<-w$focus;ok<-is.finite(ab)&is.finite(fo)
  if(sum(ok)<2||sd(ab[ok]-fo[ok])==0) return(NA)
  t.test(ab[ok],fo[ok],paired=TRUE)$p.value
})
fdr <- p.adjust(tp,method="BH")
cat(sprintf("n metrics with t_p: %d ; min raw p=%.4f ; min FDR=%.4f ; n raw<.05=%d ; n FDR<.10=%d\n",
            sum(!is.na(tp)), min(tp,na.rm=TRUE), min(fdr,na.rm=TRUE),
            sum(tp<0.05,na.rm=TRUE), sum(fdr<0.10,na.rm=TRUE)))
cat("duplicate metric pairs (identical dz => perfectly correlated, inflate the 93 count):\n")
## detect exact-duplicate columns among subject means
mm <- sapply(metrics,function(m) sm[[paste0(m,"__mean")]])
dupes <- c()
for(i in 1:(length(metrics)-1)) for(j in (i+1):length(metrics)){
  a<-mm[,i];b<-mm[,j]
  if(all(is.finite(a))&&all(is.finite(b))&&isTRUE(all.equal(a,b))) dupes<-c(dupes,paste(metrics[i],"==",metrics[j]))
}
cat(paste(" ",dupes,collapse="\n"),"\n")

cat("\n=== 4. Lens C group Wilcoxon vs 0.5, and drop-AS sensitivity ===\n")
C <- suppressMessages(read_csv(paste0(tab,"taskcmp_C_decoding.csv")))
for(fsname in c("gamma_only","gamma_airflow")){
  sub<-C%>%filter(feature_set==fsname)
  p_all<-suppressWarnings(wilcox.test(sub$auc,mu=0.5,alternative="greater")$p.value)
  subx<-sub%>%filter(subject!="AS")
  p_noAS<-suppressWarnings(wilcox.test(subx$auc,mu=0.5,alternative="greater")$p.value)
  cat(sprintf("%-13s AUCs: %s | n>0.5=%d/6 | wilcox p=%.4f | dropAS p=%.4f\n",
              fsname, paste(sprintf("%.2f",sub$auc),collapse=","),
              sum(sub$auc>0.5), p_all, p_noAS))
}
cat("Note: perm_p==0.25 for most => (1+0)/(1+NPERM) with NPERM=3 (preliminary). Fisher perm p (group)=0.228 NS.\n")

cat("\n=== 5. Lens D coupling: confirm identical values across conditions ===\n")
cp <- suppressMessages(read_csv(paste0(tab,"coupling_session.csv")))
cat("columns:", paste(names(cp),collapse=", "),"\n")
print(cp %>% arrange(participant, sessNum) %>% head(20), width=200)
