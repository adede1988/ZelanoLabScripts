# taskcmp_C2_decoding.R — rigorous multivariate decoding of audiobook vs focus,
# within each control subject, with (1) a proper within-subject label-permutation
# null and (2) a within-session TEMPORAL CONFOUND control. Audiobook and focus are
# (partly) separate blocks, so slow drift could separate them with no condition
# effect. We report: gamma-only OOF-AUC, time-only AUC (breath order as sole
# feature), and gamma-after-detrending AUC (each feature linearly residualised on
# breath order inside the training fold). Speed: the imputation/detrend/PCA are
# label-independent and precomputed once per fold; each permutation only re-fits
# the cheap LDA head.
suppressWarnings(suppressMessages({library(dplyr); library(readr); library(stringr); library(purrr)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
T <- file.path(proj,"out","tables"); Fg <- file.path(proj,"out","figs","taskcmp")
NPERM <- if (length(args)>=2) as.integer(args[[2]]) else 500
set.seed(11)

pb  <- read_csv(file.path(T,"taskcmp_perbreath.csv"), show_col_types=FALSE)
ml  <- read_csv(file.path(T,"taskcmp_metric_list.csv"), show_col_types=FALSE)
gam <- intersect(ml$metric[ml$family!="airflow_morphology"], names(pb))

auc <- function(score,y){ ok<-is.finite(score); score<-score[ok]; y<-y[ok]
  n1<-sum(y==1); n0<-sum(y==0); if(n1==0||n0==0) return(NA_real_)
  r<-rank(score); (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }

prep_folds <- function(X, ord, K=5){
  n<-nrow(X); folds<-sample(rep(1:K,length.out=n)); FP<-vector("list",K)
  for(k in 1:K){ tr<-which(folds!=k); te<-which(folds==k)
    Xtr<-X[tr,,drop=FALSE]; Xte<-X[te,,drop=FALSE]
    med<-apply(Xtr,2,median,na.rm=TRUE); med[!is.finite(med)]<-0
    fix<-function(M){ for(j in 1:ncol(M)){ v<-M[,j]; v[!is.finite(v)]<-med[j]; M[,j]<-v }; M }
    FP[[k]]<-list(tr=tr,te=te,Xtr=fix(Xtr),Xte=fix(Xte),ttr=ord[tr],tte=ord[te]) }
  FP
}
process_fold <- function(fp, detrend){          # label-independent feature processing
  Xtr<-fp$Xtr; Xte<-fp$Xte
  if(detrend){ Dtr<-cbind(1,fp$ttr); B<-solve(crossprod(Dtr)+diag(1e-8,2),crossprod(Dtr,Xtr))
    Xtr<-Xtr-Dtr%*%B; Xte<-Xte-cbind(1,fp$tte)%*%B }
  keep<-apply(Xtr,2,sd)>1e-9; if(sum(keep)<2) return(NULL)
  Xtr<-Xtr[,keep,drop=FALSE]; Xte<-Xte[,keep,drop=FALSE]
  mu<-colMeans(Xtr); sdv<-apply(Xtr,2,sd); Xtr<-scale(Xtr,mu,sdv); Xte<-scale(Xte,mu,sdv)
  pc<-prcomp(Xtr,center=FALSE,scale.=FALSE); npc<-min(10,sum(pc$sdev>1e-8),nrow(Xtr)-2); if(npc<1) return(NULL)
  list(tr=fp$tr, te=fp$te, Ztr=pc$x[,1:npc,drop=FALSE], Zte=Xte%*%pc$rotation[,1:npc,drop=FALSE])
}
score_oof <- function(PF, y){                    # cheap LDA head per fold
  s<-rep(NA_real_,length(y))
  for(P in PF){ if(is.null(P)) next; ytr<-y[P$tr]; if(length(unique(ytr))<2) next
    m1<-colMeans(P$Ztr[ytr==1,,drop=FALSE]); m0<-colMeans(P$Ztr[ytr==0,,drop=FALSE])
    S<-cov(P$Ztr)+diag(1e-3,ncol(P$Ztr)); w<-tryCatch(solve(S,m1-m0),error=function(e)rep(0,length(m1)))
    s[P$te]<-as.numeric(P$Zte%*%w) }
  s
}

subs<-sort(unique(pb$subject)); rows<-list()
for(su in subs){
  d<-pb|>filter(subject==su)|>arrange(sessNum,breathIdx)
  y<-as.integer(d$condition=="focus"); ord<-seq_len(nrow(d)); X<-as.matrix(d[,gam])
  FP<-prep_folds(X,ord); PF0<-lapply(FP,process_fold,detrend=FALSE); PFd<-lapply(FP,process_fold,detrend=TRUE)
  real <-auc(score_oof(PF0,y),y); resid<-auc(score_oof(PFd,y),y)
  timeauc<-auc(ord,y); timeauc<-max(timeauc,1-timeauc)
  nr<-numeric(NPERM); nd<-numeric(NPERM)
  for(b in 1:NPERM){ yp<-sample(y); nr[b]<-auc(score_oof(PF0,yp),yp); nd[b]<-auc(score_oof(PFd,yp),yp) }
  rows[[su]]<-tibble(subject=su,nBreath=nrow(d),nFocus=sum(y),auc_gamma=real,
    auc_gamma_detrended=resid,auc_time_only=timeauc,
    perm_p_gamma=(1+sum(nr>=real,na.rm=TRUE))/(NPERM+1),
    perm_p_detrended=(1+sum(nd>=resid,na.rm=TRUE))/(NPERM+1))
}
R<-bind_rows(rows); write_csv(R, file.path(T,"taskcmp_C2_decoding.csv"))
fisher<-function(p){ p<-pmin(pmax(p,1e-6),1); pchisq(-2*sum(log(p)),df=2*length(p),lower.tail=FALSE) }
grp<-tibble(measure=c("gamma","gamma_detrended","time_only"),
  mean_auc=c(mean(R$auc_gamma),mean(R$auc_gamma_detrended),mean(R$auc_time_only)),
  n_gt_chance=c(sum(R$auc_gamma>0.5),sum(R$auc_gamma_detrended>0.5),sum(R$auc_time_only>0.5)),
  wilcox_p=c(tryCatch(wilcox.test(R$auc_gamma-0.5)$p.value,error=function(e)NA),
             tryCatch(wilcox.test(R$auc_gamma_detrended-0.5)$p.value,error=function(e)NA),
             tryCatch(wilcox.test(R$auc_time_only-0.5)$p.value,error=function(e)NA)),
  fisher_perm_p=c(fisher(R$perm_p_gamma),fisher(R$perm_p_detrended),NA))
write_csv(grp, file.path(T,"taskcmp_C2_group.csv"))
cat(sprintf("NPERM=%d\n",NPERM)); print(as.data.frame(R),row.names=FALSE); cat("\nGROUP:\n"); print(as.data.frame(grp),row.names=FALSE)
png(file.path(Fg,"C2_decoding.png"),width=1150,height=560,res=120); op<-par(mar=c(4,5,3,1))
M<-t(as.matrix(R[,c("auc_gamma","auc_gamma_detrended","auc_time_only")])); colnames(M)<-R$subject
barplot(M,beside=TRUE,ylim=c(0,1),col=c("#1b7837","#a6dba0","#999999"),
  legend.text=c("gamma","gamma detrended","time only"),args.legend=list(x="topright",bty="n"),
  las=1,ylab="OOF AUC (focus vs audiobook)",main=sprintf("Lens C rigorous: per-subject decoding (NPERM=%d)",NPERM))
abline(h=0.5,lty=2); par(op); dev.off()
cat("\nwrote taskcmp_C2_decoding.csv, taskcmp_C2_group.csv, C2_decoding.png\n")
