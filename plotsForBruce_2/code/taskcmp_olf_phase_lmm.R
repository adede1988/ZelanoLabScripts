# Respiratory-phase-matched olfactory vs non-olfactory gamma.
# medZ ~ band*phase*taskCategory + (1|sessionID)  (single-trial myChanZscore, [3 30]).
# Per-cell olf-nonolf contrast (2 bands x 5 phases = 10 cells) + holm; 2x5 figure.
suppressWarnings(suppressMessages({library(lme4); library(lmerTest); library(emmeans)}))
setwd("C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2")
T<-"out/tables"; Fg<-"out/figs/taskcmp"
PH<-c("inhaleRise","peakToExh","exhaleFall","troughToPause","pause")
PHlab<-c("onset->peak","peak->exhOn","exhOn->trough","trough->pause","pause->end")

d<-read.csv(file.path(T,"taskcmp_olf_phase_roi.csv"),stringsAsFactors=FALSE)
d<-d[is.finite(d$medZ),]
d$band<-factor(d$band,levels=c("high","low")); d$phase<-factor(d$phase,levels=PH)
d$taskCategory<-factor(d$category,levels=c("nonolf","olf")); d$sessID<-factor(d$sessID)
cat(sprintf("rows=%d | sessions=%d | olf blocks=%d, nonolf blocks=%d\n",nrow(d),nlevels(d$sessID),
    length(unique(d$block[d$category=="olf"&d$band=="high"&d$phase=="inhaleRise"])),
    length(unique(paste(d$sessID,d$block)[d$category=="nonolf"&d$band=="high"&d$phase=="inhaleRise"]))))
print(table(d$taskCategory,d$phase,d$band))

m<-lmer(medZ~band*phase*taskCategory+(1|sessID),data=d,REML=TRUE)
ms<-lmer(medZ~band*phase*taskCategory+(1|sessID),data=d,REML=TRUE,
         contrasts=list(band=contr.sum,phase=contr.sum,taskCategory=contr.sum))
cat("\n== Type III ANOVA ==\n"); print(anova(ms,type=3))
ct<-as.data.frame(pairs(emmeans(m,~taskCategory|band*phase),reverse=TRUE))  # olf - nonolf per cell
ct<-ct[order(ct$band,match(ct$phase,PH)),]; ct$p_holm<-p.adjust(ct$p.value,"holm")
cat("\n== olf - nonolf per band x phase (holm across 10 cells) ==\n")
print(ct[,c("band","phase","estimate","SE","t.ratio","p.value","p_holm")],digits=3)
write.csv(ct,file.path(T,"taskcmp_olf_phase_contrasts.csv"),row.names=FALSE)

## 2 bands x 5 phases figure
png(file.path(Fg,"olf_phase_grid.png"),width=1650,height=760,res=140)
par(mfrow=c(2,5),mar=c(3.6,4,2.6,0.5),oma=c(0,0,2.2,0)); GREY<-"#9e9e9e"; GREEN<-"#22883a"
yr<-range(d$medZ,na.rm=TRUE)
for(bd in c("high","low")) for(pi in seq_along(PH)){
  s<-d[d$band==bd&d$phase==PH[pi],]; y<-s$medZ; xj<-as.integer(s$taskCategory)+runif(nrow(s),-0.08,0.08)
  plot(xj,y,xlim=c(0.6,2.4),ylim=yr,xaxt="n",pch=19,cex=0.9,
       col=ifelse(s$taskCategory=="olf",GREEN,GREY),
       xlab="",ylab=if(pi==1)sprintf("%s band | medZ",bd) else "",
       main=sprintf("%s Hz\n%s",if(bd=="high")"40-58" else "25-40",PHlab[pi]),cex.main=0.9)
  axis(1,at=1:2,labels=c("non","olf"))
  mu<-tapply(y,s$taskCategory,mean); segments(c(0.8,1.8),mu,c(1.2,2.2),mu,lwd=3,col=c(GREY,GREEN)); abline(h=0,lty=3,col="gray60")
  r<-ct[ct$band==bd&ct$phase==PH[pi],]
  star<-if(length(r$p_holm)&&r$p_holm<.001)"***" else if(length(r$p_holm)&&r$p_holm<.01)"**" else if(length(r$p_holm)&&r$p_holm<.05)"*" else if(length(r$p.value)&&r$p.value<.05)"(*)" else ""
  if(nrow(r)) mtext(sprintf("d=%+.2f p=%.2g%s",r$estimate,r$p.value,star),side=3,line=-1.0,cex=0.68)
}
mtext("Control: olfactory (green) vs non-olfactory=audiobook+focus (grey), phase-matched — medZ (single-trial myChanZscore)",outer=TRUE,cex=0.9,font=2)
dev.off()
cat("\nwrote taskcmp_olf_phase_contrasts.csv, olf_phase_grid.png\n")
