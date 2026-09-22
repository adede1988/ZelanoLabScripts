# Sensitivity: non-olfactory = focusedBreathing ONLY (drop audiobook), a better-matched
# active/attentional control. Same parametric 2 bands x 4 windows model on the cached grid.
suppressWarnings(suppressMessages({library(lme4); library(lmerTest); library(emmeans)}))
setwd("C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2")
T<-"out/tables"; Fg<-"out/figs/taskcmp"; WL<-c("0-500","500-1000","1000-1500","1500-2000")

d <- read.csv(file.path(T,"taskcmp_olf_roi_grid.csv"), stringsAsFactors=FALSE)
d <- d[d$category=="olf" | d$block=="focusedBreathing", ]     # drop audiobook
d$band<-factor(d$band,levels=c("high","low")); d$window<-factor(d$window,levels=WL)
d$taskCategory<-factor(d$category,levels=c("nonolf","olf")); d$sessID<-factor(d$sessID)
cat(sprintf("blocks: olf=%d, focusedBreathing=%d | sessions=%d\n",
    length(unique(d$block[d$category=="olf"])), sum(d$block=="focusedBreathing")/8, nlevels(d$sessID)))

cells <- function(resp){
  d$y<-d[[resp]]
  ms<-lmer(y~band*window*taskCategory+(1|sessID),data=d,REML=TRUE,
           contrasts=list(band=contr.sum,window=contr.sum,taskCategory=contr.sum))
  cat("\n#### ",resp," Type III ANOVA ####\n"); print(anova(ms,type=3))
  m<-lmer(y~band*window*taskCategory+(1|sessID),data=d,REML=TRUE)
  ct<-as.data.frame(pairs(emmeans(m,~taskCategory|band*window),reverse=TRUE))
  ct<-ct[order(ct$band,match(ct$window,WL)),]; ct$p_holm<-p.adjust(ct$p.value,"holm")
  cat("\n-- olf - focusedBreathing per band x window (holm across 8) --\n")
  print(ct[,c("band","window","estimate","SE","t.ratio","p.value","p_holm")],digits=3)
  ct
}
cat("\n============ PRIMARY medZ (olf vs focus only) ============\n"); P<-cells("medZ")
write.csv(P, file.path(T,"taskcmp_olf_grid_contrasts_focusonly.csv"), row.names=FALSE)
cat("\n============ ROBUSTNESS medZ_ninv (olf vs focus only) ============\n"); Q<-cells("medZ_ninv")
write.csv(Q, file.path(T,"taskcmp_olf_grid_contrasts_focusonly_ninv.csv"), row.names=FALSE)

draw<-function(resp,ctab,fn){
  png(file.path(Fg,fn),width=1500,height=760,res=140)
  par(mfrow=c(2,4),mar=c(3.4,4,2.6,0.6),oma=c(0,0,2.2,0)); GREY<-"#9e9e9e"; GREEN<-"#22883a"
  yr<-range(d[[resp]],na.rm=TRUE)
  for(bd in c("high","low")) for(w in WL){
    s<-d[d$band==bd&d$window==w,]; y<-s[[resp]]; xj<-as.integer(s$taskCategory)+runif(nrow(s),-0.08,0.08)
    plot(xj,y,xlim=c(0.6,2.4),ylim=yr,xaxt="n",pch=19,cex=0.9,
         col=ifelse(s$taskCategory=="olf",GREEN,GREY),
         xlab="",ylab=if(w=="0-500")sprintf("%s band | %s",bd,resp) else "",
         main=sprintf("%s Hz | %s ms",if(bd=="high")"40-58" else "25-40",w),cex.main=0.95)
    axis(1,at=1:2,labels=c("focus","olf"))
    mu<-tapply(y,s$taskCategory,mean); segments(c(0.8,1.8),mu,c(1.2,2.2),mu,lwd=3,col=c(GREY,GREEN)); abline(h=0,lty=3,col="gray60")
    r<-ctab[ctab$band==bd&ctab$window==w,]
    star<-if(r$p_holm<.001)"***" else if(r$p_holm<.01)"**" else if(r$p_holm<.05)"*" else if(r$p.value<.05)"(*)" else ""
    mtext(sprintf("d=%+.2f p=%.3g%s",r$estimate,r$p.value,star),side=3,line=-1.1,cex=0.72)
  }
  mtext(sprintf("Control: olfactory (green) vs FOCUSED BREATHING only (grey), %s",resp),outer=TRUE,cex=0.95,font=2)
  dev.off()
}
draw("medZ",P,"olf_grid_focusonly_medZ.png"); draw("medZ_ninv",Q,"olf_grid_focusonly_medZ_ninv.png")
cat("\nwrote taskcmp_olf_grid_contrasts_focusonly{,_ninv}.csv, olf_grid_focusonly_medZ{,_ninv}.png\n")
