# Parametric olfactory vs non-olfactory gamma: 2 bands x 4 windows.
# medZ ~ band*window*taskCategory + (1|sessionID)  (myChanZscore ROI medians, superlet [3 30]).
# Per-cell olf-nonolf contrast (8 cells) + holm correction; n-invariant robustness; 2x4 figure.
suppressWarnings(suppressMessages({library(lme4); library(lmerTest); library(emmeans)}))
setwd("C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2")
T<-"out/tables"; Fg<-"out/figs/taskcmp"
WL <- c("0-500","500-1000","1000-1500","1500-2000")

d <- read.csv(file.path(T,"taskcmp_olf_roi_grid.csv"), stringsAsFactors=FALSE)
d$band         <- factor(d$band, levels=c("high","low"))
d$window       <- factor(d$window, levels=WL)
d$taskCategory <- factor(d$category, levels=c("nonolf","olf"))
d$sessID<-factor(d$sessID); d$participant<-factor(d$participant)
d$winMid <- as.numeric(d$winMid)
cat(sprintf("rows=%d | blocks=%d | sessions=%d | participants=%d\n",
            nrow(d), nrow(d)/8, nlevels(d$sessID), nlevels(d$participant)))

cell_contrasts <- function(resp){
  d$y <- d[[resp]]
  m <- lmer(y ~ band*window*taskCategory + (1|sessID), data=d, REML=TRUE)
  cat("\n#### RESPONSE:",resp," -- Type III ANOVA ####\n")
  ms <- lmer(y ~ band*window*taskCategory + (1|sessID), data=d, REML=TRUE,
             contrasts=list(band=contr.sum, window=contr.sum, taskCategory=contr.sum))
  print(anova(ms, type=3))
  em <- emmeans(m, ~ taskCategory | band*window)
  ct <- as.data.frame(pairs(em, reverse=TRUE))          # olf - nonolf per cell
  ct <- ct[order(ct$band, match(ct$window, WL)),]
  ct$p_holm <- p.adjust(ct$p.value, "holm")
  cat("\n-- olf - nonolf per band x window (holm across 8 cells) --\n")
  print(ct[,c("band","window","estimate","SE","t.ratio","p.value","p_holm")], digits=3)
  # parametric trend across windows within band (does olf-nonolf change over time?)
  mt <- lmer(y ~ band*scale(winMid)*taskCategory + (1|sessID), data=d, REML=TRUE)
  cat("\n-- olf-nonolf linear trend across windows, by band (emtrends on winMid) --\n")
  print(emtrends(mt, ~ band, var="winMid", at=list(taskCategory=c("nonolf","olf"))) |> summary(), digits=4)
  print(pairs(emtrends(mt, ~ taskCategory|band, var="winMid"), reverse=TRUE))
  list(model=m, em=em, ct=ct)
}

cat("\n================= PRIMARY: medZ (myChanZscore) =================\n")
P <- cell_contrasts("medZ")
write.csv(P$ct, file.path(T,"taskcmp_olf_grid_contrasts.csv"), row.names=FALSE)

cat("\n================= ROBUSTNESS: medZ_ninv (trial-count invariant) =================\n")
Q <- cell_contrasts("medZ_ninv")
write.csv(Q$ct, file.path(T,"taskcmp_olf_grid_contrasts_ninv.csv"), row.names=FALSE)

## ---- 2 bands x 4 windows figure (rows=band, cols=window) ----
draw <- function(resp, ctab, fn){
  png(file.path(Fg,fn), width=1500, height=760, res=140)
  par(mfrow=c(2,4), mar=c(3.4,4,2.6,0.6), oma=c(0,0,2.2,0))
  GREY<-"#9e9e9e"; GREEN<-"#22883a"
  yr <- range(d[[resp]], na.rm=TRUE)
  for(bd in c("high","low")) for(w in WL){
    s <- d[d$band==bd & d$window==w,]; y<-s[[resp]]
    xj <- as.integer(s$taskCategory)+runif(nrow(s),-0.08,0.08)
    plot(xj,y,xlim=c(0.6,2.4),ylim=yr,xaxt="n",pch=19,cex=0.9,
         col=ifelse(s$taskCategory=="olf",GREEN,GREY),
         xlab="",ylab=if(w=="0-500") sprintf("%s band  |  %s",bd,resp) else "",
         main=sprintf("%s Hz | %s ms", if(bd=="high")"40-58" else "25-40", w), cex.main=0.95)
    axis(1,at=1:2,labels=c("non","olf"))
    mu<-tapply(y,s$taskCategory,mean); segments(c(0.8,1.8),mu,c(1.2,2.2),mu,lwd=3,col=c(GREY,GREEN))
    abline(h=0,lty=3,col="gray60")
    r <- ctab[ctab$band==bd & ctab$window==w,]
    star <- if(r$p_holm<.001)"***" else if(r$p_holm<.01)"**" else if(r$p_holm<.05)"*" else if(r$p.value<.05)"(*)" else ""
    mtext(sprintf("d=%+.2f p=%.3g%s", r$estimate, r$p.value, star), side=3, line=-1.1, cex=0.72)
  }
  mtext(sprintf("Control: olfactory (green) vs non-olfactory (grey), %s  —  * holm-sig, (*) uncorrected only", resp),
        outer=TRUE, cex=0.95, font=2)
  dev.off()
}
draw("medZ",      P$ct, "olf_grid_medZ.png")
draw("medZ_ninv", Q$ct, "olf_grid_medZ_ninv.png")
cat("\nwrote taskcmp_olf_grid_contrasts{,_ninv}.csv, olf_grid_medZ{,_ninv}.png\n")
