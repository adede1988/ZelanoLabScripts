# Olfactory vs non-olfactory gamma: mixed model on ROI median z (myChanZscore, superlet [3 30]).
# Primary model (user spec):  medZ ~ gammaType * taskCategory + (1|sessionID)
#   gammaType: earlyHigh (40-58 Hz, 0-1000 ms) vs lateLow (25-40 Hz, 1000-2000 ms)
#   taskCategory: olf (cue+thresh+O15) vs nonolf (audiobook+focusedBreathing). Control, CP dropped.
# Response variants:
#   medZ      = myChanZscore ROI median (SEM denom -> scales ~sqrt(nBreaths))  [user spec]
#   medZ_ninv = medZ/sqrt(nBreaths) = per-breath-scale z, trial-count invariant [robustness]
suppressWarnings(suppressMessages({library(lme4); library(lmerTest); library(emmeans)}))
setwd("C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2")
T <- "out/tables"; Fg <- "out/figs/taskcmp"

d <- read.csv(file.path(T,"taskcmp_olf_roi.csv"), stringsAsFactors=FALSE)
d$gammaType    <- factor(d$gammaType, levels=c("earlyHigh","lateLow"))
d$taskCategory <- factor(d$category,  levels=c("nonolf","olf"))
d$sessID <- factor(d$sessID); d$participant <- factor(d$participant)
cat(sprintf("rows=%d | blocks=%d | sessions=%d | participants=%d\n",
            nrow(d), nrow(d)/2, nlevels(d$sessID), nlevels(d$participant)))
print(table(d$taskCategory, d$gammaType))
both <- names(which(tapply(d$taskCategory, d$sessID, function(x) length(unique(x))==2)))

fit_all <- function(resp){
  cat("\n\n#################### RESPONSE:", resp, "####################\n")
  d$y <- d[[resp]]
  ## PRIMARY: (1|sessionID)
  cat("\n== PRIMARY  y ~ gammaType*taskCategory + (1|sessID) ==\n")
  m <- lmer(y ~ gammaType*taskCategory + (1|sessID), data=d, REML=TRUE)
  print(round(summary(m)$coefficients,4))
  ms <- lmer(y ~ gammaType*taskCategory + (1|sessID), data=d, REML=TRUE,
             contrasts=list(gammaType=contr.sum, taskCategory=contr.sum))
  cat("\n-- Type III ANOVA --\n"); print(anova(ms, type=3))
  em <- emmeans(m, ~ taskCategory | gammaType)
  cat("\n-- cell means --\n"); print(em)
  cat("\n-- olf - nonolf within gammaType --\n"); print(pairs(em, reverse=TRUE))
  cat("\n-- interaction (olf-nonolf: earlyHigh vs lateLow) --\n")
  print(contrast(emmeans(m, ~ taskCategory*gammaType), interaction="pairwise"))
  if(resp=="medZ"){
    co <- as.data.frame(summary(m)$coefficients); co$term <- rownames(co)
    write.csv(co, file.path(T,"taskcmp_olf_lmm_coef.csv"), row.names=FALSE)
    write.csv(as.data.frame(anova(ms,type=3)), file.path(T,"taskcmp_olf_lmm_anova_typeIII.csv"))
    write.csv(as.data.frame(em), file.path(T,"taskcmp_olf_lmm_emmeans.csv"), row.names=FALSE)
  }
  ## SENS 1: participant random effect
  cat("\n== SENS 1  (1|participant/sessID) ==\n")
  m2 <- tryCatch(lmer(y ~ gammaType*taskCategory + (1|participant/sessID), data=d, REML=TRUE),
                 error=function(e){cat("nested singular -> (1|participant)\n");
                                   lmer(y ~ gammaType*taskCategory + (1|participant), data=d, REML=TRUE)})
  print(round(summary(m2)$coefficients,4))
  ## SENS 2: within-session (both-category sessions only)
  db <- droplevels(d[d$sessID %in% both,])
  cat(sprintf("\n== SENS 2  both-category sessions only (n=%d: %s) ==\n", length(both), paste(both,collapse=", ")))
  m3 <- lmer(y ~ gammaType*taskCategory + (1|sessID), data=db, REML=TRUE)
  print(round(summary(m3)$coefficients,4))
  cat("-- within-session olf - nonolf by gammaType --\n")
  print(pairs(emmeans(m3, ~ taskCategory | gammaType), reverse=TRUE))
}

fit_all("medZ")
fit_all("medZ_ninv")

## ---- figure: 2x2 interaction for BOTH responses (block points + cell means) ----
png(file.path(Fg,"olf_lmm_interaction.png"), width=1200, height=900, res=140)
par(mfrow=c(2,2), mar=c(4,4.4,2.6,1)); GREY<-"#9e9e9e"; GREEN<-"#22883a"
for(resp in c("medZ","medZ_ninv")) for(gt in levels(d$gammaType)){
  s <- d[d$gammaType==gt,]; y <- s[[resp]]
  xj <- as.integer(s$taskCategory) + runif(nrow(s),-0.08,0.08)
  plot(xj, y, xlim=c(0.6,2.4), xaxt="n", pch=19,
       col=ifelse(s$taskCategory=="olf",GREEN,GREY),
       xlab="", ylab=resp, main=sprintf("%s | %s", resp, gt), cex.main=1.0, cex.lab=1.0)
  axis(1, at=1:2, labels=c("non-olf","olfactory"))
  mu <- tapply(y, s$taskCategory, mean); segments(c(0.8,1.8),mu,c(1.2,2.2),mu,lwd=3,col=c(GREY,GREEN))
  abline(h=0, lty=3, col="gray60")
}
dev.off()
cat("\nwrote taskcmp_olf_lmm_{coef,anova_typeIII,emmeans}.csv, olf_lmm_interaction.png\n")
