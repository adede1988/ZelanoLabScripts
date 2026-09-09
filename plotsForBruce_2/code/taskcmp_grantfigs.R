# taskcmp_grantfigs.R — final "grant figs": the two key GAMMA measures
# (peak latency ms; early 100-500 ms band-max), baseline vs ATB, with connected
# per-subject dots, baseline=translucent grey, ATB=forest green. Tall+narrow, big
# bold axes. Renames audiobook->baseline, focus->ATB. Plus two power analyses:
#   (1) raw CONTROL paired baseline-vs-ATB difference;
#   (2) group x condition INTERACTION, Control vs Dupi S1 (two-sample on the
#       within-subject ATB-baseline contrast), for each measure.
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(ggplot2)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
T <- file.path(proj,"out","tables"); Fg <- file.path(proj,"out","figs","taskcmp")
GREY <- "#9e9e9e"; GREEN <- "#228B22"
bm <- read_csv(file.path(T,"taskcmp_bandmetrics.csv"), show_col_types=FALSE) |>
  filter(band=="gamma", group %in% c("Control","DupiS1","DupiS2")) |>
  mutate(cond = recode(condition, audiobook="baseline", focus="ATB"),
         cond = factor(cond, levels=c("baseline","ATB")),
         group = factor(group, levels=c("Control","DupiS1","DupiS2")))
# SESSION-level n (each session is one observation; AD's two control sessions count separately)
subj <- bm |> group_by(group, sessID, cond) |>
  summarise(peakLatMs=mean(peakLatMs,na.rm=TRUE), earlyBandMax=mean(earlyBandMax,na.rm=TRUE), .groups="drop")

metspec <- list(
  list(col="peakLatMs",    ylab="gamma peak latency (ms)",           file="grantfig_gamma_peakLatMs.png"),
  list(col="earlyBandMax", ylab="gamma early band-max (100-500 ms, z)", file="grantfig_gamma_earlyBandMax.png"))

grantplot <- function(mc, ylab, outfile){
  d <- subj |> select(group, sessID, cond, val=all_of(mc)) |> filter(is.finite(val))
  gx <- as.numeric(d$group); d$xpos <- gx + ifelse(d$cond=="ATB", 0.19, -0.19)
  p <- ggplot(d, aes(xpos, val)) +
    geom_line(aes(group=interaction(sessID,group)), color=GREY, alpha=0.55, linewidth=0.5) +
    geom_boxplot(aes(group=interaction(group,cond), fill=cond), width=0.32,
                 outlier.shape=NA, alpha=0.55, color="black", linewidth=0.7) +
    geom_point(aes(fill=cond), shape=21, size=2.4, stroke=0.5, color="black",
               position=position_jitter(width=0.03, height=0)) +
    scale_fill_manual(values=c(baseline=GREY, ATB=GREEN), name=NULL) +
    scale_x_continuous(breaks=1:3, labels=levels(d$group), expand=expansion(add=0.45)) +
    labs(x=NULL, y=ylab) +
    theme_classic(base_size=20) +
    theme(legend.position="top", legend.text=element_text(size=18),
          axis.text=element_text(size=17, color="black"),
          axis.text.x=element_text(size=18, face="bold"),
          axis.title.y=element_text(size=19, face="bold"),
          axis.line=element_line(linewidth=1.1, color="black"),
          axis.ticks=element_line(linewidth=1.0, color="black"), axis.ticks.length=unit(4,"pt"))
  ggsave(file.path(Fg,outfile), p, width=5.2, height=7.6, dpi=200)
  invisible(p)
}
for(m in metspec) grantplot(m$col, m$ylab, m$file)
cat("wrote grantfig_gamma_peakLatMs.png, grantfig_gamma_earlyBandMax.png\n")

# ---------------- power analyses ----------------
n80  <- function(d, type){ if(!is.finite(d)||d==0) return(NA); ceiling(tryCatch(power.t.test(delta=abs(d),sd=1,sig.level=.05,power=.8,type=type)$n,error=function(e)NA)) }
pw   <- function(d, n, type){ if(!is.finite(d)||d==0||n<2) return(NA); tryCatch(power.t.test(n=n,delta=abs(d),sd=1,sig.level=.05,type=type)$power,error=function(e)NA) }
dz_paired <- function(a,f){ ok<-is.finite(a)&is.finite(f); d<-a[ok]-f[ok]; mean(d)/sd(d) }
POW <- list()
for(mc in c("peakLatMs","earlyBandMax")){
  w <- subj |> select(group,sessID,cond,val=all_of(mc)) |>
       pivot_wider(names_from=cond, values_from=val)
  ctrl <- w |> filter(group=="Control"); s1 <- w |> filter(group=="DupiS1")
  # (1) raw control paired diff
  dzc <- dz_paired(ctrl$baseline, ctrl$ATB); nc <- sum(is.finite(ctrl$baseline)&is.finite(ctrl$ATB))
  POW[[length(POW)+1]] <- tibble(measure=mc, test="control baseline-vs-ATB (paired)",
    effect=round(dzc,2), effect_type="dz", n_now=as.character(nc), power_now=round(pw(dzc,nc,"paired"),2),
    n_for_80=n80(dzc,"paired"))
  # (2) interaction: ATB-baseline contrast, Control vs DupiS1 (two-sample)
  cC <- ctrl$ATB-ctrl$baseline; cS <- s1$ATB-s1$baseline; cC<-cC[is.finite(cC)]; cS<-cS[is.finite(cS)]
  sp <- sqrt(((length(cC)-1)*var(cC)+(length(cS)-1)*var(cS))/(length(cC)+length(cS)-2))
  dint <- (mean(cC)-mean(cS))/sp
  pint <- tryCatch(t.test(cC,cS)$p.value, error=function(e)NA)
  POW[[length(POW)+1]] <- tibble(measure=mc, test="interaction Control x DupiS1 (two-sample on ATB-baseline)",
    effect=round(dint,2), effect_type="d", n_now=sprintf("%d vs %d",length(cC),length(cS)),
    power_now=round(pw(dint, 2*length(cC)*length(cS)/(length(cC)+length(cS)), "two.sample"),2),
    n_for_80=n80(dint,"two.sample"))
}
POWt <- bind_rows(POW); write_csv(POWt, file.path(T,"taskcmp_grant_power.csv"))
cat("\n=== Grant power analyses ===\n"); print(as.data.frame(POWt), row.names=FALSE)
cat("\n(n_for_80: paired = n subjects; two-sample = n PER GROUP for 80% power, alpha=.05)\n")
