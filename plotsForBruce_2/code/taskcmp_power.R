# taskcmp_power.R — sample-size / power analysis for the paired timing effects.
# Frames it prospectively (n needed to detect an effect of a given standardized
# size dz at 80/90% power, paired t-test), plus the current power at n=5 and n=6.
# Post-hoc power on the observed dz is only a planning heuristic — the observed dz
# is itself uncertain at n=5 — so we also tabulate a range of plausible effect sizes.
suppressWarnings(suppressMessages({library(readr); library(dplyr)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
T <- file.path(proj,"out","tables"); Fg <- file.path(proj,"out","figs","taskcmp")

n_for_power <- function(d, power=0.8, alpha=0.05){
  if(!is.finite(d)||d==0) return(NA)
  out <- tryCatch(power.t.test(delta=abs(d), sd=1, sig.level=alpha, power=power, type="paired")$n,
                  error=function(e) NA); ceiling(out) }
power_at <- function(d, n, alpha=0.05){
  if(!is.finite(d)||d==0||n<2) return(NA)
  tryCatch(power.t.test(n=n, delta=abs(d), sd=1, sig.level=alpha, type="paired")$power, error=function(e) NA) }

# observed timing effects (n=5, CP excluded)
st <- read_csv(file.path(T,"taskcmp_latency_stats.csv"), show_col_types=FALSE)
key <- st |> filter(metric %in% c("peakLatMs","centroidMs","xcorrLagMs")) |>
  transmute(effect=paste(band, metric), dz=abs(dz)) |> arrange(desc(dz))
key <- key |> mutate(power_n5=round(sapply(dz, power_at, 5),2), power_n6=round(sapply(dz, power_at, 6),2),
                     n_for_80=sapply(dz, n_for_power, power=0.8),
                     n_for_90=sapply(dz, n_for_power, power=0.9))
key$dz <- round(key$dz,2)
write_csv(key, file.path(T,"taskcmp_power.csv"))
cat("=== Power for observed paired timing effects (n=5, CP excluded) ===\n")
print(as.data.frame(key), row.names=FALSE)

# generic planning table across a plausible dz range
grid <- tibble(dz=c(0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2)) |>
  mutate(n_for_80=sapply(dz, n_for_power, 0.8), n_for_90=sapply(dz, n_for_power, 0.9),
         power_at_n5=round(sapply(dz, power_at, 5),2), power_at_n12=round(sapply(dz, power_at, 12),2))
write_csv(grid, file.path(T,"taskcmp_power_grid.csv"))
cat("\n=== Planning: paired t-test, alpha=0.05 ===\n"); print(as.data.frame(grid), row.names=FALSE)

# power curves for the headline effect sizes
png(file.path(Fg,"power_curves.png"), width=1050, height=560, res=120)
ns <- 3:24; cols <- c("#1b7837","#4575b4","#d73027")
ds <- c(0.90,0.72,0.50); labs <- c("gamma peak latency (dz=0.90)","theta centroid (dz=0.72)","theta peak latency (dz=0.50)")
plot(NA, xlim=range(ns), ylim=c(0,1), xlab="n (paired subjects)", ylab="power (paired t-test, alpha=0.05)",
     main="Power to detect the observed timing effects", las=1)
abline(h=0.8, lty=3, col="grey40"); abline(v=5, lty=2, col="grey60")
for(i in seq_along(ds)){ y <- sapply(ns, function(n) power_at(ds[i], n)); lines(ns, y, col=cols[i], lwd=2.5) }
legend("bottomright", legend=labs, col=cols, lwd=2.5, bty="n")
text(5, 0.02, "n=5 now", col="grey40", pos=4, cex=0.9); dev.off()
cat("\nwrote taskcmp_power.csv, taskcmp_power_grid.csv, power_curves.png\n")
