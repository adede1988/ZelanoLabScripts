# coupling_aggregate.R — session respiration-gamma coupling: table, plots, goodness.
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2); library(purrr)}))
args <- commandArgs(trailingOnly=TRUE)
proj <- if (length(args)>=1) args[[1]] else "C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2"
tdir <- file.path(proj,"out","tables"); cdir <- file.path(proj,"out","gamma","coupling"); fdir <- file.path(proj,"out","figs","gamma")
dir.create(fdir, showWarnings=FALSE, recursive=TRUE)
files <- list.files(cdir, pattern="\\.csv$", full.names=TRUE)
if(length(files)==0){ cat("no coupling files\n"); quit(save="no") }
cp <- files |> map(~suppressMessages(read_csv(.x, show_col_types=FALSE))) |> list_rbind()
resp <- read_csv(file.path(tdir,"responder_table.csv"), show_col_types=FALSE) |> transmute(participant, class)
cp <- cp |> left_join(resp, by="participant") |>
  mutate(cohort=ifelse(cohort=="OBE","Control",cohort),
         grpColor=case_when(cohort=="Control"~"control", class=="responder"~"responder",
                            class=="non-responder"~"non-responder", TRUE~"unclassified"),
         xpos=ifelse(cohort=="Control","Control", paste0("S",sessNum)),
         xpos=factor(xpos, levels=c("S1","S2","S3","Control")))
write_csv(cp, file.path(tdir,"coupling_session.csv"))

pal <- c(responder="#1b9e77", `non-responder`="#d95f02", control="#7570b3", unclassified="#999999")
mets <- c("coup_MI","coup_resultantLen","coup_inhExhRatio")
for(mc in mets){
  d <- cp |> select(taskRow, xpos, grpColor, val=all_of(mc)) |> filter(is.finite(val))
  ms <- d |> group_by(taskRow,xpos) |> summarise(m=mean(val,na.rm=TRUE), se=sd(val,na.rm=TRUE)/sqrt(sum(is.finite(val))), .groups="drop")
  p <- ggplot(d, aes(xpos,val)) + geom_jitter(aes(color=grpColor), width=.12, height=0, size=1.3, alpha=.7) +
    geom_errorbar(data=ms, aes(xpos, ymin=m-se, ymax=m+se), inherit.aes=FALSE, width=.2) +
    geom_point(data=ms, aes(xpos,m), inherit.aes=FALSE, shape=95, size=4) +
    facet_wrap(~taskRow, scales="free_y", nrow=1) + scale_color_manual(values=pal, name="") +
    labs(x=NULL, y=mc, title=paste("Respiration-gamma coupling:", mc)) + theme_bw(base_size=10) +
    theme(legend.position="bottom", axis.text.x=element_text(size=7))
  ggsave(file.path(fdir, paste0("coup_", mc, ".png")), p, width=13, height=3.2, dpi=120)
}
# goodness for coupling measures
d_between <- function(a,b){ a<-a[is.finite(a)]; b<-b[is.finite(b)]; if(length(a)<2||length(b)<2) return(NA)
  sp<-sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2)); if(!is.finite(sp)||sp==0) return(NA); (mean(a)-mean(b))/sp }
G<-list()
for(tk in unique(cp$taskRow)) for(mc in mets){
  st<-cp|>filter(taskRow==tk); ct<-st|>filter(cohort=="Control"); du<-st|>filter(cohort=="Dupi")
  rr<-du|>filter(grpColor=="responder"); nn<-du|>filter(grpColor=="non-responder")
  cv<-{v<-ct[[mc]];v<-v[is.finite(v)]; if(length(v)<2) NA else sd(v)/(abs(mean(v))+1e-9)}
  G[[length(G)+1]]<-tibble(taskRow=tk, metric=mc, control_CV=cv,
    sep_control_dupi=d_between(du[[mc]],ct[[mc]]), sep_resp_nonresp=d_between(rr[[mc]],nn[[mc]]))
}
write_csv(list_rbind(G), file.path(tdir,"coupling_goodness.csv"))
cat("wrote coupling_session.csv, coupling_goodness.csv, coup_*.png\n")
cat("\nper-xpos MI means by taskRow:\n")
print(as.data.frame(cp |> group_by(taskRow,xpos) |> summarise(MI=round(mean(coup_MI,na.rm=TRUE),4), .groups="drop") |>
  pivot_wider(names_from=xpos, values_from=MI)), row.names=FALSE)
