a <- read.csv("C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2/out/tables/taskcmp_B_distribution.csv",
              stringsAsFactors = FALSE)
sel <- c("burstLatMs","p3_aucZ","timeAboveMs_base","p4_aucZ","p2_aucZ","p4_rpowZ",
         "w1_bandDb","w1_rpowDb","w3_rpowDb","p3_bandDb","p3_rpowZ","gammaBumpDb")
cols <- c("metric","family","n_subj","n_pos","n_neg","mean_logVarRatio",
          "varRatio_p","varRatio_p_fdr","varRatio_t_p","varRatio_sign_p",
          "mean_KS","mean_KS_centered","KS_combined_p_fdr",
          "KS_centered_combined_p_fdr","direction")
b <- a[match(sel, a$metric), cols]
print(format(b, digits = 3), row.names = FALSE)

g <- a[!as.logical(a$is_manipulation), ]
cat("\n--- gamma-set FDR<0.10 survivor counts ---\n")
cat("wilcox:", sum(g$varRatio_p_fdr < 0.10, na.rm = TRUE),
    " paired-t:", sum(g$varRatio_t_p_fdr < 0.10, na.rm = TRUE),
    " rawKS:", sum(g$KS_combined_p_fdr < 0.10, na.rm = TRUE),
    " centeredKS:", sum(g$KS_centered_combined_p_fdr < 0.10, na.rm = TRUE), "\n")
cat("min raw varRatio_p:", min(g$varRatio_p, na.rm = TRUE),
    " min varRatio_t_p:", min(g$varRatio_t_p, na.rm = TRUE), "\n")
cat("gamma metrics with 6/6 sign agreement on logVarRatio:",
    sum(g$n_pos == 6 | g$n_neg == 6, na.rm = TRUE), "\n")
cat("  higher_focus (n_pos==6):", sum(g$n_pos == 6, na.rm = TRUE),
    " higher_audiobook (n_neg==6):", sum(g$n_neg == 6, na.rm = TRUE), "\n")
cat("\n6/6 higher_focus (more variable in focus):\n")
print(g$metric[which(g$n_pos == 6)])
cat("6/6 higher_audiobook (more variable in audiobook):\n")
print(g$metric[which(g$n_neg == 6)])
