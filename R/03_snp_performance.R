# =============================================================================
# Step 3 - Prioritised SNPs: distribution (Table 2), individual and combined
# diagnostic performance (Table 3), and hierarchical screening (Table 4).
# Haplotype-level allele states are converted to individual-level TP/FN/FP/TN
# using the number of individuals carrying each haplotype. All values are
# within-sample (apparent) estimates.
# =============================================================================
source("R/00_setup.R")
d <- load_haptrace_data()
message("Step 3: SNP performance (target = ", TARGET, ")")
tgt_haps <- d$cnt$Target > 0

# ---- Table 2: allele states in target-population haplotypes -----------------
t2 <- do.call(rbind, lapply(PANEL, function(s) {
  hit <- carries(d, s); col <- d$mat[, snp_pos(s)]
  dominant <- names(which.max(tapply(d$cnt$All, col, sum)))
  states <- setNames(ifelse(hit[tgt_haps], snp_base(s), "-"), d$cnt$Haplotype[tgt_haps])
  data.frame(SNP = sub("[ACGT]$", "", s), Prioritised_allele = snp_base(s),
             Variant_count = sum(col != dominant), as.list(states),
             Target_haplotypes = sprintf("%d (%d)", sum(hit & tgt_haps), sum(d$cnt$Target[hit])),
             NonTarget_haplotypes = sprintf("%d (%d)", sum(hit & d$cnt$NonTarget > 0), sum(d$cnt$NonTarget[hit])),
             EI = round(sum(hit & tgt_haps) / sum(hit), 3), check.names = FALSE)
}))
write_table(t2, "Table2_SNP_haplotype_distribution.csv")

# ---- Table 3: individual SNPs, Halda-specific subset, combined panel --------
row3 <- function(label, hit) { k <- counts_for(d, hit); cbind(SNP = label, diag_metrics(k["tp"], k["fn"], k["fp"], k["tn"])) }
t3 <- rbind(do.call(rbind, lapply(PANEL, function(s) row3(s, carries(d, s)))),
            row3(sprintf("Halda-specific SNPs (%s)", paste(SPECIFIC, collapse = ", ")), carries_any(d, SPECIFIC)),
            row3("Combined five-SNP panel", carries_any(d, PANEL)))
write_table(t3, "Table3_diagnostic_performance.csv")

# ---- Table 4: hierarchical screening ----------------------------------------
# Samples negative at all preceding steps enter the next step; metrics are conditional on them.
remaining <- rep(TRUE, nrow(d$cnt)); cum <- 0; N_T <- sum(d$cnt$Target); rows <- list()
for (i in seq_along(PANEL)) {
  s <- PANEL[i]; pos <- remaining & carries(d, s); neg <- remaining & !pos
  tp <- sum(d$cnt$Target[pos]); fp <- sum(d$cnt$NonTarget[pos])
  fn <- sum(d$cnt$Target[neg]); tn <- sum(d$cnt$NonTarget[neg]); cum <- cum + tp
  m <- diag_metrics(tp, fn, fp, tn)
  rows[[i]] <- data.frame(Step = i, SNP = s,
    Entering = sprintf("%d (%d/%d)", tp + fn + fp + tn, tp + fn, fp + tn),
    TP = tp, FN = fn, FP = fp, TN = tn, Conditional_sensitivity = m$Sensitivity,
    Specificity = m$Specificity, Accuracy = m$Accuracy,
    Cumulative_target_detected = sprintf("%d (%.1f)", cum, 100 * cum / N_T),
    Detected_target_haplotypes = paste(d$cnt$Haplotype[pos & tgt_haps], collapse = ", "))
  remaining <- neg
}
t4 <- do.call(rbind, rows)
write_table(t4, "Table4_hierarchical_screening.csv")
print(t3[, c("SNP", "TP", "FN", "FP", "TN", "Sensitivity", "Specificity")], row.names = FALSE)
