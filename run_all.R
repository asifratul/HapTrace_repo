# =============================================================================
# HapTrace - run the full analysis (from the repository root):
#   Rscript run_all.R
# Each step can also be run on its own, e.g. Rscript R/03_snp_performance.R
# =============================================================================
steps <- c("R/01_snp_matrix_ei.R", "R/02_diversity.R", "R/03_snp_performance.R",
           "R/04_batch_rule.R", "R/05_figS1_tree.R")
for (s in steps) {
  message("\n=== ", s, " ===")
  local(source(s, local = TRUE))
}
message("\nDone. Tables in results/tables, figures in results/figures.")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
