# =============================================================================
# Step 1 - SNP-haplotype matrix and Exclusivity Index (EI)
# Output: results/tables/SNP_haplotype_matrix.csv, results/tables/Exclusivity_index_all_alleles.csv
# EI = (target-population haplotypes carrying the allele) / (all haplotypes carrying the allele);
# each haplotype is counted once, irrespective of the number of individuals carrying it.
# =============================================================================
source("R/00_setup.R")
d <- load_haptrace_data()

var_sites <- which(apply(d$mat, 2, function(col) length(unique(col[col %in% c("A", "C", "G", "T")])) > 1))
message("Step 1: ", nrow(d$mat), " haplotypes, ", d$L, " bp, ", length(var_sites), " polymorphic sites")

# SNP-haplotype matrix (rows = polymorphic sites, columns = haplotypes)
snp_matrix <- data.frame(Position_manuscript = var_sites + POS_OFFSET,
                         Position_alignment = var_sites,
                         t(d$mat[, var_sites, drop = FALSE]), check.names = FALSE)
write_table(snp_matrix, "SNP_haplotype_matrix.csv")

# EI and individual counts for every allele at every polymorphic site
is_target_hap <- d$cnt$Target > 0
ei <- do.call(rbind, lapply(var_sites, function(j) {
  col <- d$mat[, j]
  w <- tapply(d$cnt$All, col, sum)                     # individuals per allele
  dominant <- names(w)[which.max(w)]
  do.call(rbind, lapply(sort(unique(col)), function(b) {
    has <- col == b
    data.frame(Allele = paste0(j + POS_OFFSET, b), Position_manuscript = j + POS_OFFSET,
               Position_alignment = j, Base = b, Majority_allele = dominant,
               Variant_count = sum(col != dominant),        # haplotypes carrying a non-majority allele
               Target_haplotypes = sum(has & is_target_hap), All_haplotypes = sum(has),
               EI = round(sum(has & is_target_hap) / sum(has), 3),
               Target_individuals = sum(d$cnt$Target[has]), NonTarget_individuals = sum(d$cnt$NonTarget[has]),
               Target_haplotype_names = paste(d$cnt$Haplotype[has & is_target_hap], collapse = ";"))
  }))
}))
ei <- ei[order(-ei$EI, -ei$Target_individuals, ei$Position_manuscript), ]
write_table(ei, "Exclusivity_index_all_alleles.csv")
message("  alleles with EI = 1: ", paste(ei$Allele[ei$EI == 1], collapse = ", "))
