# =============================================================================
# Step 2 - Control-region diversity (Table S1) and haplotype distribution (Table S2)
# Hd = n/(n-1) * (1 - sum(x_i^2)); pi = mean pairwise differences per site (Nei 1987).
# Singleton site: variable, with at most one nucleotide occurring more than once.
# Parsimony-informative site: at least two nucleotides each occurring at least twice.
# =============================================================================
source("R/00_setup.R")
d <- load_haptrace_data()
message("Step 2: diversity indices")

diversity <- function(w) {            # w = individuals per haplotype in the group
  keep <- w > 0; w <- w[keep]; m <- d$mat[keep, , drop = FALSE]; n <- sum(w); L <- d$L
  hd <- n / (n - 1) * (1 - sum((w / n)^2))
  pairs <- if (nrow(m) > 1) combn(nrow(m), 2) else matrix(0, 2, 0)
  diffs <- apply(pairs, 2, function(p) sum(m[p[1], ] != m[p[2], ]))
  pi <- sum(w[pairs[1, ]] * w[pairs[2, ]] * diffs) / choose(n, 2) / L
  site_counts <- lapply(seq_len(L), function(j) tapply(w, m[, j], sum))
  S   <- sum(vapply(site_counts, length, 1L) > 1)
  pis <- sum(vapply(site_counts, function(x) sum(x >= 2) >= 2, TRUE))
  c(`No of individuals` = n, `No of haplotypes` = length(w), `Haplotype diversity (Hd)` = round(hd, 3),
    `Nucleotide diversity (pi)` = round(pi, 5), `Monomorphic sites` = L - S, `Polymorphic sites (S)` = S,
    `Singleton sites` = S - pis, `Parsimony-informative sites` = pis)
}
groups <- c(setNames(lapply(POPS, function(p) d$cnt[[p]]), POPS), list(All = d$cnt$All))
s1 <- sapply(groups, diversity)
write_table(data.frame(Parameter = rownames(s1), s1, check.names = FALSE), "TableS1_diversity.csv")

s2 <- data.frame(Haplotype = d$cnt$Haplotype,
                 `Total %` = round(100 * d$cnt$All / sum(d$cnt$All), 2),
                 sapply(POPS, function(p) round(100 * d$cnt[[p]] / sum(d$cnt[[p]]), 2)), check.names = FALSE)
names(s2)[-(1:2)] <- paste(POPS, "%")
write_table(s2, "TableS2_haplotype_distribution.csv")
message(sprintf("  haplotypes: %d; with >=2 individuals: %d; singletons: %d", nrow(d$cnt),
                sum(d$cnt$All >= 2), sum(d$cnt$All == 1)))
