# =============================================================================
# Step 5 - Supplementary Fig. S1: maximum-likelihood tree of the haplotypes
# represented by at least two individuals.
# (a) writes mega/CR_haplotypes_n2plus.fasta (input for MEGA-CC);
# (b) plots the tree in mega/CR_19haplotypes_ML.nwk, produced with
#     megacc -a mega/ML_T92G_1000bs.mao -d mega/CR_haplotypes_n2plus.fasta -o mega/CR_19haplotypes_ML
#     (Tamura 3-parameter + G, selected by BIC; 1000 bootstrap replicates; MEGA X/MEGA-CC 10.2.6).
# Requires the Bioconductor package 'ggtree' and CRAN 'ape' for plotting.
# =============================================================================
source("R/00_setup.R")
d <- load_haptrace_data()
message("Step 5: Fig. S1")

keep <- d$cnt$All >= 2
fa <- "mega/CR_haplotypes_n2plus.fasta"
writeLines(unlist(lapply(which(keep), function(i)
  c(paste0(">", d$cnt$Haplotype[i]), substring(d$seq[i], seq(1, d$L, 60), pmin(seq(60, d$L + 59, 60), d$L))))), fa)
message(sprintf("  wrote %s (%d haplotypes with >= 2 individuals)", fa, sum(keep)))

nwk <- "mega/CR_19haplotypes_ML.nwk"
if (!file.exists(nwk) || !all(vapply(c("ape", "ggtree", "ggplot2"), requireNamespace, TRUE, quietly = TRUE))) {
  message("  tree file or ape/ggtree/ggplot2 missing - plot skipped")
} else {
suppressPackageStartupMessages({ library(ape); library(ggtree); library(ggplot2) })
tr <- read.tree(nwk)
hap <- d$cnt[keep, ]
hap$group <- factor(ifelse(hap$Target == hap$All, sprintf("%s-private", TARGET),
                    ifelse(hap$Target > 0, sprintf("Shared with %s", TARGET), sprintf("Non-%s", TARGET))),
                    levels = c(sprintf("%s-private", TARGET), sprintf("Shared with %s", TARGET), sprintf("Non-%s", TARGET)))
hap$label2 <- ifelse(hap$Target > 0, sprintf("%s (n = %d; %s %d)", hap$Haplotype, hap$All, TARGET, hap$Target),
                     sprintf("%s (n = %d)", hap$Haplotype, hap$All))
names(hap)[names(hap) == "Haplotype"] <- "label"
bs <- suppressWarnings(as.numeric(tr$node.label)) * 100           # MEGA bootstrap proportions -> %
tr$node.label <- ifelse(!is.na(bs) & bs >= 50, sprintf("%.0f", bs), "")
tr <- ladderize(tr); xmax <- max(node.depth.edgelength(tr))
cols <- setNames(c("#C0392B", "#E59866", "#7F8C8D"), levels(hap$group))
p <- ggtree(tr, linewidth = 0.5) %<+% hap +
  geom_tippoint(aes(colour = group, size = All)) +
  geom_tiplab(aes(label = label2, colour = group), size = 3.3, offset = xmax * 0.03, show.legend = FALSE) +
  geom_nodelab(size = 3, hjust = 1.25, vjust = -0.45, colour = "black") +
  geom_treescale(x = 0, y = -0.3, width = 0.002, fontsize = 3, linesize = 0.4, offset = 0.35) +
  scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
  scale_size_continuous(range = c(1.8, 6), breaks = c(2, 5, 10, 30), name = "Individuals") +
  coord_cartesian(xlim = c(0, xmax * 1.75), ylim = c(-0.8, Ntip(tr) + 0.5), clip = "off") +
  theme(legend.position = c(0.86, 0.22), legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        legend.background = element_blank(), plot.margin = margin(8, 8, 8, 8)) +
  guides(colour = guide_legend(override.aes = list(size = 4), order = 1), size = guide_legend(order = 2))
ggsave(file.path(FIG_DIR, "FigS1_ML_tree.png"), p, width = 7.5, height = 6.5, dpi = 600, bg = "white")
ggsave(file.path(FIG_DIR, "FigS1_ML_tree.pdf"), p, width = 7.5, height = 6.5, bg = "white")
message("  wrote ", file.path(FIG_DIR, "FigS1_ML_tree.png/.pdf"))
}
