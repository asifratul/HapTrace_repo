# =============================================================================
# Step 4 - Batch-level acceptance rule (Section 4.7, Fig. 3)
# Rule: accept a batch as target (Halda) origin if at least k of n tested fish
# carry a target-specific allele. The number of positive fish is binomial, so
#   P(accept | target batch)     = P(X >= k | n, pT)
#   P(accept | non-target batch) = P(X >= k | n, pN)
# pT = proportion of target individuals carrying >= 1 target-specific allele
#      (with its 95% Clopper-Pearson interval);
# pN = false-positive rate; 0 positives were observed among non-target fish, so
#      0.5%, 1% and the upper 95% confidence limit (3.0%) are evaluated.
# No prior probability is needed; fish are assumed to be sampled independently.
# =============================================================================
source("R/00_setup.R")
d <- load_haptrace_data()
message("Step 4: batch acceptance rule")

hit <- carries_any(d, SPECIFIC)
x_T <- sum(d$cnt$Target[hit]);    N_T <- sum(d$cnt$Target)
x_N <- sum(d$cnt$NonTarget[hit]); N_N <- sum(d$cnt$NonTarget)
pT  <- x_T / N_T; ciT <- binom.test(x_T, N_T)$conf.int
fpU <- binom.test(x_N, N_N)$conf.int[2]                 # upper 95% limit of the false-positive rate
FP_RATES <- c(0.005, 0.01, fpU)
accept <- function(n, k, p) 1 - pbinom(k - 1, n, p)
n_for  <- function(k, p, target = 0.95) min(which(accept(1:500, k, p) >= target))

summary_tab <- data.frame(
  Quantity = c("Target individuals carrying >= 1 target-specific allele", "Detection rate pT (%)",
               "pT 95% CI lower (%)", "pT 95% CI upper (%)",
               "Non-target individuals carrying a target-specific allele",
               "False-positive rate, upper 95% CI (%)",
               "Fish for >=95% chance of >=1 positive (point / lower CI)",
               "Fish for >=99% chance of >=1 positive (point / lower CI)"),
  Value = c(sprintf("%d/%d", x_T, N_T), round(100 * pT, 1), round(100 * ciT[1], 1), round(100 * ciT[2], 1),
            sprintf("%d/%d", x_N, N_N), round(100 * fpU, 1),
            sprintf("%d / %d", n_for(1, pT), n_for(1, ciT[1])),
            sprintf("%d / %d", n_for(1, pT, 0.99), n_for(1, ciT[1], 0.99))))
write_table(summary_tab, "Batch_detection_summary.csv")

# operating characteristics of the rule for k = 1-3
rules <- do.call(rbind, lapply(1:3, function(k) {
  n <- n_for(k, pT); n_lo <- n_for(k, ciT[1])
  data.frame(Min_positives_k = k, Fish_tested_n = n,
             Accept_target_batch_pct = round(100 * accept(n, k, pT), 1),
             Accept_nontarget_pct_at_0.5 = round(100 * accept(n, k, FP_RATES[1]), 1),
             Accept_nontarget_pct_at_1.0 = round(100 * accept(n, k, FP_RATES[2]), 1),
             Accept_nontarget_pct_at_upperCI = round(100 * accept(n, k, FP_RATES[3]), 1),
             Fish_tested_n_conservative = n_lo)
}))
write_table(rules, "Batch_acceptance_rules.csv")

grid <- expand.grid(n = 1:50, k = 1:3)
grid$Accept_target <- accept(grid$n, grid$k, pT)
grid$Accept_target_lo <- accept(grid$n, grid$k, ciT[1]); grid$Accept_target_hi <- accept(grid$n, grid$k, ciT[2])
for (i in seq_along(FP_RATES)) grid[[sprintf("Accept_nontarget_fp%.1f", 100 * FP_RATES[i])]] <- accept(grid$n, grid$k, FP_RATES[i])
write_table(cbind(grid[, 1:2], round(100 * grid[, -(1:2)], 2)), "Batch_acceptance_curves.csv")

rec <- rules[rules$Min_positives_k == 2, ]
message(sprintf("  recommended rule: test %d fish, accept if >= 2 positive (conservative: %d fish); false acceptance %.1f / %.1f / %.1f %%",
                rec$Fish_tested_n, rec$Fish_tested_n_conservative, rec$Accept_nontarget_pct_at_0.5,
                rec$Accept_nontarget_pct_at_1.0, rec$Accept_nontarget_pct_at_upperCI))

# ---- Fig. 3 -------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  k <- 2; n <- 1:50; n95 <- rec$Fish_tested_n; n95lo <- rec$Fish_tested_n_conservative
  halda <- data.frame(n = n, y = 100 * accept(n, k, pT), lo = 100 * accept(n, k, ciT[1]), hi = 100 * accept(n, k, ciT[2]))
  non <- do.call(rbind, lapply(FP_RATES, function(p)
    data.frame(n = n, y = 100 * accept(n, k, p), rate = sprintf("Non-%s batch (false-positive rate %.1f%%)", TARGET, 100 * p))))
  non$rate <- factor(non$rate, levels = unique(non$rate))
  lab_t <- sprintf("%s-origin batch", TARGET)
  p <- ggplot() +
    geom_ribbon(data = halda, aes(n, ymin = lo, ymax = hi), fill = "#0072C6", alpha = 0.15) +
    geom_line(data = halda, aes(n, y, colour = lab_t), linewidth = 1.2) +
    geom_point(data = halda, aes(n, y, colour = lab_t), size = 1.8) +
    geom_line(data = non, aes(n, y, colour = rate), linewidth = 0.9, linetype = "longdash") +
    geom_hline(yintercept = 95, linetype = "dotted", colour = "grey40") +
    geom_vline(xintercept = n95, linetype = "dashed", colour = "#D7263D", linewidth = 0.7) +
    geom_vline(xintercept = n95lo, linetype = "dotted", colour = "#D7263D", linewidth = 0.7) +
    annotate("text", x = n95 + 0.7, y = 55, hjust = 0, colour = "#D7263D", size = 3.8,
             label = sprintf("%d fish: accept if ≥2 positive", n95)) +
    annotate("text", x = n95lo + 0.7, y = 62, hjust = 0, colour = "#D7263D", size = 3.4,
             label = sprintf("%d fish (conservative)", n95lo)) +
    scale_colour_manual(NULL, values = c(setNames("#0072C6", lab_t), setNames(c("grey55", "#E69F00", "#B22222"), levels(non$rate))),
                        breaks = c(lab_t, levels(non$rate))) +
    guides(colour = guide_legend(ncol = 2, override.aes = list(linetype = c("solid", rep("longdash", 3)), shape = c(16, NA, NA, NA)))) +
    scale_x_continuous(breaks = seq(0, 50, 5), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    labs(x = "Number of fish tested", y = sprintf("Probability of accepting the batch as %s origin (%%)", TARGET)) +
    theme_bw(base_size = 13) +
    theme(axis.title = element_text(face = "bold"), panel.grid.minor = element_blank(), legend.position = "bottom",
          legend.text = element_text(size = 9.5), legend.key.width = unit(1.6, "lines"))
  ggsave(file.path(FIG_DIR, "Fig3_batch_acceptance_rule.png"), p, width = 8, height = 6, dpi = 600)
  ggsave(file.path(FIG_DIR, "Fig3_batch_acceptance_rule.pdf"), p, width = 8, height = 6)
  message("  wrote ", file.path(FIG_DIR, "Fig3_batch_acceptance_rule.png/.pdf"))
} else message("  ggplot2 not installed - Fig. 3 skipped")
