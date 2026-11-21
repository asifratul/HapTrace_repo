# =============================================================================
# HapTrace - shared settings and helper functions
# Sourced by every analysis script (run from the repository root).
# =============================================================================

# ---- settings ---------------------------------------------------------------
TARGET     <- "Halda"                                   # target population
POPS       <- c("Padma", "Jamuna", "Halda", "Culture")  # population columns in the count table
POS_OFFSET <- 4     # SNP labels use manuscript coordinates = position in data/CR_haplotypes.fasta + 4
PANEL      <- c("460C", "171A", "131T", "207C", "292C") # prioritised SNPs (screening order)
SPECIFIC   <- c("171A", "131T", "207C", "292C")         # Halda-specific subset (EI = 1)

FASTA_FILE  <- "data/CR_haplotypes.fasta"     # 50 control-region haplotypes, 594 bp
COUNTS_FILE <- "data/haplotype_counts.csv"    # individuals per haplotype and population (n = 166)
TAB_DIR     <- "results/tables"
FIG_DIR     <- "results/figures"
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- input ------------------------------------------------------------------
read_fasta <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x[nzchar(trimws(x))])
  hdr <- grepl("^>", x)
  id  <- sub("^>\\s*(\\S+).*", "\\1", x[hdr])
  seq <- vapply(split(x[!hdr], cumsum(hdr)[!hdr]), paste, "", collapse = "")
  setNames(toupper(seq), id)
}

load_haptrace_data <- function() {
  seq <- read_fasta(FASTA_FILE)
  cnt <- read.csv(COUNTS_FILE, stringsAsFactors = FALSE)
  names(cnt)[1] <- "Haplotype"
  stopifnot(setequal(names(seq), cnt$Haplotype), length(unique(nchar(seq))) == 1,
            all(rowSums(cnt[, POPS]) == cnt$All))
  cnt <- cnt[match(names(seq), cnt$Haplotype), ]          # same order as the FASTA
  mat <- do.call(rbind, strsplit(seq, ""))                  # haplotypes x sites
  rownames(mat) <- names(seq)
  cnt$Target    <- cnt[[TARGET]]
  cnt$NonTarget <- cnt$All - cnt$Target
  list(seq = seq, mat = mat, cnt = cnt, L = ncol(mat))
}

# ---- SNP helpers ------------------------------------------------------------
snp_pos  <- function(snp) as.integer(sub("[ACGT]$", "", snp)) - POS_OFFSET   # alignment column
snp_base <- function(snp) sub("^[0-9]+", "", snp)
carries  <- function(d, snp) d$mat[, snp_pos(snp)] == snp_base(snp)            # logical, per haplotype
carries_any <- function(d, snps) Reduce(`|`, lapply(snps, function(s) carries(d, s)))

# ---- diagnostic metrics -----------------------------------------------------
pct <- function(a, b) ifelse(b > 0, round(100 * a / b, 1), NA)
diag_metrics <- function(tp, fn, fp, tn) {
  data.frame(TP = tp, FN = fn, FP = fp, TN = tn,
             Sensitivity = pct(tp, tp + fn), Specificity = pct(tn, tn + fp),
             Accuracy = pct(tp + tn, tp + fn + fp + tn),
             PPV = pct(tp, tp + fp), NPV = pct(tn, tn + fn))
}
counts_for <- function(d, hit) {                    # hit = logical per haplotype
  tp <- sum(d$cnt$Target[hit]);    fn <- sum(d$cnt$Target[!hit])
  fp <- sum(d$cnt$NonTarget[hit]); tn <- sum(d$cnt$NonTarget[!hit])
  c(tp = tp, fn = fn, fp = fp, tn = tn)
}

write_table <- function(x, name) {
  write.csv(x, file.path(TAB_DIR, name), row.names = FALSE)
  message("  wrote ", file.path(TAB_DIR, name))
}
