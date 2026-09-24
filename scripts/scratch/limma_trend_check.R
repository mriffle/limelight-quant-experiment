# Reference fits for the limma-trend validation (called by validate_limma_trend.py).
#
# Usage: Rscript limma_trend_check.R <in_dir> <out_dir>
#   <in_dir>/samples.tsv            sample_id, condition, candidate_pair, batch
#   <in_dir>/<quantity>.tsv         feature x sample log2 matrix (first column = feature)
# Writes <out_dir>/<quantity>_<design>.tsv with per-feature coef, t, p, BH q, s2.prior,
# Amean, and the scalar df.prior / df.residual, for lmFit + eBayes(trend = TRUE)
# (and eBayes(trend = FALSE) for reference), designs paired / batch / unadjusted.
suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
in_dir <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- read.delim(file.path(in_dir, "samples.tsv"), quote = "", comment.char = "",
                      check.names = FALSE, stringsAsFactors = FALSE)
condition <- factor(samples$condition, levels = c("control", "raloxifene-d0"))
pair <- factor(samples$candidate_pair)
batch <- factor(samples$batch)
designs <- list(
  paired = model.matrix(~ condition + pair),
  batch = model.matrix(~ condition + batch),
  unadjusted = model.matrix(~ condition)
)
fmt <- function(x) sprintf("%.17g", x)

for (qfile in list.files(in_dir, pattern = "\\.tsv$")) {
  if (qfile == "samples.tsv") next
  quantity <- sub("\\.tsv$", "", qfile)
  m <- read.delim(file.path(in_dir, qfile), quote = "", comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE)
  y <- as.matrix(m[, samples$sample_id])
  rownames(y) <- m$feature
  for (d in names(designs)) {
    fit <- lmFit(y, designs[[d]])
    tr <- eBayes(fit, trend = TRUE)
    nt <- eBayes(fit, trend = FALSE)
    out <- data.frame(
      feature = rownames(y),
      coef = fmt(tr$coefficients[, 2]),
      t = fmt(tr$t[, 2]),
      p = fmt(tr$p.value[, 2]),
      q = fmt(p.adjust(tr$p.value[, 2], method = "BH")),
      s2_prior = fmt(tr$s2.prior),
      amean = fmt(tr$Amean),
      df_prior = fmt(tr$df.prior),
      df_residual = fmt(tr$df.residual),
      notrend_t = fmt(nt$t[, 2]),
      notrend_p = fmt(nt$p.value[, 2]),
      notrend_df_prior = fmt(nt$df.prior),
      stringsAsFactors = FALSE
    )
    write.table(out, file.path(out_dir, paste0(quantity, "_", d, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
cat("limma", as.character(packageVersion("limma")), "\n")
