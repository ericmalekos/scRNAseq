#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
})

# ---------- CLI ----------
opt <- OptionParser() |>
  add_option("--counts_rds", help="RDS file with dgCMatrix of raw counts (genes x cells).", metavar="file") |>
  add_option("--outdir",     help="Output directory.", metavar="dir") |>
  add_option("--target_sum", type="double", default=1e6,
             help="Target library size per cell for normalization (e.g. 1e6 for CPM). [default %default]") |>
  add_option("--min_genes",  type="integer", default=500,
             help="Minimum number of detected genes per cell to keep. [default %default]") |>
  add_option("--min_counts", type="double", default=50000,
             help="Minimum total counts (reads) per cell to keep. [default %default]") |>
  parse_args()

if (is.null(opt$counts_rds) || is.null(opt$outdir)) {
  stop("Required: --counts_rds and --outdir")
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Load counts ----------
cat("Loading counts matrix from: ", opt$counts_rds, "\n", sep = "")
counts <- readRDS(opt$counts_rds)

if (!inherits(counts, "dgCMatrix")) {
  counts <- Matrix(counts, sparse = TRUE)
}
cat(sprintf("Raw counts matrix dimensions (before filtering): %d genes x %d cells\n",
            nrow(counts), ncol(counts)))

# ---------- Compute per-cell QC metrics (before filtering) ----------
cat("Computing per-cell QC metrics (before filtering)...\n")
lib_sizes_full       <- Matrix::colSums(counts)          # total reads per cell
n_genes_detected_full <- Matrix::colSums(counts > 0)     # genes with count > 0 per cell
cell_ids_full        <- if (!is.null(colnames(counts))) colnames(counts) else sprintf("cell_%d", seq_len(ncol(counts)))

cell_qc_full <- data.table(
  cell_id          = cell_ids_full,
  lib_size         = lib_sizes_full,
  n_genes_detected = n_genes_detected_full
)
fwrite(cell_qc_full, file.path(opt$outdir, "cell_qc_full_before_filtering.tsv"), sep = "\t")

# ---------- Apply cell filters ----------
keep_cells <- (n_genes_detected_full >= opt$min_genes) & (lib_sizes_full >= opt$min_counts)

cat(sprintf("Filtering cells using thresholds: min_genes >= %d, min_counts >= %.0f\n",
            opt$min_genes, opt$min_counts))
cat(sprintf("Keeping %d / %d cells (%.2f%%)\n",
            sum(keep_cells), length(keep_cells),
            100 * sum(keep_cells) / length(keep_cells)))

# Save dropped cells info
dropped_cells <- cell_qc_full[!keep_cells]
if (nrow(dropped_cells) > 0) {
  fwrite(dropped_cells,
         file.path(opt$outdir, "dropped_cells_qc.tsv"),
         sep = "\t")
  cat(sprintf("Dropped %d cells; details in dropped_cells_qc.tsv\n", nrow(dropped_cells)))
} else {
  cat("No cells dropped by QC thresholds.\n")
}

# Apply filter to counts
counts <- counts[, keep_cells, drop = FALSE]
cell_ids <- if (!is.null(colnames(counts))) colnames(counts) else sprintf("cell_%d", seq_len(ncol(counts)))

cat(sprintf("Counts matrix dimensions (after filtering): %d genes x %d cells\n",
            nrow(counts), ncol(counts)))

# ---------- Recompute detection stats AFTER filtering ----------
cat("Computing gene detection stats (after filtering)...\n")
n_cells_total <- ncol(counts)

# per-gene detection: counts > 0
n_cells_expressed <- Matrix::rowSums(counts > 0)
frac_cells        <- n_cells_expressed / n_cells_total

zero_cells  <- sum(n_cells_expressed == 0)
lt_1pct     <- sum(frac_cells > 0 & frac_cells < 0.01)
lt_5pct     <- sum(frac_cells > 0 & frac_cells < 0.05)
lt_10pct    <- sum(frac_cells > 0 & frac_cells < 0.10)

gene_ids <- if (!is.null(rownames(counts))) rownames(counts) else sprintf("gene_%d", seq_len(nrow(counts)))

# per-gene table
gene_det_dt <- data.table(
  gene_id               = gene_ids,
  n_cells_expressed     = n_cells_expressed,
  frac_cells_expressed  = frac_cells
)
fwrite(gene_det_dt,
       file.path(opt$outdir, "gene_detection_per_gene.tsv"),
       sep = "\t")

# summary table
gene_summary_dt <- data.table(
  metric = c("zero_cells",
             "expressed_in_<1pct_cells",
             "expressed_in_<5pct_cells",
             "expressed_in_<10pct_cells"),
  n_genes = c(zero_cells, lt_1pct, lt_5pct, lt_10pct),
  frac_genes = c(zero_cells, lt_1pct, lt_5pct, lt_10pct) / nrow(counts)
)
fwrite(gene_summary_dt,
       file.path(opt$outdir, "gene_detection_summary.tsv"),
       sep = "\t")

cat("Gene detection summary (after filtering):\n")
print(gene_summary_dt)

# ---------- Cell detection stats AFTER filtering ----------
cat("Computing cell detection stats (after filtering)...\n")
n_genes_detected <- Matrix::colSums(counts > 0)

pct_cells_lt_50  <- sum(n_genes_detected < 50)  / n_cells_total * 100
pct_cells_lt_500 <- sum(n_genes_detected < 500) / n_cells_total * 100

cell_det_dt <- data.table(
  cell_id          = cell_ids,
  n_genes_detected = n_genes_detected
)
fwrite(cell_det_dt,
       file.path(opt$outdir, "cell_detection_per_cell.tsv"),
       sep = "\t")

cell_summary_dt <- data.table(
  metric = c("cells_with_<50_genes_detected",
             "cells_with_<500_genes_detected"),
  n_cells = c(sum(n_genes_detected < 50),
              sum(n_genes_detected < 500)),
  frac_cells = c(pct_cells_lt_50, pct_cells_lt_500) / 100,
  pct_cells  = c(pct_cells_lt_50, pct_cells_lt_500)
)
fwrite(cell_summary_dt,
       file.path(opt$outdir, "cell_detection_summary.tsv"),
       sep = "\t")

cat("Cell detection summary (after filtering):\n")
print(cell_summary_dt)

# ---------- Normalization (CPM + log1p) ----------
cat("Normalizing counts (CP", opt$target_sum, " + log1p)...\n", sep = "")
lib_sizes <- Matrix::colSums(counts)

# avoid division by zero (should not happen after filtering, but just in case)
zero_lib <- which(lib_sizes == 0)
if (length(zero_lib) > 0L) {
  warning(length(zero_lib), " cells have library size 0 after filtering; setting their library size to 1.")
  lib_sizes[zero_lib] <- 1
}

scaling_factors <- lib_sizes / opt$target_sum  # larger lib size -> larger divisor
norm_counts <- t(t(counts) / scaling_factors)  # counts per target_sum (CPM if target_sum=1e6)

# log1p transform
log_norm_counts <- log1p(norm_counts)

# ---------- Save matrices ----------
saveRDS(counts,
        file.path(opt$outdir, "gene_counts.filtered.raw.rds"))
saveRDS(norm_counts,
        file.path(opt$outdir, "gene_counts.normalized_cp.rds"))
saveRDS(log_norm_counts,
        file.path(opt$outdir, "gene_counts.normalized_logcp.rds"))

cat("Wrote:\n")
cat("  - ", file.path(opt$outdir, "cell_qc_full_before_filtering.tsv"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "dropped_cells_qc.tsv"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "gene_detection_per_gene.tsv"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "gene_detection_summary.tsv"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "cell_detection_per_cell.tsv"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "cell_detection_summary.tsv"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "gene_counts.filtered.raw.rds"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "gene_counts.normalized_cp.rds"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "gene_counts.normalized_logcp.rds"), "\n", sep = "")
cat("Done.\n")
