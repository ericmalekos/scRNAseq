#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
})

# ---------- CLI ----------
opt <- OptionParser() |>
  add_option("--counts_rds",     help="RDS with dgCMatrix of raw counts (genes x cells).", metavar="file") |>
  add_option("--logcounts_rds",  help="RDS with dgCMatrix of normalized logCP counts (same dims as counts).", metavar="file") |>
  add_option("--fastq_manifest", help="TSV/CSV with obs_names and clean_cell_id columns.", metavar="file") |>
  add_option("--metadata",       help="TSV/CSV with per-cell metadata including obs_names.", metavar="file") |>
  add_option("--outdir",         help="Output directory.", metavar="dir") |>
  add_option("--output_rds",     help="Name of output SCE RDS file [default sce_plate_seq.rds].",
             default="sce_plate_seq.rds") |>
  parse_args()

if (is.null(opt$counts_rds) || is.null(opt$logcounts_rds) ||
    is.null(opt$fastq_manifest) || is.null(opt$metadata) ||
    is.null(opt$outdir)) {
  stop("Required: --counts_rds, --logcounts_rds, --fastq_manifest, --metadata, --outdir")
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Load matrices ----------
cat("Loading counts from: ", opt$counts_rds, "\n", sep = "")
counts <- readRDS(opt$counts_rds)
if (!inherits(counts, "dgCMatrix")) {
  counts <- Matrix(counts, sparse = TRUE)
}

cat("Loading logcounts from: ", opt$logcounts_rds, "\n", sep = "")
logcounts <- readRDS(opt$logcounts_rds)
if (!inherits(logcounts, "dgCMatrix")) {
  logcounts <- Matrix(logcounts, sparse = TRUE)
}

# basic consistency checks
if (!identical(dim(counts), dim(logcounts))) {
  stop("Dimensions of counts and logcounts do not match.")
}
if (!identical(rownames(counts), rownames(logcounts))) {
  stop("Row (gene) names of counts and logcounts do not match.")
}
if (!identical(colnames(counts), colnames(logcounts))) {
  stop("Column (cell) names of counts and logcounts do not match.")
}

cat(sprintf("Matrix dimensions: %d genes x %d cells\n",
            nrow(counts), ncol(counts)))

cell_ids <- colnames(counts)

# ---------- Load fastq manifest (obs_names <-> clean_cell_id) ----------
cat("Loading fastq manifest from: ", opt$fastq_manifest, "\n", sep = "")
fq <- fread(opt$fastq_manifest)
if (!all(c("obs_names","clean_cell_id") %in% names(fq))) {
  stop("fastq_manifest must contain columns: obs_names, clean_cell_id")
}

# ensure unique mapping from clean_cell_id to obs_names
if (any(duplicated(fq$clean_cell_id))) {
  dups <- fq[duplicated(clean_cell_id), unique(clean_cell_id)]
  warning("Duplicate clean_cell_id entries in fastq_manifest; keeping first occurrence. Examples: ",
          paste(head(dups, 10), collapse = ", "))
  fq <- fq[!duplicated(clean_cell_id)]
}

# ---------- Load full metadata (with obs_names) ----------
cat("Loading metadata from: ", opt$metadata, "\n", sep = "")
meta <- fread(opt$metadata)
if (!"obs_names" %in% names(meta)) {
  stop("metadata file must contain an 'obs_names' column.")
}

# ---------- Join metadata <- fastq_manifest to get clean_cell_id ----------
meta_join <- merge(meta, fq[, .(obs_names, clean_cell_id)],
                   by = "obs_names", all.x = TRUE)

n_missing_clean <- sum(is.na(meta_join$clean_cell_id))
if (n_missing_clean > 0) {
  warning(n_missing_clean, " rows in metadata had no matching clean_cell_id from fastq_manifest.")
}

# ---------- Build colData aligned to counts ----------
# start from the vector of cell_ids (columns of counts)
cdt <- data.table(clean_cell_id = cell_ids)

# left join metadata onto this vector
cdt <- merge(cdt, meta_join, by = "clean_cell_id", all.x = TRUE, sort = FALSE)

# confirm row order matches cell_ids
if (!identical(cdt$clean_cell_id, cell_ids)) {
  # enforce order
  cdt <- cdt[match(cell_ids, clean_cell_id)]
}

# number of cells with completely missing metadata (aside from clean_cell_id)
meta_only_cols <- setdiff(names(cdt), "clean_cell_id")
n_all_na_meta <- sum(apply(as.data.frame(cdt[, ..meta_only_cols]), 1,
                           function(x) all(is.na(x))))
if (n_all_na_meta > 0) {
  warning(n_all_na_meta, " cells have no metadata (all NA except clean_cell_id).")
}

coldata_df <- as.data.frame(cdt)
rownames(coldata_df) <- coldata_df$clean_cell_id

cat("colData columns:\n")
print(colnames(coldata_df))

# ---------- Build SingleCellExperiment ----------
sce <- SingleCellExperiment(
  assays = list(
    counts    = counts,
    logcounts = logcounts
  ),
  colData = coldata_df
)

# add gene IDs explicitly to rowData (if rownames already are gene IDs, this is redundant but harmless)
if (!is.null(rownames(counts))) {
  rowData(sce)$gene_id <- rownames(counts)
}

# ---------- Save ----------
out_rds <- file.path(opt$outdir, opt$output_rds)
saveRDS(sce, out_rds)

cat("Saved SingleCellExperiment object to:\n  ", out_rds, "\n", sep = "")
cat("Done.\n")
