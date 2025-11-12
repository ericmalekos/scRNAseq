#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(readr)
  library(Matrix)
  library(tximport)
})

# ---- CLI (unchanged option names) ----
opt <- OptionParser() |>
  add_option("--tx2gene",  help="Path to tx2gene TSV/CSV (read as-is)", metavar="file") |>
  add_option("--quant_dir", help="Parent directory; subdirs (clean_cell_id) each contain quant.sf", metavar="dir") |>
  add_option("--outdir",    help="Output directory", metavar="dir") |>
  parse_args()

if (is.null(opt$tx2gene) || is.null(opt$quant_dir) || is.null(opt$outdir)) {
  stop("Required: --tx2gene, --quant_dir, --outdir")
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Read tx2gene AS-IS (do not modify input file), build 2-col map for tximport ----
tx2g_in <- data.table::fread(opt$tx2gene, header = FALSE)  # auto-detects sep
names(tx2g_in) <- c("TXNAME_RAW", "GENEID_RAW")

if (ncol(tx2g_in) < 2L) stop("tx2gene must have at least two columns (first=transcript ID, second=gene ID)")
cat(sprintf("tx2gene loaded: %d rows, %d cols. Column names: %s\n",
            nrow(tx2g_in), ncol(tx2g_in), paste(names(tx2g_in), collapse=", ")))

# Create a two-column data.frame with names tximport expects (does not touch your file on disk)
tx2g_map <- data.frame(
  TXNAME = sub("\\.\\d+$", "", trimws(tx2g_in[[1]])),  # strip version suffix if present (safe if absent)
  GENEID = trimws(tx2g_in[[2]]),
  stringsAsFactors = FALSE
)
tx2g_map <- unique(tx2g_map)
cat(sprintf("Prepared tx2gene map: %d transcript→gene mappings\n", nrow(tx2g_map)))

# Save a tiny preview for sanity
utils::write.table(head(tx2g_in, 5),
                   file = file.path(opt$outdir, "tx2gene_preview_head5.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ---- Discover quant.sf files recursively; parent dir name is the clean_cell_id ----
all_qf <- list.files(opt$quant_dir, pattern = "^quant\\.sf$", recursive = TRUE, full.names = TRUE)
if (length(all_qf) == 0L) {
  stop("No 'quant.sf' files found under --quant_dir: ", opt$quant_dir)
}
cell_ids <- basename(dirname(all_qf))
if (any(duplicated(cell_ids))) {
  dups <- unique(cell_ids[duplicated(cell_ids)])
  stop("Duplicate cell IDs (parent dir names) detected: ", paste(head(dups, 20), collapse=", "),
       ". Ensure each cell subdir name is unique.")
}
qfiles <- stats::setNames(all_qf, cell_ids)
data.table::fwrite(
  data.table::data.table(clean_cell_id = names(qfiles), quant_sf = as.character(qfiles)),
  file.path(opt$outdir, "quant_sf_manifest.tsv"), sep = "\t"
)
cat(sprintf("Discovered %d quant.sf files. Example: %s -> %s\n",
            length(qfiles), names(qfiles)[1], qfiles[1]))

# ---- Identify transcript IDs missing from your mapping (compare to FIRST quant.sf) ----
# Assumes all cells were quantified against the same Salmon index (standard practice).
first_qf <- all_qf[1]
q_tx <- readr::read_tsv(first_qf, show_col_types = FALSE, col_types = cols_only(Name = col_character()))$Name
q_tx <- sub("\\.\\d+$", "", q_tx)  # align with TXNAME processing above
missing_in_tx2g <- setdiff(unique(q_tx), unique(tx2g_map$TXNAME))
extra_in_tx2g   <- setdiff(unique(tx2g_map$TXNAME), unique(q_tx))

# Write diagnostic lists (plain text, one ID per line), and print counts
writeLines(missing_in_tx2g, con = file.path(opt$outdir, "transcripts_missing_from_tx2gene.txt"))
writeLines(extra_in_tx2g,   con = file.path(opt$outdir, "tx2gene_tx_not_in_quant.txt"))
cat(sprintf("Transcripts in quant index but missing from tx2gene: %d (see transcripts_missing_from_tx2gene.txt)\n",
            length(missing_in_tx2g)))
cat(sprintf("TX IDs present in tx2gene but not in quant index: %d (see tx2gene_tx_not_in_quant.txt)\n",
            length(extra_in_tx2g)))

# ---- tximport in batches to avoid giant transcript x sample matrix ----

library(Matrix)

all_files <- qfiles                       # named vector: names = cell IDs, values = quant.sf paths
n_samples <- length(all_files)
cat(sprintf("Total samples (quant.sf files): %d\n", n_samples))

# Choose a batch size so that n_tx * batch_size << 2e9.
# With ~278k transcripts, 1000 is very safe (278M elements).
batch_size <- 1000L
idx <- seq_len(n_samples)
batches <- split(idx, ceiling(idx / batch_size))

counts_list <- vector("list", length(batches))
gene_order  <- NULL

for (b in seq_along(batches)) {
  this_idx   <- batches[[b]]
  this_files <- all_files[this_idx]
  
  cat(sprintf("Running tximport batch %d / %d with %d samples\n",
              b, length(batches), length(this_files)))
  
  txi_b <- tximport(
    files   = this_files,
    type    = "salmon",
    tx2gene = tx2g_map,
    countsFromAbundance = "no",
    ignoreTxVersion     = TRUE,
    dropInfReps         = TRUE   # important to avoid extra blow-up
  )
  
  # Gene x samples (dense for this batch)
  counts_b <- Matrix(txi_b$counts, sparse = TRUE)
  colnames(counts_b) <- names(this_files)
  
  # Ensure consistent gene order across batches
  if (is.null(gene_order)) {
    gene_order <- rownames(counts_b)
  } else {
    if (!identical(gene_order, rownames(counts_b))) {
      # Reorder if needed, just to be safe
      counts_b <- counts_b[gene_order, , drop = FALSE]
    }
  }
  
  counts_list[[b]] <- counts_b
  
  # Free memory from this batch
  txi_b$abundance <- NULL
  txi_b$length    <- NULL
  txi_b$infReps   <- NULL
  rm(txi_b, counts_b)
  gc()
}

# Column-bind all batches into a single genes x cells matrix
if (length(counts_list) == 1L) {
  counts_mat <- counts_list[[1L]]
} else {
  counts_mat <- do.call(cbind, counts_list)  # base cbind, works with dgCMatrix
}
rownames(counts_mat) <- gene_order
cat(sprintf("Final counts matrix: %d genes x %d cells\n",
            nrow(counts_mat), ncol(counts_mat)))

# ---- Save outputs (NO GZIP) ----
saveRDS(counts_mat, file = file.path(opt$outdir, "gene_counts.sparse.rds"))

Matrix::writeMM(counts_mat, file.path(opt$outdir, "gene_counts.mtx"))
data.table::fwrite(data.table(gene_id = rownames(counts_mat)),
                   file.path(opt$outdir, "genes.tsv"), sep = "\t", col.names = FALSE)
data.table::fwrite(data.table(cell_id = colnames(counts_mat)),
                   file.path(opt$outdir, "cells.tsv"), sep = "\t", col.names = FALSE)

cat("Wrote:\n")
cat("  - ", file.path(opt$outdir, "gene_counts.sparse.rds"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "gene_counts.mtx"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "genes.tsv"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "cells.tsv"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "quant_sf_manifest.tsv"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "transcripts_missing_from_tx2gene.txt"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "tx2gene_tx_not_in_quant.txt"), "\n", sep = "")

