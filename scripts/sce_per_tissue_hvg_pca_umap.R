#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(scran)
  library(scater)
  library(igraph)
  library(BiocSingular)
})

# ---------- CLI ----------
opt <- OptionParser() |>
  add_option("--sce_rds", help="RDS file with SingleCellExperiment (counts + logcounts + metadata).", metavar="file") |>
  add_option("--outdir",  help="Output directory for results and per-tissue SCEs.", metavar="dir") |>
  
  # Optional: only process these tissues (comma-separated)
  add_option("--tissues", type = "character", default = NULL,
             help = "Optional comma-separated list of tissue names to process (default: all tissues).") |>
  
  # Global gene filter
  add_option("--min_frac_cells", type="double", default=0.001,
             help="Global gene filter: keep genes expressed in >= this fraction of cells [default %default].") |>
  add_option("--min_mean_logexpr", type="double", default=0,
             help="Global gene filter: keep genes with mean(logcounts) > this value [default %default].") |>
  
  # Per-tissue HVG / PCA / UMAP / clustering params
  add_option("--hvg_bio_threshold", type="double", default=0.5,
             help="Per-tissue HVG rule: keep genes with biological variance (bio) > this value [default %default].") |>
  add_option("--min_hvg", type="integer", default=1000,
             help="Minimum number of HVGs per tissue (fallback if threshold yields fewer) [default %default].") |>
  add_option("--max_hvg", type="integer", default=3000,
             help="Maximum number of HVGs per tissue (cap for getTopHVGs) [default %default].") |>
  add_option("--k_snn", type="integer", default=20,
             help="k for SNN graph construction [default %default].") |>
  add_option("--n_pcs", type="integer", default=50,
             help="Maximum number of PCs to compute per tissue with runPCA [default %default].") |>
  
  # UMAP params (user-controlled)
  add_option("--umap_neighbors", type="integer", default=15,
             help="Target n_neighbors for per-tissue UMAP (capped at n_cells-1) [default %default].") |>
  add_option("--umap_min_dist", type="double", default=0.3,
             help="min_dist for per-tissue UMAP [default %default].") |>
  parse_args()

if (is.null(opt$sce_rds) || is.null(opt$outdir)) {
  stop("Required: --sce_rds and --outdir")
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Load SCE ----------
cat("Loading SCE from: ", opt$sce_rds, "\n", sep = "")
sce <- readRDS(opt$sce_rds)
stopifnot(inherits(sce, "SingleCellExperiment"))

if (!"counts" %in% assayNames(sce)) {
  stop("SCE must contain a 'counts' assay.")
}
if (!"logcounts" %in% assayNames(sce)) {
  stop("SCE must contain a 'logcounts' assay (log-normalized counts).")
}
if (!"tissue" %in% colnames(colData(sce))) {
  stop("colData(sce) must contain 'tissue'.")
}
if (!"cell_ontology_class" %in% colnames(colData(sce))) {
  stop("colData(sce) must contain 'cell_ontology_class' for filtering.")
}

cat(sprintf("Initial SCE dimensions: %d genes x %d cells\n",
            nrow(sce), ncol(sce)))

# ---------- Drop NA/'nan' tissues ----------
n_before <- ncol(sce)
sce <- sce[, !is.na(sce$tissue) & sce$tissue != "nan", drop = FALSE]
n_after <- ncol(sce)
cat(sprintf("Dropped %d cells with NA/'nan' tissue. Remaining: %d cells.\n",
            n_before - n_after, n_after))

# ---------- Drop NA/'nan' cell_ontology_class ----------
n_before2 <- ncol(sce)
sce <- sce[, !is.na(sce$cell_ontology_class) & sce$cell_ontology_class != "nan", drop = FALSE]
n_after2 <- ncol(sce)
cat(sprintf("Dropped %d cells with NA/'nan' cell_ontology_class. Remaining: %d cells.\n",
            n_before2 - n_after2, n_after2))

# ---------- Global gene filter ----------
cat("Applying global gene filter...\n")
cnt <- assay(sce, "counts")
lgc <- assay(sce, "logcounts")

n_cells <- ncol(sce)
det_frac <- Matrix::rowSums(cnt > 0) / n_cells
mean_logexpr <- Matrix::rowMeans(lgc)

keep_genes <- det_frac >= opt$min_frac_cells & mean_logexpr > opt$min_mean_logexpr

cat(sprintf("Keeping %d / %d genes after global filter (>=%.3f frac cells, mean logexpr > %.3f)\n",
            sum(keep_genes), length(keep_genes),
            opt$min_frac_cells, opt$min_mean_logexpr))

sce <- sce[keep_genes, , drop = FALSE]

# ---------- Prepare per-tissue analysis ----------
all_tissues <- sort(unique(as.character(sce$tissue)))
all_tissues <- all_tissues[!is.na(all_tissues) & all_tissues != "nan"]

cat("Found ", length(all_tissues), " distinct tissues after filtering.\n", sep = "")
cat("All tissues in object:\n")
cat("  - ", paste(all_tissues, collapse = "\n  - "), "\n", sep = "")

# If user passed a comma-separated list, restrict to those
if (!is.null(opt$tissues)) {
  requested <- strsplit(opt$tissues, ",")[[1]]
  requested <- trimws(requested)
  requested <- requested[nchar(requested) > 0]
  
  cat("\nUser requested tissues:\n")
  cat("  - ", paste(requested, collapse = "\n  - "), "\n", sep = "")
  
  missing <- setdiff(requested, all_tissues)
  if (length(missing) > 0) {
    warning("The following requested tissues are not present in the data and will be ignored: ",
            paste(missing, collapse = ", "))
  }
  
  tissues <- intersect(all_tissues, requested)
  if (length(tissues) == 0) {
    stop("After filtering by requested tissues, no tissues remain to process.")
  }
  
  cat("\nTissues that will actually be processed:\n")
  cat("  - ", paste(tissues, collapse = "\n  - "), "\n", sep = "")
} else {
  tissues <- all_tissues
  cat("\nNo --tissues supplied; processing ALL tissues.\n")
}

summary_tab <- data.table(
  tissue  = character(),
  n_cells = integer(),
  n_genes = integer(),
  n_hvgs  = integer()
)

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (nchar(x) == 0L) x <- "NA"
  x
}

# ---------- Loop over tissues ----------
for (tiss in tissues) {
  cat("\n=== Processing tissue: ", tiss, " ===\n", sep = "")
  sce_sub <- sce[, sce$tissue == tiss, drop = FALSE]
  cat(sprintf("  Subset: %d genes x %d cells\n", nrow(sce_sub), ncol(sce_sub)))
  
  # skip tiny subsets
  if (ncol(sce_sub) < 10) {
    cat("  Fewer than 10 cells; skipping PCA/UMAP/clustering for this tissue.\n")
    next
  }
  
  # drop genes never expressed in this tissue
  cnt_sub <- assay(sce_sub, "counts")
  keep_genes_tiss <- Matrix::rowSums(cnt_sub > 0) > 0
  sce_sub <- sce_sub[keep_genes_tiss, , drop = FALSE]
  cat(sprintf("  After dropping never-expressed genes: %d genes x %d cells\n",
              nrow(sce_sub), ncol(sce_sub)))
  
  # HVG selection with per-tissue modelGeneVar
  cat("  Fitting mean-variance trend and selecting HVGs...\n")
  dec <- scran::modelGeneVar(sce_sub, assay.type = "logcounts")
  hvg <- rownames(dec)[dec$bio > opt$hvg_bio_threshold]
  cat(sprintf("  HVGs with bio > %.3f: %d genes\n", opt$hvg_bio_threshold, length(hvg)))
  
  if (length(hvg) < opt$min_hvg) {
    n_top <- min(opt$max_hvg, nrow(dec))
    cat(sprintf("  Fewer than %d HVGs; taking top %d HVGs by bio.\n", opt$min_hvg, n_top))
    hvg <- scran::getTopHVGs(dec, n = n_top)
  }
  
  rowData(sce_sub)$is_hvg <- rownames(sce_sub) %in% hvg
  
  # ---------- PCA per tissue ----------
  cat("  Running PCA on HVGs (runPCA, ExactParam)...\n")
  max_rank <- min(
    opt$n_pcs,
    ncol(sce_sub) - 1L,
    length(hvg) - 1L
  )
  if (max_rank < 2L) {
    cat("  Not enough rank (max_rank < 2); skipping PCA/UMAP/clustering for this tissue.\n")
    next
  }
  
  set.seed(42)
  sce_sub <- scater::runPCA(
    sce_sub,
    exprs_values = "logcounts",
    subset_row   = hvg,
    ncomponents  = max_rank,
    BSPARAM      = BiocSingular::ExactParam()
  )
  if (!"PCA" %in% names(reducedDims(sce_sub))) {
    stop("runPCA did not produce a 'PCA' reducedDim.")
  }
  pca_mat <- reducedDim(sce_sub, "PCA")
  cat(sprintf("  PCA completed: %d PCs\n", ncol(pca_mat)))
  
  # ---------- UMAP per tissue ----------
  umap_neighbors <- min(as.integer(opt$umap_neighbors), ncol(sce_sub) - 1L)
  if (umap_neighbors < 2L) {
    cat("  Too few cells for UMAP (effective n_neighbors < 2); skipping UMAP and clustering.\n")
    next
  }
  
  cat("  Running UMAP on PCA embedding (n_neighbors = ", umap_neighbors,
      ", min_dist = ", opt$umap_min_dist, ")...\n", sep = "")
  set.seed(42)
  sce_sub <- scater::runUMAP(
    sce_sub,
    dimred      = "PCA",
    name        = "UMAP",
    n_neighbors = umap_neighbors,
    min_dist    = opt$umap_min_dist
  )
  
  # ---------- SNN graph + Louvain clustering per tissue ----------
  cat("  Building SNN graph (k = ", opt$k_snn, ") and clustering (Louvain)...\n", sep = "")
  k_effective <- min(opt$k_snn, ncol(sce_sub) - 1L)
  if (k_effective < 2L) {
    cat("  Too few cells for meaningful SNN clustering (effective k < 2); skipping clustering.\n")
  } else {
    g <- scran::buildSNNGraph(sce_sub, use.dimred = "PCA", k = k_effective)
    cl <- igraph::cluster_louvain(g)$membership
    colData(sce_sub)$cluster_louvain <- factor(cl)
    cat("  Clustering completed: ", length(unique(cl)),
        " clusters for this tissue.\n", sep = "")
  }
  
  summary_tab <- rbind(
    summary_tab,
    data.table(
      tissue  = tiss,
      n_cells = ncol(sce_sub),
      n_genes = nrow(sce_sub),
      n_hvgs  = length(hvg)
    )
  )
  
  tiss_sanitized <- sanitize_name(tiss)
  out_file <- file.path(opt$outdir, paste0("sce_tissue_", tiss_sanitized, ".rds"))
  saveRDS(sce_sub, out_file)
  cat("  Saved per-tissue SCE to: ", out_file, "\n", sep = "")
}

# ---------- Save summary ----------
if (nrow(summary_tab) > 0) {
  fwrite(summary_tab,
         file.path(opt$outdir, "tissue_summary.tsv"),
         sep = "\t")
  cat("\nSummary across tissues written to tissue_summary.tsv\n")
}

cat("\nPer-tissue SCEs saved with prefix 'sce_tissue_'.\n")
