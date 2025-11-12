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

# ===================== CLI =====================
opt <- OptionParser() |>
  add_option("--sce_rds", help="Input SCE RDS (global object).", metavar = "file") |>
  add_option("--outdir",  help="Output directory.", metavar = "dir") |>
  
  # optional global gene filter (set min_frac_cells=0 to disable)
  add_option("--min_frac_cells", type = "double", default = 0.01,
             help = "Global gene filter: keep genes expressed in >= this fraction of cells [default %default].") |>
  add_option("--min_mean_logexpr", type = "double", default = 0,
             help = "Global gene filter: keep genes with mean(logcounts) > this value [default %default].") |>
  
  # HVG / PCA params
  add_option("--n_hvg",    type = "integer", default = 3000,
             help = "Number of HVGs to use for PCA/UMAP [default %default].") |>
  add_option("--n_pcs",    type = "integer", default = 50,
             help = "Maximum number of PCs for PCA [default %default].") |>
  
  # UMAP params (more separation)
  add_option("--umap_neighbors", type = "integer", default = 10,
             help = "n_neighbors for UMAP [default %default].") |>
  add_option("--umap_min_dist",  type = "double",  default = 0.1,
             help = "min_dist for UMAP [default %default].") |>
  add_option("--umap_name",      type = "character", default = "UMAP_refined",
             help = "Name of reducedDim slot for UMAP [default %default].") |>
  
  # Graph + clustering
  add_option("--k_snn", type = "integer", default = 20,
             help = "k for SNN graph [default %default].") |>
  add_option("--cluster_method", type = "character", default = "louvain",
             help = "Clustering method: 'louvain' or 'leiden' [default %default].") |>
  parse_args()

if (is.null(opt$sce_rds) || is.null(opt$outdir)) {
  stop("Required: --sce_rds and --outdir")
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ===================== Load SCE =====================
cat("Loading SCE from: ", opt$sce_rds, "\n", sep = "")
sce <- readRDS(opt$sce_rds)
stopifnot(inherits(sce, "SingleCellExperiment"))

if (!"counts" %in% assayNames(sce)) {
  stop("SCE must contain a 'counts' assay.")
}
if (!"logcounts" %in% assayNames(sce)) {
  stop("SCE must contain a 'logcounts' assay (log-normalized counts).")
}

cat(sprintf("Initial SCE dimensions: %d genes x %d cells\n",
            nrow(sce), ncol(sce)))

# ===================== Global gene filter =====================
if (opt$min_frac_cells > 0 || opt$min_mean_logexpr > 0) {
  cat("Applying global gene filter...\n")
  cnt <- assay(sce, "counts")
  lgc <- assay(sce, "logcounts")
  
  n_cells <- ncol(sce)
  det_frac <- Matrix::rowSums(cnt > 0) / n_cells
  mean_logexpr <- Matrix::rowMeans(lgc)
  
  keep_genes <- det_frac >= opt$min_frac_cells & mean_logexpr > opt$min_mean_logexpr
  
  cat(sprintf("Keeping %d / %d genes (>=%.3f frac cells, mean logexpr > %.3f)\n",
              sum(keep_genes), length(keep_genes),
              opt$min_frac_cells, opt$min_mean_logexpr))
  
  sce <- sce[keep_genes, , drop = FALSE]
} else {
  cat("Global gene filter disabled (min_frac_cells <= 0 and min_mean_logexpr <= 0).\n")
}

cat(sprintf("After filtering (if any): %d genes x %d cells\n",
            nrow(sce), ncol(sce)))

# ===================== HVGs and PCA =====================
cat("Computing HVGs...\n")
set.seed(42)
dec <- scran::modelGeneVar(sce, assay.type = "logcounts")
hvg <- scran::getTopHVGs(dec, n = opt$n_hvg)
cat("  Selected ", length(hvg), " HVGs.\n", sep = "")

max_pcs <- min(
  opt$n_pcs,
  ncol(sce) - 1L,
  length(hvg) - 1L
)
if (max_pcs < 2L) {
  stop("Not enough rank for PCA (max_pcs < 2).")
}
cat("Running PCA on HVGs (", max_pcs, " PCs)...\n", sep = "")

set.seed(42)
sce <- scater::runPCA(
  sce,
  exprs_values = "logcounts",
  subset_row   = hvg,
  ncomponents  = max_pcs,
  BSPARAM      = BiocSingular::IrlbaParam(deferred = FALSE)
)
cat("PCA completed.\n")

# ===================== UMAP =====================
umap_neighbors <- min(opt$umap_neighbors, ncol(sce) - 1L)
if (umap_neighbors < 2L) {
  stop("Too few cells for UMAP (effective n_neighbors < 2).")
}

cat("Running UMAP (n_neighbors = ", umap_neighbors,
    ", min_dist = ", opt$umap_min_dist,
    ", name = '", opt$umap_name, "')...\n", sep = "")

set.seed(42)
sce <- scater::runUMAP(
  sce,
  dimred      = "PCA",
  name        = opt$umap_name,
  n_neighbors = umap_neighbors,
  min_dist    = opt$umap_min_dist
)
cat("UMAP completed.\n")

# ===================== SNN graph + clustering =====================
cat("Building SNN graph (k = ", opt$k_snn, ")...\n", sep = "")
k_effective <- min(opt$k_snn, ncol(sce) - 1L)
if (k_effective < 2L) {
  stop("Too few cells for meaningful SNN clustering (effective k < 2).")
}

g <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = k_effective)

method <- tolower(opt$cluster_method)
cluster_col <- paste0("cluster_", method)

cat("Clustering with method: ", method, "\n", sep = "")
if (method == "leiden" && "cluster_leiden" %in% ls("package:igraph")) {
  cl <- igraph::cluster_leiden(g)$membership
} else if (method == "leiden") {
  warning("cluster_leiden not available in this igraph; falling back to louvain.")
  cl <- igraph::cluster_louvain(g)$membership
  cluster_col <- "cluster_louvain"
} else {
  # default: louvain
  cl <- igraph::cluster_louvain(g)$membership
  cluster_col <- "cluster_louvain"
}

colData(sce)[[cluster_col]] <- factor(cl)
cat("Clustering completed: ", length(unique(cl)),
    " clusters (stored in colData(sce)$", cluster_col, ").\n", sep = "")

# cluster size summary
clus_tab <- table(colData(sce)[[cluster_col]])
clus_dt <- data.table(
  cluster = names(clus_tab),
  n_cells = as.integer(clus_tab)
)[order(-n_cells)]
print(clus_dt)

fwrite(
  clus_dt,
  file = file.path(opt$outdir, "global_cluster_sizes.tsv"),
  sep = "\t"
)

# ===================== Save updated SCE =====================
out_file <- file.path(opt$outdir, "sce_global_umap_clusters.rds")
saveRDS(sce, out_file)
cat("Saved updated SCE (with PCA, ", opt$umap_name, ", and ", cluster_col,
    ") to: ", out_file, "\n", sep = "")
cat("Cluster size summary written to global_cluster_sizes.tsv\n")
