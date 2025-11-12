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
  add_option("--outdir",  help="Output directory for results and per-class SCEs.", metavar="dir") |>
  add_option("--min_frac_cells", type="double", default=0.001,
             help="Global gene filter: keep genes expressed in >= this fraction of cells [default %default].") |>
  add_option("--min_mean_logexpr", type="double", default=0,
             help="Global gene filter: keep genes with mean(logcounts) > this value [default %default].") |>
  
  # Global HVG / PCA / UMAP params
  add_option("--global_n_hvg", type="integer", default=3000,
             help="Number of HVGs to use globally for PCA/UMAP [default %default].") |>
  add_option("--global_n_pcs", type="integer", default=50,
             help="Maximum number of PCs for global PCA [default %default].") |>
  add_option("--global_umap_neighbors", type="integer", default=15,
             help="n_neighbors for global UMAP [default %default].") |>
  add_option("--global_umap_min_dist", type="double", default=0.3,
             help="min_dist for global UMAP [default %default].") |>
  
  # Per-class HVG / PCA / UMAP / clustering params
  add_option("--hvg_bio_threshold", type="double", default=0.5,
             help="Per-class HVG rule: keep genes with biological variance (bio) > this value [default %default].") |>
  add_option("--min_hvg", type="integer", default=1000,
             help="Minimum number of HVGs per class (fallback if threshold yields fewer) [default %default].") |>
  add_option("--max_hvg", type="integer", default=3000,
             help="Maximum number of HVGs per class (cap for getTopHVGs) [default %default].") |>
  add_option("--k_snn", type="integer", default=20,
             help="k for SNN graph construction [default %default].") |>
  add_option("--n_pcs", type="integer", default=50,
             help="Maximum number of PCs to compute per class with runPCA [default %default].") |>
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
if (!"cell_ontology_class" %in% colnames(colData(sce))) {
  stop("colData(sce) must contain 'cell_ontology_class'.")
}

cat(sprintf("Initial SCE dimensions: %d genes x %d cells\n",
            nrow(sce), ncol(sce)))

# ---------- Drop NA/'nan' classes ----------
n_before <- ncol(sce)
sce <- sce[, !is.na(sce$cell_ontology_class) & sce$cell_ontology_class != "nan", drop = FALSE]
n_after <- ncol(sce)
cat(sprintf("Dropped %d cells with NA/'nan' cell_ontology_class. Remaining: %d cells.\n",
            n_before - n_after, n_after))

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

sce <- sce[keep_genes, , drop = FALSE]  # globally filtered SCE

# Save globally filtered SCE (no PCA/UMAP yet)
saveRDS(sce, file.path(opt$outdir, "sce_filtered_global.rds"))
cat("Saved globally filtered SCE to sce_filtered_global.rds\n")

# ---------- Global HVG / PCA / UMAP ----------
cat("Computing global HVGs, PCA, and UMAP...\n")
set.seed(42)
dec_global <- scran::modelGeneVar(sce, assay.type = "logcounts")
hvg_global <- scran::getTopHVGs(dec_global, n = opt$global_n_hvg)
cat("  Selected ", length(hvg_global), " global HVGs.\n", sep = "")

max_pcs_global <- min(
  opt$global_n_pcs,
  ncol(sce) - 1L,
  length(hvg_global) - 1L
)
if (max_pcs_global < 2L) {
  warning("Not enough rank for global PCA (max_pcs_global < 2); skipping global PCA/UMAP.")
} else {
  cat("  Running global PCA (", max_pcs_global, " PCs)...\n", sep = "")
  set.seed(42)
  sce <- scater::runPCA(
    sce,
    exprs_values = "logcounts",
    subset_row   = hvg_global,
    ncomponents  = max_pcs_global,
    BSPARAM      = BiocSingular::IrlbaParam(deferred = FALSE)
  )
  cat("  Global PCA completed.\n")
  
  # Global UMAP
  umap_neighbors_global <- min(opt$global_umap_neighbors, ncol(sce) - 1L)
  if (umap_neighbors_global < 2L) {
    warning("Too few cells for global UMAP (effective n_neighbors < 2); skipping global UMAP.")
  } else {
    cat("  Running global UMAP (n_neighbors = ", umap_neighbors_global,
        ", min_dist = ", opt$global_umap_min_dist, ")...\n", sep = "")
    set.seed(42)
    sce <- scater::runUMAP(
      sce,
      dimred      = "PCA",
      name        = "UMAP",
      n_neighbors = umap_neighbors_global,
      min_dist    = opt$global_umap_min_dist
    )
    cat("  Global UMAP completed.\n")
  }
  
  # Save global SCE with PCA+UMAP
  saveRDS(sce, file.path(opt$outdir, "sce_filtered_global_pca_umap.rds"))
  cat("Saved global SCE with PCA+UMAP to sce_filtered_global_pca_umap.rds\n")
}

# ---------- Prepare per-class analysis ----------
cell_classes <- sort(unique(as.character(sce$cell_ontology_class)))
cell_classes <- cell_classes[!is.na(cell_classes) & cell_classes != "nan"]

cat("Found ", length(cell_classes), " distinct cell_ontology_class values after filtering.\n", sep = "")
cat("First few classes:\n")
cat("  - ", paste(head(cell_classes, 10), collapse = "\n  - "), "\n", sep = "")

summary_tab <- data.table(
  cell_ontology_class = character(),
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

# ---------- Loop over cell_ontology_class ----------
for (cls in cell_classes) {
  cat("\n=== Processing cell_ontology_class: ", cls, " ===\n", sep = "")
  sce_sub <- sce[, sce$cell_ontology_class == cls, drop = FALSE]
  cat(sprintf("  Subset: %d genes x %d cells\n", nrow(sce_sub), ncol(sce_sub)))
  
  # skip tiny subsets outright
  if (ncol(sce_sub) < 10) {
    cat("  Fewer than 10 cells; skipping PCA/UMAP/clustering for this class.\n")
    next
  }
  
  # drop genes never expressed in this class
  cnt_sub <- assay(sce_sub, "counts")
  keep_genes_cls <- Matrix::rowSums(cnt_sub > 0) > 0
  sce_sub <- sce_sub[keep_genes_cls, , drop = FALSE]
  cat(sprintf("  After dropping never-expressed genes: %d genes x %d cells\n",
              nrow(sce_sub), ncol(sce_sub)))
  
  # HVG selection with per-class modelGeneVar
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
  
  # ---------- PCA per class ----------
  cat("  Running PCA on HVGs (runPCA, ExactParam)...\n")
  max_rank <- min(
    opt$n_pcs,
    ncol(sce_sub) - 1L,
    length(hvg) - 1L
  )
  if (max_rank < 2L) {
    cat("  Not enough rank (max_rank < 2); skipping PCA/UMAP/clustering for this class.\n")
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
  
  # ---------- UMAP per class ----------
  umap_neighbors <- min(15L, ncol(sce_sub) - 1L)
  if (umap_neighbors < 2L) {
    cat("  Too few cells for UMAP (effective n_neighbors < 2); skipping UMAP and clustering.\n")
    next
  }
  
  cat("  Running UMAP on PCA embedding (n_neighbors = ", umap_neighbors, ")...\n", sep = "")
  set.seed(42)
  sce_sub <- scater::runUMAP(
    sce_sub,
    dimred      = "PCA",
    name        = "UMAP",
    n_neighbors = umap_neighbors
  )
  
  # ---------- SNN graph + Louvain clustering per class ----------
  cat("  Building SNN graph (k = ", opt$k_snn, ") and clustering (Louvain)...\n", sep = "")
  k_effective <- min(opt$k_snn, ncol(sce_sub) - 1L)
  if (k_effective < 2L) {
    cat("  Too few cells for meaningful SNN clustering (effective k < 2); skipping clustering.\n")
  } else {
    g <- scran::buildSNNGraph(sce_sub, use.dimred = "PCA", k = k_effective)
    cl <- igraph::cluster_louvain(g)$membership
    colData(sce_sub)$cluster_louvain <- factor(cl)
    cat("  Clustering completed: ", length(unique(cl)),
        " clusters for this cell_ontology_class.\n", sep = "")
  }
  
  summary_tab <- rbind(
    summary_tab,
    data.table(
      cell_ontology_class = cls,
      n_cells = ncol(sce_sub),
      n_genes = nrow(sce_sub),
      n_hvgs  = length(hvg)
    )
  )
  
  cls_sanitized <- sanitize_name(cls)
  out_file <- file.path(opt$outdir, paste0("sce_cellclass_", cls_sanitized, ".rds"))
  saveRDS(sce_sub, out_file)
  cat("  Saved per-class SCE to: ", out_file, "\n", sep = "")
}

# ---------- Save summary ----------
if (nrow(summary_tab) > 0) {
  fwrite(summary_tab,
         file.path(opt$outdir, "cellclass_summary.tsv"),
         sep = "\t")
  cat("\nSummary across cell_ontology_class subsets written to cellclass_summary.tsv\n")
}

cat("\nAll done. Global filtered SCE saved as 'sce_filtered_global.rds', ",
    "global PCA+UMAP SCE saved as 'sce_filtered_global_pca_umap.rds', ",
    "and per-class SCEs saved with prefix 'sce_cellclass_'.\n", sep = "")
