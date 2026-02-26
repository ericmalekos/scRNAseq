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

# ========================================================
# GLOBAL ANALYSIS (ALL AGES COMBINED)
# ========================================================
cat("\n========================================\n")
cat("GLOBAL ANALYSIS (ALL AGES COMBINED)\n")
cat("========================================\n\n")

sce_global <- sce  # work on a copy

# ===================== Global gene filter =====================
if (opt$min_frac_cells > 0 || opt$min_mean_logexpr > 0) {
  cat("Applying global gene filter...\n")
  cnt <- assay(sce_global, "counts")
  lgc <- assay(sce_global, "logcounts")
  
  n_cells <- ncol(sce_global)
  det_frac <- Matrix::rowSums(cnt > 0) / n_cells
  mean_logexpr <- Matrix::rowMeans(lgc)
  
  keep_genes <- det_frac >= opt$min_frac_cells & mean_logexpr > opt$min_mean_logexpr
  
  cat(sprintf("Keeping %d / %d genes (>=%.3f frac cells, mean logexpr > %.3f)\n",
              sum(keep_genes), length(keep_genes),
              opt$min_frac_cells, opt$min_mean_logexpr))
  
  sce_global <- sce_global[keep_genes, , drop = FALSE]
} else {
  cat("Global gene filter disabled (min_frac_cells <= 0 and min_mean_logexpr <= 0).\n")
}

cat(sprintf("After filtering (if any): %d genes x %d cells\n",
            nrow(sce_global), ncol(sce_global)))

# ===================== HVGs and PCA =====================
cat("Computing HVGs...\n")
set.seed(42)
dec_global <- scran::modelGeneVar(sce_global, assay.type = "logcounts")
hvg_global <- scran::getTopHVGs(dec_global, n = opt$n_hvg)
cat("  Selected ", length(hvg_global), " HVGs.\n", sep = "")

max_pcs_global <- min(
  opt$n_pcs,
  ncol(sce_global) - 1L,
  length(hvg_global) - 1L
)
if (max_pcs_global < 2L) {
  stop("Not enough rank for PCA (max_pcs < 2).")
}
cat("Running PCA on HVGs (", max_pcs_global, " PCs)...\n", sep = "")

set.seed(42)
sce_global <- scater::runPCA(
  sce_global,
  exprs_values = "logcounts",
  subset_row   = hvg_global,
  ncomponents  = max_pcs_global,
  BSPARAM      = BiocSingular::IrlbaParam(deferred = FALSE)
)
cat("PCA completed.\n")

# ===================== UMAP =====================
umap_neighbors_global <- min(opt$umap_neighbors, ncol(sce_global) - 1L)
if (umap_neighbors_global < 2L) {
  stop("Too few cells for UMAP (effective n_neighbors < 2).")
}

cat("Running UMAP (n_neighbors = ", umap_neighbors_global,
    ", min_dist = ", opt$umap_min_dist,
    ", name = '", opt$umap_name, "')...\n", sep = "")

set.seed(42)
sce_global <- scater::runUMAP(
  sce_global,
  dimred      = "PCA",
  name        = opt$umap_name,
  n_neighbors = umap_neighbors_global,
  min_dist    = opt$umap_min_dist
)
cat("UMAP completed.\n")

# ===================== SNN graph + clustering =====================
cat("Building SNN graph (k = ", opt$k_snn, ")...\n", sep = "")
k_effective_global <- min(opt$k_snn, ncol(sce_global) - 1L)
if (k_effective_global < 2L) {
  stop("Too few cells for meaningful SNN clustering (effective k < 2).")
}

g_global <- scran::buildSNNGraph(sce_global, use.dimred = "PCA", k = k_effective_global)

method <- tolower(opt$cluster_method)
cluster_col <- paste0("cluster_", method)

cat("Clustering with method: ", method, "\n", sep = "")
if (method == "leiden" && "cluster_leiden" %in% ls("package:igraph")) {
  cl_global <- igraph::cluster_leiden(g_global)$membership
} else if (method == "leiden") {
  warning("cluster_leiden not available in this igraph; falling back to louvain.")
  cl_global <- igraph::cluster_louvain(g_global)$membership
  cluster_col <- "cluster_louvain"
} else {
  # default: louvain
  cl_global <- igraph::cluster_louvain(g_global)$membership
  cluster_col <- "cluster_louvain"
}

colData(sce_global)[[cluster_col]] <- factor(cl_global)
cat("Clustering completed: ", length(unique(cl_global)),
    " clusters (stored in colData(sce)$", cluster_col, ").\n", sep = "")

# cluster size summary
clus_tab_global <- table(colData(sce_global)[[cluster_col]])
clus_dt_global <- data.table(
  cluster = names(clus_tab_global),
  n_cells = as.integer(clus_tab_global)
)[order(-n_cells)]
print(clus_dt_global)

fwrite(
  clus_dt_global,
  file = file.path(opt$outdir, "global_cluster_sizes.tsv"),
  sep = "\t"
)

# ===================== Save global SCE =====================
out_file_global <- file.path(opt$outdir, "sce_global_all_ages_umap_clusters.rds")
saveRDS(sce_global, out_file_global)
cat("Saved global (all ages) SCE to: ", out_file_global, "\n", sep = "")
cat("Cluster size summary written to global_cluster_sizes.tsv\n")

# ========================================================
# AGE-STRATIFIED ANALYSIS
# ========================================================
if ("age" %in% colnames(colData(sce))) {
  cat("\n========================================\n")
  cat("AGE-STRATIFIED ANALYSIS\n")
  cat("========================================\n\n")
  
  age_groups <- sort(unique(as.character(sce$age)))
  age_groups <- age_groups[!is.na(age_groups) & age_groups != "nan" & age_groups != ""]
  
  cat("Found ", length(age_groups), " age groups: ", paste(age_groups, collapse = ", "), "\n\n", sep = "")
  
  # Initialize summary table
  summary_dt <- data.table(
    analysis = character(),
    n_cells = integer(),
    n_genes = integer(),
    n_hvgs = integer(),
    n_pcs = integer(),
    n_clusters = integer()
  )
  
  # Add global summary
  summary_dt <- rbind(summary_dt, data.table(
    analysis = "global_all_ages",
    n_cells = ncol(sce_global),
    n_genes = nrow(sce_global),
    n_hvgs = length(hvg_global),
    n_pcs = max_pcs_global,
    n_clusters = length(unique(cl_global))
  ))
  
  for (age in age_groups) {
    cat("--- Processing age: ", age, " ---\n", sep = "")
    
    # Subset to this age
    age_mask <- !is.na(sce$age) & sce$age == age
    sce_age <- sce[, age_mask, drop = FALSE]
    cat(sprintf("  Age subset: %d genes x %d cells\n", nrow(sce_age), ncol(sce_age)))
    
    # Skip if too few cells
    if (ncol(sce_age) < 50) {
      cat("  Fewer than 50 cells; skipping analysis for this age.\n\n")
      next
    }
    
    # Gene filter for this age
    if (opt$min_frac_cells > 0 || opt$min_mean_logexpr > 0) {
      cat("  Applying gene filter...\n")
      cnt_age <- assay(sce_age, "counts")
      lgc_age <- assay(sce_age, "logcounts")
      
      n_cells_age <- ncol(sce_age)
      det_frac_age <- Matrix::rowSums(cnt_age > 0) / n_cells_age
      mean_logexpr_age <- Matrix::rowMeans(lgc_age)
      
      keep_genes_age <- det_frac_age >= opt$min_frac_cells & mean_logexpr_age > opt$min_mean_logexpr
      
      cat(sprintf("  Keeping %d / %d genes\n", sum(keep_genes_age), length(keep_genes_age)))
      sce_age <- sce_age[keep_genes_age, , drop = FALSE]
    }
    
    # HVG selection for this age
    cat("  Computing HVGs...\n")
    set.seed(42)
    dec_age <- scran::modelGeneVar(sce_age, assay.type = "logcounts")
    hvg_age <- scran::getTopHVGs(dec_age, n = opt$n_hvg)
    cat("    Selected ", length(hvg_age), " HVGs.\n", sep = "")
    
    # PCA for this age
    max_pcs_age <- min(
      opt$n_pcs,
      ncol(sce_age) - 1L,
      length(hvg_age) - 1L
    )
    if (max_pcs_age < 2L) {
      cat("  Not enough rank for PCA (max_pcs < 2); skipping this age.\n\n")
      next
    }
    
    cat("  Running PCA (", max_pcs_age, " PCs)...\n", sep = "")
    set.seed(42)
    sce_age <- scater::runPCA(
      sce_age,
      exprs_values = "logcounts",
      subset_row   = hvg_age,
      ncomponents  = max_pcs_age,
      BSPARAM      = BiocSingular::IrlbaParam(deferred = FALSE)
    )
    cat("  PCA completed.\n")
    
    # UMAP for this age
    umap_neighbors_age <- min(opt$umap_neighbors, ncol(sce_age) - 1L)
    if (umap_neighbors_age < 2L) {
      cat("  Too few cells for UMAP (effective n_neighbors < 2); skipping this age.\n\n")
      next
    }
    
    cat("  Running UMAP (n_neighbors = ", umap_neighbors_age, ")...\n", sep = "")
    set.seed(42)
    sce_age <- scater::runUMAP(
      sce_age,
      dimred      = "PCA",
      name        = opt$umap_name,
      n_neighbors = umap_neighbors_age,
      min_dist    = opt$umap_min_dist
    )
    cat("  UMAP completed.\n")
    
    # Clustering for this age
    cat("  Building SNN graph and clustering...\n")
    k_effective_age <- min(opt$k_snn, ncol(sce_age) - 1L)
    if (k_effective_age < 2L) {
      cat("  Too few cells for SNN clustering (effective k < 2); skipping clustering.\n\n")
      next
    }
    
    g_age <- scran::buildSNNGraph(sce_age, use.dimred = "PCA", k = k_effective_age)
    
    if (method == "leiden" && "cluster_leiden" %in% ls("package:igraph")) {
      cl_age <- igraph::cluster_leiden(g_age)$membership
    } else if (method == "leiden") {
      cl_age <- igraph::cluster_louvain(g_age)$membership
    } else {
      cl_age <- igraph::cluster_louvain(g_age)$membership
    }
    
    colData(sce_age)[[cluster_col]] <- factor(cl_age)
    cat("  Clustering completed: ", length(unique(cl_age)), " clusters.\n", sep = "")
    
    # Cluster size summary for this age
    clus_tab_age <- table(colData(sce_age)[[cluster_col]])
    clus_dt_age <- data.table(
      cluster = names(clus_tab_age),
      n_cells = as.integer(clus_tab_age)
    )[order(-n_cells)]
    
    # Save age-specific outputs
    age_sanitized <- gsub("[^A-Za-z0-9]+", "_", age)
    age_sanitized <- gsub("^_+|_+$", "", age_sanitized)
    
    out_file_age <- file.path(opt$outdir, paste0("sce_global_", age_sanitized, "_umap_clusters.rds"))
    saveRDS(sce_age, out_file_age)
    cat("  Saved age-specific SCE to: ", out_file_age, "\n", sep = "")
    
    clus_file_age <- file.path(opt$outdir, paste0("global_", age_sanitized, "_cluster_sizes.tsv"))
    fwrite(clus_dt_age, clus_file_age, sep = "\t")
    cat("  Cluster sizes written to: ", clus_file_age, "\n\n", sep = "")
    
    # Add to summary
    summary_dt <- rbind(summary_dt, data.table(
      analysis = age,
      n_cells = ncol(sce_age),
      n_genes = nrow(sce_age),
      n_hvgs = length(hvg_age),
      n_pcs = max_pcs_age,
      n_clusters = length(unique(cl_age))
    ))
  }
  
  # Save summary table
  summary_file <- file.path(opt$outdir, "global_analysis_summary.tsv")
  fwrite(summary_dt, summary_file, sep = "\t")
  cat("\n========================================\n")
  cat("Summary of all analyses written to: ", summary_file, "\n", sep = "")
  
} else {
  warning("'age' column not found in colData; skipping age-stratified analysis.")
}

cat("\n========================================\n")
cat("All done!\n")
cat("========================================\n")
cat("Outputs:\n")
cat("  - Global (all ages): sce_global_all_ages_umap_clusters.rds\n")
if ("age" %in% colnames(colData(sce))) {
  cat("  - Age-stratified:    sce_global_<age>_umap_clusters.rds (one per age)\n")
  cat("  - Summary table:     global_analysis_summary.tsv\n")
}