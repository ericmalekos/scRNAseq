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
  add_option("--outdir",  help="Output directory for results.", metavar="dir") |>
  add_option("--cell_classes", help="Comma-separated list of cell_ontology_class values to analyze together.", metavar="string") |>
  add_option("--list_cells", action="store_true", default=FALSE,
             help="Print unique cell classes with counts and exit (no analysis).") |>
  add_option("--min_frac_cells", type="double", default=0.00000001,
             help="Global gene filter: keep genes expressed in >= this fraction of cells [default %default].") |>
  add_option("--min_mean_logexpr", type="double", default=0,
             help="Global gene filter: keep genes with mean(logcounts) > this value [default %default].") |>
  
  # HVG / PCA / UMAP / clustering params
  add_option("--hvg_bio_threshold", type="double", default=0.5,
             help="HVG rule: keep genes with biological variance (bio) > this value [default %default].") |>
  add_option("--min_hvg", type="integer", default=1000,
             help="Minimum number of HVGs (fallback if threshold yields fewer) [default %default].") |>
  add_option("--max_hvg", type="integer", default=3000,
             help="Maximum number of HVGs (cap for getTopHVGs) [default %default].") |>
  add_option("--k_snn", type="integer", default=20,
             help="k for SNN graph construction [default %default].") |>
  add_option("--n_pcs", type="integer", default=50,
             help="Maximum number of PCs to compute with runPCA [default %default].") |>
  add_option("--umap_neighbors", type="integer", default=15,
             help="n_neighbors for UMAP [default %default].") |>
  add_option("--umap_min_dist", type="double", default=0.3,
             help="min_dist for UMAP [default %default].") |>
  parse_args()

if (is.null(opt$sce_rds)) {
  stop("Required: --sce_rds")
}

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

# ---------- Handle --list_cells short circuit ----------
if (opt$list_cells) {
  cat("\n=== Cell Class Summary ===\n")
  class_counts <- table(sce$cell_ontology_class)
  class_summary <- data.table(
    cell_ontology_class = names(class_counts),
    n_cells = as.integer(class_counts)
  )
  # Sort by count descending
  class_summary <- class_summary[order(-n_cells)]
  
  # Print formatted output
  cat(sprintf("\nTotal: %d unique cell classes\n\n", nrow(class_summary)))
  cat(sprintf("%-50s %10s\n", "Cell Class", "Count"))
  cat(sprintf("%-50s %10s\n", paste(rep("-", 50), collapse=""), paste(rep("-", 10), collapse="")))
  for (i in 1:nrow(class_summary)) {
    cat(sprintf("%-50s %10d\n", 
                class_summary$cell_ontology_class[i], 
                class_summary$n_cells[i]))
  }
  cat("\n")
  quit(save = "no", status = 0)
}

# ---------- Check required parameters for analysis ----------
if (is.null(opt$outdir) || is.null(opt$cell_classes)) {
  stop("Required for analysis: --outdir and --cell_classes")
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Parse cell class list ----------
requested_classes <- strsplit(opt$cell_classes, ",")[[1]]
requested_classes <- trimws(requested_classes)
cat("\nRequested cell classes:\n")
cat("  - ", paste(requested_classes, collapse = "\n  - "), "\n", sep = "")

# Check which requested classes exist in data
available_classes <- unique(as.character(sce$cell_ontology_class))
missing_classes <- setdiff(requested_classes, available_classes)
if (length(missing_classes) > 0) {
  warning("The following requested classes were not found in the data:\n  - ",
          paste(missing_classes, collapse = "\n  - "))
}

valid_classes <- intersect(requested_classes, available_classes)
if (length(valid_classes) == 0) {
  stop("None of the requested cell classes were found in the data.")
}

cat("\nProcessing ", length(valid_classes), " valid cell classes:\n", sep = "")
cat("  - ", paste(valid_classes, collapse = "\n  - "), "\n", sep = "")

# ---------- Subset to requested cell classes ----------
cat("\nSubsetting to requested cell classes...\n")
sce_subset <- sce[, sce$cell_ontology_class %in% valid_classes, drop = FALSE]
cat(sprintf("Subset contains: %d genes x %d cells\n", nrow(sce_subset), ncol(sce_subset)))

# Print breakdown by class
class_counts <- table(sce_subset$cell_ontology_class)
cat("\nCells per class in subset:\n")
for (cls in valid_classes) {
  cat(sprintf("  %-40s: %6d cells\n", cls, class_counts[cls]))
}

# ---------- Global gene filter on subset ----------
cat("\nApplying gene filter to subset...\n")
cnt <- assay(sce_subset, "counts")
lgc <- assay(sce_subset, "logcounts")

n_cells <- ncol(sce_subset)
det_frac <- Matrix::rowSums(cnt > 0) / n_cells
mean_logexpr <- Matrix::rowMeans(lgc)

keep_genes <- det_frac >= opt$min_frac_cells & mean_logexpr > opt$min_mean_logexpr

cat(sprintf("Keeping %d / %d genes after filter (>=%.3f frac cells, mean logexpr > %.3f)\n",
            sum(keep_genes), length(keep_genes),
            opt$min_frac_cells, opt$min_mean_logexpr))

sce_subset <- sce_subset[keep_genes, , drop = FALSE]

# Save filtered subset before PCA/UMAP
saveRDS(sce_subset, file.path(opt$outdir, "sce_subset_filtered.rds"))
cat("Saved filtered subset SCE to sce_subset_filtered.rds\n")

# ---------- Check minimum cells ----------
if (ncol(sce_subset) < 10) {
  stop("Fewer than 10 cells in subset; cannot proceed with PCA/UMAP/clustering.")
}

# ---------- HVG selection ----------
cat("\nFitting mean-variance trend and selecting HVGs...\n")
set.seed(42)
dec <- scran::modelGeneVar(sce_subset, assay.type = "logcounts")
hvg <- rownames(dec)[dec$bio > opt$hvg_bio_threshold]
cat(sprintf("HVGs with bio > %.3f: %d genes\n", opt$hvg_bio_threshold, length(hvg)))

if (length(hvg) < opt$min_hvg) {
  n_top <- min(opt$max_hvg, nrow(dec))
  cat(sprintf("Fewer than %d HVGs; taking top %d HVGs by bio.\n", opt$min_hvg, n_top))
  hvg <- scran::getTopHVGs(dec, n = n_top)
}

cat(sprintf("Using %d HVGs for PCA/UMAP\n", length(hvg)))
rowData(sce_subset)$is_hvg <- rownames(sce_subset) %in% hvg

# ---------- PCA ----------
cat("\nRunning PCA on HVGs...\n")
max_rank <- min(
  opt$n_pcs,
  ncol(sce_subset) - 1L,
  length(hvg) - 1L
)
if (max_rank < 2L) {
  stop("Not enough rank for PCA (max_rank < 2); cannot proceed.")
}

set.seed(42)
sce_subset <- scater::runPCA(
  sce_subset,
  exprs_values = "logcounts",
  subset_row   = hvg,
  ncomponents  = max_rank,
  BSPARAM      = BiocSingular::IrlbaParam(deferred = FALSE)
)
if (!"PCA" %in% names(reducedDims(sce_subset))) {
  stop("runPCA did not produce a 'PCA' reducedDim.")
}
pca_mat <- reducedDim(sce_subset, "PCA")
cat(sprintf("PCA completed: %d PCs\n", ncol(pca_mat)))

# ---------- UMAP ----------
umap_neighbors <- min(opt$umap_neighbors, ncol(sce_subset) - 1L)
if (umap_neighbors < 2L) {
  warning("Too few cells for UMAP (effective n_neighbors < 2); skipping UMAP and clustering.")
} else {
  cat("\nRunning UMAP on PCA embedding (n_neighbors = ", umap_neighbors,
      ", min_dist = ", opt$umap_min_dist, ")...\n", sep = "")
  set.seed(42)
  sce_subset <- scater::runUMAP(
    sce_subset,
    dimred      = "PCA",
    name        = "UMAP",
    n_neighbors = umap_neighbors,
    min_dist    = opt$umap_min_dist
  )
  cat("UMAP completed.\n")
  
  # ---------- SNN graph + Louvain clustering ----------
  cat("\nBuilding SNN graph (k = ", opt$k_snn, ") and clustering (Louvain)...\n", sep = "")
  k_effective <- min(opt$k_snn, ncol(sce_subset) - 1L)
  if (k_effective < 2L) {
    cat("Too few cells for meaningful SNN clustering (effective k < 2); skipping clustering.\n")
  } else {
    g <- scran::buildSNNGraph(sce_subset, use.dimred = "PCA", k = k_effective)
    cl <- igraph::cluster_louvain(g)$membership
    colData(sce_subset)$cluster_louvain <- factor(cl)
    cat("Clustering completed: ", length(unique(cl)),
        " clusters identified.\n", sep = "")
  }
}

# ---------- Save final output (all ages combined) ----------
out_file <- file.path(opt$outdir, "sce_subset_all_ages_pca_umap.rds")
saveRDS(sce_subset, out_file)
cat("\nSaved combined (all ages) subset SCE with PCA/UMAP to: ", out_file, "\n", sep = "")

# ---------- Save summary for combined ----------
summary_data <- data.table(
  analysis = "all_ages",
  cell_ontology_class = paste(valid_classes, collapse = ", "),
  n_cells = ncol(sce_subset),
  total_genes = nrow(sce_subset),
  n_hvgs = length(hvg),
  n_pcs = ncol(pca_mat)
)

# ---------- Age-stratified analysis ----------
if ("age" %in% colnames(colData(sce_subset))) {
  cat("\n========================================\n")
  cat("Starting age-stratified analysis...\n")
  cat("========================================\n")
  
  age_groups <- sort(unique(as.character(sce_subset$age)))
  age_groups <- age_groups[!is.na(age_groups) & age_groups != "nan" & age_groups != ""]
  
  cat("\nFound ", length(age_groups), " age groups: ", paste(age_groups, collapse = ", "), "\n", sep = "")
  
  for (age in age_groups) {
    cat("\n--- Processing age: ", age, " ---\n", sep = "")
    
    # Subset to this age
    sce_age <- sce_subset[, sce_subset$age == age, drop = FALSE]
    cat(sprintf("  Age subset: %d genes x %d cells\n", nrow(sce_age), ncol(sce_age)))
    
    # Skip if too few cells
    if (ncol(sce_age) < 10) {
      cat("  Fewer than 10 cells; skipping PCA/UMAP for this age.\n")
      next
    }
    
    # Drop genes never expressed in this age subset
    cnt_age <- assay(sce_age, "counts")
    keep_genes_age <- Matrix::rowSums(cnt_age > 0) > 0
    sce_age <- sce_age[keep_genes_age, , drop = FALSE]
    cat(sprintf("  After dropping never-expressed genes: %d genes x %d cells\n",
                nrow(sce_age), ncol(sce_age)))
    
    # HVG selection for this age
    cat("  Fitting mean-variance trend and selecting HVGs...\n")
    set.seed(42)
    dec_age <- scran::modelGeneVar(sce_age, assay.type = "logcounts")
    hvg_age <- rownames(dec_age)[dec_age$bio > opt$hvg_bio_threshold]
    cat(sprintf("  HVGs with bio > %.3f: %d genes\n", opt$hvg_bio_threshold, length(hvg_age)))
    
    if (length(hvg_age) < opt$min_hvg) {
      n_top <- min(opt$max_hvg, nrow(dec_age))
      cat(sprintf("  Fewer than %d HVGs; taking top %d HVGs by bio.\n", opt$min_hvg, n_top))
      hvg_age <- scran::getTopHVGs(dec_age, n = n_top)
    }
    
    cat(sprintf("  Using %d HVGs for PCA/UMAP\n", length(hvg_age)))
    rowData(sce_age)$is_hvg <- rownames(sce_age) %in% hvg_age
    
    # PCA for this age
    cat("  Running PCA...\n")
    max_rank_age <- min(
      opt$n_pcs,
      ncol(sce_age) - 1L,
      length(hvg_age) - 1L
    )
    if (max_rank_age < 2L) {
      cat("  Not enough rank for PCA (max_rank < 2); skipping this age.\n")
      next
    }
    
    set.seed(42)
    sce_age <- scater::runPCA(
      sce_age,
      exprs_values = "logcounts",
      subset_row   = hvg_age,
      ncomponents  = max_rank_age,
      BSPARAM      = BiocSingular::IrlbaParam(deferred = FALSE)
    )
    pca_mat_age <- reducedDim(sce_age, "PCA")
    cat(sprintf("  PCA completed: %d PCs\n", ncol(pca_mat_age)))
    
    # UMAP for this age
    umap_neighbors_age <- min(opt$umap_neighbors, ncol(sce_age) - 1L)
    if (umap_neighbors_age < 2L) {
      cat("  Too few cells for UMAP (effective n_neighbors < 2); skipping UMAP/clustering.\n")
    } else {
      cat("  Running UMAP (n_neighbors = ", umap_neighbors_age, ")...\n", sep = "")
      set.seed(42)
      sce_age <- scater::runUMAP(
        sce_age,
        dimred      = "PCA",
        name        = "UMAP",
        n_neighbors = umap_neighbors_age,
        min_dist    = opt$umap_min_dist
      )
      cat("  UMAP completed.\n")
      
      # Clustering for this age
      cat("  Building SNN graph and clustering...\n")
      k_effective_age <- min(opt$k_snn, ncol(sce_age) - 1L)
      if (k_effective_age < 2L) {
        cat("  Too few cells for SNN clustering (effective k < 2); skipping clustering.\n")
      } else {
        g_age <- scran::buildSNNGraph(sce_age, use.dimred = "PCA", k = k_effective_age)
        cl_age <- igraph::cluster_louvain(g_age)$membership
        colData(sce_age)$cluster_louvain <- factor(cl_age)
        cat("  Clustering completed: ", length(unique(cl_age)), " clusters.\n", sep = "")
      }
    }
    
    # Save age-specific output
    age_sanitized <- gsub("[^A-Za-z0-9]+", "_", age)
    age_sanitized <- gsub("^_+|_+$", "", age_sanitized)
    out_file_age <- file.path(opt$outdir, paste0("sce_subset_", age_sanitized, "_pca_umap.rds"))
    saveRDS(sce_age, out_file_age)
    cat("  Saved age-specific SCE to: ", out_file_age, "\n", sep = "")
    
    # Add to summary
    summary_data <- rbind(
      summary_data,
      data.table(
        analysis = age,
        cell_ontology_class = paste(valid_classes, collapse = ", "),
        n_cells = ncol(sce_age),
        total_genes = nrow(sce_age),
        n_hvgs = length(hvg_age),
        n_pcs = ncol(pca_mat_age)
      )
    )
  }
  
  cat("\nAge-stratified analysis completed.\n")
} else {
  warning("'age' column not found in colData; skipping age-stratified analysis.")
}

# ---------- Save summary ----------
fwrite(summary_data,
       file.path(opt$outdir, "subset_summary.tsv"),
       sep = "\t")
cat("\nSummary written to subset_summary.tsv\n")

cat("\n========================================\n")
cat("All done!\n")
cat("========================================\n")
cat("Outputs:\n")
cat("  - Filtered subset:       sce_subset_filtered.rds\n")
cat("  - All ages combined:     sce_subset_all_ages_pca_umap.rds\n")
if ("age" %in% colnames(colData(sce_subset))) {
  cat("  - Age-stratified:        sce_subset_<age>_pca_umap.rds (one per age)\n")
}
cat("  - Summary:               subset_summary.tsv\n")