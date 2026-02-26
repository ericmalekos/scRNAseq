# ---------- Plotting utility functions ----------
#
# Source this file to load plotting functions into your R session or notebook:
#   source("plotting_functions/plot_umap.R")
#
# Available functions:
#   plot_gene_top_violin_points_by_group() - Violin + jittered points by cell type or tissue
#   plot_gene_top_boxplot_by_group()       - Boxplot by cell type or tissue
#   add_umap_group_labels()                - Add text labels to an existing UMAP plot
#   plot_gene_by_cellgroup_in_tissue()     - Expression boxplot within a single tissue
#   plot_gene_top15_boxplot()              - Top 15 cell type x tissue boxplot
#   plot_gene_umap_size_color()            - UMAP colored/sized by gene expression
#   plot_genes_umap_batch()                - Batch plot multiple genes to files

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(ggplot2)
  library(viridis)
  library(data.table)
  library(scales)
  library(dplyr)
  library(ggrepel)
})


plot_gene_top_violin_points_by_group <- function(
    sce,
    gene_id,
    assay_name          = "logcounts",
    group_mode          = c("cell_type", "tissue"),
    cell_type_col       = "cell_ontology_class",
    tissue_col          = "tissue",
    top_n               = 5,
    min_cells_per_group = 20,
    tissue_filter       = NULL,   # INCLUDE only these tissues (comma-separated string or character vector)
    celltype_filter     = NULL,   # INCLUDE only these cell types
    tissue_exclude      = NULL,   # EXCLUDE these tissues
    celltype_exclude    = NULL,   # EXCLUDE these cell types
    top_stat            = c("median", "mean"),  # how to rank groups
    x_title             = NULL,
    y_title             = NULL,
    main_title          = NULL,
    legend_title        = NULL,
    expr_min            = 0,      # threshold for "expressed" (expr > expr_min)
    show_legend         = TRUE,   # if FALSE, hide legend entirely
    expr_cap            = NULL,   # if set, cap expression at this value
    save_prefix         = NULL,   # if set, save PDF/SVG/PNG + group_stats TSV
    fig_width           = 4,      # fallback total figure width if panel_width_in is NULL
    fig_height          = 5,   # total figure height (inches); if NULL, auto by #groups
    panel_width_in      = NULL    # NEW: desired panel width (inches) for the plotting area
) {
  group_mode <- match.arg(group_mode)
  top_stat   <- match.arg(top_stat)  # "median" or "mean"
  
  # --- basic checks ---
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment.")
  }
  if (!assay_name %in% assayNames(sce)) {
    stop("Assay '", assay_name, "' not found in assayNames(sce).")
  }
  if (!gene_id %in% rownames(sce)) {
    stop("Gene '", gene_id, "' not found in rownames(sce).")
  }
  
  meta <- as.data.table(as.data.frame(colData(sce)))
  
  if (!tissue_col %in% colnames(meta)) {
    stop("tissue_col '", tissue_col, "' not found in colData(sce).")
  }
  if (!cell_type_col %in% colnames(meta)) {
    stop("cell_type_col '", cell_type_col, "' not found in colData(sce).")
  }
  
  # choose grouping variable
  if (group_mode == "cell_type") {
    grouping_vec <- as.character(meta[[cell_type_col]])
    group_label  <- cell_type_col
  } else {  # "tissue"
    grouping_vec <- as.character(meta[[tissue_col]])
    group_label  <- tissue_col
  }
  
  tissue_vec   <- as.character(meta[[tissue_col]])
  celltype_vec <- as.character(meta[[cell_type_col]])
  
  # --- expression ---
  expr_vec <- as.numeric(assay(sce, assay_name)[gene_id, ])
  
  dt <- data.table(
    expr     = expr_vec,
    group    = grouping_vec,
    tissue   = tissue_vec,
    celltype = celltype_vec
  )
  
  # drop missing group / tissue / celltype
  dt <- dt[!is.na(group) & !is.na(tissue) & !is.na(celltype)]
  
  if (nrow(dt) == 0) {
    stop("No cells with non-NA group, tissue, and cell_ontology_class.")
  }
  
  # --- helper to parse filter/exclude vectors ---
  parse_vec <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) == 1L) {
      v <- strsplit(x, ",")[[1]]
    } else {
      v <- x
    }
    v <- trimws(v)
    v[nchar(v) > 0]
  }
  
  # --- apply tissue_filter (INCLUDE) ---
  tissue_filter_vec <- parse_vec(tissue_filter)
  if (!is.null(tissue_filter_vec)) {
    all_tissues <- sort(unique(dt$tissue))
    missing_t <- setdiff(tissue_filter_vec, all_tissues)
    if (length(missing_t) > 0) {
      warning(
        "The following tissues in tissue_filter are not present and will be ignored: ",
        paste(missing_t, collapse = ", ")
      )
    }
    
    keep_tissues <- intersect(tissue_filter_vec, all_tissues)
    if (length(keep_tissues) == 0) {
      stop("After applying tissue_filter, no tissues remain in the data.")
    }
    
    cat("Restricting to tissues in tissue_filter:\n  - ",
        paste(keep_tissues, collapse = "\n  - "), "\n", sep = "")
    
    dt <- dt[tissue %in% keep_tissues]
    if (nrow(dt) == 0) {
      stop("No cells remain after applying tissue_filter.")
    }
  }
  
  # --- apply tissue_exclude (EXCLUDE) ---
  tissue_exclude_vec <- parse_vec(tissue_exclude)
  if (!is.null(tissue_exclude_vec)) {
    all_tissues <- sort(unique(dt$tissue))
    missing_t_ex <- setdiff(tissue_exclude_vec, all_tissues)
    if (length(missing_t_ex) > 0) {
      warning(
        "The following tissues in tissue_exclude are not present and will be ignored: ",
        paste(missing_t_ex, collapse = ", ")
      )
    }
    rm_tissues <- intersect(tissue_exclude_vec, all_tissues)
    if (length(rm_tissues) > 0) {
      cat("Excluding tissues in tissue_exclude:\n  - ",
          paste(rm_tissues, collapse = "\n  - "), "\n", sep = "")
      dt <- dt[!tissue %in% rm_tissues]
    }
    if (nrow(dt) == 0) {
      stop("No cells remain after applying tissue_exclude.")
    }
  }
  
  # --- apply celltype_filter (INCLUDE) ---
  celltype_filter_vec <- parse_vec(celltype_filter)
  if (!is.null(celltype_filter_vec)) {
    all_celltypes <- sort(unique(dt$celltype))
    missing_ct <- setdiff(celltype_filter_vec, all_celltypes)
    if (length(missing_ct) > 0) {
      warning(
        "The following cell_ontology_class values in celltype_filter are not present and will be ignored: ",
        paste(missing_ct, collapse = ", ")
      )
    }
    
    keep_celltypes <- intersect(celltype_filter_vec, all_celltypes)
    if (length(keep_celltypes) == 0) {
      stop("After applying celltype_filter, no cell types remain in the data.")
    }
    
    cat("Restricting to cell_ontology_class in celltype_filter:\n  - ",
        paste(keep_celltypes, collapse = "\n  - "), "\n", sep = "")
    
    dt <- dt[celltype %in% keep_celltypes]
    if (nrow(dt) == 0) {
      stop("No cells remain after applying celltype_filter.")
    }
  }
  
  # --- apply celltype_exclude (EXCLUDE) ---
  celltype_exclude_vec <- parse_vec(celltype_exclude)
  if (!is.null(celltype_exclude_vec)) {
    all_celltypes <- sort(unique(dt$celltype))
    missing_ct_ex <- setdiff(celltype_exclude_vec, all_celltypes)
    if (length(missing_ct_ex) > 0) {
      warning(
        "The following cell_ontology_class values in celltype_exclude are not present and will be ignored: ",
        paste(missing_ct_ex, collapse = ", ")
      )
    }
    rm_celltypes <- intersect(celltype_exclude_vec, all_celltypes)
    if (length(rm_celltypes) > 0) {
      cat("Excluding cell_ontology_class in celltype_exclude:\n  - ",
          paste(rm_celltypes, collapse = "\n  - "), "\n", sep = "")
      dt <- dt[!celltype %in% rm_celltypes]
    }
    if (nrow(dt) == 0) {
      stop("No cells remain after applying celltype_exclude.")
    }
  }
  
  # --- median & mean expression + cell counts per group ---
  summary_dt <- dt[
    ,
    .(
      median_expr = median(expr, na.rm = TRUE),
      mean_expr   = mean(expr, na.rm = TRUE),
      n_cells     = .N
    ),
    by = group
  ]
  
  metric_col <- if (top_stat == "median") "median_expr" else "mean_expr"
  
  # sort by chosen metric, descending
  summary_dt <- summary_dt[order(-summary_dt[[metric_col]])]
  
  cat("Per-group ", assay_name, " expression of ", gene_id,
      " by ", group_label, " (sorted by ", top_stat, "):\n", sep = "")
  print(summary_dt)
  
  # --- filter by min_cells_per_group ---
  summary_dt_filt <- summary_dt[n_cells >= min_cells_per_group]
  
  if (nrow(summary_dt_filt) == 0) {
    warning("No groups have at least ", min_cells_per_group, " cells. Nothing to plot.")
    return(invisible(list(
      summary          = summary_dt,
      filtered_summary = summary_dt_filt,
      top_summary      = NULL,
      group_stats      = NULL,
      plot             = NULL
    )))
  }
  
  cat("\nAfter filtering groups with < ", min_cells_per_group,
      " cells, remaining groups:\n", sep = "")
  print(summary_dt_filt)
  
  # --- take top N groups by chosen statistic ---
  summary_dt_filt <- summary_dt_filt[order(-summary_dt_filt[[metric_col]])]
  top_n_eff   <- min(top_n, nrow(summary_dt_filt))
  top_summary <- summary_dt_filt[seq_len(top_n_eff)]
  
  cat("\nTop ", top_n_eff, " ", group_label,
      " groups (after min_cells filter; ranked by ", top_stat, "):\n", sep = "")
  print(top_summary)
  
  # --- determine display order for groups ---
  display_order <- as.character(top_summary$group)
  
  if (group_mode == "cell_type" && !is.null(celltype_filter_vec)) {
    candidate <- intersect(celltype_filter_vec, display_order)
    if (length(candidate) > 0) {
      display_order <- candidate
    }
  } else if (group_mode == "tissue" && !is.null(tissue_filter_vec)) {
    candidate <- intersect(tissue_filter_vec, display_order)
    if (length(candidate) > 0) {
      display_order <- candidate
    }
  }
  
  # subset dt to top groups and set factor levels in desired order
  dt_top <- merge(
    dt,
    top_summary[, .(group, median_expr, mean_expr)],
    by = "group",
    all.x = FALSE,
    all.y = TRUE
  )
  dt_top[, group := factor(group, levels = display_order)]
  
  # --- cap expression if expr_cap is set ---
  if (!is.null(expr_cap)) {
    if (!is.numeric(expr_cap) || length(expr_cap) != 1L) {
      stop("expr_cap must be a single numeric value if provided.")
    }
    cat("Capping expression at expr_cap = ", expr_cap, " for plotting and stats.\n", sep = "")
    dt_top[, expr := pmin(expr, expr_cap)]
  }
  
  # --- per-group stats for output table (only for plotted groups, using possibly capped expr) ---
  group_stats <- dt_top[
    ,
    .(
      total_cells = .N,
      n_expr_pos  = sum(expr > expr_min, na.rm = TRUE),
      q1          = as.numeric(quantile(expr, 0.25, na.rm = TRUE)),
      median      = median(expr, na.rm = TRUE),
      mean        = mean(expr, na.rm = TRUE),
      q3          = as.numeric(quantile(expr, 0.75, na.rm = TRUE)),
      max         = max(expr, na.rm = TRUE)
    ),
    by = group
  ]
  
  cat("\nPer-group stats for ", gene_id, " (only plotted groups, using ",
      ifelse(is.null(expr_cap), "raw", paste0("capped at ", expr_cap)),
      " expr):\n", sep = "")
  print(group_stats)
  
  # --- default labels if not provided ---
  if (is.null(x_title)) {
    x_title <- paste0(assay_name, "(", gene_id, ")")
  }
  if (is.null(y_title)) {
    y_title <- group_label
  }
  if (is.null(main_title)) {
    main_title <- paste0(
      "Top ", top_n_eff, " ", group_label,
      " groups for ", gene_id,
      " (violin + black points; ≥", min_cells_per_group, " cells/group; ranked by ",
      top_stat, ")"
    )
  }
  if (is.null(legend_title)) {
    legend_title <- group_label
  }
  
  # --- build ggplot (panel layout only) ---
  p <- ggplot(dt_top, aes(x = expr, y = group)) +
    geom_violin(
      aes(fill = group),
      color = NA,
      alpha = 0.8,
      scale = "width",
      trim  = TRUE
    ) +
    geom_point(
      color    = "black",
      size     = 0.4,
      alpha    = 0.6,
      position = position_jitter(height = 0.1, width = 0)
    ) +
    scale_fill_manual(values = rep("grey80", length(unique(dt_top$group))), name = legend_title) +
    labs(
      y     = y_title,
      x     = x_title,
      title = main_title
    ) +
    theme_bw(base_size = 8) +
    theme(
      legend.text  = element_text(size = 8),
      legend.title = element_text(size = 8),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8)
    )
  
  # legend handling
  if (show_legend) {
    p <- p + guides(fill = guide_legend(reverse = TRUE)) +
      theme(legend.position = "right")
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  # --- convert to grob & fix panel width if requested ---
  if (!is.null(panel_width_in)) {
    g <- ggplotGrob(p)
    panel_cols <- unique(g$layout$l[g$layout$name == "panel"])
    g$widths[panel_cols] <- grid::unit(panel_width_in, "in")
  } else {
    g <- p
  }
  
  # --- print to current device ---
  if (inherits(g, "grob")) {
    grid::grid.newpage()
    grid::grid.draw(g)
  } else {
    print(g)
  }
  
  # --- saving if requested ---
  if (!is.null(save_prefix)) {
    dir.create(dirname(save_prefix), showWarnings = FALSE, recursive = TRUE)
    
    n_groups <- length(levels(dt_top$group))
    if (is.null(fig_height)) {
      height <- max(3, 0.25 * n_groups + 1)  # auto height based on #groups
    } else {
      height <- fig_height
    }
    
    # If we adjusted panel width, compute total width from grob
    if (!is.null(panel_width_in) && inherits(g, "grob")) {
      total_width <- sum(grid::convertWidth(g$widths, "in", valueOnly = TRUE))
    } else {
      total_width <- fig_width
    }
    
    pdf_file  <- paste0(save_prefix, ".pdf")
    svg_file  <- paste0(save_prefix, ".svg")
    png_file  <- paste0(save_prefix, ".png")
    tsv_file  <- paste0(save_prefix, "_group_stats.tsv")
    
    ggsave(pdf_file, plot = g, width = total_width, height = height)
    ggsave(svg_file, plot = g, width = total_width, height = height)
    ggsave(png_file, plot = g, width = total_width, height = height, dpi = 300)
    
    data.table::fwrite(group_stats, tsv_file, sep = "\t")
    
    cat("Saved:\n  ", pdf_file, "\n  ", svg_file, "\n  ",
        png_file, "\n  ", tsv_file, "\n", sep = "")
  }
  
  invisible(list(
    summary          = summary_dt,
    filtered_summary = summary_dt_filt,
    top_summary      = top_summary,
    group_stats      = group_stats,
    plot             = g
  ))
}


plot_gene_top_boxplot_by_group <- function(
    sce,
    gene_id,
    assay_name      = "logcounts",
    group_mode      = c("cell_type", "tissue"),
    cell_type_col   = "cell_ontology_class",
    tissue_col      = "tissue",
    top_n           = 5
) {
  group_mode <- match.arg(group_mode)
  
  # --- checks ---
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment.")
  }
  if (!assay_name %in% assayNames(sce)) {
    stop("Assay '", assay_name, "' not found in assayNames(sce).")
  }
  if (!gene_id %in% rownames(sce)) {
    stop("Gene '", gene_id, "' not found in rownames(sce).")
  }
  
  meta <- as.data.table(as.data.frame(colData(sce)))
  
  if (group_mode == "cell_type") {
    if (!cell_type_col %in% colnames(meta)) {
      stop("cell_type_col '", cell_type_col, "' not found in colData(sce).")
    }
    group_vec <- as.character(meta[[cell_type_col]])
    group_label <- cell_type_col
  } else if (group_mode == "tissue") {
    if (!tissue_col %in% colnames(meta)) {
      stop("tissue_col '", tissue_col, "' not found in colData(sce).")
    }
    group_vec <- as.character(meta[[tissue_col]])
    group_label <- tissue_col
  }
  
  # --- extract expression ---
  expr_vec <- as.numeric(assay(sce, assay_name)[gene_id, ])
  
  dt <- data.table(
    expr  = expr_vec,
    group = group_vec
  )
  
  # drop missing groups
  dt <- dt[!is.na(group)]
  
  if (nrow(dt) == 0) {
    stop("No cells with non-NA '", group_label, "'.")
  }
  
  # --- compute median expression per group ---
  summary_dt <- dt[
    ,
    .(
      median_expr = median(expr, na.rm = TRUE),
      n_cells     = .N
    ),
    by = group
  ][order(-median_expr)]
  
  cat("Median ", assay_name, " expression of ", gene_id,
      " by ", group_label, " (sorted descending):\n", sep = "")
  print(summary_dt)
  
  if (nrow(summary_dt) == 0) {
    warning("No groups to summarize.")
    return(invisible(list(summary = summary_dt, plot = NULL)))
  }
  
  # --- take top N groups ---
  top_n_eff <- min(top_n, nrow(summary_dt))
  top_summary <- summary_dt[seq_len(top_n_eff)]
  
  cat("\nTop ", top_n_eff, " ", group_label, " groups:\n", sep = "")
  print(top_summary)
  
  # subset dt to only top groups
  dt_top <- merge(
    dt,
    top_summary[, .(group, median_expr)],
    by = "group",
    all.x = FALSE,
    all.y = TRUE
  )
  
  # reorder factor by median expression (descending)
  dt_top[, group := factor(group, levels = top_summary$group)]
  
  # --- boxplot ---
  p <- ggplot(dt_top, aes(x = group, y = expr)) +
    geom_boxplot(outlier.size = 0.5) +
    coord_flip() +
    labs(
      x = group_label,
      y = paste0(assay_name, "(", gene_id, ")"),
      title = paste0(
        "Top ", top_n_eff, " ", group_label, " groups for ", gene_id
      )
    ) +
    theme_bw(base_size = 12)
  
  print(p)
  
  invisible(list(summary = summary_dt, top_summary = top_summary, plot = p))
}

add_umap_group_labels <- function(
    p,
    sce,
    dimred      = "UMAP_refined",
    group_col   = "tissue",     # or "cell_ontology_class", etc.
    groups      = NULL,         # vector of group names to label; NULL = all
    repel       = FALSE,        # use ggrepel if TRUE (and installed)
    text_size   = 3,
    text_color  = "black",
    fontface    = "bold",
    override_alpha = NULL,      # set to 1 for solid colors, NULL to keep original
    x_label     = NULL,         # X axis label; NULL = use default from dimred
    y_label     = NULL,         # Y axis label; NULL = use default from dimred
    legend_title = NULL,        # Legend title; NULL = use group_col
    plot_title  = NULL,         # Plot title; NULL = no title
    add_grid    = TRUE,         # Whether to add grid lines
    grid_color  = "grey90",     # Grid line color
    grid_size   = 0.3,          # Grid line size
    output_path = NULL,         # Base path for saving (without extension)
    width       = 8,            # Width for saved plots (inches)
    height      = 6,            # Height for saved plots (inches)
    dpi         = 1200           # DPI for PNG output
) {
  # --- checks ---
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment.")
  }
  if (!dimred %in% names(reducedDims(sce))) {
    stop("dimred '", dimred, "' not found in reducedDims(sce).")
  }
  meta <- as.data.frame(colData(sce))
  if (!group_col %in% colnames(meta)) {
    stop("group_col '", group_col, "' not found in colData(sce).")
  }
  
  # --- UMAP + metadata ---
  umap <- as.data.frame(reducedDim(sce, dimred))
  colnames(umap)[1:2] <- c("UMAP1", "UMAP2")
  
  df <- cbind(
    umap[, c("UMAP1", "UMAP2"), drop = FALSE],
    group = meta[[group_col]]
  )
  
  # drop NAs
  df <- df[!is.na(df$group), , drop = FALSE]
  
  # determine which groups to label
  all_groups <- sort(unique(df$group))
  
  # Handle empty character vector as "no labels"
  if (!is.null(groups) && length(groups) == 0) {
    label_groups <- character(0)
  } else if (is.null(groups)) {
    label_groups <- all_groups
  } else {
    label_groups <- intersect(all_groups, groups)
    missing <- setdiff(groups, all_groups)
    if (length(missing) > 0) {
      warning("The following requested groups are not present and will be ignored: ",
              paste(missing, collapse = ", "))
    }
  }
  
  # --- Override point alpha if requested ---
  if (!is.null(override_alpha)) {
    # Find geom_point layer(s) and modify alpha
    for (i in seq_along(p$layers)) {
      if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
        p$layers[[i]]$aes_params$alpha <- override_alpha
      }
    }
  }
  
  # --- Only add labels if there are groups to label ---
  if (length(label_groups) > 0) {
    # --- compute median UMAP coords per group ---
    label_df <- df |>
      dplyr::filter(group %in% label_groups) |>
      dplyr::group_by(group) |>
      dplyr::summarise(
        UMAP1 = median(UMAP1, na.rm = TRUE),
        UMAP2 = median(UMAP2, na.rm = TRUE),
        .groups = "drop"
      )
    
    # --- Extract colors from the plot ---
    # Try to get the color scale from the plot
    plot_built <- ggplot2::ggplot_build(p)
    
    # Get the data with colors assigned
    plot_data <- plot_built$data[[1]]
    
    # Create a mapping of group to color
    if ("colour" %in% names(plot_data) && "group" %in% names(df)) {
      # Match groups to colors
      color_map <- df |>
        dplyr::select(group) |>
        dplyr::distinct() |>
        dplyr::mutate(row_id = dplyr::row_number())
      
      # Get unique colors from plot data
      if (length(unique(plot_data$colour)) > 1) {
        unique_colors <- unique(plot_data$colour)
        unique_groups <- unique(df$group)
        
        # Create mapping
        if (length(unique_colors) == length(unique_groups)) {
          color_mapping <- data.frame(
            group = unique_groups,
            color = unique_colors,
            stringsAsFactors = FALSE
          )
          
          # Filter to only labeled groups
          labeled_colors <- color_mapping |>
            dplyr::filter(group %in% label_groups) |>
            dplyr::arrange(group)
          
          # Print to console
          cat("\n=== Group Labels and Colors ===\n")
          for (i in seq_len(nrow(labeled_colors))) {
            cat(sprintf("%-30s : %s\n", 
                        labeled_colors$group[i], 
                        labeled_colors$color[i]))
          }
          cat("================================\n\n")
        }
      }
    }
    
    # --- build label layer ---
    if (repel && requireNamespace("ggrepel", quietly = TRUE)) {
      label_layer <- ggrepel::geom_text_repel(
        data     = label_df,
        aes(x = UMAP1, y = UMAP2, label = group),
        color    = text_color,
        fontface = fontface,
        size     = text_size,
        max.overlaps = Inf
      )
    } else {
      if (repel) {
        warning("repel=TRUE requested but ggrepel not available; using geom_text instead.")
      }
      label_layer <- geom_text(
        data     = label_df,
        aes(x = UMAP1, y = UMAP2, label = group),
        color    = text_color,
        fontface = fontface,
        size     = text_size
      )
    }
    
    # Add label layer to plot
    p <- p + label_layer
  }
  
  # --- Add axis labels ---
  if (!is.null(x_label)) {
    p <- p + xlab(x_label)
  } else {
    # Default to UMAP_1 or similar based on dimred name
    default_x <- if (grepl("UMAP", dimred, ignore.case = TRUE)) {
      "UMAP_1"
    } else if (grepl("TSNE", dimred, ignore.case = TRUE)) {
      "tSNE_1"
    } else {
      paste0(dimred, "_1")
    }
    p <- p + xlab(default_x)
  }
  
  if (!is.null(y_label)) {
    p <- p + ylab(y_label)
  } else {
    # Default to UMAP_2 or similar based on dimred name
    default_y <- if (grepl("UMAP", dimred, ignore.case = TRUE)) {
      "UMAP_2"
    } else if (grepl("TSNE", dimred, ignore.case = TRUE)) {
      "tSNE_2"
    } else {
      paste0(dimred, "_2")
    }
    p <- p + ylab(default_y)
  }
  
  # --- Add legend title ---
  if (!is.null(legend_title)) {
    p <- p + labs(color = legend_title, fill = legend_title)
  } else {
    # Use the group_col as default legend title
    p <- p + labs(color = group_col, fill = group_col)
  }
  
  # --- Add plot title ---
  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title)
  }
  
  # --- Add grid lines ---
  if (add_grid) {
    p <- p + theme(
      panel.grid.major = element_line(
        color = grid_color,
        size = grid_size,
        linetype = "solid"
      ),
      panel.grid.minor = element_line(
        color = grid_color,
        size = grid_size * 0.5,
        linetype = "dotted"
      )
    )
  } else {
    # Remove grid lines if requested
    p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }
  
  # --- Save plots if output_path is provided ---
  if (!is.null(output_path)) {
    # Ensure the directory exists
    output_dir <- dirname(output_path)
    if (!dir.exists(output_dir) && output_dir != ".") {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Remove any extension from output_path if present
    output_base <- sub("\\.[^.]*$", "", output_path)
    
    # Save PNG
    png_file <- paste0(output_base, ".png")
    ggsave(
      filename = png_file,
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      units = "in",
      bg="white"
    )
    cat("Saved PNG: ", png_file, "\n", sep = "")
    
    # Save SVG
    svg_file <- paste0(output_base, ".svg")
    ggsave(
      filename = svg_file,
      plot = p,
      width = width,
      height = height,
      units = "in",
      device = "svg",
      bg="white"
    )
    cat("Saved SVG: ", svg_file, "\n", sep = "")
    
    # Save PDF
    pdf_file <- paste0(output_base, ".pdf")
    ggsave(
      filename = pdf_file,
      plot = p,
      width = width,
      height = height,
      units = "in",
      device = "pdf"
    )
    cat("Saved PDF: ", pdf_file, "\n", sep = "")
  }
  
  # Return the modified plot
  return(p)
}

plot_gene_by_cellgroup_in_tissue <- function(
    sce,
    gene_id,
    tissue_value,
    cell_group_col = "cell_ontology_class",
    assay_name     = "logcounts"
) {
  # --- checks ---
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment.")
  }
  if (!assay_name %in% assayNames(sce)) {
    stop("Assay '", assay_name, "' not found in assayNames(sce).")
  }
  if (!"tissue" %in% colnames(colData(sce))) {
    stop("colData(sce) must contain a 'tissue' column.")
  }
  if (!cell_group_col %in% colnames(colData(sce))) {
    stop("cell_group_col '", cell_group_col, "' not found in colData(sce).")
  }
  if (!gene_id %in% rownames(sce)) {
    stop("Gene '", gene_id, "' not found in rownames(sce).")
  }
  
  # --- subset to requested tissue ---
  tissue_vec <- as.character(colData(sce)$tissue)
  keep <- tissue_vec == tissue_value
  
  if (!any(keep)) {
    stop("No cells found with tissue == '", tissue_value, "'.")
  }
  
  sce_sub <- sce[, keep, drop = FALSE]
  cat("Subsetting to tissue == '", tissue_value, "': ",
      ncol(sce_sub), " cells.\n", sep = "")
  
  # --- extract expression + grouping column ---
  expr_vec <- as.numeric(assay(sce_sub, assay_name)[gene_id, ])
  group_vec <- as.character(colData(sce_sub)[[cell_group_col]])
  
  dt <- data.table(
    expr      = expr_vec,
    group     = group_vec
  )
  
  # drop missing groups
  dt <- dt[!is.na(group)]
  
  if (nrow(dt) == 0) {
    stop("No cells with non-NA '", cell_group_col,
         "' in tissue '", tissue_value, "'.")
  }
  
  # --- summarize by group ---
  summary_dt <- dt[
    ,
    .(
      median_expr = median(expr, na.rm = TRUE),
      n_cells     = .N
    ),
    by = group
  ][order(-median_expr)]
  
  cat("\nMedian ", assay_name, " expression of ", gene_id,
      " in tissue '", tissue_value, "' by ", cell_group_col, ":\n", sep = "")
  print(summary_dt)
  
  # reorder groups in plot by median expression
  dt <- merge(
    dt,
    summary_dt[, .(group, median_expr)],
    by = "group",
    all.x = TRUE
  )
  dt[, group := factor(group, levels = summary_dt$group)]
  
  # --- boxplot ---
  p <- ggplot(dt, aes(x = group, y = expr)) +
    geom_boxplot(outlier.size = 0.5) +
    coord_flip() +
    labs(
      x = cell_group_col,
      y = paste0(assay_name, "(", gene_id, ")"),
      title = paste0(
        "Expression of ", gene_id,
        " in tissue '", tissue_value,
        "' by ", cell_group_col
      )
    ) +
    theme_bw(base_size = 12)
  
  print(p)
  
  invisible(list(summary = summary_dt, plot = p, n_cells = ncol(sce_sub)))
}

plot_gene_top15_boxplot <- function(
    sce,
    gene_id,
    assay_name    = "logcounts",
    cell_type_col = "cell_ontology_class",
    tissue_col    = "tissue"
) {
  # --- checks ---
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment.")
  }
  if (!assay_name %in% assayNames(sce)) {
    stop("Assay '", assay_name, "' not found in assayNames(sce).")
  }
  if (!gene_id %in% rownames(sce)) {
    stop("Gene '", gene_id, "' not found in rownames(sce).")
  }
  if (!cell_type_col %in% colnames(colData(sce))) {
    stop("cell_type_col '", cell_type_col, "' not found in colData(sce).")
  }
  if (!tissue_col %in% colnames(colData(sce))) {
    stop("tissue_col '", tissue_col, "' not found in colData(sce).")
  }
  
  # --- extract expression and metadata ---
  expr_vec <- as.numeric(assay(sce, assay_name)[gene_id, ])
  ct       <- as.character(colData(sce)[[cell_type_col]])
  tis      <- as.character(colData(sce)[[tissue_col]])
  
  dt <- data.table(
    expr      = expr_vec,
    cell_type = ct,
    tissue    = tis
  )
  
  # drop rows with missing annotations
  dt <- dt[!is.na(cell_type) & !is.na(tissue)]
  
  # --- compute median expression by (cell_type, tissue) ---
  summary_dt <- dt[
    ,
    .(median_expr = median(expr, na.rm = TRUE),
      n_cells     = .N),
    by = .(cell_type, tissue)
  ][order(-median_expr)]
  
  cat("Median ", assay_name, " expression of ", gene_id,
      " by cell_type x tissue (sorted descending):\n", sep = "")
  print(summary_dt)
  
  if (nrow(summary_dt) == 0) {
    warning("No valid (cell_type, tissue) groups to summarize.")
    return(invisible(list(summary = summary_dt, plot = NULL)))
  }
  
  # --- take top 5 groups by median expression ---
  top_n <- min(15L, nrow(summary_dt))
  top_summary <- summary_dt[seq_len(top_n)]
  
  cat("\nTop ", top_n, " cell_type x tissue groups:\n", sep = "")
  print(top_summary)
  
  # build a group label
  top_summary[, group := paste(cell_type, tissue, sep = " | ")]
  
  # join back to cell-level data
  dt[, group := paste(cell_type, tissue, sep = " | ")]
  dt_top <- merge(
    dt,
    top_summary[, .(group, median_expr)],
    by = "group",
    all.x = FALSE,
    all.y = TRUE
  )
  
  # reorder factor by median expression
  dt_top[, group := factor(group, levels = unique(top_summary$group[order(top_summary$median_expr)]) )]
  
  # --- boxplot for top 5 groups ---
  p <- ggplot(dt_top, aes(x = group, y = expr)) +
    geom_boxplot(outlier.size = 0.5) +
    coord_flip() +
    labs(
      x = "cell type | tissue",
      y = paste0(assay_name, "(", gene_id, ")"),
      title = paste("Top", top_n, "cell type / tissue groups for", gene_id)
    ) +
    theme_bw(base_size = 12)
  
  print(p)
  
  invisible(list(summary = summary_dt, top_summary = top_summary, plot = p))
  setorder(summary_dt, -median_expr)
  return(summary_dt)
}

plot_gene_umap_size_color <- function(
    sce,
    gene_id,
    dimred      = "UMAP_refined",
    assay_name  = "logcounts",
    size_range  = c(0.2, 1.),
    color_limits = NULL,        # c(min, max) for color scale; NULL = auto
    color_map   = "turbo",      # NEW: viridis color palette option
    nothing     = FALSE,
    title       = NULL,
    x_label     = "UMAP 1",
    y_label     = "UMAP 2",
    legend_title = "log(CPM+1)",
    xlim        = NULL,
    ylim        = NULL,
    aspect_ratio = 1,
    add_grid    = TRUE,
    grid_color  = "grey90",
    grid_size   = 0.3,
    output_path = NULL,
    width       = 8,
    height      = 6,
    dpi         = 1200
) {
  # --- checks ---
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment.")
  }
  if (!assay_name %in% assayNames(sce)) {
    stop("Assay '", assay_name, "' not found in assayNames(sce).")
  }
  if (!dimred %in% names(reducedDims(sce))) {
    stop("dimred '", dimred, "' not found in reducedDims(sce).")
  }
  if (!gene_id %in% rownames(sce)) {
    stop("Gene '", gene_id, "' not found in rownames(sce).")
  }
  
  # --- extract UMAP coordinates + expression ---
  coords <- reducedDim(sce, dimred)
  df <- data.frame(
    UMAP1 = coords[, 1],
    UMAP2 = coords[, 2],
    expr  = as.numeric(assay(sce, assay_name)[gene_id, ]),
    stringsAsFactors = FALSE
  )
  
  # clean up expression (avoid -Inf/NA)
  df$expr[!is.finite(df$expr)] <- 0
  
  # order by expression so low expr plotted first, high expr last
  df <- df[order(df$expr), ]
  
  # scale point sizes based on expression
  df$size_scaled <- rescale(df$expr, to = size_range)
  
  # --- Set default title (ignored if nothing=TRUE) ---
  if (is.null(title)) {
    title <- paste0(dimred, " - ", gene_id, " expression (", assay_name, ")")
  }
  
  # --- Build base plot ---
  p <- ggplot(df, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = expr, size = size_scaled), alpha = 0.9) +
    scale_size(range = size_range, guide = "none")
  
  # Add color scale with custom color_map
  if (!is.null(color_limits)) {
    p <- p + scale_colour_viridis_c(option = color_map, limits = color_limits)
  } else {
    p <- p + scale_colour_viridis_c(option = color_map)
  }
  
  # --- Set coordinate system with aspect ratio and limits ---
  if (!is.null(aspect_ratio)) {
    # Fixed aspect ratio
    if (!is.null(xlim) && !is.null(ylim)) {
      p <- p + coord_fixed(ratio = aspect_ratio, xlim = xlim, ylim = ylim)
    } else if (!is.null(xlim)) {
      p <- p + coord_fixed(ratio = aspect_ratio, xlim = xlim)
    } else if (!is.null(ylim)) {
      p <- p + coord_fixed(ratio = aspect_ratio, ylim = ylim)
    } else {
      p <- p + coord_fixed(ratio = aspect_ratio)
    }
  } else {
    # Free aspect ratio
    if (!is.null(xlim)) {
      p <- p + xlim(xlim)
    }
    if (!is.null(ylim)) {
      p <- p + ylim(ylim)
    }
  }
  
  # --- Add theme ---
  p <- p + theme_bw()
  
  # --- Handle "nothing" mode ---
  if (nothing) {
    # Strip everything: title, axis labels, axis text, legend
    p <- p + theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
  } else {
    # Normal mode: add title, axis labels, legend title
    
    # --- Add title ---
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    }
    
    # --- Add axis labels ---
    if (!is.null(x_label)) {
      p <- p + xlab(x_label)
    } else {
      # Default to UMAP_1 or similar based on dimred name
      default_x <- if (grepl("UMAP", dimred, ignore.case = TRUE)) {
        "UMAP_1"
      } else if (grepl("TSNE", dimred, ignore.case = TRUE)) {
        "tSNE_1"
      } else {
        paste0(dimred, "_1")
      }
      p <- p + xlab(default_x)
    }
    
    if (!is.null(y_label)) {
      p <- p + ylab(y_label)
    } else {
      # Default to UMAP_2 or similar based on dimred name
      default_y <- if (grepl("UMAP", dimred, ignore.case = TRUE)) {
        "UMAP_2"
      } else if (grepl("TSNE", dimred, ignore.case = TRUE)) {
        "tSNE_2"
      } else {
        paste0(dimred, "_2")
      }
      p <- p + ylab(default_y)
    }
    
    # --- Add legend title ---
    if (!is.null(legend_title)) {
      p <- p + labs(color = legend_title)
    } else {
      # Use the assay_name as default legend title
      p <- p + labs(color = assay_name)
    }
  }
  
  # --- Add grid lines ---
  if (add_grid) {
    p <- p + theme(
      panel.grid.major = element_line(
        color = grid_color,
        size = grid_size,
        linetype = "solid"
      ),
      panel.grid.minor = element_line(
        color = grid_color,
        size = grid_size * 0.5,
        linetype = "dotted"
      )
    )
  } else {
    # Remove grid lines if requested
    p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }
  
  # --- Save plots if output_path is provided ---
  if (!is.null(output_path)) {
    # Ensure the directory exists
    output_dir <- dirname(output_path)
    if (!dir.exists(output_dir) && output_dir != ".") {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Remove any extension from output_path if present
    output_base <- sub("\\.[^.]*$", "", output_path)
    
    # Save PNG
    png_file <- paste0(output_base, ".png")
    ggsave(
      filename = png_file,
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      units = "in"
    )
    cat("Saved PNG: ", png_file, "\n", sep = "")
    
    # Save SVG
    svg_file <- paste0(output_base, ".svg")
    ggsave(
      filename = svg_file,
      plot = p,
      width = width,
      height = height,
      units = "in",
      device = "svg"
    )
    cat("Saved SVG: ", svg_file, "\n", sep = "")
    
    # Save PDF
    pdf_file <- paste0(output_base, ".pdf")
    ggsave(
      filename = pdf_file,
      plot = p,
      width = width,
      height = height,
      units = "in",
      device = "pdf"
    )
    cat("Saved PDF: ", pdf_file, "\n", sep = "")
  }
  
  # Print and return the plot
  print(p)
  invisible(p)
}

plot_genes_umap_batch <- function(
    sce,
    genes,
    output_dir,
    dimred = "UMAP_refined",
    assay_name = "logcounts",
    nothing     = TRUE,
    size_range = c(0.1, 0.4),
    aspect_ratio = 0.75,
    color_limits = c(0, 10),
    width = 8,
    height = 6,
    dpi = 1200,
    add_grid = TRUE,
    grid_color = "grey90",
    grid_size = 0.3,
    color_map = "turbo"
) {
  # Validate inputs
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment object.")
  }
  if (!dimred %in% names(reducedDims(sce))) {
    stop("dimred '", dimred, "' not found in reducedDims(sce).")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Track results
  found_genes <- character()
  missing_genes <- character()
  error_genes <- character()
  
  cat("\n========================================\n")
  cat("Plotting", length(genes), "genes to:", output_dir, "\n")
  cat("Using dimred:", dimred, "\n")
  cat("========================================\n\n")
  
  for (gene in genes) {
    cat("--- Processing gene:", gene, "---\n")
    
    # Check if gene exists in the data
    if (!gene %in% rownames(sce)) {
      cat("  Gene '", gene, "' not found in data. Skipping.\n\n", sep = "")
      missing_genes <- c(missing_genes, gene)
      next
    }
    
    # Try to generate the plot
    tryCatch({
      # Create unique output path for this gene
      output_path <- file.path(output_dir, paste0(gene, "_", dimred))
      
      plot_gene_umap_size_color(
        sce,
        gene_id      = gene,
        dimred       = dimred,
        assay_name   = assay_name,
        title        = paste0(gene, " expression"),
        size_range   = size_range,
        aspect_ratio = aspect_ratio,
        color_limits = color_limits,
        width        = width,
        height       = height,
        dpi          = dpi,
        add_grid     = add_grid,
        grid_color   = grid_color,
        grid_size    = grid_size,
        output_path  = output_path,
        nothing = nothing,
        color_map = color_map
      )
      
      cat("  Successfully plotted:", gene, "\n\n")
      found_genes <- c(found_genes, gene)
      
    }, error = function(e) {
      cat("  Error plotting gene '", gene, "': ", e$message, "\n", sep = "")
      cat("  Continuing to next gene...\n\n")
      error_genes <- c(error_genes, gene)
    })
  }
  
  # Summary
  cat("\n========================================\n")
  cat("SUMMARY\n")
  cat("========================================\n")
  cat("Total genes requested:", length(genes), "\n")
  cat("Successfully plotted:  ", length(found_genes), "\n")
  cat("Not found in data:     ", length(missing_genes), "\n")
  cat("Errors during plotting:", length(error_genes), "\n")
  
  if (length(missing_genes) > 0) {
    cat("\nMissing genes:\n")
    cat("  ", paste(missing_genes, collapse = ", "), "\n", sep = "")
  }
  
  if (length(error_genes) > 0) {
    cat("\nGenes with errors:\n")
    cat("  ", paste(error_genes, collapse = ", "), "\n", sep = "")
  }
  
  # Return summary invisibly
  invisible(list(
    found = found_genes,
    missing = missing_genes,
    errors = error_genes
  ))
}


