# =============================================================================
# ValFunctions.R
# -----------------------------------------------------------------------------
# Utility functions for scRNA-seq pseudotime analysis:
#   1. ComplexHeatmap visualisation (refactored — 5 duplicate blocks -> 2 fns)
#   2. Mutual Information computation for GRN inference
#   3. Graph export/import (JSON / GML)
#   4. Mouse -> Human gene symbol conversion via BioMart
#   5. Expression curve visualisation along pseudotime
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin
# Project : scRNA-seq Human Pre-implantation Embryo - Lineage Specification
# Paper   : Meistermann D. et al., Cell Stem Cell, 2021
#           DOI: 10.1016/j.stem.2021.04.027
# Date    : 2018 (MSc M2 internship, CRTI UMR 1064, Nantes Universite)
# =============================================================================


# =============================================================================
# DEPENDENCIES
# =============================================================================

.load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Required package missing: '", pkg,
         "'. Install with: BiocManager::install('", pkg, "')")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


# =============================================================================
# I. COMPLEXHEATMAP UTILITIES
# =============================================================================
# REFACTORING NOTE
# The original file contained 5 near-identical functions
# (ComplexHeatMapLineage, ComplexHeatMapCluster, SmoothComplexHeatMap,
#  SmoothComplexHeatMap2, SmoothComplexHeatMap3), each repeating the same
# 5-lineage Heatmap construction block.
# Replaced by two composable functions:
#   build_lineage_heatmap()  -> builds one Heatmap object for one lineage
#   draw_multi_heatmap()     -> assembles list of Heatmaps and exports PDF
# All original use cases are reproduced via the LINEAGE_CONFIG examples below.
# =============================================================================

#' Build a single-lineage Heatmap object
#'
#' Subsets the expression matrix to cells belonging to a given lineage,
#' removes zero-variance genes, and constructs a ComplexHeatmap Heatmap.
#'
#' @param mat              Matrix or data.frame (genes x cells).
#' @param annot_df         data.frame. Sample annotation for this lineage.
#' @param annot_col        Character or integer. Column(s) to use for annotation bar.
#' @param col_annot        Named list. Colour mapping for annotation bar.
#'                         e.g. list(Lineage = c("EPI" = "red"))
#' @param panel_name       Character. Heatmap panel name (legend label).
#' @param panel_title      Character. Column title shown above the panel.
#' @param col_ht           colorRamp2 object. Colour scale for expression values.
#' @param cluster_rows     Logical or dendrogram. Row clustering specification.
#' @param cluster_dist     Character or NULL. Row clustering distance method.
#' @param cluster_method   Character or NULL. Row clustering linkage method.
#' @param split            Vector or NULL. Row split specification.
#' @param panel_width      grid::unit or NULL. Override default panel width.
#' @param filter_zero_var  Logical. Remove genes with zero variance. Default TRUE.
#'
#' @return A ComplexHeatmap Heatmap object.
#' @examples
#' \dontrun{
#' col_ht <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
#' ht_epi <- build_lineage_heatmap(
#'   mat         = expr_mat,
#'   annot_df    = sample_annot[sample_annot$Lineage == "6.Epiblast", ],
#'   annot_col   = "Lineage",
#'   col_annot   = list(Lineage = c("6.Epiblast" = "red")),
#'   panel_name  = "EPI",
#'   panel_title = "Epiblast",
#'   col_ht      = col_ht
#' )
#' }
build_lineage_heatmap <- function(mat,
                                   annot_df,
                                   annot_col,
                                   col_annot,
                                   panel_name,
                                   panel_title,
                                   col_ht,
                                   cluster_rows   = TRUE,
                                   cluster_dist   = NULL,
                                   cluster_method = NULL,
                                   split          = NULL,
                                   panel_width    = NULL,
                                   filter_zero_var = TRUE) {
  .load_pkg("ComplexHeatmap")

  # Subset to lineage cells
  cell_ids <- rownames(annot_df)
  mat_sub  <- mat[, colnames(mat) %in% cell_ids, drop = FALSE]

  # Remove zero-variance genes to avoid clustering errors
  if (filter_zero_var) {
    keep    <- apply(mat_sub, 1, sd, na.rm = TRUE) != 0
    mat_sub <- mat_sub[keep, , drop = FALSE]
  }

  # Build annotation bar
  annot_data <- as.data.frame(annot_df[, annot_col, drop = FALSE])
  if (is.integer(annot_col)) colnames(annot_data) <- names(col_annot)
  ha <- HeatmapAnnotation(df = annot_data, col = col_annot)

  # Compose Heatmap
  ht <- Heatmap(
    mat_sub,
    name                     = panel_name,
    column_title             = panel_title,
    row_title                = "Genes",
    col                      = col_ht,
    top_annotation           = ha,
    row_names_gp             = gpar(cex = 0.1),
    column_names_gp          = gpar(cex = 0.2),
    column_title_gp          = gpar(cex = 0.5),
    show_heatmap_legend      = TRUE,
    show_row_names           = TRUE,
    show_column_names        = FALSE,
    cluster_columns          = FALSE,
    column_dend_reorder      = FALSE,
    cluster_rows             = cluster_rows,
    clustering_distance_rows = cluster_dist,
    clustering_method_rows   = cluster_method,
    split                    = split,
    gap                      = unit(5, "mm")
  )

  if (!is.null(panel_width)) ht@matrix_param$width <- panel_width
  ht
}


#' Assemble and export a multi-lineage ComplexHeatmap to PDF
#'
#' @param heatmap_list List of Heatmap objects (one per lineage, in display order).
#' @param gene_group   Character. Label used in output filename and plot title.
#' @param output_dir   Character. Directory for PDF output. Default: "fig".
#' @param pdf_height   Numeric. PDF height in inches. Default: 15.
#' @param pdf_width    Numeric. PDF width in inches. Default: 10.
#'
#' @return Invisible NULL. Saves PDF to output_dir/HeatMapComplex_<gene_group>.pdf
draw_multi_heatmap <- function(heatmap_list,
                                gene_group,
                                output_dir = "fig",
                                pdf_height = 15,
                                pdf_width  = 10) {
  .load_pkg("ComplexHeatmap")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  combined <- Reduce(`+`, heatmap_list)
  pdf_path <- file.path(output_dir,
                         paste0("HeatMapComplex_", gene_group, ".pdf"))
  pdf(pdf_path, height = pdf_height, width = pdf_width)
  draw(combined,
       row_title         = paste("Gene Clustering -", gene_group),
       row_title_gp      = gpar(col = "red"),
       column_title      = "Sample Pseudotime",
       column_title_side = "bottom",
       gap               = unit(0.2, "cm"))
  dev.off()
  message("Heatmap saved: ", pdf_path)
  invisible(NULL)
}


# =============================================================================
# EXAMPLE: Reproducing original SmoothComplexHeatMap() behaviour
# =============================================================================
# col_ht <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
#
# lineage_config <- list(
#   list(
#     filter    = sample_annot$Lineage == "3.Morula4.B1",
#     col_annot = list(Lineage = c("3.Morula4.B1" = "black")),
#     name      = "MorulaB1", title = "Morula"
#   ),
#   list(
#     filter    = sample_annot$Lineage == "5. Unspecified ICM",
#     col_annot = list(Lineage = c("5. Unspecified ICM" = "yellow")),
#     name      = "ICM", title = "ICM"
#   ),
#   list(
#     filter    = sample_annot$Lineage == "6.Epiblast",
#     col_annot = list(Lineage = c("6.Epiblast" = "red")),
#     name      = "EPI", title = "EPI"
#   ),
#   list(
#     filter    = sample_annot$Lineage == "7.Primitive endoderm",
#     col_annot = list(Lineage = c("7.Primitive endoderm" = "green")),
#     name      = "PE", title = "PE",
#     width     = unit(1, "cm")
#   ),
#   list(
#     filter    = sample_annot$Lineage == "9.Trophectoderm",
#     col_annot = list(Lineage = c("9.Trophectoderm" = "blue")),
#     name      = "TE", title = "TE"
#   )
# )
#
# heatmaps <- lapply(lineage_config, function(cfg) {
#   build_lineage_heatmap(
#     mat          = expr_mat,
#     annot_df     = sample_annot[cfg$filter, ],
#     annot_col    = "Lineage",
#     col_annot    = cfg$col_annot,
#     panel_name   = cfg$name,
#     panel_title  = cfg$title,
#     col_ht       = col_ht,
#     split        = gene_split,
#     panel_width  = cfg$width
#   )
# })
# draw_multi_heatmap(heatmaps, gene_group = "Kinome", pdf_height = 10)


# =============================================================================
# II. MUTUAL INFORMATION (GRN inference)
# =============================================================================

#' Compute pairwise or one-vs-all Mutual Information between genes
#'
#' @param expr_log    data.frame (genes x cells). Log-normalised expression.
#' @param method_mi   Character. MI estimator: 'emp', 'mm', 'shrink', or 'sg'.
#' @param disc        Character. Discretisation: 'equalfreq', 'equalwidth',
#'                    or 'globalequalwidth'.
#' @param n_bins      Numeric. Number of bins. Default: n_genes^(1/3).
#' @param discretise  Logical. Discretise before MI computation.
#' @param one_gene    Logical. If TRUE, compute MI of one gene vs all others.
#' @param gene        Character. Target gene (required if one_gene = TRUE).
#'
#' @return Named numeric vector (one_gene = TRUE) or full MI matrix.
compute_mutual_information <- function(expr_log,
                                        method_mi  = "emp",
                                        disc       = "equalfreq",
                                        n_bins     = nrow(expr_log)^(1 / 3),
                                        discretise = TRUE,
                                        one_gene   = TRUE,
                                        gene       = NULL) {
  .load_pkg("infotheo")

  valid_mi   <- c("emp", "mm", "shrink", "sg")
  valid_disc <- c("equalfreq", "equalwidth", "globalequalwidth")
  if (!method_mi %in% valid_mi)
    stop("Invalid MI method. Choose from: ", paste(valid_mi, collapse = ", "))
  if (discretise && !disc %in% valid_disc)
    stop("Invalid discretisation method. Choose from: ",
         paste(valid_disc, collapse = ", "))

  expr_t <- t(expr_log)  # infotheo expects cells x genes

  if (discretise) {
    message("Discretising (method: ", disc, ", bins: ", round(n_bins), ")...")
    expr_t <- discretize(expr_t, disc = disc, nbins = n_bins)
  }

  message("Computing mutual information (method: ", method_mi, ")...")

  if (one_gene) {
    if (is.null(gene)) stop("Argument 'gene' required when one_gene = TRUE.")
    if (!gene %in% colnames(expr_t)) stop("Gene not found: ", gene)
    mi        <- vapply(seq_len(ncol(expr_t)), function(i)
      mutinformation(expr_t[, i], expr_t[, gene], method = method_mi),
      numeric(1))
    names(mi) <- colnames(expr_t)
    return(mi)
  } else {
    message("Full pairwise MI matrix (may take several minutes)...")
    return(mutinformation(expr_t, method = method_mi))
  }
}


# =============================================================================
# III. GRAPH EXPORT / IMPORT  (JSON & GML)
# =============================================================================

#' Export an igraph object to JSON
#'
#' @param g        igraph object.
#' @param filename Character or NULL. Output path (.json). If NULL, returns string.
#' @return JSON string (invisible if filename given).
export_graph_json <- function(g, filename = NULL) {
  .load_pkg("igraph")
  .load_pkg("jsonlite")

  graph_list           <- list()
  graph_list$edges     <- as_data_frame(g, what = "edges")
  graph_list$vertices  <- as_data_frame(g, what = "vertices")
  rownames(graph_list$vertices) <- NULL
  if (ncol(graph_list$vertices) == 0) graph_list$vertices <- NULL
  graph_list$directed  <- is.directed(g)
  graph_list$name      <- g$name

  json_str <- toJSON(graph_list, pretty = TRUE)

  if (!is.null(filename)) {
    writeLines(json_str, filename)
    message("Graph exported to: ", filename)
    invisible(json_str)
  } else {
    return(json_str)
  }
}


#' Import a JSON file into an igraph object
#'
#' @param filename Character. Path to JSON file or raw JSON string.
#' @return igraph object.
import_graph_json <- function(filename) {
  .load_pkg("igraph")
  .load_pkg("jsonlite")

  built <- fromJSON(filename, flatten = TRUE)
  g <- if ("vertices" %in% names(built)) {
    graph_from_data_frame(built$edges,
                          directed = built$directed,
                          vertices = built$vertices)
  } else {
    graph_from_data_frame(built$edges, directed = built$directed)
  }
  if ("name" %in% names(built)) g$name <- built$name
  return(g)
}


#' Export an igraph object to GML format
#'
#' @param graph    igraph object.
#' @param filename Character. Output file path (.gml).
#' @return Invisible NULL.
export_graph_gml <- function(graph, filename) {
  .load_pkg("igraph")

  con <- file(filename, "w")
  on.exit(close(con))  # ensure file is closed even on error

  cat("Creator \"igraph exportGML\"\nVersion 1.0\ngraph\n[\n", file = con)
  cat("  directed", as.integer(is.directed(graph)), "\n", file = con)

  for (i in seq_len(vcount(graph))) {
    cat("  node\n  [\n", file = con)
    cat("    id", i - 1, "\n", file = con)
    cat("    graphics\n    [\n", file = con)
    cat("      fill \"", V(graph)$color[i], "\"\n", sep = "", file = con)
    cat("      type \"rectangle\"\n", file = con)
    cat("      outline \"#000000\"\n", file = con)
    cat("    ]\n  ]\n", file = con)
  }

  el <- get.edgelist(graph, names = FALSE)
  for (i in seq_len(nrow(el))) {
    cat("  edge\n  [\n", file = con)
    cat("    source", el[i, 1], "\n", file = con)
    cat("    target", el[i, 2], "\n", file = con)
    cat("  ]\n", file = con)
  }
  cat("]\n", file = con)
  message("Graph exported to GML: ", filename)
  invisible(NULL)
}


# =============================================================================
# IV. MOUSE -> HUMAN GENE SYMBOL CONVERSION
# =============================================================================

#' Convert mouse MGI symbols to human HGNC symbols via Ensembl BioMart
#'
#' @param gene_list Character vector. Mouse gene symbols.
#' @return Character vector. Unique human HGNC symbols.
convert_mouse_to_human <- function(gene_list) {
  .load_pkg("biomaRt")

  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  genes_v2 <- getLDS(
    attributes  = "mgi_symbol",
    filters     = "mgi_symbol",
    values      = gene_list,
    mart        = mouse,
    attributesL = "hgnc_symbol",
    martL       = human,
    uniqueRows  = TRUE
  )
  human_genes <- unique(genes_v2[, 2])
  message("Converted ", length(gene_list), " mouse genes -> ",
          length(human_genes), " human homologs.")
  return(human_genes)
}


# =============================================================================
# V. EXPRESSION VISUALISATION ALONG PSEUDOTIME
# =============================================================================

#' Plot raw expression of a gene along pseudotime for EPI / PE / TE lineages
#'
#' Produces a three-panel cowplot figure (one panel per fate).
#'
#' @param expr_log     data.frame (genes x cells). Log-normalised expression.
#' @param gene         Character. Gene name.
#' @param annot_epi    data.frame. Sample annotation for EPI lineage cells.
#'                     Must contain 'Pseudotime' and 'BranchesVal' columns.
#' @param annot_pe     data.frame. Sample annotation for PE lineage cells.
#' @param annot_te     data.frame. Sample annotation for TE lineage cells.
#' @param colours      Named list. Colour vectors for each fate.
#'                     Default: list(EPI = c("black","orange","red"), ...)
#' @param ylim_max     Numeric. Upper y-axis limit. Default: 20.
#'
#' @return Invisible cowplot object. Also prints to active device.
plot_expr_gene <- function(expr_log,
                            gene,
                            annot_epi,
                            annot_pe,
                            annot_te,
                            colours  = list(
                              EPI = c("black", "orange", "red"),
                              PE  = c("black", "orange", "green"),
                              TE  = c("black", "cyan",   "blue")
                            ),
                            ylim_max = 20) {
  .load_pkg("ggplot2")
  .load_pkg("cowplot")

  .make_panel <- function(annot, fate, col_vals) {
    cells     <- rownames(annot)
    cells_ok  <- cells[cells %in% colnames(expr_log)]
    expr_v    <- as.numeric(expr_log[gene, cells_ok])
    df        <- data.frame(
      Pseudotime  = annot[cells_ok, "Pseudotime"],
      Expression  = expr_v,
      BranchesVal = as.factor(as.character(annot[cells_ok, "BranchesVal"]))
    )
    ggplot(df, aes(x = Pseudotime, y = Expression, colour = BranchesVal)) +
      geom_point(shape = 18) +
      ylim(0, ylim_max) +
      geom_smooth(method = "loess", se = TRUE,
                  linetype = "solid", colour = tail(col_vals, 1),
                  fullrange = TRUE, formula = y ~ x) +
      scale_colour_manual(values = col_vals) +
      labs(title    = fate,
           x        = "Pseudotime",
           y        = if (fate == "EPI") paste("Log2 Expr -", gene) else "Log2 Expr",
           colour   = "Lineage") +
      theme(legend.position  = "none",
            axis.title.x     = element_blank(),
            axis.title.y     = if (fate != "EPI") element_blank() else NULL)
  }

  p_epi <- .make_panel(annot_epi, "EPI", colours$EPI)
  p_pe  <- .make_panel(annot_pe,  "PE",  colours$PE)
  p_te  <- .make_panel(annot_te,  "TE",  colours$TE)

  p <- ggdraw() +
    draw_plot(p_epi, 0,    0, 0.33, 1) +
    draw_plot(p_pe,  0.33, 0, 0.33, 1) +
    draw_plot(p_te,  0.66, 0, 0.33, 1)

  print(p)
  invisible(p)
}


#' Plot cluster mean expression along pseudotime for EPI / PE / TE lineages
#'
#' Same layout as plot_expr_gene() but for a pre-computed mean expression vector.
#'
#' @param cluster_mean_expr data.frame (1 row x cells) or numeric vector.
#' @param annot_epi         data.frame. EPI sample annotation.
#' @param annot_pe          data.frame. PE sample annotation.
#' @param annot_te          data.frame. TE sample annotation.
#' @param colours           Named list. Colour vectors (same structure as plot_expr_gene).
#' @param ylim_max          Numeric. Upper y-axis limit. Default: 20.
#'
#' @return Invisible cowplot object.
plot_cluster_expr <- function(cluster_mean_expr,
                               annot_epi,
                               annot_pe,
                               annot_te,
                               colours  = list(
                                 EPI = c("black", "orange", "red"),
                                 PE  = c("black", "orange", "green"),
                                 TE  = c("black", "cyan",   "blue")
                               ),
                               ylim_max = 20) {
  .load_pkg("ggplot2")
  .load_pkg("cowplot")

  .make_cluster_panel <- function(annot, fate, col_vals) {
    cells    <- rownames(annot)
    cells_ok <- cells[cells %in% colnames(cluster_mean_expr)]
    expr_v   <- as.numeric(cluster_mean_expr[, cells_ok])
    df       <- data.frame(
      Pseudotime  = annot[cells_ok, "Pseudotime"],
      Expression  = expr_v,
      BranchesVal = as.factor(as.character(annot[cells_ok, "BranchesVal"]))
    )
    ggplot(df, aes(x = Pseudotime, y = Expression, colour = BranchesVal)) +
      geom_point(shape = 18) +
      ylim(0, ylim_max) +
      geom_smooth(method = "loess", se = TRUE,
                  linetype = "solid", colour = tail(col_vals, 1),
                  fullrange = TRUE, formula = y ~ x) +
      scale_colour_manual(values = col_vals) +
      labs(title  = fate,
           x      = "Pseudotime",
           y      = if (fate == "EPI") "Expr Genes Cluster" else "Expr Gene",
           colour = "Lineage") +
      theme(legend.position = "none",
            axis.title.x    = element_blank(),
            axis.title.y    = if (fate != "EPI") element_blank() else NULL)
  }

  p_epi <- .make_cluster_panel(annot_epi, "EPI", colours$EPI)
  p_pe  <- .make_cluster_panel(annot_pe,  "PE",  colours$PE)
  p_te  <- .make_cluster_panel(annot_te,  "TE",  colours$TE)

  p <- ggdraw() +
    draw_plot(p_epi, 0,    0, 0.33, 1) +
    draw_plot(p_pe,  0.33, 0, 0.33, 1) +
    draw_plot(p_te,  0.66, 0, 0.33, 1)

  print(p)
  invisible(p)
}
