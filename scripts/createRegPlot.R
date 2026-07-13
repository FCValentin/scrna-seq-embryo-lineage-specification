# =============================================================================
# createRegPlot.R
# -----------------------------------------------------------------------------
# Main script: runs the LOESS smooth regression pipeline on scRNA-seq data
# and plots fitted expression curves along pseudotime for marker genes.
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin
# Project : scRNA-seq Human Pre-implantation Embryo — Lineage Specification
# Paper   : Meistermann D. et al., Cell Stem Cell, 2021
#           DOI: 10.1016/j.stem.2021.04.027
# Date    : 2018 (MSc M2 internship, CRTI UMR 1064, Nantes Universite)
#
# Usage:
#   Rscript createRegPlot.R
#   Or source() interactively in RStudio.
#
# Input (in data/):
#   sampleAnnot.tsv     : cell x metadata (requires Lineage + Pseudotime)
#   exprDat.norm.tsv    : gene x cell normalised expression matrix
#   SegmentsPosfai.tsv  : lineage tree model (Lineages + leef columns)
#
# Output:
#   results/            : per-fate and fused expression matrices (.tsv)
#   figs/               : expression curve plots per marker gene (.svg)
# =============================================================================


# =============================================================================
# DEPENDENCIES
# =============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
library(ggplot2)

source("SmoothRegressionPipeline.R")

# =============================================================================
# PARAMETERS — edit here
# =============================================================================

DATA_DIR    <- "data"
RESULTS_DIR <- "results"
FIGS_DIR    <- "figs"

N_POINTS    <- 100      # pseudotime interpolation points per segment
SPAN        <- 0.75     # LOESS smoothing bandwidth (0-1)
AGG_FUN     <- median   # aggregation for shared segments: median, mean...
ERROR_TYPE  <- TRUE     # TRUE = standard error  /  FALSE = SD

MARKER_GENES <- c(
  "Nanog", "Sox2", "Pou5f1", "Klf2",
  "Gata6", "Sox17", "Pgdfra", "Gata4",
  "Cdx2",  "Gata2", "Gata3", "Klf6"
)

# Named colour vector — must match Lineage factor levels in sampleAnnot
LINEAGE_COLOURS <- c(
  "Epiblast"           = "red",
  "Primitive endoderm" = "green",
  "Trophectoderm"      = "blue"
)

EXPORT_SVG <- FALSE   # set TRUE to save .svg files to FIGS_DIR


# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading data...")
sample_annot <- read_tsv(file.path(DATA_DIR, "sampleAnnot.tsv"))
expr         <- read_tsv(file.path(DATA_DIR, "exprDat.norm.tsv"))
tree_model   <- read_tsv(file.path(DATA_DIR, "SegmentsPosfai.tsv"))
message(sprintf("Loaded: %d genes x %d cells | %d annotated cells",
                nrow(expr), ncol(expr), nrow(sample_annot)))


# =============================================================================
# RUN PIPELINE
# =============================================================================

smooth_pipeline(
  expr        = expr,
  sample_annot = sample_annot,
  tree_model  = tree_model,
  n           = N_POINTS,
  span        = SPAN,
  agg_fun     = AGG_FUN,
  error_type  = ERROR_TYPE,
  output_dir  = RESULTS_DIR
)
# Output structure:
#   results/<fate>/ExpressionSample.tsv
#   results/<fate>/ExpressionSdSample.tsv
#   results/<fate>/SampleAnnotation.tsv
#   results/SegmentFusion/  (fused model across all fates)


# =============================================================================
# LOAD FUSED RESULTS
# =============================================================================

message("Loading fused results...")
fusion_dir   <- file.path(RESULTS_DIR, "SegmentFusion")
sample_final <- read_tsv(file.path(fusion_dir, "SampleAnnotation.tsv"))
expr_final   <- read_tsv(file.path(fusion_dir, "ExpressionSample.tsv"))
expr_sd      <- read_tsv(file.path(fusion_dir, "ExpressionSdSample.tsv"))

# Clip negative fitted values to zero
expr_final[expr_final < 0] <- 0
sd_min <- expr_final - expr_sd
sd_max <- expr_final + expr_sd


# =============================================================================
# PLOT FUNCTION
# =============================================================================

#' Plot fitted LOESS expression curve for one gene
#'
#' @param gene         Character. Gene name (rowname in expr_final).
#' @param expr_mat     data.frame. Fitted expression matrix (genes x pts).
#' @param sd_lo        data.frame. Lower ribbon bound.
#' @param sd_hi        data.frame. Upper ribbon bound.
#' @param annot        data.frame. SampleAnnotation with Pseudotime + Lineage.
#' @param colours      Named character vector. Lineage -> colour.
#' @param export_svg   Logical. Save to FIGS_DIR as .svg.
#' @param figs_dir     Character. Output directory for SVG files.
#' @return ggplot object (invisible).
plot_expression_curve <- function(gene,
                                   expr_mat,
                                   sd_lo,
                                   sd_hi,
                                   annot,
                                   colours    = LINEAGE_COLOURS,
                                   export_svg = EXPORT_SVG,
                                   figs_dir   = FIGS_DIR) {
  if (!gene %in% rownames(expr_mat)) {
    warning("Gene not found: ", gene)
    return(invisible(NULL))
  }

  pt      <- as.numeric(as.character(annot$Pseudotime))
  lineage <- annot$Lineage
  expr_v  <- as.numeric(unlist(expr_mat[gene, ]))

  ribbon_df <- data.frame(
    pt      = pt,
    Min     = pmax(as.numeric(unlist(sd_lo[gene, ])), 0),
    Max     = as.numeric(unlist(sd_hi[gene, ])),
    Lineage = lineage
  )

  p <- ggplot(
    data.frame(Pseudotime = pt, Expression = expr_v, Lineage = lineage),
    aes(x = Pseudotime, y = Expression, colour = Lineage)
  ) +
    geom_ribbon(
      data        = ribbon_df,
      aes(x = pt, ymin = Min, ymax = Max, fill = Lineage),
      alpha       = 0.3,
      inherit.aes = FALSE,
      colour      = NA
    ) +
    geom_line(linewidth = 1.5) +
    scale_colour_manual(values = colours) +
    scale_fill_manual(values = colours) +
    guides(colour = "none", fill = "none") +
    labs(x = "Pseudotime", y = paste0(gene, " — fitted expression")) +
    theme(
      panel.background = element_rect(fill = "#EEEEEE", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background  = element_rect(fill = "transparent", colour = NA)
    )

  if (export_svg) {
    if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)
    svg_path <- file.path(figs_dir, paste0(gene, ".svg"))
    svg(filename = svg_path, width = 5, height = 3, bg = "transparent")
    print(p)
    dev.off()
    message("Saved: ", svg_path)
  } else {
    print(p)
  }
  invisible(p)
}


# =============================================================================
# PLOT MARKER GENES
# =============================================================================

message("Plotting ", length(MARKER_GENES), " marker genes...")
for (gene in MARKER_GENES) {
  plot_expression_curve(
    gene      = gene,
    expr_mat  = expr_final,
    sd_lo     = sd_min,
    sd_hi     = sd_max,
    annot     = sample_final
  )
}
message("Done.")
