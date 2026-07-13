# =============================================================================
# SmoothRegressionPipeline.R
# -----------------------------------------------------------------------------
# LOESS-based smooth regression pipeline for scRNA-seq pseudotime analysis.
# Fits expression curves along pseudotime trajectories per fate lineage,
# aggregates shared segments, and exports results.
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin
# Project : scRNA-seq Human Pre-implantation Embryo — Lineage Specification
# Paper   : Meistermann D. et al., Cell Stem Cell, 2021
#           DOI: 10.1016/j.stem.2021.04.027
# Date    : 2018 (MSc M2 internship, CRTI UMR 1064, Nantes Universite)
# =============================================================================


# =============================================================================
# I. UTILITY FUNCTIONS
# =============================================================================

#' Create a directory if it does not already exist
#' @param path Character. Directory path.
create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created directory: ", path)
  }
}


#' Read a tab-separated file
#' @param filepath Character. Path to .tsv file.
#' @param as_character Logical. Read all columns as character.
#' @return data.frame with row names.
read_tsv <- function(filepath, as_character = FALSE) {
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  read.table(
    file       = filepath,
    sep        = "\t",
    header     = TRUE,
    row.names  = 1,
    colClasses = if (as_character) "character" else NA,
    quote      = ""
  )
}

#' Write a data frame to a tab-separated file
#' @param x data.frame.
#' @param filepath Character. Output path.
#' @param row_header Character. Row names column header label.
#' @param row.names Logical.
#' @param col.names Logical.
write_tsv <- function(x, filepath = "output.tsv",
                      row_header = "Name",
                      row.names  = TRUE,
                      col.names  = TRUE) {
  if (row.names && col.names) {
    # Prepend row-name header without triggering append warning
    writeLines(paste0(row_header, "\t", paste(colnames(x), collapse = "\t")),
               con = filepath)
    suppressWarnings(
      write.table(x, file = filepath, sep = "\t",
                  row.names = TRUE, col.names = FALSE,
                  quote = FALSE, append = TRUE)
    )
  } else {
    write.table(x, file = filepath, sep = "\t",
                row.names = row.names, col.names = col.names, quote = FALSE)
  }
  invisible(NULL)
}


#' Validate and coerce sample annotation
#' Ensures Lineage is a factor and Pseudotime is numeric.
#' @param sample_annot data.frame with Lineage and Pseudotime columns.
#' @return Validated data.frame.
validate_sample_annot <- function(sample_annot) {
  required <- c("Lineage", "Pseudotime")
  missing  <- setdiff(required, colnames(sample_annot))
  if (length(missing) > 0)
    stop("Missing columns in sample annotation: ", paste(missing, collapse = ", "))
  if (!is.factor(sample_annot$Lineage))
    sample_annot$Lineage <- as.factor(sample_annot$Lineage)
  if (!is.numeric(sample_annot$Pseudotime))
    sample_annot$Pseudotime <- as.numeric(as.character(sample_annot$Pseudotime))
  return(sample_annot)
}


# =============================================================================
# II. MAIN PIPELINE
# =============================================================================

#' LOESS Smooth Regression Pipeline for Pseudotime Expression
#'
#' Fits per-gene LOESS curves along pseudotime for each fate lineage,
#' aggregates shared segments by a summary function (e.g. median),
#' and exports expression matrices + sample annotations.
#'
#' @param expr         data.frame (genes x cells). Normalised expression.
#' @param sample_annot data.frame. Must contain 'Lineage' (factor) and
#'                     'Pseudotime' (numeric) columns.
#' @param tree_model   data.frame. Lineage tree with 'Lineages' (comma-
#'                     separated segment names) and 'leef' (leaf flags).
#' @param n            Integer. Pseudotime points per segment. Default: 100.
#' @param span         Numeric. LOESS bandwidth (0-1). Default: 0.75.
#' @param agg_fun      Function. Aggregation for shared segments. Default: median.
#' @param error_type   Logical. TRUE = standard error; FALSE = SD. Default: TRUE.
#' @param output_dir   Character. Root output directory. Default: "results".
#'
#' @return Invisible NULL. All output written to output_dir/.
#'
#' @examples
#' \dontrun{
#' smooth_pipeline(expr, sample_annot, tree_model)
#' }
smooth_pipeline <- function(expr,
                             sample_annot,
                             tree_model,
                             n          = 100,
                             span       = 0.75,
                             agg_fun    = median,
                             error_type = TRUE,
                             output_dir = "results") {

  # ── Validate inputs ───────────────────────────────────────────────────────
  sample_annot <- validate_sample_annot(sample_annot)
  message("[1/5] Inputs validated.")

  # ── Parse fate model ──────────────────────────────────────────────────────
  fate_lineages        <- strsplit(as.character(tree_model$Lineages), ",")
  names(fate_lineages) <- rownames(tree_model)
  fate_lineages        <- lapply(fate_lineages, as.factor)

  # ── Select cells per fate ─────────────────────────────────────────────────
  samples_by_fate <- lapply(fate_lineages, function(lins) {
    sample_annot[sample_annot$Lineage %in% lins, ]
  })

  # ── Build pseudotime grids ────────────────────────────────────────────────
  message("[2/5] Building pseudotime grids...")
  new_samples <- lapply(names(samples_by_fate), function(fate) {
    s    <- samples_by_fate[[fate]]
    segs <- levels(droplevels(s$Lineage))
    pts  <- c()
    lins <- c()
    begin <- TRUE
    last  <- NULL
    for (seg in segs) {
      rng <- sample_annot$Pseudotime[sample_annot$Lineage == seg]
      if (begin) {
        pt_seq <- round(seq(min(rng), max(rng), length.out = n), 2)
        begin  <- FALSE
      } else {
        last_max <- max(sample_annot$Pseudotime[sample_annot$Lineage == last])
        pt_seq   <- round(seq(last_max + 0.01, max(rng) - 0.05, length.out = n), 2)
      }
      pts  <- c(pts, pt_seq)
      lins <- c(lins, rep(seg, n))
      last <- seg
    }
    df           <- data.frame(Pseudotime = pts, Lineage = lins)
    rownames(df) <- paste0(fate, "_pt_", seq_len(nrow(df)))
    df
  })
  names(new_samples) <- names(samples_by_fate)

  # ── Fit LOESS per fate ────────────────────────────────────────────────────
  message("[3/5] Fitting LOESS curves...")
  fit_gene <- function(x, pt_obs, pt_new, se = FALSE) {
    model <- loess(x ~ pt_obs, span = span, degree = 2)
    predict(model, newdata = as.numeric(as.character(pt_new)), se = se)
  }

  new_expr    <- list()
  new_expr_sd <- list()

  for (fate in names(samples_by_fate)) {
    message("  Fate: ", fate)
    cells  <- rownames(samples_by_fate[[fate]])
    mat    <- as.matrix(expr[, cells])
    pt_obs <- samples_by_fate[[fate]]$Pseudotime
    pt_new <- as.numeric(as.character(new_samples[[fate]]$Pseudotime))
    cols   <- rownames(new_samples[[fate]])

    new_expr[[fate]] <- as.data.frame(
      t(apply(mat, 1, function(x) fit_gene(x, pt_obs, pt_new, se = FALSE)))
    )
    new_expr_sd[[fate]] <- as.data.frame(
      t(apply(mat, 1, function(x) {
        res <- fit_gene(x, pt_obs, pt_new, se = TRUE)
        if (error_type) res$se.fit else res$se.fit
      }))
    )
    colnames(new_expr[[fate]])    <- cols
    colnames(new_expr_sd[[fate]]) <- cols
    rownames(new_expr[[fate]])    <- rownames(expr)
    rownames(new_expr_sd[[fate]]) <- rownames(expr)
  }

  # ── Export per-fate ───────────────────────────────────────────────────────
  message("[4/5] Exporting per-fate results...")
  create_dir(output_dir)
  for (fate in names(samples_by_fate)) {
    d <- file.path(output_dir, fate)
    create_dir(d)
    write_tsv(new_expr[[fate]],    file.path(d, "ExpressionSample.tsv"))
    write_tsv(new_expr_sd[[fate]], file.path(d, "ExpressionSdSample.tsv"))
    write_tsv(new_samples[[fate]], file.path(d, "SampleAnnotation.tsv"))
  }

  # ── Aggregate shared segments ─────────────────────────────────────────────
  message("[5/5] Aggregating shared segments...")
  leaf_lins <- tree_model$leef

  unique_sample <- do.call(rbind, lapply(new_samples, function(s)
    s[s$Lineage %in% leaf_lins, ]))
  unique_expr <- do.call(cbind, lapply(names(new_expr), function(fate) {
    cols <- rownames(new_samples[[fate]][new_samples[[fate]]$Lineage %in% leaf_lins, ])
    new_expr[[fate]][, cols, drop = FALSE]
  }))
  unique_expr_sd <- do.call(cbind, lapply(names(new_expr_sd), function(fate) {
    cols <- rownames(new_samples[[fate]][new_samples[[fate]]$Lineage %in% leaf_lins, ])
    new_expr_sd[[fate]][, cols, drop = FALSE]
  }))

  shared_lins   <- setdiff(levels(sample_annot$Lineage), leaf_lins)
  common_sample <- do.call(rbind, list())
  common_expr   <- data.frame(row.names = rownames(expr))
  common_sd     <- data.frame(row.names = rownames(expr))
  common_fit_sd <- data.frame(row.names = rownames(expr))

  for (lin in shared_lins) {
    fates_w <- names(fate_lineages)[sapply(fate_lineages, function(f) lin %in% f)]
    seg_e   <- lapply(fates_w, function(fate) {
      cols <- rownames(new_samples[[fate]][new_samples[[fate]]$Lineage == lin, ])
      new_expr[[fate]][, cols, drop = FALSE]
    })
    seg_sd <- lapply(fates_w, function(fate) {
      cols <- rownames(new_samples[[fate]][new_samples[[fate]]$Lineage == lin, ])
      new_expr_sd[[fate]][, cols, drop = FALSE]
    })
    ref_cols   <- colnames(seg_e[[1]])
    ref_sample <- new_samples[[fates_w[1]]][new_samples[[fates_w[1]]]$Lineage == lin, ]
    common_sample <- rbind(common_sample, ref_sample)

    agg_e  <- data.frame(row.names = rownames(expr))
    agg_sd <- data.frame(row.names = rownames(expr))
    agg_fs <- data.frame(row.names = rownames(expr))
    for (i in seq_len(n)) {
      vals <- do.call(cbind, lapply(seg_e,  `[[`, i))
      sds  <- do.call(cbind, lapply(seg_sd, `[[`, i))
      agg_e  <- cbind(agg_e,  apply(vals, 1, agg_fun))
      agg_sd <- cbind(agg_sd, apply(sds,  1, agg_fun))
      agg_fs <- cbind(agg_fs, apply(vals, 1, sd))
    }
    colnames(agg_e) <- colnames(agg_sd) <- colnames(agg_fs) <- ref_cols
    common_expr   <- cbind(common_expr,   agg_e)
    common_sd     <- cbind(common_sd,     agg_sd)
    common_fit_sd <- cbind(common_fit_sd, agg_fs)
  }

  unique_fit_sd      <- unique_expr
  unique_fit_sd[, ]  <- 0
  fusion_dir <- file.path(output_dir, "SegmentFusion")
  create_dir(fusion_dir)
  write_tsv(cbind(common_expr,   unique_expr),    file.path(fusion_dir, "ExpressionSample.tsv"))
  write_tsv(cbind(common_sd,     unique_expr_sd), file.path(fusion_dir, "ExpressionSdSample.tsv"))
  write_tsv(cbind(common_fit_sd, unique_fit_sd),  file.path(fusion_dir, "FittedSdSample.tsv"))
  write_tsv(rbind(common_sample, unique_sample),  file.path(fusion_dir, "SampleAnnotation.tsv"))

  message("Done. Results in: ", output_dir)
  invisible(NULL)
}
