#' Extract z-score matrix for barcode heatmap
#'
#' @description
#' Performs all data-loading and GCTX-extraction steps needed for
#' \code{\link{plot_signature_direction_tile_barcode}} without producing any
#' plot. The returned matrix can be inspected, filtered, exported, or passed
#' directly to \code{plot_signature_direction_tile_barcode(precomputed = )}
#' to avoid re-reading large files on repeated calls.
#'
#' @param results_df Data frame containing at least perturbation id and score
#'   columns.
#' @param signature_file Path to signature file with gene and log2FC columns.
#' @param reference_df Optional reference matrix with genes as rows and
#'   perturbation ids as columns. If supplied, z-scores are taken directly from
#'   this object and no GCTX file is required.
#' @param gctx_file Optional path to GCTX expression file.
#' @param geneinfo_file Optional path to geneinfo file.
#' @param siginfo_file Optional path to siginfo file for readable labels.
#' @param data_dir Optional directory containing CMap files.
#' @param selected_drug Optional drug identifier to filter \code{results_df}.
#' @param selected_drug_col Column used with \code{selected_drug}.
#' @param pert_id_col Perturbation id column (default: \code{"compound"}).
#' @param score_col Score column (default: \code{"Score"}).
#' @param max_genes Maximum number of signature genes (default: 100).
#' @param max_perts Maximum number of perturbations (default: 60).
#' @param output_zscores Optional TSV path to save the matrix. \code{NULL}
#'   to skip.
#' @param verbose Logical; print progress messages.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{z_plot}}{Numeric matrix — rows = perturbations (labelled
#'       \code{cmap_name | dose | time | cell}), cols = genes (down→up order).}
#'     \item{\code{ordered_genes}}{Character vector of gene symbols.}
#'     \item{\code{logfc_map}}{Named numeric vector of signature log2FC values.}
#'     \item{\code{sig_ids}}{Character vector of selected \code{sig_id}s.}
#'     \item{\code{sig_labels}}{Character vector of human-readable labels.}
#'   }
#'
#' @examples
#' is.function(extract_signature_zscores)
#'
#' \donttest{
#' # Requires the full CMap GCTX file (downloaded from clue.io)
#' sig_file <- system.file("extdata", "example_signature.txt",
#'                         package = "CONCERTDR")
#' ref_file <- system.file("extdata", "example_reference_df.csv",
#'                         package = "CONCERTDR")
#' ref_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
#' ref_df$gene_symbol <- rownames(ref_df)
#' zmat <- extract_signature_zscores(
#'   results_df     = data.frame(sig_id = "DEMO001", Score = -0.72),
#'   signature_file = sig_file,
#'   reference_df   = ref_df,
#'   pert_id_col    = "sig_id"
#' )
#' }
#'
#' @export
extract_signature_zscores <- function(results_df,
                                      signature_file,
                                      reference_df       = NULL,
                                      gctx_file         = NULL,
                                      geneinfo_file     = NULL,
                                      siginfo_file      = NULL,
                                      data_dir          = getOption("CONCERTDR.data_dir", NULL),
                                      selected_drug     = NULL,
                                      selected_drug_col = NULL,
                                      pert_id_col       = "compound",
                                      score_col         = "Score",
                                      max_genes         = 100,
                                      max_perts         = 60,
                                      output_zscores    = NULL,
                                      verbose           = TRUE) {
  if (!is.data.frame(results_df)) stop("results_df must be a data.frame")
  if (!pert_id_col %in% names(results_df)) stop("results_df must contain pert_id_col: ", pert_id_col)
  if (!score_col %in% names(results_df)) stop("results_df must contain score_col: ", score_col)
  if (!file.exists(signature_file)) stop("signature_file not found: ", signature_file)
  
  first_existing <- function(candidates) {
    candidates <- as.character(candidates)
    candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
    if (length(candidates) == 0) return(NULL)
    ok <- candidates[file.exists(candidates)]
    if (length(ok) == 0) return(NULL)
    ok[1]
  }

  use_reference_df <- !is.null(reference_df)

  if (use_reference_df) {
    if (is.data.frame(reference_df)) {
      if ("gene_symbol" %in% names(reference_df)) {
        reference_df$gene_symbol <- toupper(as.character(reference_df$gene_symbol))
        rownames(reference_df) <- reference_df$gene_symbol
        reference_df$gene_symbol <- NULL
      }
      reference_mat <- as.matrix(reference_df)
    } else if (is.matrix(reference_df)) {
      reference_mat <- reference_df
    } else {
      stop("reference_df must be a data.frame or matrix")
    }

    if (is.null(rownames(reference_mat)) || is.null(colnames(reference_mat))) {
      stop("reference_df must have gene identifiers as row names and perturbation ids as column names")
    }

    rownames(reference_mat) <- toupper(as.character(rownames(reference_mat)))
    reference_mat <- apply(reference_mat, 2, as.numeric)
    rownames(reference_mat) <- toupper(as.character(rownames(reference_df)))
    colnames(reference_mat) <- colnames(reference_df)
  }

  gctx_default_name <- "level5_beta_all_n1201944x12328.gctx"

  gctx_file <- if (!use_reference_df) {
    first_existing(c(
      gctx_file,
      getOption("CONCERTDR.gctx_file", NULL),
      Sys.getenv("CONCERTDR_GCTX_FILE", unset = ""),
      if (!is.null(data_dir) && nzchar(data_dir)) file.path(data_dir, gctx_default_name) else NULL,
      gctx_default_name
    ))
  } else {
    first_existing(c(
      gctx_file,
      getOption("CONCERTDR.gctx_file", NULL),
      Sys.getenv("CONCERTDR_GCTX_FILE", unset = "")
    ))
  }

  if (is.null(gctx_file) && !use_reference_df) {
    stop("Could not resolve gctx_file. Provide gctx_file explicitly, set options(CONCERTDR.gctx_file='...'), or supply reference_df")
  }

  inferred_data_dir <- if (!is.null(gctx_file)) {
    if (!is.null(data_dir) && nzchar(data_dir)) data_dir else dirname(gctx_file)
  } else {
    data_dir
  }

  gene_map <- NULL
  if (!use_reference_df) {
    geneinfo_file <- first_existing(c(
      geneinfo_file,
      getOption("CONCERTDR.geneinfo_file", NULL),
      Sys.getenv("CONCERTDR_GENEINFO_FILE", unset = ""),
      if (!is.null(inferred_data_dir) && nzchar(inferred_data_dir)) file.path(inferred_data_dir, "geneinfo_beta.txt") else NULL,
      "geneinfo_beta.txt"
    ))
    if (is.null(geneinfo_file)) {
      stop("Could not resolve geneinfo_file. Provide geneinfo_file explicitly, or set ",
           "options(CONCERTDR.geneinfo_file='...') or data_dir containing geneinfo_beta.txt")
    }
  }
  
  siginfo_file <- first_existing(c(
    siginfo_file,
    getOption("CONCERTDR.siginfo_file", NULL),
    Sys.getenv("CONCERTDR_SIGINFO_FILE", unset = ""),
    if (!is.null(inferred_data_dir) && nzchar(inferred_data_dir)) file.path(inferred_data_dir, "siginfo_beta.txt") else NULL,
    "siginfo_beta.txt"
  ))
  
  if (verbose) {
    if (use_reference_df) {
      message("Using in-memory reference_df with ", nrow(reference_mat),
              " genes and ", ncol(reference_mat), " perturbations")
    } else {
      message("Using files:")
      message(" - gctx_file: ", gctx_file)
      message(" - geneinfo_file: ", geneinfo_file)
    }
    if (!is.null(siginfo_file)) message(" - siginfo_file: ", siginfo_file)
  }

  if (!use_reference_df) {
    geneinfo <- if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::fread(geneinfo_file, data.table = FALSE)
    } else {
      utils::read.table(geneinfo_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                        quote = "", comment.char = "", fill = TRUE)
    }
    names(geneinfo) <- trimws(names(geneinfo))
    if (!all(c("gene_symbol", "gene_id") %in% names(geneinfo))) {
      stop("geneinfo_file must contain columns: gene_symbol, gene_id")
    }
    geneinfo$gene_symbol <- toupper(as.character(geneinfo$gene_symbol))
    geneinfo <- geneinfo[!is.na(geneinfo$gene_symbol) & !is.na(geneinfo$gene_id), c("gene_symbol", "gene_id")]
    gene_map <- as.character(geneinfo$gene_id)
    names(gene_map) <- geneinfo$gene_symbol
  }

  # signature genes and order (down then up)
  sig <- if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(signature_file, data.table = FALSE)
  } else {
    utils::read.table(signature_file, header = TRUE, sep = "", stringsAsFactors = FALSE,
                      quote = "", comment.char = "", fill = TRUE)
  }
  names(sig) <- trimws(names(sig))
  gene_col   <- if ("Gene"   %in% names(sig)) "Gene"   else names(sig)[1]
  log2fc_col <- if ("log2FC" %in% names(sig)) "log2FC" else names(sig)[2]
  sig[[gene_col]]   <- toupper(as.character(sig[[gene_col]]))
  sig[[log2fc_col]] <- suppressWarnings(as.numeric(sig[[log2fc_col]]))
  sig <- sig[!is.na(sig[[gene_col]]) & nzchar(sig[[gene_col]]) & !is.na(sig[[log2fc_col]]), , drop = FALSE]
  if (!is.null(max_genes)) sig <- utils::head(sig, max_genes)
  if (use_reference_df) {
    sig <- sig[sig[[gene_col]] %in% rownames(reference_mat), , drop = FALSE]
    if (nrow(sig) == 0) stop("No signature genes matched the row names of reference_df")
  } else {
    sig <- sig[sig[[gene_col]] %in% names(gene_map), , drop = FALSE]
    if (nrow(sig) == 0) stop("No signature genes mapped to geneinfo gene_id")
  }

  down <- sig[sig[[log2fc_col]] < 0, , drop = FALSE]
  down <- down[order(down[[log2fc_col]], decreasing = FALSE), , drop = FALSE]
  up   <- sig[sig[[log2fc_col]] > 0, , drop = FALSE]
  up   <- up[order(up[[log2fc_col]], decreasing = TRUE), , drop = FALSE]
  
  ordered_genes <- c(as.character(down[[gene_col]]), as.character(up[[gene_col]]))
  ordered_genes <- ordered_genes[nzchar(ordered_genes)]
  if (length(ordered_genes) == 0) stop("No ordered genes available after down/up split")
  ordered_ids <- if (use_reference_df) ordered_genes else unname(gene_map[ordered_genes])
  logfc_map   <- sig[[log2fc_col]]
  names(logfc_map) <- sig[[gene_col]]
  
  # select perturbations by score
  tech <- results_df
  tech[[score_col]] <- suppressWarnings(as.numeric(tech[[score_col]]))
  tech <- tech[!is.na(tech[[score_col]]) & !is.na(tech[[pert_id_col]]), , drop = FALSE]
  
  if (!is.null(selected_drug)) {
    drug_col <- selected_drug_col
    if (is.null(drug_col)) {
      candidates <- intersect(c("display_name", "pert_id", "cmap_name", "pert_name"), names(tech))
      if (length(candidates) == 0) stop("selected_drug given but no suitable selected_drug_col found")
      drug_col <- candidates[1]
    }
    keep <- toupper(as.character(tech[[drug_col]])) == toupper(selected_drug)
    tech <- tech[keep, , drop = FALSE]
    if (verbose) message("Filtered rows for selected_drug using ", drug_col, ": ", nrow(tech))
  }
  
  tech    <- tech[order(tech[[score_col]], decreasing = FALSE), , drop = FALSE]
  tech    <- utils::head(tech, max_perts)
  sig_ids <- as.character(tech[[pert_id_col]])
  sig_ids <- sig_ids[!is.na(sig_ids) & nzchar(sig_ids)]
  sig_ids <- unique(sig_ids)
  if (use_reference_df) {
    sig_ids <- intersect(sig_ids, colnames(reference_mat))
  }
  if (length(sig_ids) == 0) stop("No perturbation ids selected")
  if (verbose) message("Perturbations selected: ", length(sig_ids))
  
  # optional labels from siginfo
  sig_labels <- sig_ids
  if (!is.null(siginfo_file)) {
    si <- if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::fread(siginfo_file, data.table = FALSE)
    } else {
      utils::read.table(siginfo_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                        quote = "", comment.char = "", fill = TRUE)
    }
    names(si) <- trimws(names(si))
    req  <- c("sig_id", "cmap_name", "pert_idose", "pert_itime", "cell_iname")
    have <- intersect(req, names(si))
    if ("sig_id" %in% have) {
      si <- si[, have, drop = FALSE]
      si <- si[!duplicated(si$sig_id), , drop = FALSE]
      rownames(si) <- as.character(si$sig_id)
      make_label <- function(sid) {
        if (!sid %in% rownames(si)) return(sid)
        row  <- si[sid, , drop = FALSE]
        vals <- c()
        for (nm in c("cmap_name", "pert_idose", "pert_itime", "cell_iname")) {
          if (nm %in% names(row)) {
            v <- as.character(row[[nm]])[1]
            if (!is.na(v) && nzchar(v)) vals <- c(vals, v)
          }
        }
        if (length(vals) == 0) sid else paste(vals, collapse = " | ")
      }
      sig_labels <- vapply(sig_ids, make_label, character(1))
    }
  }

  if (use_reference_df) {
    z_mat <- reference_mat[ordered_ids, sig_ids, drop = FALSE]
    rownames(z_mat) <- ordered_genes
  } else {
    z_mat <- fast_parse_gctx(fname = gctx_file, rid = as.character(ordered_ids), cid = sig_ids)
    z_mat <- z_mat[as.character(ordered_ids), sig_ids, drop = FALSE]
    rownames(z_mat) <- ordered_genes
  }
  z_plot          <- t(z_mat)          # rows = perturbations, cols = genes
  rownames(z_plot) <- sig_labels
  
  if (nrow(z_plot) == 0 || ncol(z_plot) == 0) {
    stop("No data available after GCTX extraction")
  }
  
  if (!is.null(output_zscores)) {
    utils::write.table(as.data.frame(z_plot), file = output_zscores, sep = "\t",
                       quote = FALSE, row.names = TRUE, col.names = NA)
    if (verbose) message("Saved z-score matrix: ", output_zscores)
  }
  
  list(
    z_plot        = z_plot,
    ordered_genes = ordered_genes,
    logfc_map     = logfc_map,
    sig_ids       = sig_ids,
    sig_labels    = sig_labels
  )
}


#' Plot 2D Barcode Heatmap from GCTX (Single Drug Context)
#'
#' @description
#' Build a 2D barcode heatmap using z-scores extracted from a GCTX file.
#' Genes are ordered by signature direction (down then up), and perturbations
#' are chosen as top \code{max_perts} rows by best (lowest) \code{score_col}
#' from the provided technical results table.
#'
#' @param results_df Data frame (typically \code{tech_view_all}) containing at
#'   least perturbation id and score columns.
#' @param signature_file Path to signature file with gene and log2FC columns.
#' @param reference_df Optional reference matrix with genes as rows and
#'   perturbation ids as columns. If supplied, no GCTX file is needed and the
#'   heatmap is built directly from this matrix.
#' @param gctx_file Optional path to GCTX expression file. If NULL, the function
#'   tries (in order): option \code{CONCERTDR.gctx_file}, env var
#'   \code{CONCERTDR_GCTX_FILE}, \code{data_dir/level5_beta_all_n1201944x12328.gctx},
#'   then a file with that name in the working directory.
#' @param geneinfo_file Optional path to geneinfo file with \code{gene_symbol}
#'   and \code{gene_id} columns. If NULL, the function tries option
#'   \code{CONCERTDR.geneinfo_file}, env var \code{CONCERTDR_GENEINFO_FILE},
#'   \code{data_dir/geneinfo_beta.txt}, then \code{geneinfo_beta.txt} in the
#'   working directory.
#' @param siginfo_file Optional path to siginfo file for readable perturbation
#'   labels. If NULL, the function tries option \code{CONCERTDR.siginfo_file},
#'   env var \code{CONCERTDR_SIGINFO_FILE}, \code{data_dir/siginfo_beta.txt},
#'   then \code{siginfo_beta.txt} in the working directory.
#' @param data_dir Optional directory containing CMap files. Useful to avoid
#'   repeatedly passing full file paths.
#' @param selected_drug Optional drug identifier to filter \code{results_df}
#'   before selecting perturbations.
#' @param selected_drug_col Optional column used with \code{selected_drug}.
#'   If NULL, auto-detect from \code{display_name}, \code{pert_id},
#'   \code{cmap_name}, \code{pert_name}.
#' @param pert_id_col Perturbation id column in \code{results_df}
#'   (default: \code{"compound"}, i.e. sig_id in CONCERTDR results).
#' @param score_col Score column in \code{results_df} (default: \code{"Score"}).
#' @param max_genes Maximum number of signature genes to use (default: 100).
#' @param max_perts Maximum number of perturbations to show (default: 60).
#' @param cluster_rows Logical; cluster perturbation rows (default: TRUE).
#' @param cluster_method Clustering linkage method for row clustering
#'   (default: \code{"complete"}).
#' @param show_row_dendrogram Logical; whether to draw a row dendrogram when
#'   \code{cluster_rows = TRUE} (default: TRUE).
#' @param cluster_cols Logical; cluster gene columns (default: FALSE — genes
#'   are kept in signature direction order: down then up).
#' @param cluster_method_cols Clustering linkage method for column clustering
#'   (default: \code{"complete"}).
#' @param show_col_dendrogram Logical; whether to draw a column dendrogram
#'   when \code{cluster_cols = TRUE} (default: TRUE).
#' @param precomputed Optional; the output of \code{\link{extract_signature_zscores}}.
#'   If provided, skips all data-loading and GCTX-extraction steps — all
#'   data-related arguments (\code{results_df}, \code{signature_file},
#'   \code{gctx_file}, etc.) are ignored.
#' @param save_png Logical; save PNG output (default: FALSE).
#' @param output_png Output PNG path.
#' @param output_zscores Optional TSV path for exported z-score matrix.
#'   Default is \code{NULL}, so no file is written unless explicitly requested.
#' @param width,height Figure dimensions in inches for PNG. If \code{NULL}
#'   (default), dimensions are computed automatically from the number of genes
#'   and perturbations (~0.22 in/gene width, ~0.28 in/perturbation height).
#' @param dpi PNG resolution.
#' @param verbose Logical; print progress messages.
#'
#' @return Invisibly returns a list containing the plotted matrix,
#' selected perturbation ids, and file paths.
#'
#' @examples
#' is.function(plot_signature_direction_tile_barcode)
#'
#' \donttest{
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("circlize", quietly = TRUE)) {
#'   sig_file <- system.file("extdata", "example_signature.txt",
#'                           package = "CONCERTDR")
#'   plot_signature_direction_tile_barcode(
#'     results_df    = data.frame(sig_id = "DEMO001", Score = -0.72),
#'     signature_file = sig_file,
#'     reference_df  = read.csv(system.file("extdata", "example_reference_df.csv",
#'                                          package = "CONCERTDR"),
#'                              row.names = 1, check.names = FALSE),
#'     pert_id_col   = "sig_id"
#'   )
#' }
#' }
#'
#' @export
plot_signature_direction_tile_barcode <- function(results_df = NULL,
                                                  signature_file = NULL,
                                                  reference_df = NULL,
                                                  gctx_file = NULL,
                                                  geneinfo_file = NULL,
                                                  siginfo_file = NULL,
                                                  data_dir = getOption("CONCERTDR.data_dir", NULL),
                                                  selected_drug = NULL,
                                                  selected_drug_col = NULL,
                                                  pert_id_col = "compound",
                                                  score_col = "Score",
                                                  max_genes = 100,
                                                  max_perts = 60,
                                                  cluster_rows = TRUE,
                                                  cluster_method = "complete",
                                                  show_row_dendrogram = TRUE,
                                                  cluster_cols = FALSE,
                                                  cluster_method_cols = "complete",
                                                  show_col_dendrogram = TRUE,
                                                  precomputed         = NULL,
                                                  save_png = FALSE,
                                                  output_png = "barcode_heatmap.png",
                                                  output_zscores = NULL,
                                                  width = NULL,
                                                  height = NULL,
                                                  dpi = 150,
                                                  verbose = TRUE) {
  # ── data preparation ────────────────────────────────────────────────────────
  if (!is.null(precomputed)) {
    if (!is.list(precomputed) ||
        !all(c("z_plot", "ordered_genes", "logfc_map", "sig_ids") %in% names(precomputed))) {
      stop("precomputed must be the output of extract_signature_zscores()")
    }
    z_plot        <- precomputed$z_plot
    ordered_genes <- precomputed$ordered_genes
    logfc_map     <- precomputed$logfc_map
    sig_ids       <- precomputed$sig_ids
    if (verbose) message("Using precomputed matrix: ", nrow(z_plot),
                         " perturbations \u00d7 ", ncol(z_plot), " genes")
  } else {
    if (is.null(results_df))     stop("Provide results_df or precomputed")
    if (is.null(signature_file)) stop("Provide signature_file or precomputed")
    extracted <- extract_signature_zscores(
      results_df        = results_df,
      signature_file    = signature_file,
      reference_df      = reference_df,
      gctx_file         = gctx_file,
      geneinfo_file     = geneinfo_file,
      siginfo_file      = siginfo_file,
      data_dir          = data_dir,
      selected_drug     = selected_drug,
      selected_drug_col = selected_drug_col,
      pert_id_col       = pert_id_col,
      score_col         = score_col,
      max_genes         = max_genes,
      max_perts         = max_perts,
      output_zscores    = output_zscores,
      verbose           = verbose
    )
    z_plot        <- extracted$z_plot
    ordered_genes <- extracted$ordered_genes
    logfc_map     <- extracted$logfc_map
    sig_ids       <- extracted$sig_ids
  }
  
  if (nrow(z_plot) == 0 || ncol(z_plot) == 0) {
    stop("No data available for heatmap")
  }
  
  # ── ComplexHeatmap ─────────────────────────────────────────────────────────
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop(
      "Package 'ComplexHeatmap' is required. Install with:\n",
      "  BiocManager::install('ComplexHeatmap')"
    )
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Install with: install.packages('circlize')")
  }
  
  # z-score colour scale (symmetric around 0)
  zlim <- max(abs(z_plot), na.rm = TRUE)
  if (!is.finite(zlim) || zlim == 0) zlim <- 10
  col_fun <- circlize::colorRamp2(c(-zlim, 0, zlim), c("#3B4CC0", "#F7F7F7", "#B40426"))
  
  # signature log2FC colour strip — BrBG (teal → white → brown).
  # Teal/green = down-regulated (negative log2FC), brown/orange = up-regulated.
  # Completely different hue family from the coolwarm z-score palette so the
  # two scales are immediately distinguishable.
  logfc_vals    <- as.numeric(logfc_map[ordered_genes])
  lim           <- max(abs(logfc_vals), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  logfc_col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("#01665E", "#F5F5F5", "#8C510A"))
  
  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    "Signature log2FC" = logfc_vals,
    col = list("Signature log2FC" = logfc_col_fun),
    annotation_legend_param = list(
      "Signature log2FC" = list(
        title      = "Signature log2FC",
        title_gp   = grid::gpar(fontsize = 9),
        labels_gp  = grid::gpar(fontsize = 8)
      )
    ),
    show_annotation_name = TRUE,
    annotation_name_gp   = grid::gpar(fontsize = 9)
  )
  
  ttl <- if (is.null(selected_drug)) {
    paste0("Top ", nrow(z_plot), " perturbations by ", score_col)
  } else {
    paste0(selected_drug, " \u2013 Top ", nrow(z_plot), " perturbations by ", score_col)
  }
  
  ht <- ComplexHeatmap::Heatmap(
    z_plot,
    name                   = "z-score",
    col                    = col_fun,
    cluster_rows           = isTRUE(cluster_rows),
    clustering_method_rows = cluster_method,
    cluster_columns           = isTRUE(cluster_cols),
    clustering_method_columns = cluster_method_cols,
    show_row_dend          = isTRUE(show_row_dendrogram) && isTRUE(cluster_rows),
    show_column_dend       = isTRUE(show_col_dendrogram) && isTRUE(cluster_cols),
    top_annotation         = top_ann,
    show_row_names         = TRUE,
    show_column_names      = TRUE,
    row_names_gp           = grid::gpar(fontsize = 8),
    row_names_max_width    = grid::unit(7, "cm"),
    column_names_gp        = grid::gpar(fontsize = 8),
    column_names_rot       = 60,
    column_title           = ttl,
    column_title_gp        = grid::gpar(fontsize = 11, fontface = "bold"),
    row_title              = "Perturbations",
    row_title_gp           = grid::gpar(fontsize = 10),
    heatmap_legend_param   = list(
      title     = "z-score",
      title_gp  = grid::gpar(fontsize = 9),
      labels_gp = grid::gpar(fontsize = 8)
    ),
    use_raster    = TRUE,
    raster_quality = 2
  )
  
  # dynamic figure size: scale with data dimensions if not supplied
  n_cols <- ncol(z_plot)
  n_rows <- nrow(z_plot)
  auto_width  <- max(14, 4 + n_cols * 0.22 + 7)   # 4in base + ~0.22in/gene + 7in for labels/dendro/legend
  auto_height <- max(8,  2 + n_rows * 0.28)         # 2in base + ~0.28in/perturbation
  fig_width  <- if (!is.null(width))  width  else auto_width
  fig_height <- if (!is.null(height)) height else auto_height
  
  draw_fn <- function() {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side    = "right",
      annotation_legend_side = "right",
      padding                = grid::unit(c(5, 20, 8, 5), "mm")
    )
  }
  
  if (isTRUE(save_png)) {
    grDevices::png(filename = output_png, width = fig_width, height = fig_height,
                   units = "in", res = dpi)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_fn()
    if (verbose) message("Saved heatmap: ", output_png)
  } else {
    draw_fn()
  }
  
  invisible(list(
    z_plot        = z_plot,
    ordered_genes = ordered_genes,
    sig_ids       = sig_ids,
    output_png    = if (isTRUE(save_png)) output_png else NULL,
    output_zscores = output_zscores
  ))
}
