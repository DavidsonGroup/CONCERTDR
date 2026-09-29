#' Extract z-score matrix for barcode heatmap
#'
#' @description
#' Performs all data-loading and GCTX-extraction steps needed for
#' \code{\link{plot_signature_direction_tile_barcode}} without producing any
#' plot. The returned matrix can be inspected, filtered, exported, or passed
#' directly to \code{plot_signature_direction_tile_barcode(precomputed = )}
#' to avoid re-reading large files on repeated calls.
#'
#' @param direction Order in which perturbations are picked by score:
#'   \code{"reversal"} (default) lowest first, \code{"mimic"} highest first.
#'   No perturbation is excluded by the sign of its score. Precomputed
#'   matrices are plotted as supplied; choose the direction when extracting
#'   them.
#' @param results_df Data frame containing at least perturbation id and score
#'   columns.
#' @param signature_file Path to a signature file with gene and log2FC columns,
#'   or a data frame with \code{Gene} and \code{log2FC} columns.
#' @param reference_df Optional reference matrix with genes as rows and
#'   perturbation ids as columns. If supplied, z-scores are taken directly from
#'   this object and no GCTX file is required.
#' @param gctx_file Optional path to GCTX expression file.
#' @param geneinfo_file Optional path to geneinfo file.
#' @param siginfo_file Optional path to siginfo file for readable labels.
#' @param data_dir Optional directory containing CMap files.
#' @param selected_drug Optional drug identifier to filter \code{results_df}.
#' @param selected_drug_col Column used with \code{selected_drug}.
#' @param pert_id_col Perturbation id column (default: \code{"sig_id"}).
#' @param score_col Score column (default: \code{"Score"}).
#' @param max_genes Maximum number of signature genes to use (default: 100).
#'   When \code{split_direction = FALSE} (default) the first \code{max_genes}
#'   rows of the signature file are taken. When \code{split_direction = TRUE}
#'   the top \code{max_genes} up-regulated genes (by log2FC) and the top
#'   \code{max_genes} down-regulated genes (by |log2FC|) are selected
#'   independently, so the matrix may contain up to \code{2 * max_genes}
#'   gene columns in total.
#' @param max_perts Maximum number of perturbations (default: 60).
#' @param split_direction Logical; if \code{TRUE}, apply \code{max_genes}
#'   separately to the up-regulated (log2FC > 0) and down-regulated
#'   (log2FC < 0) gene sets rather than to the combined signature.
#'   The result contains up to \code{2 * max_genes} columns ordered
#'   down-regulated then up-regulated, matching the layout expected by
#'   \code{\link{plot_signature_direction_tile_barcode}} with
#'   \code{split_direction = TRUE}. Default: \code{FALSE}.
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
#'     \item{\code{sig_scores}}{Numeric vector of \code{score_col} values in
#'       the same order as \code{sig_ids} and the rows of \code{z_plot}.}
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
#' demo_sig <- colnames(ref_df)[1]
#' zmat <- extract_signature_zscores(
#'   results_df     = data.frame(sig_id = demo_sig, Score = -0.72),
#'   signature_file = sig_file,
#'   reference_df   = ref_df,
#'   pert_id_col    = "sig_id"
#' )
#' }
#'
#' @export
extract_signature_zscores <- function(results_df,
                                      signature_file,
                                      reference_df = NULL,
                                      gctx_file = NULL,
                                      geneinfo_file = NULL,
                                      siginfo_file = NULL,
                                      data_dir = getOption("CONCERTDR.data_dir", NULL),
                                      selected_drug = NULL,
                                      selected_drug_col = NULL,
                                      pert_id_col = "sig_id",
                                      score_col = "Score",
                                      max_genes = 100,
                                      max_perts = 60,
                                      split_direction = FALSE,
                                      output_zscores = NULL,
                                      verbose = TRUE,
                                      direction = c("reversal", "mimic")) {
  direction <- match.arg(direction)
  if (!is.data.frame(results_df)) stop("results_df must be a data.frame")
  if (!pert_id_col %in% names(results_df)) {
    fallback_cols <- c("sig_id", "compound")
    fallback_cols <- fallback_cols[fallback_cols %in% names(results_df)]
    if (length(fallback_cols) > 0) {
      if (verbose) {
        message("pert_id_col '", pert_id_col, "' not found; using '", fallback_cols[1], "' instead")
      }
      pert_id_col <- fallback_cols[1]
    } else {
      stop("results_df must contain pert_id_col: ", pert_id_col)
    }
  }
  if (!score_col %in% names(results_df)) stop("results_df must contain score_col: ", score_col)

  # Signature input: data.frame or file path
  if (!is.data.frame(signature_file)) {
    if (!(is.character(signature_file) && length(signature_file) == 1L)) {
      stop("signature_file must be either a file path (character) or a data.frame with 'Gene' and 'log2FC' columns.")
    }
    if (!file.exists(signature_file)) stop("signature_file not found: ", signature_file)
  }

  use_reference_df <- !is.null(reference_df)
  if (use_reference_df) reference_mat <- .reference_to_matrix(reference_df)

  # With reference_df no GCTX/geneinfo is needed, so only explicit sources count.
  gctx_default_name <- "level5_beta_all_n1201944x12328.gctx"
  gctx_file <- .resolve_cmap_file(
    gctx_file, "gctx",
    default = if (!use_reference_df) gctx_default_name,
    data_dir = data_dir
  )
  if (is.null(gctx_file) && !use_reference_df) {
    stop("Could not resolve gctx_file. Provide gctx_file explicitly, set options(CONCERTDR.gctx_file='...'), or supply reference_df")
  }

  has_data_dir <- !is.null(data_dir) && nzchar(data_dir)
  inferred_data_dir <- if (!has_data_dir && !is.null(gctx_file)) dirname(gctx_file) else data_dir

  gene_map <- NULL
  if (!use_reference_df) {
    geneinfo_file <- .resolve_cmap_file(geneinfo_file, "geneinfo", "geneinfo_beta.txt", inferred_data_dir)
    if (is.null(geneinfo_file)) {
      stop("Could not resolve geneinfo_file. Provide geneinfo_file explicitly, or set ",
           "options(CONCERTDR.geneinfo_file='...') or data_dir containing geneinfo_beta.txt")
    }
  }
  siginfo_file <- .resolve_cmap_file(siginfo_file, "siginfo", "siginfo_beta.txt", inferred_data_dir)

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

  if (!use_reference_df) gene_map <- .read_gene_map(geneinfo_file)

  # Signature genes. When split_direction = FALSE the whole signature is
  # truncated upfront; otherwise truncation is applied per direction below.
  parsed <- .read_signature(signature_file)
  sig <- parsed$sig
  gene_col <- parsed$gene_col
  log2fc_col <- parsed$log2fc_col
  if (!isTRUE(split_direction) && !is.null(max_genes)) sig <- utils::head(sig, max_genes)
  per_dir_max <- if (isTRUE(split_direction)) max_genes else NULL

  # Snapshot the full (pre-filter) signature so visualisation can show genes
  # that were dropped by the data-source intersection as gray placeholder columns.
  pre_filter_sig <- sig

  if (use_reference_df) {
    sig <- sig[sig[[gene_col]] %in% rownames(reference_mat), , drop = FALSE]
    if (nrow(sig) == 0) stop("No signature genes matched the row names of reference_df")
  } else {
    sig <- sig[sig[[gene_col]] %in% names(gene_map), , drop = FALSE]
    if (nrow(sig) == 0) stop("No signature genes mapped to geneinfo gene_id")
  }

  ordered_genes <- .direction_ordered_genes(sig, gene_col, log2fc_col, per_dir_max)
  if (length(ordered_genes) == 0) stop("No ordered genes available after down/up split")
  ordered_ids <- if (use_reference_df) ordered_genes else unname(gene_map[ordered_genes])
  logfc_map <- stats::setNames(sig[[log2fc_col]], sig[[gene_col]])

  # Full ordered gene list (pre data-source filter, same sort/truncation logic).
  all_ordered_genes <- .direction_ordered_genes(pre_filter_sig, gene_col, log2fc_col, per_dir_max)
  all_logfc_map <- stats::setNames(pre_filter_sig[[log2fc_col]], pre_filter_sig[[gene_col]])

  # Select perturbations by score
  tech <- .select_perturbations(results_df, pert_id_col, score_col, selected_drug,
                                selected_drug_col, max_perts, verbose, direction)
  sig_ids <- as.character(tech[[pert_id_col]])
  sig_ids <- unique(sig_ids[!.is_blank(sig_ids)])
  if (use_reference_df) sig_ids <- intersect(sig_ids, colnames(reference_mat))
  if (length(sig_ids) == 0) stop("No perturbation ids selected")
  if (verbose) message("Perturbations selected: ", length(sig_ids))

  # Scores in the same order as sig_ids (NA for any id not found in tech)
  score_lookup <- stats::setNames(tech[[score_col]], as.character(tech[[pert_id_col]]))
  if (use_reference_df) {
    z_mat <- reference_mat[ordered_ids, sig_ids, drop = FALSE]
  } else {
    z_mat <- fast_parse_gctx(fname = gctx_file, rid = as.character(ordered_ids), cid = sig_ids)
    present <- as.character(ordered_ids) %in% rownames(z_mat)
    ordered_genes <- ordered_genes[present]
    ordered_ids <- ordered_ids[present]
    sig_ids <- sig_ids[sig_ids %in% colnames(z_mat)]
    if (!length(ordered_ids) || !length(sig_ids)) {
      stop("No data available after GCTX extraction")
    }
    z_mat <- z_mat[as.character(ordered_ids), sig_ids, drop = FALSE]
  }
  sig_scores <- unname(score_lookup[sig_ids])
  sig_labels <- .siginfo_labels(sig_ids, siginfo_file)
  rownames(z_mat) <- ordered_genes
  z_plot <- t(z_mat) # rows = perturbations, cols = genes
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
    z_plot = z_plot,
    ordered_genes = ordered_genes,
    all_ordered_genes = all_ordered_genes,
    logfc_map = logfc_map,
    all_logfc_map = all_logfc_map,
    sig_ids = sig_ids,
    sig_labels = sig_labels,
    sig_scores = sig_scores
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
#' @param direction Order in which perturbations are picked by score:
#'   \code{"reversal"} (default) lowest first, \code{"mimic"} highest first.
#'   No perturbation is excluded by the sign of its score. Precomputed
#'   matrices are plotted as supplied; choose the direction when extracting
#'   them.
#' @param results_df Data frame (typically \code{tech_view_all}) containing at
#'   least perturbation id and score columns.
#' @param signature_file Path to a signature file with gene and log2FC columns,
#'   or a data frame with \code{Gene} and \code{log2FC} columns.
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
#'   If NULL, auto-detect from \code{perturbation_name}, \code{display_name},
#'   \code{pert_id}, \code{cmap_name}, \code{pert_name}.
#' @param pert_id_col Perturbation id column in \code{results_df}
#'   (default: \code{"sig_id"}; falls back to \code{"compound"} for older results).
#' @param score_col Score column in \code{results_df} (default: \code{"Score"}).
#' @param max_genes Maximum number of signature genes to use (default: 100).
#' @param max_perts Maximum number of perturbations to show (default: 60).
#' @param cluster_rows Logical; cluster perturbation rows (default: TRUE).
#' @param cluster_method Clustering linkage method for row clustering
#'   (default: \code{"complete"}).
#' @param show_row_dendrogram Logical; whether to draw a row dendrogram when
#'   \code{cluster_rows = TRUE} (default: TRUE).
#' @param cluster_cols Logical; cluster gene columns (default: FALSE - genes
#'   are kept in signature direction order: down then up).
#' @param cluster_method_cols Clustering linkage method for column clustering
#'   (default: \code{"complete"}).
#' @param show_col_dendrogram Logical; whether to draw a column dendrogram
#'   when \code{cluster_cols = TRUE} (default: TRUE).
#' @param split_direction Logical; if \code{TRUE}, split the heatmap into two
#'   side-by-side panels — up-regulated genes (log2FC > 0) on the left and
#'   down-regulated genes (log2FC < 0) on the right — with a small gap between
#'   them. Row order is determined by clustering the full matrix so both panels
#'   stay aligned. Default: \code{FALSE} (original single-panel behaviour).
#' @param gap_width Width of the gap between the two panels in mm when
#'   \code{split_direction = TRUE}. Default: \code{5}.
#' @param precomputed Optional; the output of \code{\link{extract_signature_zscores}}.
#'   If provided, skips all data-loading and GCTX-extraction steps - all
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
#'   ref_df <- read.csv(system.file("extdata", "example_reference_df.csv",
#'                                  package = "CONCERTDR"),
#'                      row.names = 1, check.names = FALSE)
#'   ref_df$gene_symbol <- rownames(ref_df)
#'   demo_sig <- colnames(ref_df)[1]
#'   plot_signature_direction_tile_barcode(
#'     results_df    = data.frame(sig_id = demo_sig, Score = -0.72),
#'     signature_file = sig_file,
#'     reference_df  = ref_df,
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
                                                  pert_id_col = "sig_id",
                                                  score_col = "Score",
                                                  max_genes = 100,
                                                  max_perts = 60,
                                                  cluster_rows = TRUE,
                                                  cluster_method = "complete",
                                                  show_row_dendrogram = TRUE,
                                                  cluster_cols = FALSE,
                                                  cluster_method_cols = "complete",
                                                  show_col_dendrogram = TRUE,
                                                  split_direction = FALSE,
                                                  gap_width = 5,
                                                  precomputed = NULL,
                                                  save_png = FALSE,
                                                  output_png = "barcode_heatmap.png",
                                                  output_zscores = NULL,
                                                  width = NULL,
                                                  height = NULL,
                                                  dpi = 150,
                                                  verbose = TRUE,
                                                  direction = c("reversal", "mimic")) {
  direction <- match.arg(direction)
  # ── data preparation ────────────────────────────────────────────────────────
  if (!is.null(precomputed)) {
    if (!is.list(precomputed) ||
        !all(c("z_plot", "ordered_genes", "logfc_map", "sig_ids") %in% names(precomputed))) {
      stop("precomputed must be the output of extract_signature_zscores()")
    }
    src <- precomputed
    if (verbose) message("Using precomputed matrix: ", nrow(src$z_plot),
                         " perturbations \u00d7 ", ncol(src$z_plot), " genes")
  } else {
    if (is.null(results_df)) stop("Provide results_df or precomputed")
    if (is.null(signature_file)) stop("Provide signature_file or precomputed")
    src <- extract_signature_zscores(
      results_df = results_df,
      signature_file = signature_file,
      reference_df = reference_df,
      gctx_file = gctx_file,
      geneinfo_file = geneinfo_file,
      siginfo_file = siginfo_file,
      data_dir = data_dir,
      selected_drug = selected_drug,
      selected_drug_col = selected_drug_col,
      pert_id_col = pert_id_col,
      score_col = score_col,
      max_genes = max_genes,
      max_perts = max_perts,
      split_direction = split_direction,
      direction = direction,
      output_zscores = output_zscores,
      verbose = verbose
    )
  }
  z_plot <- src$z_plot
  ordered_genes <- src$ordered_genes
  logfc_map <- src$logfc_map
  sig_ids <- src$sig_ids
  sig_scores <- src$sig_scores # NULL for old precomputed objects

  if (nrow(z_plot) == 0 || ncol(z_plot) == 0) {
    stop("No data available for heatmap")
  }

  # append score to each row label: "drug | dose | time | cell (score)"
  if (!is.null(sig_scores) && length(sig_scores) == nrow(z_plot)) {
    rownames(z_plot) <- paste0(
      rownames(z_plot), " (",
      formatC(sig_scores, format = "f", digits = 3),
      ")"
    )
  }

  # Reconstruct the full signature gene set: extract_signature_zscores stores
  # all_ordered_genes (before filtering to the data source). Genes dropped by
  # that filter are padded back as NA columns, shown as silver-gray placeholders.
  all_g <- src$all_ordered_genes
  if (!is.null(all_g)) {
    miss_g <- setdiff(all_g, ordered_genes)
    if (length(miss_g) > 0) {
      na_pad <- matrix(NA_real_, nrow = nrow(z_plot), ncol = length(miss_g),
                       dimnames = list(rownames(z_plot), miss_g))
      z_plot <- cbind(z_plot, na_pad)[, all_g, drop = FALSE]
      ordered_genes <- all_g
      if (!is.null(src$all_logfc_map)) logfc_map <- src$all_logfc_map
    }
  }

  # in_ref: TRUE = gene has actual z-score data; FALSE = NA placeholder
  in_ref <- vapply(seq_len(ncol(z_plot)),
                   function(j) any(!is.na(z_plot[, j])), logical(1))
  names(in_ref) <- ordered_genes
  any_out <- any(!in_ref)

  if (verbose && any_out) {
    message(sum(in_ref), " / ", length(ordered_genes),
            " signature genes found in data source; ",
            sum(!in_ref), " absent genes shown as gray placeholder columns")
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

  logfc_vals <- stats::setNames(as.numeric(logfc_map[ordered_genes]), ordered_genes)
  style <- .barcode_style(z_plot, logfc_vals)

  ttl <- if (is.null(selected_drug)) {
    paste0("Top ", nrow(z_plot), " perturbations by ", score_col)
  } else {
    paste0(selected_drug, " \u2013 Top ", nrow(z_plot), " perturbations by ", score_col)
  }

  up_genes <- ordered_genes[logfc_vals > 0]
  down_genes <- ordered_genes[logfc_vals < 0]
  use_split <- isTRUE(split_direction) && length(up_genes) > 0 && length(down_genes) > 0
  if (isTRUE(split_direction) && !use_split) {
    warning("split_direction = TRUE but all genes are in the same direction; drawing single heatmap")
  }

  # Options shared by every panel
  shared <- list(
    style = style,
    lfc_all = logfc_vals,
    row_method = cluster_method,
    col_method = cluster_method_cols,
    row_dend = isTRUE(show_row_dendrogram),
    col_dend = isTRUE(show_col_dendrogram)
  )
  clustered_rows <- isTRUE(cluster_rows) && nrow(z_plot) > 1L
  clustered_cols <- isTRUE(cluster_cols)
  dimmed <- function(genes) ifelse(in_ref[genes], "black", "#AAAAAA")

  if (!use_split) {
    if (clustered_cols && any_out) {
      # In-ref genes clustered, out-of-ref genes muted in signature order.
      genes_in <- ordered_genes[in_ref]
      genes_out <- ordered_genes[!in_ref]
      ht <- .barcode_panel(
        z_plot[, genes_in, drop = FALSE], genes_in, shared, "z-score", ttl,
        title_size = 11, cluster_rows = clustered_rows, show_row_names = TRUE,
        row_title = "Perturbations", cluster_cols = TRUE
      ) + .barcode_panel(
        z_plot[, genes_out, drop = FALSE], genes_out, shared, "z-score (not in ref)",
        muted = TRUE
      )
    } else {
      # Single panel in signature order; absent genes are gray with dimmed names.
      ht <- .barcode_panel(
        z_plot, ordered_genes, shared, "z-score", ttl, title_size = 11,
        na_col = style$muted_col, cluster_rows = clustered_rows,
        show_row_names = TRUE, row_title = "Perturbations",
        cluster_cols = clustered_cols, col_name_col = dimmed(ordered_genes)
      )
    }
    draw_args <- list()
  } else {
    # Row order comes from clustering the full matrix so all panels align.
    row_clust <- if (clustered_rows) {
      stats::hclust(stats::dist(z_plot), method = cluster_method)
    } else {
      FALSE
    }
    row_ord <- if (clustered_rows) row_clust$order else seq_len(nrow(z_plot))
    up_in <- in_ref[up_genes]
    down_in <- in_ref[down_genes]
    up_title <- "#8C510A"
    down_title <- "#01665E"
    # Panel by panel: a panel with no row clustering follows row_ord instead.
    rows_for <- function(first) {
      use_clust <- clustered_rows && first
      list(cluster_rows = if (use_clust) row_clust else FALSE,
           row_order = if (use_clust) NULL else row_ord)
    }

    if (clustered_cols && any_out) {
      # Up/down each split into an in-ref (clustered) and a muted out-of-ref panel.
      # The first in-ref panel carries row clustering, names and the legend.
      panels <- list()
      first <- TRUE
      if (any(up_in)) {
        g <- up_genes[up_in]
        panels$up_in <- do.call(.barcode_panel, c(
          list(z_plot[, g, drop = FALSE], g, shared, "z-score",
               paste0("Up-in-ref (", sum(up_in), ")"), title_col = up_title,
               show_row_names = FALSE, row_title = "Perturbations", cluster_cols = TRUE),
          rows_for(first)
        ))
        first <- FALSE
      }
      if (any(!up_in)) {
        g <- up_genes[!up_in]
        panels$up_out <- .barcode_panel(z_plot[, g, drop = FALSE], g, shared, "z-up-out",
                                        muted = TRUE, row_order = row_ord)
      }
      if (any(down_in)) {
        g <- down_genes[down_in]
        panels$down_in <- do.call(.barcode_panel, c(
          list(z_plot[, g, drop = FALSE], g, shared, if (first) "z-score" else "z-down-in",
               paste0("Down-in-ref (", sum(down_in), ")"), title_col = down_title,
               show_row_names = TRUE, cluster_cols = TRUE, show_legend = first,
               show_ann_name = FALSE, show_ann_legend = FALSE),
          rows_for(first)
        ))
        first <- FALSE
      }
      if (any(!down_in)) {
        g <- down_genes[!down_in]
        panels$down_out <- .barcode_panel(z_plot[, g, drop = FALSE], g, shared, "z-down-out",
                                          muted = TRUE, row_order = row_ord,
                                          show_row_names = !any(down_in))
      }
      ht <- Reduce(`+`, panels)
    } else {
      # Two panels (up | down); absent genes are gray with dimmed names.
      ht <- .barcode_panel(
        z_plot[, up_genes, drop = FALSE], up_genes, shared, "z-score",
        paste0("Up-regulated (", length(up_genes), " genes)"),
        title_col = up_title, na_col = style$muted_col,
        cluster_rows = row_clust, row_order = if (clustered_rows) NULL else row_ord,
        row_title = "Perturbations", cluster_cols = clustered_cols,
        col_name_col = dimmed(up_genes)
      ) + .barcode_panel(
        z_plot[, down_genes, drop = FALSE], down_genes, shared, "z-score-down",
        paste0("Down-regulated (", length(down_genes), " genes)"),
        title_col = down_title, na_col = style$muted_col, row_order = row_ord,
        show_row_names = TRUE, cluster_cols = clustered_cols,
        col_name_col = dimmed(down_genes), show_legend = FALSE,
        show_ann_name = FALSE, show_ann_legend = FALSE
      )
    }
    draw_args <- list(ht_gap = grid::unit(gap_width, "mm"), column_title = ttl,
                      column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"))
  }

  fig_size <- .barcode_fig_size(nrow(z_plot), ncol(z_plot), width, height)
  if (isTRUE(save_png)) {
    grDevices::png(filename = output_png, width = fig_size[1], height = fig_size[2],
                   units = "in", res = dpi)
    on.exit(grDevices::dev.off(), add = TRUE)
    .draw_barcode(ht, draw_args)
    if (verbose) message("Saved heatmap: ", output_png)
  } else {
    .draw_barcode(ht, draw_args)
  }

  invisible(list(
    z_plot = z_plot,
    ordered_genes = ordered_genes,
    sig_ids = sig_ids,
    output_png = if (isTRUE(save_png)) output_png else NULL,
    output_zscores = output_zscores
  ))
}


# ── Helpers for extract_signature_zscores() ───────────────────────────────────

#' Resolve a CMap file from explicit argument, option, env var or data_dir
#'
#' Candidates are tried in order: \code{explicit}, option
#' \code{CONCERTDR.<key>_file}, env var \code{CONCERTDR_<KEY>_FILE},
#' \code{data_dir/default}, then \code{default} in the working directory. The
#' last two are skipped when \code{default} is \code{NULL}.
#' @param explicit User-supplied path (may be \code{NULL}).
#' @param key Short file key (for example \code{"gctx"}), used to build the
#'   option and environment variable names.
#' @param default Default file name, or \code{NULL} for none.
#' @param data_dir Optional directory searched for \code{default}.
#' @return Path of the first candidate that exists, or \code{NULL}.
#' @keywords internal
.resolve_cmap_file <- function(explicit, key, default = NULL, data_dir = NULL) {
  in_dir <- if (!is.null(default) && !is.null(data_dir) && nzchar(data_dir)) {
    file.path(data_dir, default)
  }
  .first_existing_file(
    explicit,
    getOption(paste0("CONCERTDR.", key, "_file"), NULL),
    Sys.getenv(paste0("CONCERTDR_", toupper(key), "_FILE"), unset = ""),
    in_dir,
    default
  )
}

#' Coerce a reference data frame or matrix to a numeric gene x sample matrix
#' @param reference_df Data frame (optionally with a \code{gene_symbol} column)
#'   or matrix with genes as row names and perturbation ids as column names.
#' @return Numeric matrix with upper-cased gene row names.
#' @keywords internal
.reference_to_matrix <- function(reference_df) {
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
  gene_ids <- toupper(as.character(rownames(reference_df)))
  sample_ids <- colnames(reference_df)
  storage.mode(reference_mat) <- "double"
  rownames(reference_mat) <- gene_ids
  colnames(reference_mat) <- sample_ids
  reference_mat
}

#' Read a geneinfo file into a symbol -> gene_id lookup
#' @param geneinfo_file Path to a geneinfo file with \code{gene_symbol} and
#'   \code{gene_id} columns.
#' @return Named character vector of gene ids named by upper-cased symbol.
#' @keywords internal
.read_gene_map <- function(geneinfo_file) {
  geneinfo <- .read_cmap_table(geneinfo_file, "geneinfo_file")
  if (!all(c("gene_symbol", "gene_id") %in% names(geneinfo))) {
    stop("geneinfo_file must contain columns: gene_symbol, gene_id")
  }
  geneinfo$gene_symbol <- toupper(as.character(geneinfo$gene_symbol))
  geneinfo <- geneinfo[!is.na(geneinfo$gene_symbol) & !is.na(geneinfo$gene_id), , drop = FALSE]
  geneinfo <- .dedup_geneinfo(geneinfo)
  stats::setNames(as.character(geneinfo$gene_id), geneinfo$gene_symbol)
}

#' Read and clean a signature (gene and log2FC columns)
#' @param signature_file Signature file path or data frame. Falls back to the
#'   first two columns when \code{Gene}/\code{log2FC} are absent.
#' @return List with the cleaned data frame \code{sig} (upper-cased genes,
#'   numeric log2FC, incomplete rows dropped) and the column names
#'   \code{gene_col} and \code{log2fc_col}.
#' @keywords internal
.read_signature <- function(signature_file) {
  sig <- .read_cmap_table(signature_file, "signature_file")
  gene_col <- if ("Gene" %in% names(sig)) "Gene" else names(sig)[1]
  log2fc_col <- if ("log2FC" %in% names(sig)) "log2FC" else names(sig)[2]
  sig[[gene_col]] <- toupper(as.character(sig[[gene_col]]))
  sig[[log2fc_col]] <- .as_numeric(sig[[log2fc_col]])
  keep <- !.is_blank(sig[[gene_col]]) & !is.na(sig[[log2fc_col]])
  list(sig = sig[keep, , drop = FALSE], gene_col = gene_col, log2fc_col = log2fc_col)
}

#' Order signature genes down-regulated first, then up-regulated
#' @param sig Cleaned signature data frame (see \code{.read_signature}).
#' @param gene_col,log2fc_col Names of the gene and log2FC columns in \code{sig}.
#' @param max_genes Optional cap applied separately to each direction;
#'   \code{NULL} keeps all genes.
#' @return Character vector: most negative to zero, then most positive first.
#' @keywords internal
.direction_ordered_genes <- function(sig, gene_col, log2fc_col, max_genes = NULL) {
  fc <- sig[[log2fc_col]]
  down <- sig[fc < 0, , drop = FALSE]
  down <- down[order(down[[log2fc_col]], decreasing = FALSE), , drop = FALSE]
  up <- sig[fc > 0, , drop = FALSE]
  up <- up[order(up[[log2fc_col]], decreasing = TRUE), , drop = FALSE]
  if (!is.null(max_genes)) {
    down <- utils::head(down, max_genes)
    up <- utils::head(up, max_genes)
  }
  genes <- c(as.character(down[[gene_col]]), as.character(up[[gene_col]]))
  genes[nzchar(genes)]
}

#' Select the best-scoring perturbations from a results table
#' @param direction Order in which perturbations are picked by score:
#'   \code{"reversal"} (default) lowest first, \code{"mimic"} highest first.
#'   No perturbation is excluded by the sign of its score. Precomputed
#'   matrices are plotted as supplied; choose the direction when extracting
#'   them.
#' @param results_df Data frame with perturbation id and score columns.
#' @param pert_id_col,score_col Names of the id and score columns.
#' @param selected_drug Optional drug identifier used to filter rows.
#' @param selected_drug_col Column matched against \code{selected_drug}; when
#'   \code{NULL} it is auto-detected.
#' @param max_perts Maximum number of rows to keep.
#' @param verbose Logical; print progress messages.
#' @return Data frame sorted by the selected score direction, at most \code{max_perts} rows.
#' @keywords internal
.select_perturbations <- function(results_df, pert_id_col, score_col, selected_drug,
                                  selected_drug_col, max_perts, verbose, direction) {
  tech <- results_df
  tech[[score_col]] <- .as_numeric(tech[[score_col]])
  tech <- tech[!is.na(tech[[score_col]]) & !is.na(tech[[pert_id_col]]), , drop = FALSE]

  if (!is.null(selected_drug)) {
    drug_col <- selected_drug_col
    if (is.null(drug_col)) {
      candidates <- intersect(c("perturbation_name", "display_name", "pert_id", "cmap_name", "pert_name"), names(tech))
      if (length(candidates) == 0) stop("selected_drug given but no suitable selected_drug_col found")
      drug_col <- candidates[1]
    }
    keep <- toupper(as.character(tech[[drug_col]])) == toupper(selected_drug)
    tech <- tech[keep, , drop = FALSE]
    if (verbose) message("Filtered rows for selected_drug using ", drug_col, ": ", nrow(tech))
  }

  tech <- .rank_by_direction(tech, score_col, direction)
  utils::head(tech, max_perts)
}

#' Build readable perturbation labels from a siginfo file
#'
#' Labels look like \code{name | dose | time | cell}, where name is the first
#' non-empty of \code{cmap_name}, \code{pert_iname}, \code{pert_id}. Ids that are
#' not found (or when no siginfo file is available) keep their raw id.
#' @param sig_ids Character vector of \code{sig_id}s.
#' @param siginfo_file Path to a siginfo file, or \code{NULL}.
#' @return Character vector the same length as \code{sig_ids}.
#' @keywords internal
.siginfo_labels <- function(sig_ids, siginfo_file) {
  if (is.null(siginfo_file)) return(sig_ids)
  si <- .read_cmap_table(
    siginfo_file, "siginfo_file",
    select = c("sig_id", "cmap_name", "pert_iname", "pert_id", "pert_idose", "pert_itime", "cell_iname")
  )
  if (!"sig_id" %in% names(si)) return(sig_ids)
  si <- si[!duplicated(si$sig_id), , drop = FALSE]
  rownames(si) <- as.character(si$sig_id)

  label_one <- function(sid) {
    i <- match(sid, rownames(si))
    if (is.na(i)) return(sid)
    field <- function(nm) {
      if (!nm %in% names(si)) return(NULL)
      v <- as.character(si[[nm]][i])
      if (!.is_blank(v)) v
    }
    name_value <- unlist(lapply(c("cmap_name", "pert_iname", "pert_id"), field))[1]
    vals <- c(name_value, unlist(lapply(c("pert_idose", "pert_itime", "cell_iname"), field)))
    if (length(vals) == 0) sid else paste(vals, collapse = " | ")
  }
  vapply(sig_ids, label_one, character(1))
}

# ── Helpers for plot_signature_direction_tile_barcode() ───────────────────────

#' Build a grid gpar, dropping NULL entries
#' @param ... Graphical parameters; \code{NULL} values are ignored.
#' @return A \code{gpar} object.
#' @keywords internal
.barcode_gp <- function(...) {
  args <- list(...)
  do.call(grid::gpar, args[!vapply(args, is.null, logical(1))])
}

#' Colour scales shared by all barcode heatmap panels
#' @param z_plot Numeric matrix of z-scores (perturbations x genes).
#' @param logfc_vals Numeric vector of signature log2FC values.
#' @return List with \code{col_fun} (symmetric z-score scale), \code{lfc_col_fun}
#'   (BrBG-like log2FC scale, teal = down, brown = up), \code{muted_lfc_col_fun}
#'   (gray log2FC scale for absent genes) and \code{muted_col} (NA fill).
#' @keywords internal
.barcode_style <- function(z_plot, logfc_vals) {
  zlim <- max(abs(z_plot), na.rm = TRUE)
  if (!is.finite(zlim) || zlim == 0) zlim <- 10
  lim <- max(abs(logfc_vals), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  list(
    col_fun = circlize::colorRamp2(c(-zlim, 0, zlim), c("#3B4CC0", "#F7F7F7", "#B40426")),
    lfc_col_fun = circlize::colorRamp2(c(-lim, 0, lim), c("#01665E", "#F5F5F5", "#8C510A")),
    muted_lfc_col_fun = circlize::colorRamp2(c(-lim, 0, lim), c("#D0D0D0", "#E8E8E8", "#D0D0D0")),
    muted_col = "#CCCCCC"
  )
}

#' Build one panel of the barcode heatmap
#'
#' Every panel of every layout (single, clustered with absent genes, split
#' up/down) is produced here; the layouts differ only in the arguments.
#' @param z Numeric matrix (perturbations x genes) for this panel.
#' @param genes Gene names of the columns of \code{z}, used to look up log2FC.
#' @param shared List of options common to all panels: \code{style} (from
#'   \code{.barcode_style}), \code{lfc_all} (named log2FC vector),
#'   \code{row_method}, \code{col_method}, \code{row_dend}, \code{col_dend}.
#' @param name Heatmap name (must be unique within a heatmap list).
#' @param title Column title.
#' @param title_size,title_col Font size and colour of the column title.
#' @param muted Logical; draw a gray placeholder panel for absent genes (fixed
#'   title style, no legends, no clustering).
#' @param na_col Fill colour for \code{NA} cells.
#' @param cluster_rows \code{FALSE}, \code{TRUE} or an \code{hclust} object.
#' @param row_order Optional row order used when rows are not clustered.
#' @param show_row_names Logical; show perturbation labels.
#' @param row_title Row title.
#' @param cluster_cols Logical; cluster the gene columns.
#' @param col_name_col Optional column-name colours.
#' @param show_legend Logical; show the z-score legend.
#' @param show_ann_name,show_ann_legend Logical; show the annotation name and
#'   the annotation legend.
#' @return A \code{ComplexHeatmap::Heatmap}.
#' @keywords internal
.barcode_panel <- function(z, genes, shared, name, title = "(not in ref)",
                           title_size = 10, title_col = NULL, muted = FALSE,
                           na_col = "grey", cluster_rows = FALSE, row_order = NULL,
                           show_row_names = FALSE, row_title = character(0),
                           cluster_cols = FALSE, col_name_col = NULL,
                           show_legend = TRUE, show_ann_name = TRUE,
                           show_ann_legend = TRUE) {
  style <- shared$style
  if (muted) {
    na_col <- style$muted_col
    col_name_col <- "#999999"
    show_legend <- show_ann_name <- show_ann_legend <- FALSE
    cluster_cols <- FALSE
    title_gp <- .barcode_gp(fontsize = 9, col = "#999999", fontface = "italic")
  } else {
    title_gp <- .barcode_gp(fontsize = title_size, fontface = "bold", col = title_col)
  }

  annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Signature log2FC" = as.numeric(shared$lfc_all[genes]),
    col = list("Signature log2FC" = if (muted) style$muted_lfc_col_fun else style$lfc_col_fun),
    annotation_legend_param = list(
      "Signature log2FC" = list(
        title = "Signature log2FC",
        title_gp = grid::gpar(fontsize = 9),
        labels_gp = grid::gpar(fontsize = 8)
      )
    ),
    show_annotation_name = show_ann_name,
    annotation_name_gp = grid::gpar(fontsize = 9),
    show_legend = show_ann_legend
  )

  ComplexHeatmap::Heatmap(
    z,
    name = name,
    col = style$col_fun,
    na_col = na_col,
    cluster_rows = cluster_rows,
    clustering_method_rows = shared$row_method,
    row_order = row_order,
    show_row_dend = shared$row_dend && !isFALSE(cluster_rows),
    cluster_columns = cluster_cols,
    clustering_method_columns = shared$col_method,
    column_order = if (cluster_cols) NULL else seq_len(ncol(z)),
    show_column_dend = shared$col_dend && cluster_cols,
    top_annotation = annotation,
    show_row_names = show_row_names,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    row_names_max_width = grid::unit(7, "cm"),
    column_names_gp = .barcode_gp(fontsize = 8, col = col_name_col),
    column_names_rot = 60,
    column_title = title,
    column_title_gp = title_gp,
    row_title = row_title,
    row_title_gp = grid::gpar(fontsize = 10),
    heatmap_legend_param = list(
      title = "z-score",
      title_gp = grid::gpar(fontsize = 9),
      labels_gp = grid::gpar(fontsize = 8)
    ),
    show_heatmap_legend = show_legend,
    use_raster = TRUE,
    raster_quality = 2
  )
}

#' Figure size for the barcode heatmap
#' @param n_rows,n_cols Number of perturbations and genes.
#' @param width,height User-supplied size in inches, or \code{NULL} to
#'   compute from the data (about 0.22 in per gene and 0.28 in per perturbation).
#' @return Numeric vector \code{c(width, height)} in inches.
#' @keywords internal
.barcode_fig_size <- function(n_rows, n_cols, width = NULL, height = NULL) {
  c(if (is.null(width)) max(14, 4 + n_cols * 0.22 + 7) else width,
    if (is.null(height)) max(8, 2 + n_rows * 0.28) else height)
}

#' Draw a barcode heatmap (list) with the standard legend placement
#' @param ht A \code{Heatmap} or \code{HeatmapList}.
#' @param extra Named list of further arguments for \code{ComplexHeatmap::draw}.
#' @return Called for its side effect of drawing.
#' @keywords internal
.draw_barcode <- function(ht, extra = list()) {
  do.call(ComplexHeatmap::draw, c(
    list(ht,
         heatmap_legend_side = "right",
         annotation_legend_side = "right",
         padding = grid::unit(c(5, 20, 8, 5), "mm")),
    extra
  ))
}
