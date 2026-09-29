#' Extract Data from CMap GCTX Files for Specified Combinations
#'
#' @description Functions to extract expression data from CMap GCTX files based on
#' specified combinations of time points, dosages, and cell lines.
#' @return None; this is an internal documentation topic.
#'
#' @name cmap_extract
#' @keywords internal
NULL

#' Read first existing HDF5 dataset from candidate paths
#' @param fname HDF5/GCTX file path
#' @param paths Candidate dataset paths
#' @return Dataset content
#' @keywords internal
fast_read_meta <- function(fname, paths) {
  for (p in paths) {
    val <- tryCatch(rhdf5::h5read(fname, p), error = function(e) NULL)
    if (!is.null(val)) return(val)
  }
  stop("Could not read metadata from GCTX file. Checked paths: ",
       paste(paths, collapse = ", "))
}

#' Row and column ids of a GCTX file, read once per session
#'
#' Cached by file path and modification time, so repeated block reads of a
#' large file do not re-read and re-trim hundreds of thousands of ids.
#' @param fname GCTX file path.
#' @return List with \code{row} and \code{col} character vectors.
#' @keywords internal
.gctx_ids <- function(fname) {
  .camsum_cached(paste("ids", .camsum_file_key(fname)), {
    read_ids <- function(axis) {
      trimws(as.character(fast_read_meta(
        fname, paste0(c("/0/META/", "/META/"), axis, "/id"))))
    }
    list(row = read_ids("ROW"), col = read_ids("COL"))
  })
}

#' Read part of the GCTX data matrix
#'
#' Reads sorted, unique indices and explicitly restores the requested order
#' (including repeated indices), independently on both axes.
#' @param fname GCTX file path.
#' @param index List of row and column indices; \code{NULL} selects all.
#' @return Numeric matrix without dimnames.
#' @keywords internal
.gctx_read_matrix <- function(fname, index) {
  sorted <- lapply(index, function(i) if (is.null(i)) NULL else sort(unique(i)))
  mat <- tryCatch(
    rhdf5::h5read(fname, "/0/DATA/0/matrix", index = sorted, drop = FALSE),
    error = function(e) rhdf5::h5read(fname, "/DATA/0/matrix", index = sorted,
                                     drop = FALSE)
  )
  # Full-library CamSum blocks are already sorted; avoid copying those matrices.
  if (identical(index, sorted)) return(mat)
  rows <- if (is.null(index[[1L]])) seq_len(nrow(mat)) else match(index[[1L]], sorted[[1L]])
  cols <- if (is.null(index[[2L]])) seq_len(ncol(mat)) else match(index[[2L]], sorted[[2L]])
  mat[rows, cols, drop = FALSE]
}

#' Fast parser for GCTX matrices with optional row/column subsetting
#' @param fname Path to GCTX file
#' @param rid Optional character vector of row ids (gene ids)
#' @param cid Optional character vector of column ids (signature ids)
#' @return Numeric matrix with selected rows/columns
#' @keywords internal
fast_parse_gctx <- function(fname, rid = NULL, cid = NULL) {
  if (!file.exists(fname)) stop("GCTX file not found: ", fname)

  ids <- .gctx_ids(fname)
  select <- function(wanted, available) {
    idx <- if (is.null(wanted)) seq_along(available) else match(as.character(wanted), available)
    idx[!is.na(idx)]
  }
  row_idx <- select(rid, ids$row)
  col_idx <- select(cid, ids$col)

  mat <- if (length(row_idx) && length(col_idx)) {
    as.matrix(.gctx_read_matrix(fname, list(row_idx, col_idx)))
  } else {
    matrix(numeric(0), nrow = length(row_idx), ncol = length(col_idx))
  }
  dimnames(mat) <- list(ids$row[row_idx], ids$col[col_idx])
  mat
}

#' Select one feature per gene symbol
#'
#' Prefer landmark, then best inferred, then inferred, then unrecognised
#' feature spaces. Ties retain the first input row. Retained rows stay in
#' input order; a warning identifies every discarded gene ID and its replacement.
#' @param geneinfo_df Gene information with gene_id and gene_symbol columns,
#'   and optionally feature_space.
#' @return Gene information with unique gene symbols.
#' @keywords internal
.dedup_geneinfo <- function(geneinfo_df) {
  symbols <- as.character(geneinfo_df$gene_symbol)
  if (!anyDuplicated(symbols)) return(geneinfo_df)
  priority <- rep(4L, nrow(geneinfo_df))
  if ("feature_space" %in% names(geneinfo_df)) {
    priority <- match(as.character(geneinfo_df$feature_space),
                      c("landmark", "best inferred", "inferred"), nomatch = 4L)
  }
  preferred <- order(priority, seq_along(symbols))
  retained <- preferred[!duplicated(symbols[preferred])]
  dropped <- setdiff(seq_along(symbols), retained)
  replacements <- retained[match(symbols[dropped], symbols[retained])]
  warning("Duplicate gene symbols: ", paste(sprintf(
    "%s (dropped gene_id=%s; retained gene_id=%s)",
    symbols[dropped], geneinfo_df$gene_id[dropped],
    geneinfo_df$gene_id[replacements]), collapse = "; "),
    ". Priority: landmark > best inferred > inferred > other; ties keep the first row.",
    call. = FALSE)
  geneinfo_df[sort(retained), , drop = FALSE]
}

#' Get Unique Gene IDs and Gene Names
#'
#' @param geneinfo_df Data frame containing gene information
#' @param landmark Restrict to landmark genes before resolving duplicate symbols.
#' @return List with rid (gene IDs) and genenames (gene symbols)
#' @keywords internal
get_rid <- function(geneinfo_df, landmark = TRUE) {
  if (landmark) {
    geneinfo_df <- geneinfo_df[geneinfo_df$feature_space == "landmark", ]
  }
  geneinfo_df <- .dedup_geneinfo(geneinfo_df)
  list(
    rid = as.character(geneinfo_df$gene_id),
    genenames = as.character(geneinfo_df$gene_symbol)
  )
}

#' Process a Single Combination of Time, Dose, and Cell Line
#'
#' @param combination Data frame row containing time, dose, and cell line information
#' @param rid Vector of gene IDs to extract
#' @param genenames Vector of gene names corresponding to the gene IDs
#' @param sig_info Data frame containing signature information
#' @param gctx_file Path to the GCTX file
#' @param output_dir Directory to save output files
#' @return Path to the output file (invisibly)
#' @keywords internal
process_combination <- function(combination, rid, genenames, sig_info,
                                gctx_file = "level5_beta_trt_cp_n720216x12328.gctx",
                                output_dir = ".") {
  itime <- combination$itime
  idose <- combination$idose
  cell <- combination$cell

  message(sprintf("Processing: time=%s, dose=%s, cell=%s", itime, idose, cell))

  cid <- sig_info$sig_id[
    sig_info$pert_itime == itime &
      sig_info$pert_idose == idose &
      sig_info$cell_iname == cell
  ]

  # Underscores instead of spaces in the file name
  filename <- paste0("filtered_", paste(gsub(" ", "_", c(itime, idose, cell)),
                                        collapse = "_"), ".csv")
  output_file <- file.path(output_dir, filename)

  empty_reason <- NULL
  if (length(cid) == 0) {
    empty_reason <- "No matching signatures found for this combination"
  } else {
    mat <- fast_parse_gctx(fname = gctx_file, cid = cid, rid = rid)
    message(sprintf("Data dimensions: %d x %d", nrow(mat), ncol(mat)))
    if (nrow(mat) == 0 || ncol(mat) == 0) empty_reason <- "No data found for this combination"
  }

  if (is.null(empty_reason)) {
    df <- as.data.frame(mat)
    rownames(df) <- genenames[match(rownames(mat), rid)]
    utils::write.table(df, file = output_file, sep = "\t", quote = FALSE)
    message(sprintf("Successfully wrote data to %s", output_file))
  } else {
    message(empty_reason)
    writeLines(c(paste("#", empty_reason),
                 paste("# Time:", itime),
                 paste("# Dose:", idose),
                 paste("# Cell:", cell)), output_file)
  }

  invisible(output_file)
}

#' Process Multiple Combinations of Time, Dose, and Cell Line
#'
#' @param combinations_file Path to a file containing combinations to process
#' @param task_id Specific task ID to process (for SLURM array jobs)
#' @param geneinfo_file Path to the gene info file
#' @param siginfo_file Path to the signature info file
#' @param gctx_file Path to the GCTX file
#' @param output_dir Directory to save output files
#' @return List of processed files (invisibly)
#' @examples
#' combinations <- data.frame(
#'   itime = "24 h",
#'   idose = "1 uM",
#'   cell = "K562",
#'   stringsAsFactors = FALSE
#' )
#' combinations_file <- tempfile(fileext = ".tsv")
#' utils::write.table(
#'   combinations,
#'   combinations_file,
#'   sep = "\t",
#'   row.names = FALSE,
#'   quote = FALSE
#' )
#' file.exists(combinations_file)
#'
#' # With real CMap files available locally, run:
#' # process_combinations_file(
#' #   combinations_file = combinations_file,
#' #   geneinfo_file = "path/to/geneinfo_beta.txt",
#' #   siginfo_file = "path/to/siginfo_beta.txt",
#' #   gctx_file = "path/to/level5_beta_all_n1201944x12328.gctx",
#' #   output_dir = tempdir()
#' # )
#' @export
process_combinations_file <- function(combinations_file, task_id = NULL,
                                      geneinfo_file = "geneinfo_beta.txt",
                                      siginfo_file = "siginfo_beta.txt",
                                      gctx_file = "level5_beta_trt_cp_n720216x12328.gctx",
                                      output_dir = ".") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  combinations <- .read_cmap_table(combinations_file, "combinations_file", sep = "\t")
  message(sprintf("Loaded %d combinations from %s", nrow(combinations), combinations_file))

  message("Reading gene info file...")
  genes <- get_rid(.read_cmap_table(geneinfo_file, "geneinfo_file", sep = "\t"))

  message("Reading signature info file...")
  sig_info <- .read_cmap_table(siginfo_file, "siginfo_file", sep = "\t")
  sig_info <- sig_info[sig_info$pert_type == "trt_cp" & sig_info$is_hiq == 1, ]

  # An explicit task_id takes precedence over the SLURM array task ID; with
  # neither, every combination is processed sequentially.
  slurm_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  if (!is.null(task_id)) {
    task_id <- as.integer(task_id)
    message(sprintf("Using explicit task ID: %d", task_id))
  } else if (slurm_task_id != "") {
    task_id <- as.integer(slurm_task_id)
    message(sprintf("Using SLURM array task ID: %d", task_id))
  } else {
    message("No task ID specified. Processing all combinations sequentially.")
  }

  if (is.null(task_id)) {
    rows <- seq_len(nrow(combinations))
  } else if (task_id < 1 || task_id > nrow(combinations)) {
    stop("Invalid task ID: ", task_id, ". Must be between 1 and ", nrow(combinations))
  } else {
    rows <- task_id
  }

  output_files <- character()
  for (i in rows) {
    if (is.null(task_id)) {
      message(sprintf("\n--- Combination %d of %d ---\n", i, nrow(combinations)))
    }
    output_files <- c(output_files, process_combination(
      combinations[i, ], genes$rid, genes$genenames, sig_info, gctx_file, output_dir))
  }

  invisible(output_files)
}
