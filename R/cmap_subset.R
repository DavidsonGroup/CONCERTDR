#' Extract CMap Data Using Siginfo File Directly
#'
#' @description Extract expression data from CMap GCTX files using all signatures
#' present in the provided siginfo file. This function processes all signatures
#' found in the siginfo file without requiring a configuration file.
#'
#' @param siginfo_file Path to the signature info file (can be original or pre-filtered)
#' @param geneinfo_file Path to the gene info file
#' @param gctx_file Path to the GCTX file
#' @param max_signatures Integer; maximum number of signatures to process (default: NULL for all)
#' @param filter_quality Logical; whether to keep only high-quality signatures
#'   (\code{is_hiq == 1}) (default: TRUE)
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#' @param landmark Logical; whether to restrict to landmark genes only (default: TRUE)
#' @details Duplicate gene symbols are resolved by preferring landmark,
#'   then best inferred, then inferred features. Other feature spaces rank
#'   last and ties retain the first input row. A warning lists discarded and
#'   retained gene IDs. Gene symbols are not renamed or averaged. GCTX data
#'   are read in sorted index order and restored to the requested gene and
#'   signature order explicitly.
#'
#' @return A data frame with expression data for all signatures, with annotation
#'         columns indicating sample metadata. Metadata is stored as an attribute.
#'
#' @examples
#' # Build a tiny GCTX file so the example is self-contained and runnable.
#' gctx_file <- tempfile(fileext = ".gctx")
#' rhdf5::h5createFile(gctx_file)
#' for (group in c("0", "0/DATA", "0/DATA/0", "0/META",
#'                 "0/META/ROW", "0/META/COL")) {
#'   rhdf5::h5createGroup(gctx_file, group)
#' }
#' rhdf5::h5write(matrix(c(1.2, -0.4, 0.7, -1.1), nrow = 2),
#'                gctx_file, "0/DATA/0/matrix")
#' rhdf5::h5write(c("1", "2"), gctx_file, "0/META/ROW/id")
#' rhdf5::h5write(c("SIG_A", "SIG_B"), gctx_file, "0/META/COL/id")
#'
#' siginfo <- data.frame(
#'   sig_id = c("SIG_A", "SIG_B"), is_hiq = 1,
#'   pert_itime = "24 h", pert_idose = "1 uM", cell_iname = "A375"
#' )
#' geneinfo <- data.frame(
#'   gene_id = c("1", "2"), gene_symbol = c("GENE1", "GENE2"),
#'   feature_space = "landmark"
#' )
#' reference_df <- extract_cmap_data_from_siginfo(
#'   siginfo_file = siginfo,
#'   geneinfo_file = geneinfo,
#'   gctx_file = gctx_file,
#'   verbose = FALSE
#' )
#' reference_df
#' attr(reference_df, "metadata")
#' unlink(gctx_file)
#'
#' @export
extract_cmap_data_from_siginfo <- function(siginfo_file = "siginfo_beta.txt",
                                           geneinfo_file = "geneinfo_beta.txt",
                                           gctx_file = "level5_beta_trt_cp_n720216x12328.gctx",
                                           max_signatures = NULL,
                                           filter_quality = TRUE,
                                           verbose = TRUE,
                                           landmark = TRUE) {
  say <- function(...) if (verbose) message(...)

  say("Reading gene info file...")
  genes <- get_rid(.read_cmap_table(geneinfo_file, "geneinfo_file"), landmark)
  say("Found ", length(genes$rid), if (landmark) " landmark", " genes")

  say("Reading signature info file...")
  sig_info <- .read_cmap_table(siginfo_file, "siginfo_file")
  say("Loaded siginfo with ", nrow(sig_info), " signatures")
  sig_info <- .filter_siginfo(sig_info, filter_quality, max_signatures, verbose)

  say("\nProcessing ", nrow(sig_info), " signatures from siginfo file")
  if (verbose) {
    .message_column_summary(sig_info, "pert_itime", "Time points")
    .message_column_summary(sig_info, "pert_idose", "Doses", max_shown = 5L)
    .message_column_summary(sig_info, "cell_iname", "Cell lines", max_shown = 10L)
  }
  say("\n", strrep("=", 60), "\n",
      "NOTE: Extracting data from GCTX file...\n",
      "This step may take a long time depending on the number of\n",
      "signatures (", nrow(sig_info), ") and the size of the GCTX file.\n",
      "Please be patient and do not interrupt the process.\n",
      strrep("=", 60), "\n")

  mat <- tryCatch(
    fast_parse_gctx(fname = gctx_file, cid = sig_info$sig_id, rid = genes$rid),
    error = function(e) stop("Error reading GCTX file: ", e$message)
  )
  say(sprintf("Data dimensions: %d genes x %d signatures", nrow(mat), ncol(mat)))

  if (!nrow(mat) || !ncol(mat)) stop("No matching genes or signatures in GCTX file")
  expression_data <- as.data.frame(mat)
  rownames(expression_data) <- genes$genenames[match(rownames(mat), genes$rid)]

  # Metadata rows follow the column order of the expression data
  full_siginfo <- sig_info[match(colnames(expression_data), sig_info$sig_id), ]
  metadata <- .simplify_siginfo(full_siginfo)

  result <- data.frame(gene_symbol = rownames(expression_data), expression_data,
                       check.names = FALSE)
  attr(result, "metadata") <- metadata
  attr(result, "full_siginfo") <- full_siginfo

  say("\nSuccessfully extracted data:\n",
      sprintf("  - %d genes\n", nrow(result)),
      sprintf("  - %d signatures\n", ncol(result) - 1L),
      "  - Metadata includes: ", paste(names(metadata), collapse = ", "))
  result
}

#' Apply the quality filter and signature limit to a siginfo table
#' @param sig_info Siginfo data frame.
#' @param filter_quality Whether to keep only \code{is_hiq == 1} rows.
#' @param max_signatures Maximum number of signatures, or \code{NULL}.
#' @param verbose Whether to print progress messages.
#' @return The filtered data frame. Errors when \code{sig_id} is missing or
#'   no signature is left.
#' @keywords internal
.filter_siginfo <- function(sig_info, filter_quality, max_signatures, verbose) {
  say <- function(...) if (verbose) message(...)

  if (filter_quality) {
    original_count <- nrow(sig_info)
    if ("is_hiq" %in% names(sig_info)) {
      sig_info <- sig_info[sig_info$is_hiq == 1, ]
      say("Filtered to high-quality signatures: ", nrow(sig_info), " signatures")
    } else {
      warning("Column 'is_hiq' not found in siginfo file")
    }
    if (original_count > nrow(sig_info)) {
      say("Quality filtering reduced signatures from ", original_count, " to ",
          nrow(sig_info))
    }
  }

  if (!"sig_id" %in% names(sig_info)) {
    stop("Required column 'sig_id' not found in siginfo file")
  }
  if (!is.null(max_signatures) && nrow(sig_info) > max_signatures) {
    sig_info <- sig_info[seq_len(max_signatures), ]
    say("Limited to first ", max_signatures, " signatures")
  }
  if (nrow(sig_info) == 0) stop("No signatures found to process after filtering")
  sig_info
}

#' Report the distinct values of one siginfo column
#' @param sig_info Siginfo data frame.
#' @param col Column name; nothing is reported when it is absent.
#' @param label Label used in the message.
#' @param max_shown Maximum number of values listed.
#' @return NULL, invisibly.
#' @keywords internal
.message_column_summary <- function(sig_info, col, label, max_shown = Inf) {
  if (!col %in% names(sig_info)) return(invisible())
  values <- names(table(sig_info[[col]]))
  message(label, ": ", paste(utils::head(values, max_shown), collapse = ", "),
          if (length(values) > max_shown) "...")
}

#' Reduce a siginfo table to the metadata columns used downstream
#' @param sig_info Siginfo data frame.
#' @return Data frame with \code{sample_id} and, when present in
#'   \code{sig_info}, \code{time}, \code{dose}, \code{cell},
#'   \code{pert_name}, \code{pert_type} and \code{is_hiq}.
#' @keywords internal
.simplify_siginfo <- function(sig_info) {
  renamed <- c(time = "pert_itime", dose = "pert_idose", cell = "cell_iname",
               pert_name = "pert_iname", pert_type = "pert_type",
               is_hiq = "is_hiq")
  renamed <- renamed[renamed %in% names(sig_info)]
  metadata <- data.frame(sample_id = sig_info$sig_id)
  for (new_name in names(renamed)) {
    metadata[[new_name]] <- sig_info[[renamed[[new_name]]]]
  }
  metadata
}
