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
#' @param filter_quality Logical; whether to filter for pert_type="trt_cp" and is_hiq=1 (default: TRUE)
#' @param keep_all_genes Logical; whether to keep all genes (TRUE) or only common genes (FALSE) (default: TRUE)
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#' @param landmark Logical; whether to restrict to landmark genes only (default: TRUE)
#'
#' @return A data frame with expression data for all signatures, with annotation
#'         columns indicating sample metadata. Metadata is stored as an attribute.
#'
#' @examples
#' is.function(extract_cmap_data_from_siginfo)
#'
#' \donttest{
#' # Example 1: Use the original siginfo_beta.txt directly
#' reference_df <- extract_cmap_data_from_siginfo(
#'   siginfo_file = "path/to/siginfo_beta.txt",
#'   geneinfo_file = "path/to/geneinfo_beta.txt",
#'   gctx_file = "path/to/level5_beta_trt_cp_n720216x12328.gctx"
#' )
#'
#' # Example 2: Use a pre-filtered siginfo file
#' # First filter the siginfo
#' filtered_sig <- subset_siginfo_beta(
#'   "siginfo_beta.txt",
#'   interactive = FALSE,
#'   filters = list(
#'     pert_type = "trt_cp",
#'     pert_itime = c("6 h", "24 h"),
#'     pert_idose = "10 uM",
#'     cell_iname = c("A375", "MCF7")
#'   )
#' )
#' write.table(filtered_sig, "filtered_siginfo.txt", sep="\t", row.names=FALSE, quote=FALSE)
#'
#' # Then use the filtered siginfo
#' reference_df <- extract_cmap_data_from_siginfo(
#'   siginfo_file = "filtered_siginfo.txt",
#'   geneinfo_file = "path/to/geneinfo_beta.txt",
#'   gctx_file = "path/to/level5_beta_trt_cp_n720216x12328.gctx",
#'   filter_quality = FALSE  # Already filtered
#' )
#'
#' # Example 3: Process only first 5000 signatures for testing
#' reference_df <- extract_cmap_data_from_siginfo(
#'   siginfo_file = "path/to/siginfo_beta.txt",
#'   geneinfo_file = "path/to/geneinfo_beta.txt",
#'   gctx_file = "path/to/level5_beta_trt_cp_n720216x12328.gctx",
#'   max_signatures = 5000
#' )
#'
#' # Access the metadata
#' metadata <- attr(reference_df, "metadata")
#' table(metadata$cell)
#' table(metadata$time, metadata$dose)
#' }
#'
#' @export
extract_cmap_data_from_siginfo <- function(siginfo_file = "siginfo_beta.txt",
                                           geneinfo_file = "geneinfo_beta.txt",
                                           gctx_file = "level5_beta_trt_cp_n720216x12328.gctx",
                                           max_signatures = NULL,
                                           filter_quality = TRUE,
                                           keep_all_genes = TRUE,
                                           verbose = TRUE,
                                           landmark=TRUE) {

  # Handle geneinfo input
  if (is.character(geneinfo_file)) {
    if (!file.exists(geneinfo_file)) {
      stop("Gene info file not found: ", geneinfo_file)
    }
    if (verbose) message("Reading gene info file...")
    tryCatch({
      if (requireNamespace("data.table", quietly = TRUE)) {
        geneinfo_df <- data.table::fread(geneinfo_file, header = TRUE,
                                         stringsAsFactors = FALSE, data.table = FALSE)
      } else {
        geneinfo_df <- utils::read.table(geneinfo_file, sep = "\t", header = TRUE,
                                         stringsAsFactors = FALSE, quote = "",
                                         comment.char = "", fill = TRUE)
      }
    }, error = function(e) {
      stop("Error reading gene info file: ", e$message)
    })
  } else if (is.data.frame(geneinfo_file)) {
    geneinfo_df <- geneinfo_file
  } else {
    stop("geneinfo_file must be either a file path (character) or a data.frame.")
  }

  # Parse gene info
  result <- get_rid(geneinfo_df,landmark)
  rid <- result$rid
  genenames <- result$genenames
  if (verbose) message("Found ", length(rid), " landmark genes")

  # Handle siginfo input
  if (is.character(siginfo_file)) {
    if (!file.exists(siginfo_file)) {
      stop("Signature info file not found: ", siginfo_file)
    }
    if (verbose) message("Reading signature info file...")
    tryCatch({
      if (requireNamespace("data.table", quietly = TRUE)) {
        sig_info <- data.table::fread(siginfo_file, header = TRUE,
                                      stringsAsFactors = FALSE, data.table = FALSE)
      } else {
        sig_info <- utils::read.table(siginfo_file, sep = "\t", header = TRUE,
                                      stringsAsFactors = FALSE, quote = "",
                                      comment.char = "", fill = TRUE)
      }
    }, error = function(e) {
      stop("Error reading signature info file: ", e$message)
    })
  } else if (is.data.frame(siginfo_file)) {
    sig_info <- siginfo_file
  } else {
    stop("siginfo_file must be either a file path (character) or a data.frame.")
  }

  if (verbose) message("Loaded siginfo with ", nrow(sig_info), " signatures")

  # Apply quality filters if requested
  if (filter_quality) {
    original_count <- nrow(sig_info)

    # Filter for high quality
    if ("is_hiq" %in% names(sig_info)) {
      sig_info <- sig_info[sig_info$is_hiq == 1, ]
      if (verbose) message("Filtered to high-quality signatures: ", nrow(sig_info), " signatures")
    } else {
      warning("Column 'is_hiq' not found in siginfo file")
    }

    if (verbose && original_count > nrow(sig_info)) {
      message("Quality filtering reduced signatures from ", original_count, " to ", nrow(sig_info))
    }
  }

  # Check if sig_id column exists
  if (!"sig_id" %in% names(sig_info)) {
    stop("Required column 'sig_id' not found in siginfo file")
  }

  # Limit number of signatures if requested
  if (!is.null(max_signatures) && nrow(sig_info) > max_signatures) {
    sig_info <- sig_info[seq_len(max_signatures), ]
    if (verbose) message("Limited to first ", max_signatures, " signatures")
  }

  # Get all signature IDs
  all_cids <- sig_info$sig_id

  if (length(all_cids) == 0) {
    stop("No signatures found to process after filtering")
  }

  if (verbose) {
    message("\nProcessing ", length(all_cids), " signatures from siginfo file")

    # Show distribution of key parameters if available
    if ("pert_itime" %in% names(sig_info)) {
      time_table <- table(sig_info$pert_itime)
      message("Time points: ", paste(names(time_table), collapse = ", "))
    }
    if ("pert_idose" %in% names(sig_info)) {
      dose_table <- table(sig_info$pert_idose)
      message("Doses: ", paste(names(dose_table)[seq_len(min(5, length(dose_table)))], collapse = ", "),
              if(length(dose_table) > 5) "..." else "")
    }
    if ("cell_iname" %in% names(sig_info)) {
      cell_table <- table(sig_info$cell_iname)
      message("Cell lines: ", paste(names(cell_table)[seq_len(min(10, length(cell_table)))], collapse = ", "),
              if(length(cell_table) > 10) "..." else "")
    }
  }

  # Parse the gctx file using all signature IDs
  if (verbose) message("\nExtracting data from GCTX file...")

  tryCatch({
    pert_data <- cmapR::parse_gctx(
      fname = gctx_file,
      cid = all_cids,
      rid = rid
    )
  }, error = function(e) {
    stop("Error reading GCTX file: ", e$message)
  })

  # Get the data matrix
  mat <- cmapR::mat(pert_data)

  if (verbose) message(sprintf("Data dimensions: %d genes x %d signatures", nrow(mat), ncol(mat)))

  # Convert to data frame and set row names as gene names
  expression_data <- as.data.frame(mat)
  rownames(expression_data) <- genenames

  # Create metadata from sig_info
  # Match the order of signatures in the expression data
  metadata <- sig_info[match(colnames(expression_data), sig_info$sig_id), ]

  # Create a simplified metadata data frame with key columns
  metadata_df <- data.frame(
    sample_id = metadata$sig_id,
    stringsAsFactors = FALSE
  )

  # Add optional metadata columns if they exist
  optional_cols <- list(
    time = "pert_itime",
    dose = "pert_idose",
    cell = "cell_iname",
    pert_name = "pert_iname",
    pert_type = "pert_type",
    is_hiq = "is_hiq"
  )

  for (new_name in names(optional_cols)) {
    old_name <- optional_cols[[new_name]]
    if (old_name %in% names(metadata)) {
      metadata_df[[new_name]] <- metadata[[old_name]]
    }
  }

  # Add gene symbols as first column
  final_result <- data.frame(
    gene_symbol = rownames(expression_data),
    expression_data,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Add metadata as attribute
  attr(final_result, "metadata") <- metadata_df

  # Add full siginfo as attribute for reference
  attr(final_result, "full_siginfo") <- metadata

  if (verbose) {
    message(sprintf("\nSuccessfully extracted data:"))
    message(sprintf("  - %d genes", nrow(final_result)))
    message(sprintf("  - %d signatures", ncol(final_result) - 1))  # -1 for gene_symbol column
    message(sprintf("  - Metadata includes: %s", paste(names(metadata_df), collapse = ", ")))
  }

  return(final_result)
}
