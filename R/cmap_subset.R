#' Extract Combined Data from CMap GCTX Based on Config File
#'
#' @description Extract expression data from CMap GCTX files based on parameters
#' specified in a configuration file, returning a combined data frame.
#'
#' @param config_file Path to configuration file with selected parameters
#' @param geneinfo_file Path to the gene info file
#' @param siginfo_file Path to the signature info file
#' @param gctx_file Path to the GCTX file
#' @param keep_all_genes Logical; whether to keep all genes across combinations (TRUE) or only common genes (FALSE) (default: TRUE)
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return A data frame with combined results from all combinations, with annotation
#'         columns indicating time, dose, and cell line for each column
#'
#' @examples
#' \dontrun{
#' # First create and edit a configuration file
#' create_cmap_config_template()
#' # Edit the conf/cmap_options_template.conf file
#'
#' # Then extract data directly based on the config
#' results <- extract_cmap_data_from_config("conf/cmap_options_template.conf")
#'
#' # Write to CSV if needed
#' write.csv(results, file = "cmap_results.csv")
#' }
#'
#' @export
extract_cmap_data_from_config <- function(config_file,
                                          geneinfo_file = "geneinfo_beta.txt",
                                          siginfo_file = "siginfo_beta.txt",
                                          gctx_file = "level5_beta_trt_cp_n720216x12328.gctx",
                                          keep_all_genes = TRUE,
                                          verbose = TRUE,
                                          landmark=TRUE) {

  # Check if config file exists
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }

  # Generate combinations from config file
  if (verbose) message("Reading configuration from: ", config_file)
  combinations <- generate_combinations_from_config(config_file, verbose = verbose)

  if (nrow(combinations) == 0) {
    stop("No valid combinations found in config file.")
  }

  # Read gene info file or accept data.frame
  if (verbose) message("Reading gene info...")
  tryCatch({
    if (is.character(geneinfo_file)) {
      if (requireNamespace("data.table", quietly = TRUE)) {
        geneinfo_df <- data.table::fread(geneinfo_file, header = TRUE)
      } else {
        geneinfo_df <- utils::read.table(geneinfo_file, sep = "\t", header = TRUE,
                                         stringsAsFactors = FALSE, quote = "",
                                         comment.char = "", fill = TRUE)
      }
    } else if (is.data.frame(geneinfo_file)) {
      geneinfo_df <- geneinfo_file
    } else {
      stop("geneinfo_file must be a filename or a data.frame.")
    }

    result <- get_rid(geneinfo_df,landmark)
    rid <- result$rid
    genenames <- result$genenames

    if (verbose) message("Found ", length(rid), " landmark genes")
  }, error = function(e) {
    stop("Error reading gene info: ", e$message)
  })

  # Read siginfo file or accept data.frame
  if (verbose) message("Reading signature info...")
  tryCatch({
    if (is.character(siginfo_file)) {
      if (requireNamespace("data.table", quietly = TRUE)) {
        sig_info <- data.table::fread(siginfo_file, header = TRUE,
                                      stringsAsFactors = FALSE)
      } else {
        sig_info <- utils::read.table(siginfo_file, sep = "\t", header = TRUE,
                                      stringsAsFactors = FALSE, quote = "",
                                      comment.char = "", fill = TRUE)
      }
    } else if (is.data.frame(siginfo_file)) {
      sig_info <- siginfo_file
    } else {
      stop("siginfo_file must be a filename or a data.frame.")
    }

    if (verbose) message("Found ", nrow(sig_info), " treatment signatures")
  }, error = function(e) {
    stop("Error reading signature info: ", e$message)
  })


  # Initialize results list to store all dataframes
  result_list <- list()

  # Process each combination
  for (i in 1:nrow(combinations)) {
    combination <- combinations[i, ]
    itime <- combination$itime
    idose <- combination$idose
    cell <- combination$cell

    if (verbose) {
      message(sprintf("\n--- Combination %d of %d ---", i, nrow(combinations)))
      message(sprintf("Processing: time=%s, dose=%s, cell=%s", itime, idose, cell))
    }

    # Filter data
    filtered_df <- sig_info[
      sig_info$pert_itime == itime &
        sig_info$pert_idose == idose &
        sig_info$cell_iname == cell,
    ]

    cid <- filtered_df$sig_id

    if (length(cid) > 0) {
      # Parse the gctx file using CmapR
      pert_data <- cmapR::parse_gctx(
        fname = gctx_file,
        cid = cid,
        rid = rid
      )

      # Get the data matrix
      mat <- cmapR::mat(pert_data)

      if (verbose) message(sprintf("Data dimensions: %d x %d", nrow(mat), ncol(mat)))

      if (nrow(mat) > 0) {
        # Convert to data frame and set row names as gene names
        unique_idx <- !duplicated(genenames)
        df <- as.data.frame(mat[unique_idx, , drop = FALSE])
        rownames(df) <- genenames[unique_idx]

        # Add metadata columns for this combination
        metadata <- data.frame(
          sample_id = colnames(df),
          time = itime,
          dose = idose,
          cell = cell,
          stringsAsFactors = FALSE
        )

        # Add to result list
        result_list[[i]] <- list(data = df, metadata = metadata)
      } else {
        if (verbose) message("No data found for this combination")
      }
    } else {
      if (verbose) message("No matching signatures found for this combination")
    }
  }

  # Check if any results were found
  if (length(result_list) == 0) {
    warning("No data found for any of the specified combinations.")
    return(data.frame())
  }

  # Combine all results
  if (verbose) message("\nCombining results...")

  # Initialize variables to hold the combined data
  combined_data <- NULL
  combined_metadata <- NULL

  # Process each result based on the keep_all_genes parameter
  if (keep_all_genes) {
    # Strategy for keeping ALL genes across all combinations

    # First, collect all unique gene names
    all_genes <- character()
    for (result in result_list) {
      if (!is.null(result$data) && ncol(result$data) > 0) {
        all_genes <- union(all_genes, rownames(result$data))
      }
    }

    if (verbose) message("Total unique genes across all combinations: ", length(all_genes))

    # Create the combined matrix with all genes and all samples
    total_samples <- sum(sapply(result_list, function(r) {
      if (!is.null(r$data)) ncol(r$data) else 0
    }))

    # Initialize combined data with NAs
    combined_matrix <- matrix(NA, nrow = length(all_genes), ncol = total_samples)
    rownames(combined_matrix) <- all_genes

    # Collect all sample IDs
    all_sample_ids <- character()

    # Fill in the data for each combination
    col_start <- 1
    for (result in result_list) {
      if (!is.null(result$data) && ncol(result$data) > 0) {
        data_df <- result$data
        num_cols <- ncol(data_df)

        # Add sample IDs
        all_sample_ids <- c(all_sample_ids, colnames(data_df))

        # Fill in data for current combination
        for (gene in rownames(data_df)) {
          idx <- match(gene, all_genes)
          if (!is.na(idx)) {
            combined_matrix[idx, col_start:(col_start + num_cols - 1)] <-
              as.numeric(data_df[gene, ])
          }
        }

        # Combine metadata
        if (is.null(combined_metadata)) {
          combined_metadata <- result$metadata
        } else {
          combined_metadata <- rbind(combined_metadata, result$metadata)
        }

        # Update column start for next combination
        col_start <- col_start + num_cols
      }
    }

    # Set column names
    colnames(combined_matrix) <- all_sample_ids

    # Convert to data frame
    combined_data <- as.data.frame(combined_matrix)

  } else {
    # Original behavior - only keep genes common to all combinations

    for (result in result_list) {
      if (!is.null(result$data) && ncol(result$data) > 0) {
        if (is.null(combined_data)) {
          combined_data <- result$data
          combined_metadata <- result$metadata
        } else {
          # Find common genes
          common_genes <- intersect(rownames(combined_data), rownames(result$data))
          if (verbose) message("Common genes with previous combinations: ", length(common_genes))

          # Append columns to the existing dataframe, keeping only common genes
          combined_data <- cbind(combined_data[common_genes, ], result$data[common_genes, ])
          combined_metadata <- rbind(combined_metadata, result$metadata)
        }
      }
    }
  }

  # Final check
  if (is.null(combined_data)) {
    warning("No valid data found to combine.")
    return(data.frame())
  }

  # Create result metadata
  if (verbose) {
    message(sprintf("Final combined data dimensions: %d genes × %d samples",
                    nrow(combined_data), ncol(combined_data)))
  }

  # Add gene symbols as a column
  final_result <- data.frame(gene_symbol = rownames(combined_data), combined_data,
                             check.names = FALSE, stringsAsFactors = FALSE)

  # Add metadata as attributes
  attr(final_result, "metadata") <- combined_metadata

  return(final_result)
}
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
#'
#' @return A data frame with expression data for all signatures, with annotation
#'         columns indicating sample metadata. Metadata is stored as an attribute.
#'
#' @examples
#' \dontrun{
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
    sig_info <- sig_info[1:max_signatures, ]
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
      message("Doses: ", paste(names(dose_table)[1:min(5, length(dose_table))], collapse = ", "),
              if(length(dose_table) > 5) "..." else "")
    }
    if ("cell_iname" %in% names(sig_info)) {
      cell_table <- table(sig_info$cell_iname)
      message("Cell lines: ", paste(names(cell_table)[1:min(10, length(cell_table))], collapse = ", "),
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

  if (verbose) message(sprintf("Data dimensions: %d genes × %d signatures", nrow(mat), ncol(mat)))

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
