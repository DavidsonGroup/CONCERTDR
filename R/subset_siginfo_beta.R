#' Subset siginfo_beta file interactively with smart preview
#'
#' Enhanced interactive filtering that shows what options will be available
#' in subsequent parameters based on your current selection, helping you make
#' informed decisions about the filtering sequence.
#'
#' @param siginfo_file Path to siginfo_beta.txt file
#' @param output_file Path to save the filtered siginfo file (optional)
#' @param interactive Logical; whether to run in interactive mode (default: TRUE)
#' @param filters List of pre-defined filters to apply in non-interactive mode
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#' @param show_preview Logical; whether to show preview of subsequent options (default: TRUE)
#'
#' @return A filtered data frame of the siginfo file
#'
#' @examples
#' ex_sig <- data.frame(
#'   pert_type = c("trt_cp", "trt_cp", "trt_sh"),
#'   pert_itime = c("6 h", "24 h", "24 h"),
#'   pert_idose = c("10 uM", "10 uM", "5 uM"),
#'   cell_iname = c("A375", "MCF7", "A375"),
#'   stringsAsFactors = FALSE
#' )
#' ex_sig_file <- tempfile(fileext = ".txt")
#' write.table(ex_sig, ex_sig_file, sep = "\t", row.names = FALSE, quote = FALSE)
#' subset_siginfo_beta(
#'   ex_sig_file,
#'   interactive = FALSE,
#'   filters = list(pert_type = "trt_cp"),
#'   verbose = FALSE,
#'   show_preview = FALSE
#' )
#'
#' \donttest{
#' ex_file <- system.file("extdata", "example_siginfo.txt", package = "CONCERTDR")
#'
#' # Interactive mode (run only in an interactive R session)
#' if (interactive() && nzchar(ex_file)) {
#'   filtered_siginfo <- subset_siginfo_beta(ex_file)
#' }
#'
#' # Non-interactive mode with pre-defined filters
#' if (nzchar(ex_file)) {
#'   filtered_siginfo <- subset_siginfo_beta(
#'     ex_file,
#'     interactive = FALSE,
#'     filters = list(
#'       pert_type = "trt_cp",
#'       pert_itime = c("6 h", "24 h"),
#'       pert_idose = "10 uM",
#'       cell_iname = c("A375", "MCF7")
#'     ),
#'     verbose = FALSE,
#'     show_preview = FALSE
#'   )
#' }
#' }
#'
#' @export
subset_siginfo_beta <- function(siginfo_file,
                                output_file = NULL,
                                interactive = TRUE,
                                filters = NULL,
                                verbose = TRUE,
                                show_preview = TRUE) {
  
  # Check if file exists
  if (!file.exists(siginfo_file)) {
    stop("Siginfo file not found: ", siginfo_file)
  }
  
  # Read the siginfo file
  if (verbose) message("Reading siginfo file: ", siginfo_file)
  
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
    stop("Error reading siginfo file: ", e$message)
  })
  
  if (requireNamespace("data.table", quietly = TRUE)) {
    sig_info <- data.table::fread(siginfo_file, header = TRUE,
                                  stringsAsFactors = FALSE, data.table = FALSE)
  } else {
    sig_info <- utils::read.table(siginfo_file, sep = "\t", header = TRUE,
                                  stringsAsFactors = FALSE, quote = "",
                                  comment.char = "", fill = TRUE)
  }

  if (verbose) {
    message(sprintf("Loaded siginfo data with %d rows and %d columns",
                    nrow(sig_info), ncol(sig_info)))
  }
  
  # Store original count
  original_count <- nrow(sig_info)
  
  # Define the key columns to filter on
  key_columns <- c("pert_type", "pert_itime", "pert_idose", "cell_iname")
  
  # Check if all key columns exist
  missing_cols <- setdiff(key_columns, names(sig_info))
  if (length(missing_cols) > 0) {
    warning("Missing columns: ", paste(missing_cols, collapse = ", "))
    key_columns <- intersect(key_columns, names(sig_info))
  }
  
  # Apply filters
  filtered_data <- sig_info
  
  if (interactive) {
    # Interactive mode with smart preview
    message("\n", paste(rep("=", 80), collapse = ""))
    message("INTERACTIVE SIGINFO FILTERING WITH SMART PREVIEW")
    message(paste(rep("=", 80), collapse = ""))
    message(sprintf("Starting with %d signatures", original_count))

    for (col_idx in seq_along(key_columns)) {
      col <- key_columns[col_idx]
      
      # Get unique values from current filtered data
      unique_vals <- sort(unique(filtered_data[[col]]))
      
      if (length(unique_vals) == 0) {
        next
      }
      
      # Display current status
      current_count <- nrow(filtered_data)
      current_proportion <- round((current_count / original_count) * 100, 1)

      message("\n", paste(rep("=", 80), collapse = ""))
      message(sprintf("STEP %d/%d: Filtering by %s", col_idx, length(key_columns), col))
      message(sprintf("Current signatures: %d (%.1f%% of original)",
              current_count, current_proportion))

      # Show what each selection would yield and preview downstream options
      if (show_preview && col_idx < length(key_columns)) {
        message("\nPREVIEW: Impact on downstream parameters")
        message(paste(rep("-", 80), collapse = ""))
      }

      message("\nAvailable values for ", col, ":")

      # Calculate impact for each value
      value_impacts <- list()
      for (i in seq_along(unique_vals)) {
        val <- unique_vals[i]
        temp_data <- filtered_data[filtered_data[[col]] == val, ]
        count <- nrow(temp_data)
        proportion_of_current <- round((count / current_count) * 100, 1)
        proportion_of_original <- round((count / original_count) * 100, 1)
        
        value_impacts[[val]] <- list(
          data = temp_data,
          count = count,
          prop_current = proportion_of_current,
          prop_original = proportion_of_original
        )

        preview_line <- sprintf("  %d: %s [%d sigs, %.1f%% current, %.1f%% original]",
              i, val, count, proportion_of_current, proportion_of_original)

        # Show preview of option counts in next parameters
        if (show_preview && col_idx < length(key_columns)) {
          remaining_cols <- key_columns[(col_idx + 1):length(key_columns)]
          preview_info <- ""

          for (next_col in remaining_cols[seq_len(min(2, length(remaining_cols)))]) {  # Show max 2 next columns
            if (next_col %in% names(temp_data)) {
              next_count <- length(unique(temp_data[[next_col]]))
              preview_info <- paste0(preview_info, sprintf(" -> %s: %d options",
                                                           next_col, next_count))
            }
          }
          if (preview_info != "") preview_line <- paste0(preview_line, preview_info)
        }
        message(preview_line)
      }
      
      # Special preview option: show what combinations would be available
      if (show_preview && col_idx < length(key_columns)) {
        message("\n[Tip] Type 'preview X' to see detailed downstream options for choice X")
      }
      
      # Get user selection
      message("\nSelect values to keep:")
      message("  - Enter comma-separated numbers (e.g., '1,3,5')")
      message("  - Enter 'all' to keep all values")
      message("  - Enter 'skip' to skip this parameter")
      if (show_preview) {
        message("  - Enter 'preview X' to see detailed downstream options for choice X")
      }
      
      while (TRUE) {
        selection <- readline(prompt = "> ")
        
        if (tolower(selection) == "skip" || selection == "") {
          message("Skipping ", col)
          break
          
        } else if (tolower(selection) == "all") {
          message("Keeping all values for ", col)
          break
          
        } else if (grepl("^preview\\s+", tolower(selection))) {
          # Detailed preview mode
          preview_input <- gsub("^preview\\s+", "", selection, ignore.case = TRUE)
          preview_idx <- as.integer(preview_input)
          
          if (!is.na(preview_idx) && preview_idx >= 1 && preview_idx <= length(unique_vals)) {
            preview_val <- unique_vals[preview_idx]
            preview_data <- value_impacts[[preview_val]]$data

            message("\n", paste(rep("~", 60), collapse = ""))
            message(sprintf("DETAILED PREVIEW: If you select '%s'", preview_val))
            message(sprintf("This would retain %d signatures (%.1f%% of original)",
                            value_impacts[[preview_val]]$count,
                            value_impacts[[preview_val]]$prop_original))

            # Show what would be available in ALL remaining parameters
            remaining_cols <- key_columns[(col_idx + 1):length(key_columns)]
            for (next_col in remaining_cols) {
              if (next_col %in% names(preview_data)) {
                next_unique <- sort(unique(preview_data[[next_col]]))
                next_counts <- table(preview_data[[next_col]])

                message(sprintf("\nAvailable %s options (%d):", next_col, length(next_unique)))
                for (j in seq_along(next_unique)) {
                  next_val <- next_unique[j]
                  next_count <- next_counts[next_val]
                  next_prop <- round((next_count / nrow(preview_data)) * 100, 1)
                  message(sprintf("    %s: %d sigs (%.1f%%)", next_val, next_count, next_prop))

                  # Limit display to avoid overwhelming output
                  if (j >= 10) {
                    message(sprintf("    ... and %d more options", length(next_unique) - j))
                    break
                  }
                }
              }
            }
            message(paste(rep("~", 60), collapse = ""), "\n")

          } else {
            message("Invalid preview selection. Please enter a valid number.")
          }
          
        } else {
          # Regular selection
          indices <- as.integer(unlist(strsplit(selection, ",")))
          indices <- indices[!is.na(indices)]
          
          # Validate indices
          valid_indices <- indices[indices >= 1 & indices <= length(unique_vals)]
          
          if (length(valid_indices) > 0) {
            selected_vals <- unique_vals[valid_indices]
            
            # Calculate combined impact
            combined_data <- filtered_data[filtered_data[[col]] %in% selected_vals, ]
            new_count <- nrow(combined_data)
            new_proportion <- round((new_count / original_count) * 100, 1)

            message(sprintf("\nThis selection will retain %d signatures (%.1f%% of original)",
                            new_count, new_proportion))
            message("Selected values: ", paste(selected_vals, collapse = ", "))

            # Show combined preview of next option counts
            if (show_preview && col_idx < length(key_columns)) {
              remaining_cols <- key_columns[(col_idx + 1):length(key_columns)]
              next_col <- remaining_cols[1]
              if (next_col %in% names(combined_data)) {
                next_count <- length(unique(combined_data[[next_col]]))
                message(sprintf("Next parameter (%s) will have %d options",
                                next_col, next_count))
              }
            }
            
            # Confirm selection
            confirm <- readline(prompt = "Confirm this selection? (y/n): ")
            
            if (tolower(confirm) == "y" || tolower(confirm) == "yes") {
              filtered_data <- combined_data
              
              if (verbose) {
                message(sprintf("[OK] Applied filter on %s: %d values selected, %d rows remaining (%.1f%%)",
                                col, length(selected_vals), nrow(filtered_data), new_proportion))
              }
              break
            } else {
              message("Selection cancelled. You can make a new selection.")
              # Continue the loop to allow retry
            }
            
          } else {
            message("Invalid selection. Please enter valid numbers, 'all', 'skip', or 'preview X'.")
          }
        }
      }
    }
    
  } else {
    # Non-interactive mode - use provided filters
    if (!is.null(filters)) {
      for (col in names(filters)) {
        if (col %in% names(filtered_data)) {
          old_count <- nrow(filtered_data)
          filtered_data <- filtered_data[filtered_data[[col]] %in% filters[[col]], ]
          new_count <- nrow(filtered_data)
          new_proportion <- round((new_count / original_count) * 100, 1)
          
          if (verbose) {
            message(sprintf("Filtered %s: %d -> %d signatures (%.1f%% of original)",
                            col, old_count, new_count, new_proportion))
          }
        }
      }
    }
  }
  
  # Show final summary
  final_count <- nrow(filtered_data)
  final_proportion <- round((final_count / original_count) * 100, 1)
  
  if (verbose) {
    message("\n", paste(rep("=", 80), collapse = ""))
    message("FILTERING COMPLETE!")
    message(paste(rep("=", 80), collapse = ""))
    message(sprintf("Original signatures: %s", format(original_count, big.mark = ",")))
    message(sprintf("Final signatures: %s (%.1f%% retained)",
                    format(final_count, big.mark = ","), final_proportion))
    message(sprintf("Signatures removed: %s (%.1f%% reduction)",
                    format(original_count - final_count, big.mark = ","), 100 - final_proportion))

    # Show final distribution
    message("\n[Summary] Final distribution by parameter:")
    for (col in key_columns) {
      if (col %in% names(filtered_data) && nrow(filtered_data) > 0) {
        val_counts <- sort(table(filtered_data[[col]]), decreasing = TRUE)
        message(sprintf("\n%s (%d unique values):", col, length(val_counts)))

        # Show top values
        show_n <- min(length(val_counts), 8)
        for (i in seq_len(show_n)) {
          val_prop <- round((val_counts[i] / final_count) * 100, 1)
          message(sprintf("    %s: %s (%.1f%%)",
                          names(val_counts)[i],
                          format(val_counts[i], big.mark = ","),
                          val_prop))
        }
        if (length(val_counts) > 8) {
          remaining_total <- sum(val_counts[9:length(val_counts)])
          message(sprintf("    ... and %d more values: %s signatures",
                          length(val_counts) - 8,
                          format(remaining_total, big.mark = ",")))
        }
      }
    }
  }
  
  # Save filtered data if output file is specified
  if (!is.null(output_file)) {
    if (verbose) message("\n[Saving] Saving filtered data to: ", output_file)

    # Create directory if it doesn't exist
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir) && output_dir != ".") {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Write the file
    utils::write.table(filtered_data, file = output_file,
                       sep = "\t", row.names = FALSE, quote = FALSE)
    
    if (verbose) {
      message(sprintf("[OK] Saved %s signatures to %s",
                      format(nrow(filtered_data), big.mark = ","), output_file))
    }
  }
  
  # After interactive filtering loop is complete and before final return
  if (interactive) {
    message("\n[Tip] You can reuse the following filter settings in non-interactive mode:\n")

    filter_list <- list()
    for (col in key_columns) {
      if (col %in% names(filtered_data)) {
        selected_vals <- sort(unique(filtered_data[[col]]))
        if (length(selected_vals) > 0 && length(selected_vals) < length(unique(sig_info[[col]]))) {
          filter_list[[col]] <- selected_vals
        }
      }
    }
    
    # Pretty-print filter list as R code
    if (length(filter_list) > 0) {
      filter_lines <- character(length(filter_list) + 2)
      filter_lines[1] <- "filters = list("
      for (i in seq_along(filter_list)) {
        colname <- names(filter_list)[i]
        values <- filter_list[[i]]
        quoted_vals <- paste0('"', values, '"', collapse = ", ")

        if (length(values) == 1) {
          filter_lines[i + 1] <- sprintf("  %s = \"%s\"", colname, values)
        } else {
          filter_lines[i + 1] <- sprintf("  %s = c(%s)", colname, quoted_vals)
        }
        if (i < length(filter_list)) {
          filter_lines[i + 1] <- paste0(filter_lines[i + 1], ",")
        }
      }
      filter_lines[length(filter_lines)] <- ")"
      message(paste(filter_lines, collapse = "\n"), "\n")
    }
  }
  
  
  return(filtered_data)
}
