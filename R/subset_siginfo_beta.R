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
#' \dontrun{
#' # Interactive mode with smart preview of downstream options
#' filtered_siginfo <- subset_siginfo_beta("path/to/siginfo_beta.txt")
#'
#' # Non-interactive mode with pre-defined filters
#' filtered_siginfo <- subset_siginfo_beta(
#'   "path/to/siginfo_beta.txt",
#'   interactive = FALSE,
#'   filters = list(
#'     pert_type = "trt_cp",
#'     pert_itime = c("6 h", "24 h"),
#'     pert_idose = "10 uM",
#'     cell_iname = c("A375", "MCF7")
#'   )
#' )
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
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("INTERACTIVE SIGINFO FILTERING WITH SMART PREVIEW\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat(sprintf("Starting with %d signatures\n", original_count))

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

      cat("\n", paste(rep("=", 80), collapse = ""), "\n")
      cat(sprintf("STEP %d/%d: Filtering by %s\n", col_idx, length(key_columns), col))
      cat(sprintf("Current signatures: %d (%.1f%% of original)\n",
                  current_count, current_proportion))

      # Show what each selection would yield and preview downstream options
      if (show_preview && col_idx < length(key_columns)) {
        cat("\nPREVIEW: Impact on downstream parameters\n")
        cat(paste(rep("-", 80), collapse = ""), "\n")
      }

      cat("\nAvailable values for", col, ":\n")

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

        cat(sprintf("  %d: %s [%d sigs, %.1f%% current, %.1f%% original]",
                    i, val, count, proportion_of_current, proportion_of_original))

        # Show preview of option counts in next parameters
        if (show_preview && col_idx < length(key_columns)) {
          remaining_cols <- key_columns[(col_idx + 1):length(key_columns)]
          preview_info <- ""

          for (next_col in remaining_cols[1:min(2, length(remaining_cols))]) {  # Show max 2 next columns
            if (next_col %in% names(temp_data)) {
              next_count <- length(unique(temp_data[[next_col]]))
              preview_info <- paste0(preview_info, sprintf(" → %s: %d options",
                                                           next_col, next_count))
            }
          }
          if (preview_info != "") {
            cat(preview_info)
          }
        }
        cat("\n")
      }

      # Special preview option: show what combinations would be available
      if (show_preview && col_idx < length(key_columns)) {
        cat(sprintf("\n📋 Type 'preview X' to see detailed downstream options for choice X\n"))
      }

      # Get user selection
      cat("\nSelect values to keep:\n")
      cat("  - Enter comma-separated numbers (e.g., '1,3,5')\n")
      cat("  - Enter 'all' to keep all values\n")
      cat("  - Enter 'skip' to skip this parameter\n")
      if (show_preview) {
        cat("  - Enter 'preview X' to see detailed downstream options for choice X\n")
      }

      while (TRUE) {
        selection <- readline(prompt = "> ")

        if (tolower(selection) == "skip" || selection == "") {
          cat("Skipping", col, "\n")
          break

        } else if (tolower(selection) == "all") {
          cat("Keeping all values for", col, "\n")
          break

        } else if (grepl("^preview\\s+", tolower(selection))) {
          # Detailed preview mode
          preview_input <- gsub("^preview\\s+", "", selection, ignore.case = TRUE)
          preview_idx <- as.integer(preview_input)

          if (!is.na(preview_idx) && preview_idx >= 1 && preview_idx <= length(unique_vals)) {
            preview_val <- unique_vals[preview_idx]
            preview_data <- value_impacts[[preview_val]]$data

            cat("\n", paste(rep("~", 60), collapse = ""), "\n")
            cat(sprintf("DETAILED PREVIEW: If you select '%s'\n", preview_val))
            cat(sprintf("This would retain %d signatures (%.1f%% of original)\n",
                        value_impacts[[preview_val]]$count,
                        value_impacts[[preview_val]]$prop_original))

            # Show what would be available in ALL remaining parameters
            remaining_cols <- key_columns[(col_idx + 1):length(key_columns)]
            for (next_col in remaining_cols) {
              if (next_col %in% names(preview_data)) {
                next_unique <- sort(unique(preview_data[[next_col]]))
                next_counts <- table(preview_data[[next_col]])

                cat(sprintf("\nAvailable %s options (%d):\n", next_col, length(next_unique)))
                for (j in seq_along(next_unique)) {
                  next_val <- next_unique[j]
                  next_count <- next_counts[next_val]
                  next_prop <- round((next_count / nrow(preview_data)) * 100, 1)
                  cat(sprintf("    %s: %d sigs (%.1f%%)\n", next_val, next_count, next_prop))

                  # Limit display to avoid overwhelming output
                  if (j >= 10) {
                    cat(sprintf("    ... and %d more options\n", length(next_unique) - j))
                    break
                  }
                }
              }
            }
            cat(paste(rep("~", 60), collapse = ""), "\n\n")

          } else {
            cat("Invalid preview selection. Please enter a valid number.\n")
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

            cat(sprintf("\nThis selection will retain %d signatures (%.1f%% of original)\n",
                        new_count, new_proportion))
            cat("Selected values: ", paste(selected_vals, collapse = ", "), "\n")

            # Show combined preview of next option counts
            if (show_preview && col_idx < length(key_columns)) {
              remaining_cols <- key_columns[(col_idx + 1):length(key_columns)]
              next_col <- remaining_cols[1]
              if (next_col %in% names(combined_data)) {
                next_count <- length(unique(combined_data[[next_col]]))
                cat(sprintf("Next parameter (%s) will have %d options\n",
                            next_col, next_count))
              }
            }

            # Confirm selection
            confirm <- readline(prompt = "Confirm this selection? (y/n): ")

            if (tolower(confirm) == "y" || tolower(confirm) == "yes") {
              filtered_data <- combined_data

              if (verbose) {
                message(sprintf("✓ Applied filter on %s: %d values selected, %d rows remaining (%.1f%%)",
                                col, length(selected_vals), nrow(filtered_data), new_proportion))
              }
              break
            } else {
              cat("Selection cancelled. You can make a new selection.\n")
              # Continue the loop to allow retry
            }

          } else {
            cat("Invalid selection. Please enter valid numbers, 'all', 'skip', or 'preview X'.\n")
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
            message(sprintf("Filtered %s: %d → %d signatures (%.1f%% of original)",
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
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("FILTERING COMPLETE!\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat(sprintf("Original signatures: %s\n", format(original_count, big.mark = ",")))
    cat(sprintf("Final signatures: %s (%.1f%% retained)\n",
                format(final_count, big.mark = ","), final_proportion))
    cat(sprintf("Signatures removed: %s (%.1f%% reduction)\n",
                format(original_count - final_count, big.mark = ","), 100 - final_proportion))

    # Show final distribution
    cat("\n🔍 Final distribution by parameter:\n")
    for (col in key_columns) {
      if (col %in% names(filtered_data) && nrow(filtered_data) > 0) {
        val_counts <- sort(table(filtered_data[[col]]), decreasing = TRUE)
        cat(sprintf("\n%s (%d unique values):\n", col, length(val_counts)))

        # Show top values
        show_n <- min(length(val_counts), 8)
        for (i in 1:show_n) {
          val_prop <- round((val_counts[i] / final_count) * 100, 1)
          cat(sprintf("    %s: %s (%.1f%%)\n",
                      names(val_counts)[i],
                      format(val_counts[i], big.mark = ","),
                      val_prop))
        }
        if (length(val_counts) > 8) {
          remaining_total <- sum(val_counts[9:length(val_counts)])
          cat(sprintf("    ... and %d more values: %s signatures\n",
                      length(val_counts) - 8,
                      format(remaining_total, big.mark = ",")))
        }
      }
    }
  }

  # Save filtered data if output file is specified
  if (!is.null(output_file)) {
    if (verbose) message("\n💾 Saving filtered data to: ", output_file)

    # Create directory if it doesn't exist
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir) && output_dir != ".") {
      dir.create(output_dir, recursive = TRUE)
    }

    # Write the file
    utils::write.table(filtered_data, file = output_file,
                       sep = "\t", row.names = FALSE, quote = FALSE)

    if (verbose) {
      message(sprintf("✅ Saved %s signatures to %s",
                      format(nrow(filtered_data), big.mark = ","), output_file))
    }
  }

  # After interactive filtering loop is complete and before final return
  if (interactive) {
    cat("\n📋 You can reuse the following filter settings in non-interactive mode:\n\n")

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
    cat("filters = list(\n")
    for (i in seq_along(filter_list)) {
      colname <- names(filter_list)[i]
      values <- filter_list[[i]]
      quoted_vals <- paste0('"', values, '"', collapse = ", ")

      if (length(values) == 1) {
        cat(sprintf("  %s = \"%s\"", colname, values))
      } else {
        cat(sprintf("  %s = c(%s)", colname, quoted_vals))
      }

      if (i < length(filter_list)) {
        cat(",\n")
      } else {
        cat("\n")
      }
    }
    cat(")\n\n")
  }


  return(filtered_data)
}
