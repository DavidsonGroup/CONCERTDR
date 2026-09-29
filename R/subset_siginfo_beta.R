# Columns offered for filtering, in prompt order.
.siginfo_key_columns <- c("pert_type", "pert_itime", "pert_idose", "cell_iname")

#' Read a siginfo file
#' @param siginfo_file Path to the siginfo table.
#' @param verbose Logical; print progress messages.
#' @return Data frame.
#' @keywords internal
.read_siginfo <- function(siginfo_file, verbose) {
  if (!file.exists(siginfo_file)) stop("Siginfo file not found: ", siginfo_file)
  if (verbose) message("Reading siginfo file: ", siginfo_file)
  sig_info <- .read_cmap_table(siginfo_file, "siginfo_file")
  if (verbose) {
    message(sprintf("Loaded siginfo data with %d rows and %d columns",
                    nrow(sig_info), ncol(sig_info)))
  }
  sig_info
}

#' Read one line of user input
#'
#' Thin wrapper around \code{readline} so that the prompt loop can be mocked.
#' @param prompt Prompt string.
#' @return Character scalar.
#' @keywords internal
.siginfo_ask <- function(prompt) readline(prompt = prompt)

#' Percentage rounded to one decimal
#' @param n,total Numerator and denominator.
#' @return Numeric scalar.
#' @keywords internal
.siginfo_pct <- function(n, total) round((n / total) * 100, 1)

#' Columns that follow position \code{col_idx} in \code{key_columns}
#' @param key_columns Character vector of filter columns.
#' @param col_idx Current position.
#' @return Character vector (possibly empty).
#' @keywords internal
.remaining_columns <- function(key_columns, col_idx) {
  if (col_idx < length(key_columns)) key_columns[(col_idx + 1):length(key_columns)]
  else character()
}

#' Print the available values of a column with their impact
#' @param filtered_data Current data frame.
#' @param col Column being filtered.
#' @param unique_vals Sorted unique values of \code{col}.
#' @param remaining_cols Columns after \code{col}.
#' @param original_count Row count before any filtering.
#' @param show_preview Logical; append downstream option counts.
#' @return Named list (one element per value) with \code{data}, \code{count},
#'   \code{prop_current} and \code{prop_original}.
#' @keywords internal
.list_column_values <- function(filtered_data, col, unique_vals, remaining_cols,
                                original_count, show_preview) {
  current_count <- nrow(filtered_data)
  value_impacts <- list()
  for (i in seq_along(unique_vals)) {
    val <- unique_vals[i]
    temp_data <- filtered_data[filtered_data[[col]] == val, ]
    count <- nrow(temp_data)
    prop_current <- .siginfo_pct(count, current_count)
    prop_original <- .siginfo_pct(count, original_count)
    value_impacts[[val]] <- list(data = temp_data, count = count,
                                 prop_current = prop_current,
                                 prop_original = prop_original)

    preview_line <- sprintf("  %d: %s [%d sigs, %.1f%% current, %.1f%% original]",
                            i, val, count, prop_current, prop_original)
    if (show_preview && length(remaining_cols) > 0) {
      for (next_col in remaining_cols[seq_len(min(2, length(remaining_cols)))]) {
        if (next_col %in% names(temp_data)) {
          preview_line <- paste0(preview_line, sprintf(
            " -> %s: %d options", next_col, length(unique(temp_data[[next_col]]))))
        }
      }
    }
    message(preview_line)
  }
  value_impacts
}

#' Print the detailed downstream preview for one candidate value
#' @param preview_val The value being previewed.
#' @param impact Its entry of the list from \code{.list_column_values()}.
#' @param remaining_cols Columns after the current one.
#' @return NULL, invisibly.
#' @keywords internal
.show_value_preview <- function(preview_val, impact, remaining_cols) {
  preview_data <- impact$data
  message("")
  .message_rule("~", 60)
  message(sprintf("DETAILED PREVIEW: If you select '%s'", preview_val))
  message(sprintf("This would retain %d signatures (%.1f%% of original)",
                  impact$count, impact$prop_original))

  for (next_col in remaining_cols) {
    if (!next_col %in% names(preview_data)) next
    next_unique <- sort(unique(preview_data[[next_col]]))
    next_counts <- table(preview_data[[next_col]])
    message(sprintf("\nAvailable %s options (%d):", next_col, length(next_unique)))
    for (j in seq_along(next_unique)) {
      next_val <- next_unique[j]
      next_count <- next_counts[next_val]
      message(sprintf("    %s: %d sigs (%.1f%%)", next_val, next_count,
                      round((next_count / nrow(preview_data)) * 100, 1)))
      if (j >= 10) {
        message(sprintf("    ... and %d more options", length(next_unique) - j))
        break
      }
    }
  }
  .message_rule("~", 60)
  message("")
  invisible(NULL)
}

#' Interactively filter the data on one column
#' @param filtered_data Current data frame.
#' @param col Column to filter.
#' @param col_idx,key_columns Position of \code{col} within \code{key_columns}.
#' @param original_count Row count before any filtering.
#' @param verbose,show_preview Logical flags.
#' @return The (possibly reduced) data frame.
#' @keywords internal
.prompt_column_filter <- function(filtered_data, col, col_idx, key_columns,
                                  original_count, verbose, show_preview) {
  unique_vals <- sort(unique(filtered_data[[col]]))
  if (length(unique_vals) == 0) return(filtered_data)

  current_count <- nrow(filtered_data)
  remaining_cols <- .remaining_columns(key_columns, col_idx)
  has_next <- length(remaining_cols) > 0

  message("")
  .message_rule("=")
  message(sprintf("STEP %d/%d: Filtering by %s", col_idx, length(key_columns), col))
  message(sprintf("Current signatures: %d (%.1f%% of original)",
                  current_count, .siginfo_pct(current_count, original_count)))
  if (show_preview && has_next) {
    message("\nPREVIEW: Impact on downstream parameters")
    .message_rule("-")
  }

  message("\nAvailable values for ", col, ":")
  value_impacts <- .list_column_values(filtered_data, col, unique_vals,
                                       remaining_cols, original_count,
                                       show_preview)
  if (show_preview && has_next) {
    message("\n[Tip] Type 'preview X' to see detailed downstream options for choice X")
  }

  message("\nSelect values to keep:")
  message("  - Enter comma-separated numbers (e.g., '1,3,5')")
  message("  - Enter 'all' to keep all values")
  message("  - Enter 'skip' to skip this parameter")
  if (show_preview) {
    message("  - Enter 'preview X' to see detailed downstream options for choice X")
  }

  while (TRUE) {
    selection <- .siginfo_ask("> ")

    if (tolower(selection) == "skip" || selection == "") {
      message("Skipping ", col)
      break
    }
    if (tolower(selection) == "all") {
      message("Keeping all values for ", col)
      break
    }
    if (grepl("^preview\\s+", tolower(selection))) {
      preview_idx <- as.integer(gsub("^preview\\s+", "", selection, ignore.case = TRUE))
      if (!is.na(preview_idx) && preview_idx >= 1 && preview_idx <= length(unique_vals)) {
        preview_val <- unique_vals[preview_idx]
        .show_value_preview(preview_val, value_impacts[[preview_val]], remaining_cols)
      } else {
        message("Invalid preview selection. Please enter a valid number.")
      }
      next
    }

    indices <- as.integer(unlist(strsplit(selection, ",")))
    indices <- indices[!is.na(indices)]
    valid_indices <- indices[indices >= 1 & indices <= length(unique_vals)]
    if (length(valid_indices) == 0) {
      message("Invalid selection. Please enter valid numbers, 'all', 'skip', or 'preview X'.")
      next
    }

    selected_vals <- unique_vals[valid_indices]
    combined_data <- filtered_data[filtered_data[[col]] %in% selected_vals, ]
    new_count <- nrow(combined_data)
    new_proportion <- .siginfo_pct(new_count, original_count)

    message(sprintf("\nThis selection will retain %d signatures (%.1f%% of original)",
                    new_count, new_proportion))
    message("Selected values: ", paste(selected_vals, collapse = ", "))
    if (show_preview && has_next && remaining_cols[1] %in% names(combined_data)) {
      message(sprintf("Next parameter (%s) will have %d options", remaining_cols[1],
                      length(unique(combined_data[[remaining_cols[1]]]))))
    }

    confirm <- .siginfo_ask("Confirm this selection? (y/n): ")
    if (tolower(confirm) %in% c("y", "yes")) {
      filtered_data <- combined_data
      if (verbose) {
        message(sprintf("[OK] Applied filter on %s: %d values selected, %d rows remaining (%.1f%%)",
                        col, length(selected_vals), nrow(filtered_data), new_proportion))
      }
      break
    }
    message("Selection cancelled. You can make a new selection.")
  }
  filtered_data
}

#' Run the interactive filtering session over all key columns
#' @param filtered_data Data frame to filter.
#' @param key_columns Columns to filter on, in order.
#' @param original_count Row count before any filtering.
#' @param verbose,show_preview Logical flags.
#' @return The filtered data frame.
#' @keywords internal
.filter_interactively <- function(filtered_data, key_columns, original_count,
                                  verbose, show_preview) {
  message("")
  .message_rule("=")
  message("INTERACTIVE SIGINFO FILTERING WITH SMART PREVIEW")
  .message_rule("=")
  message(sprintf("Starting with %d signatures", original_count))
  for (col_idx in seq_along(key_columns)) {
    filtered_data <- .prompt_column_filter(
      filtered_data, key_columns[col_idx], col_idx, key_columns,
      original_count, verbose, show_preview
    )
  }
  filtered_data
}

#' Apply pre-defined filters (non-interactive mode)
#' @param filtered_data Data frame to filter.
#' @param filters Named list; names are columns, values the values to keep.
#'   Unknown columns are skipped.
#' @param original_count Row count before any filtering.
#' @param verbose Logical; print one line per applied filter.
#' @return The filtered data frame.
#' @keywords internal
.apply_filters <- function(filtered_data, filters, original_count, verbose) {
  for (col in names(filters)) {
    if (!col %in% names(filtered_data)) next
    old_count <- nrow(filtered_data)
    filtered_data <- filtered_data[filtered_data[[col]] %in% filters[[col]], ]
    new_count <- nrow(filtered_data)
    if (verbose) {
      message(sprintf("Filtered %s: %d -> %d signatures (%.1f%% of original)",
                      col, old_count, new_count, .siginfo_pct(new_count, original_count)))
    }
  }
  filtered_data
}

#' Print the final filtering summary
#' @param filtered_data Filtered data frame.
#' @param key_columns Columns to summarise.
#' @param original_count Row count before any filtering.
#' @return NULL, invisibly.
#' @keywords internal
.print_filter_summary <- function(filtered_data, key_columns, original_count) {
  final_count <- nrow(filtered_data)
  final_proportion <- .siginfo_pct(final_count, original_count)

  message("")
  .message_rule("=")
  message("FILTERING COMPLETE!")
  .message_rule("=")
  message(sprintf("Original signatures: %s", format(original_count, big.mark = ",")))
  message(sprintf("Final signatures: %s (%.1f%% retained)",
                  format(final_count, big.mark = ","), final_proportion))
  message(sprintf("Signatures removed: %s (%.1f%% reduction)",
                  format(original_count - final_count, big.mark = ","),
                  100 - final_proportion))

  message("\n[Summary] Final distribution by parameter:")
  for (col in key_columns) {
    if (!col %in% names(filtered_data) || nrow(filtered_data) == 0) next
    val_counts <- sort(table(filtered_data[[col]]), decreasing = TRUE)
    message(sprintf("\n%s (%d unique values):", col, length(val_counts)))
    for (i in seq_len(min(length(val_counts), 8))) {
      message(sprintf("    %s: %s (%.1f%%)", names(val_counts)[i],
                      format(val_counts[i], big.mark = ","),
                      round((val_counts[i] / final_count) * 100, 1)))
    }
    if (length(val_counts) > 8) {
      message(sprintf("    ... and %d more values: %s signatures",
                      length(val_counts) - 8,
                      format(sum(val_counts[9:length(val_counts)]), big.mark = ",")))
    }
  }
  invisible(NULL)
}

#' Write the filtered siginfo table
#' @param filtered_data Data frame to save.
#' @param output_file Destination path; missing directories are created.
#' @param verbose Logical; print progress messages.
#' @return NULL, invisibly.
#' @keywords internal
.save_siginfo <- function(filtered_data, output_file, verbose) {
  if (verbose) message("\n[Saving] Saving filtered data to: ", output_file)
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE)
  }
  utils::write.table(filtered_data, file = output_file,
                     sep = "\t", row.names = FALSE, quote = FALSE)
  if (verbose) {
    message(sprintf("[OK] Saved %s signatures to %s",
                    format(nrow(filtered_data), big.mark = ","), output_file))
  }
  invisible(NULL)
}

#' Print the filters that reproduce the interactive selection
#' @param filtered_data Filtered data frame.
#' @param sig_info Unfiltered data frame.
#' @param key_columns Filter columns.
#' @return NULL, invisibly.
#' @keywords internal
.print_reuse_filters <- function(filtered_data, sig_info, key_columns) {
  message("\n[Tip] You can reuse the following filter settings in non-interactive mode:\n")

  filter_list <- list()
  for (col in key_columns) {
    if (!col %in% names(filtered_data)) next
    selected_vals <- sort(unique(filtered_data[[col]]))
    if (length(selected_vals) > 0 &&
        length(selected_vals) < length(unique(sig_info[[col]]))) {
      filter_list[[col]] <- selected_vals
    }
  }
  if (length(filter_list) == 0) return(invisible(NULL))

  entries <- vapply(names(filter_list), function(nm) {
    values <- filter_list[[nm]]
    if (length(values) == 1) {
      sprintf("  %s = \"%s\"", nm, values)
    } else {
      sprintf("  %s = c(%s)", nm, paste0('"', values, '"', collapse = ", "))
    }
  }, character(1))
  message(paste(c("filters = list(", paste0(entries, c(rep(",", length(entries) - 1), "")), ")"),
                collapse = "\n"), "\n")
  invisible(NULL)
}

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
#'   cell_iname = c("A375", "MCF7", "A375")
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
  sig_info <- .read_siginfo(siginfo_file, verbose)
  original_count <- nrow(sig_info)

  key_columns <- .siginfo_key_columns
  missing_cols <- setdiff(key_columns, names(sig_info))
  if (length(missing_cols) > 0) {
    warning("Missing columns: ", paste(missing_cols, collapse = ", "))
    key_columns <- intersect(key_columns, names(sig_info))
  }

  filtered_data <- if (interactive) {
    .filter_interactively(sig_info, key_columns, original_count, verbose, show_preview)
  } else {
    .apply_filters(sig_info, filters, original_count, verbose)
  }

  if (verbose) .print_filter_summary(filtered_data, key_columns, original_count)
  if (!is.null(output_file)) .save_siginfo(filtered_data, output_file, verbose)
  if (interactive) .print_reuse_filters(filtered_data, sig_info, key_columns)

  filtered_data
}
