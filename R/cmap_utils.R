#' Create an interactive configuration using the console
#'
#' Provides a text-based interactive interface in the console that allows users to
#' select CMap experimental parameters (time points, dosages, cell lines).
#' This is a fallback option when the GUI version is not available.
#'
#' @param siginfo_file Path to siginfo_beta.txt file
#' @param config_dir Directory to save the configuration file (default: "conf")
#' @param config_filename Name of the configuration file (default: "cmap_options.conf")
#' @param filter_quality Logical; whether to filter for high-quality compounds only (default: TRUE)
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return A list with selected parameters (times, doses, cells) and the config file path
#'
#' @examples
#' ex_sig <- data.frame(
#'   pert_itime = c("6 h", "24 h"),
#'   pert_idose = c("1 uM", "10 uM"),
#'   cell_iname = c("A375", "MCF7"),
#'   pert_type = c("trt_cp", "trt_cp"),
#'   is_hiq = c(1, 1),
#'   stringsAsFactors = FALSE
#' )
#' ex_sig_file <- tempfile(fileext = ".txt")
#' write.table(ex_sig, ex_sig_file, sep = "\t", row.names = FALSE, quote = FALSE)
#' sel <- create_console_config(ex_sig_file, config_dir = tempdir(), verbose = FALSE)
#' names(sel)
#'
#' \donttest{
#' # Create interactive console configuration
#' siginfo_file <- system.file("extdata", "example_siginfo.txt", package = "CONCERTDR")
#' selections <- create_console_config(siginfo_file)
#'
#' # Use selections for further analysis
#' combinations <- generate_combinations_from_selections(selections)
#' }
#'
#' @export
create_console_config <- function(siginfo_file,
                                  config_dir = "conf",
                                  config_filename = "cmap_options.conf",
                                  filter_quality = TRUE,
                                  verbose = TRUE) {

  # Extract parameters from siginfo file first
  if (verbose) message("Extracting parameters from ", siginfo_file)
  params <- extract_cmap_parameters(siginfo_file,
                                    write_config = FALSE,
                                    filter_quality = filter_quality,
                                    verbose = verbose)

  # Check if parameters were extracted successfully
  if (is.null(params$times) || is.null(params$doses) || is.null(params$cells)) {
    stop("Failed to extract parameters from siginfo file")
  }

  options_list <- list(
    times = params$times,
    doses = params$doses,
    cells = params$cells
  )

  # Non-interactive sessions default to selecting all options
  if (!interactive()) {
    if (verbose) {
      message("Non-interactive session detected; selecting all available options.")
    }

    selected <- options_list
  } else {
    # Console helper for one section
    select_section <- function(section_name, items) {
      message("\n", paste(rep("=", 72), collapse = ""))
      message(sprintf("Select %s (n = %d)", section_name, length(items)))
      message(paste(rep("-", 72), collapse = ""))

      for (i in seq_along(items)) {
        message(sprintf("%4d. %s", i, items[i]))
      }

      message("\nInput options:")
      message("  - Comma-separated indices, e.g. 1,3,5")
      message("  - 'all' to select all")
      message("  - 'none' to keep this section empty")

      repeat {
        answer <- trimws(readline(prompt = paste0("Selected ", section_name, ": ")))

        if (answer == "" || tolower(answer) == "all") {
          return(items)
        }

        if (tolower(answer) == "none") {
          return(character())
        }

        idx <- as.integer(unlist(strsplit(answer, ",")))
        idx <- idx[!is.na(idx)]
        idx <- unique(idx[idx >= 1 & idx <= length(items)])

        if (length(idx) == 0) {
          message("Invalid selection. Please try again.")
          next
        }

        selected_items <- items[idx]
        message("Selected: ", paste(selected_items, collapse = ", "))

        confirm <- tolower(trimws(readline(prompt = "Confirm selection? (y/n): ")))
        if (confirm %in% c("y", "yes")) {
          return(selected_items)
        }
      }
    }

    selected <- list(
      times = select_section("times", options_list$times),
      doses = select_section("doses", options_list$doses),
      cells = select_section("cells", options_list$cells)
    )
  }

  # Create configuration directory if it doesn't exist
  if (!dir.exists(config_dir)) {
    dir.create(config_dir, recursive = TRUE)
    if (verbose) message("Created configuration directory: ", config_dir)
  }

  # Convert selected values to index representation used by config file
  selected_indices <- list(
    times = which(options_list$times %in% selected$times),
    doses = which(options_list$doses %in% selected$doses),
    cells = which(options_list$cells %in% selected$cells)
  )

  # Write configuration
  config_file <- file.path(config_dir, config_filename)

  lines <- character()
  lines <- c(lines, "# CMap experimental parameters configuration")
  lines <- c(lines, "# Generated on", as.character(Sys.time()))
  lines <- c(lines, "# This file was created using the console configuration tool")
  lines <- c(lines, "")

  for (section_name in names(options_list)) {
    lines <- c(lines, paste0("[", section_name, "]"))
    items <- options_list[[section_name]]
    for (i in seq_along(items)) {
      lines <- c(lines, paste0(i, " = ", items[i]))
    }
    lines <- c(lines, "")
  }

  lines <- c(lines, "[selected]")
  for (section_name in names(selected_indices)) {
    idx <- selected_indices[[section_name]]
    lines <- c(lines, paste0(section_name, " = ", paste(idx, collapse = ",")))
  }

  writeLines(lines, config_file)

  if (verbose) {
    message("Configuration saved to: ", config_file)
    message("Selected ", length(selected$times), " time points, ",
            length(selected$doses), " doses, and ",
            length(selected$cells), " cell lines.")
  }

  return(list(
    selected = selected,
    all_options = options_list,
    config_file = config_file
  ))
}

#' Generate combinations from selected parameters
#'
#' Creates a data frame of combinations based on the selections from the interactive config functions
#'
#' @param selections The list returned by create_console_config()
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return Data frame of combinations
#'
#' @examples
#' sel <- list(
#'   selected = list(times = c("6 h"), doses = c("10 uM"), cells = c("A375", "MCF7")),
#'   all_options = list(times = c("6 h", "24 h"), doses = c("10 uM"), cells = c("A375", "MCF7"))
#' )
#' generate_combinations_from_selections(sel, verbose = FALSE)
#'
#' \donttest{
#' siginfo_file <- system.file("extdata", "example_siginfo.txt", package = "CONCERTDR")
#' selections <- create_console_config(siginfo_file)
#'
#' # Generate combinations
#' combinations <- generate_combinations_from_selections(selections)
#' }
#'
#' @export
generate_combinations_from_selections <- function(selections, verbose = TRUE) {
  # Check input
  if (is.null(selections)) {
    stop("No selections provided. Run create_console_config() first.")
  }

  if (!is.list(selections) || is.null(selections$selected)) {
    stop("Invalid selections object. Must be the result of an interactive configuration function.")
  }

  # Extract selected values
  times <- selections$selected$times
  doses <- selections$selected$doses
  cells <- selections$selected$cells

  # Check if any selections were made
  if (length(times) == 0 || length(doses) == 0 || length(cells) == 0) {
    warning("One or more parameter types have no selections.")
    # Use all options for empty selections
    if (length(times) == 0) times <- selections$all_options$times
    if (length(doses) == 0) doses <- selections$all_options$doses
    if (length(cells) == 0) cells <- selections$all_options$cells
  }

  # Log selections
  if (verbose) {
    message("Selected time points: ", paste(times, collapse = ", "))
    message("Selected dosages: ", paste(doses, collapse = ", "))
    message("Selected cell lines: ", paste(cells, collapse = ", "))
  }

  # Generate combinations
  combinations <- expand.grid(
    itime = times,
    idose = doses,
    cell = cells,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message("Generated ", nrow(combinations), " combinations")
  }

  return(combinations)
}

#' Interactive CMap configuration and workflow
#'
#' Provides a user-friendly interface to configure and run the CMap workflow.
#' Uses the console configuration workflow for reliability and speed.
#'
#' @param siginfo_file Path to siginfo_beta.txt file (default: NULL to prompt user)
#' @param geneinfo_file Path to geneinfo_beta.txt file (default: NULL to prompt user)
#' @param gctx_file Path to GCTX file (default: NULL to prompt user)
#' @param config_dir Directory to save configuration files (default: "conf")
#' @param output_dir Directory to save output files (default: "output")
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return A list with workflow results
#'
#' @examples
#' is.function(interactive_cmap_setup)
#'
#' \donttest{
#' results <- interactive_cmap_setup()
#' }
#'
#' @export
interactive_cmap_setup <- function(siginfo_file = NULL,
                                   geneinfo_file = NULL,
                                   gctx_file = NULL,
                                   config_dir = "conf",
                                   output_dir = "output",
                                   verbose = TRUE) {

  # Check for prerequisites if not provided
  if (is.null(siginfo_file)) {
    siginfo_file <- readline(prompt = "Enter path to siginfo_beta.txt file: ")
  }

  if (is.null(geneinfo_file)) {
    geneinfo_file <- readline(prompt = "Enter path to geneinfo_beta.txt file: ")
  }

  if (is.null(gctx_file)) {
    gctx_file <- readline(prompt = "Enter path to GCTX file: ")
  }

  selections <- create_console_config(
    siginfo_file = siginfo_file,
    config_dir = config_dir,
    verbose = verbose
  )

  # If selection was cancelled, return NULL
  if (is.null(selections)) {
    message("Configuration was cancelled.")
    return(NULL)
  }

  # Generate combinations
  combinations <- generate_combinations_from_selections(selections, verbose = verbose)

  # Ask if user wants to process combinations now
  process_now <- readline(prompt = "Process combinations now? (y/n): ")

  if (tolower(process_now) == "y") {
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      if (verbose) message("Created output directory: ", output_dir)
    }

    # Process combinations
    results <- process_combinations(
      combinations = combinations,
      output_dir = output_dir,
      geneinfo_file = geneinfo_file,
      siginfo_file = siginfo_file,
      gctx_file = gctx_file,
      verbose = verbose
    )

    return(list(
      selections = selections,
      combinations = combinations,
      results = results
    ))
  } else {
    message("Combinations generated but not processed.")
    return(list(
      selections = selections,
      combinations = combinations
    ))
  }
}
