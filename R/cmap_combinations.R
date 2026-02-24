#' CMap Parameter Management and Combination Generation
#'
#' @description Functions to extract experimental parameters from CMap data,
#' create configuration files, and generate combinations for data extraction.
#'
#' @name cmap_combinations
#' @keywords internal
NULL

#' Write CMap options to configuration file
#'
#' @param options_list List with options (times, doses, cells)
#' @param file_path Path to output file
#'
#' @return Path to the written file (invisibly)
#'
#' @keywords internal
write_config_file <- function(options_list, file_path) {
  lines <- character()

  # Write header
  lines <- c(lines, "# CMap experimental parameters configuration")
  lines <- c(lines, "# Generated on", as.character(Sys.time()))
  lines <- c(lines, "# This file lists available options extracted from siginfo_beta.txt")
  lines <- c(lines, "")
  lines <- c(lines, "# HOW TO USE THIS FILE:")
  lines <- c(lines, "# 1. In the [selected] section at the bottom, specify which options to use")
  lines <- c(lines, "# 2. For each parameter, provide comma-separated indices of options you want")
  lines <- c(lines, "#    For example: times = 1,3,5 will select the 1st, 3rd and 5th time points")
  lines <- c(lines, "#    To select all options for a parameter, use 'all' or leave the value empty")
  lines <- c(lines, "")

  # Write each section
  for (section_name in names(options_list)) {
    lines <- c(lines, paste0("[", section_name, "]"))
    items <- options_list[[section_name]]
    for (i in seq_along(items)) {
      lines <- c(lines, paste0(i, " = ", items[i]))
    }
    lines <- c(lines, "")
  }

  # Write default selections section
  lines <- c(lines, "[selected]")
  for (section_name in names(options_list)) {
    # Select the first item in each section as default
    if (length(options_list[[section_name]]) > 0) {
      lines <- c(lines, paste0(section_name, " = 1  # Edit this to select specific ", section_name, " or leave empty for all"))
    } else {
      lines <- c(lines, paste0(section_name, " = "))
    }
  }

  # Write to file
  writeLines(lines, file_path)
  return(invisible(file_path))
}

#' Generate combinations of time points, dosages, and cell lines from a configuration file
#'
#' @param config_file Path to configuration file with time, dose, and cell options
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return Data frame of combinations
#'
#' @examples
#' \dontrun{
#' # First create a configuration file
#' extract_cmap_parameters("databases/siginfo_beta.txt")
#'
#' # Then edit the config file to select desired parameters
#' # Leave sections empty to use all options for that parameter
#' # Now generate combinations from the edited config file
#' combinations <- generate_combinations_from_config("conf/cmap_options.conf")
#' }
#'
#' @export
generate_combinations_from_config <- function(config_file, verbose = TRUE) {

  # Validate input
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }

  # Read configuration
  if (verbose) message("Reading configuration from: ", config_file)
  config <- read_cmap_config(config_file)

  # Extract selected values
  times <- config$selected$times
  doses <- config$selected$doses
  cells <- config$selected$cells

  # Log selections
  if (verbose) {
    message("Selected time points: ", paste(times, collapse = ", "),
            if(length(times) == length(config$all_options$times)) " (all)" else "")
    message("Selected dosages: ", paste(doses, collapse = ", "),
            if(length(doses) == length(config$all_options$doses)) " (all)" else "")
    message("Selected cell lines: ", paste(cells, collapse = ", "),
            if(length(cells) == length(config$all_options$cells)) " (all)" else "")
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

#' Create a template configuration file
#'
#' @param dest_dir Directory to write the template (default: "conf")
#' @param template_name Name of template file (default: "cmap_options_template.conf")
#' @param overwrite Logical; whether to overwrite existing template (default: FALSE)
#'
#' @return Path to the template file (invisibly)
#'
#' @export
create_cmap_config_template <- function(dest_dir = "conf",
                                        template_name = "cmap_options_template.conf",
                                        overwrite = FALSE) {

  # Create the conf directory if it doesn't exist
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
    message("Created configuration directory: ", dest_dir)
  }

  template_path <- file.path(dest_dir, template_name)

  # Check if template exists and we're not overwriting
  if (file.exists(template_path) && !overwrite) {
    message("Template file already exists. Use overwrite=TRUE to replace.")
    return(invisible(template_path))
  }

  # Create a minimal template with common values
  options_list <- list(
    times = c("6 h", "24 h", "96 h"),
    doses = c("0.04 uM", "0.37 uM", "1.11 uM", "3.33 uM", "10 uM"),
    cells = c("A375", "A549", "MCF7", "PC3", "HT29", "VCAP", "HA1E")
  )

  write_config_file(options_list, template_path)
  message("Created template configuration file at: ", template_path)
  return(invisible(template_path))
}

#' Read a CMap configuration file
#'
#' @param config_file Path to configuration file
#' @param return_selected Whether to return selected values (default: TRUE)
#'
#' @return List with parsed configuration (all options and selected options)
#'
#' @export
read_cmap_config <- function(config_file, return_selected = TRUE) {
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }

  # Read the config file
  lines <- readLines(config_file)

  # Initialize variables
  current_section <- NULL
  all_options <- list()
  selected <- list()

  # Parse each line
  for (line in lines) {
    # Skip comments and empty lines
    if (grepl("^\\s*#", line) || grepl("^\\s*$", line)) {
      next
    }

    # Check for section header
    section_match <- regexpr("^\\s*\\[(.*?)\\]\\s*$", line, perl = TRUE)
    if (section_match > 0) {
      current_section <- gsub("^\\s*\\[(.*?)\\]\\s*$", "\\1", line)
      if (current_section != "selected") {
        all_options[[current_section]] <- character()
      }
      next
    }

    # Process key-value pairs
    if (!is.null(current_section) && grepl("=", line)) {
      key_value <- strsplit(line, "=", fixed = TRUE)[[1]]
      key <- trimws(key_value[1])
      value <- trimws(key_value[2])

      # Remove comments
      value <- gsub("#.*$", "", value)
      value <- trimws(value)

      if (current_section == "selected") {
        # For the selected section, handle special "all" value, empty value, or parse indices
        if (value == "all" || value == "") {
          # "all" or empty means select all options for this section
          if (key %in% names(all_options)) {
            selected[[key]] <- seq_along(all_options[[key]])
          }
        } else {
          # Parse comma-separated indices
          selected[[key]] <- as.integer(unlist(strsplit(value, ",\\s*")))
        }
      } else {
        # For other sections, store option values
        all_options[[current_section]][as.integer(key)] <- value
      }
    }
  }

  # If return_selected is TRUE, extract the selected values
  if (return_selected) {
    result <- list(all_options = all_options, selected = list())

    # Extract selected values from each section
    for (section in names(all_options)) {
      if (section %in% names(selected) && length(selected[[section]]) > 0) {
        indices <- selected[[section]]
        valid_indices <- indices[indices <= length(all_options[[section]])]
        if (length(valid_indices) > 0) {
          result$selected[[section]] <- all_options[[section]][valid_indices]
        } else {
          # Empty selection - use all options in this section
          result$selected[[section]] <- all_options[[section]]
        }
      } else {
        # No selection - use all options in this section
        result$selected[[section]] <- all_options[[section]]
      }
    }

    return(result)
  } else {
    # Just return the raw parsed config
    return(list(all_options = all_options, selected_indices = selected))
  }
}
