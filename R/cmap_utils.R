#' Create an interactive configuration window
#'
#' Opens a GUI window that allows users to select CMap experimental parameters
#' (time points, dosages, cell lines) interactively. The selections are saved
#' to a configuration file and returned as an R object.
#'
#' This function requires the tcltk package, which is usually included with R.
#' If the GUI version fails, you can use create_console_config() as an alternative.
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
#' \dontrun{
#' # Create interactive configuration
#' selections <- create_interactive_config("databases/siginfo_beta.txt")
#'
#' # Use selections for further analysis
#' combinations <- generate_combinations_from_selections(selections)
#' }
#'
#' @export
create_interactive_config <- function(siginfo_file,
                                      config_dir = "conf",
                                      config_filename = "cmap_options.conf",
                                      filter_quality = TRUE,
                                      verbose = TRUE) {

  # Check if tcltk is available
  if (!requireNamespace("tcltk", quietly = TRUE)) {
    stop("The 'tcltk' package is required for the interactive GUI.\n",
         "It should be included with R, but seems to be unavailable.\n",
         "You may need to reinstall R with Tk support or install Tk separately.\n",
         "For console-based configuration, use create_console_config() instead.")
  }

  # Check if Tk is actually working
  tryCatch({
    tcltk::tclRequire("Tk")
  }, error = function(e) {
    stop("Tk is not properly installed or configured.\n",
         "You may need to install Tk on your system:\n",
         " - Ubuntu/Debian: sudo apt-get install tk-dev\n",
         " - RHEL/Fedora: sudo yum install tk-devel\n",
         " - macOS: brew install tcl-tk\n",
         "For console-based configuration, use create_console_config() instead.")
  })

  # Extract parameters from siginfo file first (reusing existing functionality)
  if (verbose) message("Extracting parameters from ", siginfo_file)
  params <- extract_cmap_parameters(siginfo_file,
                                    write_config = FALSE,
                                    filter_quality = filter_quality,
                                    verbose = verbose)

  # Check if parameters were extracted successfully
  if (is.null(params$times) || is.null(params$doses) || is.null(params$cells)) {
    stop("Failed to extract parameters from siginfo file")
  }

  # Initialize selected parameters (empty by default)
  selected <- list(
    times = character(),
    doses = character(),
    cells = character()
  )

  # Create tcltk window
  if (verbose) message("Creating interactive selection window...")

  tt <- tcltk::tktoplevel()
  tcltk::tktitle(tt) <- "CMap Parameter Selection"

  # Set window size and position
  tcltk::tkwm.geometry(tt, "800x600")
  tcltk::tkwm.resizable(tt, 1, 1)

  # Create frames for each parameter type
  frame_header <- tcltk::tkframe(tt, borderwidth = 2)
  frame_times <- tcltk::tkframe(tt, borderwidth = 2, relief = "groove")
  frame_doses <- tcltk::tkframe(tt, borderwidth = 2, relief = "groove")
  frame_cells <- tcltk::tkframe(tt, borderwidth = 2, relief = "groove")
  frame_buttons <- tcltk::tkframe(tt, borderwidth = 2)

  # Add header with instructions
  header_label <- tcltk::tklabel(frame_header,
                                 text = "Select the parameters to include in your analysis.\nCheck the boxes next to the parameters you want to use.\nUse the 'Select All' buttons to quickly select all options in a category.",
                                 justify = "left")
  tcltk::tkpack(header_label, side = "top", fill = "x", padx = 5, pady = 5)

  # Helper function to create a section with a title, checkboxes, and a "Select All" button
  create_selection_section <- function(frame, title, items, selected_items) {
    # Create and pack the title
    title_label <- tcltk::tklabel(frame, text = title, font = "helvetica 12 bold")
    tcltk::tkpack(title_label, side = "top", anchor = "w", padx = 5, pady = 2)

    # Create a sub-frame for the checkboxes with a scrollbar
    checkbox_frame <- tcltk::tkframe(frame)
    scrollbar <- tcltk::tkscrollbar(checkbox_frame, command = function(...) tcltk::tkyview(checkbox_canvas, ...))
    checkbox_canvas <- tcltk::tkcanvas(checkbox_frame, yscrollcommand = function(...) tcltk::tkset(scrollbar, ...),
                                       width = 750, height = 120, background = "white")

    # Pack the scrollbar and canvas
    tcltk::tkpack(scrollbar, side = "right", fill = "y")
    tcltk::tkpack(checkbox_canvas, side = "left", fill = "both", expand = TRUE)

    # Create a frame inside the canvas to hold the checkboxes
    content_frame <- tcltk::tkframe(checkbox_canvas, background = "white")
    content_window <- tcltk::tkcreate(checkbox_canvas, "window", 0, 0, anchor = "nw", window = content_frame)

    # Track variable states for checkboxes
    check_vars <- list()

    # Create checkboxes
    n_cols <- 4  # Number of columns for checkbox layout
    for (i in seq_along(items)) {
      row <- (i - 1) %/% n_cols
      col <- (i - 1) %% n_cols

      var <- tcltk::tclVar("0")  # Initialize as unchecked
      check_vars[[i]] <- var

      # Create the checkbox
      chk <- tcltk::tkcheckbutton(content_frame,
                                  text = items[i],
                                  variable = var,
                                  onvalue = "1",
                                  offvalue = "0",
                                  padx = 5)

      # Grid layout for better organization
      tcltk::tkgrid(chk, row = row, column = col, sticky = "w", padx = 5, pady = 2)
    }

    # Configure the canvas scroll region
    tcltk::tkupdate()  # Force update to get widget dimensions
    bbox <- as.numeric(tcltk::tkwinfo("reqwidth", content_frame))
    bbox2 <- as.numeric(tcltk::tkwinfo("reqheight", content_frame))
    tcltk::tkconfigure(checkbox_canvas, scrollregion = c(0, 0, bbox, bbox2))

    # Pack the checkbox frame
    tcltk::tkpack(checkbox_frame, side = "top", fill = "both", expand = TRUE, padx = 5, pady = 5)

    # Create a frame for buttons
    button_frame <- tcltk::tkframe(frame)
    tcltk::tkpack(button_frame, side = "top", fill = "x", padx = 5, pady = 5)

    # Create "Select All" button
    select_all_button <- tcltk::tkbutton(button_frame, text = "Select All",
                                         command = function() {
                                           for (var in check_vars) {
                                             tcltk::tclvalue(var) <- "1"
                                           }
                                         })
    tcltk::tkpack(select_all_button, side = "left", padx = 5)

    # Create "Deselect All" button
    deselect_all_button <- tcltk::tkbutton(button_frame, text = "Deselect All",
                                           command = function() {
                                             for (var in check_vars) {
                                               tcltk::tclvalue(var) <- "0"
                                             }
                                           })
    tcltk::tkpack(deselect_all_button, side = "left", padx = 5)

    return(list(frame = frame, check_vars = check_vars, items = items))
  }

  # Create sections for each parameter type
  times_section <- create_selection_section(frame_times, "Time Points", params$times, selected$times)
  doses_section <- create_selection_section(frame_doses, "Dosages", params$doses, selected$doses)
  cells_section <- create_selection_section(frame_cells, "Cell Lines", params$cells, selected$cells)

  # Set window title with parameter counts
  tcltk::tktitle(tt) <- sprintf("CMap Parameter Selection (%d times, %d doses, %d cell lines)",
                                length(params$times), length(params$doses), length(params$cells))

  # Variable to track if the user confirmed or cancelled
  result_var <- tcltk::tclVar("cancel")

  # Buttons for Save, Cancel
  save_button <- tcltk::tkbutton(frame_buttons, text = "Save Selections",
                                 command = function() {
                                   tcltk::tclvalue(result_var) <- "save"
                                   tcltk::tkdestroy(tt)
                                 })

  cancel_button <- tcltk::tkbutton(frame_buttons, text = "Cancel",
                                   command = function() {
                                     tcltk::tclvalue(result_var) <- "cancel"
                                     tcltk::tkdestroy(tt)
                                   })

  tcltk::tkpack(save_button, side = "left", padx = 10)
  tcltk::tkpack(cancel_button, side = "left", padx = 10)

  # Pack frames
  tcltk::tkpack(frame_header, side = "top", fill = "x", padx = 5, pady = 5)
  tcltk::tkpack(frame_times, side = "top", fill = "both", expand = TRUE, padx = 5, pady = 5)
  tcltk::tkpack(frame_doses, side = "top", fill = "both", expand = TRUE, padx = 5, pady = 5)
  tcltk::tkpack(frame_cells, side = "top", fill = "both", expand = TRUE, padx = 5, pady = 5)
  tcltk::tkpack(frame_buttons, side = "bottom", fill = "x", padx = 10, pady = 10)

  # Wait for user to interact with window and get result
  tcltk::tkwait.variable(result_var)
  result <- tcltk::tclvalue(result_var)

  if (result == "cancel") {
    if (verbose) message("Selection cancelled by user")
    return(NULL)
  }

  # Process selections
  if (verbose) message("Processing selections...")

  # Get selected times
  selected$times <- times_section$items[sapply(times_section$check_vars, function(var) tcltk::tclvalue(var) == "1")]

  # Get selected doses
  selected$doses <- doses_section$items[sapply(doses_section$check_vars, function(var) tcltk::tclvalue(var) == "1")]

  # Get selected cell lines
  selected$cells <- cells_section$items[sapply(cells_section$check_vars, function(var) tcltk::tclvalue(var) == "1")]

  # Create configuration directory if it doesn't exist
  if (!dir.exists(config_dir)) {
    dir.create(config_dir, recursive = TRUE)
    if (verbose) message("Created configuration directory: ", config_dir)
  }

  # Create options list for all available parameters
  options_list <- list(
    times = params$times,
    doses = params$doses,
    cells = params$cells
  )

  # Create a function to find indices of selected items
  find_indices <- function(all_items, selected_items) {
    sapply(selected_items, function(item) which(all_items == item))
  }

  # Get indices of selected items
  selected_indices <- list(
    times = find_indices(params$times, selected$times),
    doses = find_indices(params$doses, selected$doses),
    cells = find_indices(params$cells, selected$cells)
  )

  # Write to configuration file
  config_file <- file.path(config_dir, config_filename)

  if (verbose) message("Writing configuration to ", config_file)

  # Write header and options
  lines <- character()
  lines <- c(lines, "# CMap experimental parameters configuration")
  lines <- c(lines, "# Generated on", as.character(Sys.time()))
  lines <- c(lines, "# This file was created using the interactive configuration tool")
  lines <- c(lines, "")

  # Write each section with available options
  for (section_name in names(options_list)) {
    lines <- c(lines, paste0("[", section_name, "]"))
    items <- options_list[[section_name]]
    for (i in seq_along(items)) {
      lines <- c(lines, paste0(i, " = ", items[i]))
    }
    lines <- c(lines, "")
  }

  # Write selected section
  lines <- c(lines, "[selected]")
  for (section_name in names(selected_indices)) {
    indices <- selected_indices[[section_name]]
    if (length(indices) > 0) {
      # Format indices as comma-separated list
      selected_str <- paste(indices, collapse = ",")
      lines <- c(lines, paste0(section_name, " = ", selected_str))
    } else {
      lines <- c(lines, paste0(section_name, " = "))
    }
  }

  # Write to file
  writeLines(lines, config_file)

  if (verbose) {
    message("Configuration saved successfully.")
    message("Selected ", length(selected$times), " time points, ",
            length(selected$doses), " dosages, and ",
            length(selected$cells), " cell lines.")
  }

  # Return the result with both selected parameters and config file path
  result <- list(
    selected = selected,
    all_options = options_list,
    config_file = config_file
  )

  return(result)
}

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
#' \dontrun{
#' # Create interactive console configuration
#' selections <- create_console_config("databases/siginfo_beta.txt")
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
      cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
      cat(sprintf("Select %s (n = %d)\n", section_name, length(items)))
      cat(paste(rep("-", 72), collapse = ""), "\n", sep = "")

      for (i in seq_along(items)) {
        cat(sprintf("%4d. %s\n", i, items[i]))
      }

      cat("\nInput options:\n")
      cat("  - Comma-separated indices, e.g. 1,3,5\n")
      cat("  - 'all' to select all\n")
      cat("  - 'none' to keep this section empty\n")

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
          cat("Invalid selection. Please try again.\n")
          next
        }

        selected_items <- items[idx]
        cat("Selected: ", paste(selected_items, collapse = ", "), "\n", sep = "")

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
#' @param selections The list returned by create_interactive_config() or create_console_config()
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return Data frame of combinations
#'
#' @examples
#' \dontrun{
#' # GUI version (if available)
#' selections <- create_interactive_config("databases/siginfo_beta.txt")
#' # OR console version
#' selections <- create_console_config("databases/siginfo_beta.txt")
#'
#' # Generate combinations
#' combinations <- generate_combinations_from_selections(selections)
#' }
#'
#' @export
generate_combinations_from_selections <- function(selections, verbose = TRUE) {
  # Check input
  if (is.null(selections)) {
    stop("No selections provided. Run create_interactive_config() or create_console_config() first.")
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
#' Tries to use the GUI interface if available, falls back to console version if not.
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
#' \dontrun{
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

  # Try GUI version first, fall back to console if not available
  selections <- tryCatch({
    if (requireNamespace("tcltk", quietly = TRUE)) {
      create_interactive_config(
        siginfo_file = siginfo_file,
        config_dir = config_dir,
        verbose = verbose
      )
    } else {
      message("GUI interface not available. Using console interface instead.")
      create_console_config(
        siginfo_file = siginfo_file,
        config_dir = config_dir,
        verbose = verbose
      )
    }
  }, error = function(e) {
    message("GUI interface error: ", e$message)
    message("Falling back to console interface.")
    create_console_config(
      siginfo_file = siginfo_file,
      config_dir = config_dir,
      verbose = verbose
    )
  })

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
