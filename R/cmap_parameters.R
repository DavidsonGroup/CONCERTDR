#' Extract CMap experimental parameters from siginfo file
#'
#' Analyzes the CMap siginfo_beta.txt file to extract all available
#' experimental parameters (time points, dosages, cell lines) and saves them to
#' a configuration file for use in the DR pipeline.
#'
#' @param siginfo_file Path to siginfo_beta.txt file
#' @param write_config Logical; whether to write the configuration file (default: TRUE)
#' @param config_dir Directory to save the configuration file; if NULL uses "conf"
#' @param config_filename Name of the configuration file to write (default: "cmap_options.conf")
#' @param filter_quality Logical; whether to filter for high-quality compounds only (default: TRUE)
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return A list with extracted parameters (times, doses, cells) and the config file path
#'
#' @examples
#' \dontrun{
#' # Extract parameters and write config
#' params <- extract_cmap_parameters("databases/siginfo_beta.txt")
#' }
#'
#' @export
extract_cmap_parameters <- function(siginfo_file, write_config = TRUE,
                                    config_dir = "conf",
                                    config_filename = "cmap_options.conf",
                                    filter_quality = TRUE,
                                    verbose = TRUE) {

  # Validate inputs
  if (!file.exists(siginfo_file)) {
    stop("Siginfo file not found: ", siginfo_file)
  }

  # Set up configuration directory
  if (write_config) {
    # Create config directory if it doesn't exist
    if (!dir.exists(config_dir)) {
      dir.create(config_dir, recursive = TRUE)
      if (verbose) message("Created configuration directory: ", config_dir)
    }

    output_file <- file.path(config_dir, config_filename)
  }

  # Read siginfo file
  if (verbose) message("Reading signature info from ", siginfo_file)
  sig_info <- try(utils::read.table(siginfo_file, sep = "\t", header = TRUE,
                                    quote = "", stringsAsFactors = FALSE))

  if (inherits(sig_info, "try-error")) {
    stop("Failed to read siginfo file. Make sure it's properly formatted.")
  }

  # Filter for treatment compounds with high quality if requested
  if (filter_quality) {
    if (!"pert_type" %in% names(sig_info) || !"is_hiq" %in% names(sig_info)) {
      warning("Could not find pert_type or is_hiq columns for quality filtering")
    } else {
      sig_info <- sig_info[sig_info$pert_type == "trt_cp", ]
      sig_info <- sig_info[sig_info$is_hiq == 1, ]
    }
  }

  # Extract unique values
  if (verbose) message("Extracting unique experimental parameters...")

  # Check if required columns exist
  required_cols <- c("pert_itime", "pert_idose", "cell_iname")
  missing_cols <- setdiff(required_cols, names(sig_info))

  if (length(missing_cols) > 0) {
    stop("Missing required columns in siginfo file: ",
         paste(missing_cols, collapse = ", "))
  }

  # Extract and sort values
  times <- sort(unique(sig_info$pert_itime))
  doses <- sort(unique(sig_info$pert_idose))
  cells <- sort(unique(sig_info$cell_iname))

  # Create options list
  options_list <- list(
    times = times,
    doses = doses,
    cells = cells
  )

  # Write config if requested
  if (write_config) {
    if (verbose) message("Writing options to ", output_file)
    write_config_file(options_list, output_file)

    if (verbose) {
      message("Done! Configuration file created at: ", output_file)
      message("Available options:")
      message("  Times: ", length(times), " options")
      message("  Doses: ", length(doses), " options")
      message("  Cell lines: ", length(cells), " options")
      message("\nEdit the [selected] section in the config file to choose which options to use.")
      message("To select all options for a parameter, use 'all' or leave the value empty.")
      message("Then run generate_combinations_from_config() to create data extracts.")
    }
  }

  # Return results
  result <- list(
    times = times,
    doses = doses,
    cells = cells,
    config_file = if (write_config) output_file else NULL
  )

  return(invisible(result))
}

#' Process combinations by extracting data for each
#'
#' @param combinations Data frame of combinations to process
#' @param output_dir Directory for output files (default: "output")
#' @param geneinfo_file Path to geneinfo file (default: "geneinfo_beta.txt")
#' @param siginfo_file Path to siginfo file (default: "siginfo_beta.txt")
#' @param gctx_file Path to gctx file (default: "level5_beta_trt_cp_n720216x12328.gctx")
#' @param install_packages Logical; whether to install required packages if missing (default: FALSE)
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return Invisibly returns vector of output files
#'
#' @examples
#' \dontrun{
#' # First create a configuration file
#' extract_cmap_parameters("databases/siginfo_beta.txt")
#' # Edit the configuration file as needed
#'
#' # Then generate combinations and process them
#' combinations <- generate_combinations_from_config("conf/cmap_options.conf")
#' process_combinations(combinations,
#'                      output_dir = "output",
#'                      geneinfo_file = "databases/geneinfo_beta.txt",
#'                      siginfo_file = "databases/siginfo_beta.txt",
#'                      gctx_file = "databases/level5_beta_trt_cp_n720216x12328.gctx")
#' }
#'
#' @export
process_combinations <- function(combinations, output_dir = "output",
                                 geneinfo_file = "geneinfo_beta.txt",
                                 siginfo_file = "siginfo_beta.txt",
                                 gctx_file = "level5_beta_trt_cp_n720216x12328.gctx",
                                 install_packages = FALSE,
                                 verbose = TRUE) {
  # Check required packages
  required_packages <- c("BiocManager", "data.table", "rhdf5")
  missing_packages <- setdiff(required_packages, utils::installed.packages()[,"Package"])

  if (length(missing_packages) > 0) {
    if (install_packages && verbose) {
      warning("install_packages=TRUE is deprecated and automatic installation is disabled in package functions.")
    }

    cran_missing <- intersect(missing_packages, c("BiocManager", "data.table"))
    bioc_missing <- intersect(missing_packages, c("rhdf5"))

    install_msg <- character()
    if (length(cran_missing) > 0) {
      install_msg <- c(install_msg,
                       paste0("install.packages(c('", paste(cran_missing, collapse = "', '"), "'))"))
    }
    if (length(bioc_missing) > 0) {
      install_msg <- c(install_msg,
                       paste0("BiocManager::install(c('", paste(bioc_missing, collapse = "', '"), "'))")
      )
    }

    stop(
      "Missing required packages: ", paste(missing_packages, collapse = ", "), ". ",
      "Please install them before running this function. Suggested commands: ",
      paste(install_msg, collapse = " ; ")
    )
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Created output directory: ", output_dir)
  }

  # Read gene info file and get rid & gene names
  if (verbose) message("Reading gene info file: ", geneinfo_file)
  tryCatch({
    if (requireNamespace("data.table", quietly = TRUE)) {
      geneinfo_df <- data.table::fread(geneinfo_file, header = TRUE)
    } else {
      geneinfo_df <- utils::read.table(geneinfo_file, sep = "\t", header = TRUE)
    }
    result <- get_rid(geneinfo_df)
    rid <- result$rid
    genenames <- result$genenames

    if (verbose) message("Found ", length(rid), " landmark genes")
  }, error = function(e) {
    stop("Error reading gene info file: ", e$message)
  })

  # Read siginfo file
  if (verbose) message("Reading signature info file: ", siginfo_file)
  tryCatch({
    if (requireNamespace("data.table", quietly = TRUE)) {
      sig_info <- data.table::fread(siginfo_file, header = TRUE, stringsAsFactors = FALSE)
    } else {
      sig_info <- utils::read.table(siginfo_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    }
    #sig_info <- sig_info[sig_info$pert_type == "trt_cp", ]
    sig_info <- sig_info[sig_info$is_hiq == 1, ]

    if (verbose) message("Found ", nrow(sig_info), " high-quality treatment signatures")
  }, error = function(e) {
    stop("Error reading signature info file: ", e$message)
  })

  # Process each combination
  output_files <- character(nrow(combinations))

  for (i in 1:nrow(combinations)) {
    if (verbose) {
      message(sprintf("\nProcessing combination %d of %d:", i, nrow(combinations)))
    }

    # Process using the core function
    output_files[i] <- process_combination(
      combinations[i, ], rid, genenames, sig_info, gctx_file, output_dir
    )
  }

  if (verbose) {
    message("\nProcessing complete. Generated ", length(output_files), " output files in ", output_dir)
  }

  return(invisible(output_files))
}
