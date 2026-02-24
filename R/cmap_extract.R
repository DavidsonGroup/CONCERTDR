#' Extract Data from CMap GCTX Files for Specified Combinations
#'
#' @description Functions to extract expression data from CMap GCTX files based on
#' specified combinations of time points, dosages, and cell lines.
#'
#' @name cmap_extract
#' @keywords internal
NULL

#' Get Gene IDs and Gene Names for Landmark Genes
#'
#' @param geneinfo_df Data frame containing gene information
#' @return List with rid (gene IDs) and genenames (gene symbols)
#' @keywords internal
get_rid <- function(geneinfo_df,landmark=TRUE) {
  gene_overlapping <- geneinfo_df
  if (landmark)
  {
    gene_overlapping <- gene_overlapping[gene_overlapping$feature_space == "landmark", ]
  }

  return(list(
    rid = as.character(gene_overlapping$gene_id),
    genenames = as.character(gene_overlapping$gene_symbol)
  ))
}

#' Process a Single Combination of Time, Dose, and Cell Line
#'
#' @param combination Data frame row containing time, dose, and cell line information
#' @param rid Vector of gene IDs to extract
#' @param genenames Vector of gene names corresponding to the gene IDs
#' @param sig_info Data frame containing signature information
#' @param gctx_file Path to the GCTX file
#' @param output_dir Directory to save output files
#' @return Path to the output file (invisibly)
#' @keywords internal
process_combination <- function(combination, rid, genenames, sig_info,
                                gctx_file = "level5_beta_trt_cp_n720216x12328.gctx",
                                output_dir = ".") {
  itime <- combination$itime
  idose <- combination$idose
  cell <- combination$cell

  message(sprintf("Processing: time=%s, dose=%s, cell=%s", itime, idose, cell))

  # Filter data
  filtered_df <- sig_info[
    sig_info$pert_itime == itime &
      sig_info$pert_idose == idose &
      sig_info$cell_iname == cell,
  ]

  cid <- filtered_df$sig_id

  # Create filename with underscores instead of spaces
  filename <- paste0(
    "filtered_",
    gsub(" ", "_", itime), "_",
    gsub(" ", "_", idose), "_",
    gsub(" ", "_", cell), ".csv"
  )

  # Full path to output file
  output_file <- file.path(output_dir, filename)

  if (length(cid) > 0) {
    # Parse the gctx file using CmapR
    pert_data <- cmapR::parse_gctx(
      fname = gctx_file,
      cid = cid,
      rid = rid
    )

    # Get the data matrix
    mat <- cmapR::mat(pert_data)

    message(sprintf("Data dimensions: %d x %d", nrow(mat), ncol(mat)))

    if (nrow(mat) > 0) {
      # Convert to data frame and set row names as gene names
      df <- as.data.frame(mat)
      rownames(df) <- genenames

      # Write to file
      write.table(df, file = output_file, sep = "\t", quote = FALSE)
      message(sprintf("Successfully wrote data to %s", output_file))
    } else {
      message(sprintf("No data found for this combination"))
      # Create an empty file with header
      writeLines(paste(
        "# No data found for this combination",
        paste("# Time:", itime),
        paste("# Dose:", idose),
        paste("# Cell:", cell),
        sep = "\n"
      ), output_file)
    }
  } else {
    message(sprintf("No matching signatures found for this combination"))
    # Create an empty file with header
    writeLines(paste(
      "# No matching signatures found for this combination",
      paste("# Time:", itime),
      paste("# Dose:", idose),
      paste("# Cell:", cell),
      sep = "\n"
    ), output_file)
  }

  return(invisible(output_file))
}

#' Process Multiple Combinations of Time, Dose, and Cell Line
#'
#' @param combinations_file Path to a file containing combinations to process
#' @param task_id Specific task ID to process (for SLURM array jobs)
#' @param geneinfo_file Path to the gene info file
#' @param siginfo_file Path to the signature info file
#' @param gctx_file Path to the GCTX file
#' @param output_dir Directory to save output files
#' @return List of processed files (invisibly)
#' @export
process_combinations_file <- function(combinations_file, task_id = NULL,
                                      geneinfo_file = "geneinfo_beta.txt",
                                      siginfo_file = "siginfo_beta.txt",
                                      gctx_file = "level5_beta_trt_cp_n720216x12328.gctx",
                                      output_dir = ".") {

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Check if combinations file exists
  if (!file.exists(combinations_file)) {
    stop("Combinations file not found: ", combinations_file)
  }

  # Read combinations file
  combinations <- utils::read.table(combinations_file, header = TRUE, stringsAsFactors = FALSE)
  message(sprintf("Loaded %d combinations from %s", nrow(combinations), combinations_file))

  # Read gene info file
  message("Reading gene info file...")
  geneinfo_df <- utils::read.table(geneinfo_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  result <- get_rid(geneinfo_df)
  rid <- result$rid
  genenames <- result$genenames

  # Read siginfo file
  message("Reading signature info file...")
  sig_info <- utils::read.table(siginfo_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  sig_info <- sig_info[sig_info$pert_type == "trt_cp", ]
  sig_info <- sig_info[sig_info$is_hiq == 1, ]

  # Determine execution mode
  # Check if we're running as part of a SLURM array job
  slurm_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

  # Determine which task_id to use
  if (!is.null(task_id)) {
    # Command line argument takes precedence
    task_id <- as.integer(task_id)
    message(sprintf("Using explicit task ID: %d", task_id))
  } else if (slurm_task_id != "") {
    # Next, use SLURM array task ID if available
    task_id <- as.integer(slurm_task_id)
    message(sprintf("Using SLURM array task ID: %d", task_id))
  } else {
    # If neither is available, process all combinations
    task_id <- NULL
    message("No task ID specified. Processing all combinations sequentially.")
  }

  # Process combination(s)
  output_files <- character()

  if (!is.null(task_id)) {
    # Process just one combination
    if (task_id < 1 || task_id > nrow(combinations)) {
      stop("Invalid task ID: ", task_id, ". Must be between 1 and ", nrow(combinations))
    }

    output_file <- process_combination(combinations[task_id, ], rid, genenames, sig_info, gctx_file, output_dir)
    output_files <- c(output_files, output_file)
  } else {
    # Process all combinations
    for (i in 1:nrow(combinations)) {
      message(sprintf("\n--- Combination %d of %d ---\n", i, nrow(combinations)))
      output_file <- process_combination(combinations[i, ], rid, genenames, sig_info, gctx_file, output_dir)
      output_files <- c(output_files, output_file)
    }
  }

  return(invisible(output_files))
}
