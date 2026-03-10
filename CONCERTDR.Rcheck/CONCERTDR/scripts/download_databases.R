#!/usr/bin/env Rscript

# Script to download required CMap databases
# Can be run directly from command line: Rscript download_databases.R [output_dir]

download_cmap_data <- function(dataset = c("metadata", "perturbations", "signatures", "all"),
                               output_dir = "databases",
                               overwrite = FALSE,
                               timeout = 600,
                               verbose = TRUE) {
  dataset <- match.arg(dataset)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  file_catalog <- list(
    metadata = c(
      siginfo_beta = "https://s3.amazonaws.com/data.clue.io/builds/LINCS2020/metadata/siginfo_beta.txt",
      geneinfo_beta = "https://s3.amazonaws.com/data.clue.io/builds/LINCS2020/metadata/geneinfo_beta.txt",
      compoundinfo_beta = "https://s3.amazonaws.com/data.clue.io/builds/LINCS2020/metadata/compoundinfo_beta.txt"
    ),
    perturbations = c(
      level5_trt_cp = "https://s3.amazonaws.com/data.clue.io/builds/LINCS2020/level5/level5_beta_trt_cp_n720216x12328.gctx"
    ),
    signatures = c(
      repurposing_drugs = "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt"
    )
  )

  selected_sets <- if (dataset == "all") names(file_catalog) else dataset
  urls <- unlist(file_catalog[selected_sets], use.names = TRUE)

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = timeout)

  out_files <- character()
  for (u in urls) {
    dest <- file.path(output_dir, basename(u))

    if (file.exists(dest) && !overwrite) {
      if (verbose) message("Using existing file: ", dest)
      out_files <- c(out_files, dest)
      next
    }

    if (verbose) message("Downloading ", basename(u), " ...")
    ok <- tryCatch({
      utils::download.file(url = u, destfile = dest, mode = "wb", quiet = !verbose)
      TRUE
    }, error = function(e) {
      warning("Failed to download ", u, ": ", e$message)
      FALSE
    })

    if (ok && file.exists(dest)) {
      out_files <- c(out_files, dest)
    }
  }

  invisible(out_files)
}

if (sys.nframe() == 0L) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- if (length(args) > 0) args[1] else "databases"

  message("Starting download of CMap databases to directory: ", output_dir)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Download the different types of data
  message("Downloading signatures data...")
  sig_files <- download_cmap_data("signatures", output_dir)

  message("Downloading metadata...")
  meta_files <- download_cmap_data("metadata", output_dir)

  message("Downloading perturbation data...")
  pert_files <- download_cmap_data("perturbations", output_dir)

  message("All downloads completed successfully.")
  message("Downloaded files are available in: ", normalizePath(output_dir))

  # Print a summary of downloaded files
  all_files <- c(sig_files, meta_files, pert_files)
  file_info <- file.info(all_files)

  message("\nSummary of downloaded files:")
  for (i in seq_along(all_files)) {
    file_size <- format(file_info$size[i] / 1024^2, digits = 2)
    message(sprintf("%s (%s MB)", basename(all_files[i]), file_size))
  }
}
