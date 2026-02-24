#!/usr/bin/env Rscript

# Script to download required CMap databases
# Can be run directly from command line: Rscript download_databases.R [output_dir]

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) > 0) args[1] else "databases"

# Load download helper from sibling script
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- all_args[grepl("^--file=", all_args)]
if (length(file_arg) == 0) {
  stop("Unable to determine script path from command line")
}
script_path <- sub("^--file=", "", file_arg[1])
script_dir <- dirname(normalizePath(script_path))
source(file.path(script_dir, "download_databases.R"))

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
