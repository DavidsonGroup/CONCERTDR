#!/usr/bin/env Rscript

# Recreate the five files distributed in inst/extdata/.
#
# Four files are subsets of the Broad Institute Expanded CMap LINCS2020
# release:
#
#   example_siginfo.txt       <- siginfo_beta.txt
#   example_geneinfo.txt      <- geneinfo_beta.txt
#   example_compoundinfo.txt  <- compoundinfo_beta.txt
#   example_reference_df.csv  <- three Level 5 GCTX matrices
#
# The source files were downloaded from these public release locations:
#
# https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/
#   siginfo_beta.txt
#   geneinfo_beta.txt
#   compoundinfo_beta.txt
#   level5/level5_beta_trt_cp_n720216x12328.gctx
#   level5/level5_beta_trt_oe_n34171x12328.gctx
#   level5/level5_beta_trt_xpr_n142901x12328.gctx
#
# The checked-in example_siginfo.txt and example_geneinfo.txt files are the
# selection manifests. They make the example stable if CMap metadata ordering
# changes. The signature manifest contains 496 high-quality, QC-passing
# profiles from four well-represented cell lines. It covers 30 AML-relevant
# compounds, 22 overexpression perturbagens, and 21 CRISPR perturbagens. The
# gene manifest contains all 58 genes in the demonstration query plus 105
# additional cancer and drug-response genes available in the GCTX matrices.
#
# example_signature.txt is deliberately copied rather than derived from CMap.
# It is a manually curated AML/BCR-ABL-themed query with deterministic example
# weights, not a patient-derived differential-expression result.
#
# Usage:
#   Rscript inst/scripts/make_example_data.R RAW_CMAP_DIR [OUTPUT_DIR]
#
# RAW_CMAP_DIR must contain the six source files listed above. OUTPUT_DIR
# defaults to example-data-regenerated/ in the current directory. The default
# does not overwrite the package's checked-in data.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop(
    "Usage: Rscript inst/scripts/make_example_data.R ",
    "RAW_CMAP_DIR [OUTPUT_DIR]",
    call. = FALSE
  )
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.", call. = FALSE)
}
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  stop("Bioconductor package 'rhdf5' is required.", call. = FALSE)
}

raw_dir <- normalizePath(args[1L], mustWork = TRUE)
output_dir <- if (length(args) == 2L) {
  args[2L]
} else {
  file.path(getwd(), "example-data-regenerated")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

# Find this script so it works both from a source checkout and after package
# installation. The checked-in outputs are used only as immutable manifests.
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- all_args[grepl("^--file=", all_args)]
if (length(file_arg) != 1L) {
  stop("Unable to determine the path to this script.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
manifest_dir <- normalizePath(
  file.path(dirname(script_path), "..", "extdata"),
  mustWork = TRUE
)

source_files <- c(
  siginfo = "siginfo_beta.txt",
  geneinfo = "geneinfo_beta.txt",
  compoundinfo = "compoundinfo_beta.txt",
  trt_cp = "level5_beta_trt_cp_n720216x12328.gctx",
  trt_oe = "level5_beta_trt_oe_n34171x12328.gctx",
  trt_xpr = "level5_beta_trt_xpr_n142901x12328.gctx"
)
source_paths <- setNames(
  file.path(raw_dir, unname(source_files)),
  names(source_files)
)
missing_files <- source_paths[!file.exists(source_paths)]
if (length(missing_files)) {
  stop(
    "Missing CMap source file(s):\n  ",
    paste(missing_files, collapse = "\n  "),
    call. = FALSE
  )
}

read_tsv <- function(path, select = NULL) {
  data.table::fread(
    path,
    sep = "\t",
    select = select,
    data.table = FALSE,
    showProgress = interactive()
  )
}

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

read_first_h5_dataset <- function(path, candidates) {
  for (candidate in candidates) {
    value <- tryCatch(
      rhdf5::h5read(path, candidate),
      error = function(e) NULL
    )
    if (!is.null(value)) {
      return(value)
    }
  }
  stop(
    "None of the expected datasets was found in ", basename(path), ": ",
    paste(candidates, collapse = ", "),
    call. = FALSE
  )
}

read_gctx_subset <- function(path, gene_ids, signature_ids) {
  row_ids <- trimws(as.character(read_first_h5_dataset(
    path,
    c("/0/META/ROW/id", "/META/ROW/id")
  )))
  column_ids <- trimws(as.character(read_first_h5_dataset(
    path,
    c("/0/META/COL/id", "/META/COL/id")
  )))

  row_index <- match(gene_ids, row_ids)
  column_index <- match(signature_ids, column_ids)
  if (anyNA(row_index)) {
    stop(
      "Genes missing from ", basename(path), ": ",
      paste(gene_ids[is.na(row_index)], collapse = ", "),
      call. = FALSE
    )
  }
  if (anyNA(column_index)) {
    stop(
      "Signatures missing from ", basename(path), ": ",
      paste(signature_ids[is.na(column_index)], collapse = ", "),
      call. = FALSE
    )
  }

  matrix_value <- tryCatch(
    rhdf5::h5read(
      path,
      "/0/DATA/0/matrix",
      index = list(row_index, column_index)
    ),
    error = function(e) {
      rhdf5::h5read(
        path,
        "/DATA/0/matrix",
        index = list(row_index, column_index)
      )
    }
  )
  matrix_value <- as.matrix(matrix_value)
  rownames(matrix_value) <- gene_ids
  colnames(matrix_value) <- signature_ids
  matrix_value
}

siginfo_columns <- c(
  "sig_id", "pert_id", "pert_mfc_id", "cmap_name", "pert_type",
  "cell_iname", "pert_idose", "pert_itime", "is_hiq", "cc_q75"
)
raw_siginfo_columns <- c(siginfo_columns, "qc_pass")
geneinfo_columns <- c(
  "gene_id", "gene_symbol", "gene_title", "gene_type", "src",
  "feature_space"
)
compoundinfo_columns <- c(
  "pert_id", "cmap_name", "target", "moa", "compound_aliases",
  "inchi_key"
)

sig_manifest <- read_tsv(
  file.path(manifest_dir, "example_siginfo.txt"),
  select = siginfo_columns
)
gene_manifest <- read_tsv(
  file.path(manifest_dir, "example_geneinfo.txt"),
  select = geneinfo_columns
)

if (anyDuplicated(sig_manifest$sig_id)) {
  stop("The signature selection manifest contains duplicate IDs.")
}
if (anyDuplicated(gene_manifest$gene_id)) {
  stop("The gene selection manifest contains duplicate IDs.")
}
if (!all(sig_manifest$is_hiq == 1L)) {
  stop("Every selected CMap signature must have is_hiq = 1.")
}

message("Reading and subsetting CMap metadata ...")
raw_siginfo <- read_tsv(source_paths[["siginfo"]], raw_siginfo_columns)
sig_index <- match(sig_manifest$sig_id, raw_siginfo$sig_id)
if (anyNA(sig_index)) {
  stop(
    "Selected signatures missing from siginfo_beta.txt: ",
    paste(sig_manifest$sig_id[is.na(sig_index)], collapse = ", ")
  )
}
if (!all(raw_siginfo$qc_pass[sig_index] == 1L)) {
  stop("Every selected CMap signature must have qc_pass = 1.")
}
example_siginfo <- raw_siginfo[sig_index, siginfo_columns, drop = FALSE]
if (!identical(example_siginfo$sig_id, sig_manifest$sig_id)) {
  stop("Signature order changed unexpectedly.")
}

raw_geneinfo <- read_tsv(source_paths[["geneinfo"]], geneinfo_columns)
gene_index <- match(gene_manifest$gene_id, raw_geneinfo$gene_id)
if (anyNA(gene_index)) {
  stop(
    "Selected genes missing from geneinfo_beta.txt: ",
    paste(gene_manifest$gene_id[is.na(gene_index)], collapse = ", ")
  )
}
example_geneinfo <- raw_geneinfo[gene_index, geneinfo_columns, drop = FALSE]

raw_compoundinfo <- read_tsv(
  source_paths[["compoundinfo"]],
  compoundinfo_columns
)
compound_ids <- unique(
  example_siginfo$pert_id[example_siginfo$pert_type == "trt_cp"]
)
example_compoundinfo <- raw_compoundinfo[
  raw_compoundinfo$pert_id %in% compound_ids,
  compoundinfo_columns,
  drop = FALSE
]
missing_compounds <- setdiff(compound_ids, example_compoundinfo$pert_id)
if (length(missing_compounds)) {
  stop(
    "Selected compounds missing from compoundinfo_beta.txt: ",
    paste(missing_compounds, collapse = ", ")
  )
}

message("Extracting Level 5 z-scores from the three GCTX matrices ...")
gene_ids <- as.character(example_geneinfo$gene_id)
matrices <- lapply(c("trt_cp", "trt_oe", "trt_xpr"), function(type) {
  signature_ids <- example_siginfo$sig_id[
    example_siginfo$pert_type == type
  ]
  read_gctx_subset(source_paths[[type]], gene_ids, signature_ids)
})
reference_matrix <- do.call(cbind, matrices)
reference_matrix <- reference_matrix[
  gene_ids,
  example_siginfo$sig_id,
  drop = FALSE
]
rownames(reference_matrix) <- example_geneinfo$gene_symbol
reference_matrix <- round(reference_matrix, digits = 2L)

write_tsv(
  example_siginfo,
  file.path(output_dir, "example_siginfo.txt")
)
write_tsv(
  example_geneinfo,
  file.path(output_dir, "example_geneinfo.txt")
)
write_tsv(
  example_compoundinfo,
  file.path(output_dir, "example_compoundinfo.txt")
)
utils::write.csv(
  as.data.frame(reference_matrix, check.names = FALSE),
  file = file.path(output_dir, "example_reference_df.csv"),
  row.names = TRUE
)
file.copy(
  file.path(manifest_dir, "example_signature.txt"),
  file.path(output_dir, "example_signature.txt"),
  overwrite = TRUE
)

created <- file.path(
  output_dir,
  c(
    "example_signature.txt", "example_reference_df.csv",
    "example_siginfo.txt", "example_geneinfo.txt",
    "example_compoundinfo.txt"
  )
)
message("Created example files in ", output_dir, ":")
print(data.frame(
  file = basename(created),
  bytes = file.info(created)$size,
  md5 = unname(tools::md5sum(created)),
  row.names = NULL
))
