pkgname <- "CONCERTDR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "CONCERTDR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('CONCERTDR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("annotate_drug_results")
### * annotate_drug_results

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: annotate_drug_results
### Title: Annotate CONCERTDR Results with Drug Information
### Aliases: annotate_drug_results

### ** Examples

ex_results <- data.frame(
  compound = "CVD001_HEPG2_6H:BRD-K03652504-001-01-9:10.0497",
  Score = -0.72,
  pValue = 0.002,
  pAdjValue = 0.02,
  stringsAsFactors = FALSE
)
ex_siginfo <- data.frame(
  sig_id = ex_results$compound,
  pert_type = "trt_cp",
  pert_id = "BRD-K03652504-001-01-9",
  pert_iname = "imatinib",
  stringsAsFactors = FALSE
)
ex_compinfo <- data.frame(
  pert_id = "BRD-K03652504-001-01-9",
  pert_name = "imatinib",
  cmap_name = "imatinib",
  stringsAsFactors = FALSE
)
views <- annotate_drug_results(
  results_df = ex_results,
  sig_info_file = ex_siginfo,
  comp_info_file = ex_compinfo,
  write_outputs = FALSE,
  verbose = FALSE
)
names(views)

## No test: 
if (file.exists("sig_match_xsum_results.csv") &&
    file.exists("siginfo_beta.txt") &&
    file.exists("compoundinfo_beta.txt")) {
  results <- read.csv("sig_match_xsum_results.csv")

  views <- annotate_drug_results(
    results_df = results,
    sig_info_file = "siginfo_beta.txt",
    comp_info_file = "compoundinfo_beta.txt",
    output_dir = "results",
    write_outputs = TRUE
  )

  head(views$wetlab_drug_view)
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("annotate_drug_results", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("extract_cmap_data_from_siginfo")
### * extract_cmap_data_from_siginfo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: extract_cmap_data_from_siginfo
### Title: Extract CMap Data Using Siginfo File Directly
### Aliases: extract_cmap_data_from_siginfo

### ** Examples

is.function(extract_cmap_data_from_siginfo)

## Not run: 
##D # Example 1: Use the original siginfo_beta.txt directly
##D reference_df <- extract_cmap_data_from_siginfo(
##D   siginfo_file = "path/to/siginfo_beta.txt",
##D   geneinfo_file = "path/to/geneinfo_beta.txt",
##D   gctx_file = "path/to/level5_beta_trt_cp_n720216x12328.gctx"
##D )
##D 
##D # Example 2: Use a pre-filtered siginfo file
##D # First filter the siginfo
##D filtered_sig <- subset_siginfo_beta(
##D   "siginfo_beta.txt",
##D   interactive = FALSE,
##D   filters = list(
##D     pert_type = "trt_cp",
##D     pert_itime = c("6 h", "24 h"),
##D     pert_idose = "10 uM",
##D     cell_iname = c("A375", "MCF7")
##D   )
##D )
##D write.table(filtered_sig, "filtered_siginfo.txt", sep="\t", row.names=FALSE, quote=FALSE)
##D 
##D # Then use the filtered siginfo
##D reference_df <- extract_cmap_data_from_siginfo(
##D   siginfo_file = "filtered_siginfo.txt",
##D   geneinfo_file = "path/to/geneinfo_beta.txt",
##D   gctx_file = "path/to/level5_beta_trt_cp_n720216x12328.gctx",
##D   filter_quality = FALSE  # Already filtered
##D )
##D 
##D # Example 3: Process only first 5000 signatures for testing
##D reference_df <- extract_cmap_data_from_siginfo(
##D   siginfo_file = "path/to/siginfo_beta.txt",
##D   geneinfo_file = "path/to/geneinfo_beta.txt",
##D   gctx_file = "path/to/level5_beta_trt_cp_n720216x12328.gctx",
##D   max_signatures = 5000
##D )
##D 
##D # Access the metadata
##D metadata <- attr(reference_df, "metadata")
##D table(metadata$cell)
##D table(metadata$time, metadata$dose)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("extract_cmap_data_from_siginfo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("extract_compound_id")
### * extract_compound_id

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: extract_compound_id
### Title: Extract Compound ID from Signature Matching Results
### Aliases: extract_compound_id

### ** Examples

extract_compound_id("A:BRD-K03652504-001-01-9:10")

## No test: 
# Example compound strings
compounds <- c("CVD001_HEPG2_6H:BRD-K03652504-001-01-9:10.0497",
               "CVD001_HEPG2_6H:BRD-A37828317-001-03-0:10")

# Extract using colon splitting (default)
ids <- extract_compound_id(compounds)

# Extract using custom regex
ids <- extract_compound_id(compounds, method = "regex",
                          regex_pattern = "BRD-[A-Z0-9-]+")
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("extract_compound_id", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("extract_signature_zscores")
### * extract_signature_zscores

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: extract_signature_zscores
### Title: Extract z-score matrix for barcode heatmap
### Aliases: extract_signature_zscores

### ** Examples

is.function(extract_signature_zscores)

## No test: 
# Requires the full CMap GCTX file (downloaded from clue.io)
sig_file <- system.file("extdata", "example_signature.txt",
                        package = "CONCERTDR")
ref_file <- system.file("extdata", "example_reference_df.csv",
                        package = "CONCERTDR")
ref_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
ref_df$gene_symbol <- rownames(ref_df)
zmat <- extract_signature_zscores(
  results_df     = data.frame(sig_id = "DEMO001", Score = -0.72),
  signature_file = sig_file,
  reference_df   = ref_df,
  pert_id_col    = "sig_id"
)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("extract_signature_zscores", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fuzzy_drug_match")
### * fuzzy_drug_match

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fuzzy_drug_match
### Title: Fuzzy String Matching for Drug Names
### Aliases: fuzzy_drug_match

### ** Examples

if (requireNamespace("RecordLinkage", quietly = TRUE)) {
  fuzzy_drug_match(c("asprin"), c("aspirin", "ibuprofen"), threshold = 70)
}

## No test: 
if (requireNamespace("RecordLinkage", quietly = TRUE)) {
  query <- c("aspirin", "ibuprofen", "acetaminophen")
  reference <- c("aspirin", "ibuprofen", "acetaminophen", "naproxen", "diclofenac")
  matches <- fuzzy_drug_match(query, reference, threshold = 80)
  print(matches)
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fuzzy_drug_match", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.cmap_signature_result")
### * plot.cmap_signature_result

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.cmap_signature_result
### Title: Plot method for cmap_signature_result objects
### Aliases: plot.cmap_signature_result

### ** Examples

mock_res <- structure(
  list(
    results = list(
      ks = data.frame(compound = c("imatinib", "dasatinib", "doxorubicin"),
                      Score = c(-0.72, -0.45, 0.31),
                      pValue = c(0.002, 0.041, 0.280),
                      pAdjValue = c(0.020, 0.205, 0.560),
                      stringsAsFactors = FALSE)
    ),
    summary = data.frame(compound = "imatinib", method = "ks",
                         Score = -0.72, pValue = 0.002, global_rank = 1L,
                         stringsAsFactors = FALSE),
    gene_data    = data.frame(Gene = c("TP53", "MYC"),
                              log2FC = c(1.5, -2.1), stringsAsFactors = FALSE),
    settings     = list(signature_file = "ex.txt",
                        time_completed = Sys.time(), time_taken_mins = 0.1,
                        methods = "ks", permutations = 10),
    common_genes = list(up   = list(found = "TP53", count = 1L, percent = 50),
                        down = list(found = "MYC",  count = 1L, percent = 50))
  ), class = "cmap_signature_result"
)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(mock_res, method = "ks", plot_type = "scores")
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.cmap_signature_result", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_signature_direction_tile_barcode")
### * plot_signature_direction_tile_barcode

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_signature_direction_tile_barcode
### Title: Plot 2D Barcode Heatmap from GCTX (Single Drug Context)
### Aliases: plot_signature_direction_tile_barcode

### ** Examples

is.function(plot_signature_direction_tile_barcode)

## No test: 
if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {
  sig_file <- system.file("extdata", "example_signature.txt",
                          package = "CONCERTDR")
  plot_signature_direction_tile_barcode(
    results_df    = data.frame(sig_id = "DEMO001", Score = -0.72),
    signature_file = sig_file,
    reference_df  = read.csv(system.file("extdata", "example_reference_df.csv",
                                         package = "CONCERTDR"),
                             row.names = 1, check.names = FALSE),
    pert_id_col   = "sig_id"
  )
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_signature_direction_tile_barcode", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.cmap_signature_result")
### * print.cmap_signature_result

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.cmap_signature_result
### Title: Print method for cmap_signature_result objects
### Aliases: print.cmap_signature_result

### ** Examples

mock_res <- structure(
  list(
    results = list(
      ks = data.frame(compound = c("imatinib", "dasatinib"),
                      Score = c(-0.72, -0.45), pValue = c(0.002, 0.041),
                      pAdjValue = c(0.020, 0.205), stringsAsFactors = FALSE)
    ),
    summary = data.frame(compound = "imatinib", method = "ks",
                         Score = -0.72, pValue = 0.002, global_rank = 1L,
                         stringsAsFactors = FALSE),
    gene_data    = data.frame(Gene = c("TP53", "MYC"),
                              log2FC = c(1.5, -2.1), stringsAsFactors = FALSE),
    settings     = list(signature_file = "ex.txt",
                        time_completed = Sys.time(), time_taken_mins = 0.1,
                        methods = "ks", permutations = 10),
    common_genes = list(up   = list(found = "TP53", count = 1L, percent = 50),
                        down = list(found = "MYC",  count = 1L, percent = 50))
  ), class = "cmap_signature_result"
)
print(mock_res)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.cmap_signature_result", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("process_signature_with_df")
### * process_signature_with_df

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: process_signature_with_df
### Title: Process drug response signatures against reference data in a
###   dataframe
### Aliases: process_signature_with_df

### ** Examples

sig_file <- system.file("extdata", "example_signature.txt",
                        package = "CONCERTDR")
ref_file <- system.file("extdata", "example_reference_df.csv",
                        package = "CONCERTDR")
ref_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
ref_df$gene_symbol <- rownames(ref_df)
results <- process_signature_with_df(
  signature_file = sig_file,
  reference_df = ref_df,
  output_dir = tempdir(),
  permutations = 10,
  methods = "ks",
  save_files = FALSE
)
summary(results, top_n = 3)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("process_signature_with_df", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("subset_siginfo_beta")
### * subset_siginfo_beta

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: subset_siginfo_beta
### Title: Subset siginfo_beta file interactively with smart preview
### Aliases: subset_siginfo_beta

### ** Examples

ex_sig <- data.frame(
  pert_type = c("trt_cp", "trt_cp", "trt_sh"),
  pert_itime = c("6 h", "24 h", "24 h"),
  pert_idose = c("10 uM", "10 uM", "5 uM"),
  cell_iname = c("A375", "MCF7", "A375"),
  stringsAsFactors = FALSE
)
ex_sig_file <- tempfile(fileext = ".txt")
write.table(ex_sig, ex_sig_file, sep = "\t", row.names = FALSE, quote = FALSE)
subset_siginfo_beta(
  ex_sig_file,
  interactive = FALSE,
  filters = list(pert_type = "trt_cp"),
  verbose = FALSE,
  show_preview = FALSE
)

## No test: 
ex_file <- system.file("extdata", "example_siginfo.txt", package = "CONCERTDR")

# Interactive mode (run only in an interactive R session)
if (interactive() && nzchar(ex_file)) {
  filtered_siginfo <- subset_siginfo_beta(ex_file)
}

# Non-interactive mode with pre-defined filters
if (nzchar(ex_file)) {
  filtered_siginfo <- subset_siginfo_beta(
    ex_file,
    interactive = FALSE,
    filters = list(
      pert_type = "trt_cp",
      pert_itime = c("6 h", "24 h"),
      pert_idose = "10 uM",
      cell_iname = c("A375", "MCF7")
    ),
    verbose = FALSE,
    show_preview = FALSE
  )
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("subset_siginfo_beta", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.cmap_signature_result")
### * summary.cmap_signature_result

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.cmap_signature_result
### Title: Summary method for cmap_signature_result objects
### Aliases: summary.cmap_signature_result

### ** Examples

mock_res <- structure(
  list(
    results = list(
      ks = data.frame(compound = c("imatinib", "dasatinib"),
                      Score = c(-0.72, -0.45), pValue = c(0.002, 0.041),
                      pAdjValue = c(0.020, 0.205), stringsAsFactors = FALSE)
    ),
    summary = data.frame(compound = c("imatinib", "dasatinib"),
                         method = c("ks", "ks"),
                         Score = c(-0.72, -0.45), pValue = c(0.002, 0.041),
                         global_rank = 1:2, stringsAsFactors = FALSE),
    gene_data    = data.frame(Gene = c("TP53", "MYC"),
                              log2FC = c(1.5, -2.1), stringsAsFactors = FALSE),
    settings     = list(signature_file = "ex.txt",
                        time_completed = Sys.time(), time_taken_mins = 0.1,
                        methods = "ks", permutations = 10),
    common_genes = list(up   = list(found = "TP53", count = 1L, percent = 50),
                        down = list(found = "MYC",  count = 1L, percent = 50))
  ), class = "cmap_signature_result"
)
summary(mock_res, top_n = 2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.cmap_signature_result", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
