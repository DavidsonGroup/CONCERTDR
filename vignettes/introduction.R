## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  message   = FALSE,
  warning   = FALSE
)


## ----install, eval=FALSE------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("CONCERTDR")
# 
# # Optional but recommended for visualisation
# BiocManager::install("ComplexHeatmap")
# install.packages("circlize")


## ----load---------------------------------------------------------------------
library(CONCERTDR)


## ----bundled-files------------------------------------------------------------
sig_file  <- system.file("extdata", "example_signature.txt",
                          package = "CONCERTDR")
ref_file  <- system.file("extdata", "example_reference_df.csv",
                          package = "CONCERTDR")
sig_info_file <- system.file("extdata", "example_siginfo.txt",
                              package = "CONCERTDR")
gene_info_file <- system.file("extdata", "example_geneinfo.txt",
                               package = "CONCERTDR")


## ----read-signature-----------------------------------------------------------
signature <- read.delim(sig_file)
head(signature, 10)


## ----signature-summary--------------------------------------------------------
cat("Up-regulated genes:  ", sum(signature$log2FC > 0), "\n")
cat("Down-regulated genes:", sum(signature$log2FC < 0), "\n")


## ----use-df-signature, eval=FALSE---------------------------------------------
# results <- process_signature_with_df(
#   signature_file = signature,    # data.frame, not a file path
#   reference_df   = reference_df,
#   methods        = "ks"
# )


## ----gene-list-signature------------------------------------------------------
sig_from_lists <- create_signature_from_gene_lists(
  up_genes   = c("TP53", "MYC", "BRCA1"),
  down_genes = c("EGFR", "VEGFA", "KRAS")
)
sig_from_lists


## ----read-reference-----------------------------------------------------------
reference_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
# The gene_symbol column is required by process_signature_with_df()
reference_df$gene_symbol <- rownames(reference_df)

cat("Reference dimensions:", nrow(reference_df), "genes ×",
    ncol(reference_df) - 1, "signatures\n")
reference_df[1:5, 1:6]


## ----run-scoring--------------------------------------------------------------
results <- process_signature_with_df(
  signature_file = sig_file,
  reference_df   = reference_df,
  methods        = c("ks", "xsum"),
  topN           = 4,          # genes used by XSum; increase for real data
  permutations   = 100,        # increase to ≥1000 for publication-quality p-values
  save_files     = FALSE
)


## ----print-results------------------------------------------------------------
print(results)


## ----summary-results----------------------------------------------------------
summary(results, top_n = 5)


## ----inspect-ks---------------------------------------------------------------
ks_df <- results$results$ks
ks_df[order(ks_df$Score), ]


## ----plot-scores, fig.width=7, fig.height=4-----------------------------------
plot(results, method = "ks", plot_type = "scores", top_n = 10)


## ----plot-volcano, fig.width=6, fig.height=4----------------------------------
plot(results, method = "ks", plot_type = "volcano")


## ----set-options, eval=FALSE--------------------------------------------------
# data_dir <- "/path/to/your/cmap_data"
# 
# options(
#   CONCERTDR.data_dir          = data_dir,
#   CONCERTDR.gctx_file         = file.path(data_dir, "level5_beta_all_n1201944x12328.gctx"),
#   CONCERTDR.siginfo_file      = file.path(data_dir, "siginfo_beta.txt"),
#   CONCERTDR.geneinfo_file     = file.path(data_dir, "geneinfo_beta.txt"),
#   CONCERTDR.compoundinfo_file = file.path(data_dir, "compoundinfo_beta.txt")
# )
# 
# # Verify paths before running
# stopifnot(file.exists(getOption("CONCERTDR.gctx_file")))
# stopifnot(file.exists(getOption("CONCERTDR.siginfo_file")))


## ----session-info-------------------------------------------------------------
sessionInfo()

