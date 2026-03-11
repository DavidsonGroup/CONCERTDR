## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  message  = FALSE,
  warning  = FALSE
)


## ----load---------------------------------------------------------------------
library(CONCERTDR)


## ----score--------------------------------------------------------------------
sig_file     <- system.file("extdata", "example_signature.txt",
                             package = "CONCERTDR")
ref_file     <- system.file("extdata", "example_reference_df.csv",
                             package = "CONCERTDR")
reference_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
reference_df$gene_symbol <- rownames(reference_df)

scoring_results <- process_signature_with_df(
  signature_file = sig_file,
  reference_df   = reference_df,
  methods        = "ks",
  permutations   = 100,
  save_files     = FALSE
)
results_df <- scoring_results$results$ks
head(results_df[order(results_df$Score), ])


## ----extract-zscores----------------------------------------------------------
z <- extract_signature_zscores(
  results_df     = results_df,
  signature_file = sig_file,
  reference_df   = reference_df,
  max_genes      = 20,    # use all 20 genes in the example
  max_perts      = 10,    # show all 10 example signatures
  verbose        = TRUE
)


## ----inspect-z----------------------------------------------------------------
names(z)


## ----inspect-dimensions-------------------------------------------------------
cat("Matrix dimensions:", nrow(z$z_plot), "perturbations ×",
    ncol(z$z_plot), "genes\n")

# Genes are ordered: down-regulated (left) → up-regulated (right)
cat("First 5 genes (down-regulated in disease):\n",
    head(z$ordered_genes, 5), "\n")
cat("Last 5 genes (up-regulated in disease):\n",
    tail(z$ordered_genes, 5), "\n")

# Row labels combine drug name, dose, time, and cell line
head(z$sig_labels, 4)


## ----logfc-distribution-------------------------------------------------------
# Confirm the ordering: negative logFC on left, positive on right
lfc_ordered <- z$logfc_map[z$ordered_genes]
cat("Range of logFC for left half (down genes):",
    round(range(lfc_ordered[seq_len(sum(lfc_ordered < 0))]), 2), "\n")
cat("Range of logFC for right half (up genes):",
    round(range(lfc_ordered[seq_len(sum(lfc_ordered > 0)) + sum(lfc_ordered < 0)]), 2), "\n")


## ----save-zscores, eval=FALSE-------------------------------------------------
# z <- extract_signature_zscores(
#   results_df     = results_df,
#   signature_file = sig_file,
#   reference_df   = reference_df,
#   output_zscores = "results/barcode_zscores.tsv"
# )


## ----heatmap, fig.width=10, fig.height=6--------------------------------------
if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {

  plot_signature_direction_tile_barcode(
    precomputed  = z,
    cluster_rows = FALSE,   # keep score-order rows for now
    cluster_cols = FALSE,   # keep signature gene order
    verbose      = FALSE
  )

} else {
  message("Install ComplexHeatmap and circlize for heatmap rendering:\n",
          "  BiocManager::install('ComplexHeatmap')\n",
          "  install.packages('circlize')")
}


## ----heatmap-clustered, fig.width=10, fig.height=6----------------------------
if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {

  plot_signature_direction_tile_barcode(
    precomputed         = z,
    cluster_rows        = TRUE,
    cluster_method      = "complete",
    show_row_dendrogram = TRUE,
    cluster_cols        = FALSE,
    verbose             = FALSE
  )
}


## ----heatmap-png, eval=FALSE--------------------------------------------------
# plot_signature_direction_tile_barcode(
#   precomputed  = z,
#   cluster_rows = TRUE,
#   save_png     = TRUE,
#   output_png   = "results/barcode_heatmap.png",
#   width        = 14,     # inches; NULL = auto-compute from data size
#   height       = 8,
#   dpi          = 150
# )


## ----single-drug, eval=FALSE--------------------------------------------------
# # With full CMap data after annotate_drug_results():
# z_imatinib <- extract_signature_zscores(
#   results_df        = views$tech_view_all,
#   signature_file    = sig_file,
#   selected_drug     = "imatinib",
#   selected_drug_col = "display_name",
#   max_perts         = 30
# )
# 
# plot_signature_direction_tile_barcode(
#   precomputed = z_imatinib,
#   cluster_rows = TRUE,
#   save_png     = FALSE
# )


## ----gctx-workflow, eval=FALSE------------------------------------------------
# # Set global paths once
# options(
#   CONCERTDR.gctx_file     = "/data/cmap/level5_beta_all_n1201944x12328.gctx",
#   CONCERTDR.geneinfo_file = "/data/cmap/geneinfo_beta.txt",
#   CONCERTDR.siginfo_file  = "/data/cmap/siginfo_beta.txt"
# )
# 
# # extract_signature_zscores picks up paths from options automatically
# z <- extract_signature_zscores(
#   results_df     = views$tech_view_all,
#   signature_file = "my_signature.txt",
#   max_genes      = 100,
#   max_perts      = 60,
#   output_zscores = "results/barcode_zscores.tsv"
# )
# 
# # Plot — reuse the precomputed object for all iterations
# plot_signature_direction_tile_barcode(
#   precomputed  = z,
#   cluster_rows = TRUE,
#   save_png     = TRUE,
#   output_png   = "results/barcode_heatmap.png"
# )


## ----session-info-------------------------------------------------------------
sessionInfo()

