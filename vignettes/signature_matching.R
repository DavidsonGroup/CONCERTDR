## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  message  = FALSE,
  warning  = FALSE
)


## ----load---------------------------------------------------------------------
library(CONCERTDR)


## ----filter-noninteractive, eval=FALSE----------------------------------------
# filtered_siginfo <- subset_siginfo_beta(
#   siginfo_file = getOption("CONCERTDR.siginfo_file"),
#   interactive  = FALSE,
#   filters = list(
#     pert_type  = c("trt_cp"),          # small-molecule compounds only
#     pert_itime = c("6 h", "24 h"),     # standard time points
#     cell_iname = c("HL60", "THP1",     # AML-relevant cell lines
#                    "K562", "JURKAT",
#                    "U937", "HAP1")
#   )
# )
# cat("Retained signatures:", nrow(filtered_siginfo), "\n")


## ----filter-interactive, eval=FALSE-------------------------------------------
# filtered_siginfo <- subset_siginfo_beta(
#   siginfo_file = getOption("CONCERTDR.siginfo_file"),
#   interactive  = TRUE,
#   show_preview = TRUE
# )


## ----filter-demo--------------------------------------------------------------
example_siginfo_file <- system.file("extdata", "example_siginfo.txt",
                                     package = "CONCERTDR")
# Filter to K562 only
filtered_demo <- subset_siginfo_beta(
  siginfo_file = example_siginfo_file,
  interactive  = FALSE,
  filters      = list(cell_iname = "K562"),
  verbose      = TRUE,
  show_preview = FALSE
)
filtered_demo[, c("sig_id", "cmap_name", "cell_iname", "pert_itime", "pert_idose")]


## ----load-signature-----------------------------------------------------------
sig_file  <- system.file("extdata", "example_signature.txt",
                          package = "CONCERTDR")
signature <- read.delim(sig_file)

# Separate up- and down-regulated genes
up_genes   <- signature$Gene[signature$log2FC > 0]
down_genes <- signature$Gene[signature$log2FC < 0]

cat("Up-regulated genes:  ", length(up_genes),   "\n")
cat("Down-regulated genes:", length(down_genes),  "\n")

# Preview the most-extreme genes
head(signature[order(-abs(signature$log2FC)), ], 8)


## ----sig-df-demo--------------------------------------------------------------
# 'signature' is already a data.frame from the previous chunk
head(signature)
# We will pass it directly in Step 4 below


## ----gene-list-signature------------------------------------------------------
sig_from_lists <- create_signature_from_gene_lists(
  up_genes   = c("BRCA1", "TP53", "MYC"),
  down_genes = c("RB1", "PTEN", "APC")
)
sig_from_lists


## ----build-reference, eval=FALSE----------------------------------------------
# # Pass the in-memory filtered_siginfo directly — no need to write it to disk
# reference_df <- extract_cmap_data_from_siginfo(
#   siginfo_file   = filtered_siginfo,    # data.frame from Step 1
#   geneinfo_file  = getOption("CONCERTDR.geneinfo_file"),
#   gctx_file      = getOption("CONCERTDR.gctx_file"),
#   filter_quality = FALSE,   # already filtered above
#   landmark       = TRUE     # restrict to 978 landmark genes
# )
# 
# cat("Reference matrix:", nrow(reference_df), "genes ×",
#     ncol(reference_df) - 1, "signatures\n")
# 
# # Siginfo metadata is stored as an attribute
# metadata <- attr(reference_df, "metadata")
# head(metadata)


## ----load-reference-----------------------------------------------------------
ref_file     <- system.file("extdata", "example_reference_df.csv",
                             package = "CONCERTDR")
reference_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
reference_df$gene_symbol <- rownames(reference_df)

cat("Reference:", nrow(reference_df), "genes ×",
    ncol(reference_df) - 1, "signatures\n")


## ----run-single-method--------------------------------------------------------
results_ks <- process_signature_with_df(
  signature_file = sig_file,
  reference_df   = reference_df,
  methods        = "ks",
  permutations   = 100,
  save_files     = FALSE
)
print(results_ks)


## ----top-hits-ks--------------------------------------------------------------
ks_df <- results_ks$results$ks
ks_df[order(ks_df$Score), ]


## ----run-all-methods----------------------------------------------------------
results_all <- process_signature_with_df(
  signature_file = sig_file,
  reference_df   = reference_df,
  methods        = c("ks", "xsum", "xcos", "zhang", "gsea0", "gsea1", "gsea2"),
  topN           = 4,          # for real data use topN = 400
  permutations   = 100,
  save_files     = FALSE
)


## ----compare-methods----------------------------------------------------------
# Extract the top-3 ranked compound per method
lapply(results_all$results, function(df) {
  head(df[order(df$Score), c("compound", "Score", "pAdjValue")], 3)
})


## ----summary-all-methods------------------------------------------------------
summary(results_all, top_n = 5)


## ----score-interpretation-----------------------------------------------------
ks_df <- results_all$results$ks
ks_df$direction <- ifelse(
  ks_df$Score < 0, "Reversal (potentially therapeutic)",
  ifelse(ks_df$Score > 0, "Mimic / Aggravating", "Neutral")
)
ks_df[, c("compound", "Score", "direction")]


## ----plot-scores, fig.width=7, fig.height=4-----------------------------------
plot(results_all, method = "ks", plot_type = "scores", top_n = 10)


## ----plot-volcano, fig.width=6, fig.height=4----------------------------------
plot(results_all, method = "xsum", plot_type = "volcano")


## ----save-results, eval=FALSE-------------------------------------------------
# results_saved <- process_signature_with_df(
#   signature_file = sig_file,
#   reference_df   = reference_df,
#   methods        = c("ks", "xsum"),
#   permutations   = 1000,
#   output_dir     = "results",
#   save_files     = TRUE     # writes sig_match_ks_results.csv, etc.
# )


## ----topn-demo----------------------------------------------------------------
# Compare xsum at different topN values
scores_topn <- lapply(c(2, 4, 6), function(n) {
  res <- process_signature_with_df(
    signature_file = sig_file,
    reference_df   = reference_df,
    methods        = "xsum",
    topN           = n,
    permutations   = 50,
    save_files     = FALSE
  )
  df <- res$results$xsum[order(res$results$xsum$Score), ]
  df$topN <- n
  head(df, 3)
})
do.call(rbind, scores_topn)


## ----session-info-------------------------------------------------------------
sessionInfo()

