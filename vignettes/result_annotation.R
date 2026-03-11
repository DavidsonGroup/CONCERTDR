## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  message  = FALSE,
  warning  = FALSE
)


## ----load---------------------------------------------------------------------
library(CONCERTDR)


## ----run-scoring--------------------------------------------------------------
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


## ----load-metadata------------------------------------------------------------
sig_info_file <- system.file("extdata", "example_siginfo.txt",
                              package = "CONCERTDR")
siginfo       <- read.delim(sig_info_file, stringsAsFactors = FALSE)
head(siginfo)


## ----build-compinfo-----------------------------------------------------------
compinfo <- data.frame(
  pert_id           = c("BRD-K00001", "BRD-K00002", "BRD-K00003",
                         "BRD-K00004", "BRD-K00005"),
  pert_name         = c("imatinib",    "dasatinib",   "doxorubicin",
                         "methotrexate", "vincristine"),
  cmap_name         = c("imatinib",    "dasatinib",   "doxorubicin",
                         "methotrexate", "vincristine"),
  target            = c("BCR-ABL|KIT", "BCR-ABL|SRC|LCK",
                         "TOP2A",        "DHFR",   "TUBB"),
  moa               = c("kinase inhibitor", "kinase inhibitor",
                         "DNA synthesis inhibitor",
                         "antimetabolite", "microtubule inhibitor"),
  phase             = c("Launched", "Launched", "Launched",
                         "Launched", "Launched"),
  stringsAsFactors  = FALSE
)

# Patch the bundled siginfo to include matching pert_id values
siginfo$pert_id <- compinfo$pert_id[
  match(siginfo$cmap_name, compinfo$cmap_name)
]
head(siginfo[, c("sig_id", "cmap_name", "pert_id", "pert_type")])


## ----annotate-----------------------------------------------------------------
views <- annotate_drug_results(
  results_df     = results_df,
  sig_info_file  = siginfo,     # data.frame, not a file path
  comp_info_file = compinfo,
  write_outputs  = FALSE,
  verbose        = TRUE
)


## ----view-names---------------------------------------------------------------
names(views)


## ----tech-view----------------------------------------------------------------
tv <- views$tech_view_all
# Key columns for drug repurposing interpretation
tv[order(tv$Score),
   c("sig_id", "perturbation_name", "Score", "pAdjValue",
     "effect_direction", "cell_line", "time_h", "dose_uM",
     "moa", "target")]


## ----effect-direction---------------------------------------------------------
table(tv$effect_direction)


## ----pert-kind----------------------------------------------------------------
table(tv$pert_kind, tv$pert_type)


## ----wetlab-drug--------------------------------------------------------------
views$wetlab_drug_view[,
  c("perturbation_name", "Score", "effect_direction",
    "target", "moa")]


## ----wetlab-drug-dose---------------------------------------------------------
views_with_dose <- annotate_drug_results(
  results_df             = results_df,
  sig_info_file          = siginfo,
  comp_info_file         = compinfo,
  keep_dose_in_drug_view = TRUE,
  write_outputs          = FALSE,
  verbose                = FALSE
)
views_with_dose$wetlab_drug_view[, c("perturbation_name", "Score", "dose_uM")]


## ----wetlab-gene--------------------------------------------------------------
nrow(views$wetlab_gene_view)
# In a full analysis this would list top shRNA/CRISPR hits:
# views$wetlab_gene_view[, c("pert_name", "mode", "Score", "cell_line")]


## ----drug-context-------------------------------------------------------------
views$drug_context_summary[,
  c("perturbation_name", "best_score", "n_contexts",
    "n_cell_lines", "moa_status")]


## ----write-outputs, eval=FALSE------------------------------------------------
# views_saved <- annotate_drug_results(
#   results_df     = results_df,
#   sig_info_file  = "path/to/siginfo_beta.txt",
#   comp_info_file = "path/to/compoundinfo_beta.txt",
#   output_dir     = "results",
#   write_outputs  = TRUE
# )
# # Written files:
# # results/wetlab_drug_view.tsv
# # results/wetlab_gene_view.tsv
# # results/tech_view_all.tsv
# # results/drug_context_summary.tsv


## ----show-output-files, eval=FALSE--------------------------------------------
# views_saved$output_files


## ----extract-compound-id------------------------------------------------------
example_ids <- c(
  "CVD001_K562_6H:BRD-K52492843-001-01-8:10:6",
  "CVD001_HL60_24H:BRD-A00619745-001-03-9:1:24"
)

# Extract the BRD identifier (field 2)
extract_compound_id(example_ids, method = "split_colon", part_index = 2)

# Extract using a regular expression
extract_compound_id(example_ids, method = "regex",
                    regex_pattern = "BRD-[A-Z0-9-]+")

# Extract the cell line from the first field
extract_compound_id(
  vapply(strsplit(example_ids, ":"), `[[`, character(1), 1),
  method = "split_underscore",
  part_index = 2
)


## ----fuzzy-match--------------------------------------------------------------
if (requireNamespace("RecordLinkage", quietly = TRUE)) {
  query_names <- c("imatinb", "doxorubicn", "vincristin")  # deliberate typos
  ref_names   <- c("imatinib", "dasatinib", "doxorubicin",
                   "methotrexate", "vincristine")

  matches <- fuzzy_drug_match(
    query_names     = query_names,
    reference_names = ref_names,
    method          = "levenshtein",
    threshold       = 70
  )
  matches
} else {
  message("Install RecordLinkage for fuzzy name matching: ",
          "install.packages('RecordLinkage')")
}


## ----full-workflow, eval=FALSE------------------------------------------------
# # 1. Set paths once
# options(
#   CONCERTDR.gctx_file         = "/data/cmap/level5_beta_all_n1201944x12328.gctx",
#   CONCERTDR.siginfo_file      = "/data/cmap/siginfo_beta.txt",
#   CONCERTDR.geneinfo_file     = "/data/cmap/geneinfo_beta.txt",
#   CONCERTDR.compoundinfo_file = "/data/cmap/compoundinfo_beta.txt"
# )
# 
# # 2. Filter siginfo to disease-relevant contexts
# filtered_siginfo <- subset_siginfo_beta(
#   getOption("CONCERTDR.siginfo_file"),
#   interactive = FALSE,
#   filters = list(
#     pert_type  = "trt_cp",
#     pert_itime = c("6 h", "24 h"),
#     cell_iname = c("K562", "HL60", "THP1", "JURKAT")
#   )
# )
# 
# # 3. Prepare your disease signature
# # Option A: read from a file
# signature <- read.delim("my_disease_signature.txt")
# # Option B: create from gene lists
# # signature <- create_signature_from_gene_lists(
# #   up_genes   = c("TP53", "MYC"),
# #   down_genes = c("EGFR", "KRAS")
# # )
# 
# # 4. Build reference matrix (this step may take a long time)
# reference_df <- extract_cmap_data_from_siginfo(
#   siginfo_file  = filtered_siginfo,
#   geneinfo_file = getOption("CONCERTDR.geneinfo_file"),
#   gctx_file     = getOption("CONCERTDR.gctx_file"),
#   filter_quality = FALSE,
#   landmark       = TRUE
# )
# 
# # 5. Score
# results <- process_signature_with_df(
#   signature_file = signature,   # data.frame or file path
#   reference_df   = reference_df,
#   methods        = c("xsum", "ks"),
#   topN           = 400,
#   permutations   = 1000,
#   save_files     = FALSE
# )
# 
# # 6. Annotate
# views <- annotate_drug_results(
#   results_df     = results$results$xsum,
#   sig_info_file  = filtered_siginfo,
#   comp_info_file = getOption("CONCERTDR.compoundinfo_file"),
#   output_dir     = "results",
#   write_outputs  = TRUE
# )
# 
# head(views$wetlab_drug_view)


## ----session-info-------------------------------------------------------------
sessionInfo()

