# tests/testthat/test-signature-matching-extended.R
#
# Extended tests for process_signature_with_df(), create_summary_from_results(),
# and the cmap_signature_result S3 methods (print / summary / plot).
# Complements the existing test-signature-matching.R.

# ── Helpers ──────────────────────────────────────────────────────────────────

sig_file <- function() {
  system.file("extdata", "example_signature.txt", package = "CONCERTDR")
}

ref_df <- function() {
  f <- system.file("extdata", "example_reference_df.csv", package = "CONCERTDR")
  df <- read.csv(f, row.names = 1, check.names = FALSE)
  df$gene_symbol <- rownames(df)
  df
}

skip_if_no_example <- function() {
  skip_if(!nzchar(sig_file()), "Example signature file not found")
}

# ── Input validation ─────────────────────────────────────────────────────────

test_that("process_signature_with_df errors if signature_file does not exist", {
  expect_error(
    process_signature_with_df(
      signature_file = "no_such_file.txt",
      reference_df   = ref_df(),
      permutations   = 5,
      methods        = "ks",
      save_files     = FALSE
    ),
    regexp = NULL   # any error is acceptable
  )
})

test_that("process_signature_with_df errors if gene_symbol column is absent", {
  skip_if_no_example()
  bad_ref <- ref_df()
  bad_ref$gene_symbol <- NULL

  expect_error(
    process_signature_with_df(
      signature_file = sig_file(),
      reference_df   = bad_ref,
      permutations   = 5,
      methods        = "ks",
      save_files     = FALSE
    ),
    "gene_symbol"
  )
})

test_that("process_signature_with_df errors if no genes overlap", {
  skip_if_no_example()
  # Reference with totally different gene names
  bad_ref <- ref_df()
  rownames(bad_ref) <- paste0("FAKE_", rownames(bad_ref))
  bad_ref$gene_symbol <- rownames(bad_ref)

  expect_error(
    process_signature_with_df(
      signature_file = sig_file(),
      reference_df   = bad_ref,
      permutations   = 5,
      methods        = "ks",
      save_files     = FALSE
    ),
    "No matching genes"
  )
})

test_that("process_signature_with_df messages and skips invalid method names", {
  skip_if_no_example()
  # The function uses message(), not warning(), to report invalid methods
  expect_message(
    process_signature_with_df(
      signature_file = sig_file(),
      reference_df   = ref_df(),
      permutations   = 5,
      methods        = c("ks", "INVALID_METHOD"),
      save_files     = FALSE
    ),
    "Invalid methods"
  )
})

test_that("process_signature_with_df errors if only invalid methods are given", {
  skip_if_no_example()
  expect_error(
    suppressWarnings(
      process_signature_with_df(
        signature_file = sig_file(),
        reference_df   = ref_df(),
        permutations   = 5,
        methods        = "INVALID_METHOD",
        save_files     = FALSE
      )
    ),
    "No valid methods"
  )
})

# ── Return object structure ───────────────────────────────────────────────────

test_that("process_signature_with_df returns a cmap_signature_result", {
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )
  expect_s3_class(res, "cmap_signature_result")
  expect_named(res, c("results", "summary", "gene_data", "settings", "common_genes"))
})

test_that("result list has all expected components", {
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = c("ks", "xsum"),
    save_files     = FALSE
  )
  # results contains both methods
  expect_true(all(c("ks", "xsum") %in% names(res$results)))

  # each method result has required columns
  for (method in c("ks", "xsum")) {
    df <- res$results[[method]]
    expect_true(all(c("compound", "Score", "pValue", "pAdjValue", "rank") %in% names(df)))
    expect_equal(nrow(df), 10)  # example_reference_df has 10 signatures
  }
})

test_that("settings record correct method list and permutation count", {
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 7,
    methods        = c("ks", "zhang"),
    save_files     = FALSE
  )
  expect_equal(res$settings$permutations, 7)
  expect_setequal(res$settings$methods, c("ks", "zhang"))
})

test_that("common_genes counts are positive and leq total signature genes", {
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )
  expect_gt(res$common_genes$up$count, 0)
  expect_gt(res$common_genes$down$count, 0)
  expect_lte(res$common_genes$up$percent,   100)
  expect_lte(res$common_genes$down$percent, 100)
})

# ── All seven methods run without error ───────────────────────────────────────

test_that("all seven scoring methods complete successfully on example data", {
  skip_if_no_example()
  all_methods <- c("ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang")

  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 10,
    methods        = all_methods,
    topN           = 4,
    save_files     = FALSE
  )

  expect_setequal(names(res$results), all_methods)
  for (method in all_methods) {
    expect_gt(nrow(res$results[[method]]), 0,
              label = paste("method", method, "returned empty data frame"))
  }
})

# ── save_files = TRUE ─────────────────────────────────────────────────────────

test_that("save_files = TRUE writes CSVs to output_dir", {
  skip_if_no_example()
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    output_dir     = tmp,
    permutations   = 5,
    methods        = "ks",
    save_files     = TRUE
  )

  written <- list.files(tmp)
  expect_true(any(grepl("sig_match_ks", written)))
  expect_true(any(grepl("summary_results", written)))
})

test_that("save_files = FALSE writes nothing", {
  skip_if_no_example()
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    output_dir     = tmp,
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )

  expect_equal(length(list.files(tmp)), 0)
})

# ── topN parameter ────────────────────────────────────────────────────────────

test_that("topN = 2 still produces valid xsum results", {
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = c("xsum", "xcos"),
    topN           = 2,
    save_files     = FALSE
  )
  expect_gt(nrow(res$results$xsum), 0)
  expect_gt(nrow(res$results$xcos), 0)
})

# ── create_summary_from_results (internal) ───────────────────────────────────

test_that("create_summary_from_results returns empty df when results_list is empty", {
  out <- CONCERTDR:::create_summary_from_results(list())
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})

test_that("create_summary_from_results skips methods with errors", {
  bad_result <- data.frame(compound = character(0), Score = numeric(0),
                           error = character(0))
  good_result <- data.frame(
    compound = c("drugA","drugB"),
    Score    = c(0.8, 0.3),
    pValue   = c(0.01, 0.2),
    rank     = 1:2
  )
  out <- CONCERTDR:::create_summary_from_results(
    list(bad_method = bad_result, good_method = good_result),
    top_n = 2
  )
  expect_true(all(out$method == "good_method"))
})

test_that("create_summary_from_results includes global_rank column", {
  r1 <- data.frame(compound = c("A","B","C"),
                   Score = c(0.9, 0.5, 0.1),
                   pValue = c(0.01, 0.1, 0.5),
                   rank = 1:3)
  out <- CONCERTDR:::create_summary_from_results(list(ks = r1), top_n = 3)
  expect_true("global_rank" %in% names(out))
  expect_true(all(!is.na(out$global_rank)))
})

test_that("create_summary_from_results respects top_n", {
  r1 <- data.frame(
    compound = paste0("drug", 1:20),
    Score    = runif(20),
    pValue   = runif(20),
    rank     = 1:20
  )
  out <- CONCERTDR:::create_summary_from_results(list(ks = r1), top_n = 5)
  expect_lte(nrow(out), 5)
})

# ── S3 methods ────────────────────────────────────────────────────────────────

test_that("print.cmap_signature_result is invisible and prints to console", {
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )
  # capture.output() itself makes the return value visible, so use
  # withVisible() to check invisibility, and capture output separately.
  output <- capture.output(wv <- withVisible(print(res)))
  expect_false(wv$visible)
  expect_true(any(grepl("CMap Signature Matching", output)))
  expect_true(any(grepl("Methods used", output)))
  expect_true(any(grepl("Gene coverage", output)))
})

test_that("summary.cmap_signature_result returns top_n rows", {
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )
  out <- summary(res, top_n = 3)
  expect_lte(nrow(out), 3)
  expect_true("global_rank" %in% names(out))
})

test_that("summary returns NULL and messages when summary is empty", {
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )
  res$summary <- res$summary[0, ]
  expect_null(summary(res))
})

test_that("plot.cmap_signature_result returns ggplot for 'scores' type", {
  skip_if_not_installed("ggplot2")
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )
  p <- plot(res, method = "ks", plot_type = "scores", top_n = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot.cmap_signature_result returns ggplot for 'volcano' type", {
  skip_if_not_installed("ggplot2")
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )
  p <- plot(res, method = "ks", plot_type = "volcano")
  expect_s3_class(p, "ggplot")
})

test_that("plot.cmap_signature_result errors for 'heatmap' type", {
  skip_if_not_installed("ggplot2")
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )
  expect_error(plot(res, method = "ks", plot_type = "heatmap"),
               "not currently supported")
})

test_that("plot uses first available method when method = NULL", {
  skip_if_not_installed("ggplot2")
  skip_if_no_example()
  res <- process_signature_with_df(
    signature_file = sig_file(),
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = "xsum",
    save_files     = FALSE
  )
  # Should use "xsum" automatically
  expect_no_error(plot(res))
})

test_that("plot errors when no results exist", {
  skip_if_not_installed("ggplot2")
  res <- structure(list(results = list()), class = "cmap_signature_result")
  expect_error(plot(res), "No results available")
})

# ── create_signature_from_gene_lists ─────────────────────────────────────────

test_that("create_signature_from_gene_lists returns correct data frame", {
  sig <- create_signature_from_gene_lists(
    up_genes   = c("TP53", "MYC"),
    down_genes = c("EGFR", "KRAS")
  )
  expect_s3_class(sig, "data.frame")
  expect_named(sig, c("Gene", "log2FC"))
  expect_equal(nrow(sig), 4)
  expect_true(all(sig$log2FC[sig$Gene %in% c("TP53", "MYC")] == 1))
  expect_true(all(sig$log2FC[sig$Gene %in% c("EGFR", "KRAS")] == -1))
})

test_that("create_signature_from_gene_lists accepts custom values", {
  sig <- create_signature_from_gene_lists(
    up_genes   = c("A", "B"),
    down_genes = c("C"),
    up_value   = 2.5,
    down_value = -3
  )
  expect_equal(sig$log2FC[sig$Gene == "A"], 2.5)
  expect_equal(sig$log2FC[sig$Gene == "C"], -3)
})

test_that("create_signature_from_gene_lists errors when both lists are empty", {
  expect_error(
    create_signature_from_gene_lists(up_genes = character(0), down_genes = character(0)),
    "non-empty"
  )
})

test_that("create_signature_from_gene_lists warns on overlap", {
  expect_warning(
    create_signature_from_gene_lists(
      up_genes   = c("TP53", "MYC"),
      down_genes = c("MYC", "EGFR")
    ),
    "both up and down"
  )
})

test_that("create_signature_from_gene_lists warns when only up_genes provided", {
  expect_warning(
    create_signature_from_gene_lists(
      up_genes   = c("TP53"),
      down_genes = character(0)
    ),
    "No down-regulated"
  )
})

test_that("create_signature_from_gene_lists removes duplicates", {
  sig <- create_signature_from_gene_lists(
    up_genes   = c("TP53", "TP53", "MYC"),
    down_genes = c("EGFR", "EGFR")
  )
  expect_equal(nrow(sig), 3)
})

# ── process_signature_with_df with data.frame input ──────────────────────────

test_that("process_signature_with_df accepts a data.frame instead of file path", {
  skip_if_no_example()
  sig_data <- read.delim(sig_file())
  res <- process_signature_with_df(
    signature_file = sig_data,
    reference_df   = ref_df(),
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )
  expect_s3_class(res, "cmap_signature_result")
  expect_equal(res$settings$signature_file, "(data.frame)")
})

test_that("process_signature_with_df with gene-list-derived data.frame works", {
  skip_if_no_example()
  # Get gene names from the reference to ensure overlap
  rdf <- ref_df()
  gene_names <- rownames(rdf)
  sig <- create_signature_from_gene_lists(
    up_genes   = gene_names[1:5],
    down_genes = gene_names[6:10]
  )
  res <- process_signature_with_df(
    signature_file = sig,
    reference_df   = rdf,
    permutations   = 5,
    methods        = "ks",
    save_files     = FALSE
  )
  expect_s3_class(res, "cmap_signature_result")
  expect_gt(nrow(res$results$ks), 0)
})

test_that("process_signature_with_df errors on invalid signature_file type", {
  expect_error(
    process_signature_with_df(
      signature_file = 42,
      reference_df   = ref_df(),
      permutations   = 5,
      methods        = "ks",
      save_files     = FALSE
    ),
    "file path.*data.frame"
  )
})
