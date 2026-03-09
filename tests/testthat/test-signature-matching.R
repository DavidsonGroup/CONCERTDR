make_mock_signature_result <- function() {
  structure(
    list(
      results = list(
        ks = data.frame(
          compound = c("imatinib", "dasatinib", "doxorubicin"),
          Score = c(-0.72, -0.45, 0.31),
          pValue = c(0.002, 0.041, 0.280),
          pAdjValue = c(0.020, 0.205, 0.560),
          stringsAsFactors = FALSE
        )
      ),
      summary = data.frame(
        compound = c("imatinib", "dasatinib", "doxorubicin"),
        method = c("ks", "ks", "ks"),
        Score = c(-0.72, -0.45, 0.31),
        pValue = c(0.002, 0.041, 0.280),
        rank = 1:3,
        global_rank = 1:3,
        stringsAsFactors = FALSE
      ),
      gene_data = data.frame(
        Gene = c("TP53", "MYC", "EGFR", "STAT1"),
        log2FC = c(1.5, -2.1, 1.2, -1.4),
        stringsAsFactors = FALSE
      ),
      settings = list(
        signature_file = "example_signature.txt",
        time_completed = Sys.time(),
        time_taken_mins = 0.1,
        methods = "ks",
        permutations = 10
      ),
      common_genes = list(
        up = list(found = c("TP53", "EGFR"), count = 2L, percent = 100),
        down = list(found = c("MYC", "STAT1"), count = 2L, percent = 100)
      )
    ),
    class = "cmap_signature_result"
  )
}

test_that("create_summary_from_results aggregates and ranks methods", {
  res1 <- data.frame(
    compound = c("A", "B"),
    Score = c(0.9, 0.2),
    pValue = c(0.01, 0.20),
    stringsAsFactors = FALSE
  )
  res2 <- data.frame(
    compound = c("C", "D"),
    Score = c(0.8, 0.1),
    pValue = c(0.02, 0.50),
    stringsAsFactors = FALSE
  )

  out <- CONCERTDR:::create_summary_from_results(list(ks = res1, xsum = res2), top_n = 1)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
  expect_true(all(c("compound", "method", "Score", "pValue", "global_rank") %in% names(out)))
})

test_that("print method works and returns invisibly", {
  obj <- make_mock_signature_result()
  expect_invisible(print(obj))
})

test_that("summary method returns top compounds invisibly", {
  obj <- make_mock_signature_result()
  out <- summary(obj, top_n = 2)
  expect_equal(nrow(out), 2)
  expect_true(all(c("compound", "method", "Score", "pValue", "global_rank") %in% names(out)))
})

test_that("summary method returns NULL when no summary data", {
  obj <- make_mock_signature_result()
  obj$summary <- obj$summary[0, ]
  expect_null(summary(obj))
})

test_that("plot method returns ggplot object for scores", {
  skip_if_not_installed("ggplot2")
  obj <- make_mock_signature_result()
  p <- plot(obj, method = "ks", plot_type = "scores", top_n = 2)
  expect_s3_class(p, "ggplot")
})

test_that("plot method errors for invalid method and plot type", {
  obj <- make_mock_signature_result()
  expect_error(plot(obj, method = "bad_method"), "not found")
  expect_error(plot(obj, method = "ks", plot_type = "bad_plot"), "Invalid plot_type")
})

test_that("run_cmap_workflow composes extraction and matching", {
  fake_ref <- data.frame(gene_symbol = c("TP53", "MYC"), sig1 = c(1, -1), stringsAsFactors = FALSE)
  fake_result <- make_mock_signature_result()

  testthat::local_mocked_bindings(
    extract_cmap_data_from_config = function(...) fake_ref,
    process_signature_with_df = function(...) fake_result,
    .package = "CONCERTDR"
  )

  out <- run_cmap_workflow(
    config_file = "dummy_config.txt",
    signature_file = "dummy_signature.txt",
    verbose = FALSE
  )

  expect_equal(class(out), "cmap_signature_result")
})

test_that("process_signature_with_df reports missing RCSM dependency", {
  if (requireNamespace("RCSM", quietly = TRUE)) {
    skip("RCSM installed; missing-dependency branch not applicable")
  }

  sig_file <- system.file("extdata", "example_signature.txt", package = "CONCERTDR")
  ref_file <- system.file("extdata", "example_reference_df.csv", package = "CONCERTDR")
  ref_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
  ref_df$gene_symbol <- rownames(ref_df)

  expect_error(
    process_signature_with_df(
      signature_file = sig_file,
      reference_df = ref_df,
      save_files = FALSE
    ),
    "Package 'RCSM' is required"
  )
})
