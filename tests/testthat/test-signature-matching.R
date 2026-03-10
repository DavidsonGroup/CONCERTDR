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

# ── Internal scoring method tests ────────────────────────────────────────────

test_that("score_ks returns expected structure", {
  set.seed(1234)
  ref <- matrix(rnorm(1000), nrow = 10,
                dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
  Up <- c("gene1", "gene2")
  Down <- c("gene9", "gene10")

  results <- CONCERTDR:::score_ks(ref, Up, Down, permuteNum = 100)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
  expect_true(all(c("Score", "pValue", "pAdjValue") %in% names(results)))
  expect_true(all(results$pValue >= 0 & results$pValue <= 1))
  expect_true(all(results$pAdjValue >= 0 & results$pAdjValue <= 1))
})

test_that("score_xcos returns expected structure", {
  set.seed(1234)
  ref <- matrix(rnorm(1000), nrow = 10,
                dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
  query <- rnorm(4)
  names(query) <- c("gene1", "gene2", "gene9", "gene10")

  results <- CONCERTDR:::score_xcos(ref, query, topN = 4, permuteNum = 100)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
  # XCos scores should be bounded by [-1, 1] (cosine similarity)
  expect_true(all(results$Score >= -1 & results$Score <= 1, na.rm = TRUE))
})

test_that("score_xsum returns expected structure", {
  set.seed(1234)
  ref <- matrix(rnorm(1000), nrow = 10,
                dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
  Up <- c("gene1", "gene2")
  Down <- c("gene9", "gene10")

  results <- CONCERTDR:::score_xsum(ref, Up, Down, topN = 4, permuteNum = 100)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
})

test_that("score_gsea0 returns expected structure", {
  set.seed(1234)
  ref <- matrix(rnorm(1000), nrow = 10,
                dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
  Up <- c("gene1", "gene2")
  Down <- c("gene9", "gene10")

  results <- CONCERTDR:::score_gsea0(ref, Up, Down, permuteNum = 100)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
})

test_that("score_gsea1 returns expected structure", {
  set.seed(1234)
  ref <- matrix(rnorm(1000), nrow = 10,
                dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
  Up <- c("gene1", "gene2")
  Down <- c("gene9", "gene10")

  results <- CONCERTDR:::score_gsea1(ref, Up, Down, permuteNum = 100)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
})

test_that("score_gsea2 returns expected structure", {
  set.seed(1234)
  ref <- matrix(rnorm(1000), nrow = 10,
                dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
  Up <- c("gene1", "gene2")
  Down <- c("gene9", "gene10")

  results <- CONCERTDR:::score_gsea2(ref, Up, Down, permuteNum = 100)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
})

test_that("score_zhang returns expected structure", {
  set.seed(1234)
  ref <- matrix(rnorm(1000), nrow = 10,
                dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
  Up <- c("gene1", "gene2")
  Down <- c("gene9", "gene10")

  results <- CONCERTDR:::score_zhang(ref, Up, Down, permuteNum = 100)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
  # Zhang scores should be bounded by [-1, 1]
  expect_true(all(abs(results$Score) <= 1, na.rm = TRUE))
})

test_that("score_xcos errors on bad input", {
  ref <- matrix(rnorm(100), nrow = 10,
                dimnames = list(paste0("gene", 1:10), paste0("drug", 1:10)))
  expect_error(CONCERTDR:::score_xcos(ref, "not_numeric", topN = 4),
               "query must be a numeric")
  expect_error(CONCERTDR:::score_xcos(ref, c(1, 2), topN = 4),
               "query must have names")
  q <- c(1, 2)
  names(q) <- c("gene1", "gene2")
  expect_error(CONCERTDR:::score_xcos(ref, q, topN = 6),
               "topN is larger")
})

test_that("score_ks handles missing genes gracefully", {
  set.seed(42)
  ref <- matrix(rnorm(100), nrow = 10,
                dimnames = list(paste0("gene", 1:10), paste0("drug", 1:10)))
  Up <- c("gene1", "missing_gene_A")
  Down <- c("gene10", "missing_gene_B")

  results <- CONCERTDR:::score_ks(ref, Up, Down, permuteNum = 10)
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 10)
})

# ── Summary / print / plot tests ─────────────────────────────────────────────

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

# ── End-to-end integration test ──────────────────────────────────────────────

test_that("process_signature_with_df runs end-to-end with internal methods", {
  sig_file <- system.file("extdata", "example_signature.txt", package = "CONCERTDR")
  ref_file <- system.file("extdata", "example_reference_df.csv", package = "CONCERTDR")

  skip_if(!nzchar(sig_file), "Example signature file not found")
  skip_if(!nzchar(ref_file), "Example reference file not found")

  ref_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
  ref_df$gene_symbol <- rownames(ref_df)

  results <- process_signature_with_df(
    signature_file = sig_file,
    reference_df = ref_df,
    output_dir = tempdir(),
    permutations = 10,
    methods = c("ks", "xsum"),
    save_files = FALSE
  )

  expect_s3_class(results, "cmap_signature_result")
  expect_true("ks" %in% names(results$results))
  expect_true("xsum" %in% names(results$results))
  expect_true(nrow(results$results$ks) > 0)
  expect_true(nrow(results$results$xsum) > 0)
})