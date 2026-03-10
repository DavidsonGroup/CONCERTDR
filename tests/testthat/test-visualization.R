test_that("extract_signature_zscores works with reference_df fallback", {
  sig_file <- system.file("extdata", "example_signature.txt", package = "CONCERTDR")
  ref_file <- system.file("extdata", "example_reference_df.csv", package = "CONCERTDR")

  if (!nzchar(sig_file)) {
    sig_file <- testthat::test_path("..", "..", "inst", "extdata", "example_signature.txt")
  }
  if (!nzchar(ref_file)) {
    ref_file <- testthat::test_path("..", "..", "inst", "extdata", "example_reference_df.csv")
  }

  skip_if(!nzchar(sig_file), "Example signature file not found")
  skip_if(!nzchar(ref_file), "Example reference file not found")

  ref_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
  ref_df$gene_symbol <- rownames(ref_df)

  results_df <- data.frame(
    compound = c("DEMO001", "DEMO002", "DEMO003"),
    Score = c(-0.72, -0.45, 0.31),
    stringsAsFactors = FALSE
  )

  z <- extract_signature_zscores(
    results_df = results_df,
    signature_file = sig_file,
    reference_df = ref_df,
    max_genes = 10,
    max_perts = 3,
    verbose = FALSE
  )

  expect_true(is.list(z))
  expect_true(all(c("z_plot", "ordered_genes", "logfc_map", "sig_ids") %in% names(z)))
  expect_equal(nrow(z$z_plot), 3)
  expect_equal(colnames(z$z_plot), z$ordered_genes)
  expect_equal(unname(z$sig_ids), c("DEMO001", "DEMO002", "DEMO003"))
})
