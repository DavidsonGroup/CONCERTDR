# tests/testthat/test-cmap-extract.R
#
# Tests for get_rid() (internal) and extract_cmap_data_from_siginfo() using
# in-memory data frames (no GCTX file required).

# ── get_rid ───────────────────────────────────────────────────────────────────

make_geneinfo_df <- function() {
  data.frame(
    gene_id       = c("7157", "4609", "6774", "1956", "3845",
                       "207",  "596",  "1017",  "595",  "5594"),
    gene_symbol   = c("TP53", "MYC", "STAT3", "EGFR", "KRAS",
                       "AKT1", "BCL2", "CDK2",  "CCND1", "MAPK1"),
    feature_space = c("landmark", "landmark", "landmark", "landmark", "landmark",
                       "landmark", "best inferred", "best inferred", "landmark", "landmark"),
    stringsAsFactors = FALSE
  )
}

test_that("get_rid with landmark = TRUE returns only landmark genes", {
  gi  <- make_geneinfo_df()
  out <- CONCERTDR:::get_rid(gi, landmark = TRUE)

  expect_named(out, c("rid", "genenames"))
  # Only 8 landmark genes in our fixture
  expect_equal(length(out$rid), 8)
  expect_true(all(gi$gene_id[gi$feature_space == "landmark"] %in% out$rid))
})

test_that("get_rid with landmark = FALSE returns all genes", {
  gi  <- make_geneinfo_df()
  out <- CONCERTDR:::get_rid(gi, landmark = FALSE)
  expect_equal(length(out$rid), nrow(gi))
})

test_that("get_rid returns character vectors", {
  gi  <- make_geneinfo_df()
  out <- CONCERTDR:::get_rid(gi)
  expect_type(out$rid,       "character")
  expect_type(out$genenames, "character")
})

test_that("get_rid rid and genenames have equal length", {
  gi  <- make_geneinfo_df()
  out <- CONCERTDR:::get_rid(gi)
  expect_equal(length(out$rid), length(out$genenames))
})

test_that("get_rid returns empty vectors when no landmark genes exist", {
  gi  <- make_geneinfo_df()
  gi$feature_space <- "best inferred"   # no landmark genes
  out <- CONCERTDR:::get_rid(gi, landmark = TRUE)
  expect_equal(length(out$rid), 0)
})

# ── extract_cmap_data_from_siginfo: using df inputs (no GCTX required) ───────
#
# extract_cmap_data_from_siginfo always calls cmapR::parse_gctx internally, so
# we cannot fully test it without a GCTX file. However, we can test all the
# validation and pre-processing logic by checking error messages triggered
# before the GCTX call.

make_minimal_siginfo <- function() {
  data.frame(
    sig_id     = c("SIG001", "SIG002"),
    pert_type  = c("trt_cp", "trt_cp"),
    is_hiq     = c(1L, 1L),
    cell_iname = c("K562", "HL60"),
    pert_itime = c("6 h", "24 h"),
    pert_idose = c("10 uM", "1 uM"),
    stringsAsFactors = FALSE
  )
}

test_that("extract_cmap_data_from_siginfo errors if geneinfo_file path is missing", {
  expect_error(
    extract_cmap_data_from_siginfo(
      siginfo_file  = make_minimal_siginfo(),
      geneinfo_file = "no_such_geneinfo.txt",
      gctx_file     = "dummy.gctx",
      verbose       = FALSE
    ),
    "not found"
  )
})

test_that("extract_cmap_data_from_siginfo errors if siginfo_file path is missing", {
  tmp_gi <- tempfile(fileext = ".txt")
  write.table(make_geneinfo_df(), tmp_gi, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(tmp_gi))

  expect_error(
    extract_cmap_data_from_siginfo(
      siginfo_file  = "no_such_siginfo.txt",
      geneinfo_file = tmp_gi,
      gctx_file     = "dummy.gctx",
      verbose       = FALSE
    ),
    "not found"
  )
})

test_that("extract_cmap_data_from_siginfo errors if sig_id column is absent", {
  si_no_id <- make_minimal_siginfo()
  si_no_id$sig_id <- NULL

  tmp_gi <- tempfile(fileext = ".txt")
  write.table(make_geneinfo_df(), tmp_gi, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(tmp_gi))

  expect_error(
    extract_cmap_data_from_siginfo(
      siginfo_file  = si_no_id,
      geneinfo_file = tmp_gi,
      gctx_file     = "dummy.gctx",
      verbose       = FALSE
    ),
    "sig_id"
  )
})

test_that("extract_cmap_data_from_siginfo honours max_signatures", {
  # We expect the function to error trying to read a dummy GCTX,
  # but only after it has trimmed siginfo to max_signatures rows.
  # We intercept the error and check the message doesn't say "No signatures"
  tmp_gi <- tempfile(fileext = ".txt")
  write.table(make_geneinfo_df(), tmp_gi, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(tmp_gi))

  err <- tryCatch(
    extract_cmap_data_from_siginfo(
      siginfo_file    = make_minimal_siginfo(),
      geneinfo_file   = tmp_gi,
      gctx_file       = "dummy.gctx",
      max_signatures  = 1,
      filter_quality  = FALSE,
      verbose         = FALSE
    ),
    error = function(e) e
  )
  # The error should be about reading the GCTX, not about "No signatures"
  expect_false(grepl("No signatures found", conditionMessage(err)))
  expect_true(grepl("GCTX|gctx|Error reading", conditionMessage(err),
                    ignore.case = TRUE))
})

test_that("extract_cmap_data_from_siginfo filter_quality removes low-quality sigs", {
  # Provide a siginfo where only one sig has is_hiq = 1
  si <- make_minimal_siginfo()
  si$is_hiq[2] <- 0

  tmp_gi <- tempfile(fileext = ".txt")
  write.table(make_geneinfo_df(), tmp_gi, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(tmp_gi))

  # Expect GCTX error, but only after filtering to 1 signature
  err <- tryCatch(
    extract_cmap_data_from_siginfo(
      siginfo_file   = si,
      geneinfo_file  = tmp_gi,
      gctx_file      = "dummy.gctx",
      filter_quality = TRUE,
      verbose        = FALSE
    ),
    error = function(e) e
  )
  expect_false(grepl("No signatures found", conditionMessage(err)))
})

test_that("extract_cmap_data_from_siginfo errors when all sigs filtered out", {
  si <- make_minimal_siginfo()
  si$is_hiq <- 0L   # all low quality

  tmp_gi <- tempfile(fileext = ".txt")
  write.table(make_geneinfo_df(), tmp_gi, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(tmp_gi))

  expect_error(
    extract_cmap_data_from_siginfo(
      siginfo_file   = si,
      geneinfo_file  = tmp_gi,
      gctx_file      = "dummy.gctx",
      filter_quality = TRUE,
      verbose        = FALSE
    ),
    "No signatures found"
  )
})

test_that("extract_cmap_data_from_siginfo accepts geneinfo as a data.frame", {
  gi_df <- make_geneinfo_df()

  # Should reach the GCTX reading step (and fail there), not error on geneinfo
  err <- tryCatch(
    extract_cmap_data_from_siginfo(
      siginfo_file   = make_minimal_siginfo(),
      geneinfo_file  = gi_df,   # data.frame, not file path
      gctx_file      = "dummy.gctx",
      filter_quality = FALSE,
      verbose        = FALSE
    ),
    error = function(e) e
  )
  expect_true(grepl("GCTX|gctx|Error reading", conditionMessage(err),
                    ignore.case = TRUE))
})

test_that("extract_cmap_data_from_siginfo accepts siginfo as a data.frame", {
  tmp_gi <- tempfile(fileext = ".txt")
  write.table(make_geneinfo_df(), tmp_gi, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(tmp_gi))

  # Should reach the GCTX step
  err <- tryCatch(
    extract_cmap_data_from_siginfo(
      siginfo_file   = make_minimal_siginfo(),  # data.frame
      geneinfo_file  = tmp_gi,
      gctx_file      = "dummy.gctx",
      filter_quality = FALSE,
      verbose        = FALSE
    ),
    error = function(e) e
  )
  expect_true(grepl("GCTX|gctx|Error reading", conditionMessage(err),
                    ignore.case = TRUE))
})

test_that("extract_cmap_data_from_siginfo errors on wrong geneinfo_file type", {
  expect_error(
    extract_cmap_data_from_siginfo(
      siginfo_file   = make_minimal_siginfo(),
      geneinfo_file  = 123,   # not a string or df
      gctx_file      = "dummy.gctx",
      verbose        = FALSE
    ),
    "geneinfo_file must be"
  )
})

test_that("extract_cmap_data_from_siginfo errors on wrong siginfo_file type", {
  expect_error(
    extract_cmap_data_from_siginfo(
      siginfo_file  = 123,   # not a string or df
      geneinfo_file = make_geneinfo_df(),
      gctx_file     = "dummy.gctx",
      verbose       = FALSE
    ),
    "siginfo_file must be"
  )
})
