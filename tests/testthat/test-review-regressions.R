# Behavioural regressions found during the 0.99.4 review.
review_scores <- function() {
  data.frame(compound = c("rev", "weak_rev", "mim", "weak_mim", "zero", "missing"),
             Score = c(-2, -1, 2, 1, 0, NA_real_), pValue = .01)
}

test_that("summary only orders by direction and compares relative method ranks", {
  x <- review_scores()
  y <- x
  y$Score <- y$Score * 100
  rev <- CONCERTDR:::create_summary_from_results(list(a = x, b = y))
  mim <- CONCERTDR:::create_summary_from_results(list(a = x), direction = "mimic")
  # Every finite score is kept whatever its sign; only the order changes
  expect_equal(rev$compound, rep(c("rev", "weak_rev", "zero", "weak_mim", "mim"), 2))
  expect_equal(mim$compound, c("mim", "weak_mim", "zero", "weak_rev", "rev"))
  # Scaling the scores of one method does not change its relative ranks
  expect_equal(rev$rank_percentile[rev$method == "a"], rev$rank_percentile[rev$method == "b"])
  expect_equal(CONCERTDR:::create_summary_from_results(list(a = x), top_n = 1)$compound, "rev")
  only_pos <- x[x$Score > 0 & !is.na(x$Score), ]
  expect_equal(nrow(CONCERTDR:::create_summary_from_results(list(a = only_pos))), 2L)
  none <- CONCERTDR:::create_summary_from_results(list(a = x[0, ]))
  expect_equal(nrow(none), 0L)
  expect_named(none, names(rev))
})

test_that("direction changes ranking and bar order without changing scores or p-values", {
  sig <- data.frame(Gene = c("A", "B"), log2FC = c(1, -1))
  ref <- data.frame(gene_symbol = c("A", "B", "C", "D"),
                    rev = c(-3, 3, 1, -1), mim = c(3, -3, -1, 1))
  set.seed(42)
  rev <- suppressMessages(process_signature_with_df(sig, ref, methods = "ks", permutations = 5))
  set.seed(42)
  mim <- suppressMessages(process_signature_with_df(sig, ref, methods = "ks", permutations = 5,
                                                    direction = "mimic"))
  expect_equal(rev$results$ks[, c("Score", "pValue", "pAdjValue")],
               mim$results$ks[, c("Score", "pValue", "pAdjValue")])
  expect_equal(rev$summary$compound, c("rev", "mim"))
  expect_equal(mim$summary$compound, c("mim", "rev"))
  expect_equal(rev$results$ks$rank, c(1, 2))
  expect_equal(mim$results$ks$rank, c(2, 1))
  expect_equal(rev$settings$direction, "reversal")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_equal(plot(rev)$data$compound, c("rev", "mim"))
    expect_equal(plot(mim)$data$compound, c("mim", "rev"))
  }
})

test_that("CamSum metadata errors do not break other methods", {
  bad <- tempfile()
  writeLines("wrong_column\nvalue", bad)
  on.exit(unlink(bad))
  old_options <- options(CONCERTDR.siginfo_file = bad)
  on.exit(options(old_options), add = TRUE)
  sig <- data.frame(Gene = c("A", "B"), log2FC = c(1, -1))
  ref <- data.frame(gene_symbol = c("A", "B"), s = c(-1, 1))
  expect_no_warning(out <- suppressMessages(process_signature_with_df(
    sig, ref, methods = "ks", permutations = 2)))
  expect_false("error" %in% names(out$results$ks))
  out <- suppressWarnings(suppressMessages(process_signature_with_df(
    sig, ref, methods = c("camsum", "ks"), permutations = 2)))
  expect_true("error" %in% names(out$results$camsum))
  expect_false("error" %in% names(out$results$ks))
})

test_that("drug views keep the best context, ordered by direction, without dropping any", {
  ids <- c("L_A_6H:D:1", "L_B_24H:D:10", "L_C_48H:D:20")
  res <- data.frame(compound = ids, Score = c(-.2, -.9, .7))
  si <- data.frame(sig_id = ids, pert_id = "D", pert_type = "trt_cp", pert_iname = "drug")
  cp <- data.frame(pert_id = "D")
  rev <- annotate_drug_results(res, si, cp, keep_dose_in_drug_view = TRUE, verbose = FALSE)
  mim <- annotate_drug_results(res, si, cp, keep_dose_in_drug_view = TRUE,
                               verbose = FALSE, direction = "mimic")
  expect_equal(rev$wetlab_drug_view$Score, -.9)
  expect_equal(rev$wetlab_drug_view$cell_line, "B")
  expect_equal(rev$wetlab_drug_view$time_h, 24)
  expect_equal(rev$wetlab_drug_view$dose_uM, 10)
  expect_equal(rev$drug_context_summary$n_contexts, 3L)
  expect_equal(mim$wetlab_drug_view$Score, .7)
  expect_equal(mim$wetlab_drug_view$cell_line, "C")
  expect_equal(mim$drug_context_summary$n_contexts, 3L)
  expect_equal(nrow(rev$tech_view_all), 3L)
  # Positive-only scores are still reported, not silently dropped
  res$Score <- abs(res$Score)
  pos <- annotate_drug_results(res, si, cp, verbose = FALSE)
  expect_equal(nrow(pos$wetlab_drug_view), 1L)
  expect_equal(pos$wetlab_drug_view$Score, .2)
})

test_that("regex extraction preserves unmatched positions and missing values", {
  x <- c("A:BRD-ABC:10", "NO_MATCH", "B:BRD-DEF:10", NA_character_, "")
  expect_equal(extract_compound_id(x, method = "regex", regex_pattern = "BRD-[A-Z]+"),
               c("BRD-ABC", "NO_MATCH", "BRD-DEF", NA_character_, ""))
})

test_that("table selection is passed to fread and accepts trimmed headers", {
  f <- tempfile()
  writeLines(c(" id \tunused\tvalue", "a\tx\t1", "b\ty\t2"), f)
  on.exit(unlink(f))
  reader <- data.table::fread
  selections <- list()
  local_mocked_bindings(fread = function(...) {
    args <- list(...)
    if (!is.null(args$select)) selections[[length(selections) + 1L]] <<- args$select
    reader(...)
  }, .package = "data.table")
  out <- CONCERTDR:::.read_cmap_table(f, "fixture", select = c("value", "id", "absent"))
  expect_named(out, c("value", "id"))
  expect_equal(out$id, c("a", "b"))
  expect_equal(out$value, 1:2)
  expect_length(selections, 1L)
  expect_length(selections[[1]], 2L)
})

review_gctx <- function() {
  f <- tempfile(fileext = ".gctx")
  rhdf5::h5createFile(f)
  for (g in c("0", "0/META", "0/META/ROW", "0/META/COL", "0/DATA", "0/DATA/0")) {
    rhdf5::h5createGroup(f, g)
  }
  rhdf5::h5write(c("1", "2"), f, "0/META/ROW/id")
  rhdf5::h5write(c("s1", "s2"), f, "0/META/COL/id")
  rhdf5::h5write(matrix(c(11, 21, 12, 22), 2), f, "0/DATA/0/matrix")
  f
}

test_that("GCTX extraction aligns gene names to actual returned IDs", {
  f <- review_gctx()
  on.exit(unlink(f))
  gi <- data.frame(gene_id = c("2", "3", "1"), gene_symbol = c("B", "MISSING", "A"),
                   feature_space = "landmark")
  si <- data.frame(sig_id = c("s2", "s1"), is_hiq = 1,
                   pert_itime = "24 h", pert_idose = "1 uM", cell_iname = "A")
  out <- extract_cmap_data_from_siginfo(si, gi, f, verbose = FALSE)
  expect_equal(out$gene_symbol, c("B", "A"))
  expect_equal(out$s2, c(22, 12))
  expect_equal(attr(out, "metadata")$sample_id, c("s2", "s1"))
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  path <- suppressMessages(CONCERTDR:::process_combination(
    data.frame(itime = "24 h", idose = "1 uM", cell = "A"),
    gi$gene_id, gi$gene_symbol, si, f, d))
  tab <- read.delim(path, row.names = 1)
  expect_equal(rownames(tab), c("B", "A"))
  expect_equal(tab$s2, c(22, 12))
  gf <- tempfile(); write.table(gi, gf, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(gf), add = TRUE)
  z <- extract_signature_zscores(data.frame(sig_id = c("s1", "absent"), Score = c(-1, -2)),
    data.frame(Gene = c("A", "MISSING", "B"), log2FC = c(1, 2, -1)),
    gctx_file = f, geneinfo_file = gf, verbose = FALSE)
  expect_equal(z$sig_ids, "s1")
  expect_equal(z$ordered_genes, c("B", "A"))
  expect_equal(unname(z$z_plot), matrix(c(21, 11), 1))
})

test_that("one-gene matrices keep dimensions and direction only orders extraction", {
  m <- matrix(c(1, 2), 1, dimnames = list("G", c("rev", "mim")))
  sig <- data.frame(Gene = "G", log2FC = 1)
  res <- data.frame(sig_id = c("rev", "mim"), Score = c(-1, 1))
  rev <- extract_signature_zscores(res, sig, reference_df = m, verbose = FALSE)
  mim <- extract_signature_zscores(res, sig, reference_df = m, verbose = FALSE, direction = "mimic")
  expect_equal(dim(rev$z_plot), c(2L, 1L))
  expect_equal(rev$sig_ids, c("rev", "mim"))
  expect_equal(unname(rev$z_plot[, 1]), c(1, 2))
  expect_equal(mim$sig_ids, c("mim", "rev"))
})

test_that("direct split plots select both directions and allow a single profile", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  sig <- data.frame(Gene = c("A", "B", "C", "D"), log2FC = c(4, 3, -4, -3))
  ref <- matrix(c(1, 2, -2, -1), 4, dimnames = list(sig$Gene, "s"))
  pdf <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf)
  on.exit({grDevices::dev.off(); unlink(pdf)})
  direct <- suppressMessages(plot_signature_direction_tile_barcode(
    results_df = data.frame(sig_id = "s", Score = -1), signature_file = sig,
    reference_df = ref, split_direction = TRUE, max_genes = 2, verbose = FALSE))
  expect_equal(direct$ordered_genes, c("C", "D", "A", "B"))
  z <- extract_signature_zscores(data.frame(sig_id = "s", Score = -1), sig,
                                 reference_df = ref, split_direction = TRUE,
                                 max_genes = 2, verbose = FALSE)
  pre <- suppressMessages(plot_signature_direction_tile_barcode(
    precomputed = z, split_direction = TRUE, verbose = FALSE))
  expect_equal(direct$z_plot, pre$z_plot)
})

test_that("reported query coverage uses the complete query", {
  sig <- data.frame(Gene = c("MISSING_UP", "A", "MISSING_DOWN", "B"),
                     log2FC = c(10, 1, -10, -1))
  ref <- data.frame(gene_symbol = c("A", "B"), s = c(-1, 1))
  out <- suppressMessages(process_signature_with_df(sig, ref, methods = "ks",
                                                     topN = 1, permutations = 2))
  expect_equal(out$common_genes$up$percent, 50)
  expect_equal(out$common_genes$down$percent, 50)
})
