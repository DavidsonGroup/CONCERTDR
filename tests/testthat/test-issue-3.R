# Regression coverage for GitHub issue #3.
issue3_fixture <- function(root = "0") {
  f <- tempfile(fileext = ".gctx")
  rhdf5::h5createFile(f)
  groups <- c("DATA", "DATA/0", "META", "META/ROW", "META/COL")
  if (nzchar(root)) {
    rhdf5::h5createGroup(f, root)
    groups <- paste0(root, "/", groups)
  }
  for (g in groups) rhdf5::h5createGroup(f, g)
  prefix <- if (nzchar(root)) paste0(root, "/") else ""
  m <- matrix(c(10, 2, -1, 4, 20, -2, 3, 1, 30, 1, -3, 2,
                40, 3, 1, -2, 50, -1, 2, 3), 4,
              dimnames = list(c("117153", "4253", "1", "2"), paste0("S", 1:5)))
  rhdf5::h5write(unname(m), f, paste0(prefix, "DATA/0/matrix"))
  rhdf5::h5write(rownames(m), f, paste0(prefix, "META/ROW/id"))
  rhdf5::h5write(colnames(m), f, paste0(prefix, "META/COL/id"))
  gi <- data.frame(gene_id = c("117153", "1", "4253", "2"),
                   gene_symbol = c("MIA2", "DOWN", "MIA2", "OTHER"),
                   feature_space = c("inferred", "landmark", "best inferred", "landmark"))
  si <- data.frame(sig_id = rev(colnames(m)), is_hiq = 1, pert_type = "trt_oe")
  list(file = f, matrix = m, gi = gi, si = si)
}

test_that("full extraction resolves duplicate symbols to the preferred feature", {
  f <- issue3_fixture()
  on.exit(unlink(f$file))
  expect_warning(out <- extract_cmap_data_from_siginfo(f$si, f$gi, f$file,
    landmark = FALSE, verbose = FALSE), "MIA2.*dropped gene_id=117153; retained gene_id=4253")
  expect_equal(out$gene_symbol, c("DOWN", "MIA2", "OTHER"))
  expected <- f$matrix[c("1", "4253", "2"), f$si$sig_id, drop = FALSE]
  expect_equal(unname(as.matrix(out[, -1])), unname(expected))
  expect_equal(attr(out, "metadata")$sample_id, f$si$sig_id)
  expect_no_warning(lm <- extract_cmap_data_from_siginfo(f$si, f$gi, f$file,
    landmark = TRUE, verbose = FALSE))
  expect_equal(lm$gene_symbol, c("DOWN", "OTHER"))
  sig <- data.frame(Gene = c("MIA2", "DOWN"), log2FC = c(1, -1))
  res <- suppressMessages(process_signature_with_df(sig, out, methods = "ks", permutations = 2))
  expect_false("error" %in% names(res$results$ks))
  expect_equal(nrow(res$results$ks), 5)
})

test_that("duplicate priority and tie handling preserve retained input order", {
  gi <- data.frame(gene_id = 1:7, gene_symbol = c("A", "B", "A", "B", "A", "C", "C"),
    feature_space = c("inferred", "other", "best inferred", "inferred", "landmark", NA, "other"))
  expect_warning(x <- CONCERTDR:::get_rid(gi, landmark = FALSE), "Duplicate gene symbols")
  expect_equal(x$rid, c("4", "5", "6"))
  expect_equal(x$genenames, c("B", "A", "C"))
  gi$feature_space <- NULL
  expect_warning(x <- CONCERTDR:::get_rid(gi, landmark = FALSE), "ties keep the first row")
  expect_equal(x$rid, c("1", "2", "6"))
})

test_that("direct heatmap extraction uses the same preferred feature", {
  f <- issue3_fixture()
  gf <- tempfile(fileext = ".tsv")
  write.table(f$gi, gf, sep = "\t", quote = FALSE, row.names = FALSE)
  on.exit(unlink(c(f$file, gf)))
  expect_warning(z <- extract_signature_zscores(
    data.frame(sig_id = "S1", Score = -1),
    data.frame(Gene = c("MIA2", "DOWN"), log2FC = c(1, -1)),
    gctx_file = f$file, geneinfo_file = gf, verbose = FALSE), "MIA2.*4253")
  expect_equal(unname(z$z_plot[1, "MIA2"]), f$matrix["4253", "S1"])
})

test_that("CamSum library rho uses the same duplicate resolution as extraction", {
  f <- issue3_fixture()
  gf <- tempfile(); sf <- tempfile()
  write.table(f$gi, gf, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(f$si, sf, sep = "\t", quote = FALSE, row.names = FALSE)
  old <- options(CONCERTDR.gctx_file = f$file, CONCERTDR.geneinfo_file = gf,
                 CONCERTDR.siginfo_file = sf)
  on.exit({options(old); unlink(c(f$file, gf, sf))})
  ref <- f$matrix[c("1", "4253", "2"), , drop = FALSE]
  rownames(ref) <- c("DOWN", "MIA2", "OTHER")
  expected <- compute_camsum_rho(ref, "MIA2", "DOWN", chunk = 2)
  expect_warning(actual <- suppressMessages(CONCERTDR:::score_camsum(
    ref, "MIA2", "DOWN", chunk = 2)), "MIA2.*4253")
  expect_equal(unname(attr(actual, "rho_bar")), as.numeric(expected), tolerance = 1e-12)
  expect_equal(unname(attr(actual, "rho_source")), "library_gctx")
})

test_that("GCTX reads sorted unique indices and restores both requested axes", {
  f <- issue3_fixture()
  on.exit(unlink(f$file))
  reader <- rhdf5::h5read
  local_mocked_bindings(h5read = function(file, name, ...) {
    args <- list(...)
    if (grepl("matrix$", name)) {
      for (idx in args$index) {
        if (!is.null(idx)) expect_identical(idx, sort(unique(idx)))
      }
    }
    reader(file, name, ...)
  }, .package = "rhdf5")
  rid <- c("2", "117153", "4253", "2", "not_found")
  cid <- c("S4", "S1", "S4", "absent")
  out <- CONCERTDR:::fast_parse_gctx(f$file, rid, cid)
  expect_equal(out, f$matrix[rid[1:4], cid[1:3], drop = FALSE])
  one <- CONCERTDR:::fast_parse_gctx(f$file, "4253", "S3")
  expect_equal(one, f$matrix["4253", "S3", drop = FALSE])
  all_rows <- CONCERTDR:::.gctx_read_matrix(f$file, list(NULL, c(5L, 1L, 5L)))
  expect_equal(all_rows, unname(f$matrix[, c(5, 1, 5), drop = FALSE]))
  expect_equal(dim(CONCERTDR:::fast_parse_gctx(f$file, "absent", "S1")), c(0L, 1L))
})

test_that("the alternate GCTX layout also restores requested order", {
  f <- issue3_fixture(root = "")
  on.exit(unlink(f$file))
  out <- suppressMessages(CONCERTDR:::fast_parse_gctx(f$file, c("2", "4253"), c("S5", "S1")))
  expect_equal(out, f$matrix[c("2", "4253"), c("S5", "S1"), drop = FALSE])
})

test_that("compound targets and MOAs are merged without multiplying result rows", {
  ids <- c("L_A_6H:D:1", "L_B_24H:D:10", "L_C_24H:E:1")
  res <- data.frame(compound = ids, Score = c(-.2, -.9, -.3))
  si <- data.frame(sig_id = ids, pert_id = c("D", "D", "E"),
                   pert_type = "trt_cp", pert_iname = c("drugD", "drugD", "drugE"))
  cp <- data.frame(pert_id = c(rep("D", 5), "E", "E"),
    target = c("A", "B", "A", NA, " ", NA, ""),
    moa = c("inhibitor", "agonist", "inhibitor", "", NA, NA, " "))
  for (input in list(cp, cp[, c("pert_id", "target")])) {
    v <- annotate_drug_results(res, si, input, verbose = FALSE)
    expect_equal(nrow(v$tech_view_all), nrow(res))
    expect_equal(v$tech_view_all$target[v$tech_view_all$pert_id == "D"], rep("A; B", 2))
    expect_true(is.na(v$tech_view_all$target[v$tech_view_all$pert_id == "E"]))
    expect_equal(v$drug_context_summary$n_contexts[v$drug_context_summary$pert_id == "D"], 2L)
    expect_equal(v$wetlab_drug_view$target[v$wetlab_drug_view$pert_id == "D"], "A; B")
    expect_equal(v$wetlab_drug_view$Score[v$wetlab_drug_view$pert_id == "D"], -.9)
    if ("moa" %in% names(input)) {
      expect_equal(v$tech_view_all$moa[v$tech_view_all$pert_id == "D"], rep("inhibitor; agonist", 2))
    }
  }
  f <- tempfile(); on.exit(unlink(f))
  write.table(cp, f, sep = "\t", quote = FALSE, row.names = FALSE)
  from_file <- annotate_drug_results(res, si, f, verbose = FALSE)
  from_df <- annotate_drug_results(res, si, cp, verbose = FALSE)
  expect_equal(from_file, from_df)
})
