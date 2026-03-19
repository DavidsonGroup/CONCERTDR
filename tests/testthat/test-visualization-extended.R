# tests/testthat/test-visualization-extended.R
#
# Extended tests for extract_signature_zscores() and
# plot_signature_direction_tile_barcode(), plus subset_siginfo_beta().

# ── Helpers ──────────────────────────────────────────────────────────────────

sig_file <- function() {
  system.file("extdata", "example_signature.txt", package = "CONCERTDR")
}

ref_df <- function() {
  f  <- system.file("extdata", "example_reference_df.csv", package = "CONCERTDR")
  df <- read.csv(f, row.names = 1, check.names = FALSE)
  df$gene_symbol <- rownames(df)
  df
}

results_df <- function() {
  # Derive compound IDs directly from the reference matrix so that this
  # fixture remains valid regardless of which example dataset is installed.
  sig_ids <- setdiff(colnames(ref_df()), "gene_symbol")
  n       <- length(sig_ids)
  data.frame(
    compound = sig_ids,
    Score    = seq(-1, 1, length.out = n),
    stringsAsFactors = FALSE
  )
}

skip_if_no_example <- function() {
  skip_if(!nzchar(sig_file()), "Example signature file not found")
}

# ── extract_signature_zscores: input validation ───────────────────────────────

test_that("extract_signature_zscores errors if results_df is not a data.frame", {
  skip_if_no_example()
  expect_error(
    extract_signature_zscores(
      results_df     = "not_a_df",
      signature_file = sig_file(),
      reference_df   = ref_df()
    ),
    "data.frame"
  )
})


test_that("extract_signature_zscores errors if score_col is absent", {
  skip_if_no_example()
  expect_error(
    extract_signature_zscores(
      results_df     = results_df(),
      signature_file = sig_file(),
      reference_df   = ref_df(),
      score_col      = "MISSING_SCORE"
    ),
    "score_col"
  )
})

test_that("extract_signature_zscores errors if signature_file does not exist", {
  skip_if_no_example()
  expect_error(
    extract_signature_zscores(
      results_df     = results_df(),
      signature_file = "no_such_file.txt",
      reference_df   = ref_df()
    )
  )
})

test_that("extract_signature_zscores errors when no signature genes match reference", {
  skip_if_no_example()
  # Replace all row names in reference with fake gene names
  bad_ref <- ref_df()
  rownames(bad_ref) <- paste0("FAKE_", rownames(bad_ref))
  bad_ref$gene_symbol <- rownames(bad_ref)

  expect_error(
    extract_signature_zscores(
      results_df     = results_df(),
      signature_file = sig_file(),
      reference_df   = bad_ref,
      verbose        = FALSE
    ),
    "No signature genes matched"
  )
})

# ── extract_signature_zscores: data.frame signature input ─────────────────────

test_that("extract_signature_zscores accepts a data.frame for signature_file", {
  skip_if_no_example()
  sig_data <- read.delim(sig_file())
  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_data,
    reference_df   = ref_df(),
    verbose        = FALSE
  )
  expect_type(z, "list")
  expect_true("z_plot" %in% names(z))
  expect_gt(ncol(z$z_plot), 0)
})

test_that("extract_signature_zscores errors on invalid signature_file type", {
  skip_if_no_example()
  expect_error(
    extract_signature_zscores(
      results_df     = results_df(),
      signature_file = 42,
      reference_df   = ref_df(),
      verbose        = FALSE
    ),
    "file path.*data.frame"
  )
})

# ── extract_signature_zscores: return structure ───────────────────────────────

test_that("extract_signature_zscores returns a list with five named elements", {
  skip_if_no_example()
  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    verbose        = FALSE
  )
  expect_type(z, "list")
  expect_named(z, c("z_plot", "ordered_genes", "logfc_map", "sig_ids", "sig_labels"),
               ignore.order = TRUE)
})

test_that("z_plot dimensions match perturbation count and gene count", {
  skip_if_no_example()
  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    max_genes      = 15,
    max_perts      = 5,
    verbose        = FALSE
  )
  expect_lte(nrow(z$z_plot), 5)
  expect_lte(ncol(z$z_plot), 15)
})

test_that("column names of z_plot match ordered_genes", {
  skip_if_no_example()
  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    verbose        = FALSE
  )
  expect_equal(colnames(z$z_plot), z$ordered_genes)
})

test_that("logfc_map is a named numeric vector", {
  skip_if_no_example()
  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    verbose        = FALSE
  )
  expect_true(is.numeric(z$logfc_map))
  expect_false(is.null(names(z$logfc_map)))
})

test_that("ordered_genes starts with down-regulated then up-regulated genes", {
  skip_if_no_example()
  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    verbose        = FALSE
  )
  # logfc_map for first half should be <= 0, last half >= 0
  lfc <- z$logfc_map[z$ordered_genes]
  n_down <- sum(lfc < 0)
  n_up   <- sum(lfc > 0)
  # first n_down genes should all have negative logFC
  expect_true(all(lfc[seq_len(n_down)] < 0))
  # last n_up genes should all have positive logFC
  expect_true(all(lfc[(n_down + 1):(n_down + n_up)] > 0))
})

test_that("sig_ids are a subset of input compound IDs", {
  skip_if_no_example()
  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    verbose        = FALSE
  )
  expect_true(all(z$sig_ids %in% results_df()$compound))
})

test_that("perturbations are ordered by score ascending (worst to best)", {
  skip_if_no_example()
  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    max_perts      = 10,
    verbose        = FALSE
  )
  # Row order should correspond to ascending Score
  rd  <- results_df()
  rd  <- rd[order(rd$Score), ]
  expected_ids <- rd$compound[rd$compound %in% z$sig_ids]
  expect_equal(z$sig_ids, expected_ids)
})

# ── extract_signature_zscores: selected_drug filter ──────────────────────────

test_that("selected_drug filters to matching rows", {
  skip_if_no_example()
  # Create a results_df with a drug name column
  rd_with_drug <- results_df()
  rd_with_drug$drug_name <- rep(c("imatinib", "dasatinib"),
                                length.out = nrow(rd_with_drug))

  z_all    <- extract_signature_zscores(
    results_df     = rd_with_drug,
    signature_file = sig_file(),
    reference_df   = ref_df(),
    max_perts      = 10,
    verbose        = FALSE
  )
  z_subset <- extract_signature_zscores(
    results_df        = rd_with_drug,
    signature_file    = sig_file(),
    reference_df      = ref_df(),
    selected_drug     = "imatinib",
    selected_drug_col = "drug_name",
    max_perts         = 10,
    verbose           = FALSE
  )

  expect_lte(nrow(z_subset$z_plot), nrow(z_all$z_plot))
})

# ── extract_signature_zscores: output_zscores ────────────────────────────────

test_that("output_zscores writes a readable TSV when path is provided", {
  skip_if_no_example()
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp))

  extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    output_zscores = tmp,
    verbose        = FALSE
  )

  expect_true(file.exists(tmp))
  saved <- read.delim(tmp, row.names = 1, check.names = FALSE)
  expect_gt(nrow(saved), 0)
  expect_gt(ncol(saved), 0)
})

# ── extract_signature_zscores: matrix accepts both df and matrix input ────────

test_that("reference_df accepts a plain matrix (no gene_symbol column)", {
  skip_if_no_example()
  rd  <- ref_df()
  rd$gene_symbol <- NULL
  mat <- as.matrix(rd)

  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = mat,
    verbose        = FALSE
  )
  expect_gt(ncol(z$z_plot), 0)
})

# ── plot_signature_direction_tile_barcode ─────────────────────────────────────

test_that("plot_signature_direction_tile_barcode errors without ComplexHeatmap", {
  skip_if_no_example()
  skip_if(requireNamespace("ComplexHeatmap", quietly = TRUE),
          "ComplexHeatmap is installed — skipping absence test")

  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    verbose        = FALSE
  )
  expect_error(
    plot_signature_direction_tile_barcode(precomputed = z),
    "ComplexHeatmap"
  )
})

test_that("plot_signature_direction_tile_barcode errors without results_df when no precomputed", {
  skip_if_no_example()
  expect_error(
    plot_signature_direction_tile_barcode(
      signature_file = sig_file(),
      reference_df   = ref_df()
    ),
    "results_df"
  )
})

test_that("plot_signature_direction_tile_barcode errors without signature_file when no precomputed", {
  skip_if_no_example()
  expect_error(
    plot_signature_direction_tile_barcode(
      results_df   = results_df(),
      reference_df = ref_df()
    ),
    "signature_file"
  )
})

test_that("plot_signature_direction_tile_barcode returns invisible list", {
  skip_if_no_example()
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  z   <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    max_genes      = 10,
    max_perts      = 5,
    verbose        = FALSE
  )
  out <- plot_signature_direction_tile_barcode(
    precomputed  = z,
    save_png     = FALSE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    verbose      = FALSE
  )
  expect_type(out, "list")
  expect_named(out, c("z_plot", "ordered_genes", "sig_ids", "output_png", "output_zscores"),
               ignore.order = TRUE)
})

test_that("plot_signature_direction_tile_barcode saves PNG when save_png = TRUE", {
  skip_if_no_example()
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))

  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    max_genes      = 10,
    max_perts      = 5,
    verbose        = FALSE
  )
  plot_signature_direction_tile_barcode(
    precomputed  = z,
    save_png     = TRUE,
    output_png   = tmp,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    verbose      = FALSE
  )
  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0)
})

test_that("precomputed argument skips data loading", {
  skip_if_no_example()
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  z <- extract_signature_zscores(
    results_df     = results_df(),
    signature_file = sig_file(),
    reference_df   = ref_df(),
    verbose        = FALSE
  )
  # Passing precomputed should work even when other data args are NULL
  expect_no_error(
    plot_signature_direction_tile_barcode(
      precomputed  = z,
      save_png     = FALSE,
      verbose      = FALSE
    )
  )
})

test_that("plot_signature_direction_tile_barcode errors on malformed precomputed", {
  expect_error(
    plot_signature_direction_tile_barcode(
      precomputed = list(wrong_name = 1)
    ),
    "precomputed must be the output"
  )
})

# ── subset_siginfo_beta ───────────────────────────────────────────────────────

make_siginfo_file <- function() {
  df <- data.frame(
    sig_id     = paste0("SIG", 1:12),
    pert_type  = rep(c("trt_cp", "trt_cp", "trt_sh", "ctl_vehicle"), 3),
    pert_itime = rep(c("6 h", "24 h"), 6),
    pert_idose = rep(c("1 uM", "10 uM", "0 uM"), 4),
    cell_iname = rep(c("A375", "MCF7", "K562", "HL60"), 3),
    stringsAsFactors = FALSE
  )
  f <- tempfile(fileext = ".txt")
  write.table(df, f, sep = "\t", row.names = FALSE, quote = FALSE)
  f
}

test_that("subset_siginfo_beta returns a data.frame in non-interactive mode", {
  f   <- make_siginfo_file()
  on.exit(unlink(f))
  out <- subset_siginfo_beta(f, interactive = FALSE, verbose = FALSE,
                             show_preview = FALSE)
  expect_s3_class(out, "data.frame")
})

test_that("subset_siginfo_beta filters correctly with a single filter", {
  f   <- make_siginfo_file()
  on.exit(unlink(f))
  out <- subset_siginfo_beta(
    f,
    interactive = FALSE,
    filters     = list(pert_type = "trt_cp"),
    verbose     = FALSE,
    show_preview = FALSE
  )
  expect_true(all(out$pert_type == "trt_cp"))
})

test_that("subset_siginfo_beta filters correctly with multiple filters", {
  f   <- make_siginfo_file()
  on.exit(unlink(f))
  out <- subset_siginfo_beta(
    f,
    interactive  = FALSE,
    filters      = list(pert_type  = "trt_cp",
                        pert_itime = "6 h"),
    verbose      = FALSE,
    show_preview = FALSE
  )
  expect_true(all(out$pert_type == "trt_cp"))
  expect_true(all(out$pert_itime == "6 h"))
})

test_that("subset_siginfo_beta returns all rows when filters = NULL", {
  f    <- make_siginfo_file()
  on.exit(unlink(f))
  orig <- read.delim(f, stringsAsFactors = FALSE)
  out  <- subset_siginfo_beta(
    f,
    interactive  = FALSE,
    filters      = NULL,
    verbose      = FALSE,
    show_preview = FALSE
  )
  expect_equal(nrow(out), nrow(orig))
})

test_that("subset_siginfo_beta accepts multi-value filters", {
  f   <- make_siginfo_file()
  on.exit(unlink(f))
  out <- subset_siginfo_beta(
    f,
    interactive  = FALSE,
    filters      = list(cell_iname = c("A375", "MCF7")),
    verbose      = FALSE,
    show_preview = FALSE
  )
  expect_true(all(out$cell_iname %in% c("A375", "MCF7")))
})

test_that("subset_siginfo_beta returns zero rows when filter matches nothing", {
  f   <- make_siginfo_file()
  on.exit(unlink(f))
  out <- subset_siginfo_beta(
    f,
    interactive  = FALSE,
    filters      = list(cell_iname = "NONEXISTENT_CELL"),
    verbose      = FALSE,
    show_preview = FALSE
  )
  expect_equal(nrow(out), 0)
})

test_that("subset_siginfo_beta saves output file when output_file is specified", {
  f   <- make_siginfo_file()
  out_f <- tempfile(fileext = ".txt")
  on.exit({ unlink(f); unlink(out_f) })

  subset_siginfo_beta(
    f,
    output_file  = out_f,
    interactive  = FALSE,
    filters      = list(pert_type = "trt_cp"),
    verbose      = FALSE,
    show_preview = FALSE
  )
  expect_true(file.exists(out_f))
  saved <- read.delim(out_f, stringsAsFactors = FALSE)
  expect_gt(nrow(saved), 0)
  expect_true(all(saved$pert_type == "trt_cp"))
})

test_that("subset_siginfo_beta errors when siginfo_file does not exist", {
  expect_error(
    subset_siginfo_beta("no_such_file.txt", interactive = FALSE,
                        verbose = FALSE, show_preview = FALSE),
    "not found"
  )
})

test_that("subset_siginfo_beta silently skips unrecognised filter columns", {
  f <- make_siginfo_file()
  on.exit(unlink(f))
  out <- subset_siginfo_beta(
    f,
    interactive  = FALSE,
    filters      = list(nonexistent_col = "value"),
    verbose      = FALSE,
    show_preview = FALSE
  )
  expect_s3_class(out, "data.frame")
  orig <- read.delim(f, stringsAsFactors = FALSE)
  expect_equal(nrow(out), nrow(orig))
})

test_that("subset_siginfo_beta column names include key metadata columns", {
  f   <- make_siginfo_file()
  on.exit(unlink(f))
  out <- subset_siginfo_beta(f, interactive = FALSE, verbose = FALSE,
                             show_preview = FALSE)
  expect_true(all(c("sig_id", "pert_type", "cell_iname") %in% names(out)))
})