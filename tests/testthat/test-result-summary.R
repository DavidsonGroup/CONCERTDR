# tests/testthat/test-result-summary.R
#
# Tests for annotate_drug_results(), extract_compound_id(), and
# fuzzy_drug_match().

# ── Fixtures ──────────────────────────────────────────────────────────────────

# A minimal results_df that mimics the output of process_signature_with_df
make_results_df <- function() {
  data.frame(
    compound   = c(
      "CVD001_K562_6H:BRD-K52492843-001-01-8:10:6",
      "CVD001_K562_24H:BRD-K52492843-001-01-8:10:24",
      "CVD001_HL60_24H:BRD-A00619745-001-03-9:1:24",
      "CVD001_HL60_6H:BRD-A00619745-001-03-9:1:6"
    ),
    Score      = c(-0.72, -0.65, 0.31, 0.42),
    pValue     = c(0.002, 0.005, 0.280, 0.350),
    pAdjValue  = c(0.020, 0.050, 0.560, 0.700),
    stringsAsFactors = FALSE
  )
}

make_siginfo_df <- function() {
  data.frame(
    sig_id     = c(
      "CVD001_K562_6H:BRD-K52492843-001-01-8:10:6",
      "CVD001_K562_24H:BRD-K52492843-001-01-8:10:24",
      "CVD001_HL60_24H:BRD-A00619745-001-03-9:1:24",
      "CVD001_HL60_6H:BRD-A00619745-001-03-9:1:6"
    ),
    pert_type  = rep("trt_cp", 4),
    pert_id    = c("BRD-K52492843-001-01-8", "BRD-K52492843-001-01-8",
                   "BRD-A00619745-001-03-9", "BRD-A00619745-001-03-9"),
    pert_iname = c("imatinib", "imatinib", "dasatinib", "dasatinib"),
    stringsAsFactors = FALSE
  )
}

make_compinfo_df <- function() {
  data.frame(
    pert_id  = c("BRD-K52492843-001-01-8", "BRD-A00619745-001-03-9"),
    pert_name = c("imatinib", "dasatinib"),
    cmap_name = c("imatinib", "dasatinib"),
    target    = c("BCR-ABL", "BCR-ABL;SRC"),
    moa       = c("kinase inhibitor", "kinase inhibitor"),
    phase     = c("Launched", "Launched"),
    stringsAsFactors = FALSE
  )
}

run_annotate <- function(...) {
  annotate_drug_results(
    results_df     = make_results_df(),
    sig_info_file  = make_siginfo_df(),
    comp_info_file = make_compinfo_df(),
    write_outputs  = FALSE,
    verbose        = FALSE,
    ...
  )
}

# ── annotate_drug_results: return structure ───────────────────────────────────

test_that("annotate_drug_results returns a list with four named tables", {
  views <- run_annotate()
  expect_type(views, "list")
  expect_named(views,
    c("wetlab_drug_view", "wetlab_gene_view", "tech_view_all",
      "drug_context_summary", "output_files"),
    ignore.order = TRUE
  )
})

test_that("tech_view_all has same row count as input results_df", {
  views <- run_annotate()
  expect_equal(nrow(views$tech_view_all), nrow(make_results_df()))
})

test_that("tech_view_all contains score, effect_direction and display_name", {
  views <- run_annotate()
  tv    <- views$tech_view_all
  expect_true("Score"            %in% names(tv))
  expect_true("effect_direction" %in% names(tv))
  expect_true("display_name"     %in% names(tv))
})

test_that("wetlab_drug_view is de-duplicated to one row per drug", {
  views <- run_annotate()
  # Two drugs (imatinib, dasatinib), each in multiple contexts
  expect_equal(nrow(views$wetlab_drug_view), 2)
})

test_that("wetlab_drug_view is ordered by score ascending", {
  views <- run_annotate()
  scores <- as.numeric(views$wetlab_drug_view$Score)
  expect_equal(scores, sort(scores))
})

test_that("wetlab_gene_view is empty when all perturbations are drugs", {
  views <- run_annotate()
  # All rows have pert_type = "trt_cp" => pert_kind = "Drug", no Gene rows
  expect_equal(nrow(views$wetlab_gene_view), 0)
})

test_that("drug_context_summary has one row per unique pert_id", {
  views <- run_annotate()
  dcs <- views$drug_context_summary
  expect_equal(nrow(dcs), 2)
  expect_true("n_contexts" %in% names(dcs))
  # Each drug has 2 contexts in our fixture
  expect_true(all(dcs$n_contexts == 2))
})

test_that("drug_context_summary has best_score column", {
  views <- run_annotate()
  expect_true("best_score" %in% names(views$drug_context_summary))
})

# ── effect_direction logic ────────────────────────────────────────────────────

test_that("negative scores produce 'Reversal' effect_direction", {
  views <- run_annotate()
  tv    <- views$tech_view_all
  neg_rows <- tv[as.numeric(tv$Score) < 0, ]
  expect_true(all(grepl("Reversal", neg_rows$effect_direction)))
})

test_that("positive scores produce 'Mimic' effect_direction", {
  views <- run_annotate()
  tv    <- views$tech_view_all
  pos_rows <- tv[as.numeric(tv$Score) > 0, ]
  expect_true(all(grepl("Mimic", pos_rows$effect_direction)))
})

# ── display_name fallback chain ───────────────────────────────────────────────

test_that("display_name uses cmap_name when available", {
  views <- run_annotate()
  tv    <- views$tech_view_all
  # cmap_name is "imatinib" / "dasatinib" in compinfo
  expect_true(all(tv$display_name %in% c("imatinib", "dasatinib")))
})

test_that("display_name falls back to pert_id when cmap_name and pert_name are NA", {
  comp_no_name <- make_compinfo_df()
  comp_no_name$cmap_name <- NA
  comp_no_name$pert_name <- NA

  views <- annotate_drug_results(
    results_df     = make_results_df(),
    sig_info_file  = make_siginfo_df(),
    comp_info_file = comp_no_name,
    write_outputs  = FALSE,
    verbose        = FALSE
  )
  tv <- views$tech_view_all
  expect_true(all(tv$display_name %in% comp_no_name$pert_id))
})

# ── keep_dose_in_drug_view ────────────────────────────────────────────────────

test_that("keep_dose_in_drug_view = TRUE adds dose_uM to wetlab_drug_view", {
  views <- run_annotate(keep_dose_in_drug_view = TRUE)
  expect_true("dose_uM" %in% names(views$wetlab_drug_view))
})

test_that("keep_dose_in_drug_view = FALSE (default) omits dose_uM", {
  views <- run_annotate(keep_dose_in_drug_view = FALSE)
  expect_false("dose_uM" %in% names(views$wetlab_drug_view))
})

# ── compound string parsing ───────────────────────────────────────────────────

test_that("dose_uM is correctly parsed from compound string", {
  views <- run_annotate()
  tv    <- views$tech_view_all
  # Compound strings end in :10:... or :1:...
  expect_true(all(tv$dose_uM %in% c(1, 10)))
})

test_that("time_h is correctly parsed from compound string", {
  views <- run_annotate()
  tv    <- views$tech_view_all
  expect_true(all(tv$time_h %in% c(6, 24)))
})

test_that("cell_line is extracted from compound string", {
  views <- run_annotate()
  tv    <- views$tech_view_all
  expect_true(all(tv$cell_line %in% c("K562", "HL60")))
})

# ── write_outputs = TRUE ──────────────────────────────────────────────────────

test_that("write_outputs = TRUE creates four TSV files", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  annotate_drug_results(
    results_df     = make_results_df(),
    sig_info_file  = make_siginfo_df(),
    comp_info_file = make_compinfo_df(),
    output_dir     = tmp,
    write_outputs  = TRUE,
    verbose        = FALSE
  )

  written <- list.files(tmp)
  expect_true("wetlab_drug_view.tsv"    %in% written)
  expect_true("wetlab_gene_view.tsv"    %in% written)
  expect_true("tech_view_all.tsv"       %in% written)
  expect_true("drug_context_summary.tsv" %in% written)
})

test_that("written TSV files are non-empty and readable", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  annotate_drug_results(
    results_df     = make_results_df(),
    sig_info_file  = make_siginfo_df(),
    comp_info_file = make_compinfo_df(),
    output_dir     = tmp,
    write_outputs  = TRUE,
    verbose        = FALSE
  )

  tv <- read.delim(file.path(tmp, "tech_view_all.tsv"), stringsAsFactors = FALSE)
  expect_gt(nrow(tv), 0)
})

# ── gene perturbation rows ────────────────────────────────────────────────────

test_that("wetlab_gene_view captures gene-knockdown rows", {
  # Add a trt_sh row to results and siginfo
  res_with_gene <- rbind(
    make_results_df(),
    data.frame(
      compound  = "GEN001_K562_96H:BRDA:0:96",
      Score     = -0.5,
      pValue    = 0.03,
      pAdjValue = 0.15,
      stringsAsFactors = FALSE
    )
  )
  sig_with_gene <- rbind(
    make_siginfo_df(),
    data.frame(
      sig_id     = "GEN001_K562_96H:BRDA:0:96",
      pert_type  = "trt_sh",
      pert_id    = "BRDA",
      pert_iname = "TP53",
      stringsAsFactors = FALSE
    )
  )

  views <- annotate_drug_results(
    results_df     = res_with_gene,
    sig_info_file  = sig_with_gene,
    comp_info_file = make_compinfo_df(),
    write_outputs  = FALSE,
    verbose        = FALSE
  )

  expect_gt(nrow(views$wetlab_gene_view), 0)
  expect_true(all(views$wetlab_gene_view$pert_kind == "Gene"))
})

# ── error handling ────────────────────────────────────────────────────────────

test_that("annotate_drug_results errors if compound_col is absent", {
  res <- make_results_df()
  names(res)[names(res) == "compound"] <- "sig_id"   # wrong column name

  expect_error(
    annotate_drug_results(
      results_df     = res,
      sig_info_file  = make_siginfo_df(),
      comp_info_file = make_compinfo_df(),
      compound_col   = "compound",   # still the old name
      write_outputs  = FALSE,
      verbose        = FALSE
    ),
    "compound"
  )
})

test_that("annotate_drug_results errors if score_col is absent", {
  res <- make_results_df()
  names(res)[names(res) == "Score"] <- "score_bad"

  expect_error(
    annotate_drug_results(
      results_df     = res,
      sig_info_file  = make_siginfo_df(),
      comp_info_file = make_compinfo_df(),
      write_outputs  = FALSE,
      verbose        = FALSE
    ),
    "Score"
  )
})

test_that("annotate_drug_results errors if sig_info_file path is missing", {
  expect_error(
    annotate_drug_results(
      results_df     = make_results_df(),
      sig_info_file  = "no_such_file.txt",
      comp_info_file = make_compinfo_df(),
      write_outputs  = FALSE,
      verbose        = FALSE
    ),
    "not found"
  )
})

# ── extract_compound_id ──────────────────────────────────────────────────────

test_that("extract_compound_id split_colon extracts the second field", {
  x   <- c("LIB_CELL_6H:BRD-K12345:10", "LIB_CELL_24H:BRD-K99999:1")
  out <- extract_compound_id(x, method = "split_colon", part_index = 2)
  expect_equal(out, c("BRD-K12345", "BRD-K99999"))
})

test_that("extract_compound_id split_colon extracts first field when part_index = 1", {
  x   <- "CVD001_K562_6H:BRD-K52492843:10"
  out <- extract_compound_id(x, method = "split_colon", part_index = 1)
  expect_equal(out, "CVD001_K562_6H")
})

test_that("extract_compound_id split_colon returns input if only one field", {
  x   <- "NO_COLON_HERE"
  out <- extract_compound_id(x, method = "split_colon", part_index = 2)
  expect_equal(out, "NO_COLON_HERE")
})

test_that("extract_compound_id split_underscore extracts the second token", {
  x   <- c("CVD001_K562_6H", "CVD002_HL60_24H")
  out <- extract_compound_id(x, method = "split_underscore", part_index = 2)
  expect_equal(out, c("K562", "HL60"))
})

test_that("extract_compound_id regex method extracts matching substring", {
  x   <- c("CVD001_K562_6H:BRD-K52492843-001-01-8:10")
  out <- extract_compound_id(x, method = "regex",
                              regex_pattern = "BRD-[A-Z0-9-]+")
  expect_equal(out, "BRD-K52492843-001-01-8")
})

test_that("extract_compound_id regex errors without pattern", {
  expect_error(
    extract_compound_id("A:B:C", method = "regex"),
    "regex_pattern must be provided"
  )
})

test_that("extract_compound_id errors on unknown method", {
  expect_error(
    extract_compound_id("A:B", method = "bad_method"),
    "Invalid method"
  )
})

test_that("extract_compound_id is vectorised", {
  x   <- paste0("A:BRD-K0000", 1:5, ":10")
  out <- extract_compound_id(x)
  expect_length(out, 5)
  expect_equal(out, paste0("BRD-K0000", 1:5))
})

# ── fuzzy_drug_match ─────────────────────────────────────────────────────────

test_that("fuzzy_drug_match requires RecordLinkage", {
  skip_if_not_installed("RecordLinkage")
  # If package is present the function should run, not throw
  expect_no_error(
    fuzzy_drug_match("aspirin", c("aspirin", "ibuprofen"), threshold = 90)
  )
})

test_that("fuzzy_drug_match returns data.frame with required columns", {
  skip_if_not_installed("RecordLinkage")
  out <- fuzzy_drug_match(c("imatinb"), c("imatinib", "dasatinib"),
                          method = "levenshtein", threshold = 70)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("query_name", "matched_name", "similarity_score") %in% names(out)))
})

test_that("fuzzy_drug_match returns empty df when no match exceeds threshold", {
  skip_if_not_installed("RecordLinkage")
  out <- fuzzy_drug_match("ZZZZZ", c("aspirin", "ibuprofen"), threshold = 99)
  expect_equal(nrow(out), 0)
})

test_that("fuzzy_drug_match returns top_n matches per query", {
  skip_if_not_installed("RecordLinkage")
  ref <- c("imatinib", "imatinol", "imatinine", "dasatinib")
  out <- fuzzy_drug_match("imatinib", ref, threshold = 50, top_n = 2)
  expect_lte(nrow(out), 2)
})

test_that("fuzzy_drug_match skips NA query entries", {
  skip_if_not_installed("RecordLinkage")
  out <- fuzzy_drug_match(c(NA, "aspirin"), c("aspirin"), threshold = 80)
  expect_equal(nrow(out), 1)
})

test_that("fuzzy_drug_match errors on invalid method", {
  skip_if_not_installed("RecordLinkage")
  expect_error(
    fuzzy_drug_match("aspirin", "aspirin", method = "bad_method"),
    "Invalid method"
  )
})

test_that("fuzzy_drug_match jaro method works", {
  skip_if_not_installed("RecordLinkage")
  out <- fuzzy_drug_match("aspirin", c("aspirin", "ibuprofen"),
                          method = "jaro", threshold = 80)
  expect_gte(nrow(out), 1)
  expect_true("aspirin" %in% out$matched_name)
})
