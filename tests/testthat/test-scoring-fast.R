# Vectorised scoring engine (R/scoring_fast.R) against the row-by-row code

fast_fixture <- function(seed = 1, digits = 1, n_genes = 400, n_samples = 60) {
  set.seed(seed)
  M <- round(matrix(rnorm(n_genes * n_samples), n_genes, n_samples,
                    dimnames = list(paste0("g", seq_len(n_genes)),
                                    paste0("s", seq_len(n_samples)))), digits)
  up <- sample(rownames(M), 20)
  dn <- sample(setdiff(rownames(M), up), 15)
  list(M = M, up = up, dn = dn,
       query = stats::setNames(c(rnorm(20, 1, .5), rnorm(15, -1, .5)), c(up, dn)))
}

with_seed <- function(seed, expr) {
  set.seed(seed)
  expr
}

B <- 100

test_that("methods with per-sample nulls reproduce the row-by-row results", {
  f <- fast_fixture()
  cases <- list(
    gsea1 = function(v) score_gsea1(f$M, f$up, f$dn, permuteNum = B, vectorized = v),
    gsea2 = function(v) score_gsea2(f$M, f$up, f$dn, permuteNum = B, vectorized = v),
    zhang = function(v) score_zhang(f$M, f$up, f$dn, permuteNum = B, vectorized = v),
    xsum  = function(v) score_xsum(f$M, f$up, f$dn, topN = 50, permuteNum = B, vectorized = v),
    xcos  = function(v) score_xcos(f$M, f$query, topN = 50, permuteNum = B, vectorized = v)
  )
  for (name in names(cases)) {
    fast <- with_seed(7, cases[[name]](TRUE))
    slow <- with_seed(7, cases[[name]](FALSE))
    expect_equal(rownames(fast), rownames(slow), info = name)
    expect_equal(fast$Score, slow$Score, tolerance = 1e-8, info = name)
    # Same permuted gene sets for the same seed; a p-value may move by one
    # permutation where a permuted score equals the observed one up to rounding
    expect_lte(max(abs(fast$pValue - slow$pValue)), 1 / B + 1e-12)
  }
})

test_that("KS scores match and its shared null gives the same p-value distribution", {
  f <- fast_fixture(digits = 2, n_samples = 80)
  fast <- with_seed(3, score_ks(f$M, f$up, f$dn, permuteNum = 2000))
  slow <- with_seed(3, score_ks(f$M, f$up, f$dn, permuteNum = 2000, vectorized = FALSE))
  expect_equal(fast$Score, slow$Score, tolerance = 1e-12)
  expect_lt(mean(abs(fast$pValue - slow$pValue)), 0.03)
  expect_gt(stats::cor(fast$pValue, slow$pValue), 0.99)
})

test_that("GSEA weight 0 scores match except where the maximum equals the minimum", {
  f <- fast_fixture(seed = 2, digits = 3, n_genes = 600, n_samples = 100)
  fast <- with_seed(3, score_gsea0(f$M, f$up, f$dn, permuteNum = 50))
  slow <- with_seed(3, score_gsea0(f$M, f$up, f$dn, permuteNum = 50, vectorized = FALSE))
  expect_lte(mean(abs(fast$Score - slow$Score) > 1e-8), 0.02)
})

test_that("an exact tie between the maximum and the absolute minimum is resolved to the minimum", {
  # One hit at position 3 of 5: the running sum peaks at +2/4 and dips to -2/4
  expect_equal(.fast_gsea0_es(matrix(3, 1, 1), 1, 5), -0.5)
})

test_that("references with missing values use the row-by-row code", {
  f <- fast_fixture()
  f$M[3, 4] <- NA
  for (fun in list(score_ks, score_gsea1, score_zhang)) {
    a <- with_seed(5, suppressWarnings(fun(f$M, f$up, f$dn, permuteNum = 20)))
    b <- with_seed(5, suppressWarnings(fun(f$M, f$up, f$dn, permuteNum = 20, vectorized = FALSE)))
    expect_identical(a, b)
  }
})

test_that("block size does not change the results", {
  f <- fast_fixture()
  run <- function() {
    list(gsea2 = with_seed(4, score_gsea2(f$M, f$up, f$dn, permuteNum = 30)),
         xsum = with_seed(4, score_xsum(f$M, f$up, f$dn, topN = 50, permuteNum = 30)),
         zhang = with_seed(4, score_zhang(f$M, f$up, f$dn, permuteNum = 30)))
  }
  whole <- run()
  old <- options(CONCERTDR.chunk_size = 7L)
  on.exit(options(old))
  expect_equal(run(), whole)
})

test_that("process_signature_with_df can switch the vectorised engine off", {
  f <- fast_fixture(n_genes = 300)
  sig <- data.frame(Gene = c(f$up, f$dn), log2FC = c(rep(1, 20), rep(-1, 15)))
  ref <- data.frame(gene_symbol = rownames(f$M), f$M, check.names = FALSE)
  methods <- c("gsea1", "gsea2", "zhang", "xsum", "xcos")
  fast <- with_seed(1, suppressMessages(process_signature_with_df(
    sig, ref, methods = methods, permutations = 20, topN = 20)))
  slow <- with_seed(1, suppressMessages(process_signature_with_df(
    sig, ref, methods = methods, permutations = 20, topN = 20, vectorized = FALSE)))
  expect_named(fast$results, methods)
  for (m in methods) expect_equal(fast$results[[m]]$Score, slow$results[[m]]$Score,
                                  tolerance = 1e-8)
})
