# Tests for the CamSum scoring method

make_ref <- function(G = 300, n = 80, seed = 1) {
  set.seed(seed)
  f <- rnorm(n)   # shared factor -> correlated genes
  load <- c(rep(0.6, 30), rep(0, G - 30))
  X <- matrix(rnorm(G * n), G, n) + outer(load, f)
  X[, 5] <- X[, 5] * 3                            # one high-amplitude profile
  dimnames(X) <- list(paste0("g", seq_len(G)), paste0("s", seq_len(n)))
  X
}
up   <- paste0("g", 1:20)
down <- paste0("g", 21:35)

brute_rho <- function(X, up, down) {
  s <- apply(X, 2, stats::sd)
  Y <- sweep(X[c(up, down), ], 2, s, "/")
  Y[down, ] <- -Y[down, ]
  R <- stats::cor(t(Y))
  mean(R[upper.tri(R)])
}

# Score with rho_bar estimated on X itself (warns)
score_local <- function(X, ...) {
  suppressWarnings(
    CONCERTDR:::score_camsum(X, ..., pert_type = rep(NA, ncol(X))))
}

# Write a minimal GCTX + siginfo + geneinfo into dir; returns file paths
write_fake_cmap <- function(dir, X, pert_type) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  gctx <- file.path(dir, "fake.gctx")
  rhdf5::h5createFile(gctx)
  groups <- c("0", "0/DATA", "0/DATA/0", "0/META", "0/META/ROW", "0/META/COL")
  for (g in groups) {
    rhdf5::h5createGroup(gctx, g)
  }
  gene_id <- as.character(seq_len(nrow(X)) + 1000L)
  rhdf5::h5write(unname(X), gctx, "0/DATA/0/matrix")
  rhdf5::h5write(gene_id, gctx, "0/META/ROW/id")
  rhdf5::h5write(colnames(X), gctx, "0/META/COL/id")
  rhdf5::h5closeAll()
  geneinfo <- file.path(dir, "geneinfo.txt")
  utils::write.table(data.frame(gene_id = gene_id, gene_symbol = rownames(X)),
                     geneinfo, sep = "\t", quote = FALSE, row.names = FALSE)
  siginfo <- file.path(dir, "siginfo.txt")
  utils::write.table(data.frame(sig_id = colnames(X), pert_type = pert_type),
                     siginfo, sep = "\t", quote = FALSE, row.names = FALSE)
  list(gctx = gctx, geneinfo = geneinfo, siginfo = siginfo)
}

test_that("streaming rho_bar equals the brute-force k x k correlation", {
  X <- make_ref()
  r <- compute_camsum_rho(X, up, down, chunk = 7L)
  expect_equal(as.numeric(r), brute_rho(X, up, down), tolerance = 1e-12)
  expect_equal(attr(r, "k"), 35L)
  expect_equal(attr(r, "kU"), 20L)
  expect_equal(attr(r, "kD"), 15L)
})

test_that("score follows T = (sum_U z - sum_D z) / (s_i sqrt(k VIF))", {
  X <- make_ref()
  res <- score_local(X, up, down)
  rho <- brute_rho(X, up, down)
  vif <- max(1, 1 + (35 - 1) * rho)
  s <- apply(X, 2, stats::sd)
  expect_equal(res$Score,
               unname((colSums(X[up, ]) - colSums(X[down, ])) /
                        (s * sqrt(35 * vif))),
               tolerance = 1e-12)
  expect_equal(unname(attr(res, "VIF")), vif)
  expect_equal(unname(attr(res, "rho_source")), "refMatrix")
  expect_equal(rownames(res), colnames(X))
  expect_equal(colnames(res), c("Score", "pValue", "pAdjValue"))
  # scale invariance per profile: multiplying a column does not change its score
  X2 <- X; X2[, 3] <- X2[, 3] * 10
  expect_equal(score_local(X2, up, down)$Score[3], res$Score[3],
               tolerance = 1e-10)
})

test_that("estimating rho_bar on the supplied reference warns", {
  X <- make_ref()
  expect_warning(CONCERTDR:::score_camsum(X, up, down, pert_type = rep(NA, 80)),
                 "estimated on the supplied reference")
  expect_warning(CONCERTDR:::score_camsum(X, up, down, pert_type = rep(NA, 80)),
                 "Fewer than 500 profiles")
})

test_that("alternatives are consistent and default is two-sided", {
  X <- make_ref()
  two  <- score_local(X, up, down)
  gr   <- score_local(X, up, down, alternative = "greater")
  less <- score_local(X, up, down, alternative = "less")
  expect_equal(attr(two, "alternative"), "two.sided")
  expect_equal(two$pValue, 2 * pmin(gr$pValue, less$pValue), tolerance = 1e-12)
  expect_equal(gr$pValue + less$pValue, rep(1, ncol(X)), tolerance = 1e-12)
  expect_equal(two$pAdjValue, stats::p.adjust(two$pValue, "BH"))
})

test_that("VIF is floored at 1 and a supplied rho_bar is used as given", {
  X <- make_ref()
  neg <- CONCERTDR:::score_camsum(X, up, down, rho_bar = -0.2)
  expect_equal(unname(attr(neg, "VIF")), 1)
  expect_equal(unname(attr(neg, "rho_source")), "supplied")
  pos <- CONCERTDR:::score_camsum(X, up, down, rho_bar = 0.05)
  expect_equal(unname(attr(pos, "VIF")), 1 + 34 * 0.05)
  expect_error(CONCERTDR:::score_camsum(X, up, down, rho_bar = c(0.1, 0.2)),
               "single number")
})

test_that("genes in both directions are dropped; up-only queries work", {
  X <- make_ref()
  a <- score_local(X, c(up, "g21"), down)
  expect_equal(attr(a, "kU"), 20L)
  expect_equal(attr(a, "kD"), 14L)
  b <- score_local(X, up, character(0))
  expect_equal(attr(b, "kD"), 0L)
  expect_true(all(is.finite(b$Score)))
  expect_error(score_local(X, "nope", "none"), "no genes in common")
})

test_that("constant profiles get NA with a warning", {
  X <- make_ref(); X[, 2] <- 0
  expect_warning(
    res <- CONCERTDR:::score_camsum(X, up, down, rho_bar = 0.01),
    "zero or non-finite standard deviation")
  expect_true(is.na(res$Score[2]))
  expect_true(all(is.finite(res$Score[-2])))
})

test_that("the session cache accepts keys longer than 10,000 bytes", {
  env <- new.env()
  key <- strrep("GENE\t", 5000)
  expect_equal(CONCERTDR:::.camsum_cached(key, 1, env), 1)
  expect_equal(CONCERTDR:::.camsum_cached(key, stop("not evaluated"), env), 1)
})

test_that("rho_bar comes from the full library in the GCTX for any subset", {
  skip_if_not_installed("rhdf5")
  dir <- tempfile("camsum_gctx"); dir.create(dir)
  A <- make_ref(n = 120, seed = 2)                       # library "trt_oe"
  B <- make_ref(n = 90, seed = 3)                        # library "trt_xpr"
  B[1:30, ] <- B[1:30, ] + outer(rep(0.8, 30), rnorm(90))
  colnames(B) <- paste0("x", seq_len(ncol(B)))
  X <- cbind(A, B)
  f <- write_fake_cmap(dir, X, rep(c("trt_oe", "trt_xpr"), c(ncol(A), ncol(B))))
  old <- options(CONCERTDR.gctx_file = f$gctx,
                 CONCERTDR.siginfo_file = f$siginfo,
                 CONCERTDR.geneinfo_file = f$geneinfo)
  on.exit({ options(old); unlink(dir, recursive = TRUE) }, add = TRUE)

  # a subset of trt_oe gets the full-library rho_bar, looked up via siginfo
  expect_message(res <- CONCERTDR:::score_camsum(A[, 1:15], up, down,
                                                 chunk = 25L),
                 "all 120 trt_oe profiles")
  expect_equal(unname(attr(res, "rho_bar")), brute_rho(A, up, down),
               tolerance = 1e-12)
  expect_equal(unname(attr(res, "rho_source")), "library_gctx")
  expect_equal(unname(attr(res, "rho_n_profiles")), ncol(A))
  # same query again: cached, no second pass over the GCTX
  expect_silent(CONCERTDR:::score_camsum(A[, 20:30], up, down, chunk = 25L))

  # mixed types: each type gets its own library rho_bar and VIF
  mix <- cbind(A[, 1:10], B[, 1:10])
  res2 <- suppressMessages(CONCERTDR:::score_camsum(mix, up, down))
  expect_equal(attr(res2, "rho_bar")[["trt_oe"]], brute_rho(A, up, down),
               tolerance = 1e-12)
  expect_equal(attr(res2, "rho_bar")[["trt_xpr"]], brute_rho(B, up, down),
               tolerance = 1e-12)
  s <- apply(mix, 2, stats::sd)
  vif <- attr(res2, "VIF")[rep(c("trt_oe", "trt_xpr"), each = 10)]
  expect_equal(res2$Score,
               unname((colSums(mix[up, ]) - colSums(mix[down, ])) /
                        (s * sqrt(35 * vif))), tolerance = 1e-12)

  # a smaller gene space: s_i is taken over the reference's own genes
  sub_genes <- paste0("g", seq(1, 300, by = 3))
  lm_up <- intersect(up, sub_genes); lm_down <- intersect(down, sub_genes)
  res3 <- suppressMessages(
    CONCERTDR:::score_camsum(A[sub_genes, 1:15], lm_up, lm_down))
  expect_equal(unname(attr(res3, "rho_bar")),
               brute_rho(A[sub_genes, ], lm_up, lm_down), tolerance = 1e-12)

  # a reference gene missing from the GCTX -> falls back with a warning
  odd <- rbind(A[, 1:15], zzz = rnorm(15))
  expect_warning(CONCERTDR:::score_camsum(odd, up, down),
                 "estimated on the supplied reference")

  # a GCTX that lacks part of the library -> falls back with a warning
  f2 <- write_fake_cmap(file.path(dir, "part"), X[, -(1:5)],
                        rep(c("trt_oe", "trt_xpr"), c(ncol(A) - 5, ncol(B))))
  options(CONCERTDR.gctx_file = f2$gctx)
  expect_warning(CONCERTDR:::score_camsum(A[, 6:20], up, down),
                 "does not hold every profile")
  options(CONCERTDR.gctx_file = f$gctx)

  # process_signature_with_df uses the same route
  sig <- data.frame(Gene = c(up, down),
                    log2FC = c(rep(1, length(up)), rep(-1, length(down))))
  ref_df <- data.frame(gene_symbol = rownames(A), A[, 1:15],
                       check.names = FALSE)
  out <- suppressMessages(process_signature_with_df(
    signature_file = sig, reference_df = ref_df, methods = "camsum"))
  expect_equal(unname(out$settings$camsum$rho_bar), brute_rho(A, up, down),
               tolerance = 1e-12)
  expect_equal(unname(out$settings$camsum$rho_source), "library_gctx")
})

test_that("process_signature_with_df runs camsum and keeps its settings", {
  sig_file <- system.file("extdata", "example_signature.txt",
                          package = "CONCERTDR")
  ref_file <- system.file("extdata", "example_reference_df.csv",
                          package = "CONCERTDR")
  ref_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
  ref_df$gene_symbol <- rownames(ref_df)
  old <- options(CONCERTDR.gctx_file = NULL, CONCERTDR.siginfo_file = NULL,
                 CONCERTDR.geneinfo_file = NULL)
  on.exit(options(old), add = TRUE)
  res <- suppressWarnings(suppressMessages(process_signature_with_df(
    signature_file = sig_file, reference_df = ref_df,
    methods = "camsum", camsum_alternative = "less", save_files = FALSE)))
  cs <- res$results$camsum
  expect_true(all(c("compound", "Score", "pValue", "pAdjValue", "rank") %in%
                  colnames(cs)))
  expect_equal(nrow(cs), ncol(ref_df) - 1L)
  expect_equal(res$settings$camsum$alternative, "less")
  expect_true(is.numeric(res$settings$camsum$rho_bar))
  expect_true(all(res$settings$camsum$VIF >= 1))
  expect_equal(unname(res$settings$camsum$rho_source), "refMatrix")
  # default methods are unchanged and do not include camsum
  default_methods <- eval(formals(process_signature_with_df)$methods)
  expect_false("camsum" %in% default_methods)
})
