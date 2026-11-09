# ──────────────────────────────────────────────────────────────────────────────
# Vectorised scoring engine
#
# The scoring methods in scoring_methods.R score one reference profile at a
# time from R. The functions below compute exactly the same statistics for a
# block of profiles at once with matrixStats, so the number of R-level calls no
# longer grows with the number of profiles times the number of permutations.
#
# Formulas are unchanged. Three ideas carry the speed:
#   * profiles are processed in column blocks; one \code{colRanks()} per block
#     gives every profile's descending order (ties broken by original row, as
#     in \code{order(x, decreasing = TRUE)});
#   * the permuted gene sets are drawn up front in the same order as the
#     row-by-row code, so a given seed yields the same permuted sets;
#   * for KS and GSEA weight 0 the score depends on the gene set only through
#     its ranks, whose null distribution is identical for every profile, so one
#     shared null replaces one null per profile.
# ──────────────────────────────────────────────────────────────────────────────

#' Number of profiles scored at a time by the vectorised engine
#' @return Integer; option \code{CONCERTDR.chunk_size} (default: 2000).
#' @keywords internal
.fast_chunk_size <- function() {
  as.integer(getOption("CONCERTDR.chunk_size", 2000L))
}

#' Whether the vectorised engine can score this reference matrix
#'
#' Matrices with missing values or duplicated gene names take the row-by-row
#' path, which defines the behaviour for those cases.
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @param vectorized Whether the user asked for the vectorised engine.
#' @return Logical scalar.
#' @keywords internal
.use_fast <- function(refMatrix, vectorized) {
  isTRUE(vectorized) && is.numeric(refMatrix) && !anyNA(refMatrix) &&
    !anyDuplicated(rownames(refMatrix))
}

#' Turn permutation exceedance counts into a result data frame
#' @param score Observed scores.
#' @param count Number of permutations with \code{|permuted| >= |observed|}.
#' @param permuteNum Number of permutations.
#' @param pAdjMethod P-value adjustment method.
#' @param sample_names Row names of the result.
#' @return Data frame with Score, pValue, pAdjValue columns.
#' @keywords internal
.count_pvalues <- function(score, count, permuteNum, pAdjMethod, sample_names) {
  pValue <- count / permuteNum
  out <- data.frame(Score = score, pValue = pValue,
                    pAdjValue = stats::p.adjust(pValue, method = pAdjMethod))
  rownames(out) <- sample_names
  out
}

#' Positions (1 = highest value) of all genes in every profile of a block
#' @param M Numeric matrix, genes x profiles.
#' @return Integer matrix of the same shape.
#' @keywords internal
.fast_positions <- function(M) {
  matrixStats::colRanks(-M, ties.method = "first", preserveShape = TRUE)
}

#' Positions of a gene set sorted within each profile
#' @param pos Position matrix from \code{.fast_positions}.
#' @param S Integer row indices of the gene set.
#' @param values Optional matrix (genes x profiles) whose rows \code{S} are
#'   reordered the same way.
#' @return List with \code{SP} (k x profiles sorted positions) and, when
#'   \code{values} is given, \code{W} (the values in the same order).
#' @keywords internal
.fast_sorted_set <- function(pos, S, values = NULL) {
  P <- pos[S, , drop = FALSE]
  k <- nrow(P)
  n <- ncol(P)
  p <- as.vector(P)
  o <- order(p + rep(seq_len(n) - 1, each = k) * nrow(pos), method = "radix")
  out <- list(SP = matrix(p[o], k, n))
  if (!is.null(values)) out$W <- matrix(as.vector(values[S, , drop = FALSE])[o], k, n)
  out
}

#' KS enrichment of sorted positions
#' @param SP k x profiles matrix of sorted positions.
#' @param k Size of the gene set.
#' @param G Number of genes.
#' @return Numeric vector, one value per profile.
#' @keywords internal
.fast_ks_es <- function(SP, k, G) {
  if (k == 0L) return(numeric(ncol(SP)))
  D <- seq_len(k) / k - SP / G
  a <- matrixStats::colMaxs(D)
  b <- -matrixStats::colMins(D) + 1 / k
  ifelse(a > b, a, -b)
}

#' Unweighted GSEA enrichment of sorted positions
#'
#' The running sum ends at 0, so its maximum sits at a hit (or is 0) and its
#' minimum just before a hit (or is 0). The comparison is done on integers
#' (the sums scaled by \eqn{k (G - k)}), so a maximum that equals the absolute
#' minimum is resolved the same way for every profile.
#' @inheritParams .fast_ks_es
#' @return Numeric vector, one value per profile.
#' @keywords internal
.fast_gsea0_es <- function(SP, k, G) {
  if (k == 0L) return(numeric(ncol(SP)))
  k <- as.numeric(k)
  Nm <- G - k
  i <- seq_len(k)
  peak   <- i * Nm - (SP - i) * k
  trough <- (i - 1) * Nm - (SP - i) * k
  maxES <- pmax(matrixStats::colMaxs(peak), 0)
  minES <- pmin(matrixStats::colMins(trough), 0)
  ifelse(maxES > -minES, maxES, minES) / (k * Nm)
}

#' Weighted GSEA enrichment of sorted positions
#' @inheritParams .fast_ks_es
#' @param W k x profiles matrix of the weights, sorted like \code{SP}.
#' @return Numeric vector, one value per profile.
#' @keywords internal
.fast_gseaW_es <- function(SP, W, k, G) {
  if (k == 0L) return(numeric(ncol(SP)))
  Nm <- G - k
  i <- seq_len(k)
  cw <- matrixStats::colCumsums(W)
  sw <- cw[k, ]
  offset <- (SP - i) / Nm
  peaks   <- cw / rep(sw, each = k) - offset
  cw0 <- rbind(0, cw[-k, , drop = FALSE])
  troughs <- cw0 / rep(sw, each = k) - offset
  maxES <- pmax(matrixStats::colMaxs(peaks), 0)
  minES <- pmin(matrixStats::colMins(troughs), 0)
  es <- ifelse(maxES > -minES, maxES, minES)
  es[!is.finite(sw) | sw == 0] <- 0
  es
}

#' Score every profile of a reference matrix, block by block
#'
#' @param refMatrix Numeric matrix, genes x profiles.
#' @param make_ctx Function of a block matrix returning the per-block context.
#' @param score Function \code{(ctx, draw)} returning one score per profile of
#'   the block; \code{draw} is \code{NULL} for the observed query and one
#'   element of \code{draws} for a permutation.
#' @param draws List of permuted queries.
#' @param permuteNum Number of permutations.
#' @param pAdjMethod P-value adjustment method.
#' @return Data frame with Score, pValue, pAdjValue per profile.
#' @keywords internal
.fast_permutation_test <- function(refMatrix, make_ctx, score, draws,
                                   permuteNum, pAdjMethod) {
  n <- ncol(refMatrix)
  observed <- numeric(n)
  count <- numeric(n)
  for (cols in .camsum_blocks(n, .fast_chunk_size())) {
    ctx <- make_ctx(refMatrix[, cols, drop = FALSE])
    obs <- score(ctx, NULL)
    observed[cols] <- obs
    hits <- numeric(length(cols))
    for (draw in draws) {
      perm <- score(ctx, draw)
      perm[is.na(perm)] <- 0
      hits <- hits + (abs(perm) >= abs(obs))
    }
    count[cols] <- hits
  }
  .count_pvalues(observed, count, permuteNum, pAdjMethod, colnames(refMatrix))
}

#' Draw the permuted queries in the order the row-by-row code draws them
#' @param permuteNum Number of permutations.
#' @param draw_one Function of no arguments drawing one permuted query.
#' @return List of \code{permuteNum} draws.
#' @keywords internal
.fast_draws <- function(permuteNum, draw_one) {
  lapply(seq_len(permuteNum), function(b) draw_one())
}

#' Vectorised test for methods scoring an up and a down gene set
#'
#' @param refMatrix Numeric matrix, genes x profiles.
#' @param queryUp,queryDown Character vectors of gene symbols.
#' @param make_ctx,score See \code{.fast_permutation_test}; \code{score} is
#'   called as \code{score(ctx, list(U =, D =))} where \code{U} and \code{D}
#'   are row indices, or with \code{NULL} for the observed query.
#' @param permuteNum Number of permutations.
#' @param pAdjMethod P-value adjustment method.
#' @return Data frame with Score, pValue, pAdjValue per profile.
#' @keywords internal
.fast_up_down_test <- function(refMatrix, queryUp, queryDown, make_ctx, score,
                               permuteNum, pAdjMethod) {
  genes <- rownames(refMatrix)
  G <- length(genes)
  U <- match(intersect(as.character(queryUp), genes), genes)
  D <- match(intersect(as.character(queryDown), genes), genes)
  draws <- .fast_draws(permuteNum, function() {
    list(U = sample.int(G, length(U)), D = sample.int(G, length(D)))
  })
  .fast_permutation_test(
    refMatrix, make_ctx,
    function(ctx, draw) {
      if (is.null(draw)) draw <- list(U = U, D = D)
      score(ctx, draw)
    },
    draws, permuteNum, pAdjMethod)
}

#' KS or unweighted GSEA with one null shared by all profiles
#'
#' The score depends on the gene set only through its positions in a profile,
#' and a random gene set has the same distribution of positions in every
#' profile, so the permutation null does not depend on the profile.
#' @param refMatrix Numeric matrix, genes x profiles.
#' @param queryUp,queryDown Character vectors of gene symbols.
#' @param es Enrichment function of sorted positions (\code{.fast_ks_es} or
#'   \code{.fast_gsea0_es}).
#' @param permuteNum Number of permutations.
#' @param pAdjMethod P-value adjustment method.
#' @return Data frame with Score, pValue, pAdjValue per profile.
#' @keywords internal
.fast_shared_null_test <- function(refMatrix, queryUp, queryDown, es,
                                   permuteNum, pAdjMethod) {
  genes <- rownames(refMatrix)
  G <- length(genes)
  U <- match(intersect(as.character(queryUp), genes), genes)
  D <- match(intersect(as.character(queryDown), genes), genes)
  kU <- length(U)
  kD <- length(D)
  combined <- function(SPu, SPd) .combine_up_down(es(SPu, kU, G), es(SPd, kD, G))

  n <- ncol(refMatrix)
  observed <- numeric(n)
  for (cols in .camsum_blocks(n, .fast_chunk_size())) {
    pos <- .fast_positions(refMatrix[, cols, drop = FALSE])
    observed[cols] <- combined(.fast_sorted_set(pos, U)$SP,
                               .fast_sorted_set(pos, D)$SP)
  }

  # Same B draws for every profile; positions are sorted within each draw
  spu <- matrix(0, kU, permuteNum)
  spd <- matrix(0, kD, permuteNum)
  for (b in seq_len(permuteNum)) {
    spu[, b] <- sort.int(sample.int(G, kU))
    spd[, b] <- sort.int(sample.int(G, kD))
  }
  null <- combined(spu, spd)
  null[is.na(null)] <- 0
  sorted_null <- sort(abs(null))
  count <- permuteNum - findInterval(abs(observed), sorted_null, left.open = TRUE)
  .count_pvalues(observed, count, permuteNum, pAdjMethod, colnames(refMatrix))
}

# ── Per-method entry points ──────────────────────────────────────────────────

#' Vectorised KS score
#' @inheritParams .fast_shared_null_test
#' @keywords internal
.fast_ks <- function(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod) {
  .fast_shared_null_test(refMatrix, queryUp, queryDown, .fast_ks_es,
                         permuteNum, pAdjMethod)
}

#' Vectorised GSEA weight-0 score
#' @inheritParams .fast_shared_null_test
#' @keywords internal
.fast_gsea0 <- function(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod) {
  .fast_shared_null_test(refMatrix, queryUp, queryDown, .fast_gsea0_es,
                         permuteNum, pAdjMethod)
}

#' Vectorised weighted GSEA score
#' @inheritParams .fast_shared_null_test
#' @param weight Exponent applied to the absolute reference values.
#' @keywords internal
.fast_gseaW <- function(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod,
                        weight) {
  G <- nrow(refMatrix)
  .fast_up_down_test(
    refMatrix, queryUp, queryDown,
    make_ctx = function(M) list(pos = .fast_positions(M), W = abs(M)^weight),
    score = function(ctx, draw) {
      es <- function(S) {
        s <- .fast_sorted_set(ctx$pos, S, ctx$W)
        .fast_gseaW_es(s$SP, s$W, length(S), G)
      }
      .combine_up_down(es(draw$U), es(draw$D))
    },
    permuteNum, pAdjMethod)
}

#' Vectorised XSum score
#' @inheritParams .fast_shared_null_test
#' @param topN Number of top/bottom genes kept per profile.
#' @keywords internal
.fast_xsum <- function(refMatrix, queryUp, queryDown, topN, permuteNum, pAdjMethod) {
  G <- nrow(refMatrix)
  .fast_up_down_test(
    refMatrix, queryUp, queryDown,
    make_ctx = function(M) {
      pos <- .fast_positions(M)
      M[pos > topN & pos <= G - topN] <- 0
      list(ext = M)
    },
    score = function(ctx, draw) {
      matrixStats::colSums2(ctx$ext, rows = draw$U) -
        matrixStats::colSums2(ctx$ext, rows = draw$D)
    },
    permuteNum, pAdjMethod)
}

#' Vectorised XCos score
#' @inheritParams .fast_shared_null_test
#' @param query Named numeric vector.
#' @param topN Number of top/bottom genes kept per profile.
#' @keywords internal
.fast_xcos <- function(refMatrix, query, topN, permuteNum, pAdjMethod) {
  genes <- rownames(refMatrix)
  G <- length(genes)
  Sq <- match(names(query), genes)
  qv <- unname(query)
  draws <- .fast_draws(permuteNum, function() sample.int(G, length(qv)))
  .fast_permutation_test(
    refMatrix,
    make_ctx = function(M) {
      pos <- .fast_positions(M)
      keep <- pos <= topN | pos > G - topN
      M[!keep] <- 0
      list(ext = M, keep = keep)
    },
    score = function(ctx, draw) {
      S <- if (is.null(draw)) Sq else draw
      E <- ctx$ext[S, , drop = FALSE]
      member <- ctx$keep[S, , drop = FALSE]
      # Members are the genes in the extreme list, including those whose
      # value is exactly 0
      sc <- matrixStats::colSums2(E * qv) /
        sqrt(matrixStats::colSums2(E * E) * matrixStats::colSums2(member * qv^2))
      sc[!is.finite(sc)] <- 0
      sc[matrixStats::colSums2(member) == 0] <- NA_real_
      sc
    },
    draws, permuteNum, pAdjMethod)
}

#' Vectorised Zhang score
#' @inheritParams .fast_shared_null_test
#' @keywords internal
.fast_zhang <- function(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod) {
  genes <- rownames(refMatrix)
  G <- length(genes)
  S <- match(c(as.character(queryUp), as.character(queryDown)), genes)
  n <- length(S)
  signs <- c(rep(1, length(queryUp)), rep(-1, length(queryDown)))
  draws <- .fast_draws(permuteNum, function() {
    list(signs = sample(c(-1, 1), replace = TRUE, size = n),
         genes = sample.int(G, n))
  })
  .fast_permutation_test(
    refMatrix,
    make_ctx = function(M) {
      ranks <- matrixStats::colRanks(abs(M), ties.method = "average",
                                     preserveShape = TRUE)
      # Sum of the n largest ranks of each profile (ties share their mean rank)
      nth <- matrixStats::colOrderStats(ranks, which = G - n + 1L)
      above <- ranks > rep(nth, each = G)
      max_score <- matrixStats::colSums2(ranks * above) +
        (n - matrixStats::colSums2(above)) * nth
      list(signed_rank = ranks * sign(M), max_score = max_score)
    },
    score = function(ctx, draw) {
      s <- if (is.null(draw)) list(signs = signs, genes = S) else draw
      sc <- matrixStats::colSums2(ctx$signed_rank[s$genes, , drop = FALSE] * s$signs) /
        ctx$max_score
      sc[ctx$max_score == 0] <- 0
      sc
    },
    draws, permuteNum, pAdjMethod)
}
