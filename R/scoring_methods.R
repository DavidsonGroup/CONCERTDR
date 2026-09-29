#' Internal Connectivity Scoring Methods
#'
#' @description
#' Self-contained implementations of connectivity scoring methods for matching
#' disease gene expression signatures to compound-induced gene expression
#' profiles.
#'
#' The seven permutation-based methods share one framework
#' (\code{.permutation_test}) for computing p-values and adjusted p-values.
#' CamSum uses an analytic null instead; see \code{score_camsum}.
#' @return None; this is an internal documentation topic.
#'
#' @references
#' Lamb J et al. Science, 2006, 313(5795): 1929-1935 (KS method).
#' Subramanian A et al. PNAS, 2005, 102(43): 15545-15550 (GSEA methods).
#' Cheng J et al. Genome Medicine, 2014, 6(12): 95 (XCos, XSum methods).
#' Zhang S D et al. BMC Bioinformatics, 2008, 9(1): 258 (Zhang method).
#'
#' @name scoring_methods
#' @keywords internal
NULL

# ──────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ──────────────────────────────────────────────────────────────────────────────

#' Validate and coerce reference matrix
#' @param refMatrix Numeric matrix or data frame with genes as rows and samples as columns.
#' @return Numeric matrix with preserved row and column names.
#' @keywords internal
.validate_ref_matrix <- function(refMatrix) {
  if (is.data.frame(refMatrix)) refMatrix <- as.matrix(refMatrix)
  if (is.null(colnames(refMatrix)) || is.null(rownames(refMatrix))) {
    stop("refMatrix must have both rownames and colnames")
  }
  refMatrix
}

#' Convert expression matrix to sorted gene-name lists (for KS)
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @return List of character vectors ordered decreasingly by expression for each sample.
#' @keywords internal
.matrix_to_name_ranked_list <- function(refMatrix) {
  lapply(seq_len(ncol(refMatrix)), function(i) {
    names(refMatrix[order(refMatrix[, i], decreasing = TRUE), i])
  })
}

#' Convert expression matrix to sorted named-value lists (for GSEA)
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @return List of named numeric vectors ordered decreasingly by expression for each sample.
#' @keywords internal
.matrix_to_value_ranked_list <- function(refMatrix) {
  lapply(seq_len(ncol(refMatrix)), function(i) {
    refMatrix[order(refMatrix[, i], decreasing = TRUE), i]
  })
}

#' Convert expression matrix to top-N / bottom-N named-value lists (for XCos/XSum)
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @param topN Number of most positive and most negative genes to retain.
#' @return List of named numeric vectors containing the top and bottom genes per sample.
#' @keywords internal
.matrix_to_extreme_list <- function(refMatrix, topN) {
  lapply(seq_len(ncol(refMatrix)), function(i) {
    sorted <- refMatrix[order(refMatrix[, i], decreasing = TRUE), i]
    c(utils::head(sorted, topN), utils::tail(sorted, topN))
  })
}

#' Convert expression matrix to signed-rank lists (for Zhang)
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @return List of named numeric vectors containing signed ranks per sample.
#' @keywords internal
.matrix_to_signed_rank_list <- function(refMatrix) {
  lapply(seq_len(ncol(refMatrix)), function(i) {
    sorted <- refMatrix[order(abs(refMatrix[, i]), decreasing = TRUE), i]
    rank(abs(sorted)) * sign(sorted)
  })
}

#' Compute p-values and adjusted p-values via permutation
#'
#' @param score Numeric vector of observed scores (one per sample).
#' @param permuteScoreMat Matrix of permuted scores (nSamples x nPerms).
#' @param pAdjMethod Adjustment method passed to \code{stats::p.adjust}.
#' @param sample_names Optional sample names used as row names of the result.
#' @return Data frame with Score, pValue, pAdjValue columns.
#' @keywords internal
.permutation_pvalues <- function(score, permuteScoreMat, pAdjMethod = "BH",
                                 sample_names = NULL) {
  permuteScoreMat[is.na(permuteScoreMat)] <- 0
  pValue <- rowSums(abs(permuteScoreMat) >= abs(score)) / ncol(permuteScoreMat)
  out <- data.frame(Score = score, pValue = pValue,
                    pAdjValue = stats::p.adjust(pValue, method = pAdjMethod))
  if (!is.null(sample_names)) rownames(out) <- sample_names
  out
}

#' Apply a scoring function across reference lists
#' @param refList List of per-sample reference objects.
#' @param scoreFun Function used to score one reference object.
#' @param ... Additional arguments passed to \code{scoreFun}.
#' @return Named or unnamed numeric vector of scores, one per element of \code{refList}.
#' @keywords internal
.apply_score <- function(refList, scoreFun, ...) {
  vapply(refList, scoreFun, numeric(1), ...)
}

#' Score a query against every reference profile and against its permutations
#'
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @param refList Per-sample reference objects built from \code{refMatrix}.
#' @param scoreFun Function scoring one reference object; called as
#'   \code{scoreFun(ref, ...)} with the arguments in \code{queryArgs} or
#'   returned by \code{draw}.
#' @param queryArgs Named list of query arguments for the observed score.
#' @param draw Function of no arguments returning the named list of query
#'   arguments for one permutation.
#' @param permuteNum Number of permutations.
#' @param pAdjMethod P-value adjustment method.
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
.permutation_test <- function(refMatrix, refList, scoreFun, queryArgs, draw,
                              permuteNum, pAdjMethod) {
  score_with <- function(args) {
    do.call(.apply_score, c(list(refList, scoreFun), args))
  }
  score <- score_with(queryArgs)
  permMat <- matrix(0, nrow = ncol(refMatrix), ncol = permuteNum)
  for (n in seq_len(permuteNum)) permMat[, n] <- score_with(draw())
  .permutation_pvalues(score, permMat, pAdjMethod, colnames(refMatrix))
}

#' Permutation test for methods scoring an up and a down gene set
#'
#' Permuted queries draw random gene sets of the same sizes from the
#' reference genes.
#' @inheritParams .permutation_test
#' @param queryUp,queryDown Character vectors of gene symbols.
#' @param scoreFun Function \code{(ref, queryUp, queryDown)} scoring one
#'   reference object.
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
.permutation_test_up_down <- function(refMatrix, refList, scoreFun,
                                      queryUp, queryDown,
                                      permuteNum, pAdjMethod) {
  genes     <- rownames(refMatrix)
  queryUp   <- intersect(as.character(queryUp), genes)
  queryDown <- intersect(as.character(queryDown), genes)
  .permutation_test(
    refMatrix, refList, scoreFun,
    queryArgs = list(queryUp = queryUp, queryDown = queryDown),
    draw = function() {
      list(queryUp   = sample(genes, length(queryUp)),
           queryDown = sample(genes, length(queryDown)))
    },
    permuteNum, pAdjMethod)
}

#' Stop when \code{topN} exceeds half of the reference genes
#' @param topN Number of top/bottom genes per reference profile.
#' @param refMatrix Numeric matrix with genes as rows.
#' @return NULL, invisibly.
#' @keywords internal
.check_topN <- function(topN, refMatrix) {
  if (topN > nrow(refMatrix) / 2) {
    stop("topN is larger than half the length of the gene list")
  }
}

#' Combine up and down enrichment scores, zeroing them when they agree in sign
#' @param scoreUp,scoreDown Enrichment scores of the up and down gene sets.
#' @return Numeric.
#' @keywords internal
.combine_up_down <- function(scoreUp, scoreDown) {
  ifelse(scoreUp * scoreDown <= 0, scoreUp - scoreDown, 0)
}

# ──────────────────────────────────────────────────────────────────────────────
# KS Score  (Lamb et al. 2006)
# ──────────────────────────────────────────────────────────────────────────────

#' Kolmogorov-Smirnov connectivity score
#'
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @param queryUp Character vector of up-regulated gene symbols.
#' @param queryDown Character vector of down-regulated gene symbols.
#' @param permuteNum Number of permutations (default: 10000).
#' @param pAdjMethod P-value adjustment method (default: "BH").
#' @param vectorized Whether to use the vectorised engine (default: TRUE),
#'   which scores blocks of samples at once. It is skipped automatically for
#'   reference matrices with missing values or duplicated gene names. Scores
#'   are the same as the row-by-row code; p-values of KS and GSEA weight 0 use
#'   one permutation null shared by all samples (same distribution, different
#'   random numbers), and the other methods reproduce the row-by-row p-values
#'   for the same seed. Use \code{FALSE} for the row-by-row code.
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_ks <- function(refMatrix, queryUp, queryDown,
                     permuteNum = 10000, pAdjMethod = "BH", vectorized = TRUE) {
  refMatrix <- .validate_ref_matrix(refMatrix)
  if (.use_fast(refMatrix, vectorized)) {
    return(.fast_ks(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod))
  }

  ks_enrichment <- function(refList, query) {
    lenRef <- length(refList)
    queryRank <- match(query, refList)
    queryRank <- sort(queryRank[!is.na(queryRank)])
    lenQuery <- length(queryRank)
    if (lenQuery == 0) return(0)

    d <- seq_len(lenQuery) / lenQuery - queryRank / lenRef
    a <- max(d)
    b <- -min(d) + 1 / lenQuery
    ifelse(a > b, a, -b)
  }

  ks_combined <- function(refList, queryUp, queryDown) {
    .combine_up_down(ks_enrichment(refList, queryUp),
                     ks_enrichment(refList, queryDown))
  }

  .permutation_test_up_down(refMatrix, .matrix_to_name_ranked_list(refMatrix),
                            ks_combined, queryUp, queryDown,
                            permuteNum, pAdjMethod)
}

# ──────────────────────────────────────────────────────────────────────────────
# GSEA weights 0, 1 and 2  (Subramanian et al. 2005)
# ──────────────────────────────────────────────────────────────────────────────

#' Weighted GSEA connectivity score
#'
#' Genes are ranked by decreasing value; each gene in the query contributes
#' \eqn{|v|^{weight}} to the running enrichment sum (weight 0 gives the
#' classical unweighted KS-like statistic).
#' @inheritParams score_ks
#' @param weight Exponent applied to the absolute reference values.
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
.score_gsea <- function(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod,
                        weight, vectorized = TRUE) {
  refMatrix <- .validate_ref_matrix(refMatrix)
  if (.use_fast(refMatrix, vectorized)) {
    fast <- if (weight == 0) .fast_gsea0 else
      function(...) .fast_gseaW(..., weight = weight)
    return(fast(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod))
  }

  gsea_enrichment <- function(refList, query) {
    tagIndicator   <- sign(match(names(refList), query, nomatch = 0))
    noTagIndicator <- 1 - tagIndicator
    Nm <- length(refList) - length(query)
    correlVector <- abs(refList)^weight
    normTag   <- 1.0 / sum(correlVector[tagIndicator == 1])
    normNoTag <- 1.0 / Nm
    RES <- cumsum(tagIndicator * correlVector * normTag -
                    noTagIndicator * normNoTag)
    maxES <- max(RES)
    minES <- min(RES)
    maxES <- ifelse(is.na(maxES), 0, maxES)
    minES <- ifelse(is.na(minES), 0, minES)
    ifelse(maxES > -minES, maxES, minES)
  }

  gsea_combined <- function(refList, queryUp, queryDown) {
    .combine_up_down(gsea_enrichment(refList, queryUp),
                     gsea_enrichment(refList, queryDown))
  }

  .permutation_test_up_down(refMatrix, .matrix_to_value_ranked_list(refMatrix),
                            gsea_combined, queryUp, queryDown,
                            permuteNum, pAdjMethod)
}

#' GSEA weight-0 connectivity score
#' @inheritParams score_ks
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_gsea0 <- function(refMatrix, queryUp, queryDown,
                        permuteNum = 10000, pAdjMethod = "BH", vectorized = TRUE) {
  .score_gsea(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod, weight = 0,
              vectorized = vectorized)
}

#' GSEA weight-1 connectivity score
#' @inheritParams score_ks
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_gsea1 <- function(refMatrix, queryUp, queryDown,
                        permuteNum = 10000, pAdjMethod = "BH", vectorized = TRUE) {
  .score_gsea(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod, weight = 1,
              vectorized = vectorized)
}

#' GSEA weight-2 connectivity score
#' @inheritParams score_ks
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_gsea2 <- function(refMatrix, queryUp, queryDown,
                        permuteNum = 10000, pAdjMethod = "BH", vectorized = TRUE) {
  .score_gsea(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod, weight = 2,
              vectorized = vectorized)
}

# ──────────────────────────────────────────────────────────────────────────────
# XCos Score  (Cheng et al. 2014)
# ──────────────────────────────────────────────────────────────────────────────

#' Extreme cosine similarity score
#'
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @param query Named numeric vector (gene symbols → fold-change or rank).
#' @param topN Number of top/bottom genes per reference profile (default: 500).
#' @param permuteNum Number of permutations (default: 10000).
#' @param pAdjMethod P-value adjustment method (default: "BH").
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_xcos <- function(refMatrix, query, topN = 500,
                       permuteNum = 10000, pAdjMethod = "BH", vectorized = TRUE) {
  refMatrix <- .validate_ref_matrix(refMatrix)
  if (!is.numeric(query)) stop("query must be a numeric vector")
  if (is.null(names(query))) stop("query must have names")
  .check_topN(topN, refMatrix)
  if (.use_fast(refMatrix, vectorized) && all(names(query) %in% rownames(refMatrix)) &&
      !anyDuplicated(names(query)) && length(query) > 0) {
    return(.fast_xcos(refMatrix, query, topN, permuteNum, pAdjMethod))
  }

  xcos_single <- function(refList, query) {
    common <- intersect(names(refList), names(query))
    if (length(common) == 0) return(NA_real_)
    r <- refList[common]
    q <- query[common]
    denom <- sqrt(crossprod(r) * crossprod(q))
    if (denom == 0) return(0)
    (crossprod(r, q) / denom)[1, 1]
  }

  .permutation_test(
    refMatrix, .matrix_to_extreme_list(refMatrix, topN), xcos_single,
    queryArgs = list(query = query),
    draw = function() {
      perm_query <- query
      names(perm_query) <- sample(rownames(refMatrix), length(query))
      list(query = perm_query)
    },
    permuteNum, pAdjMethod)
}

# ──────────────────────────────────────────────────────────────────────────────
# XSum Score  (Cheng et al. 2014)
# ──────────────────────────────────────────────────────────────────────────────

#' Extreme sum score
#'
#' @inheritParams score_ks
#' @param topN Number of top/bottom genes per reference profile (default: 500).
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_xsum <- function(refMatrix, queryUp, queryDown, topN = 500,
                       permuteNum = 10000, pAdjMethod = "BH", vectorized = TRUE) {
  refMatrix <- .validate_ref_matrix(refMatrix)
  .check_topN(topN, refMatrix)
  if (.use_fast(refMatrix, vectorized)) {
    return(.fast_xsum(refMatrix, queryUp, queryDown, topN, permuteNum, pAdjMethod))
  }

  xsum_single <- function(refList, queryUp, queryDown) {
    scoreUp   <- sum(refList[match(queryUp,   names(refList))], na.rm = TRUE)
    scoreDown <- sum(refList[match(queryDown, names(refList))], na.rm = TRUE)
    scoreUp - scoreDown
  }

  .permutation_test_up_down(refMatrix, .matrix_to_extreme_list(refMatrix, topN),
                            xsum_single, queryUp, queryDown,
                            permuteNum, pAdjMethod)
}

# ──────────────────────────────────────────────────────────────────────────────
# Zhang Score  (Zhang & Gant 2008)
# ──────────────────────────────────────────────────────────────────────────────

#' Zhang connectivity score
#' @inheritParams score_ks
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_zhang <- function(refMatrix, queryUp, queryDown,
                        permuteNum = 10000, pAdjMethod = "BH", vectorized = TRUE) {
  refMatrix <- .validate_ref_matrix(refMatrix)
  queryUp   <- as.character(queryUp)
  queryDown <- as.character(queryDown)
  nQuery    <- length(queryUp) + length(queryDown)
  if (.use_fast(refMatrix, vectorized) && nQuery > 0 && nQuery <= nrow(refMatrix) &&
      all(c(queryUp, queryDown) %in% rownames(refMatrix)) &&
      !anyDuplicated(c(queryUp, queryDown))) {
    return(.fast_zhang(refMatrix, queryUp, queryDown, permuteNum, pAdjMethod))
  }

  zhang_single <- function(refRank, queryRank) {
    common <- intersect(names(refRank), names(queryRank))
    if (length(common) == 0) return(NA_real_)
    maxScore <- sum(abs(refRank)[seq_along(queryRank)] * abs(queryRank))
    if (maxScore == 0) return(0)
    sum(queryRank * refRank[names(queryRank)], na.rm = TRUE) / maxScore
  }

  queryVector <- c(rep(1, length(queryUp)), rep(-1, length(queryDown)))
  names(queryVector) <- c(queryUp, queryDown)

  .permutation_test(
    refMatrix, .matrix_to_signed_rank_list(refMatrix), zhang_single,
    queryArgs = list(queryRank = queryVector),
    draw = function() {
      bootSample <- sample(c(-1, 1), replace = TRUE, size = nQuery)
      names(bootSample) <- sample(rownames(refMatrix), replace = FALSE,
                                  size = nQuery)
      list(queryRank = bootSample)
    },
    permuteNum, pAdjMethod)
}
