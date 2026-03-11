#' Internal Connectivity Scoring Methods
#'
#' @description
#' Self-contained implementations of connectivity scoring methods for matching
#' disease gene expression signatures to compound-induced gene expression
#' profiles. 
#'
#' All seven methods share a common permutation-based framework for computing
#' p-values and BH-adjusted p-values.
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

#' Convert expression matrix to sorted gene-name lists (for KS / GSEA-w0)
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @return List of character vectors ordered decreasingly by expression for each sample.
#' @keywords internal
.matrix_to_name_ranked_list <- function(refMatrix) {
  lapply(seq_len(ncol(refMatrix)), function(i) {
    names(refMatrix[order(refMatrix[, i], decreasing = TRUE), i])
  })
}

#' Convert expression matrix to sorted named-value lists (for GSEA-w1/w2)
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
#' @return Data frame with Score, pValue, pAdjValue columns.
#' @keywords internal
.permutation_pvalues <- function(score, permuteScoreMat, pAdjMethod = "BH",
                                  sample_names = NULL) {
  permuteScoreMat[is.na(permuteScoreMat)] <- 0
  pValue <- rowSums(abs(permuteScoreMat) >= abs(score)) / ncol(permuteScoreMat)
  pAdjust <- stats::p.adjust(pValue, method = pAdjMethod)
  out <- data.frame(Score = score, pValue = pValue, pAdjValue = pAdjust)
  if (!is.null(sample_names)) rownames(out) <- sample_names
  out
}

#' Apply a scoring function across reference lists, optionally in parallel
#' @param refList List of per-sample reference objects.
#' @param scoreFun Function used to score one reference object.
#' @param ... Additional arguments passed to \code{scoreFun}.
#' @return Named or unnamed numeric vector of scores, one per element of \code{refList}.
#' @keywords internal
.apply_score <- function(refList, scoreFun, ...) {
  vapply(refList, scoreFun, numeric(1), ...)
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
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_ks <- function(refMatrix, queryUp, queryDown,
                     permuteNum = 10000, pAdjMethod = "BH") {

  refMatrix <- .validate_ref_matrix(refMatrix)
  queryUp   <- as.character(queryUp)
  queryDown <- as.character(queryDown)

  ks_enrichment <- function(refList, query) {
    lenRef <- length(refList)
    queryRank <- match(query, refList)
    queryRank <- sort(queryRank[!is.na(queryRank)])
    lenQuery <- length(queryRank)

    if (lenQuery == 0) {
      return(0)
    }

    d <- seq_len(lenQuery) / lenQuery - queryRank / lenRef
    a <- max(d)
    b <- -min(d) + 1 / lenQuery
    ifelse(a > b, a, -b)
  }

  ks_combined <- function(refList, queryUp, queryDown) {
    scoreUp   <- ks_enrichment(refList, queryUp)
    scoreDown <- ks_enrichment(refList, queryDown)
    ifelse(scoreUp * scoreDown <= 0, scoreUp - scoreDown, 0)
  }

  refList  <- .matrix_to_name_ranked_list(refMatrix)
  queryUp  <- intersect(queryUp, rownames(refMatrix))
  queryDown <- intersect(queryDown, rownames(refMatrix))

  score <- .apply_score(refList, ks_combined,
                        queryUp = queryUp, queryDown = queryDown)

  permMat <- matrix(0, nrow = ncol(refMatrix), ncol = permuteNum)
  for (n in seq_len(permuteNum)) {
    bootUp   <- sample(rownames(refMatrix), length(queryUp))
    bootDown <- sample(rownames(refMatrix), length(queryDown))
    permMat[, n] <- .apply_score(refList, ks_combined,
                                 queryUp = bootUp, queryDown = bootDown)
  }

  .permutation_pvalues(score, permMat, pAdjMethod, colnames(refMatrix))
}

# ──────────────────────────────────────────────────────────────────────────────
# GSEA weight 0  (Subramanian et al. 2005)
# ──────────────────────────────────────────────────────────────────────────────

#' GSEA weight-0 connectivity score
#' @inheritParams score_ks
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_gsea0 <- function(refMatrix, queryUp, queryDown,
                        permuteNum = 10000, pAdjMethod = "BH") {

  refMatrix <- .validate_ref_matrix(refMatrix)
  queryUp   <- as.character(queryUp)
  queryDown <- as.character(queryDown)

  w0_enrichment <- function(refList, query) {
    tagIndicator   <- sign(match(refList, query, nomatch = 0))
    noTagIndicator <- 1 - tagIndicator
    N  <- length(refList)
    Nh <- length(query)
    Nm <- N - Nh
    correlVector   <- rep(1, N)
    sumCorrelTag   <- sum(correlVector[tagIndicator == 1])
    normTag   <- 1.0 / sumCorrelTag
    normNoTag <- 1.0 / Nm
    RES <- cumsum(tagIndicator * correlVector * normTag -
                    noTagIndicator * normNoTag)
    maxES <- max(RES);  minES <- min(RES)
    maxES <- ifelse(is.na(maxES), 0, maxES)
    minES <- ifelse(is.na(minES), 0, minES)
    ifelse(maxES > -minES, maxES, minES)
  }

  w0_combined <- function(refList, queryUp, queryDown) {
    scoreUp   <- w0_enrichment(refList, queryUp)
    scoreDown <- w0_enrichment(refList, queryDown)
    ifelse(scoreUp * scoreDown <= 0, scoreUp - scoreDown, 0)
  }

  refList   <- .matrix_to_name_ranked_list(refMatrix)
  queryUp   <- intersect(queryUp, rownames(refMatrix))
  queryDown <- intersect(queryDown, rownames(refMatrix))

  score <- .apply_score(refList, w0_combined,
                        queryUp = queryUp, queryDown = queryDown)

  permMat <- matrix(0, nrow = ncol(refMatrix), ncol = permuteNum)
  for (n in seq_len(permuteNum)) {
    bootUp   <- sample(rownames(refMatrix), length(queryUp))
    bootDown <- sample(rownames(refMatrix), length(queryDown))
    permMat[, n] <- .apply_score(refList, w0_combined,
                                 queryUp = bootUp, queryDown = bootDown)
  }

  .permutation_pvalues(score, permMat, pAdjMethod, colnames(refMatrix))
}

# ──────────────────────────────────────────────────────────────────────────────
# GSEA weight 1  (Subramanian et al. 2005)
# ──────────────────────────────────────────────────────────────────────────────

#' GSEA weight-1 connectivity score
#' @inheritParams score_ks
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_gsea1 <- function(refMatrix, queryUp, queryDown,
                        permuteNum = 10000, pAdjMethod = "BH") {

  refMatrix <- .validate_ref_matrix(refMatrix)
  queryUp   <- as.character(queryUp)
  queryDown <- as.character(queryDown)

  w1_enrichment <- function(refList, query) {
    tagIndicator   <- sign(match(names(refList), query, nomatch = 0))
    noTagIndicator <- 1 - tagIndicator
    N  <- length(refList)
    Nh <- length(query)
    Nm <- N - Nh
    correlVector   <- abs(refList)
    sumCorrelTag   <- sum(correlVector[tagIndicator == 1])
    normTag   <- 1.0 / sumCorrelTag
    normNoTag <- 1.0 / Nm
    RES <- cumsum(tagIndicator * correlVector * normTag -
                    noTagIndicator * normNoTag)
    maxES <- max(RES);  minES <- min(RES)
    maxES <- ifelse(is.na(maxES), 0, maxES)
    minES <- ifelse(is.na(minES), 0, minES)
    ifelse(maxES > -minES, maxES, minES)
  }

  w1_combined <- function(refList, queryUp, queryDown) {
    scoreUp   <- w1_enrichment(refList, queryUp)
    scoreDown <- w1_enrichment(refList, queryDown)
    ifelse(scoreUp * scoreDown <= 0, scoreUp - scoreDown, 0)
  }

  refList   <- .matrix_to_value_ranked_list(refMatrix)
  queryUp   <- intersect(queryUp, rownames(refMatrix))
  queryDown <- intersect(queryDown, rownames(refMatrix))

  score <- .apply_score(refList, w1_combined,
                        queryUp = queryUp, queryDown = queryDown)

  permMat <- matrix(0, nrow = ncol(refMatrix), ncol = permuteNum)
  for (n in seq_len(permuteNum)) {
    bootUp   <- sample(rownames(refMatrix), length(queryUp))
    bootDown <- sample(rownames(refMatrix), length(queryDown))
    permMat[, n] <- .apply_score(refList, w1_combined,
                                 queryUp = bootUp, queryDown = bootDown)
  }

  .permutation_pvalues(score, permMat, pAdjMethod, colnames(refMatrix))
}

# ──────────────────────────────────────────────────────────────────────────────
# GSEA weight 2  (Subramanian et al. 2005)
# ──────────────────────────────────────────────────────────────────────────────

#' GSEA weight-2 connectivity score
#' @inheritParams score_ks
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_gsea2 <- function(refMatrix, queryUp, queryDown,
                        permuteNum = 10000, pAdjMethod = "BH") {

  refMatrix <- .validate_ref_matrix(refMatrix)
  queryUp   <- as.character(queryUp)
  queryDown <- as.character(queryDown)

  w2_enrichment <- function(refList, query) {
    tagIndicator   <- sign(match(names(refList), query, nomatch = 0))
    noTagIndicator <- 1 - tagIndicator
    N  <- length(refList)
    Nh <- length(query)
    Nm <- N - Nh
    correlVector   <- abs(refList)^2
    sumCorrelTag   <- sum(correlVector[tagIndicator == 1])
    normTag   <- 1.0 / sumCorrelTag
    normNoTag <- 1.0 / Nm
    RES <- cumsum(tagIndicator * correlVector * normTag -
                    noTagIndicator * normNoTag)
    maxES <- max(RES);  minES <- min(RES)
    maxES <- ifelse(is.na(maxES), 0, maxES)
    minES <- ifelse(is.na(minES), 0, minES)
    ifelse(maxES > -minES, maxES, minES)
  }

  w2_combined <- function(refList, queryUp, queryDown) {
    scoreUp   <- w2_enrichment(refList, queryUp)
    scoreDown <- w2_enrichment(refList, queryDown)
    ifelse(scoreUp * scoreDown <= 0, scoreUp - scoreDown, 0)
  }

  refList   <- .matrix_to_value_ranked_list(refMatrix)
  queryUp   <- intersect(queryUp, rownames(refMatrix))
  queryDown <- intersect(queryDown, rownames(refMatrix))

  score <- .apply_score(refList, w2_combined,
                        queryUp = queryUp, queryDown = queryDown)

  permMat <- matrix(0, nrow = ncol(refMatrix), ncol = permuteNum)
  for (n in seq_len(permuteNum)) {
    bootUp   <- sample(rownames(refMatrix), length(queryUp))
    bootDown <- sample(rownames(refMatrix), length(queryDown))
    permMat[, n] <- .apply_score(refList, w2_combined,
                                 queryUp = bootUp, queryDown = bootDown)
  }

  .permutation_pvalues(score, permMat, pAdjMethod, colnames(refMatrix))
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
                       permuteNum = 10000, pAdjMethod = "BH") {

  refMatrix <- .validate_ref_matrix(refMatrix)
  if (!is.numeric(query)) stop("query must be a numeric vector")
  if (is.null(names(query))) stop("query must have names")
  if (topN > nrow(refMatrix) / 2) {
    stop("topN is larger than half the length of the gene list")
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

  refList <- .matrix_to_extreme_list(refMatrix, topN)

  score <- .apply_score(refList, xcos_single, query = query)

  permMat <- matrix(0, nrow = ncol(refMatrix), ncol = permuteNum)
  for (n in seq_len(permuteNum)) {
    perm_query <- query
    names(perm_query) <- sample(rownames(refMatrix), length(query))
    permMat[, n] <- .apply_score(refList, xcos_single, query = perm_query)
  }

  .permutation_pvalues(score, permMat, pAdjMethod, colnames(refMatrix))
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
                       permuteNum = 10000, pAdjMethod = "BH") {

  refMatrix <- .validate_ref_matrix(refMatrix)
  queryUp   <- as.character(queryUp)
  queryDown <- as.character(queryDown)
  if (topN > nrow(refMatrix) / 2) {
    stop("topN is larger than half the length of the gene list")
  }

  xsum_single <- function(refList, queryUp, queryDown) {
    scoreUp   <- sum(refList[match(queryUp,   names(refList))], na.rm = TRUE)
    scoreDown <- sum(refList[match(queryDown, names(refList))], na.rm = TRUE)
    scoreUp - scoreDown
  }

  refList   <- .matrix_to_extreme_list(refMatrix, topN)
  queryUp   <- intersect(queryUp, rownames(refMatrix))
  queryDown <- intersect(queryDown, rownames(refMatrix))

  score <- .apply_score(refList, xsum_single,
                        queryUp = queryUp, queryDown = queryDown)

  permMat <- matrix(0, nrow = ncol(refMatrix), ncol = permuteNum)
  for (n in seq_len(permuteNum)) {
    bootUp   <- sample(rownames(refMatrix), length(queryUp))
    bootDown <- sample(rownames(refMatrix), length(queryDown))
    permMat[, n] <- .apply_score(refList, xsum_single,
                                 queryUp = bootUp, queryDown = bootDown)
  }

  .permutation_pvalues(score, permMat, pAdjMethod, colnames(refMatrix))
}

# ──────────────────────────────────────────────────────────────────────────────
# Zhang Score  (Zhang & Gant 2008)
# ──────────────────────────────────────────────────────────────────────────────

#' Zhang connectivity score
#' @inheritParams score_ks
#' @return Data frame with Score, pValue, pAdjValue per sample.
#' @keywords internal
score_zhang <- function(refMatrix, queryUp, queryDown,
                        permuteNum = 10000, pAdjMethod = "BH") {

  refMatrix <- .validate_ref_matrix(refMatrix)
  if (is.null(queryUp))   queryUp   <- character(0)
  if (is.null(queryDown)) queryDown <- character(0)
  queryUp   <- as.character(queryUp)
  queryDown <- as.character(queryDown)

  zhang_single <- function(refRank, queryRank) {
    common <- intersect(names(refRank), names(queryRank))
    if (length(common) == 0) return(NA_real_)
    maxScore <- sum(abs(refRank)[seq_along(queryRank)] * abs(queryRank))
    if (maxScore == 0) return(0)
    sum(queryRank * refRank[names(queryRank)], na.rm = TRUE) / maxScore
  }

  refList <- .matrix_to_signed_rank_list(refMatrix)

  queryVector <- c(rep(1, length(queryUp)), rep(-1, length(queryDown)))
  names(queryVector) <- c(queryUp, queryDown)

  score <- .apply_score(refList, zhang_single, queryRank = queryVector)

  permMat <- matrix(0, nrow = ncol(refMatrix), ncol = permuteNum)
  for (n in seq_len(permuteNum)) {
    bootSample <- sample(c(-1, 1), replace = TRUE,
                         size = length(queryUp) + length(queryDown))
    names(bootSample) <- sample(rownames(refMatrix), replace = FALSE,
                                size = length(queryUp) + length(queryDown))
    permMat[, n] <- .apply_score(refList, zhang_single,
                                 queryRank = bootSample)
  }

  .permutation_pvalues(score, permMat, pAdjMethod, colnames(refMatrix))
}