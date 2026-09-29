# ──────────────────────────────────────────────────────────────────────────────
# CamSum  (correlation-adjusted sum score with an analytic null)
# ──────────────────────────────────────────────────────────────────────────────

#' Match query genes to reference rows for CamSum
#'
#' Genes listed in both directions are dropped from both.
#' @param genes Character vector of reference row names.
#' @param queryUp,queryDown Character vectors of gene symbols.
#' @return List with integer row indices \code{U} and \code{D}.
#' @keywords internal
.camsum_rows <- function(genes, queryUp, queryDown) {
  queryUp   <- unique(as.character(queryUp))
  queryDown <- unique(as.character(queryDown))
  both <- intersect(queryUp, queryDown)
  list(U = which(genes %in% setdiff(queryUp, both)),
       D = which(genes %in% setdiff(queryDown, both)))
}

#' Split 1..n into consecutive blocks of at most \code{chunk}
#' @param n Number of items.
#' @param chunk Block size.
#' @return List of integer vectors (empty when \code{n = 0}).
#' @keywords internal
.camsum_blocks <- function(n, chunk) {
  split(seq_len(n), ceiling(seq_len(n) / chunk))
}

#' Per-profile standard deviation over all genes, computed in column blocks
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @param chunk Number of columns processed at a time.
#' @return Numeric vector, one standard deviation per column.
#' @keywords internal
.camsum_col_sd <- function(refMatrix, chunk = 2000L) {
  G <- nrow(refMatrix)
  s <- numeric(ncol(refMatrix))
  for (cols in .camsum_blocks(ncol(refMatrix), chunk)) {
    X <- refMatrix[, cols, drop = FALSE]
    X <- X - rep(colMeans(X), each = G)
    s[cols] <- sqrt(colSums(X * X) / (G - 1L))
  }
  s
}

#' Scale query rows by s_i, negate down rows, drop unusable profiles
#' @param Z k x m block of query rows (up rows first, then down rows).
#' @param s Standard deviations of the m profiles.
#' @param kU Number of up rows.
#' @return k x m' matrix of \eqn{y = \pm z / s_i}.
#' @keywords internal
.camsum_signed <- function(Z, s, kU) {
  ok <- is.finite(s) & s > 0 & colSums(!is.finite(Z)) == 0L
  Y <- Z[, ok, drop = FALSE] / rep(s[ok], each = nrow(Z))
  down <- seq_len(nrow(Z)) > kU
  Y[down, ] <- -Y[down, , drop = FALSE]
  Y
}

#' Two-pass streaming estimate of rho_bar
#'
#' \code{read_block(b)} returns block \code{b} of \code{.camsum_signed()}
#' output. Uses \eqn{1'R1 = \mathrm{Var}_n(w'(y_n - \bar y))} with
#' \eqn{w_g = 1/\mathrm{sd}_g}, so the k x k matrix is never formed.
#' @param read_block Function of a block number returning a k x m matrix.
#' @param n_blocks Number of blocks.
#' @param k Number of query genes.
#' @return List with \code{rho_bar} and \code{n} (profiles used).
#' @keywords internal
.camsum_rho_stream <- function(read_block, n_blocks, k) {
  s1 <- numeric(k); s2 <- numeric(k); n <- 0L
  for (b in seq_len(n_blocks)) {
    Y <- read_block(b)
    s1 <- s1 + rowSums(Y); s2 <- s2 + rowSums(Y * Y); n <- n + ncol(Y)
  }
  if (k < 2L || n < 2L) return(list(rho_bar = NA_real_, n = n))
  mu  <- s1 / n
  sdv <- sqrt(pmax((s2 - n * mu^2) / (n - 1), 0))
  good <- is.finite(sdv) & sdv > 0
  k_good <- sum(good)
  if (k_good < 2L) return(list(rho_bar = NA_real_, n = n))
  w <- ifelse(good, 1 / sdv, 0)
  t1 <- 0; t2 <- 0
  for (b in seq_len(n_blocks)) {
    tt <- as.numeric(crossprod(w, read_block(b))) - sum(w * mu)
    t1 <- t1 + sum(tt); t2 <- t2 + sum(tt * tt)
  }
  one_R_one <- (t2 - t1^2 / n) / (n - 1)
  list(rho_bar = (one_R_one - k_good) / (k_good * (k_good - 1)), n = n)
}

#' Mean inter-gene correlation of a query across a reference library
#'
#' Computes \eqn{\bar\rho}, the mean off-diagonal Pearson correlation of the
#' query genes measured across the profiles of \code{refMatrix}, after dividing
#' each profile by its standard deviation over all genes and negating the
#' down-regulated rows. The \eqn{k \times k} correlation matrix is never formed:
#' two passes over column blocks use the identity
#' \eqn{1'R1 = \mathrm{Var}_n(w'(y_n - \bar y))} with
#' \eqn{w_g = 1/\mathrm{sd}_g}.
#'
#' CamSum's p-values were calibrated with \eqn{\bar\rho} estimated on the whole
#' treatment library of one perturbation type (for example all \code{trt_oe}
#' or all \code{trt_xpr} profiles). \code{\link{process_signature_with_df}}
#' does this automatically from the GCTX file when the CMap file options are
#' set; this function estimates \eqn{\bar\rho} on the matrix it is given.
#'
#' @param refMatrix Numeric matrix (or data frame) with genes as rows and
#'   profiles as columns, with row and column names.
#' @param queryUp,queryDown Character vectors of up- and down-regulated gene
#'   symbols. Genes present in both are dropped from both.
#' @param chunk Number of columns processed at a time (default: 2000).
#'
#' @return A single numeric value, \eqn{\bar\rho}, with attributes \code{k},
#'   \code{kU}, \code{kD} (query genes found in \code{refMatrix}) and
#'   \code{n_profiles} (profiles used). \code{NA} when fewer than two query
#'   genes vary across the library.
#'
#' @examples
#' set.seed(1)
#' ref <- matrix(rnorm(200 * 50), 200, 50,
#'               dimnames = list(paste0("g", 1:200), paste0("s", 1:50)))
#' compute_camsum_rho(ref, queryUp = paste0("g", 1:10),
#'                    queryDown = paste0("g", 11:20))
#'
#' @export
compute_camsum_rho <- function(refMatrix, queryUp, queryDown, chunk = 2000L) {
  refMatrix <- .validate_ref_matrix(refMatrix)
  idx <- .camsum_rows(rownames(refMatrix), queryUp, queryDown)
  s_i <- .camsum_col_sd(refMatrix, chunk)
  .camsum_rho(refMatrix, idx$U, idx$D, s_i, chunk)
}

#' Streaming rho_bar on an in-memory matrix
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @param U,D Integer row indices of up- and down-regulated query genes.
#' @param s_i Per-profile standard deviations.
#' @param chunk Number of columns processed at a time.
#' @return Numeric scalar with attributes k, kU, kD, n_profiles.
#' @keywords internal
.camsum_rho <- function(refMatrix, U, D, s_i, chunk = 2000L) {
  blocks <- .camsum_blocks(ncol(refMatrix), chunk)
  read_block <- function(b) {
    cols <- blocks[[b]]
    .camsum_signed(refMatrix[c(U, D), cols, drop = FALSE], s_i[cols],
                   length(U))
  }
  k <- length(U) + length(D)
  r <- .camsum_rho_stream(read_block, length(blocks), k)
  structure(r$rho_bar, k = k, kU = length(U), kD = length(D),
            n_profiles = r$n)
}

# ── Full-library rho_bar from the GCTX file ─────────────────────────────────

.camsum_cache <- new.env(parent = emptyenv())

#' Evaluate \code{value} once per key (by default, once per session)
#'
#' Values are kept in a list inside \code{env} because keys can be longer
#' than the 10,000-byte limit on variable names.
#' @param key Character cache key.
#' @param value Expression, evaluated only on a cache miss; must not be NULL.
#' @param env Environment used as the cache.
#' @return The cached value.
#' @keywords internal
.camsum_cached <- function(key, value, env = .camsum_cache) {
  if (is.null(env$store[[key]])) env$store[[key]] <- value
  env$store[[key]]
}

#' Cache key for a file: path plus modification time
#' @param path File path.
#' @return Character string.
#' @keywords internal
.camsum_file_key <- function(path) {
  paste(normalizePath(path), file.mtime(path))
}

#' Read columns of a tab-separated CMap metadata file once per session
#' @param path File path.
#' @param select Columns to read.
#' @return data.frame.
#' @keywords internal
.camsum_read_meta_file <- function(path, select) {
  .camsum_cached(
    paste("meta", .camsum_file_key(path), paste(select, collapse = ",")),
    data.table::fread(path, select = select, quote = "", data.table = FALSE,
                      showProgress = FALSE))
}

#' Row and column ids of a GCTX file, read once per session
#' @param gctx_file GCTX file path.
#' @return List with \code{row} and \code{col} character vectors.
#' @keywords internal
.camsum_gctx_ids <- function(gctx_file) {
  .camsum_cached(paste("ids", .camsum_file_key(gctx_file)), list(
    row = trimws(as.character(fast_read_meta(
      gctx_file, c("/0/META/ROW/id", "/META/ROW/id")))),
    col = trimws(as.character(fast_read_meta(
      gctx_file, c("/0/META/COL/id", "/META/COL/id"))))))
}

#' Read whole columns of the GCTX data matrix
#' @param gctx_file GCTX file path.
#' @param cols Integer column indices.
#' @return Numeric matrix, all genes x \code{cols}.
#' @keywords internal
.camsum_gctx_columns <- function(gctx_file, cols) {
  X <- tryCatch(
    rhdf5::h5read(gctx_file, "/0/DATA/0/matrix", index = list(NULL, cols)),
    error = function(e) {
      rhdf5::h5read(gctx_file, "/DATA/0/matrix", index = list(NULL, cols))
    })
  storage.mode(X) <- "double"
  X
}

#' Path from a CONCERTDR file option, or NULL when unset or missing
#' @param name Option name.
#' @return Character path or NULL.
#' @keywords internal
.camsum_option_file <- function(name) {
  f <- getOption(name, NULL)
  if (is.character(f) && length(f) == 1L && file.exists(f)) f else NULL
}

#' Look up the perturbation type of reference profiles
#'
#' Uses \code{metadata} (the \code{"metadata"} attribute written by
#' \code{\link{extract_cmap_data_from_siginfo}}) when it covers the profiles,
#' otherwise the siginfo file in option \code{CONCERTDR.siginfo_file}.
#' @param sig_ids Character vector of profile identifiers.
#' @param metadata Optional data frame with \code{sample_id} and
#'   \code{pert_type} columns.
#' @return Character vector of perturbation types, \code{NA} where unknown.
#' @keywords internal
.camsum_pert_types <- function(sig_ids, metadata = NULL) {
  out <- rep(NA_character_, length(sig_ids))
  if (is.data.frame(metadata) &&
      all(c("sample_id", "pert_type") %in% names(metadata))) {
    out <- as.character(metadata$pert_type[match(sig_ids, metadata$sample_id)])
  }
  siginfo_file <- .camsum_option_file("CONCERTDR.siginfo_file")
  if (anyNA(out) && !is.null(siginfo_file)) {
    si <- .camsum_read_meta_file(siginfo_file, c("sig_id", "pert_type"))
    miss <- is.na(out)
    out[miss] <- as.character(si$pert_type[match(sig_ids[miss], si$sig_id)])
  }
  out
}

#' rho_bar over every profile of one perturbation type in the GCTX file
#'
#' Profiles are all siginfo rows of \code{pert_type}; every one of them must
#' be in the GCTX file. Each profile's \eqn{s_i} is taken over the GCTX rows
#' whose gene symbols are the row names of the reference, so the gene space
#' matches the reference. The result is cached for the session.
#' @param pert_type Perturbation type.
#' @param genes Row names of the reference matrix.
#' @param U,D Integer row indices (into \code{genes}) of the query genes.
#' @param chunk Number of profiles read at a time.
#' @return Numeric scalar with attributes k, kU, kD, n_profiles; \code{NULL}
#'   when the CMap file options are not set, a reference gene is not in the
#'   GCTX file, or the GCTX file lacks part of the library.
#' @keywords internal
.camsum_rho_gctx <- function(pert_type, genes, U, D, chunk = 2000L) {
  gctx_file     <- .camsum_option_file("CONCERTDR.gctx_file")
  siginfo_file  <- .camsum_option_file("CONCERTDR.siginfo_file")
  geneinfo_file <- .camsum_option_file("CONCERTDR.geneinfo_file")
  if (is.null(gctx_file) || is.null(siginfo_file) || is.null(geneinfo_file)) {
    return(NULL)
  }

  ids <- .camsum_gctx_ids(gctx_file)
  gi <- .camsum_read_meta_file(geneinfo_file, c("gene_id", "gene_symbol"))
  row_symbols <- gi$gene_symbol[match(ids$row, as.character(gi$gene_id))]
  rows <- match(genes, row_symbols)
  si <- .camsum_read_meta_file(siginfo_file, c("sig_id", "pert_type"))
  cols <- match(si$sig_id[si$pert_type == pert_type], ids$col)
  if (anyNA(rows) || length(cols) < 2L || anyNA(cols)) return(NULL)
  cols <- sort(cols)

  key <- paste("rho", .camsum_file_key(gctx_file), pert_type,
               paste(sort(genes), collapse = "\t"),
               paste(sort(genes[U]), collapse = "\t"),
               paste(sort(genes[D]), collapse = "\t"), sep = "\n")
  .camsum_cached(key, {
    message("CamSum: estimating rho_bar on all ", length(cols), " ",
            pert_type, " profiles in ", basename(gctx_file),
            " (two passes over the library)...")
    blocks <- .camsum_blocks(length(cols), chunk)
    s_cache <- new.env(parent = emptyenv())    # s_i per block, reused in pass 2
    read_block <- function(b) {
      X <- .camsum_gctx_columns(gctx_file, cols[blocks[[b]]])
      X <- X[rows, , drop = FALSE]
      s <- .camsum_cached(as.character(b), .camsum_col_sd(X, chunk), s_cache)
      .camsum_signed(X[c(U, D), , drop = FALSE], s, length(U))
    }
    r <- .camsum_rho_stream(read_block, length(blocks), length(U) + length(D))
    rhdf5::h5closeAll()
    structure(r$rho_bar, k = length(U) + length(D), kU = length(U),
              kD = length(D), n_profiles = r$n)
  })
}


# ── Scoring ──────────────────────────────────────────────────────────────────

#' CamSum connectivity score
#'
#' Correlation-adjusted sum of reference z-scores over the query genes, with an
#' analytic standard-normal p-value (no permutations). For profile \eqn{i}
#' with standard deviation \eqn{s_i} over all genes of \code{refMatrix} and a
#' query of \eqn{k = k_U + k_D} genes found in \code{refMatrix},
#' \deqn{T_i = \frac{\sum_{g \in U} z_{gi} - \sum_{g \in D} z_{gi}}
#'                  {s_i \sqrt{k \cdot \mathrm{VIF}}}, \qquad
#'       \mathrm{VIF} = \max\{1,\ 1 + (k - 1)\bar\rho\}.}
#' The floor at 1 matches \code{limma::camera} with
#' \code{allow.neg.cor = FALSE}. Positive scores mean the profile mimics the
#' query, negative scores mean it reverses it.
#'
#' \eqn{\bar\rho} is taken, in order of preference, from \code{rho_bar}; from
#' every profile of the same perturbation type in the GCTX file named by
#' option \code{CONCERTDR.gctx_file}, which must hold the whole library of
#' that type (options \code{CONCERTDR.siginfo_file} and
#' \code{CONCERTDR.geneinfo_file} are needed too; profiles of different types
#' get their own \eqn{\bar\rho}); or, with a warning, the columns of
#' \code{refMatrix} (see \code{\link{compute_camsum_rho}}). Only the p-values
#' depend on this choice; the ranking of profiles within one perturbation
#' type does not.
#'
#' Calibration was checked on LINCS \code{trt_oe} and \code{trt_xpr} with
#' library-level \eqn{\bar\rho}; on \code{trt_cp} the score distribution is
#' left-skewed and p-values are anti-conservative, so use the score for
#' ranking only there.
#'
#' @param refMatrix Numeric matrix with genes as rows and samples as columns.
#' @param queryUp,queryDown Character vectors of gene symbols. Genes present in
#'   both are dropped from both.
#' @param alternative \code{"two.sided"} (default, same convention as the
#'   permutation methods), \code{"greater"} (mimic) or \code{"less"}
#'   (reversal).
#' @param rho_bar Optional precomputed \eqn{\bar\rho}, used for all profiles.
#' @param pAdjMethod P-value adjustment method (default: "BH").
#' @param chunk Number of columns processed at a time (default: 2000).
#' @param pert_type Optional perturbation type of each column; looked up from
#'   the siginfo file in option \code{CONCERTDR.siginfo_file} when
#'   \code{NULL}.
#' @return Data frame with Score, pValue, pAdjValue per sample. Attributes
#'   \code{k}, \code{kU}, \code{kD} and \code{alternative}, plus
#'   \code{rho_bar}, \code{VIF}, \code{rho_source} and \code{rho_n_profiles}
#'   as vectors named by perturbation type.
#' @keywords internal
score_camsum <- function(refMatrix, queryUp, queryDown,
                         alternative = c("two.sided", "greater", "less"),
                         rho_bar = NULL, pAdjMethod = "BH", chunk = 2000L,
                         pert_type = NULL) {
  alternative <- match.arg(alternative)
  refMatrix <- .validate_ref_matrix(refMatrix)
  if (!is.null(rho_bar) && (!is.numeric(rho_bar) || length(rho_bar) != 1L)) {
    stop("rho_bar must be a single number")
  }
  genes <- rownames(refMatrix)
  idx <- .camsum_rows(genes, queryUp, queryDown)
  U <- idx$U; D <- idx$D
  k <- length(U) + length(D)
  if (k < 1L) stop("query has no genes in common with refMatrix")

  s_i <- .camsum_col_sd(refMatrix, chunk)
  bad <- !is.finite(s_i) | s_i <= 0
  if (any(bad)) {
    warning(sum(bad), " profile(s) have zero or non-finite standard ",
            "deviation; their CamSum score is NA")
    s_i[bad] <- NA_real_
  }

  # One rho_bar per perturbation type; "unknown" when the type is not known
  if (!is.null(rho_bar)) {
    group <- rep("all", ncol(refMatrix))
  } else {
    if (is.null(pert_type)) pert_type <- .camsum_pert_types(colnames(refMatrix))
    if (length(pert_type) != ncol(refMatrix)) {
      stop("pert_type must have one entry per column of refMatrix")
    }
    group <- ifelse(is.na(pert_type), "unknown", as.character(pert_type))
  }
  rho_for_group <- function(g) {
    if (!is.null(rho_bar)) {
      return(list(rho = as.numeric(rho_bar), source = "supplied", n = NA_real_))
    }
    r <- if (g != "unknown") .camsum_rho_gctx(g, genes, U, D, chunk)
    if (!is.null(r)) {
      return(list(rho = as.numeric(r), source = "library_gctx",
                  n = attr(r, "n_profiles")))
    }
    cols <- which(group == g)
    r <- .camsum_rho(refMatrix[, cols, drop = FALSE], U, D, s_i[cols], chunk)
    list(rho = as.numeric(r), source = "refMatrix", n = attr(r, "n_profiles"))
  }
  per_group <- lapply(stats::setNames(nm = unique(group)), rho_for_group)
  rho   <- vapply(per_group, `[[`, numeric(1), "rho")
  src   <- vapply(per_group, `[[`, character(1), "source")
  n_rho <- vapply(per_group, function(x) as.numeric(x$n), numeric(1))
  vif   <- ifelse(is.finite(rho), pmax(1, 1 + (k - 1) * rho), 1)

  local <- names(src)[src == "refMatrix"]
  if (length(local)) {
    small <- local[n_rho[local] < 500]
    warning("CamSum rho_bar was estimated on the supplied reference (",
            paste(sprintf("%s: %d profiles", local, as.integer(n_rho[local])),
                  collapse = ", "),
            ") instead of the full library, because the perturbation type ",
            "was unknown, options CONCERTDR.gctx_file / siginfo_file / ",
            "geneinfo_file were not set, a reference gene was missing from ",
            "the GCTX file, or the GCTX file does not hold every profile of ",
            "that type. Rankings are not affected, but p-values can ",
            "depart from nominal levels when the reference is a filtered ",
            "subset (e.g. one cell line) and the query has many genes; pass ",
            "camsum_rho_bar to override.",
            if (length(small)) {
              paste0(" Fewer than 500 profiles were available for: ",
                     paste(small, collapse = ", "))
            },
            call. = FALSE)
  }

  raw <- colSums(refMatrix[U, , drop = FALSE]) -
         colSums(refMatrix[D, , drop = FALSE])
  score <- raw / (s_i * sqrt(k * unname(vif[group])))
  pValue <- switch(alternative,
    two.sided = 2 * stats::pnorm(-abs(score)),
    greater   = stats::pnorm(score, lower.tail = FALSE),
    less      = stats::pnorm(score, lower.tail = TRUE))

  out <- data.frame(Score = score, pValue = pValue,
                    pAdjValue = stats::p.adjust(pValue, method = pAdjMethod))
  rownames(out) <- colnames(refMatrix)
  attributes(out)[c("k", "kU", "kD", "rho_bar", "VIF", "rho_source",
                    "rho_n_profiles", "alternative")] <-
    list(k, length(U), length(D), rho, vif, src, n_rho, alternative)
  out
}
