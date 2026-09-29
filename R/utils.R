#' Read a CMap table from a path or pass a data frame through
#'
#' Column names are trimmed of surrounding whitespace.
#' @param x File path or data frame.
#' @param what Name of the calling function's argument, used in error messages
#'   (for example \code{"geneinfo_file"}).
#' @param select Optional character vector of columns to keep. Columns that
#'   are absent from the table are dropped; errors if none are present.
#' @param sep Field separator passed to \code{data.table::fread}
#'   (default: auto-detect).
#' @return Data frame.
#' @keywords internal
.read_cmap_table <- function(x, what, select = NULL, sep = "auto") {
  if (is.data.frame(x)) {
    df <- x
  } else if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop(what, " not found: ", x)
    selected <- NULL
    if (!is.null(select)) {
      header <- names(data.table::fread(x, sep = sep, header = TRUE,
                                        nrows = 0, showProgress = FALSE))
      selected <- which(trimws(header) %in% select)
      if (!length(selected)) {
        stop(what, " contains none of the requested columns: ",
             paste(select, collapse = ", "))
      }
    }
    df <- data.table::fread(x, sep = sep, header = TRUE, data.table = FALSE,
                            select = selected, showProgress = FALSE)
  } else {
    stop(what, " must be either a file path (character) or a data.frame.")
  }
  names(df) <- trimws(names(df))
  if (!is.null(select)) df <- df[, intersect(select, names(df)), drop = FALSE]
  df
}

#' Print a horizontal rule of \code{n} copies of \code{char} as a message
#' @param char Single character.
#' @param n Width.
#' @return NULL, invisibly.
#' @keywords internal
.message_rule <- function(char = "=", n = 80L) {
  message(strrep(char, n))
}

#' First existing file among candidates
#' @param ... Candidate paths; \code{NULL}, \code{NA} and empty strings are
#'   skipped.
#' @return A path, or \code{NULL} when none of the candidates exists.
#' @keywords internal
.first_existing_file <- function(...) {
  candidates <- as.character(unlist(list(...)))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) candidates[[1L]] else NULL
}

#' Coerce to numeric without warnings
#' @param x Vector.
#' @return Numeric vector; \code{NA} where \code{x} is not numeric.
#' @keywords internal
.as_numeric <- function(x) suppressWarnings(as.numeric(x))

#' Whether strings are missing or empty
#' @param x Vector coercible to character.
#' @return Logical vector.
#' @keywords internal
.is_blank <- function(x) is.na(x) | !nzchar(as.character(x))

#' Order rows by score in the requested connectivity direction
#'
#' Rows are only ordered, never dropped by the sign of the score; rows with a
#' non-finite score are removed because they cannot be ranked.
#' @param df Data frame containing scores.
#' @param score_col Score column name.
#' @param direction \code{"reversal"} (lowest score first) or \code{"mimic"}
#'   (highest score first).
#' @return Data frame in best-score-first order.
#' @keywords internal
.rank_by_direction <- function(df, score_col = "Score", direction = c("reversal", "mimic")) {
  direction <- match.arg(direction)
  score <- .as_numeric(df[[score_col]])
  df <- df[is.finite(score), , drop = FALSE]
  df[order(.as_numeric(df[[score_col]]), decreasing = direction == "mimic"), , drop = FALSE]
}
