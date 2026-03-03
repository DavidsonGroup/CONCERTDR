#' Fast GCTX File I/O Utilities
#'
#' @description
#' Lightweight functions for reading CMap GCTX files using \pkg{rhdf5} directly,
#' without depending on \pkg{cmapR}.  \code{fast_parse_gctx} returns a plain
#' numeric matrix whose row-names and column-names are the HDF5 row/column ids.
#'
#' @name gctx_io
#' @keywords internal
NULL

#' Read Column or Row Metadata from a GCTX File
#'
#' @param fname Path to the GCTX file.
#' @param dim   Either \code{"col"} (default) or \code{"row"}.
#' @param ids   Optional character vector of ids to keep (preserving order).
#' @return A \code{data.frame} with an \code{id} column plus all other metadata
#'   fields stored in the HDF5 file.
#' @keywords internal
fast_read_meta <- function(fname, dim = "col", ids = NULL) {
  path    <- if (dim == "col") "0/META/COL" else "0/META/ROW"
  all_ids <- as.character(rhdf5::h5read(fname, paste0(path, "/id")))

  fields <- rhdf5::h5ls(fname)
  fields <- fields[fields$group == paste0("/", path) &
                     fields$name != "id", "name"]

  df <- data.frame(id = all_ids, stringsAsFactors = FALSE)
  for (f in fields) {
    v      <- rhdf5::h5read(fname, paste0(path, "/", f))
    df[[f]] <- as.character(v)
  }

  if (!is.null(ids)) {
    df <- df[match(ids, df$id), ]
    rownames(df) <- NULL
  }
  df
}

#' Parse a GCTX File and Return a Numeric Matrix
#'
#' Reads a subset of rows and/or columns from a GCTX (HDF5) file and returns
#' a plain numeric matrix.  Column subsetting is performed at the HDF5 level
#' (fast axis), while row subsetting is performed in memory (slow axis).
#'
#' @param fname Path to the GCTX file.
#' @param rid   Character vector of row ids, integer vector of row indices, or
#'   \code{NULL} (all rows).
#' @param cid   Character vector of column ids, integer vector of column
#'   indices, or \code{NULL} (all columns).
#' @return A numeric matrix with \code{rownames} set to row ids and
#'   \code{colnames} set to column ids.
#' @keywords internal
fast_parse_gctx <- function(fname, rid = NULL, cid = NULL) {
  rhdf5::h5closeAll()
  on.exit(rhdf5::h5closeAll())

  all_rid <- as.character(rhdf5::h5read(fname, "0/META/ROW/id"))
  all_cid <- as.character(rhdf5::h5read(fname, "0/META/COL/id"))

  # Column indices (fast axis – pushed into HDF5)
  cidx <- if (is.null(cid)) {
    seq_along(all_cid)
  } else {
    if (is.character(cid)) match(cid, all_cid) else cid
  }
  cidx <- sort(cidx)

  # Row indices (slow axis – filtered in memory)
  ridx <- if (is.null(rid)) {
    seq_along(all_rid)
  } else {
    if (is.character(rid)) match(rid, all_rid) else rid
  }

  # Read all rows × target columns
  m <- rhdf5::h5read(fname, "0/DATA/0/matrix",
                     index = list(seq_along(all_rid), cidx))

  # Subset rows, then label
  m <- m[ridx, , drop = FALSE]
  rownames(m) <- all_rid[ridx]
  colnames(m) <- all_cid[cidx]

  m
}
