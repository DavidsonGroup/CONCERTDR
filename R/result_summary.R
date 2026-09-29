# pert_type -> broad perturbation kind / mode lookups
.PERT_TYPE_KIND <- c(
  "trt_cp" = "Drug",
  "trt_lig" = "Drug",
  "trt_sh" = "Gene",
  "trt_sh.cgs" = "Gene",
  "trt_sh.css" = "Gene",
  "trt_xpr" = "Gene",
  "trt_oe" = "Gene",
  "trt_oe.mut" = "Gene",
  "ctl_vehicle" = "Control",
  "ctl_vector" = "Control",
  "ctl_vehicle.cns" = "Control",
  "ctl_vector.cns" = "Control",
  "ctl_untrt.cns" = "Control",
  "ctl_untrt" = "Control"
)

.PERT_TYPE_MODE <- c(
  "trt_xpr" = "XPR",
  "trt_sh" = "XPR",
  "trt_sh.cgs" = "XPR",
  "trt_sh.css" = "XPR",
  "trt_oe" = "OE",
  "trt_oe.mut" = "OE",
  "trt_cp" = "Drug",
  "trt_lig" = "Drug"
)

#' Parse library, cell line, dose and time from CMap signature ids
#'
#' Ids look like \code{LIB_CELL_6H:BRD-XXXX:dose:time}.
#' @param df Data frame.
#' @param col Column holding the signature ids.
#' @return \code{df} with added \code{library}, \code{cell_line},
#'   \code{broad_id}, \code{dose_uM} and \code{time_h} columns.
#' @keywords internal
.parse_compound_context <- function(df, col = "compound") {
  if (!col %in% names(df)) stop("Missing column for compound parsing: ", col)

  x <- as.character(df[[col]])
  x[is.na(x)] <- ""
  nth <- function(parts, n) {
    vapply(parts, function(v) if (length(v) >= n) v[[n]] else "", character(1))
  }
  blank_to_na <- function(v) ifelse(v == "", NA_character_, v)

  parts <- strsplit(x, ":", fixed = TRUE)
  left_parts <- strsplit(nth(parts, 1), "_", fixed = TRUE)
  time_token <- nth(left_parts, 3)

  time_token_num <- .as_numeric(sub(".*?(\\d+).*", "\\1", time_token))
  time_token_num[!grepl("\\d+", time_token)] <- NA_real_
  time_from_4th <- .as_numeric(nth(parts, 4))

  out <- df
  out$library <- blank_to_na(nth(left_parts, 1))
  out$cell_line <- blank_to_na(nth(left_parts, 2))
  out$broad_id <- blank_to_na(nth(parts, 2))
  out$dose_uM <- .as_numeric(nth(parts, 3))
  out$time_h <- ifelse(!is.na(time_from_4th), time_from_4th, time_token_num)
  out
}

#' Add perturbation kind, mode and dose unit from \code{pert_type}
#' @param df Data frame with a \code{pert_type} column.
#' @return \code{df} with added \code{pert_kind}, \code{mode} and
#'   \code{dose_unit} columns.
#' @keywords internal
.add_pert_kind_and_mode <- function(df) {
  pt <- tolower(trimws(as.character(df$pert_type)))
  out <- df
  out$pert_kind <- unname(.PERT_TYPE_KIND[pt])
  out$pert_kind[is.na(out$pert_kind)] <- "Unknown"
  out$mode <- unname(.PERT_TYPE_MODE[pt])
  out$mode[is.na(out$mode)] <- "Other"
  out$dose_unit <- ifelse(out$pert_kind == "Drug", "uM", "construct/NA")
  out
}

#' Collapse drug rows to one per drug
#'
#' Keeps the complete row with the best score in the requested direction.
#' @param df Data frame.
#' @param score_col Name of the score column.
#' @param group_keys Columns defining a drug (those present in \code{df}).
#' @return Data frame with one row per group.
#' @param direction \code{"reversal"} (lowest score first) or \code{"mimic"}.
#' @keywords internal
.dedup_drug_rows <- function(df, score_col, group_keys = c("pert_id", "display_name"),
                             direction = "reversal") {
  if (nrow(df) == 0) return(df)
  keys <- intersect(group_keys, names(df))
  if (length(keys) == 0) return(df)

  # Select the complete best context; never combine its score with other doses.
  df <- .order_by_score(df, score_col, direction)
  out <- df[!duplicated(df[, keys, drop = FALSE]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Whether the \code{pert_id} column is absent or entirely blank
#' @param df Data frame.
#' @return Logical scalar.
#' @keywords internal
.pert_id_missing <- function(df) {
  !"pert_id" %in% names(df) || all(.is_blank(df$pert_id))
}

#' Merge signature results with siginfo annotations
#' @param res Results data frame.
#' @param sig_info_file Siginfo path or data frame.
#' @param compound_col Signature id column in \code{res}.
#' @return Merged data frame with a guaranteed \code{pert_id} column.
#' @keywords internal
.merge_siginfo <- function(res, sig_info_file, compound_col) {
  sig_ids <- unique(as.character(res[[compound_col]]))
  sig_ids <- sig_ids[!.is_blank(sig_ids)]

  sig_cols <- c(
    "sig_id", "pert_type", "pert_id", "pert_iname", "pert_name",
    "cmap_name", "phase", "compound_aliases", "canonical_smiles",
    "inchi_key"
  )
  siginfo <- .read_cmap_table(sig_info_file, "sig_info_file",
                              select = sig_cols, sep = "\t")
  if (!all(c("sig_id", "pert_type") %in% names(siginfo))) {
    stop("sig_info_file must contain columns: sig_id, pert_type")
  }
  siginfo <- siginfo[siginfo$sig_id %in% sig_ids, , drop = FALSE]

  res2 <- merge(res, siginfo, by.x = compound_col, by.y = "sig_id", all.x = TRUE)
  if (!"pert_type" %in% names(res2) || all(is.na(res2$pert_type))) {
    stop("pert_type missing after merge; check sig_info_file")
  }
  if (.pert_id_missing(res2)) {
    res2$pert_id <- extract_compound_id(as.character(res2[[compound_col]]),
                                        method = "split_colon", part_index = 2)
  }
  if (.pert_id_missing(res2)) {
    stop("Could not determine 'pert_id'. Ensure sig_info_file contains 'pert_id' or compound strings include a BRD ID in colon-delimited format.")
  }
  res2
}

#' Merge compound annotations (target, MOA) by \code{pert_id}
#'
#' Multiple annotation rows are collapsed per perturbation: distinct nonblank
#' values in each column are joined with \code{"; "} in input order. Target
#' and MOA columns are independent lists, not positional target-MOA pairs.
#' @param res2 Data frame with a \code{pert_id} column.
#' @param comp_info_file Compound info path or data frame.
#' @return Merged data frame.
#' @keywords internal
.merge_compoundinfo <- function(res2, comp_info_file) {
  pert_ids <- unique(as.character(res2$pert_id))
  pert_ids <- pert_ids[!.is_blank(pert_ids)]

  compinfo <- .read_cmap_table(comp_info_file, "comp_info_file",
                               select = c("pert_id", "target", "moa"), sep = "\t")
  if (!"pert_id" %in% names(compinfo)) {
    stop("comp_info_file must contain column 'pert_id'")
  }
  compinfo <- compinfo[compinfo$pert_id %in% pert_ids, , drop = FALSE]
  if (anyDuplicated(compinfo$pert_id)) {
    groups <- split(seq_len(nrow(compinfo)), as.character(compinfo$pert_id))
    compinfo <- do.call(rbind, lapply(groups, function(idx) {
      row <- compinfo[idx[1L], , drop = FALSE]
      for (col in intersect(c("target", "moa"), names(compinfo))) {
        values <- trimws(as.character(compinfo[[col]][idx]))
        values <- unique(values[!.is_blank(values)])
        row[[col]] <- if (length(values)) paste(values, collapse = "; ") else NA_character_
      }
      row
    }))
    rownames(compinfo) <- NULL
  }

  merge(res2, compinfo, by = "pert_id", all.x = TRUE, suffixes = c("", "_comp"))
}

#' Fill blank entries of a name column from \code{pert_iname}
#' @param df Data frame.
#' @param col Name column; created as \code{NA} when absent.
#' @return \code{df}.
#' @keywords internal
.fill_name_from_iname <- function(df, col) {
  if (!col %in% names(df)) df[[col]] <- NA_character_
  if ("pert_iname" %in% names(df)) {
    blank <- .is_blank(df[[col]])
    df[[col]][blank] <- as.character(df$pert_iname[blank])
  }
  df
}

#' Derive names, context, perturbation kind, effect direction and MOA status
#' @param res2 Merged data frame.
#' @param compound_col Signature id column.
#' @param score_col Score column.
#' @return \code{res2} with the derived columns added.
#' @keywords internal
.derive_annotations <- function(res2, compound_col, score_col) {
  res2 <- .fill_name_from_iname(res2, "pert_name")
  res2 <- .fill_name_from_iname(res2, "cmap_name")
  res2 <- .parse_compound_context(res2, col = compound_col)
  res2 <- .add_pert_kind_and_mode(res2)

  # display_name: cmap_name -> pert_name -> pert_id
  cmap <- as.character(res2$cmap_name)
  pname <- as.character(res2$pert_name)
  res2$display_name <- ifelse(!.is_blank(cmap), cmap,
                              ifelse(!.is_blank(pname), pname,
                                     as.character(res2$pert_id)))
  res2$sig_id <- as.character(res2[[compound_col]])
  res2$perturbation_name <- as.character(res2$display_name)

  score_num <- .as_numeric(res2[[score_col]])
  res2$effect_direction <- ifelse(
    score_num < 0,
    "Reversal (potentially therapeutic)",
    ifelse(score_num > 0, "Mimic/Aggravating", "Neutral")
  )

  moa <- if ("moa" %in% names(res2)) as.character(res2$moa) else rep(NA_character_, nrow(res2))
  res2$moa_status <- ifelse(!.is_blank(moa), "Known", "Unknown")
  res2
}

#' Full technical view with preferred columns first
#' @param res2 Annotated data frame.
#' @param compound_col,score_col,p_col,padj_col Column names.
#' @return Data frame.
#' @keywords internal
.build_tech_view <- function(res2, compound_col, score_col, p_col, padj_col) {
  hidden <- c("display_name", "cmap_name", "pert_name", "pert_iname", "broad_id")
  if (!identical(compound_col, "sig_id")) hidden <- c(hidden, compound_col)
  preferred <- c(
    "sig_id", "perturbation_name", "pert_id", "pert_type", "pert_kind", "mode",
    "compound_aliases", "target", "moa", "moa_status", "phase",
    "library", "cell_line", "time_h", "dose_uM", "dose_unit",
    score_col, "effect_direction", p_col, padj_col, "rank",
    "canonical_smiles", "inchi_key"
  )
  tech_cols <- c(intersect(preferred, names(res2)),
                 setdiff(names(res2), c(preferred, hidden)))
  res2[, unique(tech_cols), drop = FALSE]
}

#' Order rows by ascending numeric score
#' @param df Data frame.
#' @param score_col Score column.
#' @return Data frame.
#' @param direction \code{"reversal"} (lowest score first) or \code{"mimic"}.
#' @keywords internal
.order_by_score <- function(df, score_col, direction = "reversal") {
  if (nrow(df) > 0) df <- df[order(.as_numeric(df[[score_col]]), decreasing = direction == "mimic"), , drop = FALSE]
  df
}

#' Drug-focused wetlab view (one row per drug, best score first)
#' @param res2 Annotated data frame.
#' @param score_col Score column.
#' @param keep_dose Logical; keep \code{dose_uM}.
#' @return Data frame.
#' @param direction \code{"reversal"} (lowest score first) or \code{"mimic"}.
#' @keywords internal
.build_wetlab_drug_view <- function(res2, score_col, keep_dose, direction = "reversal") {
  drug_cols <- c(
    "perturbation_name", "pert_id", "pert_type", "pert_kind",
    "compound_aliases", "target", "moa", "phase",
    score_col, "effect_direction", "cell_line", "time_h"
  )
  if (keep_dose) drug_cols <- c(drug_cols, "dose_uM")
  drug_cols <- intersect(drug_cols, names(res2))

  drugs <- res2[res2$pert_kind == "Drug", drug_cols, drop = FALSE]
  drugs <- .dedup_drug_rows(drugs, score_col = score_col,
                            group_keys = c("pert_id", "perturbation_name"), direction = direction)
  .order_by_score(drugs, score_col, direction)
}

#' Gene-focused wetlab view (best score first)
#' @param res2 Annotated data frame.
#' @param score_col Score column.
#' @return Data frame.
#' @param direction \code{"reversal"} (lowest score first) or \code{"mimic"}.
#' @keywords internal
.build_wetlab_gene_view <- function(res2, score_col, direction = "reversal") {
  gene_cols <- intersect(c(
    "perturbation_name", "pert_id", "pert_type", "pert_kind", "mode",
    score_col, "effect_direction", "cell_line", "time_h", "library"
  ), names(res2))
  .order_by_score(res2[res2$pert_kind == "Gene", gene_cols, drop = FALSE], score_col, direction)
}

#' Per-drug context summary
#' @param res2 Annotated data frame.
#' @param score_col Score column.
#' @return Data frame (empty when there are no drug rows).
#' @param direction \code{"reversal"} (lowest score first) or \code{"mimic"}.
#' @keywords internal
.build_drug_context_summary <- function(res2, score_col, direction = "reversal") {
  drug_only <- res2[res2$pert_kind == "Drug", , drop = FALSE]
  if (nrow(drug_only) == 0) return(data.frame())

  groups <- split(seq_len(nrow(drug_only)), as.character(drug_only$pert_id))
  idx_best <- unlist(lapply(groups, function(idx) {
    sc <- .as_numeric(drug_only[[score_col]][idx])
    idx[order(sc, decreasing = direction == "mimic", na.last = TRUE)][1]
  }), use.names = FALSE)

  summary_cols <- intersect(
    c("perturbation_name", "pert_id", "pert_type", "pert_kind", "phase",
      "moa_status", score_col),
    names(drug_only)
  )
  drug_summary <- drug_only[idx_best, summary_cols, drop = FALSE]
  names(drug_summary)[names(drug_summary) == score_col] <- "best_score"

  pid <- as.character(drug_summary$pert_id)
  drug_summary$n_contexts <- as.integer(table(as.character(drug_only$pert_id))[pid])
  n_unique <- function(col) {
    n <- tapply(drug_only[[col]], drug_only$pert_id, function(x) length(unique(x[!is.na(x)])))
    as.integer(n[pid])
  }
  if ("cell_line" %in% names(drug_only)) drug_summary$n_cell_lines <- n_unique("cell_line")
  if ("time_h" %in% names(drug_only)) drug_summary$n_time <- n_unique("time_h")

  ord_score <- .as_numeric(drug_summary$best_score)
  drug_summary[order(-drug_summary$n_contexts, if (direction == "mimic") -ord_score else ord_score), , drop = FALSE]
}

#' Write the four result tables as TSV files
#' @param tables Named list with \code{wetlab_drug_view},
#'   \code{wetlab_gene_view}, \code{tech_view_all} and
#'   \code{drug_context_summary}.
#' @param output_dir Output directory (created when missing).
#' @param verbose Logical; list written files.
#' @return Named character vector of file paths.
#' @keywords internal
.write_result_tables <- function(tables, output_dir, verbose) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  paths <- file.path(output_dir, paste0(names(tables), ".tsv"))
  names(paths) <- names(tables)
  for (nm in names(tables)) {
    utils::write.table(tables[[nm]], paths[[nm]], sep = "\t",
                       row.names = FALSE, quote = FALSE, na = "")
  }
  if (verbose) {
    message("Wrote output files:")
    for (p in paths) message(" - ", p)
  }
  paths
}


#' Annotate CONCERTDR Results with Drug Information
#'
#' @description
#' Integrate signature matching results with annotations from
#' \code{siginfo_beta.txt} and \code{compoundinfo_beta.txt}. It produces
#' three analysis-ready tables and one drug context summary table.
#' @details When compoundinfo contains multiple rows for one \code{pert_id},
#'   all distinct nonblank target and MOA values are retained, separated by
#'   \code{"; "}. This does not duplicate signature results or inflate context
#'   counts. Targets and MOAs are independent lists, not positional pairs.
#'
#' @param direction \code{"reversal"} (default) puts the lowest scores first
#'   in the wetlab views and picks each drug's lowest-scoring context as its
#'   best; \code{"mimic"} uses the highest. It only orders and picks; no row
#'   is dropped by the sign of its score.
#' @param results_df Data frame containing signature matching results with a 'compound' column
#' @param sig_info_file Path to siginfo_beta.txt file or data frame with signature information
#' @param comp_info_file Path to compound information file or data frame
#' @param output_file Deprecated single-file output path; when provided,
#'   \code{tech_view_all} is also written as CSV for compatibility.
#' @param score_col Score column in \code{results_df} (default: \code{"Score"})
#' @param padj_col Adjusted p-value column name (default: \code{"pAdjValue"})
#' @param p_col Raw p-value column name (default: \code{"pValue"})
#' @param compound_col Compound/sig-id column in \code{results_df} (default: \code{"compound"})
#' @param keep_dose_in_drug_view Logical; whether to keep \code{dose_uM} in
#'   \code{wetlab_drug_view} (default: \code{FALSE})
#' @param output_dir Optional directory to write four TSV outputs. If NULL,
#'   files are not written unless \code{write_outputs = TRUE}.
#' @param write_outputs Logical; write four TSV outputs (default: \code{FALSE})
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return A named list with:
#'   \item{wetlab_drug_view}{Drug-focused wetlab view}
#'   \item{wetlab_gene_view}{Gene-focused wetlab view}
#'   \item{tech_view_all}{Full technical table with all integrated fields}
#'   \item{drug_context_summary}{Per-drug context summary}
#'   \item{output_files}{Named character vector of written files (if any)}
#'
#' @examples
#' ex_results <- data.frame(
#'   compound = "CVD001_HEPG2_6H:BRD-K03652504-001-01-9:10.0497",
#'   Score = -0.72,
#'   pValue = 0.002,
#'   pAdjValue = 0.02
#' )
#' ex_siginfo <- data.frame(
#'   sig_id = ex_results$compound,
#'   pert_type = "trt_cp",
#'   pert_id = "BRD-K03652504-001-01-9",
#'   pert_iname = "imatinib"
#' )
#' ex_compinfo <- data.frame(
#'   pert_id = "BRD-K03652504-001-01-9",
#'   pert_name = "imatinib",
#'   cmap_name = "imatinib"
#' )
#' views <- annotate_drug_results(
#'   results_df = ex_results,
#'   sig_info_file = ex_siginfo,
#'   comp_info_file = ex_compinfo,
#'   write_outputs = FALSE,
#'   verbose = FALSE
#' )
#' names(views)
#'
#' \donttest{
#' if (file.exists("sig_match_xsum_results.csv") &&
#'     file.exists("siginfo_beta.txt") &&
#'     file.exists("compoundinfo_beta.txt")) {
#'   results <- read.csv("sig_match_xsum_results.csv")
#'
#'   views <- annotate_drug_results(
#'     results_df = results,
#'     sig_info_file = "siginfo_beta.txt",
#'     comp_info_file = "compoundinfo_beta.txt",
#'     output_dir = "results",
#'     write_outputs = TRUE
#'   )
#'
#'   head(views$wetlab_drug_view)
#' }
#' }
#'
#' @export
annotate_drug_results <- function(results_df,
                                  sig_info_file,
                                  comp_info_file,
                                  output_file = NULL,
                                  score_col = "Score",
                                  padj_col = "pAdjValue",
                                  p_col = "pValue",
                                  compound_col = "compound",
                                  keep_dose_in_drug_view = FALSE,
                                  output_dir = NULL,
                                  write_outputs = FALSE,
                                  verbose = TRUE,
                                  direction = c("reversal", "mimic")) {
  direction <- match.arg(direction)
  if (!is.data.frame(results_df)) stop("results_df must be a data.frame")
  if (!compound_col %in% names(results_df)) {
    stop("results_df must contain column: ", compound_col)
  }
  if (!score_col %in% names(results_df)) {
    stop("results_df must contain score column: ", score_col)
  }
  if (verbose) message("Integrating signature results with siginfo and compoundinfo...")

  res2 <- .merge_siginfo(results_df, sig_info_file, compound_col)
  res2 <- .merge_compoundinfo(res2, comp_info_file)
  res2 <- .derive_annotations(res2, compound_col, score_col)

  tech_view <- .build_tech_view(res2, compound_col, score_col, p_col, padj_col)
  tables <- list(
    wetlab_drug_view = .build_wetlab_drug_view(res2, score_col, keep_dose_in_drug_view, direction),
    wetlab_gene_view = .build_wetlab_gene_view(res2, score_col, direction),
    tech_view_all = tech_view,
    drug_context_summary = .build_drug_context_summary(res2, score_col, direction)
  )

  output_files <- character()
  if (isTRUE(write_outputs) || !is.null(output_dir)) {
    if (is.null(output_dir)) output_dir <- "."
    output_files <- .write_result_tables(tables, output_dir, verbose)
  }

  # Backward compatibility: legacy single output file
  if (!is.null(output_file)) {
    utils::write.csv(tech_view, file = output_file, row.names = FALSE)
    if (verbose) message("Wrote legacy output_file (tech_view_all as CSV): ", output_file)
  }

  c(tables, list(output_files = output_files))
}

#' Extract Compound ID from Signature Matching Results
#'
#' @description
#' Helper function to extract compound identifiers from complex compound strings
#' in signature matching results. This function can be customized based on the
#' specific format of compound identifiers in your data.
#'
#' @param compound_strings Vector of compound identifier strings
#' @param method Method for extraction: "split_colon" (default), "split_underscore", or "regex"
#' @param regex_pattern Regular expression pattern for extraction (used when method="regex")
#' @param part_index Which part to extract when splitting (default: 2)
#'
#' @return Vector of extracted compound identifiers
#'
#' @examples
#' extract_compound_id("A:BRD-K03652504-001-01-9:10")
#'
#' \donttest{
#' # Example compound strings
#' compounds <- c("CVD001_HEPG2_6H:BRD-K03652504-001-01-9:10.0497",
#'                "CVD001_HEPG2_6H:BRD-A37828317-001-03-0:10")
#'
#' # Extract using colon splitting (default)
#' ids <- extract_compound_id(compounds)
#'
#' # Extract using custom regex
#' ids <- extract_compound_id(compounds, method = "regex",
#'                           regex_pattern = "BRD-[A-Z0-9-]+")
#' }
#'
#' @export
extract_compound_id <- function(compound_strings,
                                method = "split_colon",
                                regex_pattern = NULL,
                                part_index = 2) {
  split_at <- function(sep) {
    vapply(strsplit(compound_strings, sep, fixed = TRUE), function(x) {
      if (length(x) >= part_index) x[part_index] else x[1]
    }, character(1))
  }
  switch(method,
    split_colon = split_at(":"),
    split_underscore = split_at("_"),
    regex = {
      if (is.null(regex_pattern)) {
        stop("regex_pattern must be provided when method='regex'")
      }
      locations <- regexpr(regex_pattern, compound_strings)
      matched <- !is.na(locations) & locations > 0L
      out <- compound_strings
      out[matched] <- regmatches(compound_strings, locations)
      out
    },
    stop("Invalid method. Must be 'split_colon', 'split_underscore', or 'regex'")
  )
}

#' Fuzzy String Matching for Drug Names
#'
#' @description
#' Performs fuzzy string matching to find the best matches between query drug names
#' and a reference list of drug names. Uses Levenshtein distance for similarity calculation.
#'
#' @param query_names Vector of query drug names
#' @param reference_names Vector of reference drug names to match against
#' @param method Similarity method: "levenshtein" (default), "jaro", or "jarowinkler"
#' @param threshold Minimum similarity threshold (0-100)
#' @param top_n Number of top matches to return for each query (default: 1)
#'
#' @return Data frame with query names, matched names, and similarity scores
#'
#' @examples
#' if (requireNamespace("RecordLinkage", quietly = TRUE)) {
#'   fuzzy_drug_match(c("asprin"), c("aspirin", "ibuprofen"), threshold = 70)
#' }
#'
#' \donttest{
#' if (requireNamespace("RecordLinkage", quietly = TRUE)) {
#'   query <- c("aspirin", "ibuprofen", "acetaminophen")
#'   reference <- c("aspirin", "ibuprofen", "acetaminophen", "naproxen", "diclofenac")
#'   matches <- fuzzy_drug_match(query, reference, threshold = 80)
#'   print(matches)
#' }
#' }
#'
#' @export
fuzzy_drug_match <- function(query_names,
                             reference_names,
                             method = "levenshtein",
                             threshold = 80,
                             top_n = 1) {
  if (!requireNamespace("RecordLinkage", quietly = TRUE)) {
    stop("Package 'RecordLinkage' is required for fuzzy matching. Please install it with: install.packages('RecordLinkage')")
  }

  matches <- list()
  for (query in query_names) {
    if (is.na(query) || query == "") next

    similarities <- switch(method,
      levenshtein = RecordLinkage::levenshteinSim(query, reference_names),
      jaro = RecordLinkage::jarowinkler(query, reference_names),
      jarowinkler = RecordLinkage::jarowinkler(query, reference_names, W = 0.1),
      stop("Invalid method. Must be 'levenshtein', 'jaro', or 'jarowinkler'")
    ) * 100

    above <- which(similarities >= threshold)
    if (length(above) == 0) next
    top <- utils::head(above[order(similarities[above], decreasing = TRUE)], top_n)
    matches[[length(matches) + 1L]] <- data.frame(
      query_name = query,
      matched_name = reference_names[top],
      similarity_score = similarities[top]
    )
  }

  if (length(matches) == 0) {
    return(data.frame(query_name = character(),
                      matched_name = character(),
                      similarity_score = numeric()))
  }
  out <- do.call(rbind, matches)
  rownames(out) <- NULL
  out
}
