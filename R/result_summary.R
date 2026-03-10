#' Annotate CONCERTDR Results with Drug Information
#'
#' @description
#' Integrate signature matching results with annotations from
#' \code{siginfo_beta.txt} and \code{compoundinfo_beta.txt}. It produces
#' three analysis-ready tables and one drug context summary table.
#'
#' @param results_df Data frame containing signature matching results with a 'compound' column
#' @param sig_info_file Path to siginfo_beta.txt file or data frame with signature information
#' @param comp_info_file Path to compound information file or data frame
#' @param drug_info_file Deprecated. Kept for backward compatibility; not used.
#' @param output_file Deprecated single-file output path; when provided,
#'   \code{tech_view_all} is also written as CSV for compatibility.
#' @param fuzzy_threshold Deprecated. Not used.
#' @param perfect_match_only Deprecated. Not used.
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
#'   pAdjValue = 0.02,
#'   stringsAsFactors = FALSE
#' )
#' ex_siginfo <- data.frame(
#'   sig_id = ex_results$compound,
#'   pert_type = "trt_cp",
#'   pert_id = "BRD-K03652504-001-01-9",
#'   pert_iname = "imatinib",
#'   stringsAsFactors = FALSE
#' )
#' ex_compinfo <- data.frame(
#'   pert_id = "BRD-K03652504-001-01-9",
#'   pert_name = "imatinib",
#'   cmap_name = "imatinib",
#'   stringsAsFactors = FALSE
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
                                  drug_info_file = NULL,
                                  output_file = NULL,
                                  fuzzy_threshold = 90,
                                  perfect_match_only = FALSE,
                                  score_col = "Score",
                                  padj_col = "pAdjValue",
                                  p_col = "pValue",
                                  compound_col = "compound",
                                  keep_dose_in_drug_view = FALSE,
                                  output_dir = NULL,
                                  write_outputs = FALSE,
                                  verbose = TRUE) {
  if (!is.null(drug_info_file) && verbose) {
    message("'drug_info_file' is deprecated in the new integration workflow and will be ignored.")
  }
  if (!missing(fuzzy_threshold) && verbose) {
    message("'fuzzy_threshold' is deprecated and ignored.")
  }
  if (!missing(perfect_match_only) && verbose) {
    message("'perfect_match_only' is deprecated and ignored.")
  }

  # Local helpers
  read_or_use_df <- function(file_or_df, file_description, select_cols = NULL) {
    if (is.data.frame(file_or_df)) {
      df <- file_or_df
    } else if (is.character(file_or_df) && length(file_or_df) == 1L) {
      if (!file.exists(file_or_df)) {
        stop(file_description, " file not found: ", file_or_df)
      }

      if (requireNamespace("data.table", quietly = TRUE)) {
        if (is.null(select_cols)) {
          df <- data.table::fread(file_or_df, sep = "\t", header = TRUE,
                                  data.table = FALSE, stringsAsFactors = FALSE)
        } else {
          header_names <- names(data.table::fread(
            file_or_df,
            sep = "\t",
            header = TRUE,
            nrows = 0,
            data.table = FALSE,
            stringsAsFactors = FALSE,
            showProgress = FALSE
          ))
          keep_cols <- intersect(select_cols, trimws(header_names))
          df <- data.table::fread(file_or_df, sep = "\t", header = TRUE,
                                  data.table = FALSE, stringsAsFactors = FALSE,
                                  select = keep_cols)
        }
      } else {
        df <- utils::read.table(file_or_df, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE, quote = "",
                                comment.char = "", fill = TRUE)
        if (!is.null(select_cols)) {
          keep <- intersect(select_cols, names(df))
          df <- df[, keep, drop = FALSE]
        }
      }
    } else {
      stop(file_description, " must be either a data.frame or a file path")
    }

    names(df) <- trimws(names(df))
    df
  }

  parse_compound_context <- function(df, col = "compound") {
    if (!col %in% names(df)) {
      stop("Missing column for compound parsing: ", col)
    }

    x <- as.character(df[[col]])
    x[is.na(x)] <- ""

    parts <- strsplit(x, ":", fixed = TRUE)
    get_nth <- function(v, n) if (length(v) >= n) v[[n]] else ""

    left <- vapply(parts, get_nth, character(1), n = 1)
    broad_id <- vapply(parts, get_nth, character(1), n = 2)
    dose_raw <- vapply(parts, get_nth, character(1), n = 3)
    time_raw <- vapply(parts, get_nth, character(1), n = 4)

    left_parts <- strsplit(left, "_", fixed = TRUE)
    library <- vapply(left_parts, get_nth, character(1), n = 1)
    cell_line <- vapply(left_parts, get_nth, character(1), n = 2)
    time_token <- vapply(left_parts, get_nth, character(1), n = 3)

    time_token_num <- suppressWarnings(as.numeric(sub(".*?(\\d+).*", "\\1", time_token)))
    no_digit <- !grepl("\\d+", time_token)
    time_token_num[no_digit] <- NA_real_

    dose_uM <- suppressWarnings(as.numeric(dose_raw))
    time_h_from_4th <- suppressWarnings(as.numeric(time_raw))
    time_h <- ifelse(!is.na(time_h_from_4th), time_h_from_4th, time_token_num)

    out <- df
    out$library <- ifelse(library == "", NA_character_, library)
    out$cell_line <- ifelse(cell_line == "", NA_character_, cell_line)
    out$broad_id <- ifelse(broad_id == "", NA_character_, broad_id)
    out$dose_uM <- dose_uM
    out$time_h <- time_h
    out
  }

  add_pert_kind_and_mode <- function(df) {
    pert_type_to_kind <- c(
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

    pert_type_to_mode <- c(
      "trt_xpr" = "XPR",
      "trt_sh" = "XPR",
      "trt_sh.cgs" = "XPR",
      "trt_sh.css" = "XPR",
      "trt_oe" = "OE",
      "trt_oe.mut" = "OE",
      "trt_cp" = "Drug",
      "trt_lig" = "Drug"
    )

    pt <- tolower(trimws(as.character(df$pert_type)))
    out <- df
    out$pert_kind <- unname(pert_type_to_kind[pt])
    out$pert_kind[is.na(out$pert_kind)] <- "Unknown"
    out$mode <- unname(pert_type_to_mode[pt])
    out$mode[is.na(out$mode)] <- "Other"
    out$dose_unit <- ifelse(out$pert_kind == "Drug", "uM", "construct/NA")
    out
  }

  dedup_drug_rows <- function(df, score_col, group_keys = c("pert_id", "display_name")) {
    if (nrow(df) == 0) return(df)

    keys <- intersect(group_keys, names(df))
    if (length(keys) == 0) return(df)

    group_id <- interaction(df[, keys, drop = FALSE], drop = TRUE, lex.order = TRUE)
    groups <- split(df, group_id)

    collapse_unique <- function(x) {
      xv <- unique(as.character(x[!is.na(x) & as.character(x) != ""]))
      if (length(xv) == 0) NA_character_ else paste(sort(xv), collapse = "; ")
    }

    rows <- lapply(groups, function(g) {
      out <- g[1, , drop = FALSE]
      for (nm in names(g)) {
        v <- g[[nm]]
        if (nm %in% keys) {
          nz <- which(!is.na(v))[1]
          out[[nm]] <- if (is.na(nz)) NA else v[nz]
        } else if (nm == score_col) {
          num <- suppressWarnings(as.numeric(v))
          out[[nm]] <- if (all(is.na(num))) NA_real_ else min(num, na.rm = TRUE)
        } else if (is.numeric(v)) {
          nz <- which(!is.na(v))[1]
          out[[nm]] <- if (is.na(nz)) NA_real_ else v[nz]
        } else {
          out[[nm]] <- collapse_unique(v)
        }
      }
      out
    })

    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }

  if (!is.data.frame(results_df)) {
    stop("results_df must be a data.frame")
  }
  if (!compound_col %in% names(results_df)) {
    stop("results_df must contain column: ", compound_col)
  }
  if (!score_col %in% names(results_df)) {
    stop("results_df must contain score column: ", score_col)
  }

  if (verbose) message("Integrating signature results with siginfo and compoundinfo...")

  res <- results_df
  sig_ids <- unique(as.character(res[[compound_col]]))
  sig_ids <- sig_ids[!is.na(sig_ids) & nzchar(sig_ids)]

  sig_cols <- c(
    "sig_id", "pert_type", "pert_id", "pert_iname", "pert_name",
    "cmap_name", "phase", "compound_aliases", "canonical_smiles",
    "inchi_key"
  )
  siginfo <- read_or_use_df(sig_info_file, "Signature info", select_cols = sig_cols)
  siginfo <- siginfo[, intersect(sig_cols, names(siginfo)), drop = FALSE]

  if (!all(c("sig_id", "pert_type") %in% names(siginfo))) {
    stop("sig_info_file must contain columns: sig_id, pert_type")
  }

  siginfo <- siginfo[siginfo$sig_id %in% sig_ids, , drop = FALSE]

  res2 <- merge(res, siginfo, by.x = compound_col, by.y = "sig_id", all.x = TRUE)
  if (!"pert_type" %in% names(res2) || all(is.na(res2$pert_type))) {
    stop("pert_type missing after merge; check sig_info_file")
  }

  if (!"pert_id" %in% names(res2) || all(is.na(res2$pert_id) | !nzchar(as.character(res2$pert_id)))) {
    if (compound_col %in% names(res2)) {
      res2$pert_id <- extract_compound_id(as.character(res2[[compound_col]]), method = "split_colon", part_index = 2)
    }
  }

  if (!"pert_id" %in% names(res2) || all(is.na(res2$pert_id) | !nzchar(as.character(res2$pert_id)))) {
    stop("Could not determine 'pert_id'. Ensure sig_info_file contains 'pert_id' or compound strings include a BRD ID in colon-delimited format.")
  }

  pert_ids <- unique(as.character(res2$pert_id))
  pert_ids <- pert_ids[!is.na(pert_ids) & nzchar(pert_ids)]

  compound_cols <- c(
    "pert_id", "pert_name", "cmap_name", "target", "moa", "phase",
    "compound_aliases", "canonical_smiles", "inchi_key"
  )

  compinfo <- read_or_use_df(comp_info_file, "Compound info", select_cols = compound_cols)
  compinfo <- compinfo[, intersect(compound_cols, names(compinfo)), drop = FALSE]

  if (!"pert_id" %in% names(compinfo)) {
    stop("comp_info_file must contain column 'pert_id'")
  }

  compinfo <- compinfo[compinfo$pert_id %in% pert_ids, , drop = FALSE]
  compinfo <- compinfo[!duplicated(compinfo$pert_id), , drop = FALSE]

  res2 <- merge(res2, compinfo, by = "pert_id", all.x = TRUE, suffixes = c("", "_comp"))

  res2 <- parse_compound_context(res2, col = compound_col)
  res2 <- add_pert_kind_and_mode(res2)

  # display_name: cmap_name -> pert_name -> pert_id
  cmap <- if ("cmap_name" %in% names(res2)) as.character(res2$cmap_name) else rep(NA_character_, nrow(res2))
  pname <- if ("pert_name" %in% names(res2)) as.character(res2$pert_name) else rep(NA_character_, nrow(res2))
  pid <- as.character(res2$pert_id)
  res2$display_name <- ifelse(!is.na(cmap) & nzchar(cmap), cmap,
                              ifelse(!is.na(pname) & nzchar(pname), pname, pid))

  score_num <- suppressWarnings(as.numeric(res2[[score_col]]))
  res2$effect_direction <- ifelse(
    score_num < 0,
    "Reversal (potentially therapeutic)",
    ifelse(score_num > 0, "Mimic/Aggravating", "Neutral")
  )

  moa <- if ("moa" %in% names(res2)) as.character(res2$moa) else rep(NA_character_, nrow(res2))
  res2$moa_status <- ifelse(!is.na(moa) & nzchar(moa), "Known", "Unknown")

  # Tech view
  preferred <- c(
    "pert_type", "pert_kind", "mode", "display_name", "pert_id", "pert_name",
    "compound_aliases", "target", "moa", "moa_status", "phase",
    "library", "cell_line", "time_h", "dose_uM", "dose_unit",
    score_col, "effect_direction", p_col, padj_col, "rank",
    compound_col, "broad_id", "cmap_name", "canonical_smiles", "inchi_key"
  )
  tech_cols <- c(intersect(preferred, names(res2)), setdiff(names(res2), preferred))
  tech_view <- res2[, tech_cols, drop = FALSE]

  # Wetlab drug view
  drug_cols <- c(
    "pert_type", "pert_kind", "display_name", "pert_id",
    "compound_aliases", "target", "moa", "phase",
    score_col, "effect_direction", "cell_line", "time_h"
  )
  if (keep_dose_in_drug_view) {
    drug_cols <- c(drug_cols, "dose_uM")
  }
  drug_cols <- intersect(drug_cols, names(res2))
  wetlab_drug <- res2[res2$pert_kind == "Drug", drug_cols, drop = FALSE]
  wetlab_drug <- dedup_drug_rows(wetlab_drug, score_col = score_col)
  if (nrow(wetlab_drug) > 0) {
    wetlab_drug <- wetlab_drug[order(suppressWarnings(as.numeric(wetlab_drug[[score_col]]))), , drop = FALSE]
  }

  # Wetlab gene view
  gene_cols <- c(
    "pert_type", "pert_kind", "mode", "pert_name", "pert_id",
    score_col, "effect_direction", "cell_line", "time_h", "library"
  )
  gene_cols <- intersect(gene_cols, names(res2))
  wetlab_gene <- res2[res2$pert_kind == "Gene", gene_cols, drop = FALSE]
  if (nrow(wetlab_gene) > 0) {
    wetlab_gene <- wetlab_gene[order(suppressWarnings(as.numeric(wetlab_gene[[score_col]]))), , drop = FALSE]
  }

  # Drug context summary
  drug_only <- res2[res2$pert_kind == "Drug", , drop = FALSE]
  if (nrow(drug_only) > 0) {
    groups <- split(seq_len(nrow(drug_only)), as.character(drug_only$pert_id))
    idx_best <- unlist(lapply(groups, function(idx) {
      sc <- suppressWarnings(as.numeric(drug_only[[score_col]][idx]))
      idx[order(sc, na.last = TRUE)][1]
    }), use.names = FALSE)

    summary_cols <- intersect(c("pert_type", "pert_kind", "pert_id", "display_name", "phase", "moa_status", score_col),
                              names(drug_only))
    drug_summary <- drug_only[idx_best, summary_cols, drop = FALSE]
    names(drug_summary)[names(drug_summary) == score_col] <- "best_score"

    ctx_count <- table(as.character(drug_only$pert_id))
    drug_summary$n_contexts <- as.integer(ctx_count[as.character(drug_summary$pert_id)])

    if ("cell_line" %in% names(drug_only)) {
      n_cells <- tapply(drug_only$cell_line, drug_only$pert_id, function(x) length(unique(x[!is.na(x)])))
      drug_summary$n_cell_lines <- as.integer(n_cells[as.character(drug_summary$pert_id)])
    }
    if ("time_h" %in% names(drug_only)) {
      n_time <- tapply(drug_only$time_h, drug_only$pert_id, function(x) length(unique(x[!is.na(x)])))
      drug_summary$n_time <- as.integer(n_time[as.character(drug_summary$pert_id)])
    }

    ord_score <- suppressWarnings(as.numeric(drug_summary$best_score))
    drug_summary <- drug_summary[order(-drug_summary$n_contexts, ord_score), , drop = FALSE]
  } else {
    drug_summary <- data.frame()
  }

  # Write outputs if requested
  output_files <- character()
  if (isTRUE(write_outputs) || !is.null(output_dir)) {
    if (is.null(output_dir)) output_dir <- "."
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    fp1 <- file.path(output_dir, "wetlab_drug_view.tsv")
    fp2 <- file.path(output_dir, "wetlab_gene_view.tsv")
    fp3 <- file.path(output_dir, "tech_view_all.tsv")
    fp4 <- file.path(output_dir, "drug_context_summary.tsv")

    utils::write.table(wetlab_drug, fp1, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
    utils::write.table(wetlab_gene, fp2, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
    utils::write.table(tech_view, fp3, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
    utils::write.table(drug_summary, fp4, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

    output_files <- c(
      wetlab_drug_view = fp1,
      wetlab_gene_view = fp2,
      tech_view_all = fp3,
      drug_context_summary = fp4
    )

    if (verbose) {
      message("Wrote output files:")
      message(" - ", fp1)
      message(" - ", fp2)
      message(" - ", fp3)
      message(" - ", fp4)
    }
  }

  # Backward compatibility: legacy single output file
  if (!is.null(output_file)) {
    utils::write.csv(tech_view, file = output_file, row.names = FALSE)
    if (verbose) message("Wrote legacy output_file (tech_view_all as CSV): ", output_file)
  }

  return(list(
    wetlab_drug_view = wetlab_drug,
    wetlab_gene_view = wetlab_gene,
    tech_view_all = tech_view,
    drug_context_summary = drug_summary,
    output_files = output_files
  ))
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

  if (method == "split_colon") {
    return(vapply(strsplit(compound_strings, ":"), function(x) {
      if (length(x) >= part_index) x[part_index] else x[1]
    }, character(1)))
  } else if (method == "split_underscore") {
    return(vapply(strsplit(compound_strings, "_"), function(x) {
      if (length(x) >= part_index) x[part_index] else x[1]
    }, character(1)))
  } else if (method == "regex") {
    if (is.null(regex_pattern)) {
      stop("regex_pattern must be provided when method='regex'")
    }
    matches <- regmatches(compound_strings, regexpr(regex_pattern, compound_strings))
    return(ifelse(matches == "", compound_strings, matches))
  } else {
    stop("Invalid method. Must be 'split_colon', 'split_underscore', or 'regex'")
  }
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

  results <- data.frame(
    query_name = character(),
    matched_name = character(),
    similarity_score = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(query_names)) {
    query <- query_names[i]

    if (is.na(query) || query == "") {
      next
    }

    # Calculate similarities
    if (method == "levenshtein") {
      similarities <- RecordLinkage::levenshteinSim(query, reference_names) * 100
    } else if (method == "jaro") {
      similarities <- RecordLinkage::jarowinkler(query, reference_names) * 100
    } else if (method == "jarowinkler") {
      similarities <- RecordLinkage::jarowinkler(query, reference_names, W = 0.1) * 100
    } else {
      stop("Invalid method. Must be 'levenshtein', 'jaro', or 'jarowinkler'")
    }

    # Find top matches above threshold
    above_threshold <- which(similarities >= threshold)
    if (length(above_threshold) > 0) {
      sorted_indices <- above_threshold[order(similarities[above_threshold], decreasing = TRUE)]
      top_indices <- utils::head(sorted_indices, top_n)

      for (idx in top_indices) {
        results <- rbind(results, data.frame(
          query_name = query,
          matched_name = reference_names[idx],
          similarity_score = similarities[idx],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(results)
}
