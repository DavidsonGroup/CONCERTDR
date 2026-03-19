# Suppress R CMD check notes for ggplot2 NSE column names
utils::globalVariables(c("compound", "Score", "pValue"))

#' Create a Signature Data Frame from Gene Lists
#'
#' Converts separate lists of up-regulated and down-regulated genes into a
#' signature data frame compatible with \code{\link{process_signature_with_df}}
#' and other CONCERTDR functions.  This is a convenience wrapper for users who
#' only have gene lists (e.g.\ from a pathway database) rather than continuous
#' fold-change values.
#'
#' @param up_genes Character vector of up-regulated gene symbols.
#' @param down_genes Character vector of down-regulated gene symbols.
#' @param up_value Numeric value to assign to up-regulated genes (default: 1).
#' @param down_value Numeric value to assign to down-regulated genes (default: -1).
#'
#' @return A data frame with columns \code{Gene} and \code{log2FC}, suitable for
#'   use as the \code{signature_file} argument in
#'   \code{\link{process_signature_with_df}} and
#'   \code{\link{extract_signature_zscores}}.
#'
#' @examples
#' sig_df <- create_signature_from_gene_lists(
#'   up_genes   = c("TP53", "MYC", "BRCA1"),
#'   down_genes = c("EGFR", "VEGFA", "KRAS")
#' )
#' head(sig_df)
#'
#' @export
create_signature_from_gene_lists <- function(up_genes, down_genes,
                                             up_value = 1, down_value = -1) {
  if (missing(up_genes))   up_genes   <- character(0)
  if (missing(down_genes)) down_genes <- character(0)

  if (length(up_genes) == 0 && length(down_genes) == 0) {
    stop("At least one of 'up_genes' or 'down_genes' must be non-empty")
  }
  if (length(up_genes) == 0) {
    warning("No up-regulated genes provided. The signature will only contain down-regulated genes.")
  }
  if (length(down_genes) == 0) {
    warning("No down-regulated genes provided. The signature will only contain up-regulated genes.")
  }

  # Remove duplicates
  up_genes   <- unique(as.character(up_genes))
  down_genes <- unique(as.character(down_genes))

  # Check for overlap
  overlap <- intersect(up_genes, down_genes)
  if (length(overlap) > 0) {
    warning("Found ", length(overlap), " gene(s) in both up and down lists: ",
            paste(utils::head(overlap, 5), collapse = ", "),
            if (length(overlap) > 5) "..." else "",
            ". These genes will be removed from the down-regulated list.")
    down_genes <- setdiff(down_genes, overlap)
  }

  sig_df <- data.frame(
    Gene   = c(up_genes, down_genes),
    log2FC = c(rep(up_value, length(up_genes)),
               rep(down_value, length(down_genes))),
    stringsAsFactors = FALSE
  )

  # Sort by absolute log2FC descending
  sig_df <- sig_df[order(-abs(sig_df$log2FC)), , drop = FALSE]
  rownames(sig_df) <- NULL

  message(sprintf("Created signature with %d up-regulated and %d down-regulated genes",
                  length(up_genes), length(down_genes)))

  return(sig_df)
}

#' Process drug response signatures against reference data in a dataframe
#'
#' @param signature_file Path to a signature gene list file with log2FC values,
#'   or a data frame with columns \code{Gene} and \code{log2FC}.
#'   When a data frame is supplied the \code{read_method} argument is ignored.
#'   You can create a suitable data frame from gene lists using
#'   \code{\link{create_signature_from_gene_lists}}.
#' @param reference_df Dataframe containing reference data (combined across all conditions)
#' @param output_dir Directory for output files (default: "results")
#' @param permutations Number of permutations for statistical testing (default: 100)
#' @param methods Vector of method names to run (default: all)
#'        Options: "ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang"
#' @param topN Integer; number of top-ranked genes to use for XCos and XSum methods (default: 4)
#' @param read_method Character; method to use for reading signature file ("auto", "fread", or "read.table") (default: "auto")
#' @param save_files Logical; whether to save results to files (default: FALSE)
#'
#' @return A structured list containing:
#'   \item{results}{List containing results for each method}
#'   \item{summary}{Data frame of top hits across all methods}
#'   \item{gene_data}{Original signature gene data}
#'   \item{settings}{List of parameters used for the analysis}
#'   \item{common_genes}{Counts of genes found in reference data}
#'   Use print() or summary() methods to get a quick overview of results
#'
#' @examples
#' sig_file <- system.file("extdata", "example_signature.txt",
#'                         package = "CONCERTDR")
#' ref_file <- system.file("extdata", "example_reference_df.csv",
#'                         package = "CONCERTDR")
#' ref_df <- read.csv(ref_file, row.names = 1, check.names = FALSE)
#' ref_df$gene_symbol <- rownames(ref_df)
#' results <- process_signature_with_df(
#'   signature_file = sig_file,
#'   reference_df = ref_df,
#'   output_dir = tempdir(),
#'   permutations = 10,
#'   methods = "ks",
#'   save_files = FALSE
#' )
#' summary(results, top_n = 3)
#'
#' # You can also pass a data.frame directly:
#' sig_df <- read.delim(sig_file)
#' results2 <- process_signature_with_df(
#'   signature_file = sig_df,
#'   reference_df = ref_df,
#'   permutations = 10,
#'   methods = "ks",
#'   save_files = FALSE
#' )
#'
#' @export
process_signature_with_df <- function(signature_file, reference_df, output_dir = "results",
                                      permutations = 100, methods = c("ks", "xcos", "xsum", "gsea0",
                                                                      "gsea1", "gsea2", "zhang"),
                                      topN = 4, read_method = "auto", save_files = FALSE) {

  # Create output directory if it doesn't exist and files will be saved
  if (save_files && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Record start time for performance tracking
  start_time <- Sys.time()
  
  # Handle signature input: data.frame or file path
  if (is.data.frame(signature_file)) {
    gene_data <- signature_file
    signature_label <- "(data.frame)"
    message("Using signature data from data.frame")
  } else if (is.character(signature_file) && length(signature_file) == 1L) {
    signature_label <- signature_file
    message("Reading signature data from ", signature_file)
    if (!file.exists(signature_file)) {
      stop("Signature file not found: ", signature_file)
    }
    tryCatch({
      if (read_method == "auto") {
        # Auto-detect best method
        if (requireNamespace("data.table", quietly = TRUE)) {
          gene_data <- data.table::fread(signature_file, header = TRUE,
                                         stringsAsFactors = FALSE, data.table = FALSE)
        } else {
          gene_data <- utils::read.delim(signature_file, header = TRUE, sep = "\t",
                                         stringsAsFactors = FALSE, comment.char = "")
        }
      } else if (read_method == "fread" && requireNamespace("data.table", quietly = TRUE)) {
        gene_data <- data.table::fread(signature_file, header = TRUE,
                                       stringsAsFactors = FALSE, data.table = FALSE)
      } else {
        # Default to read.table
        gene_data <- utils::read.delim(signature_file, header = TRUE, sep = "\t",
                                       stringsAsFactors = FALSE, comment.char = "")
      }
    }, error = function(e) {
      stop("Error reading signature file: ", e$message)
    })
  } else {
    stop("signature_file must be either a file path (character) or a data.frame with 'Gene' and 'log2FC' columns.")
  }
  
  # Check if required columns exist
  if (!all(c("Gene", "log2FC") %in% colnames(gene_data))) {
    stop("Signature file must contain 'Gene' and 'log2FC' columns")
  }
  
  # Separate up and down regulated genes based on log2FC values
  Up <- gene_data$Gene[gene_data$log2FC > 0]
  Down <- gene_data$Gene[gene_data$log2FC < 0]
  
  if (length(Up) == 0 || length(Down) == 0) {
    stop("Signature file must contain both up-regulated (log2FC > 0) and down-regulated (log2FC < 0) genes")
  }
  
  # Define query for XCos (signed values with gene names)
  query_genes <- gene_data$Gene
  query_values <- suppressWarnings(as.numeric(gene_data$log2FC))
  query <- stats::setNames(query_values, query_genes)
  
  # Prepare reference data for processing
  message("Preparing reference data for analysis...")
  
  # Check if reference_df has gene_symbol column
  if (!("gene_symbol" %in% colnames(reference_df))) {
    stop("Reference dataframe must contain a 'gene_symbol' column")
  }
  
  # Create reference matrix with genes as row names
  ref_data <- reference_df
  rownames(ref_data) <- ref_data$gene_symbol
  ref_data$gene_symbol <- NULL

  # Convert to matrix format
  ref <- as.matrix(ref_data)
  
  # Verify genes exist in reference
  common_up <- intersect(Up, rownames(ref))
  common_down <- intersect(Down, rownames(ref))
  
  # Calculate percent overlap to assess gene coverage
  pct_up <- round(length(common_up) / length(Up) * 100, 1)
  pct_down <- round(length(common_down) / length(Down) * 100, 1)
  
  message(sprintf("Found %d/%d up-regulated genes (%g%%) and %d/%d down-regulated genes (%g%%) in reference",
                  length(common_up), length(Up), pct_up,
                  length(common_down), length(Down), pct_down))
  
  if (length(common_up) == 0 || length(common_down) == 0) {
    stop("No matching genes found in reference data. Please check your signature genes.")
  }
  
  # Initialize results list
  all_results <- list()
  
  # Store metadata about the analysis
  settings <- list(
    signature_file = signature_label,
    permutations = permutations,
    methods = methods,
    topN = topN,
    time_started = start_time,
    reference_dimensions = dim(ref),
    reference_genes = length(rownames(ref))
  )
  
  # Available methods
  all_methods <- list(
    ks = function() {
      message("Running KS score...")
      score_ks(refMatrix = ref, queryUp = common_up, queryDown = common_down,
               permuteNum = permutations)
    },
    
    xcos = function() {
      message(sprintf("Running XCos score with topN = %d...", topN))
      score_xcos(refMatrix = ref, query = query[names(query) %in% rownames(ref)],
                 topN = topN, permuteNum = permutations)
    },
    
    xsum = function() {
      message(sprintf("Running XSum score with topN = %d...", topN))
      score_xsum(refMatrix = ref, queryUp = common_up, queryDown = common_down,
                 topN = topN, permuteNum = permutations)
    },
    
    gsea0 = function() {
      message("Running GSEA weight 0 score...")
      score_gsea0(refMatrix = ref, queryUp = common_up, queryDown = common_down,
                  permuteNum = permutations)
    },
    
    gsea1 = function() {
      message("Running GSEA weight 1 score...")
      score_gsea1(refMatrix = ref, queryUp = common_up, queryDown = common_down,
                  permuteNum = permutations)
    },
    
    gsea2 = function() {
      message("Running GSEA weight 2 score...")
      score_gsea2(refMatrix = ref, queryUp = common_up, queryDown = common_down,
                  permuteNum = permutations)
    },
    
    zhang = function() {
      message("Running Zhang score...")
      score_zhang(refMatrix = ref, queryUp = common_up, queryDown = common_down,
                  permuteNum = permutations)
    }
  )
  
  # Validate selected methods
  invalid_methods <- setdiff(methods, names(all_methods))
  if (length(invalid_methods) > 0) {
    message("Invalid methods specified: ", paste(invalid_methods, collapse = ", "),
            ". Will be ignored.")
    methods <- intersect(methods, names(all_methods))
  }
  
  if (length(methods) == 0) {
    stop("No valid methods to run. Valid methods are: ", paste(names(all_methods), collapse = ", "))
  }
  
  # Run selected scoring methods
  for (method in methods) {
    tryCatch({
      result_df <- all_methods[[method]]()
      
      # Add compound names if not already included
      if (!("compound" %in% colnames(result_df))) {
        result_df <- cbind(compound = rownames(result_df), result_df)
      }
      
      # Add rank column for easier interpretation
      result_df$rank <- rank(-result_df$Score)
      
      # Store in results list
      all_results[[method]] <- result_df
      
      # Save individual results if requested
      if (save_files) {
        output_file <- file.path(output_dir, paste0("sig_match_", method, "_results.csv"))
        utils::write.csv(result_df, file = output_file, row.names = FALSE)
        message("Saved ", method, " results to ", output_file)
      }
    }, error = function(e) {
      warning("Error running ", method, " method: ", e$message)
      all_results[[method]] <- data.frame(
        compound = character(0),
        Score = numeric(0),
        pValue = numeric(0),
        error = character(0)
      )
      all_results[[method]]$error <- e$message
    })
  }
  
  # Create summary of top hits across all methods
  summary_df <- create_summary_from_results(all_results)
  
  # Save summary file if requested
  if (save_files) {
    summary_file <- file.path(output_dir, "summary_results.csv")
    utils::write.csv(summary_df, file = summary_file, row.names = FALSE)
    message("Saved summary report to ", summary_file)
  }
  
  # Record end time
  end_time <- Sys.time()
  time_taken <- difftime(end_time, start_time, units = "mins")
  
  # Update settings with completion info
  settings$time_completed <- end_time
  settings$time_taken_mins <- as.numeric(time_taken)
  
  # Create a structured result object
  result_object <- structure(
    list(
      results = all_results,
      summary = summary_df,
      gene_data = gene_data,
      settings = settings,
      common_genes = list(
        up = list(found = common_up, count = length(common_up), percent = pct_up),
        down = list(found = common_down, count = length(common_down), percent = pct_down)
      )
    ),
    class = "cmap_signature_result"
  )
  
  # Return the complete result object
  return(result_object)
}

#' Create a summary from signature matching results
#'
#' @param results_list List of results from different methods
#' @param top_n Number of top hits to include (default: 20)
#'
#' @return Data frame with summary of top hits across methods
#'
#' @keywords internal
create_summary_from_results <- function(results_list, top_n = 20) {
  # Initialize summary dataframe
  summary_df <- data.frame(
    compound = character(),
    method = character(),
    Score = numeric(),
    pValue = numeric(),
    rank = integer(),
    stringsAsFactors = FALSE
  )
  
  # Extract top hits from each method
  for (method_name in names(results_list)) {
    result <- results_list[[method_name]]
    
    # Skip if error or empty
    if (nrow(result) == 0 || "error" %in% colnames(result)) {
      next
    }
    
    # Make sure rank column exists
    if (!"rank" %in% colnames(result)) {
      result$rank <- rank(-result$Score)
    }
    
    # Get top hits
    top_hits <- result[result$rank <= top_n, ]
    
    if (nrow(top_hits) > 0) {
      top_hits$method <- method_name
      summary_df <- rbind(
        summary_df,
        top_hits[, c("compound", "method", "Score", "pValue", "rank")]
      )
    }
  }
  
  # Sort by method and rank
  summary_df <- summary_df[order(summary_df$method, summary_df$rank), ]
  
  # Add global rank across methods
  if (nrow(summary_df) > 0) {
    # Calculate a weighted score considering both Score and p-value
    summary_df$weighted_score <- summary_df$Score * (1 - summary_df$pValue)
    
    # Rank across all methods
    summary_df$global_rank <- rank(-summary_df$weighted_score)
    
    # Remove the temporary weighted score column
    summary_df$weighted_score <- NULL
  }
  
  return(summary_df)
}

#' Print method for cmap_signature_result objects
#'
#' @param x A cmap_signature_result object
#' @param ... Additional arguments (not used)
#'
#' @return x invisibly
#'
#' @examples
#' mock_res <- structure(
#'   list(
#'     results = list(
#'       ks = data.frame(compound = c("imatinib", "dasatinib"),
#'                       Score = c(-0.72, -0.45), pValue = c(0.002, 0.041),
#'                       pAdjValue = c(0.020, 0.205), stringsAsFactors = FALSE)
#'     ),
#'     summary = data.frame(compound = "imatinib", method = "ks",
#'                          Score = -0.72, pValue = 0.002, global_rank = 1L,
#'                          stringsAsFactors = FALSE),
#'     gene_data    = data.frame(Gene = c("TP53", "MYC"),
#'                               log2FC = c(1.5, -2.1), stringsAsFactors = FALSE),
#'     settings     = list(signature_file = "ex.txt",
#'                         time_completed = Sys.time(), time_taken_mins = 0.1,
#'                         methods = "ks", permutations = 10),
#'     common_genes = list(up   = list(found = "TP53", count = 1L, percent = 50),
#'                         down = list(found = "MYC",  count = 1L, percent = 50))
#'   ), class = "cmap_signature_result"
#' )
#' print(mock_res)
#'
#' @export
print.cmap_signature_result <- function(x, ...) {
  cat("CMap Signature Matching Results\n")
  cat("===============================\n\n")
  
  # Print basic information
  cat(sprintf("Signature file: %s\n", x$settings$signature_file))
  cat(sprintf("Analysis completed: %s\n", format(x$settings$time_completed)))
  cat(sprintf("Time taken: %.2f minutes\n\n", x$settings$time_taken_mins))
  
  # Print methods used
  cat("Methods used:\n")
  for (method in names(x$results)) {
    result <- x$results[[method]]
    if ("error" %in% colnames(result)) {
      cat(sprintf("  - %s: ERROR - %s\n", method, result$error[1]))
    } else {
      cat(sprintf("  - %s: %d results\n", method, nrow(result)))
    }
  }
  
  # Print gene coverage
  cat("\nGene coverage:\n")
  cat(sprintf("  Up-regulated: %d/%d genes (%.1f%%)\n",
              x$common_genes$up$count, length(x$gene_data$Gene[x$gene_data$log2FC > 0]),
              x$common_genes$up$percent))
  cat(sprintf("  Down-regulated: %d/%d genes (%.1f%%)\n",
              x$common_genes$down$count, length(x$gene_data$Gene[x$gene_data$log2FC < 0]),
              x$common_genes$down$percent))
  
  # Print top compounds from summary (if available)
  if (nrow(x$summary) > 0) {
    top_n <- min(5, nrow(x$summary))
    top_compounds <- x$summary[order(x$summary$global_rank)[seq_len(top_n)], ]
    
    cat("\nTop compounds across all methods:\n")
    for (i in seq_len(nrow(top_compounds))) {
      cat(sprintf("  %d. %s (method: %s, score: %.4f, p-value: %.4f)\n",
                  i, top_compounds$compound[i], top_compounds$method[i],
                  top_compounds$Score[i], top_compounds$pValue[i]))
    }
  }
  
  cat("\nUse the following to explore the results:\n")
  cat("  - $results: List of result data frames for each method\n")
  cat("  - $summary: Summary of top hits across all methods\n")
  cat("  - $gene_data: Original signature gene data\n")
  cat("  - $settings: Analysis settings and metadata\n")
  cat("  - $common_genes: Genes found in the reference data\n")
  cat("\nExample: result_obj$results$ks to view KS score results\n")
  
  invisible(x)
}

#' Summary method for cmap_signature_result objects
#'
#' @param object A cmap_signature_result object
#' @param top_n Number of top compounds to show (default: 10)
#' @param ... Additional arguments (not used)
#'
#' @return A data frame with the top compounds
#'
#' @examples
#' mock_res <- structure(
#'   list(
#'     results = list(
#'       ks = data.frame(compound = c("imatinib", "dasatinib"),
#'                       Score = c(-0.72, -0.45), pValue = c(0.002, 0.041),
#'                       pAdjValue = c(0.020, 0.205), stringsAsFactors = FALSE)
#'     ),
#'     summary = data.frame(compound = c("imatinib", "dasatinib"),
#'                          method = c("ks", "ks"),
#'                          Score = c(-0.72, -0.45), pValue = c(0.002, 0.041),
#'                          global_rank = 1:2, stringsAsFactors = FALSE),
#'     gene_data    = data.frame(Gene = c("TP53", "MYC"),
#'                               log2FC = c(1.5, -2.1), stringsAsFactors = FALSE),
#'     settings     = list(signature_file = "ex.txt",
#'                         time_completed = Sys.time(), time_taken_mins = 0.1,
#'                         methods = "ks", permutations = 10),
#'     common_genes = list(up   = list(found = "TP53", count = 1L, percent = 50),
#'                         down = list(found = "MYC",  count = 1L, percent = 50))
#'   ), class = "cmap_signature_result"
#' )
#' summary(mock_res, top_n = 2)
#'
#' @export
summary.cmap_signature_result <- function(object, top_n = 10, ...) {
  if (nrow(object$summary) == 0) {
    message("No summary data available.")
    return(NULL)
  }
  
  # Get top compounds across all methods
  top_compounds <- object$summary[order(object$summary$global_rank), ]
  
  # Limit to requested number
  top_compounds <- utils::head(top_compounds, top_n)
  
  # Print a summary
    summary_tbl <- top_compounds[, c("compound", "method", "Score", "pValue", "global_rank")]
    message("Top ", top_n, " compounds across all methods:\n",
      paste(utils::capture.output(summary_tbl), collapse = "\n"))

  # Return the data invisibly
  invisible(top_compounds)
}

#' Plot method for cmap_signature_result objects
#'
#' @param x A cmap_signature_result object
#' @param method Which method to plot (default: first method in results)
#' @param plot_type Type of plot: "scores", "volcano", or "heatmap" (default: "scores")
#' @param top_n Number of compounds to include (default: 20)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return A ggplot object
#'
#' @examples
#' mock_res <- structure(
#'   list(
#'     results = list(
#'       ks = data.frame(compound = c("imatinib", "dasatinib", "doxorubicin"),
#'                       Score = c(-0.72, -0.45, 0.31),
#'                       pValue = c(0.002, 0.041, 0.280),
#'                       pAdjValue = c(0.020, 0.205, 0.560),
#'                       stringsAsFactors = FALSE)
#'     ),
#'     summary = data.frame(compound = "imatinib", method = "ks",
#'                          Score = -0.72, pValue = 0.002, global_rank = 1L,
#'                          stringsAsFactors = FALSE),
#'     gene_data    = data.frame(Gene = c("TP53", "MYC"),
#'                               log2FC = c(1.5, -2.1), stringsAsFactors = FALSE),
#'     settings     = list(signature_file = "ex.txt",
#'                         time_completed = Sys.time(), time_taken_mins = 0.1,
#'                         methods = "ks", permutations = 10),
#'     common_genes = list(up   = list(found = "TP53", count = 1L, percent = 50),
#'                         down = list(found = "MYC",  count = 1L, percent = 50))
#'   ), class = "cmap_signature_result"
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(mock_res, method = "ks", plot_type = "scores")
#' }
#'
#' @export
plot.cmap_signature_result <- function(x, method = NULL, plot_type = "scores", top_n = 20, ...) {
  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it with: install.packages('ggplot2')")
  }
  
  # If method is not specified, use the first available method
  if (is.null(method)) {
    if (length(x$results) == 0) {
      stop("No results available to plot.")
    }
    method <- names(x$results)[1]
  }
  
  # Check if the specified method exists
  if (!method %in% names(x$results)) {
    stop("Method '", method, "' not found in results. Available methods: ",
         paste(names(x$results), collapse = ", "))
  }
  
  # Get the data for the specified method
  result_df <- x$results[[method]]
  
  # Check if there was an error with this method
  if ("error" %in% colnames(result_df)) {
    stop("Cannot plot results for method '", method, "' due to error: ", result_df$error[1])
  }
  
  # Different plot types
  if (plot_type == "scores") {
    # Bar plot of top scores
    top_df <- result_df[order(-result_df$Score), ][seq_len(min(top_n, nrow(result_df))), ]
    
    p <- ggplot2::ggplot(top_df, ggplot2::aes(x = stats::reorder(compound, Score), y = Score)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste("Top", top_n, "Compounds by", method, "Score"),
        x = "Compound",
        y = paste(method, "Score")
      ) +
      ggplot2::theme_minimal()
    
  } else if (plot_type == "volcano") {
    # Volcano plot of scores vs p-values
    p <- ggplot2::ggplot(result_df, ggplot2::aes(x = Score, y = -log10(pValue))) +
      ggplot2::geom_point(ggplot2::aes(color = abs(Score) > 0.5 & pValue < 0.05)) +
      ggplot2::scale_color_manual(values = c("grey", "red"),
                                  name = "Significant",
                                  labels = c("No", "Yes")) +
      ggplot2::labs(
        title = paste("Volcano Plot of", method, "Results"),
        x = paste(method, "Score"),
        y = "-log10(p-value)"
      ) +
      ggplot2::theme_minimal()
    
  } else if (plot_type == "heatmap") {
    stop(
      "plot_type='heatmap' is not currently supported for cmap_signature_result objects ",
      "because the object does not store expression matrices required for heatmap rendering."
    )
    
  } else {
    stop("Invalid plot_type: '", plot_type, "'. Must be one of: 'scores', 'volcano', 'heatmap'")
  }
  
  return(p)
}
#' Complete workflow from configuration to signature matching
#'
#' @param config_file Path to configuration file with selected parameters
#' @param signature_file Path to signature gene list with log2FC values,
#'   or a data frame with \code{Gene} and \code{log2FC} columns.
#' @param geneinfo_file Path to the gene info file
#' @param siginfo_file Path to the signature info file
#' @param gctx_file Path to the GCTX file
#' @param output_dir Directory for output files (default: "results")
#' @param methods Vector of method names to run (default: all)
#'        Options: "ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang"
#' @param topN Integer; number of top-ranked genes to use for XCos and XSum methods (default: 4)
#' @param permutations Number of permutations for statistical testing (default: 100)
#' @param save_files Logical; whether to write per-method and summary result
#'   files to \code{output_dir} (default: FALSE)
#' @param keep_all_genes Logical; whether to keep all genes when extracting CMap data (default: TRUE)
#' @param read_method Character; method to use for reading signature file ("auto", "fread", or "read.table") (default: "auto")
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return List containing results from all methods
#'