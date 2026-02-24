#' @title Demonstrate CONCERTDR Workflow
#' @description
#' This function demonstrates the full workflow of the CONCERTDR package
#' with example data. Users can run this to see how the package works
#' in practice. It extracts a small subset of data from a GCTX file
#' to create a realistic reference for signature matching.
#'
#' @param demo_dir Directory to store demonstration files
#' @param use_minimal Logical; use minimal settings for faster demo
#' @param gctx_file Path to GCTX file (if available)
#' @param geneinfo_file Path to geneinfo file (if available)
#' @param siginfo_file Path to siginfo file (if available)
#' @export
demonstrate_workflow <- function(demo_dir = "CONCERTDR_demo",
                                 use_minimal = TRUE,
                                 gctx_file = NULL,
                                 geneinfo_file = NULL,
                                 siginfo_file = NULL) {
  # Create demo directory if it doesn't exist
  if (!dir.exists(demo_dir)) {
    dir.create(demo_dir, recursive = TRUE)
    message("Created demonstration directory: ", demo_dir)
  }

  # Create subdirectories
  conf_dir <- file.path(demo_dir, "conf")
  output_dir <- file.path(demo_dir, "output")
  results_dir <- file.path(demo_dir, "results")

  for (dir in c(conf_dir, output_dir, results_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
  }


  message("\n--- STEP 1: Creating Configuration Template ---")
  # Create a template config file
  config_file <- create_cmap_config_template(
    dest_dir = conf_dir,
    template_name = "cmap_options_template.conf",
    overwrite = TRUE
  )

  message("\n--- STEP 3: Creating Example Signature ---")
  # Load bundled example signature (or generate a fallback signature)
  example_sig_file <- file.path(demo_dir, "signature.txt")

  bundled_sig <- system.file("extdata", "example_signature.txt", package = "CONCERTDR")
  if (nzchar(bundled_sig) && file.exists(bundled_sig)) {
    file.copy(bundled_sig, example_sig_file, overwrite = TRUE)
  } else {
    fallback_sig <- data.frame(
      Gene = c("STAT3", "TP53", "MYC", "BRCA1", "EGFR", "BCL2"),
      log2FC = c(2.5, -1.8, 3.2, -2.1, 1.4, -1.2),
      stringsAsFactors = FALSE
    )
    utils::write.table(fallback_sig, example_sig_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  message("\n--- STEP 4: Creating Reference Data ---")
  ref_file <- file.path(demo_dir, "reference_data.csv")

  resolve_input_file <- function(user_path, extdata_name) {
    if (!is.null(user_path) && nzchar(user_path) && file.exists(user_path)) {
      return(user_path)
    }

    bundled <- system.file("extdata", extdata_name, package = "CONCERTDR")
    if (nzchar(bundled) && file.exists(bundled)) {
      return(bundled)
    }

    return(NULL)
  }

  gctx_path <- resolve_input_file(gctx_file, "level5_beta_trt_cp_n720216x12328.gctx")
  geneinfo_path <- resolve_input_file(geneinfo_file, "geneinfo_beta.txt")
  siginfo_path <- resolve_input_file(siginfo_file, "siginfo_beta.txt")

  # Check if we have real CMap data files
  have_real_data <- !is.null(gctx_path) && !is.null(geneinfo_path) && !is.null(siginfo_path)

  ref_df <- NULL
  using_synthetic_reference <- FALSE

  if (have_real_data) {
    message("Using real CMap data files to extract reference data")

    ref_df <- tryCatch({
      ref_df <- extract_cmap_data_from_config(
        config_file = config_file,
        geneinfo_file = geneinfo_path,
        siginfo_file = siginfo_path,
        gctx_file = gctx_path
      )
      ref_df
    }, error = function(e) {
      message("Failed to extract reference data from real files: ", e$message)
      NULL
    })
  }

  if (is.null(ref_df) || !is.data.frame(ref_df) || nrow(ref_df) == 0 || ncol(ref_df) <= 1) {
    message("Falling back to synthetic reference matrix for demonstration.")

    sig_df <- utils::read.delim(example_sig_file, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("Gene", "log2FC") %in% names(sig_df))) {
      stop("Example signature file must contain 'Gene' and 'log2FC' columns")
    }

    genes <- unique(sig_df$Gene)
    n_samples <- if (use_minimal) 30 else 100
    set.seed(42)
    synth_mat <- matrix(stats::rnorm(length(genes) * n_samples),
                        nrow = length(genes), ncol = n_samples)
    colnames(synth_mat) <- paste0("SYNTH_", seq_len(n_samples))

    ref_df <- data.frame(gene_symbol = genes, synth_mat,
                         stringsAsFactors = FALSE, check.names = FALSE)
    utils::write.csv(ref_df, ref_file, row.names = FALSE)
    using_synthetic_reference <- TRUE
  }

  message("\n--- STEP 5: Signature Matching ---")
  # Run signature matching with the reference data
  methods_to_use <- if (use_minimal) c("ks", "xsum") else c("ks", "xcos", "xsum", "zhang")
  perm_num <- if (use_minimal) 10 else 100

  result <- tryCatch({
    process_signature_with_df(
      signature_file = example_sig_file,
      reference_df = ref_df,
      output_dir = results_dir,
      permutations = perm_num,
      methods = methods_to_use
    )
  }, error = function(e) {
    e
  })

  if (inherits(result, "error")) {
    message("Signature matching could not be completed: ", result$message)
    message("Install dependency 'RCSM' first (for example via remotes::install_github('Jasonlinchina/RCSM')).")
  } else {
    message("Successfully ran signature matching with methods: ", paste(methods_to_use, collapse = ", "))
    message("Results saved to: ", file.path(results_dir))
  }

  message("\n--- WORKFLOW DEMONSTRATION COMPLETE ---")
  message("Demonstration files are in: ", demo_dir)

  invisible(list(
    config_file = config_file,
    signature_file = example_sig_file,
    reference_file = ref_file,
    results_dir = results_dir,
    using_synthetic_reference = using_synthetic_reference,
    result = if (inherits(result, "error")) NULL else result
  ))
}
