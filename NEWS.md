# CONCERTDR 0.99.1

## Bug fixes and improvements

* **`process_signature_with_df()` — gene coverage message now reflects topN**:
  The "Found X/Y up-regulated genes" progress message previously reported
  overlap against the full signature gene list. It now sorts each direction by
  effect size, truncates to the `topN` genes actually passed to the scoring
  functions, and reports overlap against that subset. The denominator and
  percentage therefore match the genes that scoring methods such as XSum and
  XCos will use. The message also appends `(top N per direction)` to make the
  basis of the calculation explicit.

## New features

* **Split-direction mode for `extract_signature_zscores()` and
  `plot_signature_direction_tile_barcode()`**: Both functions now accept a
  `split_direction` argument (default `FALSE`).

  - In `extract_signature_zscores()`, setting `split_direction = TRUE` applies
    `max_genes` independently to the up-regulated (log2FC > 0) and
    down-regulated (log2FC < 0) gene sets, so each direction contributes up to
    `max_genes` columns (up to `2 * max_genes` total). Genes within each group
    are ordered by effect size (most-positive first for up, most-negative first
    for down) before truncation.

  - In `plot_signature_direction_tile_barcode()`, setting `split_direction =
    TRUE` renders the heatmap as two side-by-side panels: up-regulated genes on
    the left and down-regulated genes on the right. A companion `gap_width`
    argument (default `5` mm) controls the gap between panels. Row clustering is
    performed on the full combined matrix so both panels display perturbations in
    the same order, enabling direct visual comparison of activation and
    suppression patterns within a single figure.

# CONCERTDR 0.99.0

## Bioconductor submission

* Fixed mailing list registration

# CONCERTDR 0.5.0

## New features

* **`create_signature_from_gene_lists()`**: New convenience function to convert
  separate up-regulated and down-regulated gene lists into a signature data
  frame compatible with all CONCERTDR scoring functions.  Up-regulated genes
  receive a default value of +1 and down-regulated genes receive -1.

* **Data frame input for signatures**: `process_signature_with_df()`,
  `extract_signature_zscores()`, and `plot_signature_direction_tile_barcode()`
  now accept a `data.frame` (with `Gene` and `log2FC` columns) for the
  `signature_file` argument, in addition to the existing file-path input.
  This allows users to pass in-memory signature objects throughout the
  entire workflow without writing intermediate files.

* **Progress message for GCTX extraction**: `extract_cmap_data_from_siginfo()`
  now prints a prominent notice before reading the GCTX file, letting users
  know the step may take a long time and should not be interrupted.

## Documentation

* Reorganised the workflow steps in all vignettes:
  - Steps 2 (Build reference matrix) and 3 (Prepare query signature) have
    been swapped so that users prepare their signature first.
  - Steps 4 (Annotate) and 5 (Visualise) have been merged into a single step.
* Added examples for data-frame input and gene-list conversion in
  `introduction.Rmd` and `signature_matching.Rmd`.

# CONCERTDR 0.4.0

## New features

* Added `extract_signature_zscores()` for pre-computing z-score matrices
  without rendering a plot, enabling reuse across repeated calls.
* Added `plot_signature_direction_tile_barcode()` for publication-quality
  directional barcode heatmaps via ComplexHeatmap.
* Added `annotate_drug_results()` producing four analysis-ready views:
  wetlab drug view, wetlab gene view, full technical table, and drug context
  summary.
* Added `subset_siginfo_beta()` for interactive and non-interactive filtering
  of the CMap siginfo file.

## Bug fixes and improvements

* Replaced `installed.packages()` in `.onLoad` with `requireNamespace()` to
  comply with CRAN/Bioconductor startup function guidelines.
* Removed non-ASCII characters from R source files.
* Replaced `sapply()` with `vapply()` throughout for type-safe output.
* Replaced `1:nrow(x)` idiom with `seq_len(nrow(x))` throughout.
* Replaced `\dontrun{}` with `\donttest{}` in all man page examples.
