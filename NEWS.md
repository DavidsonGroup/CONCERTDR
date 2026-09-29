# CONCERTDR 0.99.4

## Performance

* The seven permutation methods (`ks`, `xcos`, `xsum`, `gsea0`, `gsea1`,
  `gsea2`, `zhang`) now score blocks of profiles at once with `matrixStats`
  instead of one profile at a time, controlled by the new `vectorized`
  argument (default `TRUE`) of `process_signature_with_df()` and the
  `score_*()` functions. Formulas are unchanged. On a 978-gene and a
  12328-gene reference with 100 permutations the seven methods together ran
  22 and 35 times faster, with the largest gains for KS, GSEA and Zhang.
* Scores equal the row-by-row values. For `xcos`, `xsum`, `gsea1`, `gsea2` and
  `zhang`, p-values equal the row-by-row values for the same seed (up to one
  permutation where a permuted score equals the observed one up to rounding),
  because the permuted gene sets are drawn in the same order.
* `ks` and `gsea0` scores depend on a gene set only through its ranks, so the
  permutation null is the same for every profile and is now computed once and
  shared. The p-value distribution is the same, but the random numbers differ,
  so p-values change for a fixed seed. P-values are still `b / B`.
* `gsea0` now resolves a running sum whose maximum exactly equals its absolute
  minimum (about 0.3% of random gene sets) deterministically, in favour of the
  minimum. The row-by-row code decided this tie by floating-point rounding.
* Reference matrices with missing values or duplicated gene names use the
  row-by-row code automatically; `vectorized = FALSE` selects it explicitly.
  The block size can be set with `options(CONCERTDR.chunk_size = )`
  (default 2000 profiles).
* New import: `matrixStats`.

## Breaking changes

* New `direction` argument (`"reversal"` by default, or `"mimic"`) in
  `process_signature_with_df()`, `annotate_drug_results()`,
  `extract_signature_zscores()` and `plot_signature_direction_tile_barcode()`.
  It only orders results: `"reversal"` ranks the lowest scores first,
  `"mimic"` the highest first. No result is dropped by the sign of its score,
  and p-values are unchanged. The default now ranks the lowest scores first
  everywhere; the cross-method summary previously put the highest scores
  first.
* Cross-method `global_rank` now ranks each hit's rank percentile among all
  profiles of its method (`rank_percentile`), avoiding comparisons of
  incompatible raw score scales. It is a relative ordering of method-specific
  hits, not a combined p-value.
* Wetlab drug rows now retain the complete best-scoring context rather than
  combining its score with another context's dose or time.
* `annotate_drug_results()` no longer has the `drug_info_file`,
  `fuzzy_threshold` and `perfect_match_only` arguments; they were deprecated
  and ignored.
* `extract_cmap_data_from_siginfo()` no longer has the `keep_all_genes`
  argument, which was never used.
* `process_signature_with_df()` now validates `read_method` with `match.arg()`
  instead of silently falling back to `read.delim()` for unknown values.
* File-reading errors now name the offending argument (for example
  `geneinfo_file not found: ...`) consistently across the package.

## Bug fixes

* Fix GitHub issue #3: duplicate gene symbols in geneinfo (including LINCS2020
  MIA2) are resolved by preferring `landmark > best inferred > inferred > other`,
  with input order breaking ties. A warning names the discarded and retained
  gene IDs. Extraction, direct heatmap lookup and CamSum library rho estimation
  use the same rule, so `landmark = FALSE` can use the full feature space.
* GCTX reads now use sorted, unique indices on both axes and explicitly restore
  the requested order, including repeats, before assigning labels.
* Compound annotation retains all distinct nonblank target and MOA values
  across rows of the same `pert_id`, joined with `; `, instead of keeping only
  the first row. Signature row counts and context counts are preserved.
* CamSum metadata is read only when CamSum runs, inside its method-level
  error handler, so malformed CamSum configuration cannot stop other methods.
* The shared table reader selects columns during `fread()` rather than
  allocating the entire metadata table first.
* GCTX extraction maps gene symbols to actual returned gene IDs and handles
  absent genes/signatures. Empty extractions produce explicit errors.
* Regex compound-ID extraction preserves unmatched entries and vector length.
* Direct barcode plotting forwards `split_direction` to extraction; one-gene
  matrices preserve dimensions and single-profile plots skip row clustering.
* Query gene coverage now consistently reports the full query overlap:
  `topN` in XCos/XSum limits reference-profile extremes, not query genes.

* `process_signature_with_df()`: when a scoring method fails, its error is now
  recorded in `results` (a one-row data frame with an `error` column, which
  `print()` and `plot()` already expected). Before, the error handler only
  modified a local copy, so the method silently disappeared from the results.
* `process_combinations_file()` could not read tab-separated combination files
  whose values contain spaces (such as `24 h`), including the file built in its
  own example; it now reads them with an explicit tab separator.
* `subset_siginfo_beta()`: `preview N` on the last filter column no longer
  lists that column's own options as "downstream" options.

## Internal changes

* Scoring methods share one permutation framework (`.permutation_test()`);
  `score_gsea0/1/2` are one weighted implementation, with results and random
  number streams identical to before.
* GCTX reading is shared between `fast_parse_gctx()` and CamSum
  (`.gctx_ids()`, `.gctx_read_matrix()`).
* `plot_signature_direction_tile_barcode()` builds every panel with a single
  `.barcode_panel()` helper, and
  `extract_signature_zscores()`, `annotate_drug_results()` and
  `subset_siginfo_beta()` are split into small documented helpers.
* New `R/utils.R` with the shared table reader `.read_cmap_table()`; removed
  the `requireNamespace("data.table")` fallbacks (`data.table` is an Import),
  the `.onLoad()` dependency check, the empty `R/cmap_download.R`, an orphaned
  roxygen block in `R/signature_matching.R` and a committed build tarball.

# CONCERTDR 0.99.3

## New features

* **New matching method `camsum`** (CamSum). Scores each reference profile as
  `T = (sum_U z - sum_D z) / (s_i * sqrt(k * VIF))`, where `s_i` is the
  profile's standard deviation over all genes and
  `VIF = max(1, 1 + (k - 1) * rho_bar)`; `rho_bar` is the mean inter-gene
  correlation of the query across the reference library (down genes negated).
  p-values are analytic (standard normal), so no permutations are needed.
  Run it with `process_signature_with_df(..., methods = "camsum")`; it is not
  part of the default method set, so existing calls give the same results.
  New arguments `camsum_alternative` ("two.sided", "greater", "less") and
  `camsum_rho_bar`, and the new exported helper `compute_camsum_rho()`.
* CamSum estimates `rho_bar` on the whole library of each profile's
  perturbation type, not on the (possibly filtered) reference passed in: when
  options `CONCERTDR.gctx_file`, `CONCERTDR.siginfo_file` and
  `CONCERTDR.geneinfo_file` are set, every profile of that type is streamed
  from the GCTX file (two passes; the value is cached for the session).
  Profiles of different types get their own `rho_bar`. Without these options
  `rho_bar` is estimated on the reference with a warning. Where it came from is
  recorded in `settings$camsum$rho_source`.

## Other changes

* Minimum R version lowered from 4.6.0 back to 4.5.0. No code needs R 4.6;
  this lets the package be installed from GitHub on R 4.5.

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

* **Reference-intersection highlighting in `plot_signature_direction_tile_barcode()`**:
  When a `reference_df` is supplied, the function now computes the intersection
  of the signature genes with genes present in `reference_df` and visually
  distinguishes the two groups:

  - Genes **not** found in `reference_df` are shown in a muted silver-gray
    colour (`#CCCCCC`) with dimmed column labels, making the overlap immediately
    apparent without removing any genes from the display.

  - When `cluster_cols = FALSE` (default), all genes are kept in their original
    signature order (down-regulated first, then up-regulated); a `cell_fun`
    overlay is used to gray out non-intersecting columns without disturbing the
    column arrangement.

  - When `cluster_cols = TRUE`, the intersecting genes form a clustered panel
    and the non-intersecting genes are appended as a separate un-clustered panel
    in their original signature order (muted colour, labeled *not in ref*).
    When `split_direction = TRUE` is also active, this logic is applied within
    each direction sub-panel, yielding up to four panels in total
    (up-in-ref | up-not-in-ref | down-in-ref | down-not-in-ref).

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
