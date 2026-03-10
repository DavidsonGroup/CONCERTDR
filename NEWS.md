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
