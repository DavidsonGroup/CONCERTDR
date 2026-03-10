![image](https://github.com/user-attachments/assets/9b07b382-d65e-4ec1-a761-b71dd6318989)

# CONCERTDR: CONtext-aware Cellular and Tissue-specific Expression for Drug Repurposing

CONCERTDR is an R package for drug repurposing against the CMap L1000 database.
It covers CMap condition filtering, signature matching (KS, XCos, XSum, GSEA, Zhang),
result annotation, z-score extraction, and heatmap visualization.

The recommended pattern is to set all CMap file paths as global options once at
the top of your script. `extract_signature_zscores()` and
`plot_signature_direction_tile_barcode()` read these automatically; for other
functions that need explicit paths, pass `getOption(...)`.

---

## Installation

```r
install.packages("remotes")
remotes::install_github("DavidsonGroup/CONCERT_DR")
```

Visualization dependencies:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages("circlize")
```

---

## Setup

Download the CMap data files from https://clue.io/releases/data-dashboard:

- `level5_beta_all_n1201944x12328.gctx`
- `siginfo_beta.txt`
- `geneinfo_beta.txt`
- `compoundinfo_beta.txt`

Also, the repurposing_drugs.txt is available at
[<https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt>]
on Broad Institute website for drug compound reference.

Set these once at the top of any script before calling other functions. Use plain
paths or `file.path()` — **do not use `system.file()` here**; that function only
looks inside installed R packages and will return `""` for external files.

```r
library(CONCERTDR)

data_dir <- "/path/to/cmap_data"

options(
  CONCERTDR.data_dir          = data_dir,
  CONCERTDR.gctx_file         = file.path(data_dir, "level5_beta_all_n1201944x12328.gctx"),
  CONCERTDR.siginfo_file      = file.path(data_dir, "siginfo_beta.txt"),
  CONCERTDR.geneinfo_file     = file.path(data_dir, "geneinfo_beta.txt"),
  CONCERTDR.compoundinfo_file = file.path(data_dir, "compoundinfo_beta.txt")
)
```

You can sanity-check the paths before running:

```r
stopifnot(file.exists(getOption("CONCERTDR.gctx_file")))
stopifnot(file.exists(getOption("CONCERTDR.geneinfo_file")))
```

---

## Quick-start (bundled test data)

The package ships a small example dataset in `inst/extdata/` (20 genes × 10 signatures),
so you can run the full pipeline without the multi-GB GCTX.
`system.file()` is used here because the files live *inside the package* — do not
use it for your own data files on disk.

```r
library(CONCERTDR)

sig_file <- system.file("extdata", "example_signature.txt",    package = "CONCERTDR")
ref_csv  <- system.file("extdata", "example_reference_df.csv", package = "CONCERTDR")

reference_df <- read.csv(ref_csv, row.names = 1, check.names = FALSE)
reference_df$gene_symbol <- rownames(reference_df)

res <- process_signature_with_df(
  reference_df   = reference_df,
  signature_file = sig_file,
  methods        = c("ks","xsum"),
  topN           = 4,
  permutations   = 10,
  save_files     = FALSE
)

summary(res)
head(res$results$ks)
```

The bundled example data also support z-score extraction and plotting directly
from `example_reference_df.csv`:

```r
z <- extract_signature_zscores(
  results_df     = res$results$ks,
  signature_file = sig_file,
  reference_df   = reference_df,
  max_genes      = 20,
  max_perts      = 10
)

plot_signature_direction_tile_barcode(
  precomputed  = z,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)
```

If you have the original CMap files, you can still pass `gctx_file`,
`geneinfo_file`, and `siginfo_file` to extract the heatmap matrix from GCTX
instead of using the precomputed reference matrix.

By default, the core analysis and plotting functions do not write files. Files
are only created when you explicitly request them, for example with
`save_files = TRUE`, `write_outputs = TRUE`, `output_zscores = ...`, or
`save_png = TRUE`.

---

## Workflow

### Step 1 — Filter siginfo

Narrow down to the cell lines and perturbation types relevant to your analysis.
The result is an in-memory data frame that gets passed directly to later steps.

```r
filtered_siginfo <- subset_siginfo_beta(
  getOption("CONCERTDR.siginfo_file"),
  interactive = FALSE,
  filters = list(
    pert_type  = c("trt_xpr", "trt_cp", "trt_oe"),
    pert_itime = c("6 h", "24 h"),
    cell_iname = c("HL60", "THP1", "K562", "HAP1", "JURKAT", "U937")
  )
)
```

Use `interactive = TRUE` to select parameters interactively from the console.
Filterable columns: `pert_type`, `pert_itime`, `pert_idose`, `cell_iname`.

The `pert_type` values are defined by LINCS2020 as follows:

![image](https://github.com/user-attachments/assets/02ef148d-736b-4c02-92b5-5fde0935db17)


### Step 2 — Signature file

A tab-delimited file with `Gene` (HGNC symbol) and `log2FC`. Positive = up-regulated,
negative = down-regulated. Rows should be sorted by `|log2FC|` descending; `max_genes`
in later steps controls how many are used.

```text
Gene    log2FC
STAT3    2.5
TP53    -1.8
MYC      3.2
BRCA1   -2.1
```

### Step 3 — Build the reference matrix

```r
reference_df <- extract_cmap_data_from_siginfo(
  geneinfo_file  = getOption("CONCERTDR.geneinfo_file"),
  siginfo_file   = filtered_siginfo,   # data frame from Step 1, not re-read from disk
  gctx_file      = getOption("CONCERTDR.gctx_file"),
  filter_quality = FALSE,              # TRUE restricts to is_hiq == 1
  landmark       = TRUE               # TRUE restricts to 978 landmark genes
)
```

### Step 4 — Score the signature

```r
results <- process_signature_with_df(
  reference_df   = reference_df,
  signature_file = "signature.txt",
  output_dir     = "results",
  methods        = c("xsum", "xcos", "zhang", "gsea0", "gsea1", "gsea2", "ks"),
  topN           = 400,    # genes used by XCos/XSum
  permutations   = 1,      # increase for p-values
  save_files     = FALSE   # set TRUE only if you want CSV outputs written
)

print(results)            # brief overview
summary(results)          # top hits across methods
head(results$results$ks)
```

Available methods:

| ID | Description |
|---|---|
| `ks` | Kolmogorov–Smirnov |
| `xcos` | Extreme cosine similarity |
| `xsum` | Extreme sum |
| `gsea0` / `gsea1` / `gsea2` | GSEA (weight 0 / 1 / 2) |
| `zhang` | Zhang et al. |

### Step 5 — Annotate

Joins `siginfo` and `compoundinfo` to add drug names, MoA, dose, cell line, etc.
`filtered_siginfo` is already in memory from Step 1 — no file re-read.

```r
views <- annotate_drug_results(
  results_df     = results$results$ks,
  sig_info_file  = filtered_siginfo,
  comp_info_file = getOption("CONCERTDR.compoundinfo_file"),
  output_dir     = "results",
  write_outputs  = TRUE
)

head(views$tech_view_all)        # per-signature scores — pass to Step 6
head(views$wetlab_drug_view)     # drug-level summary
head(views$drug_context_summary) # drug × cell line breakdown
```

### Step 6 — Extract z-scores

Reads the GCTX once. Save the result and reuse it rather than repeating this call.

```r
z <- extract_signature_zscores(
  results_df     = views$tech_view_all,
  signature_file = "signature.txt",
  max_genes      = 100,
  max_perts      = 60,
  output_zscores = NULL                     # set a path if you want to save TSV
)
```

The returned list contains:

- `z_plot` — matrix (perturbations × genes), rows labelled `drug | dose | time | cell`
- `ordered_genes` — gene order (down-regulated → up-regulated)
- `logfc_map` — named vector of log2FC values
- `sig_ids` / `sig_labels` — CMap identifiers and readable labels

### Step 7 — Plot

Pass `precomputed = z` to skip any file I/O on repeat calls.

```r
plot_signature_direction_tile_barcode(
  precomputed         = z,
  cluster_rows        = TRUE,
  cluster_cols        = FALSE,   # keeps down → up gene order from signature
  show_row_dendrogram = TRUE,
  save_png            = TRUE,    # default is FALSE
  output_png          = "results/barcode_heatmap.png"
)
```

Colour scheme: coolwarm (`#3B4CC0` → white → `#B40426`) for z-scores; BrBG
(`#01665E` → white → `#8C510A`) for the log2FC annotation bar. Figure dimensions
are computed automatically from data size (~0.22 in/gene, ~0.28 in/perturbation),
or set explicitly with `width`, `height`, `dpi`.

---


## Function reference

### Database subsetting

| Function | Description |
|---|---|
| `subset_siginfo_beta()` | Filter `siginfo_beta.txt` to a cell-line/pert-type subset |

### Data extraction

| Function | Description |
|---|---|
| `extract_cmap_data_from_siginfo()` | Build reference matrix from filtered siginfo |

### Signature matching

| Function | Description |
|---|---|
| `process_signature_with_df()` | Score a signature against an in-memory reference |

### Result annotation

| Function | Description |
|---|---|
| `annotate_drug_results()` | Join CMap metadata; produce wetlab/tech view tables |
| `fuzzy_drug_match()` | Fuzzy-match drug names across databases |
| `extract_compound_id()` | Resolve compound identifiers |

### Visualization

| Function | Description |
|---|---|
| `extract_signature_zscores()` | Extract GCTX z-score matrix as an R object |
| `plot_signature_direction_tile_barcode()` | ComplexHeatmap barcode plot (z-score × perturbation) |

## Developing / packaging

All documentation is written as `#'` roxygen2 blocks in the `.R` source files.
Do not edit `NAMESPACE` or `man/*.Rd` directly.

```r
devtools::document()  # regenerate man/ and NAMESPACE
devtools::check()     # R CMD check
devtools::install()   # install locally for testing
```

---

## License

MIT — see the `LICENSE` file for details.
