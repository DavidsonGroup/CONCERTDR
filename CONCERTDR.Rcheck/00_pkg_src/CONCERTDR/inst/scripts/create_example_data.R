## inst/scripts/create_example_data.R
##
## Run ONCE from the package root on a machine that has access to the real
## CMap release files (e.g., the lab HPC).  The script produces five small
## files in inst/extdata/ that are committed to the repository so that users
## can run the Quick-start vignette without downloading the full GCTX.
##
## Usage:
##   cd /path/to/CONCERT_DR_updated
##   Rscript inst/scripts/create_example_data.R
##
## Optionally override default paths with environment variables:
##   CMAP_GCTX     path to level5_beta_all_n1201944x12328.gctx
##   CMAP_GENEINFO path to geneinfo_beta.txt
##   CMAP_SIGINFO  path to siginfo_beta.txt

suppressPackageStartupMessages({
  library(cmapR)
  library(data.table)
})

# ── 0. Configurable paths ────────────────────────────────────────────────────
# Set this to the directory containing your CMap release files,
# or override with environment variables CMAP_GCTX / CMAP_GENEINFO / CMAP_SIGINFO.
default_dir <- "/path/to/cmap_data"

gctx_file     <- Sys.getenv("CMAP_GCTX",     file.path(default_dir, "level5_beta_all_n1201944x12328.gctx"))
geneinfo_file <- Sys.getenv("CMAP_GENEINFO",  file.path(default_dir, "geneinfo_beta.txt"))
siginfo_file  <- Sys.getenv("CMAP_SIGINFO",   file.path(default_dir, "siginfo_beta.txt"))
out_dir       <- "inst/extdata"

for (f in c(gctx_file, geneinfo_file, siginfo_file)) {
  if (!file.exists(f)) stop("File not found: ", f)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Gene selection (landmark genes only) ───────────────────────────────────
target_genes <- c(
  "TP53",   "MYC",    "STAT3",  "EGFR",   "KRAS",
  "AKT1",   "BCL2",   "CDK2",   "CCND1",  "MAPK1",
  "PIK3CA", "PTEN",   "JAK2",   "HDAC1",  "SRC",
  "ABL1",   "NOTCH1", "VEGFA",  "FOS",    "JUN"
)

message("Reading geneinfo...")
geneinfo <- fread(geneinfo_file)

# Prefer landmark; fall back to best_inferred if a gene isn't landmark
gene_sub_lm <- geneinfo[gene_symbol %in% target_genes & feature_space == "landmark"]
missing_lm  <- setdiff(target_genes, gene_sub_lm$gene_symbol)
if (length(missing_lm) > 0) {
  gene_sub_inf <- geneinfo[gene_symbol %in% missing_lm & feature_space == "best inferred"]
  gene_sub     <- rbind(gene_sub_lm, gene_sub_inf)
  message(sprintf("  %d landmark  +  %d best-inferred genes selected",
                  nrow(gene_sub_lm), nrow(gene_sub_inf)))
} else {
  gene_sub <- gene_sub_lm
  message(sprintf("  All %d target genes found in landmark space", nrow(gene_sub)))
}

# ── 2. Signature selection (2 per drug, prefer K562 then any cell) ─────────────
target_drugs <- c("imatinib", "dasatinib", "doxorubicin", "methotrexate", "vincristine")

message("Reading siginfo...")
siginfo <- fread(siginfo_file)

pick_sigs <- function(drug) {
  pool <- siginfo[cmap_name == drug & pert_type %in% c("trt_cp", "trt_xpr")]
  if (nrow(pool) == 0) return(NULL)
  # Prefer K562 (ABL1-positive, well-characterised for BCR-ABL inhibitors)
  k562 <- pool[cell_iname == "K562"]
  chosen <- if (nrow(k562) >= 2) k562 else pool
  head(chosen[order(sig_id)], 2)
}

sig_sub <- rbindlist(lapply(target_drugs, pick_sigs))
message(sprintf("  Selected %d signatures across %d drugs", nrow(sig_sub), length(target_drugs)))
print(sig_sub[, .(sig_id, cmap_name, cell_iname, pert_dose, pert_time)])

# ── 3. Parse GCTX subset ──────────────────────────────────────────────────────
rid <- as.character(gene_sub$gene_id)
cid <- sig_sub$sig_id

message(sprintf("Subsetting GCTX (%d genes × %d sigs)...", length(rid), length(cid)))
ds <- parse_gctx(gctx_file, rid = rid, cid = cid)
message("  Matrix dimensions: ", paste(dim(ds@mat), collapse = " x "))

# ── 4. Write example_data.gctx ───────────────────────────────────────────────
out_gctx <- file.path(out_dir, "example_data")
write_gctx(ds, out_gctx)
message("Written: ", out_gctx, ".gctx  (",
        round(file.info(paste0(out_gctx, ".gctx"))$size / 1024, 1), " KB)")

# ── 5. Write example_geneinfo.txt ─────────────────────────────────────────────
keep_gene_cols <- intersect(
  c("gene_id", "gene_symbol", "gene_title", "gene_type", "src", "feature_space"),
  colnames(geneinfo)
)
fwrite(gene_sub[, ..keep_gene_cols],
       file.path(out_dir, "example_geneinfo.txt"), sep = "\t")
message("Written: example_geneinfo.txt")

# ── 6. Write example_siginfo.txt ──────────────────────────────────────────────
keep_sig_cols <- intersect(
  c("sig_id", "pert_id", "cmap_name", "pert_type", "cell_iname",
    "pert_dose", "pert_dose_unit", "pert_time", "pert_time_unit",
    "distil_cc_q75", "distil_ss"),
  colnames(siginfo)
)
fwrite(sig_sub[, ..keep_sig_cols],
       file.path(out_dir, "example_siginfo.txt"), sep = "\t")
message("Written: example_siginfo.txt")

# ── 7. Write example_reference_df.csv ─────────────────────────────────────────
mat      <- ds@mat
sym_map  <- setNames(gene_sub$gene_symbol, as.character(gene_sub$gene_id))
row_syms <- sym_map[rownames(mat)]

ref_df <- data.frame(
  gene_symbol = row_syms,
  as.data.frame(mat, stringsAsFactors = FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.csv(ref_df, file.path(out_dir, "example_reference_df.csv"), row.names = FALSE)
message("Written: example_reference_df.csv")

# ── 8. Write example_signature.txt (20 genes, simulated log2FC) ───────────────
# Keep the gene set but use reproducible synthetic fold-changes so the
# quick-start demo is fully self-contained and doesn't require biological truth.
set.seed(42)
sig_genes <- gene_sub$gene_symbol
sig_lfc   <- round(rnorm(length(sig_genes), mean = 0, sd = 1.5), 3)
# Amplify a few genes to make the up/down split clear
sig_lfc[sig_genes %in% c("MYC","STAT3","KRAS","CDK2","CCND1","PIK3CA","JAK2","NOTCH1","VEGFA","JUN")] <-
  abs(sig_lfc[sig_genes %in% c("MYC","STAT3","KRAS","CDK2","CCND1","PIK3CA","JAK2","NOTCH1","VEGFA","JUN")])
sig_lfc[sig_genes %in% c("TP53","EGFR","AKT1","BCL2","MAPK1","PTEN","HDAC1","SRC","ABL1","FOS")] <-
  -abs(sig_lfc[sig_genes %in% c("TP53","EGFR","AKT1","BCL2","MAPK1","PTEN","HDAC1","SRC","ABL1","FOS")])

sig_df <- data.frame(Gene = sig_genes, log2FC = sig_lfc)
write.table(sig_df,
            file.path(out_dir, "example_signature.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Written: example_signature.txt")

message("\n── Done ─────────────────────────────────────────────────────────────────────")
message("Files in ", out_dir, ":")
for (f in list.files(out_dir, full.names = TRUE)) {
  message(sprintf("  %-40s  %s KB",
                  basename(f),
                  round(file.info(f)$size / 1024, 1)))
}
message("\nCommit inst/extdata/ to Git so users don't need the full GCTX for Quick-start.")
