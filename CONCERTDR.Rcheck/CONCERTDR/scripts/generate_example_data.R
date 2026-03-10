## Script to regenerate inst/extdata/example_reference_df.csv
## Run once from the package root: source("inst/scripts/generate_example_data.R")

set.seed(42)

genes <- c("TP53","MYC","STAT3","EGFR","KRAS","AKT1","BCL2","CDK2","CCND1","MAPK1",
           "PIK3CA","PTEN","JAK2","HDAC1","SRC","ABL1","NOTCH1","VEGFA","FOS","JUN")

sig_ids <- paste0("DEMO", sprintf("%03d", 1:10))

# Simulate z-score-like expression (mean 0, sd ~2)
# DEMO001/002 (imatinib in K562): anti-correlated with the query signature
# to simulate a "drug reversal" hit
mat <- matrix(rnorm(length(genes) * length(sig_ids), mean = 0, sd = 2),
              nrow = length(genes), ncol = length(sig_ids),
              dimnames = list(genes, sig_ids))

# Signature log2FC direction: negative genes = TP53,EGFR,AKT1,BCL2,MAPK1,PTEN,HDAC1,ABL1,FOS
down_genes <- c("TP53","EGFR","AKT1","BCL2","MAPK1","PTEN","HDAC1","ABL1","FOS")
up_genes   <- setdiff(genes, down_genes)

# imatinib: invert the pattern (high score = reversal drug)
mat[down_genes, "DEMO001"] <- mat[down_genes, "DEMO001"] + 3
mat[up_genes,   "DEMO001"] <- mat[up_genes,   "DEMO001"] - 3
mat[down_genes, "DEMO002"] <- mat[down_genes, "DEMO002"] + 2
mat[up_genes,   "DEMO002"] <- mat[up_genes,   "DEMO002"] - 2

ref_df <- data.frame(gene_symbol = rownames(mat), mat,
                     stringsAsFactors = FALSE, check.names = FALSE)

out <- file.path("inst", "extdata", "example_reference_df.csv")
write.csv(ref_df, out, row.names = FALSE)
message("Written: ", out)
