#' Extract Data from CMap GCTX Files for Specified Combinations
#'
#' @description Functions to extract expression data from CMap GCTX files based on
#' specified combinations of time points, dosages, and cell lines.
#' @return None; this is an internal documentation topic.
#'
#' @name cmap_extract
#' @keywords internal
NULL

#' Get Gene IDs and Gene Names for Landmark Genes
#'
#' @param geneinfo_df Data frame containing gene information
#' @return List with rid (gene IDs) and genenames (gene symbols)
#' @keywords internal
get_rid <- function(geneinfo_df,landmark=TRUE) {
  gene_overlapping <- geneinfo_df
  if (landmark)
  {
    gene_overlapping <- gene_overlapping[gene_overlapping$feature_space == "landmark", ]
  }

  return(list(
    rid = as.character(gene_overlapping$gene_id),
    genenames = as.character(gene_overlapping$gene_symbol)
  ))
}
