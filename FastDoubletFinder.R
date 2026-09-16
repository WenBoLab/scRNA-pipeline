#' Fast doublet detection using DoubletFinder
#'
#' This function wraps the DoubletFinder workflow for identifying doublets
#' in single-cell RNA-seq data. It follows the standard DoubletFinder protocol:
#' PCA -> neighbors -> UMAP -> paramSweep -> model selection -> doublet classification.
#'
#' @param obj A Seurat object that has already been processed through PCA and clustering.
#' @param pcSelect Number of principal components to use. Default 30.
#' @param doublet.rate Expected doublet rate. Default 0.076 (7.6%).
#' @param annotation Cluster annotation column in meta.data. Default "seurat_clusters".
#' @param pN_value Number of doublets to synthesize (proportion). Default 0.25.
#' @param GT Whether ground truth labels are available. Default FALSE.
#' @param sct Whether SCTransform was used. Default TRUE.
#'
#' @return A Seurat object with doublet classifications added in the "Doublet" meta.data column.
#' @export
#'
#' @examples
#' \dontrun{
#' sce <- FastDoubletFinder(sce, pcSelect = 30, doublet.rate = 0.076)
#' sce <- subset(sce, Doublet == "Singlet")
#' }
FastDoubletFinder <- function(obj = NULL,
                              pcSelect = 30,
                              doublet.rate = 0.076,
                              annotation = "seurat_clusters",
                              pN_value = 0.25,
                              GT = FALSE,
                              sct = TRUE) {
  # Load required packages
  nds <- c("Seurat", "DoubletFinder")
  Plus.library(nds)

  sce <- obj

  # Estimate number of expected doublets
  nExp <- round(doublet.rate * ncol(sce))

  # Run parameter sweep
  sweep.res.list <- paramSweep(
    sce,
    PCs = 1:pcSelect,
    sct = sct
  )

  # Summarize sweep results
  sweep.stats <- summarizeSweep(sweep.res.list, GT = GT)

  # Find optimal pK
  bcmvn <- find.pK(sweep.stats)
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  # Adjust expected doublets by homotypic proportion
  annotations <- sce@meta.data[[annotation]]
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(nExp * (1 - homotypic.prop))

  # Run DoubletFinder
  sce <- doubletFinder(
    sce,
    PCs = 1:pcSelect,
    pN = pN_value,
    pK = pK,
    nExp = nExp_poi,
    reuse.pANN = FALSE,
    sct = sct
  )

  # Rename the classification column to "Doublet"
  df_cols <- grep("^DF.classifications", colnames(sce@meta.data), value = TRUE)
  if (length(df_cols) > 0) {
    sce@meta.data$Doublet <- sce@meta.data[[df_cols[1]]]
  }

  return(sce)
}
