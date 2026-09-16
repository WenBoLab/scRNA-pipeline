#' Cell quality control for single-cell RNA-seq
#'
#' This function calculates mitochondrial and ribosomal gene percentages,
#' then filters cells based on feature count and mitochondrial threshold.
#'
#' @param obj A Seurat object.
#' @param species Either "human" or "mouse". Default "human".
#' @param min.features Minimum number of detected genes per cell. Default 0.
#' @param max.features Maximum number of detected genes per cell. Default 30000000.
#' @param percent.mt.num Maximum mitochondrial percentage allowed. Default 100.
#'
#' @return A list with two elements:
#'   \item{initial_sce}{The unfiltered Seurat object with QC metrics added.}
#'   \item{sce}{The filtered Seurat object.}
#' @export
#'
#' @examples
#' \dontrun{
#' qc_result <- FastSeuratCellQuality(sce, species = "human",
#'                                    min.features = 200, max.features = 7500,
#'                                    percent.mt.num = 15)
#' sce_filtered <- qc_result$sce
#' }
FastSeuratCellQuality <- function(obj = NULL,
                                  species = c("human", "mouse")[1],
                                  min.features = 0,
                                  max.features = 30000000,
                                  percent.mt.num = 100) {
  # Load required packages
  nds <- c("Seurat", "dplyr", "cowplot", "patchwork", "harmony", "SoupX", "presto")
  Plus.library(nds)

  if (!is.null(obj)) {
    sce <- obj
  }

  if (species == "human") {
    # Add mitochondrial gene percentage
    sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")

    # Add ribosomal gene percentage
    rb.genes <- rownames(sce)[grep("^RP[SL]", rownames(sce))]
    C <- GetAssayData(object = sce, layer = "counts")
    percent.ribo <- Matrix::colSums(C[rb.genes, ]) / Matrix::colSums(C) * 100
    sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
  }

  if (species == "mouse") {
    # Add mitochondrial gene percentage
    sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")

    # Add ribosomal gene percentage
    rb.genes <- rownames(sce)[grep("^Rp[sl]", rownames(sce))]
    C <- GetAssayData(object = sce, layer = "counts")
    percent.ribo <- Matrix::colSums(C[rb.genes, ]) / Matrix::colSums(C) * 100
    sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
  }

  # Store initial (unfiltered) object
  initial_sce <- sce

  # Filter cells
  sce <- subset(
    sce,
    subset = nFeature_RNA > min.features &
      nFeature_RNA < max.features &
      percent.mt < percent.mt.num
  )

  result <- list(
    initial_sce = initial_sce,
    sce = sce
  )

  return(result)
}
