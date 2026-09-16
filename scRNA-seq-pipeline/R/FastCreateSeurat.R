#' Create a Seurat object from 10X data or a count matrix
#'
#' This function wraps Seurat::CreateSeuratObject with standard defaults for
#' single-cell RNA-seq analysis. It accepts either a count matrix or a path
#' to a 10X Genomics output directory.
#'
#' @param count A count matrix (genes x cells). Default NULL.
#' @param dir.name Path to a 10X output directory containing barcodes/features/matrix files.
#'   Default NULL.
#' @param project Project name stored in the Seurat object. Default "project".
#'
#' @return A Seurat object.
#' @export
#'
#' @examples
#' \dontrun{
#' sce <- FastCreateSeurat(dir.name = "data/GSE178481", project = "GSE178481")
#' }
FastCreateSeurat <- function(count = NULL,
                             dir.name = NULL,
                             project = "project") {
  # Load required packages
  nds <- c("Seurat", "dplyr", "cowplot", "patchwork", "harmony", "SoupX", "presto")
  Plus.library(nds)

  # Create Seurat object
  if (!is.null(count)) {
    sce <- CreateSeuratObject(
      counts = count,
      project = project
    )
  }

  if (!is.null(dir.name)) {
    sce <- CreateSeuratObject(
      counts = Read10X(dir.name),
      project = project
    )
  }

  return(sce)
}
