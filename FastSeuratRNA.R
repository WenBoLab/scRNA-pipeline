#' Fast single-cell RNA-seq analysis pipeline
#'
#' This function performs the standard Seurat workflow: normalization, variable
#' feature selection, scaling, PCA, batch correction (Harmony), clustering,
#' dimensional reduction (UMAP/tSNE), cell cycle scoring, and marker gene detection.
#'
#' @param obj A Seurat object.
#' @param species Either "human" or "mouse". Default "human".
#' @param plot Whether to generate PCA/tSNE/UMAP plots. Default FALSE.
#' @param pcSelect Number of principal components to use downstream. Default 30.
#' @param nfeatures Number of highly variable features to select. Default 2000.
#' @param sctransform Whether to use SCTransform instead of standard normalization. Default FALSE.
#' @param vars.to.regress Variables to regress out during scaling. Default NULL.
#' @param all.scale Whether to scale all genes (not just variable features). Default FALSE.
#' @param npcs Number of PCs to compute. Default 50.
#' @param resolution Clustering resolution. Default 0.5.
#' @param harmony Whether to use Harmony batch correction. Default FALSE.
#' @param harmony_by Meta.data column to group by for Harmony correction. Default "orig.ident".
#' @param doublet Whether to run doublet detection. Default FALSE.
#' @param perplexity Perplexity parameter for tSNE/UMAP. Default 30.
#' @param cellCycle Whether to calculate cell cycle scores. Default TRUE.
#' @param test.use Marker gene detection test. Default "wilcox".
#' @param isMarkers Whether to find cluster markers. Default TRUE.
#' @param algorithm Clustering algorithm. 1=Louvin, 2=Louvin multilevel, 3=SLM, 4=Leiden. Default 4.
#' @param rmOtherGene Whether to remove MT/ribo/IG genes before HVG selection. Default TRUE.
#' @param features User-specified features to use. Default NULL (auto-select).
#' @param filepath Path to save the final Seurat object. Default NULL.
#' @param outdir Output directory for doublet detection results. Default "Results".
#' @param names Prefix for output plot files. Default "love".
#'
#' @return A list with the following elements:
#'   \item{double_sce}{Doublet detection result (NULL if doublet=FALSE).}
#'   \item{sce}{The processed Seurat object.}
#'   \item{gene.markers}{Differentially expressed genes per cluster.}
#'   \item{features}{The variable features used.}
#'   \item{e}{ElbowPlot object.}
#' @export
#'
#' @examples
#' \dontrun{
#' result <- FastSeuratRNA(sce, species = "human", harmony = TRUE,
#'                         resolution = 0.8, pcSelect = 30)
#' sce <- result$sce
#' markers <- result$gene.markers
#' }
FastSeuratRNA <- function(obj = NULL,
                          species = c("human", "mouse")[1],
                          plot = FALSE,
                          pcSelect = 30,
                          nfeatures = 2000,
                          sctransform = FALSE,
                          vars.to.regress = NULL,
                          all.scale = FALSE,
                          npcs = 50,
                          resolution = 0.5,
                          harmony = FALSE,
                          harmony_by = "orig.ident",
                          doublet = FALSE,
                          perplexity = 30,
                          cellCycle = TRUE,
                          test.use = c("wilcox", "LR", "MAST")[1],
                          isMarkers = TRUE,
                          algorithm = c(1, 2, 3, 4)[4],
                          rmOtherGene = TRUE,
                          features = NULL,
                          filepath = NULL,
                          outdir = "Results",
                          names = "love") {
  # Load required packages
  nds <- c("Seurat", "dplyr", "cowplot", "patchwork", "harmony", "SoupX", "clustree")
  Plus.library(nds)

  if (!is.null(obj)) {
    sce <- obj
  }

  # Add QC metrics if not already present
  if (!"percent.mt" %in% colnames(sce@meta.data)) {
    if (species == "human") {
      sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
      rb.genes <- rownames(sce)[grep("^RP[SL]", rownames(sce))]
      C <- GetAssayData(object = sce, layer = "counts")
      percent.ribo <- Matrix::colSums(C[rb.genes, ]) / Matrix::colSums(C) * 100
      sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
    }
    if (species == "mouse") {
      sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")
      rb.genes <- rownames(sce)[grep("^Rp[sl]", rownames(sce))]
      C <- GetAssayData(object = sce, layer = "counts")
      percent.ribo <- Matrix::colSums(C[rb.genes, ]) / Matrix::colSums(C) * 100
      sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
    }
  }

  # Normalization and variable feature selection
  if (sctransform == TRUE) {
    sce <- SCTransform(sce, vars.to.regress = vars.to.regress, verbose = FALSE)
    if (is.null(features)) {
      features <- VariableFeatures(object = sce)
    } else {
      features <- features[features %in% rownames(sce)]
    }
  } else {
    # Normalize data
    sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)

    # Find variable features
    if (is.null(features)) {
      sce1 <- sce

      if (rmOtherGene == TRUE) {
        if (species == "human") {
          # Remove mitochondrial genes
          sce1 <- sce1[!grepl("^MT-", rownames(sce1)), ]
          # Remove ribosomal genes
          sce1 <- sce1[!grepl("^RP[SL]", rownames(sce1)), ]
          # Remove immunoglobulin/TCR genes
          sce1 <- sce1[!grepl("^HLA*|^IGHV*|^IGHJ*|^IGHD*|^IGKV*|^IGLV*|^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*", rownames(sce1)), ]
        }
        if (species == "mouse") {
          sce1 <- sce1[!grepl("^mt-", rownames(sce1)), ]
          sce1 <- sce1[!grepl("^Rp[sl]", rownames(sce1)), ]
        }
      }

      sce1 <- FindVariableFeatures(sce1, selection.method = "vst", nfeatures = nfeatures)
      features <- VariableFeatures(object = sce1)
    } else {
      features <- features[features %in% rownames(sce)]
    }

    # Scale data
    if (all.scale) {
      all.genes <- rownames(sce)
      if (is.null(vars.to.regress)) {
        sce <- ScaleData(sce, features = all.genes)
      } else {
        sce <- ScaleData(sce, features = all.genes, vars.to.regress = vars.to.regress)
      }
    } else {
      if (is.null(vars.to.regress)) {
        sce <- ScaleData(sce, features = features)
      } else {
        sce <- ScaleData(sce, features = features, vars.to.regress = vars.to.regress)
      }
    }
  }

  # PCA
  sce <- RunPCA(sce, npcs = npcs, features = features)
  e <- ElbowPlot(sce, ndims = ncol(Embeddings(sce, "pca")))

  # Cell cycle scoring
  if (cellCycle) {
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    sce <- CellCycleScoring(sce, s.features = s.genes, g2m.features = g2m.genes)
  }

  # Batch correction, clustering, and dimensional reduction
  if (harmony) {
    sce <- RunHarmony(sce, group.by.vars = harmony_by, dims.use = 1:pcSelect, max_iter = 50)
    sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:pcSelect)
    sce <- FindClusters(sce, resolution = resolution, algorithm = algorithm)
    sce <- RunTSNE(object = sce, reduction = "harmony", dims = 1:pcSelect, do.fast = TRUE, perplexity = perplexity)
    sce <- RunUMAP(sce, reduction = "harmony", dims = 1:pcSelect, perplexity = perplexity)
  } else {
    sce <- FindNeighbors(sce, dims = 1:pcSelect)
    sce <- FindClusters(sce, resolution = resolution, algorithm = algorithm)
    sce <- RunTSNE(object = sce, dims = 1:pcSelect, do.fast = TRUE, perplexity = perplexity)
    sce <- RunUMAP(sce, dims = 1:pcSelect, perplexity = perplexity)
  }

  # Doublet detection
  if (doublet) {
    if (!file.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }
    sce_double <- FastDoubletFinder(
      sce,
      pcSelect = pcSelect,
      doublet.rate = 0.076,
      annotation = "seurat_clusters",
      pN_value = 0.25,
      GT = FALSE,
      sct = sctransform
    )
    sce <- subset(sce_double, Doublet == "Singlet")
  } else {
    sce_double <- NULL
  }

  # Find marker genes
  if (isMarkers) {
    sce.markers <- FindAllMarkers(object = sce, test.use = test.use)
  } else {
    sce.markers <- NULL
  }

  # Generate plots
  if (plot) {
    pcaplot <- DimPlot(sce, reduction = "pca", label.size = 4, repel = TRUE, label = TRUE)
    print(pcaplot)
    ggsave(pcaplot, filename = paste0(names, "_pcaplot.pdf"), height = 10, width = 10)

    tsneplot <- DimPlot(object = sce, reduction = "tsne", label.size = 4, repel = TRUE, label = TRUE)
    print(tsneplot)
    ggsave(tsneplot, filename = paste0(names, "_tsneplot.pdf"), height = 10, width = 10)

    umapplot <- DimPlot(object = sce, reduction = "umap", label.size = 4, repel = TRUE, label = TRUE)
    print(umapplot)
    ggsave(umapplot, filename = paste0(names, "_umapplot.pdf"), height = 10, width = 10)
  }

  # Return results
  result <- list(
    double_sce = sce_double,
    sce = sce,
    gene.markers = sce.markers,
    features = features,
    e = e
  )

  if (!is.null(filepath)) {
    save(sce, file = filepath)
  }

  return(result)
}
