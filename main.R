# =============================================================================
# scRNA-seq Pipeline for GSE178481 (ccRCC) — Self-contained main script
#
# This single file contains:
#   1. All reusable wrapper functions (defined below)
#   2. The actual analysis calls on GSE178481 RCC-PR6-PTumor data
#
# Usage:
#   1. Place the count matrix at: data/RCC-PR6-PTumor.count.csv
#   2. Run:  Rscript main.R
# =============================================================================


# =============================================================================
# Part 1: Wrapper functions
# =============================================================================

# ---- Batch load / install packages ----
Plus.library <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
}

# ---- Step 1: Create a Seurat object ----
FastCreateSeurat <- function(count = NULL,
                             dir.name = NULL,
                             project = "project") {
  if (!is.null(count)) {
    sce <- CreateSeuratObject(counts = count, project = project)
  }
  if (!is.null(dir.name)) {
    sce <- CreateSeuratObject(counts = Read10X(dir.name), project = project)
  }
  return(sce)
}

# ---- Step 2: Cell quality control ----
FastSeuratCellQuality <- function(obj = NULL,
                                  species = c("human", "mouse")[1],
                                  min.features = 0,
                                  max.features = 30000000,
                                  percent.mt.num = 100) {
  sce <- obj

  if (species == "human") {
    sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
    rb.genes <- rownames(sce)[grep("^RP[SL]", rownames(sce))]
    if (length(rb.genes) > 0) {
      C <- GetAssayData(object = sce, layer = "counts")
      percent.ribo <- Matrix::colSums(C[rb.genes, , drop = FALSE]) /
                      Matrix::colSums(C) * 100
      sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
    } else {
      sce[["percent.ribo"]] <- 0
    }
  }
  if (species == "mouse") {
    sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")
    rb.genes <- rownames(sce)[grep("^Rp[sl]", rownames(sce))]
    if (length(rb.genes) > 0) {
      C <- GetAssayData(object = sce, layer = "counts")
      percent.ribo <- Matrix::colSums(C[rb.genes, , drop = FALSE]) /
                      Matrix::colSums(C) * 100
      sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
    } else {
      sce[["percent.ribo"]] <- 0
    }
  }

  initial_sce <- sce
  sce <- subset(sce,
    subset = nFeature_RNA > min.features &
             nFeature_RNA < max.features &
             percent.mt < percent.mt.num)

  return(list(initial_sce = initial_sce, sce = sce))
}

# ---- Step 3: Doublet detection (DoubletFinder) ----
FastDoubletFinder <- function(obj = NULL,
                              pcSelect = 30,
                              doublet.rate = 0.076,
                              annotation = "seurat_clusters",
                              pN_value = 0.25,
                              GT = FALSE,
                              sct = TRUE) {
  sce <- obj
  nExp <- round(doublet.rate * ncol(sce))

  sweep.res.list <- paramSweep(sce, PCs = 1:pcSelect, sct = sct)
  sweep.stats    <- summarizeSweep(sweep.res.list, GT = GT)
  bcmvn          <- find.pK(sweep.stats)
  pK             <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  annotations      <- sce@meta.data[[annotation]]
  homotypic.prop   <- modelHomotypic(annotations)
  nExp_poi         <- round(nExp * (1 - homotypic.prop))

  sce <- doubletFinder(sce, PCs = 1:pcSelect, pN = pN_value, pK = pK,
                       nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)

  df_cols <- grep("^DF.classifications", colnames(sce@meta.data), value = TRUE)
  if (length(df_cols) > 0) sce@meta.data$Doublet <- sce@meta.data[[df_cols[1]]]

  return(sce)
}

# ---- Step 4: Main RNA-seq analysis pipeline ----
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
                          algorithm = 3,
                          rmOtherGene = TRUE,
                          features = NULL,
                          outdir = "Results",
                          names = "love") {
  sce <- obj

  # QC metrics
  if (!"percent.mt" %in% colnames(sce@meta.data)) {
    if (species == "human") {
      sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
      rb.genes <- rownames(sce)[grep("^RP[SL]", rownames(sce))]
      if (length(rb.genes) > 0) {
        C <- GetAssayData(object = sce, layer = "counts")
        percent.ribo <- Matrix::colSums(C[rb.genes, , drop = FALSE]) /
                        Matrix::colSums(C) * 100
        sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
      }
    }
    if (species == "mouse") {
      sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")
      rb.genes <- rownames(sce)[grep("^Rp[sl]", rownames(sce))]
      if (length(rb.genes) > 0) {
        C <- GetAssayData(object = sce, layer = "counts")
        percent.ribo <- Matrix::colSums(C[rb.genes, , drop = FALSE]) /
                        Matrix::colSums(C) * 100
        sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
      }
    }
  }

  if (sctransform == TRUE) {
    sce <- SCTransform(sce, vars.to.regress = vars.to.regress, verbose = FALSE)
    if (is.null(features)) features <- VariableFeatures(object = sce)
  } else {
    sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)

    if (is.null(features)) {
      sce1 <- sce
      if (rmOtherGene == TRUE) {
        if (species == "human") {
          sce1 <- sce1[!grepl("^MT-", rownames(sce1)), ]
          sce1 <- sce1[!grepl("^RP[SL]", rownames(sce1)), ]
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

    if (all.scale) {
      all.genes <- rownames(sce)
      if (is.null(vars.to.regress)) sce <- ScaleData(sce, features = all.genes)
      else                          sce <- ScaleData(sce, features = all.genes, vars.to.regress = vars.to.regress)
    } else {
      if (is.null(vars.to.regress)) sce <- ScaleData(sce, features = features)
      else                          sce <- ScaleData(sce, features = features, vars.to.regress = vars.to.regress)
    }
  }

  sce <- RunPCA(sce, npcs = npcs, features = features)
  e <- ElbowPlot(sce, ndims = ncol(Embeddings(sce, "pca")))

  if (cellCycle) {
    s.genes   <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    sce <- CellCycleScoring(sce, s.features = s.genes, g2m.features = g2m.genes)
  }

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

  if (doublet) {
    if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
    sce_double <- FastDoubletFinder(sce, pcSelect = pcSelect, doublet.rate = 0.076,
                                    annotation = "seurat_clusters", pN_value = 0.25,
                                    GT = FALSE, sct = sctransform)
    sce <- subset(sce_double, Doublet == "Singlet")
  } else {
    sce_double <- NULL
  }

  if (isMarkers) {
    sce.markers <- FindAllMarkers(object = sce, test.use = test.use)
  } else {
    sce.markers <- NULL
  }

  if (plot) {
    pcaplot  <- DimPlot(sce, reduction = "pca",  label.size = 4, repel = TRUE, label = TRUE)
    tsneplot <- DimPlot(sce, reduction = "tsne", label.size = 4, repel = TRUE, label = TRUE)
    umapplot <- DimPlot(sce, reduction = "umap",  label.size = 4, repel = TRUE, label = TRUE)
    print(pcaplot);  ggsave(pcaplot,  filename = paste0(names, "_pcaplot.pdf"),  height = 10, width = 10)
    print(tsneplot); ggsave(tsneplot, filename = paste0(names, "_tsneplot.pdf"), height = 10, width = 10)
    print(umapplot); ggsave(umapplot, filename = paste0(names, "_umapplot.pdf"), height = 10, width = 10)
  }

  return(list(double_sce = sce_double, sce = sce, gene.markers = sce.markers,
              features = features, e = e))
}

# ---- Step 5: Combined QC violin plot ----
FastPlotVlnPlot <- function(object = NULL) {
  p4 <- VlnPlot(object, features = 'nCount_RNA', ncol = 1, pt.size = 0, log = TRUE, same.y.lims = TRUE) +
    geom_boxplot(outlier.size = 0, width = 0.3, show.legend = FALSE, notchwidth = 0.1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), legend.position = "none")

  p5 <- VlnPlot(object, features = 'nFeature_RNA', ncol = 1, pt.size = 0, log = TRUE) +
    geom_boxplot(outlier.size = 0, width = 0.3, show.legend = FALSE, notchwidth = 0.1) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), legend.position = "none")

  p6 <- VlnPlot(object, features = 'percent.mt', ncol = 1, pt.size = 0, y.max = 100) +
    geom_boxplot(outlier.size = 0, width = 0.3, show.legend = FALSE, notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1), legend.position = "none")

  p7 <- VlnPlot(object, features = 'percent.ribo', ncol = 1, pt.size = 0, y.max = 100) +
    geom_boxplot(outlier.size = 0, width = 0.3, show.legend = FALSE, notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1), legend.position = "none")

  plot <- cowplot::plot_grid(p4, p5, p6, p7, nrow = 4)
  return(plot)
}


# =============================================================================
# Part 2: Actual run on GSE178481 RCC-PR6-PTumor
# =============================================================================

Plus.library(c("Seurat", "dplyr", "cowplot", "patchwork", "ggplot2", "Matrix"))

outdir <- "results/gse178481"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1: Read count matrix and build Seurat object ----
cat("Step 1: Reading count matrix...\n")
counts <- read.csv("data/RCC-PR6-PTumor.count.csv", row.names = 1, check.names = FALSE)
counts <- as(as.matrix(counts), "dgCMatrix")
sce <- FastCreateSeurat(count = counts, project = "GSE178481_PR6_PTumor")
cat("  Raw:", ncol(sce), "cells x", nrow(sce), "genes\n")

# ---- Step 2: Cell quality control ----
cat("Step 2: QC filtering...\n")
qc_result <- FastSeuratCellQuality(sce, species = "human",
                                    min.features = 200, max.features = 7500, percent.mt.num = 15)
sce <- qc_result$sce
cat("  After QC:", ncol(sce), "cells\n")

# QC violin before filtering (on the unfiltered object)
p_qc_before <- VlnPlot(qc_result$initial_sce,
                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                       pt.size = 0.1, ncol = 2)
ggsave(file.path(outdir, "01_qc_violin_before_filtering.pdf"), p_qc_before, width = 12, height = 10)

# ---- Step 3: Normalize, HVG, scale, PCA, cluster, UMAP ----
cat("Step 3: Main RNA pipeline...\n")
result <- FastSeuratRNA(
  obj = sce,
  species = "human",
  pcSelect = 30,
  nfeatures = 2000,
  resolution = 0.8,
  harmony = FALSE,       # single sample -> no batch correction needed
  doublet = FALSE,       # set TRUE to run DoubletFinder
  cellCycle = TRUE,
  isMarkers = TRUE,
  algorithm = 3,         # SLM
  plot = TRUE,
  names = file.path(outdir, "GSE178481")
)
sce         <- result$sce
sce.markers <- result$gene.markers
cat("  Clusters:", length(unique(Idents(sce))), "\n")

# ---- Step 4: QC violin plot across clusters ----
cat("Step 4: QC violin by cluster...\n")
p_qc <- FastPlotVlnPlot(sce)
ggsave(file.path(outdir, "07_qc_violin_by_cluster.pdf"), p_qc, width = 8, height = 12)

# ---- Step 5: Marker heatmap ----
cat("Step 5: Marker heatmap...\n")
top10 <- sce.markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)
dh <- DoHeatmap(sce, features = top10$gene, size = 3) +
  scale_fill_gradientn(colours = c("#045a8d", "white", "#a50f15")) + NoLegend()
ggsave(file.path(outdir, "09_marker_heatmap.pdf"), dh, width = 12, height = 16)

# ---- Step 6: Cell type feature plots ----
cat("Step 6: Feature plots...\n")
markers_list <- list(
  Epithelial  = c("EPCAM", "KRT19", "KRT7", "KRT18"),
  T_cells     = c("CD3D", "CD3E", "CD8A", "CD4"),
  NK_cells    = c("GNLY", "NKG7", "FGFBP2"),
  B_cells     = c("CD19", "CD79A", "MS4A1"),
  Macrophage  = c("CD68", "CD14", "APOE"),
  Endothelial = c("ENG", "VWF"),
  Fibroblast  = c("ACTA2", "COL1A2"),
  Plasma      = c("SDC1", "CD38", "MZB1")
)
for (name in names(markers_list)) {
  genes   <- markers_list[[name]]
  present <- genes[genes %in% rownames(sce)]
  if (length(present) > 0) {
    p <- FeaturePlot(sce, features = present, ncol = 2, order = TRUE)
    ggsave(file.path(outdir, paste0("10_", name, "_feature_plot.pdf")), p, width = 12, height = 8)
  }
}

# ---- Step 7: Save ----
cat("Step 7: Saving...\n")
saveRDS(sce, file.path(outdir, "GSE178481_seurat_final.rds"))
write.csv(sce.markers, file.path(outdir, "all_markers.csv"), row.names = FALSE)

cat("\nDone. Results in:", outdir, "\n")
