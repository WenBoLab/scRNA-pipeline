# =============================================================================
# Demo script: run the pipeline on pbmc_small to generate example figures
# This demonstrates the full workflow on a small dataset
# =============================================================================

library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(harmony)

# Load demo data
data("pbmc_small")
sce <- pbmc_small

cat("Initial object:", ncol(sce), "cells and", nrow(sce), "genes\n")

# =============================================================================
# Step 1: Add QC metrics (simulate the FastSeuratCellQuality step)
# =============================================================================
cat("Step 1: Adding QC metrics...\n")

sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
rb.genes <- rownames(sce)[grep("^RP[SL]", rownames(sce))]
if (length(rb.genes) > 0) {
  C <- GetAssayData(object = sce, layer = "counts")
  percent.ribo <- Matrix::colSums(C[rb.genes, , drop = FALSE]) / Matrix::colSums(C) * 100
  sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
} else {
  sce[["percent.ribo"]] <- 0
}

# =============================================================================
# Step 2: Normalization and variable features
# =============================================================================
cat("Step 2: Normalization and variable feature selection...\n")

sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 100)
features <- VariableFeatures(object = sce)
sce <- ScaleData(sce, features = features)

# =============================================================================
# Step 3: PCA
# =============================================================================
cat("Step 3: Running PCA...\n")

sce <- RunPCA(sce, npcs = 20, features = features)
pca_plot <- DimPlot(sce, reduction = "pca", label = TRUE) + ggtitle("PCA")
ggsave("results/demo_pca_plot.pdf", pca_plot, width = 8, height = 6)

# Elbow plot
elbow_plot <- ElbowPlot(sce, ndims = 20)
ggsave("results/demo_elbow_plot.pdf", elbow_plot, width = 8, height = 6)

# =============================================================================
# Step 4: Cell cycle scoring (skipped in demo due to small dataset)
# =============================================================================
cat("Step 4: Cell cycle scoring (skipped in demo dataset)...\n")

# s.genes <- cc.genes.updated.2019$s.genes
# g2m.genes <- cc.genes.updated.2019$g2m.genes
# sce <- CellCycleScoring(sce, s.features = s.genes, g2m.features = g2m.genes)

# =============================================================================
# Step 5: Clustering and UMAP
# =============================================================================
cat("Step 5: Clustering and UMAP...\n")

sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 0.8, algorithm = 3)
sce <- RunUMAP(sce, dims = 1:10)

umap_plot <- DimPlot(sce, reduction = "umap", label = TRUE) + ggtitle("UMAP Clusters")
ggsave("results/demo_umap_clusters.pdf", umap_plot, width = 8, height = 6)

# =============================================================================
# Step 6: QC violin plot
# =============================================================================
cat("Step 6: QC violin plot...\n")

p_ncount <- VlnPlot(sce, features = 'nCount_RNA', pt.size = 0.5, log = TRUE) +
  ggtitle("nCount_RNA (log)") + NoLegend()

p_nfeature <- VlnPlot(sce, features = 'nFeature_RNA', pt.size = 0.5, log = TRUE) +
  ggtitle("nFeature_RNA (log)") + NoLegend()

p_pctmt <- VlnPlot(sce, features = 'percent.mt', pt.size = 0.5) +
  ggtitle("Mitochondrial %") + NoLegend()

p_pctribo <- VlnPlot(sce, features = 'percent.ribo', pt.size = 0.5) +
  ggtitle("Ribosomal %") + NoLegend()

qc_plot <- plot_grid(p_ncount, p_nfeature, p_pctmt, p_pctribo, nrow = 2)
ggsave("results/demo_qc_violin.pdf", qc_plot, width = 12, height = 10)

# =============================================================================
# Step 7: Find markers
# =============================================================================
cat("Step 7: Finding marker genes...\n")

sce.markers <- FindAllMarkers(object = sce, test.use = "wilcox")
write.csv(sce.markers, "results/demo_markers.csv", row.names = FALSE)

# Top 10 markers per cluster
top10 <- sce.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)

# =============================================================================
# Step 8: Marker heatmap
# =============================================================================
cat("Step 8: Marker heatmap...\n")

dh <- DoHeatmap(sce, features = top10$gene, size = 3) +
  scale_fill_gradientn(colours = c("#045a8d", "white", "#a50f15")) +
  ggtitle("Top Marker Genes Heatmap")

ggsave("results/demo_marker_heatmap.pdf", dh, width = 10, height = 12)

# =============================================================================
# Step 9: Feature plots for cell type markers
# =============================================================================
cat("Step 9: Cell type feature plots...\n")

# Example markers for PBMC
cd3d <- FeaturePlot(sce, features = "CD3D") + ggtitle("CD3D (T cells)")
ms4a1 <- FeaturePlot(sce, features = "MS4A1") + ggtitle("MS4A1 (B cells)")
gnly <- FeaturePlot(sce, features = "GNLY") + ggtitle("GNLY (NK cells)")
cd14 <- FeaturePlot(sce, features = "CD14") + ggtitle("CD14 (Monocytes)")

feature_plot <- plot_grid(cd3d, ms4a1, gnly, cd14, nrow = 2)
ggsave("results/demo_feature_plots.pdf", feature_plot, width = 12, height = 10)

# =============================================================================
# Step 10: Save final object
# =============================================================================
cat("Step 10: Saving results...\n")

saveRDS(sce, "results/demo_seurat_object.rds")

cat("\nDemo pipeline completed!\n")
cat("All output files are in the 'results/' directory.\n")
