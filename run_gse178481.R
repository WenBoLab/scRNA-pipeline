# =============================================================================
# Run scRNA-seq pipeline on real GSE178481 data (RCC-PR6-PTumor)
# =============================================================================

library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(harmony)

# Set paths
data_file <- "data/RCC-PR6-PTumor.count.csv"
outdir <- "results/gse178481"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Step 1: Read count matrix and create Seurat object
# =============================================================================
cat("Step 1: Reading count matrix...\n")

# Read CSV - first column is gene names
counts <- read.csv(data_file, row.names = 1, check.names = FALSE)
cat("Data dimensions:", nrow(counts), "genes x", ncol(counts), "cells\n")

# Convert to sparse matrix
counts <- as(as.matrix(counts), "dgCMatrix")

# Create Seurat object
sce <- CreateSeuratObject(counts = counts, project = "GSE178481_PR6_PTumor")
cat("Seurat object created:", ncol(sce), "cells,", nrow(sce), "genes\n")

# =============================================================================
# Step 2: Quality control metrics
# =============================================================================
cat("Step 2: Calculating QC metrics...\n")

sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")

rb.genes <- rownames(sce)[grep("^RP[SL]", rownames(sce))]
if (length(rb.genes) > 0) {
  C <- GetAssayData(object = sce, layer = "counts")
  percent.ribo <- Matrix::colSums(C[rb.genes, , drop = FALSE]) / Matrix::colSums(C) * 100
  sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
} else {
  sce[["percent.ribo"]] <- 0
}

# QC violin plot before filtering
p_qc_before <- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                       pt.size = 0.1, ncol = 2)
ggsave(file.path(outdir, "01_qc_violin_before_filtering.pdf"), p_qc_before, width = 12, height = 10)
cat("QC before filtering: saved\n")

# Filter cells
sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
cat("After QC filtering:", ncol(sce), "cells remain\n")

# =============================================================================
# Step 3: Normalization and variable features
# =============================================================================
cat("Step 3: Normalization and variable feature selection...\n")

sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)

# Top 10 variable features plot
top10 <- head(VariableFeatures(sce), 10)
p_hvg <- VariableFeaturePlot(sce)
p_hvg <- LabelPoints(plot = p_hvg, points = top10, repel = TRUE)
ggsave(file.path(outdir, "02_variable_features.pdf"), p_hvg, width = 10, height = 6)

# =============================================================================
# Step 4: Scaling and PCA
# =============================================================================
cat("Step 4: Scaling data and running PCA...\n")

all.genes <- rownames(sce)
sce <- ScaleData(sce, features = VariableFeatures(sce))
sce <- RunPCA(sce, npcs = 50, features = VariableFeatures(sce))

# Elbow plot
p_elbow <- ElbowPlot(sce, ndims = 50)
ggsave(file.path(outdir, "03_elbow_plot.pdf"), p_elbow, width = 10, height = 6)

# PCA plot
p_pca <- DimPlot(sce, reduction = "pca", dims = c(1, 2)) + ggtitle("PCA")
ggsave(file.path(outdir, "04_pca_plot.pdf"), p_pca, width = 8, height = 6)

# =============================================================================
# Step 5: Cell cycle scoring
# =============================================================================
cat("Step 5: Cell cycle scoring...\n")

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
sce <- CellCycleScoring(sce, s.features = s.genes, g2m.features = g2m.genes)

# =============================================================================
# Step 6: Clustering and UMAP
# =============================================================================
cat("Step 6: Clustering and UMAP...\n")

pcSelect <- 30

sce <- FindNeighbors(sce, dims = 1:pcSelect)
sce <- FindClusters(sce, resolution = 0.8, algorithm = 3)
sce <- RunUMAP(sce, dims = 1:pcSelect)
sce <- RunTSNE(sce, dims = 1:pcSelect)

cat("Number of clusters:", length(unique(Idents(sce))), "\n")

# UMAP cluster plot
p_umap <- DimPlot(sce, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("UMAP - GSE178481 RCC-PR6-PTumor")
ggsave(file.path(outdir, "05_umap_clusters.pdf"), p_umap, width = 10, height = 8)

# tSNE cluster plot
p_tsne <- DimPlot(sce, reduction = "tsne", label = TRUE, repel = TRUE) +
  ggtitle("tSNE - GSE178481 RCC-PR6-PTumor")
ggsave(file.path(outdir, "06_tsne_clusters.pdf"), p_tsne, width = 10, height = 8)

# =============================================================================
# Step 7: QC violin plot after clustering
# =============================================================================
cat("Step 7: QC violin plot across clusters...\n")

p_ncount <- VlnPlot(sce, features = 'nCount_RNA', pt.size = 0, log = TRUE) +
  ggtitle("nCount_RNA") + NoLegend()
p_nfeature <- VlnPlot(sce, features = 'nFeature_RNA', pt.size = 0, log = TRUE) +
  ggtitle("nFeature_RNA") + NoLegend()
p_pctmt <- VlnPlot(sce, features = 'percent.mt', pt.size = 0) +
  ggtitle("Mitochondrial %") + NoLegend()
p_pctribo <- VlnPlot(sce, features = 'percent.ribo', pt.size = 0) +
  ggtitle("Ribosomal %") + NoLegend()

p_qc_combined <- plot_grid(p_ncount, p_nfeature, p_pctmt, p_pctribo, nrow = 2)
ggsave(file.path(outdir, "07_qc_violin_by_cluster.pdf"), p_qc_combined, width = 14, height = 10)

# =============================================================================
# Step 8: Find marker genes
# =============================================================================
cat("Step 8: Finding marker genes...\n")

sce.markers <- FindAllMarkers(object = sce, test.use = "wilcox", only.pos = TRUE)
write.csv(sce.markers, file.path(outdir, "08_all_markers.csv"), row.names = FALSE)

# Top 10 markers per cluster
top10 <- sce.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

cat("Top markers found for", length(unique(sce.markers$cluster)), "clusters\n")

# =============================================================================
# Step 9: Marker heatmap
# =============================================================================
cat("Step 9: Marker heatmap...\n")

p_heatmap <- DoHeatmap(sce, features = top10$gene, size = 3) +
  scale_fill_gradientn(colours = c("#045a8d", "white", "#a50f15")) +
  ggtitle("Top 10 Marker Genes per Cluster")
ggsave(file.path(outdir, "09_marker_heatmap.pdf"), p_heatmap, width = 12, height = 16)

# =============================================================================
# Step 10: Cell type marker feature plots
# =============================================================================
cat("Step 10: Cell type feature plots...\n")

# Define cell type markers for ccRCC
markers_list <- list(
  Epithelial = c("EPCAM", "KRT19", "KRT7", "KRT18"),
  T_cells = c("CD3D", "CD3E", "CD8A", "CD4"),
  NK_cells = c("GNLY", "NKG7", "FGFBP2"),
  B_cells = c("CD19", "CD79A", "MS4A1"),
  Macrophage = c("CD68", "CD14", "APOE"),
  Endothelial = c("ENG", "VWF"),
  Fibroblast = c("ACTA2", "COL1A2"),
  Plasma = c("SDC1", "CD38", "MZB1")
)

# Check which markers are present
for (name in names(markers_list)) {
  genes <- markers_list[[name]]
  present <- genes[genes %in% rownames(sce)]
  if (length(present) > 0) {
    p <- FeaturePlot(sce, features = present, ncol = 2, order = TRUE)
    ggsave(file.path(outdir, paste0("10_", name, "_feature_plot.pdf")), p, width = 12, height = 8)
    cat("  Plotted:", name, "-", paste(present, collapse = ", "), "\n")
  }
}

# =============================================================================
# Step 11: Save Seurat object
# =============================================================================
cat("Step 11: Saving Seurat object...\n")

saveRDS(sce, file.path(outdir, "GSE178481_PR6_PTumor_seurat.rds"))

cat("\n========================================\n")
cat("Pipeline completed successfully!\n")
cat("Results saved in:", outdir, "\n")
cat("Final object:", ncol(sce), "cells,", nrow(sce), "genes\n")
cat("Number of clusters:", length(unique(Idents(sce))), "\n")
cat("========================================\n")
