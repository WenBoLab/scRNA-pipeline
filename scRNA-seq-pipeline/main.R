# =============================================================================
# scRNA-seq Pipeline for GSE178481
# Renal Clear Cell Carcinoma (ccRCC) Single-Cell Analysis
# =============================================================================

# Set working directory (adjust as needed)
# setwd("/path/to/scRNA-seq-pipeline")

# Load utility functions
source("R/Plus.library.R")
source("R/FastCreateSeurat.R")
source("R/FastSeuratCellQuality.R")
source("R/FastDoubletFinder.R")
source("R/FastSeuratRNA.R")
source("R/FastPlotVlnPlot.R")

# Load required packages
Plus.library(c(
  "Seurat", "dplyr", "cowplot", "patchwork", "harmony",
  "scCustomize", "dittoSeq", "ggplot2", "Matrix"
))

# =============================================================================
# Step 1: Download and prepare GSE178481 data
# =============================================================================
# GSE178481 is a ccRCC single-cell RNA-seq dataset from GEO.
#
# Option 1: Download from GEO using GEOquery (requires GEOquery)
#   BiocManager::install("GEOquery")
#   library(GEOquery)
#   getGEO(filename = "data/GSE178481_series_matrix.txt.gz")
#
# Option 2: Download 10X data directly from GEO supplementary files
#   Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178481
#   Download the filtered_feature_bc_matrix.tar.gz file(s) and extract to data/
#
# Expected directory structure:
#   data/
#   └── GSE178481/
#       ├── barcodes.tsv.gz
#       ├── features.tsv.gz
#       └── matrix.mtx.gz

# Set data path
data_dir <- "data/GSE178481"

# Create data directory if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# =============================================================================
# Step 2: Create Seurat object
# =============================================================================
cat("Step 2: Creating Seurat object...\n")

sce <- FastCreateSeurat(
  dir.name = data_dir,
  project = "GSE178481"
)

cat("Initial object:", ncol(sce), "cells and", nrow(sce), "genes\n")

# =============================================================================
# Step 3: Cell quality control
# =============================================================================
cat("Step 3: Running cell quality control...\n")

qc_result <- FastSeuratCellQuality(
  obj = sce,
  species = "human",
  min.features = 200,
  max.features = 7500,
  percent.mt.num = 15
)

sce <- qc_result$sce
cat("After QC:", ncol(sce), "cells remain\n")

# =============================================================================
# Step 4: Run full Seurat RNA analysis pipeline
# =============================================================================
cat("Step 4: Running Seurat RNA analysis pipeline...\n")

pcSelect <- 30

pipeline_result <- FastSeuratRNA(
  obj = sce,
  species = "human",
  plot = TRUE,
  pcSelect = pcSelect,
  nfeatures = 2000,
  resolution = 0.8,
  harmony = TRUE,
  harmony_by = "orig.ident",
  doublet = TRUE,
  cellCycle = TRUE,
  isMarkers = TRUE,
  algorithm = 3,
  names = "GSE178481",
  outdir = "results/doublet"
)

sce <- pipeline_result$sce
sce.markers <- pipeline_result$gene.markers

cat("After doublet removal:", ncol(sce), "cells remain\n")

# =============================================================================
# Step 5: QC violin plot across clusters
# =============================================================================
cat("Step 5: Generating QC violin plot...\n")

p_qc <- FastPlotVlnPlot(sce)
ggsave("results/qc_violin_plot.pdf", p_qc, width = 8, height = 12)

# =============================================================================
# Step 6: Marker gene heatmap
# =============================================================================
cat("Step 6: Generating marker gene heatmap...\n")

top10_cl_markers <- sce.markers %>%
  dplyr::group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

dh <- DoHeatmap(sce, features = top10_cl_markers$gene) +
  scale_fill_gradientn(colours = c("#045a8d", "white", "#a50f15")) +
  NoLegend()

ggsave("results/marker_heatmap.pdf", dh, width = 12, height = 16)

# =============================================================================
# Step 7: Cell type marker gene feature plots
# =============================================================================
cat("Step 7: Generating cell type marker feature plots...\n")

# Define cell type marker genes
EPCAM <- c("EPCAM", "KRT19", "KRT7", "KRT18")          # Epithelial cells
fibroblasts <- c("ACTA2", "COL1A2")                      # Fibroblasts
endothelial <- c("ENG", "VWF")                           # Endothelial cells
B.cell <- c("BANK1", "CD19", "CD79A", "MS4A1")           # B cells
T.cell <- c("CD3D", "CD3E", "CD8A", "CD4")               # T cells
NK <- c("GNLY", "NKG7", "FGFBP2")                        # NK cells
Macrophage <- c("CD68", "CD14", "APOE", "GPNMB")         # Macrophages
Mast <- c("KIT", "TPSAB1", "CPA3")                       # Mast cells
Monocyte <- c("CD14", "FCN1", "S100A9")                  # Monocytes
plasma <- c("SDC1", "CD38", "MZB1")                      # Plasma cells
proliferative <- c("MKI67", "PCNA", "TOP2A")             # Proliferating cells

# Feature plot for epithelial marker
p_epcam <- FeaturePlot_scCustom(object = sce, features = EPCAM, order = FALSE)
ggsave("results/epcam_feature_plot.pdf", p_epcam, width = 10, height = 8)

# =============================================================================
# Step 8: UMAP cluster plot with labels
# =============================================================================
cat("Step 8: Generating UMAP cluster plot...\n")

p_umap <- dittoDimPlot(
  sce,
  reduction.use = "umap",
  var = "seurat_clusters",
  size = 1,
  do.label = TRUE,
  do.ellipse = TRUE,
  legend.size = 9,
  shape.legend.size = 9,
  labels.size = 7,
  do.raster = TRUE,
  raster.dpi = 500
) +
  theme(legend.text = element_text(face = "bold", size = 16))

ggsave("results/umap_clusters.pdf", p_umap, width = 10, height = 8)

# =============================================================================
# Step 9: Save final Seurat object
# =============================================================================
cat("Step 9: Saving final Seurat object...\n")

saveRDS(sce, file = "results/GSE178481_seurat_final.rds")
write.csv(sce.markers, file = "results/GSE178481_markers.csv", row.names = FALSE)

cat("\nPipeline completed successfully!\n")
cat("Results saved in the 'results/' directory.\n")
