# scRNA-pipeline
This is a reproducible and modular pipeline for single-cell RNA-seq (scRNA-seq) data processing and downstream analysis, from raw count matrices to biological interpretation and visualization.
## ⚠️Reminders
> **Note 1:** This analysis is designed to run in R/RStudio environment.
> **Note 2:** Make sure R (≥4.0) and RStudio are installed before running the pipeline.
> **Note 3:** It is recommended to set up a project-specific R environment for reproducibility.
## Overview
## This pipeline performs:
### 1. Data Loading and Seurat Object Creation

This step loads raw count data and creates a Seurat object, the primary data structure for scRNA-seq analysis.

```bash
sce <- CreateSeuratObject(
  counts = Read10X("/public/home/fanyuxuan"),  
  project = "ccRCC01"                          
)
```
### 2. Quality Control Metrics Calculation

This step calculates quality metrics to identify and filter low-quality cells.

```bash
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
rb.genes <- rownames(sce)[grep("^RP[SL]", rownames(sce))]
C <- GetAssayData(object = sce, layer = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,]) / Matrix::colSums(C) * 100
sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
v <- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0)
print(v)
```

<img width="886" height="386" alt="image" src="https://github.com/user-attachments/assets/1ded03c9-a0eb-44bf-850b-98c78cc7f5dc" />

### 3. Cell Filtering

This step removes low-quality cells based on QC metrics.

```bash
sce <- subset(sce, subset = 
  nFeature_RNA > 200 &     
  nFeature_RNA < 7500 &    
  percent.mt < 15          
)
```
### 4. Data Normalization and Feature Selection

This step normalizes counts and identifies highly variable genes for downstream analysis.

```bash
sce <- NormalizeData(sce, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)
sce <- FindVariableFeatures(sce, 
  selection.method = "vst", 
  nfeatures = 2000
)
features <- VariableFeatures(object = sce)
sce <- ScaleData(sce, features = features)
```
### 5. Dimensionality Reduction (PCA)

This step performs linear dimensionality reduction to capture major sources of variation.

```bash
sce <- RunPCA(sce, 
  npcs = 50, 
  features = features
)
e <- ElbowPlot(sce, ndims = ncol(Embeddings(sce, "pca")))
e
```

<img width="889" height="398" alt="image" src="https://github.com/user-attachments/assets/85215b43-7cbc-48af-b126-5b6bc23172c3" />

### 6. Clustering and Non-linear Dimensionality Reduction

Constructs a shared nearest neighbor (SNN) graph from the dimensionality-reduced data.

```bash
sce <- FindNeighbors(sce, dims = 1:30)
```

Performs graph-based clustering to identify distinct cell populations.

```bash
sce <- FindClusters(sce, resolution = 0.8, algorithm = 3)
```

Performs t-distributed Stochastic Neighbor Embedding (t-SNE) for 2D visualization.

```bash
sce <- RunTSNE(object = sce, dims = 1:30, do.fast = TRUE)
```

Performs Uniform Manifold Approximation and Projection (UMAP) for 2D visualization.

```bash
sce <- RunUMAP(sce, dims = 1:30)
```

### 7. Doublet Detection and Removal

This step identifies and removes potential doublets (two cells in one droplet).

```bash
pcSelect =30
sce_double <- FastDoubletFinder(sce, 
  pcSelect = pcSelect,
  doublet.rate = 0.076,
  annotation = "seurat_clusters",
  pN_value = 0.25
)
sce <- subset(sce_double, Doublet == "Singlet")
```

### 8. Cell Cycle Scoring

This step assigns cell cycle phase scores based on expression of phase-specific genes.

```bash
s.genes <- cc.genes.updated.2019$s.genes    
g2m.genes <- cc.genes.updated.2019$g2m.genes 
sce <- CellCycleScoring(sce, 
  s.features = s.genes, 
  g2m.features = g2m.genes
)
```

### 9.  Marker Gene Identification

This step identifies genes that are differentially expressed in each cluster.

```bash
sce.markers <- FindAllMarkers(object = sce)
top10_cl_markers <- sce.markers %>% dplyr::group_by(cluster) %>% slice_max(avg_log2FC,n=10)
```

### 10.  Visualization and Custom Plotting Functions

Creates multi-panel violin plots for QC metrics

```bash
FastPlotVlnPlot <- function(object){
    p4 <- VlnPlot(object,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
      geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
    p5 <- VlnPlot(object,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
      geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
    p6 <- VlnPlot(object,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
      geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
      theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
    p7 <- VlnPlot(object,features='percent.ribo',ncol=1,pt.size=0,y.max = 100) +
      geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
      theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")  
    plot <- p4 + p5 + p6 + p7 + plot_layout(nrow = 4)
    return(plot)
  }
qc_plot <- FastPlotVlnPlot(sce)
  print(qc_plot) 
```

<img width="775" height="450" alt="image" src="https://github.com/user-attachments/assets/698edae8-649b-4e98-9b6d-413812e279d7" />

### 11. Marker Gene Heatmap Visualization

This step visualizes expression patterns of top marker genes across clusters.

```bash
dh <- DoHeatmap(sce, features = top10_cl_markers$gene, size = 3) + 
  scale_fill_gradientn(colours=c("#045a8d","white","#a50f15")) + NoLegend()
print(dh)
```

<img width="761" height="371" alt="image" src="https://github.com/user-attachments/assets/e987fb5f-007f-4335-9c05-6f3640cd7219" />

### 12. Cell Type Marker Gene Definition

This section defines canonical marker genes for major cell types.

```bash
cell_type_genes <- list(
  EPCAM = c("EPCAM","KRT19","KRT7","KRT18"), 
  cholangiocytes = c("FYXD2","TM4SF4", "ANXA4"), 
  hepatocytes = c("APOC3","FABP1","APOA1"), 
  fibroblasts = c("ACTA2" ,"COL1A2"), 
  endothelial = c("ENG","VWF"), 
  All_Immune = "PTPRC", 
  B_cell = c("BANK1","CD19","CD79A","CD79B","MS4A1"), 
  T_cell = c("CD3D","CD3E","CD8A","CD4","CD2"), 
  NK = c("GNLY","NKG7","FGFBP2", "KLRF1"), 
  Dendritic_cell = c("CLEC9A","CD1C","THBD","ITGAX","FLT3","IDO1","FCER1A","HLA-DQA1","LAMP3","CCR7","FSCN1"), 
  proliferative = c("MKI67","PCNA","TK1","TOP2A","TUBB","TYMS"), 
  Macrophage = c("INHBA","IL1RN","CCL4","NLRP3","EREG","IL1B","LYVE1","PLTP","SEPP1","C1QC","C1QA","APOE","GPNMB","SPP1","ISG15","VCAN","CX3CR1","PPARG"), 
  Macrophage1 = c("CD68","CD14","HLA-DRA","HLA-DRB1","APOE","MMP9","CTSK"), 
  Mast = c("KIT","TPSAB1","CPA3"), 
  Monocyte = c("CD14","FCN1","S100A9","S100A8","FCGR3A","LST1","LILRB2"), 
  pDC = c("LILRA4","UGCG","CLEC4C","GZMB","IL3RA"), 
  plasma = c("SDC1","CD38","SLAMF7","DERL3","MZB1","IGHA1","IGHA2") 
)
```

### 13. Gene Expression Feature Plots

This step visualizes expression of marker genes on UMAP plots.

```bash
for (cell_type in names(cell_type_genes)) {
  current_genes <- cell_type_genes[[cell_type]]
  if (length(current_genes) <= 1) {
    message(paste("⚠️  跳过单基因细胞类型：", cell_type))
    next
  }
  p <- FeaturePlot(
    object = sce,         
    features = current_genes, 
    ncol = 3,             
    cols = c("lightblue", "red"), 
    pt.size = 0.2,         
  ) + theme_bw()
  plot_height <- 4 * ceiling(length(current_genes) / 3)
   ggsave(
    filename = paste0("~/pipeline/", cell_type, "_FeaturePlot.pdf"),
    plot = p,
    width = 12,
    height = plot_height,
    dpi = 300,
    limitsize = FALSE
  )}
```

<img width="611" height="357" alt="image" src="https://github.com/user-attachments/assets/927063ed-28f9-44a3-940b-804cebf5e69b" />

<img width="645" height="380" alt="image" src="https://github.com/user-attachments/assets/78b0e39c-f3cb-4e7e-90b9-b2199d046917" />

<img width="658" height="329" alt="image" src="https://github.com/user-attachments/assets/329e077d-aecf-4b0b-b7f6-3e9cc7b36672" />

<img width="563" height="282" alt="image" src="https://github.com/user-attachments/assets/51d7c713-4508-4352-aa0b-2e5888dd7f99" />

<img width="546" height="273" alt="image" src="https://github.com/user-attachments/assets/aa308ac7-6abf-42b3-a52a-47d2ac59e077" />

<img width="625" height="314" alt="image" src="https://github.com/user-attachments/assets/b9fed2cb-ab0b-41f5-a885-18f5352d50c7" />

<img width="551" height="273" alt="image" src="https://github.com/user-attachments/assets/4913070b-2d14-42bb-8472-819c9416b061" />

### 14.  UMAP Clustering Visualization

This step creates publication-quality UMAP plots with cluster annotations.

```bash
p_umap <- dittoDimPlot(
  sce,
  reduction.use = "umap",
  var = "seurat_clusters",
  size = 1,
  do.label = TRUE,
  do.ellipse = TRUE,
  legend.size = 9,
  labels.size = 7,
  do.raster = T,
  raster.dpi = 300
) + theme(legend.text = element_text(face="bold",size = 16))
```

<img width="443" height="358" alt="image" src="https://github.com/user-attachments/assets/ba6d65dc-8b85-4586-b443-02b9710508c0" />












