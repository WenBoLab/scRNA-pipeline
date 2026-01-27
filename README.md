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

This step identifies cell clusters and visualizes them in 2D space.

```bash
sce <- FindNeighbors(sce, 
  reduction = "harmony", 
  dims = 1:30
)
sce <- FindClusters(sce, 
  resolution = 0.8,    
  algorithm = 3        
)
sce <- RunTSNE(sce, 
  reduction = "harmony", 
  dims = 1:30
)

sce <- RunUMAP(sce, 
  reduction = "harmony", 
  dims = 1:30
)
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

















