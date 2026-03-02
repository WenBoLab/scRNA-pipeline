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

## When analyzing multiple cell groups with RNA-seq
### 1.Data Loading and Seurat Object Creation

Load raw expression matrix and create a Seurat object, the primary data structure for single-cell RNA-seq analysis.

```bash
suppressPackageStartupMessages({
  library(Seurat)
  library(stringr) 
  library(readr)
  library(ggplot2)
  library(harmony)
  library(dplyr)
  library(clustree)
})
path <- file.path("~", "test_sclc", "GSE281740_sclc_cellranger (1).csv.gz")
df <- read_csv(file = path)
expr_matrix <- as.data.frame(df)
rownames(expr_matrix) <- expr_matrix$feature_name
expr_matrix <- expr_matrix[, -1]          
expr_matrix <- as.matrix(expr_matrix)
sclc <- CreateSeuratObject(
  counts = expr_matrix, 
  project = "SCLC"
)
sclc@meta.data$orig.ident <- stringr::str_split(
  rownames(sclc@meta.data), "-", simplify = TRUE
)[, 3]
```

### 2.Quality Control Metrics Calculation

Calculate quality metrics such as mitochondrial and ribosomal gene percentages to identify low-quality cells.

```bash
sclc[["percent.mt"]] <- PercentageFeatureSet(sclc, pattern = "^MT-")
rb.genes <- rownames(sclc)[grep("^RP[SL]", rownames(sclc))]
C <- GetAssayData(object = sclc, layer = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes, ]) / Matrix::colSums(C) * 100
sclc <- AddMetaData(sclc, percent.ribo, col.name = "percent.ribo")
v1 <- VlnPlot(sclc, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
              pt.size = 0)
print(v1)
v2 <- VlnPlot(sclc, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
              group.by = 'orig.ident', 
              pt.size = 0)
print(v2)

<img width="883" height="388" alt="image" src="https://github.com/user-attachments/assets/51d6df16-ee19-498a-af35-53a3a43417c8" />

<img width="855" height="385" alt="image" src="https://github.com/user-attachments/assets/3742b4fc-bee7-4ba0-be90-8fe35c5c1a2e" />

```

### 3.Cell Filtering

Filter low-quality cells based on QC metrics.

```bash
sclc <- subset(sclc, 
               subset = nFeature_RNA > 200 &    
                        nFeature_RNA < 7500 &   
                        percent.mt < 15)       
```

### 4.Data Normalization and Feature Selection

Normalize count data and identify highly variable genes for downstream analysis.

```bash
sclc <- NormalizeData(sclc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)
sclc <- FindVariableFeatures(sclc, 
                             selection.method = "vst", 
                             nfeatures = 2000)
features <- VariableFeatures(object = sclc)
sclc <- ScaleData(sclc, features = features)
```

### 5.Dimensionality Reduction (PCA)

Perform linear dimensionality reduction to capture major sources of variation.

```bash
sclc <- RunPCA(sclc, npcs = 50, features = features)
e <- ElbowPlot(sclc, ndims = ncol(Embeddings(sclc, "pca")))
print(e)
```

<img width="867" height="399" alt="image" src="https://github.com/user-attachments/assets/c0ec746e-8fac-4bbb-a7ce-872359f9be69" />

### 6.Batch Effect Correction (Harmony)

Correct batch effects between samples using Harmony algorithm.

```bash
sclc <- RunHarmony(sclc, 
                   group.by.vars = "orig.ident",  # Correct by original sample ID
                   dims.use = 1:30, 
                   max_iter = 50)
```

### 7.Clustering and Non-linear Dimensionality Reduction

Build SNN graph

```bash
sclc <- FindNeighbors(sclc, reduction = "harmony", dims = 1:30)
```

Graph-based clustering

```bash
sclc <- FindClusters(sclc, resolution = 1, algorithm = 3)
```

t-SNE dimensionality reduction

```bash
sclc <- RunTSNE(object = sclc, reduction = "harmony", dims = 1:30)
```

UMAP dimensionality reduction

```bash
sclc <- RunUMAP(sclc, reduction = "harmony", dims = 1:30)
```

### 8.Clustering Results Visualization

Visualize UMAP dimensionality reduction results, colored by cluster and sample.

```bash
p1 <- DimPlot(sclc, 
              reduction = "umap", 
              group.by = "seurat_clusters", 
              label = TRUE, 
              repel = TRUE) + 
  ggtitle("Single Cell Clustering UMAP")
ggsave("~/test_sclc/Single_Cell_Clustering_UMAP.pdf", p1, width = 10, height = 8)
p2 <- DimPlot(sclc, 
              reduction = "umap", 
              group.by = "orig.ident", 
              label = TRUE, 
              repel = TRUE) + 
  ggtitle("Single Cell Clustering UMAP (Colored by Sample)")
ggsave("~/test_sclc/Single_Cell_Clustering_Sample_UMAP.pdf", p2, width = 10, height = 8)
```

<img width="1003" height="788" alt="image" src="https://github.com/user-attachments/assets/12acdb64-d278-4948-a875-8d045fa1247b" />

<img width="935" height="736" alt="image" src="https://github.com/user-attachments/assets/afdc58b5-6d44-40c5-bee3-d419691eaf7b" />

### 9. Marker Gene Definition

Define canonical marker gene sets for major cell types.

```bash
sclc_markers <- list(
  EPCAM = c("EPCAM", "KRT19", "KRT7", "KRT18"),          
  fibroblasts = c("ACTA2", "COL1A2"),                      
  endothelial = c("ENG", "VWF"),                            
  All_Immune = "PTPRC",                                    
  B_cell = c("BANK1", "CD19", "CD79A", "CD79B", "MS4A1"),   
  T_cell = c("CD3D", "CD3E", "CD8A", "CD4", "CD2"),         
  NK = c("GNLY", "NKG7", "FGFBP2", "KLRF1"),                
  Dendritic_cell = c("CLEC9A", "CD1C", "THBD", "ITGAX", "FLT3", 
                      "IDO1", "FCER1A", "HLA-DQA1", "LAMP3", "CCR7", "FSCN1"), 
  proliferative = c("MKI67", "PCNA", "TK1", "TOP2A", "TUBB", "TYMS"),
  Myeloid = c("CD68", "CD14", "HLA-DRA", "HLA-DRB1", "APOE", "MMP9", 
               "CTSK", "FCN1", "S100A9", "S100A8", "FCGR3A", "LST1", "LILRB2"),
  Mast = c("KIT", "TPSAB1", "CPA3"),                        
  pDC = c("LILRA4", "UGCG", "CLEC4C", "GZMB", "IL3RA"),      
  plasma = c("SDC1", "CD38", "SLAMF7", "DERL3", "MZB1", "IGHA1", "IGHA2") 
```

### 10.Marker Gene DotPlot Visualization

Create DotPlot to visualize expression patterns of marker genes across clusters.

```bash
cat(" -> Generating SCLC marker gene DotPlot...\n")
p_dot <- DotPlot(sclc, 
                 features = unique(unlist(sclc_markers)), 
                 assay = "RNA",                             
                 group.by = "seurat_clusters") + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  coord_flip()                                                
ggsave("~/test_sclc/SCLC_Detailed_Marker_Gene_DotPlot.pdf", 
       p_dot, width = 16, height = 20)
```

<img width="624" height="778" alt="image" src="https://github.com/user-attachments/assets/9cc403db-5255-4a23-a399-7ab616adfb76" />

### 11. Gene Expression Feature Plots

Visualize expression of marker genes on UMAP plots.

```bash
reduction_type <- "umap"
for (cell_type in names(sclc_markers)) {
  # Skip single-gene cell types to avoid too few subplots
  if (length(sclc_markers[[cell_type]]) > 1) {
    p <- FeaturePlot(sclc, 
                     features = sclc_markers[[cell_type]],
                     reduction = reduction_type,
                     ncol = 3,                              
                     cols = c("lightblue", "red"))
    ggsave(paste0("~/test/SCLC_", cell_type, "_FeaturePlot.pdf"),
           p, width = 12, height = 4 * ceiling(length(sclc_markers[[cell_type]])/3))
  }
}
```
<img width="1169" height="549" alt="image" src="https://github.com/user-attachments/assets/12e340d2-f99b-4127-8a95-a215f8e1b71e" />

<img width="1103" height="723" alt="image" src="https://github.com/user-attachments/assets/7ca1a2a0-4f1b-45dd-993e-86bf3be33b6a" />

<img width="1086" height="320" alt="image" src="https://github.com/user-attachments/assets/51d391b5-39aa-440c-bec2-6723ad114ad2" />

<img width="1026" height="732" alt="image" src="https://github.com/user-attachments/assets/82fc7f92-f93b-4cbc-94bd-2a4e9421542d" />

### 12.Cell Type Annotation

Annotate clusters based on marker gene expression patterns.

```bash
cluster_annotation <- c(
  "0" = "EPCAM", "1" = "EPCAM", "2" = "T_cell", "3" = "proliferative",
  "4" = "proliferative", "5" = "T_cell", "6" = "proliferative", "7" = "Myeloid",
  "8" = "proliferative", "9" = "EPCAM", "10" = "EPCAM", "11" = "Myeloid",
  "12" = "NK", "13" = "proliferative", "14" = "EPCAM", "15" = "proliferative",
  "16" = "plasma", "17" = "fibroblasts", "18" = "B_cell", "19" = "EPCAM",
  "20" = "EPCAM", "21" = "endothelial", "22" = "EPCAM", "23" = "EPCAM",
  "24" = "T_cell", "25" = "proliferative", "26" = "EPCAM", "27" = "proliferative",
  "28" = "EPCAM", "29" = "EPCAM", "30" = "Mast", "31" = "EPCAM",
  "32" = "EPCAM", "33" = "proliferative"
)
sclc$celltype <- as.character(sclc$seurat_clusters)
for (cluster in names(cluster_annotation)) {
  if (cluster %in% sclc$seurat_clusters) {
    sclc$celltype[sclc$seurat_clusters == cluster] <- cluster_annotation[cluster]
  }
}
sclc$celltype[!sclc$celltype %in% cluster_annotation] <- "Unknown"
```

### 13.Annotated UMAP Visualization

Create publication-quality UMAP plots with cell type annotations.

```bash
Biocols <- c('#F1BB72', '#F3B1A0', '#E95C59', '#E59CC4', '#F7F398', '#E5D2DD',
             '#57C3F3', '#BD956A', '#8C549C', '#476D87', '#23452F', '#585658',
             '#D6E7A3', '#CCC9E6', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA',
             '#58A4C3', '#E4C755', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3',
             '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',
             '#625D9E', '#68A180', '#3A6963', '#968175', '#AB3282', '#53A85F')
p <- DimPlot(sclc, 
             reduction = reduction_type,
             label = TRUE, 
             repel = TRUE,                       
             group.by = "celltype", 
             cols = Biocols) +
  ggtitle("SCLC Annotated Cell Types (UMAP)") +
  theme(plot.title = element_text(hjust = 0.5))   # Center title
ggsave("~/test_sclc/SCLC_annotated_umap.pdf", p, width = 8, height = 6)
```

<img width="1012" height="754" alt="image" src="https://github.com/user-attachments/assets/084e4680-5e9d-4153-b95f-c698edd68063" />
























































