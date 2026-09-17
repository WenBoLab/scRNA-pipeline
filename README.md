# scRNA-seq Pipeline for GSE178481 (ccRCC)

A reproducible single-cell RNA-seq analysis pipeline for renal clear cell carcinoma (ccRCC) dataset **GSE178481**, built on [Seurat](https://satijalab.org/seurat/) and following standard best practices.

---

## Table of Contents

- [Overview](#overview)
- [Dataset](#dataset)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Steps (with actual output figures)](#pipeline-steps-with-actual-output-figures)
  - [Step 1: Create Seurat Object](#step-1-create-seurat-object)
  - [Step 2: Cell Quality Control](#step-2-cell-quality-control)
  - [Step 3: Normalization & Variable Features](#step-3-normalization--variable-features)
  - [Step 4: Scaling & PCA](#step-4-scaling--pca)
  - [Step 5: Clustering, UMAP & tSNE](#step-5-clustering-umap--tsne)
  - [Step 6: QC Violin Plot by Cluster](#step-6-qc-violin-plot-by-cluster)
  - [Step 7: Marker Gene Detection & Heatmap](#step-7-marker-gene-detection--heatmap)
  - [Step 8: Cell Type Marker Feature Plots](#step-8-cell-type-marker-feature-plots)
- [Function Reference](#function-reference)
- [Output Files](#output-files)
- [Files Overview (demo_run.R etc.)](#files-overview-demo_runr-etc)
- [Citation](#citation)

---

## Overview

This pipeline performs a complete single-cell RNA-seq analysis workflow:

1. **Data Import** — Load count matrix into a Seurat object
2. **Quality Control** — Filter low-quality cells using mitochondrial/ribosomal percentages and feature counts
3. **Normalization & Variable Features** — LogNormalize + top 2000 highly variable genes
4. **Dimensional Reduction** — PCA with ElbowPlot to select PCs
5. **Clustering** — KNN graph construction + community detection (SLM algorithm)
6. **Dimensional Reduction** — UMAP and tSNE embeddings
7. **Cell Cycle Scoring** — S/G2M phase scoring
8. **Marker Gene Detection** — FindAllMarkers to identify cluster-specific genes
9. **Visualization** — QC violin plots, marker heatmap, feature plots, UMAP

---

## Dataset

**GSE178481** — Single-cell RNA-seq of renal clear cell carcinoma

- **Source**: [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178481)
- **Species**: Human (*Homo sapiens*)
- **Tissue**: Kidney / Renal tumor
- **Assay**: 10x Genomics 3' scRNA-seq

### Download Instructions

1. Visit the [GSE178481 GEO page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178481)
2. Download the supplementary count matrix for your sample of interest (e.g., `GSM5392399_RCC-PR6-PTumor.count.csv.gz`)
3. Place it in `data/` and decompress:
   ```
   data/
   └── RCC-PR6-PTumor.count.csv
   ```
   The CSV is a gene × cell count matrix (first column = gene names, header row = cell barcodes).

> **Note:** The run in this README uses sample **RCC-PR6-PTumor** (GSM5392399). You can swap in any other sample from the dataset by changing the file path in `main.R`.

---

## Project Structure

```
scRNA-seq-pipeline/
├── R/
│   ├── Plus.library.R          # (optional reference) wrapper functions, split per file
│   ├── FastCreateSeurat.R
│   ├── FastSeuratCellQuality.R
│   ├── FastDoubletFinder.R
│   ├── FastSeuratRNA.R
│   └── FastPlotVlnPlot.R
├── data/                        # Raw data (not tracked in Git)
├── results/
│   └── gse178481/               # Example output figures (committed to Git)
├── main.R                       # ★ Self-contained main script: wrapper functions + actual run
├── demo_run.R                   # Tiny self-test on Seurat's pbmc_small (no download needed)
├── README.md                    # This file
└── .gitignore                   # Git ignore rules
```

> **`main.R` is fully self-contained.** All wrapper functions (`FastCreateSeurat`, `FastSeuratCellQuality`, `FastSeuratRNA`, `FastDoubletFinder`, `FastPlotVlnPlot`) are defined inside it, then called on GSE178481 data at the bottom. You do **not** need to `source()` anything from `R/` — that folder is kept only as a reference showing the functions split into one file each.

---

## Installation

### Requirements

- R >= 4.0
- R packages:
  - Seurat >= 5.0
  - Harmony
  - DoubletFinder
  - dplyr, cowplot, patchwork, ggplot2
  - scCustomize, dittoSeq
  - Matrix, presto

### Install Packages

```r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("multtest", "SingleCellExperiment"))

# Install CRAN packages
install.packages(c(
  "Seurat", "dplyr", "cowplot", "patchwork",
  "ggplot2", "Matrix", "harmony"
))

# Install DoubletFinder from GitHub
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
}

# Install scCustomize and dittoSeq
devtools::install_github("samuel-marsh/scCustomize")
BiocManager::install("dittoSeq")
```

---

## Quick Start

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/scRNA-seq-pipeline.git
   cd scRNA-seq-pipeline
   ```

2. Download a GSE178481 sample count matrix and place it in `data/`

3. Run the analysis:
   ```bash
   Rscript main.R
   ```

   Or run interactively in RStudio:
   ```r
   source("main.R")
   ```

---

## Pipeline Steps (with actual output figures)

All figures below were produced by running `run_gse178481.R` on the **RCC-PR6-PTumor** sample (3,577 raw cells → 3,555 filtered cells → 15 clusters).

### Step 1: Create Seurat Object

```r
# Read CSV count matrix (genes x cells)
counts <- read.csv("data/RCC-PR6-PTumor.count.csv", row.names = 1, check.names = FALSE)
counts <- as(as.matrix(counts), "dgCMatrix")

sce <- FastCreateSeurat(
  count = counts,
  project = "GSE178481_PR6_PTumor"
)
```

- Reads the gene × cell count matrix
- Converts to a sparse `dgCMatrix` to save memory
- Wraps in `CreateSeuratObject()`
- **Result**: 3,577 cells × 32,738 genes

---

### Step 2: Cell Quality Control

```r
qc_result <- FastSeuratCellQuality(
  obj = sce,
  species = "human",
  min.features = 200,
  max.features = 7500,
  percent.mt.num = 15
)
sce <- qc_result$sce
```

- Calculates mitochondrial gene percentage (`^MT-`)
- Calculates ribosomal gene percentage (`^RP[SL]`)
- Filters cells by: `nFeature_RNA` between 200–7500, `percent.mt < 15`
- **Result**: 3,555 cells retained (99.4%)

QC metrics before filtering (to inspect outliers):

![QC violin before filtering](results/gse178481/01_qc_violin_before_filtering.png)

---

### Step 3: Normalization & Variable Features

```r
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(sce), 10)
p_hvg <- VariableFeaturePlot(sce)
p_hvg <- LabelPoints(plot = p_hvg, points = top10, repel = TRUE)
```

- Log-normalizes counts with a scale factor of 10,000
- Selects the top 2,000 highly variable genes using the VST method

![Highly variable features](results/gse178481/02_variable_features.png)

---

### Step 4: Scaling & PCA

```r
sce <- ScaleData(sce, features = VariableFeatures(sce))
sce <- RunPCA(sce, npcs = 50, features = VariableFeatures(sce))

p_elbow <- ElbowPlot(sce, ndims = 50)
p_pca   <- DimPlot(sce, reduction = "pca", dims = c(1, 2))
```

- Scales the variable features (mean-centering + scaling to unit variance)
- Computes 50 principal components
- The Elbow plot helps choose how many PCs to keep (we use 30 downstream)

Elbow plot:

![Elbow plot](results/gse178481/03_elbow_plot.png)

PCA projection (PC1 vs PC2):

![PCA plot](results/gse178481/04_pca_plot.png)

---

### Step 5: Clustering, UMAP & tSNE

```r
pcSelect <- 30

sce <- FindNeighbors(sce, dims = 1:pcSelect)
sce <- FindClusters(sce, resolution = 0.8, algorithm = 3)   # SLM
sce <- RunUMAP(sce, dims = 1:pcSelect)
sce <- RunTSNE(sce, dims = 1:pcSelect)
```

- Builds a KNN graph on PCs 1–30
- Communities detected with the SLM algorithm at resolution 0.8
- Projects onto UMAP and tSNE embeddings
- **Result**: 15 clusters

UMAP (colored by cluster):

![UMAP clusters](results/gse178481/05_umap_clusters.png)

tSNE (colored by cluster):

![tSNE clusters](results/gse178481/06_tsne_clusters.png)

---

### Step 6: QC Violin Plot by Cluster

```r
p_ncount  <- VlnPlot(sce, features = 'nCount_RNA', pt.size = 0, log = TRUE)
p_nfeature<- VlnPlot(sce, features = 'nFeature_RNA', pt.size = 0, log = TRUE)
p_pctmt   <- VlnPlot(sce, features = 'percent.mt', pt.size = 0)
p_pctribo <- VlnPlot(sce, features = 'percent.ribo', pt.size = 0)

p_qc_combined <- plot_grid(p_ncount, p_nfeature, p_pctmt, p_pctribo, nrow = 2)
```

Checks that no cluster has suspiciously high mitochondrial content or extreme feature counts:

![QC violin by cluster](results/gse178481/07_qc_violin_by_cluster.png)

---

### Step 7: Marker Gene Detection & Heatmap

```r
sce.markers <- FindAllMarkers(object = sce, test.use = "wilcox", only.pos = TRUE)

top10 <- sce.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

p_heatmap <- DoHeatmap(sce, features = top10$gene, size = 3)
```

- Finds positive marker genes for every cluster with the Wilcoxon rank-sum test
- Takes the top 10 markers per cluster and plots a heatmap of scaled expression

![Marker heatmap](results/gse178481/09_marker_heatmap.png)

---

### Step 8: Cell Type Marker Feature Plots

Overlay canonical cell-type marker genes on the UMAP to interpret the clusters:

```r
FeaturePlot(sce, features = c("EPCAM", "KRT19", "KRT7", "KRT18"))  # epithelial
FeaturePlot(sce, features = c("CD3D", "CD3E", "CD8A", "CD4"))      # T cells
FeaturePlot(sce, features = c("GNLY", "NKG7", "FGFBP2"))           # NK cells
FeaturePlot(sce, features = c("CD19", "CD79A", "MS4A1"))           # B cells
FeaturePlot(sce, features = c("CD68", "CD14", "APOE"))             # macrophages
FeaturePlot(sce, features = c("ENG", "VWF"))                       # endothelial
FeaturePlot(sce, features = c("ACTA2", "COL1A2"))                   # fibroblasts
FeaturePlot(sce, features = c("SDC1", "CD38", "MZB1"))             # plasma cells
```

**Epithelial cells (EPCAM, KRT19, KRT7, KRT18):**

![Epithelial feature plot](results/gse178481/10_Epithelial_feature_plot.png)

**T cells (CD3D, CD3E, CD8A, CD4):**

![T cells feature plot](results/gse178481/10_T_cells_feature_plot.png)

**Macrophages (CD68, CD14, APOE):**

![Macrophage feature plot](results/gse178481/10_Macrophage_feature_plot.png)

**Endothelial cells (ENG, VWF):**

![Endothelial feature plot](results/gse178481/10_Endothelial_feature_plot.png)

**NK cells (GNLY, NKG7, FGFBP2):**

![NK cells feature plot](results/gse178481/10_NK_cells_feature_plot.png)

**B cells (CD19, CD79A, MS4A1):**

![B cells feature plot](results/gse178481/10_B_cells_feature_plot.png)

**Fibroblasts (ACTA2, COL1A2):**

![Fibroblast feature plot](results/gse178481/10_Fibroblast_feature_plot.png)

**Plasma cells (SDC1, CD38, MZB1):**

![Plasma cells feature plot](results/gse178481/10_Plasma_feature_plot.png)

---

## Function Reference

### `Plus.library(pkgs)`

Batch load and install R packages.

| Parameter | Type | Description |
|-----------|------|-------------|
| `pkgs` | character vector | Package names to load |

---

### `FastCreateSeurat(count, dir.name, project)`

Create a Seurat object from count matrix or 10X directory.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `count` | matrix | NULL | Count matrix (genes × cells) |
| `dir.name` | character | NULL | Path to 10X output directory |
| `project` | character | "project" | Project name |

**Returns**: Seurat object

---

### `FastSeuratCellQuality(obj, species, min.features, max.features, percent.mt.num)`

Calculate QC metrics and filter low-quality cells.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `obj` | Seurat | NULL | Input Seurat object |
| `species` | character | "human" | "human" or "mouse" |
| `min.features` | numeric | 0 | Min genes per cell |
| `max.features` | numeric | 3e7 | Max genes per cell |
| `percent.mt.num` | numeric | 100 | Max mitochondrial % |

**Returns**: List with `initial_sce` (unfiltered) and `sce` (filtered)

---

### `FastSeuratRNA(obj, species, ...)`

Main single-cell RNA-seq analysis pipeline.

| Key Parameter | Type | Default | Description |
|---------------|------|---------|-------------|
| `pcSelect` | numeric | 30 | Number of PCs to use |
| `nfeatures` | numeric | 2000 | Number of HVGs |
| `resolution` | numeric | 0.5 | Clustering resolution |
| `harmony` | logical | FALSE | Use Harmony batch correction |
| `doublet` | logical | FALSE | Run doublet detection |
| `cellCycle` | logical | TRUE | Score cell cycle phases |
| `isMarkers` | logical | TRUE | Find cluster markers |
| `algorithm` | numeric | 4 | Clustering algorithm (1=Louvain, 3=SLM, 4=Leiden) |

**Returns**: List with:
- `double_sce` — Doublet detection result
- `sce` — Processed Seurat object
- `gene.markers` — Marker gene data frame
- `features` — Variable features used
- `e` — ElbowPlot ggplot object

---

### `FastDoubletFinder(obj, pcSelect, doublet.rate, ...)`

Detect doublets using DoubletFinder.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pcSelect` | numeric | 30 | Number of PCs |
| `doublet.rate` | numeric | 0.076 | Expected doublet rate |
| `pN_value` | numeric | 0.25 | Synthetic doublet proportion |
| `sct` | logical | TRUE | Whether SCTransform was used |

**Returns**: Seurat object with `Doublet` column in meta.data

---

### `FastPlotVlnPlot(object)`

Combined 4-panel QC violin plot.

| Parameter | Type | Description |
|-----------|------|-------------|
| `object` | Seurat | Seurat object with clusters |

**Returns**: Combined ggplot object

---

## Output Files

After running `main.R`, the `results/gse178481/` directory contains:

| File | Description |
|------|-------------|
| `01_qc_violin_before_filtering.pdf/png` | QC metrics before cell filtering |
| `02_variable_features.pdf/png` | Top 10 highly variable genes |
| `03_elbow_plot.pdf/png` | Elbow plot for PC selection |
| `04_pca_plot.pdf/png` | PCA (PC1 vs PC2) |
| `05_umap_clusters.pdf/png` | UMAP cluster plot |
| `06_tsne_clusters.pdf/png` | tSNE cluster plot |
| `07_qc_violin_by_cluster.pdf/png` | QC violin plots per cluster |
| `08_all_markers.csv` | All cluster marker genes (not in Git) |
| `09_marker_heatmap.pdf/png` | Top 10 markers per cluster heatmap |
| `10_*.pdf/png` | Cell type marker feature plots |
| `GSE178481_PR6_PTumor_seurat.rds` | Final Seurat object (not in Git) |

---

## Files Overview

There are two R scripts at the top level:

| File | What it is | When to use it |
|------|------------|----------------|
| **`main.R`** | **Self-contained.** Part 1 defines all wrapper functions; Part 2 calls them on GSE178481 data to produce every figure in this README. No `source()` needed. | **Use this one.** It is both the library and the runnable script. |
| **`demo_run.R`** | A tiny self-contained demo that uses Seurat's built-in `pbmc_small` (80 cells, no download required). It skips data download and is only meant to verify that Seurat installs and runs correctly. | Quick smoke test after installing packages. |

The `R/` folder contains the same wrapper functions split into one file each, for reference. `main.R` writes them out inline so the code is self-contained and easy to run.

---

## Reproducibility Notes

- R version: >= 4.0
- Seurat version: >= 5.0 (uses `layer` syntax for count data)
- Random seed: Set by default in most Seurat functions
- Hardware: Minimum 8 GB RAM recommended; 16 GB+ for large datasets (multiple samples)

---

## Citation

If you use this pipeline in your research, please cite:

1. Stuart T, Butler A, Hoffman P, et al. **Seurat: an open-source R package for fast and easy analysis of single-cell RNA-seq data.** *Genome Biology*. 2021.
2. Korsunsky I, Millard N, Fan J, et al. **Fast, sensitive and accurate integration of single-cell data with Harmony.** *Nature Methods*. 2019.
3. McGinnis CS, Murrow LM, Gartner ZJ. **DoubletFinder: Doublet detection in single-cell RNA sequencing data using artificial nearest neighbors.** *Cell Systems*. 2019.
4. **GSE178481** dataset — Alchahin AM, Mei S, et al. *Nat Commun*. 2022. PMID: 36180422.

---

## License

MIT License — feel free to use and modify for your own projects.
