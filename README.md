# scRNA-seq Pipeline for GSE178481 (ccRCC)

A reproducible single-cell RNA-seq analysis pipeline for renal clear cell carcinoma (ccRCC) dataset **GSE178481**, built on [Seurat](https://satijalab.org/seurat/) and following standard best practices.

---

## Table of Contents

- [Overview](#overview)
- [Dataset](#dataset)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Steps](#pipeline-steps)
- [Function Reference](#function-reference)
- [Output Files](#output-files)
- [Expected Results](#expected-results)
- [Citation](#citation)

---

## Overview

This pipeline performs a complete single-cell RNA-seq analysis workflow:

1. **Data Import** — Load 10X Genomics data into a Seurat object
2. **Quality Control** — Filter low-quality cells using mitochondrial/ribosomal percentages and feature counts
3. **Normalization & Variable Features** — LogNormalize + top 2000 highly variable genes
4. **Dimensional Reduction** — PCA with ElbowPlot to select PCs
5. **Batch Correction** — Harmony integration across samples
6. **Clustering** — KNN graph construction + community detection (SLM algorithm)
7. **Dimensional Reduction** — UMAP and tSNE embeddings
8. **Doublet Detection** — DoubletFinder to remove doublets
9. **Cell Cycle Scoring** — S/G2M phase scoring
10. **Marker Gene Detection** — FindAllMarkers to identify cluster-specific genes
11. **Visualization** — QC violin plots, marker heatmap, feature plots, UMAP

---

## Dataset

**GSE178481** — Single-cell RNA-seq of renal clear cell carcinoma

- **Source**: [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178481)
- **Species**: Human (*Homo sapiens*)
- **Tissue**: Kidney / Renal tumor
- **Assay**: 10x Genomics 3' scRNA-seq

### Download Instructions

1. Visit the [GSE178481 GEO page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178481)
2. Download the supplementary 10X matrix files (e.g., `filtered_feature_bc_matrix.tar.gz`)
3. Extract to `data/GSE178481/` so the directory contains:
   ```
   data/GSE178481/
   ├── barcodes.tsv.gz
   ├── features.tsv.gz
   └── matrix.mtx.gz
   ```

---

## Project Structure

```
scRNA-seq-pipeline/
├── R/
│   ├── Plus.library.R          # Utility: batch package loading
│   ├── FastCreateSeurat.R      # Step 1: Create Seurat object
│   ├── FastSeuratCellQuality.R  # Step 2: Cell QC and filtering
│   ├── FastDoubletFinder.R      # Step 3: Doublet detection
│   ├── FastSeuratRNA.R          # Step 4: Main analysis (PCA/clustering/UMAP)
│   └── FastPlotVlnPlot.R        # Step 5: Combined QC violin plots
├── data/                        # Raw data (not tracked in Git)
├── results/                     # Output files (not tracked in Git)
├── docs/                        # Additional documentation
├── main.R                       # Main script: run the full pipeline
├── README.md                    # This file
└── .gitignore                   # Git ignore rules
```

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

2. Download GSE178481 data and place it in `data/GSE178481/`

3. Run the main script:
   ```r
   Rscript main.R
   ```

   Or run interactively in RStudio:
   ```r
   source("main.R")
   ```

---

## Pipeline Steps

### Step 1: Create Seurat Object

```r
sce <- FastCreateSeurat(
  dir.name = "data/GSE178481",
  project = "GSE178481"
)
```

- Reads 10X output using `Seurat::Read10X()`
- Wraps in `CreateSeuratObject()`
- Stores project name in object metadata

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
- Filters cells by:
  - `nFeature_RNA > 200` (detected genes)
  - `nFeature_RNA < 7500`
  - `percent.mt < 15` (mitochondrial content < 15%)

### Step 3: Main Analysis Pipeline

```r
result <- FastSeuratRNA(
  obj = sce,
  species = "human",
  pcSelect = 30,
  nfeatures = 2000,
  resolution = 0.8,
  harmony = TRUE,
  doublet = TRUE,
  cellCycle = TRUE,
  isMarkers = TRUE,
  algorithm = 3
)
sce <- result$sce
```

This function runs:
1. **Normalization** — `LogNormalize` with scale factor 10,000
2. **Variable Features** — Top 2000 HVGs using VST method
3. **Scaling** — Scale data for PCA
4. **PCA** — Compute 50 PCs, ElbowPlot for visualization
5. **Cell Cycle Scoring** — S and G2M phase scores using Seurat's built-in gene sets
6. **Harmony Batch Correction** — Correct for sample-level batch effects
7. **Nearest Neighbor Graph** — Find neighbors on Harmony embeddings (PCs 1:30)
8. **Clustering** — SLM algorithm at resolution 0.8
9. **tSNE & UMAP** — Non-linear dimensional reduction
10. **Doublet Removal** — DoubletFinder with 7.6% expected doublet rate
11. **Marker Detection** — FindAllMarkers with Wilcoxon rank-sum test

### Step 4: QC Visualization

```r
p_qc <- FastPlotVlnPlot(sce)
ggsave("results/qc_violin_plot.pdf", p_qc, width = 8, height = 12)
```

Generates a 4-panel violin plot showing:
- `nCount_RNA` (log scale)
- `nFeature_RNA` (log scale)
- `percent.mt`
- `percent.ribo`

Across all clusters, with boxplot overlays.

### Step 5: Marker Heatmap

```r
top10_cl_markers <- sce.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

dh <- DoHeatmap(sce, features = top10_cl_markers$gene)
```

Shows the top 10 marker genes per cluster, colored by scaled expression.

### Step 6: UMAP Visualization

```r
p_umap <- dittoDimPlot(
  sce,
  reduction.use = "umap",
  var = "seurat_clusters",
  do.label = TRUE,
  do.ellipse = TRUE
)
```

UMAP plot colored by cluster, with labels and ellipses.

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

After running `main.R`, the `results/` directory will contain:

| File | Description |
|------|-------------|
| `GSE178481_seurat_final.rds` | Final processed Seurat object |
| `GSE178481_markers.csv` | All cluster marker genes |
| `GSE178481_pcaplot.pdf` | PCA plot |
| `GSE178481_tsneplot.pdf` | tSNE cluster plot |
| `GSE178481_umapplot.pdf` | UMAP cluster plot |
| `qc_violin_plot.pdf` | QC violin plots across clusters |
| `marker_heatmap.pdf` | Top 10 markers per cluster heatmap |
| `epcam_feature_plot.pdf` | EPCAM epithelial marker feature plot |
| `umap_clusters.pdf` | Labeled UMAP with ellipses |

---

## Actual Results (GSE178481 RCC-PR6-PTumor)

### Cell Count Summary

| Stage | Cells |
|-------|-------|
| Raw input | 3,577 cells |
| After QC filtering (nFeature 200-7500, mt < 15%) | 3,555 cells (99.4% retained) |
| Final clusters | 15 clusters |

### UMAP Clustering

**UMAP embedding showing 15 clusters:**

<img width="1139" height="855" alt="image" src="https://github.com/user-attachments/assets/23f9f054-ff77-40f6-8b36-da3914666092" />


### tSNE Clustering

<img width="1131" height="863" alt="image" src="https://github.com/user-attachments/assets/9613606e-b611-487a-bdeb-06bf93aebc9c" />


### QC Violin Plot (by cluster)

4-panel QC metrics across 15 clusters:
- nCount_RNA (log scale)
- nFeature_RNA (log scale)
- Mitochondrial percentage
- Ribosomal percentage

<img width="1138" height="805" alt="image" src="https://github.com/user-attachments/assets/b2c2556a-e954-4fdd-bb2e-87cb2f71deb1" />


### Marker Gene Heatmap

Top 10 differentially expressed genes per cluster:

<img width="1126" height="1508" alt="image" src="https://github.com/user-attachments/assets/12ba615d-3826-4ba5-aea8-713d22a57e51" />

### Cell Type Marker Feature Plots

**Epithelial cells (EPCAM, KRT19, KRT7, KRT18):**

<img width="1228" height="818" alt="image" src="https://github.com/user-attachments/assets/28879586-34c9-4281-861d-6665d9888caf" />

**T cells (CD3D, CD3E, CD8A, CD4):**

<img width="1230" height="820" alt="image" src="https://github.com/user-attachments/assets/82240c7b-9350-4155-ae69-450295f35dc6" />


**Macrophages (CD68, CD14, APOE):**

<img width="1231" height="819" alt="image" src="https://github.com/user-attachments/assets/e8070a99-2f7c-4e1a-94cf-a70812eb0465" />

**Endothelial cells (ENG, VWF):**

<img width="1230" height="819" alt="image" src="https://github.com/user-attachments/assets/c56aab06-d973-4e98-8492-42011649c939" />

**NK cells (GNLY, NKG7, FGFBP2):**

<img width="1229" height="820" alt="image" src="https://github.com/user-attachments/assets/1a48386c-b089-4049-8c10-60b320618533" />


**B cells (CD19, CD79A, MS4A1):**

<img width="1230" height="819" alt="image" src="https://github.com/user-attachments/assets/6ac9386b-21c5-4916-82f0-111fe7e913a7" />

**Fibroblasts (ACTA2, COL1A2):**

<img width="1228" height="819" alt="image" src="https://github.com/user-attachments/assets/9b5030a6-ed48-4450-943d-8ed3c57f0535" />

**Plasma cells (SDC1, CD38, MZB1):**

<img width="1228" height="820" alt="image" src="https://github.com/user-attachments/assets/2d7a57e4-d85f-41d2-ad99-9ebbbaa13556" />


### PCA & Dimensionality Reduction

<img width="763" height="574" alt="image" src="https://github.com/user-attachments/assets/4b07a938-4236-4d26-9ee9-e639757b5f50" />


<img width="1015" height="609" alt="image" src="https://github.com/user-attachments/assets/26bff905-23a0-49a8-a9ba-ad8ad84c2a1a" />


### Variable Features

<img width="1023" height="613" alt="image" src="https://github.com/user-attachments/assets/34f8e8be-34c3-40e4-b624-4cc10a5e3017" />


---

## Reproducibility Notes

- R version: >= 4.0
- Seurat version: >= 5.0 (uses `layer` syntax for count data)
- Random seed: Set by default in most Seurat functions
- Hardware: Minimum 8 GB RAM recommended; 16 GB+ for large datasets

---

## Citation

If you use this pipeline in your research, please cite:

1. Stuart T, Butler A, Hoffman P, et al. **Seurat: an open-source R package for fast and easy analysis of single-cell RNA-seq data.** *Genome Biology*. 2021.
2. Korsunsky I, Millard N, Fan J, et al. **Fast, sensitive and accurate integration of single-cell data with Harmony.** *Nature Methods*. 2019.
3. McGinnis CS, Murrow LM, Gartner ZJ. **DoubletFinder: Doublet detection in single-cell RNA sequencing data using artificial nearest neighbors.** *Cell Systems*. 2019.
4. **GSE178481** dataset — refer to the original GEO publication.

---

## License

MIT License — feel free to use and modify for your own projects.
