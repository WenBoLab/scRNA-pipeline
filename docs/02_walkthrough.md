# Run the stages in R

Load the functions and settings from the project root:

```r
source("scripts/bootstrap.R")
config <- load_config("config/config.R")
check_packages()
```

Run these functions in order on a new output directory. `run_pipeline()` calls
the same functions and adds checksum-based resume behavior.

```r
stage_ingest(config)
stage_qc(config)
stage_normalize(config)
stage_reduce(config)
stage_markers(config)
stage_annotate(config)
stage_enrich(config)
stage_report(config)
```

1. Import validates nonnegative, finite integer counts and capture metadata.
   Inspect `provenance/inputs.json` and `tables/gene_mapping.csv`.
2. QC computes genes, UMIs and mitochondrial percentage. A minimum-gene screen
   precedes scDblFinder; final upper-gene and mitochondrial filters follow it.
   Doublets are predicted separately per capture. Every cell receives an audit
   row, including cells not tested for doublets. Thresholds are teaching settings,
   not universal biological boundaries.
3. LogNormalize scales each cell to 10,000 and applies log1p. VST selects 2,000
   variable genes. Multiple captures use per-capture variable features followed
   by SelectIntegrationFeatures. ScaleData acts on variable genes; raw counts
   are preserved.
4. PCA summarizes scaled variable-gene expression. Neighbors and Louvain define
   clusters. UMAP uses the R uwot backend. Resolution and PCA dimension choices
   can change membership. UMAP separation does not measure a biological distance.
5. Wilcoxon tests compare each cluster to all other cells. Genes expressed in
   at least 10% of either group are eligible. Fold-change filtering follows the
   tests. Seurat uses Bonferroni adjustment over assay genes. These exploratory
   markers are not independent donor-level disease tests.
6. Marker panels produce candidate labels with evidence tables. Scores are
   heuristic and ambiguous groups remain unresolved. Review the dot plot and
   feature plots before accepting labels.
7. Reactome over-representation uses the positive markers and each cluster's
   actually tested genes as background. BH correction includes every eligible
   pathway within the cluster, even when overlap is zero. Inspect overlap genes
   and redundant pathways before interpreting themes.
8. The R report embeds the plots for offline viewing and records configuration,
   input hashes, R session information and package versions.

Inspect saved objects and tables:

```r
pbmc <- readRDS(object_path(config, "annotate"))
SeuratObject::Layers(pbmc[["RNA"]])
head(pbmc[[]])
Seurat::DimPlot(pbmc, group.by = "cluster", label = TRUE)
head(read.csv(file.path(config$output_dir, "tables/markers_top10.csv")))
```

The executable lesson repeats these steps with questions and saves its own
results in `results/notebook_demo/`.
