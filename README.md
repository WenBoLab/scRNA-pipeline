# scRNA-seq Learning Pipeline in R

A complete count-matrix-to-report teaching workflow built with **R and Seurat v5**.
Analysis code, installation scripts, tests, and the executable lesson all use R.
All README files and project guides are English-only ASCII text.

Repository: [flee70973-coder/F](https://github.com/flee70973-coder/F)

## Start in RStudio

1. Install R 4.4.3 and RStudio. Download and extract this project.
2. Open `scrna-seq-learning.Rproj`.
3. Run these two lines in the RStudio Console:

```r
source("scripts/install.R")
source("run_pipeline.R")
```

The first command restores the versions in `renv.lock` and needs internet access.
The second runs the bundled PBMC3k demo. After packages are installed, the default
analysis needs no data download: both the count matrix and Reactome gene sets are
included.

Open `results/pbmc3k/report.html` when the run finishes. Inspect the final object:

```r
pbmc <- readRDS("results/pbmc3k/objects/06_annotated.rds")
Seurat::DimPlot(pbmc, group.by = "cell_type")
```

Terminal equivalent: `Rscript scripts/install.R`, then `Rscript run_pipeline.R`.
See [setup and troubleshooting](docs/00_setup.md) for system requirements.

## Analysis stages

| Stage | R method | Output |
|---|---|---|
| Import | Seurat Read10X / Read10X_h5 / readRDS | Raw counts and input provenance |
| Quality control | Seurat metrics and scDblFinder per capture | Every-cell audit and doublet calls |
| Normalize | LogNormalize, VST, ScaleData | Normalized values and variable genes |
| Reduce and cluster | PCA, shared-neighbor graph, Louvain, uwot | Clusters and R-native UMAP |
| Find markers | Seurat Wilcoxon tests | All tested and positive marker tables |
| Annotate | Marker panels and optional manual review | Candidate labels and evidence |
| Enrich | Hypergeometric tests with BH correction | Reactome over-representation |
| Report | R-generated HTML with embedded plots | Standalone offline report |

Optional R extensions provide Harmony integration and donor-level pseudobulk
differential expression with DESeq2. The R FASTQ launcher calls an external Cell
Ranger installation. Actual FASTQ alignment was not executed for this demo.

## Learn step by step

Open [the executable R Markdown lesson](notebooks/01_pbmc3k_walkthrough.Rmd).
It runs the same analysis functions in eight sections and includes questions
for interpreting the results. Render it from the project root:

```r
source("scripts/render_notebook.R")
```

| Guide | Topic |
|---|---|
| [00 - Setup](docs/00_setup.md) | RStudio, package restoration, troubleshooting |
| [01 - Dataset](docs/01_dataset.md) | PBMC3k, input format, checksums, attribution |
| [02 - Walkthrough](docs/02_walkthrough.md) | Stage functions and their interpretation |
| [03 - Annotation](docs/03_annotation.md) | Marker evidence and manual labels |
| [04 - Multiple samples](docs/04_multisample.md) | Capture metadata, Harmony, pseudobulk |
| [05 - FASTQ](docs/05_fastq.md) | Optional upstream Cell Ranger entry point |
| [06 - GitHub](docs/06_github.md) | Which files to upload and how to inspect CI |
| [07 - Methods](docs/07_methods.md) | Parameters, limitations, primary references |

## Dataset and example results

The bundled 10x Genomics PBMC3k matrix contains 2,700 cells, 32,738 measured genes
and 6,390,631 UMI counts from one healthy donor. It is stored in a standard 10x
sparse matrix directory that R reads directly. Source records and SHA-256 hashes
are in `resources/pbmc3k/source.json`.

The verified R run retained **2,542 cells**, **13,522 genes** and **9 clusters**.
**13 R tests and all 9 lesson code chunks passed.**

The [validation record](VALIDATION.md) reports the observed results and exact
tested R package versions. `examples/pbmc3k/` contains real R-generated plots,
tables, provenance and an HTML report. Results depend on versions and parameters.

![Candidate cell labels generated with R](examples/pbmc3k/figures/umap_celltypes.png)

Automatic annotations remain candidates until reviewed. Inspect several lineage
markers together: NKG7 and PRF1 alone cannot distinguish NK from cytotoxic T cells.
PBMC3k has no independent disease/control donor groups, so it cannot demonstrate
donor-level disease differential expression. Synthetic extension tests are
explicitly identified as synthetic.

## Run selected stages or resume

```r
source("scripts/bootstrap.R")
run_pipeline("config/config.R", until = "normalize")
run_pipeline("config/config.R")
run_pipeline("config/config.R", resume = FALSE)
source("tests/run_tests.R")
```

Stages are reused only when input hashes, configuration, R source, package
versions and saved output hashes match. Changes conservatively rerun downstream
analysis. Run one process per output directory; use different output paths for
concurrent experiments. The lesson uses its own output directory.

## Project files

| Path | Purpose |
|---|---|
| `run_pipeline.R` | One-command R entry point |
| `scrna-seq-learning.Rproj` | RStudio project with UTF-8 encoding |
| `R/` | Analysis, workflow and reporting functions |
| `config/` | R settings, marker panels and sample sheets |
| `scripts/` | Installation, rendering and optional extensions |
| `notebooks/` | R Markdown lesson |
| `resources/` | Verified counts and Reactome snapshot |
| `examples/pbmc3k/` | Compact real R results |
| `tests/` | R validation suite |
| `renv.lock` | Exact R package versions |
| `.github/workflows/` | R tests and manual real-data demonstration |

Code: MIT. PBMC3k data and derived figures: CC BY 4.0 with 10x Genomics attribution.
Reactome data: CC0. Resource source records contain the original links.
