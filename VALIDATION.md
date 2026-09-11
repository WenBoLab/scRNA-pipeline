# R validation record

This release contains R analysis, an R installer, an R Markdown lesson and R tests.
All README files are English-only ASCII and no Python scripts, Python notebooks,
Snakemake files or previous analysis outputs are included in the release archive.

## Real PBMC3k analysis

The complete count-matrix-to-report workflow executed successfully on Linux
with R 4.4.3 and Bioconductor 3.20.

| Measurement | Observed value |
|---|---:|
| Input cells | 2,700 |
| Input measured genes | 32,738 |
| Input UMI counts | 6,390,631 |
| Retained cells | 2,542 |
| Retained genes | 13,522 |
| Louvain clusters | 9 |
| Predicted doublets before final QC | 96 |
| Retained UMI counts | 5,867,141 |

Verified invariants: the final raw counts exactly equal the corresponding subset
of the imported raw matrix; the audit accounts for every input cell; the retained
cell IDs match the audit; UMAP values are finite; a second run reuses all eight
stages after checking file hashes and the environment fingerprint. The R-native
analysis did not initialize a Python runtime. Input repackaging was lossless,
including all count entries and cell/gene identifiers.

Reactome over-representation ran on the fixed bundled collection using the
tested-gene universe for each cluster. All annotation labels remain candidates.
Any unresolved groups are listed explicitly in the annotation evidence table.

## Tests and executable lesson

All 13 R tests passed with zero failed expectations and zero errored tests.
Checks include invalid-count rejection, exact bundled counts, raw-count retention,
ORA background and multiple testing, T/NK ambiguity, required annotation evidence,
technical-capture aggregation, independent-donor design, cache invalidation,
ASCII README validation, and an end-to-end synthetic study.
Harmony and DESeq2 were executed on explicitly synthetic multiple-donor data.
The DESeq2 test recovered the direction of a known synthetic effect. These
synthetic tests do not establish biological validity in a real disease study.

All 9 code chunks of the R Markdown lesson executed successfully.
The rendered lesson embeds 6 PNG images; no image path is unresolved.
Its retained counts, dimensions and cluster count match the one-command run.
The standalone analysis report embeds all 9 figures.

## Package specification

The lockfile records 232 R packages with explicit CRAN or Bioconductor
sources. `renv::restore()` confirmed that the tested installed library matches
the lockfile. No local-only package source is required by the lockfile.

| Package | Tested version |
|---|---|
| Seurat | 5.2.1 |
| SeuratObject | 5.0.2 |
| Matrix | 1.7-6 |
| scDblFinder | 1.20.2 |
| xgboost | 1.7.11.1 |
| harmony | 1.2.4 |
| DESeq2 | 1.46.0 |
| ggplot2 | 3.5.2 |
| uwot | 0.2.5 |
| renv | 1.2.4 |

The validation environment used conda-forge/bioconda packages and a CRAN source
build of xgboost. A completely fresh CRAN/Bioconductor source restoration of every
package, interactive RStudio, Windows and macOS were not tested. The prebuilt
SeuratObject binary reported an older Matrix build version; raw-count, normalization,
clustering and integration checks nevertheless passed. Fresh source restoration
builds SeuratObject against the restored Matrix version.

Seurat printed an informational UMAP-backend notice and the synthetic three-group
dot plot warned about scaling with few groups. These were not analysis failures;
UMAP explicitly used R uwot. Full test counts and real-data validation records
are included in `examples/pbmc3k/validation` and `examples/pbmc3k/provenance`.

## Scope and publication status

Actual FASTQ alignment, empty-droplet calling, ambient-RNA correction, trajectory
inference and RNA velocity were not executed. PBMC3k has one donor and cannot
support a disease/control donor comparison. The optional R Cell Ranger launcher
requires an external installation and suitable FASTQs/reference data.

GitHub Actions configurations are supplied but remote CI was not executed.
Publication to flee70973-coder/F remains blocked by GitHub integration error
403: Resource not accessible by integration. The replacement upload archive is
complete and must be extracted into the repository root.
