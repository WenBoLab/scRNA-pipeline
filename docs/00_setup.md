# Setup and troubleshooting

Use R 4.4.3 and open `scrna-seq-learning.Rproj` in RStudio. Run:

```r
source("scripts/install.R")
source("run_pipeline.R")
```

The analysis uses R packages, including Seurat 5.2.1, SeuratObject 5.0.2,
scDblFinder 1.20.2 and Bioconductor 3.20. Restore `renv.lock` instead of installing
the newest packages independently. Seurat 5.2.1 is incompatible with the removed
`slot` interface in newer SeuratObject releases; the lockfile keeps the tested pair.

Package restoration needs internet access and can take time. Windows and macOS
may use binary packages; older pinned versions may require source compilation.
Use Rtools44 on Windows or the appropriate Xcode command-line tools and Fortran
compiler on macOS when source compilation is required. Linux needs compilers and
development libraries for curl, SSL, XML, HDF5, fonts and graphics.

`environment.yml` is an optional conda-forge/bioconda bootstrap for Linux.
Run `micromamba create -f environment.yml`, activate `scrna-learning-r`, then
run `Rscript scripts/install.R` to restore the exact R lockfile before analysis.
The bootstrap includes compilers and many prebuilt packages; the lockfile remains
the authoritative R package specification. RStudio users can use renv directly.

At least 8 GB RAM is a practical starting point for this small dataset. Larger
studies require more resources. The pipeline uses one analysis worker.

Run all commands at the project root, not inside `R/` or `scripts/`.

| Problem | Action |
|---|---|
| Missing packages | Run the installation script from the project root. |
| Missing input/checksum mismatch | Restore the original bundled resources. |
| No mitochondrial genes | Check gene identifiers and the species-specific prefix. |
| Too few retained cells | Review `config/config.R` and the QC distributions. |
| No relevant marker genes | Check species and gene symbols; edit `config/markers.R`. |
| Notebook rendering fails | Install R Markdown dependencies and make Pandoc available. RStudio includes Pandoc. |
| Pseudobulk is refused | Check reviewed labels, independent donors and complete conditions. |

The README is English-only ASCII Markdown. Hash signs, brackets and backticks
in a text editor are Markdown formatting; GitHub renders them on the repository
home page. Save edited text as UTF-8. Supplied README files are checked for
non-ASCII bytes before packaging.
