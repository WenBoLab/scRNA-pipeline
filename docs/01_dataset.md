# PBMC3k and input provenance

The demonstration uses the 10x Genomics PBMC3k cell-called matrix: 2,700 cells
from peripheral blood mononuclear cells of one healthy donor. It is small enough
for teaching QC, clustering and lineage annotation. It has no independent
disease/control comparison.

The repository contains `resources/pbmc3k/matrix.mtx.gz`, `features.tsv.gz` and
`barcodes.tsv.gz`. `Seurat::Read10X()` reads this directory. There are 32,738
measured genes, 6,390,631 UMI counts and 2,286,884 nonzero entries before filtering.
The standard matrix uses genes as rows and cells as columns.

The files were losslessly repackaged from the previously verified raw count
object, with original cell barcodes restored. R's Matrix and hdf5r packages were
used for the conversion. All matrix entries and identifiers were checked after
reading the new files back into R. No QC, normalization or imputation was applied
to the bundled input. `source.json` records hashes and the conversion check.

The workflow checks the resource SHA-256 values before running. It also checks
that raw counts remain unchanged through normalization and in the final object.
The final object contains only retained cells and genes; the complete input is
saved separately as `01_counts.rds`.

For your data, use a cell-called matrix, not an unfiltered all-droplet matrix.
Inputs can be a 10x matrix directory, a 10x H5 file, or an RDS file holding a named
sparse integer matrix. Gene symbols must match the configured marker panel.
For multiple captures, measured genes are intersected and cell names receive
a capture prefix. Gene identifier changes are recorded in `gene_mapping.csv`.

Source and attribution:

- [10x Genomics PBMC3k dataset](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0)
- [Seurat PBMC3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
- [Original raw-count mirror](https://falexwolf.de/data/pbmc3k_raw.h5ad)
- [CC BY 4.0 terms](https://creativecommons.org/licenses/by/4.0/)

Reactome_2022 gene sets are a fixed Enrichr snapshot. The original download URL,
retrieval date and SHA-256 are in `resources/reactome_source.json`. Reactome data
are [CC0](https://reactome.org/license). Enrichment uses no live database request.
