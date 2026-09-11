# Optional FASTQ entry point

The default demo starts from a cell-called count matrix. Actual FASTQ alignment
has not been run or validated in this project. `scripts/cellranger_count.R` is
an R launcher for an external 10x Genomics Cell Ranger installation.

Install Cell Ranger separately, obtain the appropriate reference and FASTQs,
then edit `config/cellranger.csv`. From the project root run:

```r
source("scripts/cellranger_count.R")
```

The launcher invokes `cellranger count` with sample, FASTQ directory, reference,
local CPU/RAM settings and `--create-bam=false`. It refuses to overwrite an
existing sample output directory and fails on a nonzero external exit status.
Cell Ranger needs substantially more resources than the bundled teaching demo.
Check the documentation for your installed Cell Ranger version.

After reviewing the upstream web summary, point `config/samples.csv` at each
`filtered_feature_bc_matrix` directory and set the input mode to `samples`.
Preserve physical capture and donor metadata. The count-matrix workflow does not
perform empty-droplet calling, full ambient-RNA correction or read alignment.

[Official Cell Ranger count documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-gex-count)
