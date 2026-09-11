# Multiple captures and donor-level differential expression

Set `input$mode = "samples"` and edit `config/samples.csv`. Each row describes
one physical capture/library. `sample_id` identifies the capture, `donor` the
independent biological donor, `condition` the experimental group, and `batch`
the technical batch. Repeated captures from the same donor are not independent
biological replicates. Do not merge distinct captures into one row before
doublet detection.

```csv
sample_id,path,format,donor,condition,batch
capture1,data/capture1,10x_mtx,donor1,control,batchA
capture2,data/capture2,10x_mtx,donor2,control,batchB
capture3,data/capture3,10x_mtx,donor3,control,batchA
capture4,data/capture4,10x_mtx,donor4,case,batchB
capture5,data/capture5,10x_mtx,donor5,case,batchA
capture6,data/capture6,10x_mtx,donor6,case,batchB
```

Supported formats are `10x_mtx`, `10x_h5`, and `counts_rds`. Use gene symbols
that match the marker panel and mitochondrial prefix. For capture-specific QC,
add an override such as `qc$overrides = list(capture1 = list(max_pct_mt = 10))`.
The demo's expected doublet rate 0.04 is explicit; choose a suitable rate for
your capture loading and chemistry.

Harmony is optional. Set `embedding$integration = "harmony"` and
`embedding$batch_key = "batch"`. It adjusts the PCA embedding used for neighbors,
clustering and UMAP. Raw counts and normalized RNA data remain available for
markers and pseudobulk. Do not integrate away a real condition effect. Perfect
batch/condition confounding cannot be resolved by integration; review the
experimental design and batch mixing by cell type.

After reviewing labels, run donor-level pseudobulk:

```r
source("scripts/bootstrap.R")
pseudobulk_de(
  object_path = "results/my_study/objects/06_annotated.rds",
  output_dir = "results/my_study/pseudobulk_T",
  cell_type = "T cells", case = "case", control = "control",
  paired = FALSE, min_cells = 20L
)
```

Counts are summed within donor, condition and cell type, pooling technical
captures. Groups with fewer than 20 cells are omitted. The implementation
requires at least three independent donors per condition, a full-rank model,
and complete donor pairs when `paired = TRUE`. Unpaired models use `~condition`;
paired models use `~donor + condition`. More complex covariates require adapting
and validating the design, not assuming this basic contrast is sufficient.

Genes need at least 10 summed counts in at least three donor profiles. DESeq2
uses positive-count size factors and reports the specified case/control
contrast. Review effect sizes, uncertainty, counts and adjusted p-values.
PBMC3k has one donor: it cannot run this disease comparison. Software tests use
explicitly synthetic donors and a known synthetic effect.
