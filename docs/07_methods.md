# Methods, defaults and primary references

| Setting | Default | Interpretation |
|---|---|---|
| Minimum detected genes | 200 | Applied before doublet detection |
| Maximum detected genes | less than 2,500 | Final QC |
| Mitochondrial percentage | less than 5% | Final QC, human MT- prefix |
| Gene retention | expressed in at least 3 retained cells | Preserves raw integer counts |
| Doublets | scDblFinder, expected rate 0.04 | Explicit demo setting; per capture |
| Normalize | LogNormalize, scale factor 10,000 | RNA normalized data layer |
| Variable genes | 2,000, VST | Per-capture consensus for multiple captures |
| PCA / graph | 30 PCs, dimensions 1:20, k = 20 | Limited by available matrix dimensions |
| Clustering | Louvain, resolution 0.5 | Exploratory partition |
| UMAP | R uwot, cosine metric | Visualization, seed 42 |
| Markers | Wilcoxon, prevalence at least 0.1 | Bonferroni-adjusted p < 0.05 and log2FC >= 0.25 after testing |
| Annotation | Panel score >= 0.25, margin >= 0.10 | Heuristic candidate labels with T/NK checks |
| Enrichment | Top 200 positive markers; set sizes 5-500 | Tested-gene background, BH within cluster |

All thresholds require reconsideration for other tissues, chemistries or study
designs. One small healthy-donor dataset cannot validate every biological use.
Doublets are uncertain predictions. Automatic labels need review. Cell-level
marker statistics are exploratory because clusters were learned from these
same expression values. Pathway over-representation does not establish pathway
activation or causality. No disease comparison, trajectory, RNA velocity,
ambient-RNA correction or FASTQ alignment is claimed for the executed demo.

Primary documentation:

- [Seurat PBMC3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
- [Seurat FindMarkers](https://satijalab.org/seurat/reference/findmarkers)
- [Seurat differential expression and pseudobulk](https://satijalab.org/seurat/articles/de_vignette.html)
- [scDblFinder vignette](https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html)
- [Harmony project](https://github.com/immunogenomics/harmony)
- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [renv introduction](https://rstudio.github.io/renv/articles/renv.html)
- [Reactome licensing](https://reactome.org/license)

Online documentation can describe newer versions; the supplied lockfile records
the exact versions used for this project's validation.
