# Review cell annotations

Automatic labels are candidates, not accepted ground truth. Compare
`annotation_evidence.csv`, `annotation_expression.csv`, `markers_top10.csv`, the
marker dot plot and the feature plots. Evaluate multiple markers and incompatible
lineages. Do not accept a label solely because it has the highest panel score.

| Lineage | Examples of informative markers |
|---|---|
| T cells | CD3D, CD3E, TRAC; IL7R supports some T-cell states |
| B cells | MS4A1, CD79A, CD79B |
| NK cells | GNLY, KLRD1 together with low T-cell receptor markers |
| CD14 monocytes | LYZ, S100A8, S100A9, CD14, FCN1 |
| FCGR3A monocytes | FCGR3A, MS4A7, LST1 |
| Dendritic cells | FCER1A, CD1C, CLEC10A |
| Platelets | PPBP, PF4 |

NKG7 and PRF1 also occur in cytotoxic T cells. The candidate rule checks T-cell
markers alongside NK-specific evidence. Mixed T/NK evidence remains unresolved.
Panel means are scaled across clusters and averaged; missing genes contribute
zero. The score is a relative heuristic, not a probability or reference mapping.

To accept or change labels:

1. Copy `results/pbmc3k/tables/manual_labels_template.csv` to
   `config/manual_labels.csv`.
2. Fill `cell_type` and `evidence` for every cluster; retain the cluster IDs.
3. Set `annotation$manual_labels = "config/manual_labels.csv"` in the config.
4. Run `source("run_pipeline.R")` again.

The pipeline requires exactly one nonempty label and evidence statement per
cluster. Reviewed labels receive `annotation_status = "reviewed"`. Donor-level
DE requires reviewed labels. If upstream parameters change, review the new
clusters before reusing a manual label sheet: the same number need not identify
the same cell population.
