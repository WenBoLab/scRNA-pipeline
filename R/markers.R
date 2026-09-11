stage_markers <- function(cfg) {
  object <- readRDS(object_path(cfg, "reduce")); SeuratObject::Idents(object) <- object$cluster
  data <- normalized_data(object); tables <- list(); universes <- list(); audit <- list()
  for (cluster in sort(unique(object$cluster))) {
    inside <- which(object$cluster == cluster); outside <- which(object$cluster != cluster)
    if (length(inside) < 3 || length(outside) < 3) {
      audit[[cluster]] <- data.frame(cluster = cluster, status = "skipped_too_few_cells", tested_genes = 0)
      universes[[cluster]] <- character(); next
    }
    prevalence <- pmax(Matrix::rowMeans(data[, inside, drop = FALSE] > 0), Matrix::rowMeans(data[, outside, drop = FALSE] > 0))
    eligible <- rownames(data)[prevalence >= cfg$markers$min_pct]
    result <- Seurat::FindMarkers(object, ident.1 = cluster, assay = "RNA", test.use = "wilcox",
      features = eligible, min.pct = 0, logfc.threshold = 0, min.diff.pct = -Inf, only.pos = FALSE, verbose = FALSE)
    result$gene <- rownames(result); result$cluster <- cluster; rownames(result) <- NULL
    universes[[cluster]] <- result$gene[is.finite(result$p_val)]
    tables[[cluster]] <- result
    audit[[cluster]] <- data.frame(cluster = cluster, status = "tested", tested_genes = length(universes[[cluster]]))
  }
  empty <- data.frame(p_val = numeric(), avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric(), gene = character(), cluster = character())
  all <- if (length(tables)) do.call(rbind, tables) else empty
  positive <- all[is.finite(all$p_val_adj) & all$p_val_adj < cfg$markers$max_adjusted_p & all$avg_log2FC >= cfg$markers$min_log2fc, ]
  positive <- positive[order(positive$cluster, -positive$avg_log2FC), ]; rownames(positive) <- NULL
  top <- if (nrow(positive)) do.call(rbind, lapply(split(positive, positive$cluster), head, 10)) else empty
  features <- intersect(unique(unlist(load_panels(cfg$annotation$panel))), rownames(object))
  if (!length(features)) features <- head(SeuratObject::VariableFeatures(object), 12)
  dot <- Seurat::DotPlot(object, features = features, group.by = "cluster") + Seurat::RotatedAxis() +
    ggplot2::labs(title = "Lineage marker evidence", x = NULL, y = "Cluster")
  feature <- Seurat::FeaturePlot(object, features = head(features, 8), ncol = 4)
  c(save_rds(object, object_path(cfg, "markers")),
    save_rds(universes, file.path(cfg$output_dir, "objects/marker_universes.rds")),
    write_csv(all, file.path(cfg$output_dir, "tables/markers_all.csv")),
    write_csv(positive, file.path(cfg$output_dir, "tables/markers_positive.csv")),
    write_csv(top, file.path(cfg$output_dir, "tables/markers_top10.csv")),
    write_csv(do.call(rbind, audit), file.path(cfg$output_dir, "tables/marker_test_audit.csv")),
    save_plot(dot, file.path(cfg$output_dir, "figures/marker_dotplot.png"), 12, 5),
    save_plot(feature, file.path(cfg$output_dir, "figures/feature_markers.png"), 14, 7))
}
