candidate_label <- function(best, score, margin, t_core, nk_core, cytotoxic, min_score, min_margin) {
  if (best %in% c("T cells", "NK cells")) {
    if (t_core >= 0.20) return(if (cytotoxic >= 0.25) "Cytotoxic T cells (candidate)" else "T cells (candidate)")
    if (t_core < 0.10 && nk_core >= 0.15 && score >= min_score && margin >= min_margin) return("NK cells (candidate)")
    return("T/NK unresolved")
  }
  if (score < min_score || margin < min_margin) return("Unresolved")
  paste(best, "(candidate)")
}
validate_manual_labels <- function(manual, clusters) {
  assert(all(c("cluster", "cell_type", "evidence") %in% names(manual)), "Manual labels need cluster,cell_type,evidence.")
  assert(!anyDuplicated(manual$cluster) && setequal(as.character(manual$cluster), clusters), "Supply exactly one manual label per cluster.")
  assert(!anyNA(manual) && all(nzchar(trimws(manual$cell_type))) && all(nzchar(trimws(manual$evidence))), "Every manual label needs evidence.")
  invisible(TRUE)
}
stage_annotate <- function(cfg) {
  object <- readRDS(object_path(cfg, "markers")); data <- normalized_data(object)
  panels <- load_panels(cfg$annotation$panel); clusters <- sort(unique(object$cluster))
  genes <- unique(unlist(panels)); present <- intersect(genes, rownames(data))
  means <- fractions <- matrix(0, nrow = length(genes), ncol = length(clusters), dimnames = list(genes, clusters))
  for (cluster in clusters) {
    cells <- which(object$cluster == cluster)
    means[present, cluster] <- Matrix::rowMeans(data[present, cells, drop = FALSE])
    fractions[present, cluster] <- Matrix::rowMeans(data[present, cells, drop = FALSE] > 0)
  }
  scaled <- means / pmax(apply(means, 1, max), 1e-8)
  scores <- vapply(panels, function(p) colMeans(scaled[p, , drop = FALSE]), numeric(length(clusters)))
  rownames(scores) <- clusters
  frac <- function(gs, cl) {
    available <- intersect(gs, rownames(fractions))
    if (length(available)) mean(fractions[available, cl, drop = TRUE]) else 0
  }
  evidence <- lapply(clusters, function(cl) {
    s <- sort(scores[cl, ], decreasing = TRUE); best <- names(s)[1]; score <- unname(s[1]); margin <- unname(s[1] - s[2])
    tc <- frac(c("CD3D", "CD3E", "TRAC"), cl); nk <- frac(c("GNLY", "KLRD1"), cl); cy <- frac(c("NKG7", "PRF1", "GNLY"), cl)
    data.frame(cluster = cl, panel = best, panel_score = score, score_margin = margin,
      t_core_fraction = tc, nk_core_fraction = nk, cytotoxic_fraction = cy,
      cell_type = candidate_label(best, score, margin, tc, nk, cy, cfg$annotation$min_score, cfg$annotation$min_margin),
      annotation_status = "candidate", evidence = "Heuristic panel support; inspect marker plots before accepting")
  })
  evidence <- do.call(rbind, evidence)
  if (!is.null(cfg$annotation$manual_labels)) {
    manual <- utils::read.csv(cfg$annotation$manual_labels, colClasses = "character")
    validate_manual_labels(manual, clusters); index <- match(evidence$cluster, manual$cluster)
    evidence$cell_type <- manual$cell_type[index]; evidence$evidence <- manual$evidence[index]
    evidence$annotation_status <- "reviewed"
  }
  index <- match(object$cluster, evidence$cluster)
  object$cell_type <- evidence$cell_type[index]; object$annotation_status <- evidence$annotation_status[index]
  expression <- do.call(rbind, lapply(clusters, function(cl) data.frame(cluster = cl, gene = genes,
    present = genes %in% present, mean_log_expression = means[, cl], fraction_expressing = fractions[, cl])))
  counts <- as.data.frame(table(object$cell_type)); names(counts) <- c("cell_type", "cells")
  template <- data.frame(cluster = clusters, cell_type = "", evidence = "")
  c(save_rds(object, object_path(cfg, "annotate")),
    write_csv(evidence, file.path(cfg$output_dir, "tables/annotation_evidence.csv")),
    write_csv(data.frame(cluster = clusters, scores, check.names = FALSE), file.path(cfg$output_dir, "tables/annotation_scores.csv")),
    write_csv(expression, file.path(cfg$output_dir, "tables/annotation_expression.csv")),
    write_csv(counts, file.path(cfg$output_dir, "tables/celltype_counts.csv")),
    write_csv(template, file.path(cfg$output_dir, "tables/manual_labels_template.csv")),
    save_plot(Seurat::DimPlot(object, group.by = "cell_type") + ggplot2::ggtitle("Candidate cell annotations"),
      file.path(cfg$output_dir, "figures/umap_celltypes.png"), 12, 6))
}
