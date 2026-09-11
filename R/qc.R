stage_qc <- function(cfg) {
  root <- ensure_output(cfg); object <- readRDS(object_path(cfg, "ingest")); x <- raw_counts(object)
  mt <- grepl(cfg$qc$mitochondrial_pattern, rownames(x))
  assert(any(mt), "No mitochondrial genes matched; check gene identifiers and the configured prefix.")
  audit <- object[[]]
  audit$cell <- rownames(audit)
  audit$percent_mt <- 100 * Matrix::colSums(x[mt, , drop = FALSE]) / pmax(Matrix::colSums(x), 1)
  audit$doublet_score <- NA_real_; audit$doublet_class <- "not_tested"
  audit$qc_keep <- FALSE; audit$removal_reason <- ""
  for (i in seq_along(unique(audit$sample_id))) {
    sid <- unique(audit$sample_id)[i]; idx <- which(audit$sample_id == sid)
    limits <- utils::modifyList(cfg$qc, cfg$qc$overrides[[sid]] %||% list())
    preliminary <- idx[audit$nFeature_RNA[idx] >= limits$min_genes]
    if (cfg$doublets$enabled) {
      assert(length(preliminary) >= 100, paste("scDblFinder needs at least 100 screened cells for capture", sid))
      local <- x[, preliminary, drop = FALSE]
      local <- local[Matrix::rowSums(local > 0) >= 3, , drop = FALSE]
      sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = local))
      set.seed(cfg$seed + i)
      sce <- scDblFinder::scDblFinder(sce, dbr = cfg$doublets$expected_rate,
        BPPARAM = BiocParallel::SerialParam(RNGseed = cfg$seed + i))
      cd <- as.data.frame(SummarizedExperiment::colData(sce))
      audit$doublet_score[preliminary] <- cd$scDblFinder.score
      audit$doublet_class[preliminary] <- as.character(cd$scDblFinder.class)
    } else audit$doublet_class[preliminary] <- "disabled"
    for (j in idx) {
      reasons <- c(if (audit$nFeature_RNA[j] < limits$min_genes) "low_genes",
        if (audit$nFeature_RNA[j] >= limits$max_genes) "high_genes",
        if (audit$percent_mt[j] >= limits$max_pct_mt) "high_mitochondrial",
        if (cfg$doublets$enabled && cfg$doublets$remove && audit$doublet_class[j] == "doublet") "predicted_doublet")
      audit$removal_reason[j] <- if (length(reasons)) paste(reasons, collapse = ";") else "retained"
      audit$qc_keep[j] <- !length(reasons)
    }
  }
  keep <- rownames(audit)[audit$qc_keep]
  assert(length(keep) >= 20, "Too few cells remain; review cell_qc thresholds.")
  genes <- rownames(x)[Matrix::rowSums(x[, keep, drop = FALSE] > 0) >= cfg$qc$min_cells_per_gene]
  object <- subset(object, cells = keep, features = genes)
  object <- SeuratObject::AddMetaData(object, audit[colnames(object), c("percent_mt", "doublet_score", "doublet_class")])
  metrics <- do.call(rbind, lapply(c("nFeature_RNA", "nCount_RNA", "percent_mt"), function(m)
    data.frame(sample_id = audit$sample_id, kept = ifelse(audit$qc_keep, "Retained", "Removed"), metric = m, value = audit[[m]])))
  p <- ggplot2::ggplot(metrics, ggplot2::aes(kept, value, fill = kept)) +
    ggplot2::geom_violin(scale = "width") + ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = NULL, title = "Quality control of all input cells") +
    ggplot2::theme(legend.position = "none")
  d <- audit[is.finite(audit$doublet_score), ]
  p2 <- if (nrow(d)) ggplot2::ggplot(d, ggplot2::aes(doublet_score, fill = doublet_class)) +
    ggplot2::geom_histogram(bins = 50, position = "identity", alpha = 0.6) + ggplot2::theme_bw() +
    ggplot2::labs(title = "scDblFinder predictions", x = "Doublet score", y = "Cells") else
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "Doublet detection disabled") + ggplot2::theme_void()
  summary <- do.call(rbind, lapply(split(audit, audit$sample_id), function(a) data.frame(sample_id = a$sample_id[1],
    input_cells = nrow(a), kept_cells = sum(a$qc_keep), predicted_doublets = sum(a$doublet_class == "doublet"))))
  c(save_rds(object, object_path(cfg, "qc")), write_csv(audit, file.path(root, "tables/cell_qc.csv")),
    write_csv(summary, file.path(root, "tables/qc_summary.csv")),
    save_plot(p, file.path(root, "figures/qc.png")), save_plot(p2, file.path(root, "figures/doublets.png")))
}
