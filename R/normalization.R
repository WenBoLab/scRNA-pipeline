stage_normalize <- function(cfg) {
  object <- readRDS(object_path(cfg, "qc")); before <- digest::digest(raw_counts(object), algo = "sha256")
  object <- Seurat::NormalizeData(object, scale.factor = cfg$normalization$scale_factor, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = cfg$normalization$n_hvg, verbose = FALSE)
  if (length(unique(object$sample_id)) > 1) {
    parts <- Seurat::SplitObject(object, split.by = "sample_id")
    parts <- lapply(parts, function(x) Seurat::FindVariableFeatures(x, nfeatures = cfg$normalization$n_hvg, verbose = FALSE))
    SeuratObject::VariableFeatures(object) <- Seurat::SelectIntegrationFeatures(parts, nfeatures = cfg$normalization$n_hvg)
  }
  object <- Seurat::ScaleData(object, features = SeuratObject::VariableFeatures(object), verbose = FALSE)
  assert(identical(before, digest::digest(raw_counts(object), algo = "sha256")), "Normalization modified raw counts.")
  c(save_rds(object, object_path(cfg, "normalize")),
    write_csv(data.frame(gene = SeuratObject::VariableFeatures(object)), file.path(cfg$output_dir, "tables/variable_genes.csv")),
    save_plot(Seurat::VariableFeaturePlot(object), file.path(cfg$output_dir, "figures/variable_genes.png")))
}
stage_reduce <- function(cfg) {
  object <- readRDS(object_path(cfg, "normalize")); set.seed(cfg$seed)
  npcs <- min(cfg$embedding$n_pcs, ncol(object) - 1L, length(SeuratObject::VariableFeatures(object)) - 1L)
  assert(npcs >= 2, "Too few dimensions for PCA.")
  object <- Seurat::RunPCA(object, npcs = npcs, seed.use = cfg$seed, verbose = FALSE)
  dims <- intersect(cfg$embedding$dims, seq_len(npcs)); reduction <- "pca"
  if (cfg$embedding$integration == "harmony") {
    assert(cfg$embedding$batch_key %in% colnames(object[[]]), "Missing Harmony batch metadata.")
    assert(length(unique(object[[]][[cfg$embedding$batch_key]])) >= 2, "Harmony needs at least two batches.")
    object <- harmony::RunHarmony(object, group.by.vars = cfg$embedding$batch_key, reduction.use = "pca",
      dims.use = dims, project.dim = FALSE, verbose = FALSE)
    reduction <- "harmony"; dims <- seq_along(dims)
  }
  object <- Seurat::FindNeighbors(object, reduction = reduction, dims = dims,
    k.param = min(cfg$embedding$k, ncol(object) - 1L), verbose = FALSE)
  object <- Seurat::FindClusters(object, resolution = cfg$embedding$resolution, algorithm = 1,
    random.seed = cfg$seed, verbose = FALSE)
  object$cluster <- as.character(SeuratObject::Idents(object))
  object <- Seurat::RunUMAP(object, reduction = reduction, dims = dims, umap.method = "uwot",
    metric = "cosine", n.neighbors = min(30L, ncol(object) - 1L), seed.use = cfg$seed, verbose = FALSE)
  object@misc$embedding <- list(integration = cfg$embedding$integration, reduction = reduction, dims = dims)
  counts <- as.data.frame(table(object$cluster)); names(counts) <- c("cluster", "cells")
  coordinates <- data.frame(cell = colnames(object), SeuratObject::Embeddings(object, "umap"), cluster = object$cluster)
  c(save_rds(object, object_path(cfg, "reduce")), write_csv(counts, file.path(cfg$output_dir, "tables/cluster_counts.csv")),
    write_csv(coordinates, file.path(cfg$output_dir, "tables/umap_coordinates.csv")),
    save_plot(Seurat::ElbowPlot(object, ndims = npcs), file.path(cfg$output_dir, "figures/pca_elbow.png")),
    save_plot(Seurat::DimPlot(object, group.by = "cluster", label = TRUE) + ggplot2::ggtitle("Seurat Louvain clusters"),
      file.path(cfg$output_dir, "figures/umap_clusters.png")))
}
