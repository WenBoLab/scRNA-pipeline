read_capture <- function(path, format) {
  x <- switch(format, "10x_mtx" = Seurat::Read10X(path), "10x_h5" = Seurat::Read10X_h5(path),
              "counts_rds" = readRDS(path))
  if (is.list(x)) {
    assert("Gene Expression" %in% names(x), "Multimodal input requires a Gene Expression matrix.")
    x <- x[["Gene Expression"]]
  }
  validate_counts(x)
}
verify_demo <- function() {
  provenance <- jsonlite::read_json("resources/pbmc3k/source.json", simplifyVector = TRUE)
  for (name in names(provenance$files)) {
    assert(identical(file_sha256(file.path("resources/pbmc3k", name)), provenance$files[[name]]),
      paste("Bundled count-matrix checksum mismatch:", name))
  }
  invisible(provenance)
}
stage_ingest <- function(cfg) {
  root <- ensure_output(cfg)
  if (cfg$input$mode == "demo") verify_demo()
  samples <- sample_table(cfg); captures <- list(); metadata <- list(); mapping <- list()
  for (i in seq_len(nrow(samples))) {
    s <- samples[i, ]; x <- read_capture(s$path, s$format)
    genes <- rownames(x); rownames(x) <- make.unique(gsub("_", "-", genes, fixed = TRUE))
    mapping[[i]] <- data.frame(sample_id = s$sample_id, input_gene = genes, analysis_gene = rownames(x))
    barcode <- colnames(x); colnames(x) <- paste0(s$sample_id, "__", barcode)
    metadata[[i]] <- data.frame(sample_id = s$sample_id, donor = s$donor, condition = s$condition,
      batch = s$batch, original_barcode = barcode, row.names = colnames(x))
    captures[[i]] <- x
  }
  shared <- Reduce(intersect, lapply(captures, rownames))
  assert(length(shared) >= 20, "Too few shared genes across captures.")
  counts <- do.call(cbind, lapply(captures, function(x) x[shared, , drop = FALSE]))
  object <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = do.call(rbind, metadata),
    project = cfg$project, min.cells = 0, min.features = 0)
  if (cfg$input$mode == "demo") {
    assert(identical(dim(counts), c(32738L, 2700L)) && sum(counts) == 6390631, "Unexpected bundled input.")
  }
  inputs <- unlist(lapply(samples$path, function(p) if (dir.exists(p)) list.files(p, full.names = TRUE) else p))
  c(save_rds(object, object_path(cfg, "ingest")),
    write_csv(do.call(rbind, mapping), file.path(root, "tables/gene_mapping.csv")),
    write_csv(samples, file.path(root, "tables/samples.csv")),
    write_json(list(implementation = "R / Seurat", cells = ncol(counts), genes = nrow(counts),
      total_umi = sum(counts), shared_gene_policy = "Intersection across captures",
      files = as.list(file_hashes(inputs))), file.path(root, "provenance/inputs.json")))
}
