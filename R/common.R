`%||%` <- function(x, y) if (is.null(x)) y else x
assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}
required_packages <- function() {
  c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "scDblFinder",
    "SingleCellExperiment", "SummarizedExperiment", "BiocParallel", "harmony",
    "DESeq2", "digest", "jsonlite", "base64enc", "hdf5r", "xgboost", "uwot")
}
check_packages <- function(packages = required_packages()) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  assert(!length(missing), paste("Missing R packages:", paste(missing, collapse = ", "),
    "Run source('scripts/install.R') from the project root."))
  invisible(TRUE)
}
load_config <- function(path = "config/config.R") {
  env <- new.env(parent = baseenv()); sys.source(path, envir = env)
  cfg <- env$config
  assert(is.list(cfg), "Config must define a list named config.")
  assert(cfg$input$mode %in% c("demo", "samples"), "input$mode must be demo or samples.")
  assert(cfg$embedding$integration %in% c("none", "harmony"), "Unknown integration method.")
  assert(cfg$qc$min_genes >= 1 && cfg$qc$max_genes > cfg$qc$min_genes, "Invalid gene QC limits.")
  assert(cfg$qc$max_pct_mt > 0 && cfg$qc$max_pct_mt <= 100, "Invalid mitochondrial limit.")
  assert(cfg$doublets$expected_rate > 0 && cfg$doublets$expected_rate < 1, "Invalid doublet rate.")
  cfg
}
ensure_output <- function(cfg) {
  for (d in c("objects", "tables", "figures", "provenance", "logs", "cache")) {
    dir.create(file.path(cfg$output_dir, d), recursive = TRUE, showWarnings = FALSE)
  }
  cfg$output_dir
}
object_path <- function(cfg, stage) {
  files <- c(ingest = "01_counts.rds", qc = "02_qc.rds", normalize = "03_normalized.rds",
             reduce = "04_clustered.rds", markers = "05_markers.rds", annotate = "06_annotated.rds")
  assert(stage %in% names(files), paste("Unknown object stage:", stage))
  file.path(cfg$output_dir, "objects", unname(files[[stage]]))
}
save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path, compress = "gzip")
  invisible(path)
}
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}
write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = 12)
  invisible(path)
}
file_sha256 <- function(path) digest::digest(file = path, algo = "sha256")
file_hashes <- function(paths) {
  paths <- sort(unique(paths))
  assert(all(file.exists(paths)), "An input or output file is missing.")
  stats::setNames(vapply(paths, file_sha256, character(1)), paths)
}
validate_counts <- function(x) {
  assert(inherits(x, "Matrix") || is.matrix(x), "Input must be a count matrix.")
  x <- methods::as(x, "dgCMatrix")
  assert(nrow(x) > 0 && ncol(x) > 0, "The matrix is empty.")
  assert(all(is.finite(x@x)) && all(x@x >= 0), "Counts must be finite and nonnegative.")
  assert(all(abs(x@x - round(x@x)) < 1e-8), "Use raw integer UMI counts.")
  assert(!is.null(rownames(x)) && !is.null(colnames(x)), "Gene and cell names are required.")
  assert(!anyNA(rownames(x)) && !anyNA(colnames(x)) &&
    all(nzchar(rownames(x))) && all(nzchar(colnames(x))), "Names cannot be empty.")
  assert(!anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x)), "Names must be unique.")
  Matrix::drop0(x)
}
raw_counts <- function(object) SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
normalized_data <- function(object) SeuratObject::LayerData(object, assay = "RNA", layer = "data")
load_panels <- function(path) {
  env <- new.env(parent = baseenv()); sys.source(path, envir = env)
  assert(is.list(env$marker_panels) && length(env$marker_panels) >= 2, "Invalid marker panels.")
  lapply(env$marker_panels, unique)
}
save_plot <- function(plot, path, width = 9, height = 6) {
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in",
    dpi = 160, bg = "white", limitsize = FALSE)
  invisible(path)
}
sample_table <- function(cfg) {
  samples <- if (cfg$input$mode == "demo") data.frame(sample_id = "pbmc3k",
    path = "resources/pbmc3k", format = "10x_mtx", donor = "donor1", condition = "healthy", batch = "batch1") else
    utils::read.csv(cfg$input$samples, colClasses = "character", check.names = FALSE)
  fields <- c("sample_id", "path", "format", "donor", "condition", "batch")
  assert(all(fields %in% names(samples)) && nrow(samples) > 0, "Sample sheet requires sample_id,path,format,donor,condition,batch.")
  assert(!anyNA(samples[fields]) && all(vapply(samples[fields], function(x) all(nzchar(x)), logical(1))), "Sample metadata cannot be empty.")
  assert(!anyDuplicated(samples$sample_id), "One unique sample_id is required per physical capture.")
  assert(all(grepl("^[A-Za-z][A-Za-z0-9.-]*$", samples$sample_id)), "Use safe capture IDs: letters, digits, dots or hyphens.")
  assert(all(samples$format %in% c("10x_mtx", "10x_h5", "counts_rds")), "Unsupported input format.")
  samples
}
