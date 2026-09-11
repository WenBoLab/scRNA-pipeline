source("scripts/bootstrap.R")
cfg <- load_config(); input <- readRDS(object_path(cfg, "ingest")); final <- readRDS(object_path(cfg, "annotate"))
audit <- read.csv(file.path(cfg$output_dir, "tables/cell_qc.csv"))
expected <- raw_counts(input)[rownames(final), colnames(final), drop = FALSE]
stopifnot(identical(dim(raw_counts(input)), c(32738L, 2700L)), sum(raw_counts(input)) == 6390631,
  Matrix::nnzero(raw_counts(final) - expected) == 0,
  setequal(colnames(final), audit$cell[audit$qc_keep]),
  all(is.finite(SeuratObject::Embeddings(final, "umap"))),
  all(final$doublet_class == "singlet"))
resumed <- run_pipeline()
stopifnot(all(resumed$stages$status == "cached"))
python_initialized <- if (requireNamespace("reticulate", quietly = TRUE)) reticulate::py_available(initialize = FALSE) else FALSE
summary <- jsonlite::read_json(file.path(cfg$output_dir, "run_summary.json"))
stopifnot(!python_initialized, identical(summary$python_runtime_initialized, FALSE))
write_json(list(implementation = "R / Seurat", input_cells = ncol(input), retained_cells = ncol(final),
  retained_genes = nrow(final), clusters = length(unique(final$cluster)), raw_count_equality = TRUE,
  complete_cell_audit = TRUE, resumed_stages = resumed$stages, python_runtime_initialized = python_initialized,
  r_version = as.character(getRversion()), validated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
  file.path(cfg$output_dir, "provenance/real_data_validation.json"))
message("Real R demo validation passed, including exact raw counts and cached resume.")
