PROJECT_ROOT <- normalizePath(getwd())
while (!file.exists(file.path(PROJECT_ROOT, "scrna-seq-learning.Rproj"))) {
  parent <- dirname(PROJECT_ROOT)
  if (identical(parent, PROJECT_ROOT)) stop("Cannot find the project root.")
  PROJECT_ROOT <- parent
}
make_synthetic_study <- function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  set.seed(741)
  panels <- load_panels("config/markers.R")
  genes <- unique(c("MT-CO1", "MT-ND1", unlist(panels), sprintf("GENE%03d", 1:170)))
  rows <- list(); known <- character(); captures <- 0L
  for (donor in 1:6) for (capture in 1:2) {
    captures <- captures + 1L; id <- paste0("capture", captures)
    condition <- if (donor <= 3) "control" else "case"
    types <- rep(c("T cells", "B cells", "CD14 monocytes"), each = 24)
    mu <- matrix(0.6, nrow = length(genes), ncol = length(types), dimnames = list(genes, paste0("cell", seq_along(types))))
    mu[c("MT-CO1", "MT-ND1"), ] <- 0.05
    for (type in unique(types)) mu[panels[[type]], types == type] <- 12
    if (condition == "case") mu[sprintf("GENE%03d", 1:6), types == "T cells"] <- 7
    mu <- mu * exp(stats::rnorm(1, 0, 0.08))
    counts <- matrix(stats::rnbinom(length(mu), mu = as.vector(mu), size = 12), nrow = nrow(mu), dimnames = dimnames(mu))
    path <- file.path(directory, paste0(id, ".rds")); saveRDS(methods::as(counts, "dgCMatrix"), path)
    rows[[captures]] <- data.frame(sample_id = id, path = path, format = "counts_rds", donor = paste0("donor", donor),
      condition = condition, batch = paste0("batch", capture))
    known <- c(known, stats::setNames(types, paste0(id, "__", colnames(counts))))
  }
  samples <- file.path(directory, "samples.csv"); write_csv(do.call(rbind, rows), samples)
  cfg <- load_config("config/config.R"); cfg$project <- "explicitly_synthetic"
  cfg$output_dir <- file.path(directory, "results"); cfg$input <- list(mode = "samples", samples = samples)
  cfg$qc$min_genes <- 20L; cfg$qc$max_genes <- 500L; cfg$qc$max_pct_mt <- 20
  cfg$doublets$enabled <- FALSE; cfg$normalization$n_hvg <- 100L
  cfg$embedding$n_pcs <- 12L; cfg$embedding$dims <- 1:8; cfg$embedding$resolution <- 0.3
  cfg$enrichment$gmt <- NULL
  write_config <- function(config, path) writeLines(c("config <-", capture.output(dput(config))), path)
  config_path <- file.path(directory, "config.R"); write_config(cfg, config_path)
  list(config = cfg, config_path = config_path, known_types = known, write_config = write_config)
}
