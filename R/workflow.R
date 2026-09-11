pipeline_fingerprint <- function(cfg) {
  paths <- c(list.files("R", pattern = "\\.R$", full.names = TRUE), "scripts/bootstrap.R", "run_pipeline.R", cfg$annotation$panel)
  if (cfg$input$mode == "samples") paths <- c(paths, cfg$input$samples)
  for (path in sample_table(cfg)$path) paths <- c(paths, if (dir.exists(path)) list.files(path, full.names = TRUE) else path)
  if (!is.null(cfg$enrichment$gmt)) paths <- c(paths, cfg$enrichment$gmt)
  if (!is.null(cfg$annotation$manual_labels)) paths <- c(paths, cfg$annotation$manual_labels)
  installed <- utils::installed.packages()
  versions <- stats::setNames(installed[, "Version"], installed[, "Package"])
  versions <- versions[order(names(versions))]
  digest::digest(list(config = cfg, files = file_hashes(paths[!dir.exists(paths)]), versions = versions,
    R = as.character(getRversion())), algo = "sha256")
}
stage_is_current <- function(path, fingerprint) {
  if (!file.exists(path)) return(FALSE)
  state <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(state) || !identical(state$fingerprint, fingerprint) || !length(state$outputs)) return(FALSE)
  if (!all(file.exists(names(state$outputs)))) return(FALSE)
  identical(file_hashes(names(state$outputs)), state$outputs)
}
run_pipeline <- function(config_path = "config/config.R", resume = TRUE, until = NULL) {
  check_packages(); cfg <- load_config(config_path); root <- ensure_output(cfg)
  stages <- c(ingest = "stage_ingest", qc = "stage_qc", normalize = "stage_normalize", reduce = "stage_reduce",
    markers = "stage_markers", annotate = "stage_annotate", enrich = "stage_enrich", report = "stage_report")
  if (!is.null(until)) { assert(until %in% names(stages), "Unknown stage."); stages <- stages[seq_len(match(until, names(stages))) ] }
  fingerprint <- pipeline_fingerprint(cfg); upstream_changed <- FALSE; status <- list()
  write_json(cfg, file.path(root, "provenance/config_used.json"))
  for (stage in names(stages)) {
    state_path <- file.path(root, "cache", paste0(stage, ".rds"))
    current <- resume && !upstream_changed && stage_is_current(state_path, fingerprint)
    started <- Sys.time()
    if (current) message("[cached] ", stage) else {
      message("[run] ", stage); upstream_changed <- TRUE; set.seed(cfg$seed)
      outputs <- get(stages[[stage]], mode = "function")(cfg)
      save_rds(list(fingerprint = fingerprint, outputs = file_hashes(outputs)), state_path)
    }
    status[[stage]] <- data.frame(stage = stage, status = if (current) "cached" else "executed",
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
  }
  status <- do.call(rbind, status); rownames(status) <- NULL
  write_csv(status, file.path(root, "logs/last_run.csv"))
  write_json(list(implementation = "R / Seurat", completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    fingerprint = fingerprint, source_config = config_path, stages = status), file.path(root, "provenance/run_manifest.json"))
  if ("report" %in% status$stage) message("Report: ", file.path(root, "report.html"))
  invisible(list(config = cfg, stages = status, fingerprint = fingerprint))
}
