read_gmt <- function(path) {
  lines <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  lines <- lines[lengths(lines) >= 3]
  stats::setNames(lapply(lines, function(x) unique(x[-c(1, 2)][nzchar(x[-c(1, 2)])])), vapply(lines, `[[`, character(1), 1))
}
ora <- function(query, universe, sets, min_size = 5L, max_size = 500L) {
  universe <- unique(universe); query <- intersect(unique(query), universe)
  sets <- lapply(sets, intersect, y = universe)
  sets <- sets[lengths(sets) >= min_size & lengths(sets) <= max_size]
  empty <- data.frame(term = character(), overlap = integer(), set_size = integer(), query_size = integer(),
    universe_size = integer(), p_value = numeric(), adjusted_p = numeric(), overlap_genes = character())
  if (!length(sets) || !length(universe)) return(empty)
  result <- do.call(rbind, lapply(names(sets), function(term) {
    overlap <- intersect(query, sets[[term]])
    data.frame(term = term, overlap = length(overlap), set_size = length(sets[[term]]),
      query_size = length(query), universe_size = length(universe),
      p_value = if (length(overlap)) stats::phyper(length(overlap) - 1, length(sets[[term]]),
        length(universe) - length(sets[[term]]), length(query), lower.tail = FALSE) else 1,
      overlap_genes = paste(overlap, collapse = ";"))
  }))
  result$adjusted_p <- stats::p.adjust(result$p_value, method = "BH")
  result[order(result$adjusted_p, result$p_value), ]
}
stage_enrich <- function(cfg) {
  universes <- readRDS(file.path(cfg$output_dir, "objects/marker_universes.rds"))
  markers <- utils::read.csv(file.path(cfg$output_dir, "tables/markers_positive.csv"), colClasses = c(cluster = "character"))
  sets <- if (is.null(cfg$enrichment$gmt)) list() else read_gmt(cfg$enrichment$gmt)
  tables <- list(); audit <- list()
  for (cluster in names(universes)) {
    m <- markers[markers$cluster == cluster, ]; m <- m[order(-m$avg_log2FC), ]
    query <- head(m$gene, cfg$enrichment$top_genes)
    result <- ora(query, universes[[cluster]], sets, cfg$enrichment$min_size, cfg$enrichment$max_size)
    if (nrow(result)) { result$cluster <- cluster; tables[[cluster]] <- result }
    audit[[cluster]] <- list(query = query, tested_gene_universe = universes[[cluster]], eligible_pathways = nrow(result))
  }
  all <- if (length(tables)) do.call(rbind, tables) else data.frame(cluster = character(), term = character(), adjusted_p = numeric())
  top <- all[all$adjusted_p < 0.05, , drop = FALSE]
  if (nrow(top)) top <- do.call(rbind, lapply(split(top, top$cluster), head, 5))
  plot <- if (nrow(top)) {
    top$display <- paste0("C", top$cluster, ": ", substr(top$term, 1, 64))
    ggplot2::ggplot(top, ggplot2::aes(-log10(pmax(adjusted_p, 1e-300)), stats::reorder(display, -log10(pmax(adjusted_p, 1e-300))))) +
      ggplot2::geom_point(ggplot2::aes(size = overlap, color = cluster)) + ggplot2::theme_bw() +
      ggplot2::labs(x = "-log10(BH-adjusted p)", y = NULL, title = "Reactome marker over-representation")
  } else ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
    label = if (length(sets)) "No pathways passed BH < 0.05" else "Enrichment disabled") + ggplot2::theme_void()
  c(write_csv(all, file.path(cfg$output_dir, "tables/enrichment_all.csv")),
    write_csv(top, file.path(cfg$output_dir, "tables/enrichment_top.csv")),
    write_json(list(status = if (length(sets)) "completed" else "disabled", tested_pathways = nrow(all),
      significant_pathways = sum(all$adjusted_p < 0.05), multiple_testing = "BH over every eligible pathway within each cluster"),
      file.path(cfg$output_dir, "tables/enrichment_status.json")),
    save_rds(audit, file.path(cfg$output_dir, "objects/enrichment_universes.rds")),
    save_plot(plot, file.path(cfg$output_dir, "figures/enrichment.png"), 13, max(5, nrow(top) * 0.24)))
}
