aggregate_pseudobulk <- function(counts, metadata, min_cells = 20L) {
  fields <- c("donor", "condition", "cell_type")
  assert(all(fields %in% names(metadata)) && identical(colnames(counts), rownames(metadata)), "Pseudobulk metadata are misaligned.")
  assert(!anyNA(metadata[fields]) && all(vapply(metadata[fields], function(x) all(nzchar(as.character(x))), logical(1))), "Pseudobulk metadata cannot be empty.")
  # Length prefixes make the group key unambiguous even when labels contain separators.
  key <- do.call(paste, c(lapply(metadata[fields], function(x) paste0(nchar(as.character(x)), ":", x)), sep = "|"))
  groups <- unique(key); indices <- split(seq_len(ncol(counts)), factor(key, levels = groups))
  indices <- indices[lengths(indices) >= min_cells]
  assert(length(indices) > 0, "No donor/cell-type group passes min_cells.")
  pb <- do.call(cbind, lapply(indices, function(i) Matrix::rowSums(counts[, i, drop = FALSE])))
  rownames(pb) <- rownames(counts); colnames(pb) <- paste0("PB", seq_along(indices))
  meta <- do.call(rbind, lapply(indices, function(i) data.frame(metadata[i[1], fields, drop = FALSE], n_cells = length(i))))
  rownames(meta) <- colnames(pb)
  list(counts = pb, metadata = meta)
}
validate_donor_design <- function(meta, case, control, paired = FALSE) {
  assert(case != control && all(c(case, control) %in% meta$condition), "Both distinct conditions are required.")
  assert(!anyDuplicated(meta[c("donor", "condition")]), "Technical captures must be pooled before modeling.")
  for (condition in c(control, case)) assert(length(unique(meta$donor[meta$condition == condition])) >= 3,
    paste("At least three independent donors are required in", condition))
  if (paired) assert(all(table(meta$donor, meta$condition) == 1), "Paired analysis needs complete donor pairs.") else
    assert(!anyDuplicated(meta$donor), "Repeated donors across conditions require paired = TRUE.")
  invisible(TRUE)
}
pseudobulk_de <- function(object_path, output_dir, cell_type, case, control, paired = FALSE, min_cells = 20L) {
  object <- readRDS(object_path); meta <- object[[]]
  assert(all(c("cell_type", "annotation_status", "donor", "condition") %in% names(meta)), "Missing reviewed labels or donor metadata.")
  keep <- meta$cell_type == cell_type & meta$condition %in% c(case, control)
  assert(any(keep) && all(meta$annotation_status[keep] == "reviewed"), "Pseudobulk DE requires reviewed labels for the selected cells.")
  pb <- aggregate_pseudobulk(raw_counts(object)[, keep, drop = FALSE], meta[keep, , drop = FALSE], min_cells)
  meta <- pb$metadata; validate_donor_design(meta, case, control, paired)
  meta$condition <- factor(meta$condition, levels = c(control, case)); meta$donor <- factor(meta$donor)
  design <- if (paired) ~donor + condition else ~condition
  matrix <- stats::model.matrix(design, meta)
  assert(qr(matrix)$rank == ncol(matrix), "The donor design is not full rank.")
  keep_genes <- rowSums(pb$counts >= 10) >= 3
  assert(sum(keep_genes) >= 20, "Too few genes pass the pseudobulk expression filter.")
  counts <- round(pb$counts[keep_genes, , drop = FALSE])
  assert(all(counts <= .Machine$integer.max), "Pseudobulk counts exceed the supported integer range.")
  storage.mode(counts) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design)
  dds <- DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE)
  result <- as.data.frame(DESeq2::results(dds, contrast = c("condition", case, control)))
  result$gene <- rownames(result); result <- result[order(result$padj, na.last = TRUE), ]; rownames(result) <- NULL
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  save_rds(pb, file.path(output_dir, "pseudobulk_counts.rds"))
  save_rds(dds, file.path(output_dir, "deseq2_model.rds"))
  write_csv(data.frame(profile = rownames(meta), meta), file.path(output_dir, "design.csv"))
  write_csv(result, file.path(output_dir, "differential_expression.csv"))
  write_json(list(case = case, control = control, cell_type = cell_type, paired = paired,
    design = paste(deparse(design), collapse = " "), donors = length(unique(meta$donor)),
    expression_filter = "At least 10 summed counts in at least 3 donor profiles",
    technical_captures = "Pooled within donor, condition and cell type"), file.path(output_dir, "analysis.json"))
  result
}
