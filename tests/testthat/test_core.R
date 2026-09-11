testthat::test_that("invalid count values and duplicate identifiers are rejected", {
  x <- matrix(c(1, 2, 3, 4), 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))
  testthat::expect_s4_class(validate_counts(x), "dgCMatrix")
  fractional <- x; fractional[1, 1] <- 0.5
  testthat::expect_error(validate_counts(fractional), "raw integer")
  negative <- x; negative[1, 1] <- -1
  testthat::expect_error(validate_counts(negative), "nonnegative")
  colnames(x) <- c("c", "c")
  testthat::expect_error(validate_counts(x), "unique")
})
testthat::test_that("bundled counts have their verified dimensions and exact totals", {
  withr::local_dir(PROJECT_ROOT); verify_demo()
  x <- read_capture("resources/pbmc3k", "10x_mtx")
  testthat::expect_identical(dim(x), c(32738L, 2700L))
  testthat::expect_equal(sum(x), 6390631)
  testthat::expect_equal(Matrix::nnzero(x), 2286884)
})
testthat::test_that("normalization preserves raw UMI values", {
  x <- matrix(1:60, nrow = 10, dimnames = list(paste0("g", 1:10), paste0("c", 1:6)))
  object <- SeuratObject::CreateSeuratObject(methods::as(x, "dgCMatrix"))
  before <- raw_counts(object)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  testthat::expect_equal(Matrix::nnzero(raw_counts(object) - before), 0)
  testthat::expect_true("data" %in% SeuratObject::Layers(object[["RNA"]]))
})
testthat::test_that("ORA uses the actual background and adjusts zero-overlap pathways", {
  result <- ora(c("a", "b", "outside"), letters[1:6], list(A = c("a", "b"), B = c("c", "d")), 1, 10)
  a <- result[result$term == "A", ]; b <- result[result$term == "B", ]
  testthat::expect_equal(a$p_value, 1/15)
  testthat::expect_equal(a$adjusted_p, 2/15)
  testthat::expect_equal(a$query_size, 2)
  testthat::expect_equal(b$p_value, 1)
  testthat::expect_equal(b$overlap, 0)
})
testthat::test_that("cytotoxic T evidence cannot be labeled NK solely from cytotoxic markers", {
  testthat::expect_identical(candidate_label("NK cells", 0.8, 0.2, 0.7, 0.3, 0.6, 0.25, 0.1), "Cytotoxic T cells (candidate)")
  testthat::expect_identical(candidate_label("NK cells", 0.8, 0.2, 0.01, 0.5, 0.6, 0.25, 0.1), "NK cells (candidate)")
  testthat::expect_identical(candidate_label("NK cells", 0.8, 0.2, 0.15, 0.4, 0.6, 0.25, 0.1), "T/NK unresolved")
})
testthat::test_that("manual review requires one label and evidence per cluster", {
  manual <- data.frame(cluster = c("0", "1"), cell_type = c("T cells", "B cells"), evidence = c("CD3D/TRAC", "MS4A1/CD79A"))
  testthat::expect_silent(validate_manual_labels(manual, c("0", "1")))
  testthat::expect_error(validate_manual_labels(manual[1, ], c("0", "1")), "exactly one")
  manual$evidence[1] <- ""
  testthat::expect_error(validate_manual_labels(manual, c("0", "1")), "evidence")
})
testthat::test_that("technical captures pool into donor sums without losing counts", {
  x <- matrix(1:12, nrow = 2, dimnames = list(c("g1", "g2"), paste0("c", 1:6)))
  meta <- data.frame(donor = rep(c("d1", "d2"), each = 3), condition = "control", cell_type = "T", row.names = colnames(x))
  pb <- aggregate_pseudobulk(methods::as(x, "dgCMatrix"), meta, min_cells = 1)
  testthat::expect_equal(ncol(pb$counts), 2)
  testthat::expect_equal(unname(pb$counts[, 1]), unname(rowSums(x[, 1:3])))
  testthat::expect_equal(sum(pb$counts), sum(x))
})
testthat::test_that("donor design prevents pseudoreplication and incomplete pairs", {
  meta <- data.frame(donor = paste0("d", 1:6), condition = rep(c("control", "case"), each = 3))
  testthat::expect_silent(validate_donor_design(meta, "case", "control"))
  testthat::expect_error(validate_donor_design(meta[-1, ], "case", "control"), "three independent")
  paired <- data.frame(donor = rep(paste0("d", 1:3), 2), condition = rep(c("control", "case"), each = 3))
  testthat::expect_silent(validate_donor_design(paired, "case", "control", paired = TRUE))
  testthat::expect_error(validate_donor_design(paired, "case", "control"), "paired")
})
testthat::test_that("resume validates output bytes and configuration fingerprint", {
  directory <- withr::local_tempdir(); output <- file.path(directory, "output.txt"); state <- file.path(directory, "cache.rds")
  writeLines("original", output); saveRDS(list(fingerprint = "a", outputs = file_hashes(output)), state)
  testthat::expect_true(stage_is_current(state, "a"))
  testthat::expect_false(stage_is_current(state, "b"))
  writeLines("changed", output)
  testthat::expect_false(stage_is_current(state, "a"))
})
testthat::test_that("README files are ASCII and the project ships R analysis sources", {
  paths <- list.files(PROJECT_ROOT, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  paths <- paths[!grepl("/(renv/library|renv/staging|results|\\.git)/", paths)]
  readmes <- paths[tolower(basename(paths)) == "readme.md"]
  testthat::expect_gt(length(readmes), 0)
  for (path in readmes) testthat::expect_true(all(as.integer(readBin(path, "raw", n = file.info(path)$size)) < 128))
  testthat::expect_length(paths[grepl("\\.(py|ipynb)$|/Snakefile$|/pyproject\\.toml$|/uv\\.lock$", paths)], 0)
})
