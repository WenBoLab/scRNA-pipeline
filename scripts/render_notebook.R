if (!file.exists("scrna-seq-learning.Rproj")) stop("Run from the project root.")
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Restore the R packages first.")
if (!rmarkdown::pandoc_available()) stop("Pandoc is required; RStudio includes it.")
dir.create("results/notebook", recursive = TRUE, showWarnings = FALSE)
rmarkdown::render("notebooks/01_pbmc3k_walkthrough.Rmd",
  output_file = "pbmc3k_walkthrough.html", output_dir = normalizePath("results/notebook"),
  knit_root_dir = getwd(), envir = new.env(parent = globalenv()), quiet = TRUE)
