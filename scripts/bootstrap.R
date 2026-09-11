# Load definitions only; this file never starts the analysis.
if (!file.exists("scrna-seq-learning.Rproj") || !dir.exists("R")) stop("Run from the project root or open scrna-seq-learning.Rproj.")
source("R/common.R", local = FALSE)
for (path in setdiff(list.files("R", pattern = "\\.R$", full.names = TRUE), "R/common.R")) source(path, local = FALSE)
options(stringsAsFactors = FALSE, Seurat.object.assay.version = "v5")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
