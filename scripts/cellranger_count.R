# Optional external FASTQ alignment. Not executed in the bundled matrix demo.
if (!file.exists("scrna-seq-learning.Rproj")) stop("Run from the project root.")
executable <- Sys.which("cellranger")
if (!nzchar(executable)) stop("Install Cell Ranger separately and add it to PATH.")
samples <- read.csv("config/cellranger.csv", colClasses = "character")
required <- c("sample_id", "fastqs", "transcriptome", "sample_name")
stopifnot(all(required %in% names(samples)), nrow(samples) > 0,
  !anyNA(samples[required]), !anyDuplicated(samples$sample_id))
dir.create("results/cellranger", recursive = TRUE, showWarnings = FALSE)
project_root <- getwd(); setwd("results/cellranger")
tryCatch({
  for (i in seq_len(nrow(samples))) {
    s <- samples[i, ]
    if (!grepl("^[A-Za-z][A-Za-z0-9_-]*$", s$sample_id)) stop("Use a safe Cell Ranger sample ID.")
    if (file.exists(s$sample_id) || dir.exists(s$sample_id)) stop("Refusing to overwrite an existing Cell Ranger run.")
    if (!dir.exists(s$fastqs) || !dir.exists(s$transcriptome)) stop("Use existing absolute FASTQ/reference directory paths.")
    arguments <- c("count", paste0("--id=", s$sample_id),
      paste0("--fastqs=", shQuote(s$fastqs)), paste0("--transcriptome=", shQuote(s$transcriptome)),
      paste0("--sample=", shQuote(s$sample_name)), "--localcores=4", "--localmem=32", "--create-bam=false")
    status <- system2(executable, arguments, stdout = paste0(s$sample_id, ".log"), stderr = paste0(s$sample_id, ".err"))
    if (status != 0L) stop("Cell Ranger failed; inspect the sample logs.")
  }
}, finally = setwd(project_root))
