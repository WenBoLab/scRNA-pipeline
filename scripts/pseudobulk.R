# Edit these values for a real study with independent donors and reviewed labels.
source("scripts/bootstrap.R")
pseudobulk_de(
  object_path = "results/my_study/objects/06_annotated.rds",
  output_dir = "results/my_study/pseudobulk_T",
  cell_type = "T cells", case = "case", control = "control",
  paired = FALSE, min_cells = 20L
)
