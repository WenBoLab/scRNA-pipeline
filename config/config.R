config <- list(
  project = "pbmc3k", seed = 42L, output_dir = "results/pbmc3k",
  input = list(mode = "demo", samples = "config/samples.csv"),
  qc = list(min_genes = 200L, max_genes = 2500L, max_pct_mt = 5,
            min_cells_per_gene = 3L, mitochondrial_pattern = "^MT-", overrides = list()),
  doublets = list(enabled = TRUE, remove = TRUE, expected_rate = 0.04),
  normalization = list(scale_factor = 10000, n_hvg = 2000L),
  embedding = list(n_pcs = 30L, dims = 1:20, k = 20L, resolution = 0.5,
                   integration = "none", batch_key = "batch"),
  markers = list(min_pct = 0.1, min_log2fc = 0.25, max_adjusted_p = 0.05),
  annotation = list(panel = "config/markers.R", manual_labels = NULL,
                    min_score = 0.25, min_margin = 0.10),
  enrichment = list(gmt = "resources/reactome_2022.gmt", min_size = 5L,
                    max_size = 500L, top_genes = 200L)
)
