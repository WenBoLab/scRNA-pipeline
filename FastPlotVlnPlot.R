#' Plot combined QC violin plots
#'
#' This function creates a combined violin plot showing QC metrics (nCount_RNA,
#' nFeature_RNA, percent.mt, percent.ribo) across clusters.
#'
#' @param object A Seurat object with seurat_clusters in meta.data.
#'
#' @return A combined patchwork/ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- FastPlotVlnPlot(sce)
#' ggsave("qc_violin.pdf", p, width = 8, height = 12)
#' }
FastPlotVlnPlot <- function(object = NULL) {
  # Load required packages
  nds <- c("Seurat", "cowplot", "ggplot2")
  Plus.library(nds)

  p4 <- VlnPlot(object, features = 'nCount_RNA', ncol = 1, pt.size = 0, log = TRUE, same.y.lims = TRUE) +
    geom_boxplot(outlier.size = 0, width = 0.3, show.legend = FALSE, notchwidth = 0.1) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )

  p5 <- VlnPlot(object, features = 'nFeature_RNA', ncol = 1, pt.size = 0, log = TRUE) +
    geom_boxplot(outlier.size = 0, width = 0.3, show.legend = FALSE, notchwidth = 0.1) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )

  p6 <- VlnPlot(object, features = 'percent.mt', ncol = 1, pt.size = 0, y.max = 100) +
    geom_boxplot(outlier.size = 0, width = 0.3, show.legend = FALSE, notchwidth = 0.1) +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 1),
      legend.position = "none"
    )

  p7 <- VlnPlot(object, features = 'percent.ribo', ncol = 1, pt.size = 0, y.max = 100) +
    geom_boxplot(outlier.size = 0, width = 0.3, show.legend = FALSE, notchwidth = 0.1) +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 1),
      legend.position = "none"
    )

  plot <- cowplot::plot_grid(p4, p5, p6, p7, nrow = 4)

  return(plot)
}
