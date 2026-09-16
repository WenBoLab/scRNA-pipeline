#' Load multiple R packages at once
#'
#' This function loads a vector of R packages, installing any that are missing
#' from the local library.
#'
#' @param pkgs A character vector of package names to load.
#'
#' @return Invisible. The packages are loaded (and installed if necessary).
#' @export
#'
#' @examples
#' \dontrun{
#' Plus.library(c("Seurat", "dplyr", "ggplot2"))
#' }
Plus.library <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
}
