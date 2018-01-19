#' Plotting a heatmap using pheatmap
#'
#' @param data A data frame containing values of interest. Genes in rows, samples across columns.
#' @param color A string for what colours you want. e.g. c("blue","white,"red")
#' @param scale see ?pheatmap
#' @param show_colnames see ?pheatmap
#' @param show_rownames see ?pheatmap
#' @param cluster_rows see ?pheatmap
#' @param cluster_cols see ?pheatmap
#' @param fontsize see ?pheatmap
#' @param n see ?viridis
#' @param option see ?viridis
#' @param direction see ?viridis
#' @return A wrapper for pheatmap.
#' @examples
#' data(iris)
#' head(iris)
#' #   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#' #1          5.1         3.5          1.4         0.2  setosa
#' #2          4.9         3.0          1.4         0.2  setosa
#' #3          4.7         3.2          1.3         0.2  setosa
#' #4          4.6         3.1          1.5         0.2  setosa
#' #5          5.0         3.6          1.4         0.2  setosa
#' #6          5.4         3.9          1.7         0.4  setosa
#' df <- iris[,c(1:4)]
#' plotHeat(df)
#' @import pheatmap
#' @import RColorBrewer
#' @import viridis
#' @export
plotHeat <- function(dm , color="RdWhBlu", scale="row", show_colnames=TRUE, show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=TRUE, fontsize=9, n=50, option="D", direction=1) {
  if(dm == null)
  if(class(dm[,1]) != "numeric") {
  rnames <- dm[,1] # assign labels in column 1 to "rnames"
  mat_data <- data.matrix(dm[,2:ncol(m)])    # transform column 2 to last column into a matrix
  rownames(mat_data) <- rnames
  } else {
  mat_data <- data.matrix(dm)
  }

  if (color == "RdWhBlu") {
    col_scale = c(brewer.pal(9, "RdBu"))
    col_palette <- colorRampPalette(col_scale)(n = n)
  } else if (color == "viridis") {
    col_scale = viridis(n, option=option, direction=direction)
    col_palette = col_scale
    } else if (color == color) {
    col_scale = color
    col_palette <- colorRampPalette(col_scale)(n = n)
  }
  p <- pheatmap::pheatmap(
    mat               = mat_data,
    color             = col_palette,
    border_color      = NA,
    scale = scale,
    legend = TRUE,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    drop_levels       = TRUE,
    fontsize          = fontsize,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols
  )
  return(p)
}
