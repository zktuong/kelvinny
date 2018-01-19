#' Plotting a heatmap using pheatmap
#'
#' @param data A data frame containing values of interest. Genes in rows, samples across columns.
#' @param color A string for what colours you want. e.g. c("blue","white,"red")
#' @param scale see ?pheatmap
#' @param show_colnames see ?pheatmap
#' @param show_rownames see ?pheatmap
#' @param cluster_rows see ?pheatmap
#' @param clsuter_cols see ?pheatmap
#' @param fontsize see ?pheatmap
#' @param n see ?viridis
#' @param option see ?viridis
#' @param direction see ?viridis
#' @return A wrapper for pheatmap.
#' @examples
#' data(iris)
#' data <- iris[,c(1:4)]
#' pheatmapping(data)
#' @import pheatmap
#' @import RColorBrewer
#' @import viridis
#' @export
heatmap.k <- function(data , color="RdWhBlu", scale="row", show_colnames=TRUE, show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=TRUE, fontsize=9, n=50, option="A", direction=1) {
  if(class(data[,1]) != "numeric") {
  rnames <- data[,1] # assign labels in column 1 to "rnames"
  mat_data <- data.matrix(data[,2:ncol(data)])    # transform column 2 to last column into a matrix
  rownames(mat_data) <- rnames
  } else {
  mat_data <- data.matrix(data)
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
