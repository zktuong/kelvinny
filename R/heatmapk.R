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
#' @param ... Passes on arguments to pheatmap or viridis
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
#'
#' # using example from pheatmap
#' # Create test matrix
#' test = matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#' # Generate annotations for rows and columns
#' annotation_col = data.frame(
#'   CellType = factor(rep(c("CT1", "CT2"), 5)),
#'   Time = 1:5
#' )
#' rownames(annotation_col) = paste("Test", 1:10, sep = "")
#'
#' annotation_row = data.frame(
#'   GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
#' )
#' rownames(annotation_row) = paste("Gene", 1:20, sep = "")
#' #plot!
#' plotHeat(test, annotation_col=annotation_col, annotation_row=annotation_row)
#' @import pheatmap
#' @import RColorBrewer
#' @import viridis
#' @export
plotHeat <- function(d , color = "RdWhBlu", scale = "row", n = 50, alpha = 1, begin = 0, end = 1,  option = "D", direction = 1, show_colnames=TRUE, show_rownames=TRUE, cluster_rows=TRUE, cluster_cols=TRUE, fontsize=9, ...) {
  if(class(d[,1]) != "numeric") {
  rnames <- d[,1] # assign labels in column 1 to "rnames"
  mat_data <- data.matrix(d[,2:ncol(d)])    # transform column 2 to last column into a matrix
  rownames(mat_data) <- rnames
  } else {
  mat_data <- data.matrix(d)
  }

  if (color == "RdWhBlu") {
    col_scale = c(brewer.pal(9, "RdBu"))
    col_palette <- colorRampPalette(col_scale)(n = n)
  } else if (color == "viridis") {
    col_scale = viridis(n, alpha = alpha, begin = begin, end = end, direction = direction, option = option)
    col_palette = col_scale
    } else if (color == color) {
    col_scale = color
    col_palette <- colorRampPalette(col_scale)(n = n)
  }
  p <- pheatmap(
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
    cluster_cols = cluster_cols,
    ...
  )
  return(p)
}
