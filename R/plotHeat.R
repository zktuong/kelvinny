#' Plotting a heatmap using pheatmap
#'
#' @param data A data frame containing values of interest. Genes in rows, samples across columns. if supplied as a matrix with rownames as genes and the 1st column is exprs values, it will still go through.
#' @param color Default has "BluWhRd", "RdWhBlu", "viridis" or a string of the colours you want. e.g. c("blue","yellow","orange"). 
#' @param scale Default is "row". Other options include "column" or "none". see ?pheatmap.
#' @param genes if provided, will highlight the vector provided in the rows
#' @param n see ?viridis
#' @param alpha see?viridis
#' @param begin see ?viridis
#' @param end see ?viridis
#' @param option see ?viridis
#' @param direction see ?viridis
#' @param show_colnames TRUE/FALSE. see ?pheatmap.
#' @param show_rownames TRUE/FALSE. see ?pheatmap.
#' @param cluster_rows TRUE/FALSE. see ?pheatmap.
#' @param cluster_cols TRUE/FALSE. see ?pheatmap.
#' @param fontsize Default is 9. see ?pheatmap.
#' @param drop_levels TRUE/FALSE. see ?pheatmap.
#' @param legend TRUE/FALSE. see ?pheatmap.
#' @param border_color see ?pheatmap.
#' @param repel_degree amount of repel for flag. Only used if genes is specified
#' @param ... Passes on arguments to pheatmap.
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
#' @import grDevices
#' @import viridisLite
#' @export
plotHeat <- function(d , color = "BluWhRd", scale = "row", genes = NULL,  n = 50, alpha = 1, begin = 0, end = 1,  option = "D", direction = 1, show_colnames = TRUE, show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, drop_levels = TRUE, fontsize = 9, legend = TRUE, border_color = NA, repel_degree = 0, ...) {
  if(class(d[,1]) != "numeric") {
  rnames <- d[,1] # assign labels in column 1 to "rnames"
  mat_data <- data.matrix(d[,2:ncol(d)])    # transform column 2 to last column into a matrix
  rownames(mat_data) <- rnames
  } else {
  mat_data <- data.matrix(d)
  }

  if (length(color) > 1){
    col_scale = color
    col_palette <- grDevices::colorRampPalette(col_scale)(n = n)
  } else if (color == "RdWhBlu") {
    col_scale = RColorBrewer::brewer.pal(9, "RdBu")
    col_palette <- grDevices::colorRampPalette(col_scale)(n = n)
  } else if (color == "BluWhRd") {
    col_scale = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    col_palette <- grDevices::colorRampPalette(col_scale)(n = n)
  } else if (color == "viridis") {
    col_scale = viridisLite::viridis(n, alpha = alpha, begin = begin, end = end, direction = direction, option = option)
    col_palette = col_scale
  }
  p <- pheatmap::pheatmap(
    mat               = mat_data,
    color             = col_palette,
    border_color      = border_color,
    scale = scale,
    legend = legend,
    show_colnames     = show_colnames,
    show_rownames     = show_rownames,
    drop_levels       = drop_levels,
    fontsize          = fontsize,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    ...
  )

  add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {

  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  require(grid)
  heatmap <- pheatmap$gtable

  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 

  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")

  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant

    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }

      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }

    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))

    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)

  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions

  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                   grobs = new.flag,
                                   t = 4, 
                                   l = 4
  )

  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label

  # plot result
  grid.newpage()
  grid.draw(heatmap)

  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

  if (!is.null(genes)){
  p <- add.flag(p,
         kept.labels = genes,
         repel.degree = repel_degree)
                       }

  return(p)
}
