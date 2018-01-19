#' Plotting and coloring Monocle 2 trajectory plot.
#' 
#' @param gene Gene name.
#' @param data Monocle CellDataSet for the experiment.
#' @param exprs Expression matrix with genes in rows and cells in columns.
#' @param manual_col Color option. Default is FALSE. Replace with string of colours or choose preset options. Options are: FALSE = viridis colour scale, TRUE = ggplot default palette, OrPU = Orange and purple, GrRd=Green and Red.
#' @param cutOff Expression value cut off threshold. Default is "0".
#' @param cell_size Size of geom_point. Default is "0.5".
#' @param show_branch_points See ?plot_cell_trajectory. Default is FALSE.
#' @param option see ?viridis. Default is "D".
#' @param begin see ?viridis. Default is "0.2".
#' @param end see ?viridis. Default is "0.5".
#' @param direction see ?viridis. Default is "1".
#' @param scale_to_expr TRUE/FALSE. True applies a log2 transformation to the counts. Default is FALSE.
#' @return Plotting and coloring Monocle 2 trajectory plot.
#' @examples
#' colorMyTrajectory('Il33', manual_col="OrPu")
#' colorMyTrajectory('Il33', manual_col="GrRd", show_branch_points=TRUE)
#' @import ggplot2
#' @import viridis
#' @import monocle
#' @import ascend
#' @export
colorMyTrajectory <- function(gene, data = HSMM, exprs = exp.unlog.cln, manual_col=FALSE, cutOff= 0, cell_size = 0.5, show_branch_points = FALSE, option="D", begin=0.2, end=0.5, direction=1, log2=FALSE){
if (log2 == TRUE) {use_cell_size <- log2(exp.unlog.cln[which(rownames(exp.unlog.cln) == gene),])
} else {use_cell_size <- cell_size}

rownames(exprs) <- gsub('_.*', '', rownames(exprs))
gene_idx <- which(rownames(exprs) == gene)
pos_idx <- which(exprs[gene_idx, ] > cutOff)
new_column <- ncol(pData(data))+1

WT <- grep('^1_*|^3_*', colnames(exprs), value = TRUE)
TG <- grep('^2_*|^4_*', colnames(exprs), value = TRUE)
WT_pos_idx <- pos_idx[which(names(pos_idx) %in% WT)]
TG_pos_idx <- pos_idx[which(names(pos_idx) %in% TG)]

if(length(pos_idx) != 0) {
    pData(data)[,new_column] <- 'WT_positive'
    colnames(pData(data))[new_column] <- gene
    pData(data)[-WT_pos_idx, new_column] <- 'zero'
	  pData(data)[+TG_pos_idx, new_column] <- 'TG_positive'    
  } else {
    pData(data)[,new_column] <- 'zero'
  }

p <- monocle::plot_cell_trajectory(data, color_by=gene, cell_size = use_cell_size, show_branch_points = show_branch_points)

n <- length(unique(pData(data)[,new_column]))

p <- p + 
  if(manual_col == "OrPu" && n == 3 ) {
      colour_scale = c("#ff7f00", "#984ea3", "#e0e0e0")
      ggplot2::scale_color_manual(values = colour_scale)
    } else if(manual_col == "GrRd" && n == 3) {
      colour_scale = c("#c91212","#0e9122", "#e0e0e0")
      ggplot2::scale_color_manual(values = colour_scale)
    } else if(manual_col == TRUE && n == 3) {
      colour_scale = gg_color_hue(n-1)
      ggplot2::scale_color_manual(values = c(colour_scale,"#e0e0e0"))
    } else if (manual_col == FALSE && n > 1) {
      colour_scale = viridis::viridis(n-1, option=option, begin=begin, end=end, direction=direction)
      ggplot2::scale_color_manual(values = c(colour_scale, "#e0e0e0"))
    } else if (manual_col != FALSE && manual_col != TRUE && n > 1) {
      colour_scale = manual_col
      ggplot2::scale_color_manual(values = c(colour_scale,"#e0e0e0"))
    } else if(n == 2) {
      colour_scale = c("#A020F0", "#e0e0e0")
      ggplot2::scale_color_manual(values = colour_scale)
    } else if (n == 1) {
      colour_scale =c("#B3B3B3","#B3B3B3","#B3B3B3","#B3B3B3","#B3B3B3")
      ggplot2::scale_color_manual(values = colour_scale) } 

p$data[,ncol(p$data)] <- as.factor(p$data[,ncol(p$data)])
p$data <- p$data[rev(order(p$data[,ncol(p$data)])),]  

p <- p + ggplot2::labs(title = gene) + ggplot2::guides(color = 'none')
return(p)
}


#monocle.factor <- function(x, data=HSMM, cell_size = 0.5, show_branch_points = TRUE){
#p.t <- plot_cell_trajectory(data, color_by = x, cell_size = cell_size, show_branch_points = show_branch_points)
#p.t <- p.t + ggplot2::labs(title = x) + ggplot2::guides(color = 'none') 
#return(p.t)
#}

#monocle.trajectory.split <- function(x, data=HSMM, cell_size = 0.5, show_branch_points = TRUE, nrow=5){
#class(pData(data)[,x]) <- 'factor'
#color1 <- brewer.pal(9,"Set3")
#color2 <- brewer.pal(9,"Paired")
#color <- c(color1,color2)
#p.ts <- plot_cell_trajectory(data, color_by= x, cell_size=cell_size , show_branch_points = show_branch_points) + facet_wrap(x, nrow=nrow)
#p.ts <- p.ts + ggplot2::guides(color = 'none') + ggplot2::scale_color_manual(values = color[1:length(x)])
#return(p.ts)
#}