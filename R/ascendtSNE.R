#' Plotting and coloring ascend tSNE plots when samples/batch 1 and 2 are WT and TG respectively..
#'
#' @param gene Gene name.
#' @param data tSNE data table containing coordinate and cell information.
#' @param exprs Expression matrix with genes in rows and cells in columns.
#' @param cutOff Expression value cut off threshold. Default is "0".
#' @param alpha Transparency for geom_point. Default is "0.5".
#' @param cell_size Size of geom_point. Default is "0.5".
#' @param show_branch_points See ?plot_cell_trajectory. Default is FALSE.
#' @param option see ?viridis. Default is "D".
#' @param begin see ?viridis. Default is "0.2".
#' @param end see ?viridis. Default is "0.5".
#' @param direction see ?viridis. Default is "1".
#' @param lab_size Size of plot labels. Default is "4".
#' @param heat Plot as heatmap? TRUE/FALSE
#' @param scale_to_expr TRUE/FALSE. True applies a log2 transformation to the counts. Default is FALSE.
#' @param manual_col Color option. Default is FALSE. Replace with string of colours or choose preset options. Options are: FALSE = viridis colour scale, TRUE = ggplot default palette, OrPU = Orange and purple, GrRd=Green and Red.
#' @return Plotting and coloring ascend tSNE plots.
#' @examples
#' em.set <- readRDS('em.set.RDS')
#' tsne.plot <- PlotTSNE(em.set, PCA = TRUE, condition = "cluster", seed = 122, perplexity = 30, theta = 0.5) 
#' # expression matrix
#' exp.unlog.cln <- as.matrix(GetExpressionMatrix(em.set))
#' # tSNE coordinates and cell information from PlotTSNE and em.set from ascend
#' dat3d <- cbind(tsne.plot$data,em.set@CellInformation)
#' colnames(dat3d) <- c('x','y','conditions','cell_barcode','batch','cluster')
#' ascendtSNE.gene('Cd207')            # plots which cells express > 0 transcript counts
#' ascendtSNE.gene('Cd207', heat=TRUE) # overlay heatmap of expression level over tSNE plot
#' ascendtSNE.info('cluster')          # plots cluster information
#' ascendtSNE.info('batch')            # plots batch/sample information
#' @import ggplot2
#' @import viridis
#' @import ascend
#' @import grid
#' @export
ascendtSNE.gene <- function(gene, data = dat3d, exprs = exp.unlog.cln, cutOff= 0, alpha = 0.5, size = 1.5, option="D", begin=0.2, end=0.5, direction=1, lab_size=4, heat=FALSE, manual_col=FALSE){
  gene_idx <- which(rownames(exprs) == gene)
  pos_idx <- which(exprs[gene_idx, ] > cutOff)

  WT <- grep('*-1', colnames(exprs), value = TRUE)
  TG <- grep('*-2', colnames(exprs), value = TRUE)
  WT_pos_idx <- pos_idx[which(names(pos_idx) %in% WT)]
  TG_pos_idx <- pos_idx[which(names(pos_idx) %in% TG)]

  WTpercent <- length(WT_pos_idx)/length(WT)*100
  TGpercent <- length(TG_pos_idx)/length(TG)*100

  cat('Total WT cells:', ' ', length(WT), "\n")
  cat('Total TG cells:', ' ', length(TG), "\n")
  cat('Positive WT cells:', ' ',length(WT_pos_idx),' ')
  cat('Percentage of Total WT:', ' ',WTpercent, '%',"\n")
  cat('Positive TG cells:', ' ',length(TG_pos_idx),' ')
  cat('Percentage of Total TG:', ' ',TGpercent, '%',"\n")

if (heat == FALSE) {
    if(length(pos_idx) != 0) {
    col_vec <- rep('WT_positive', length(exprs[1, ]))
    col_vec[-WT_pos_idx] <- 'zero'
    col_vec[+TG_pos_idx] <- 'TG_positive'
  } else {
    col_vec <- rep('zero', length(exprs[1, ]))
  }
  
  col_vec <- as.factor(col_vec)
  data$col_vec <- col_vec
  data <- data[rev(order(data$col_vec)),]
 # levels(col_vec) <- c("WT_positive","TG_positive","zero")
     
  n = length(unique(col_vec))
  gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
  
  p <- ggplot(data, aes(x, y)) +
    geom_point(aes(color = col_vec), alpha = alpha, size = size) + 
    if(manual_col == "OrPu" && length(unique(col_vec)) == 3) {
		colour_scale = c("#ff7f00", "#984ea3", "#e0e0e0")
		scale_color_manual(values = colour_scale)
	} else if(manual_col == "GrRd" && length(unique(col_vec)) == 3) {
		colour_scale = c("#c91212","#0e9122", "#e0e0e0")
		scale_color_manual(values = colour_scale)
	} else if(manual_col == TRUE && length(unique(col_vec)) == 3) {
		colour_scale = gg_color_hue(n-1)
		scale_color_manual(values = c(colour_scale,"#e0e0e0"))
	} else if (manual_col == FALSE && length(unique(col_vec)) > 1) {
		colour_scale = viridis(n-1, option=option, begin=begin, end=end, direction=direction)
		scale_color_manual(values = c(colour_scale, "#e0e0e0"))
	} else if (manual_col != FALSE && manual_col != TRUE && length(unique(col_vec)) > 1) {
		colour_scale = manual_col
		scale_color_manual(values = c(colour_scale,"#e0e0e0"))
	} else if (length(unique(col_vec)) == 2) {
		colour_scale =c("#A020F0", "#e0e0e0")
		scale_color_manual(values = colour_scale)  
   	} else if (length(unique(col_vec)) == 1) {
		colour_scale =c("#B3B3B3","#B3B3B3","#B3B3B3","#B3B3B3","#B3B3B3")
    scale_color_manual(values = colour_scale) } 

  p <- p + labs(title = gene) + 
    xlab('tSNE 1') +
    ylab('tSNE 2') +
    theme(plot.title = element_text(face = "bold", colour = "#990000", size = 16)) +
    theme(axis.text = element_text(size = 8)) +
    theme(axis.title = element_text(size = 12)) +
    #theme_classic() +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom") + 
    guides(color = 'none') +
    annotate("text", hjust = 1, x=max(data$x), y=min(data$y)*0.9, label=paste('WT: ', length(WT_pos_idx), ' ',formatC(WTpercent, digits = 2,format="f"), '%',"\n"), size = lab_size, color=colour_scale[2])+
    annotate("text", hjust = 1, x=max(data$x), y=min(data$y)*0.999, label=paste('TG: ', length(TG_pos_idx), ' ',formatC(TGpercent,digits = 2,format="f"), '%',"\n"), size = lab_size, color=colour_scale[1])

    return(p)
} else {
    data$gene.exp <- exprs[gene_idx,]
    data <- data[order(data$gene.exp),]
    p <- ggplot(data, aes(x, y)) +
    geom_point(aes(color=data$gene.exp), alpha = alpha, size = size) + 
    scale_color_gradientn(colours = c("#e0e0e0",viridis(50, option=option, begin=begin*0, end=end/end, direction=direction))) +
    labs(title = gene) + 
    xlab('tSNE 1') +
    ylab('tSNE 2') +
    theme(plot.title = element_text(face = "bold", colour = "#990000", size = 16)) +
    theme(axis.text = element_text(size = 8)) +
    theme(axis.title = element_text(size = 12)) +
    #theme_classic() +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=7), legend.key.size=unit(2,"pt")) +
    guides(colour= guide_legend(title = "Transcript Count", title.position = "top", title.hjust=0.5,override.aes = list(alpha = 1)))
}
  return(p)
} 

#' @export
ascendtSNE.info <- function(condition="cluster", data = dat3d, alpha = 0.5, size = 1.5, option="D", begin=0.2, end=0.5, direction=1, manual_col=FALSE){
  options(warn=-1)
  library(ggplot2)
  library(viridis)
  library(grid)

  data$batch <- as.character(data$batch)
  data$cluster <- as.character(data$cluster)
  
  if(condition == "batch") {
  n = length(unique(data$batch))
  } 
  if(condition == "cluster") {
  n = length(unique(data$cluster))
  } 

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

p <- ggplot(data, aes(x, y)) + 
	xlab('tSNE 1') +
    ylab('tSNE 2') +
    theme(axis.text = element_text(size = 8)) +
    theme(axis.title = element_text(size = 12)) +
    #theme_classic() +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
if(manual_col == FALSE & condition == "batch") {
		colour_scale = viridis(n, option=option, begin=begin, end=end, direction=direction)
		p <- p + geom_point(aes(color=data$batch), alpha = alpha, size = size) + scale_color_manual(values = colour_scale, guide = guide_legend(title = "batch"))
	} else if(manual_col== FALSE & condition == "cluster") {
		colour_scale = viridis(n)
		p <- p + geom_point(aes(color=data$cluster), alpha = alpha, size = size) + scale_color_manual(values = colour_scale, guide = guide_legend(title = "cluster"))
	} else if(manual_col == TRUE & condition == "batch") {
		colour_scale = gg_color_hue(n)
		p <- p + geom_point(aes(color=data$batch), alpha = alpha, size = size) + scale_color_manual(values = colour_scale, guide = guide_legend(title = "batch"))
	} else if(manual_col != TRUE & manual_col != FALSE & condition == "batch") {
		colour_scale = manual_col
		p <- p + geom_point(aes(color=data$batch), alpha = alpha, size = size) + scale_color_manual(values = colour_scale, guide = guide_legend(title = "batch"))
	} else if(manual_col == TRUE & condition == "cluster") {
		colour_scale = gg_color_hue(n)
		p <- p + geom_point(aes(color=data$cluster), alpha = alpha, size = size) + scale_color_manual(values = colour_scale, guide = guide_legend(title = "cluster"))
	} else if(manual_col != TRUE & manual_col != FALSE & condition == "cluster") {
		colour_scale = manual_col
		p <- p + geom_point(aes(color=data$cluster), alpha = alpha, size = size) + scale_color_manual(values = colour_scale, guide = guide_legend(title = "cluster"))
	}
	
if(condition == "batch") {
      p <- p + labs(title = "Labelled by batch") + theme(plot.title = element_text(face = "bold", colour = "#000000", size = 16))
	} else if(condition == "cluster") {
	p <- p + labs(title = "Labelled by cluster") + theme(plot.title = element_text(face = "bold", colour = "#000000", size = 16))
	} 

return(p)
}
