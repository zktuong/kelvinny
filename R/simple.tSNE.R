#' Plotting and coloring tSNE plot.
#' 
#' @param gene Gene name.
#' @param data tSNE data table containing coordinate and cell information.
#' @param exprs Expression matrix with genes in rows and cells in columns.
#' @param alpha Transparency for geom_point. Default is "0.5".
#' @param size Size of geom_point. Default is "0.5".
#' @param orange TRUE/FALSE. TRUE is orange/purple palette. FALSE is red/green palette
#' @return Plotting and coloring tSNE plot.
#' @examples
#' simple.tSNE('Il33', orange=TRUE) # plots which cells express > 0 transcript counts using orange and purple palette
#' @export
simple.tSNE <- function(gene, data = dat3d, exprs = exp.unlog.cln, alpha = 0.5, size = 1.5, orange = FALSE){
  if(orange == FALSE ) {color_scale <- c("#c91212","#0e9122", "#D3D3D3")
} else {color_scale <- c("#ff7f00", "#984ea3", "#D3D3D3")}

  gene_idx <- which(rownames(exprs) == gene)
  cln.gene <- gsub('_.*', '', gene)
  pos_idx <- which(exprs[gene_idx, ] > 0 )

  WT <- grep('^1_*|^3_*', colnames(exprs), value = TRUE)
  TG <- grep('^2_*|^4_*', colnames(exprs), value = TRUE)
  WT_pos_idx <- pos_idx[which(names(pos_idx) %in% WT)]
  TG_pos_idx <- pos_idx[which(names(pos_idx) %in% TG)]

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
  
  p <- ggplot(data, aes(x, y))
  p <- p + 
    geom_point(aes(color = col_vec), alpha = alpha, size = size) + 
    if(length(unique(col_vec)) == 3) {
      scale_color_manual(values = color_scale)
    } else {
      scale_color_manual(values = c("#B3B3B3")) }  
  p <- p + labs(title = cln.gene) + 
    xlab('tSNE 1') +
    ylab('tSNE 2') +
    theme(plot.title = element_text(face = "bold", colour = "#990000", size = 16)) +
    theme(axis.text = element_text(size = 8)) +
    theme(axis.title = element_text(size = 12)) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    guides(color = 'none')
    
  return(p)
}