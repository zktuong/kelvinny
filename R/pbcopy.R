#' a quick copy to/paste from clipboard. function
#' 
#' @param object anything object within R.
#' @return Copy-pasta.
#' @export
pbcopy <- function(object){
clip <- pipe("pbcopy", "w")
write.table(object, file = clip, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
close(clip)
}

#' @export
pbpaste <- function(unlist = TRUE){
  if(unlist == TRUE) {
data_string <- unlist(read.table(pipe("pbpaste"), sep="\t", header=F), use.names = FALSE, recursive = FALSE)
  } else {
data <- read.table(pipe("pbpaste"), sep="\t", header=F)
  }
  if(is.numeric(data_string) == FALSE | is.integer(data_string) == FALSE) {
data <- levels(droplevels(data_string))
 } else {
data <- data_string    
  }
return(data)
}
