#' shortcut to initiate building of kelvinny
#' 
#' @return runs init kelvinny
#' @examples
#' kelvinny()
#' @export
#'      
            
kelvinny <- function() 
{
	setwd("~/Documents/GitHub/kelvinny")
	library(kelvinny)
	library(roxygen2)
	init_Kelvinny()
}