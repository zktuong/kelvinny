#' shortcut to create directory
#' 
#' @param directory Directory to create
#' @return creates directory if not found
#' @examples
#' createDir('./plots/test/test2/test3/test4')
#' @export
#'      
            
createDir <- function(directory) 
{
	if (!dir.exists(directory))
	{
	print(paste0("creating path to ", file.path(directory)))
    dir.create(directory, recursive = T)
	}
}