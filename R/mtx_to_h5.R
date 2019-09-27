#' Saves matrix to .h5 format quickly
#' 
#' @param object input matrix (or maybe even data frame) - can be sparse or dense, but this is to sort the problem with dense matrices
#' @param filename name of the h5 object to save to
#' @param datasetname name of dataset
#' @param ... pass to rhdf5 associated functions
#' @return Save large (dense) matrix to .h5 format
#' @examples
#' mtx_to_h5(matrixA, "matrix.h5")
#' mtx_to_h5totxt(matrixA, "matrix.h5")
#' @export
#'   

mtx_to_h5 <- function(object, filename, datasetname = "counts", ...){
	if(!grepl(".h5$", filename)){
		filename <- paste0(filename, ".h5")
	}
	cat(paste0("Creating ", crayon::magenta(filename)), sep = "\n")
	rhdf5::h5createFile(filename, ...)
	cat(crayon::green("Saving file "), sep = "\n")
	rhdf5::h5write(object, filename, datasetname, ...)
}

#' @export
mtx_to_h5totxt <- function(object, filename, datasetname = "counts",...){
	if(!grepl(".h5$", filename)){
		filename <- paste0(filename, ".h5")
	}
	cat(paste0("Creating ", crayon::magenta(filename)), sep = "\n")
	rhdf5::h5createFile(filename, ...)
	cat(crayon::green("Saving file "), sep = "\n")
	rhdf5::h5write(object, filename, datasetname, ...)
	
	require(reticulate)
	cat(crayon::green("Converting .h5 file to .txt"), sep = "\n")
	test <- reticulate::py_module_available('kelvinnypy')
	if(!test) {
		stop("python can't seem to find the module in this session of R")
	} else {
		kp <- reticulate::import("kelvinnypy")
		kp$h5totxt(filename)
	}
}