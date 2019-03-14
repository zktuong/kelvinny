#' Random Forest Classifcation for single-cell expression data or just general matrix
#' 
#' @param training training expression matrix.
#' @param test test expression matrix.
#' @param training.genes training genes to build classifier.
#' @param training.classes training classes to build classifier.
#' @param verbose verbosity.
#' @param probability Default = FALSE, which returns class when predicting. Otherwise TRUE returns probability estimates.
#' @param ... pass to ranger::ranger.
#' @return Training random forest classification model.
#' @examples
#' classifier <- RFclassifier(park.seurat, training.classes = park.seurat@ident, importance = "impurity")
#' prediction <- RFpredictor(classifier, test@data)
#' @import Seurat
#' @import ranger
#' @export
RFclassifier <- function(training, test, training.genes = NULL, training.classes = NULL, verbose = TRUE, probability = FALSE, ...) {
    SetIfNull <- function(x, default) {
  if(is.null(x = x)){
    return(default)
  } else {
    return(x)
  }
}
if(class(training) == "seurat"){
    training.classes <- as.vector(x = training.classes)
    training.genes <- SetIfNull(x = training.genes, default = rownames(x = training@data))
    training.data <- as.data.frame(x = as.matrix(x = t(x = training@data[training.genes, ])))
    training.data$class <- factor(x = training.classes)
    if (verbose) {
        message("Training Classifier ...")
    }
    classifier <- ranger::ranger(data = training.data, dependent.variable.name = "class", 
        classification = TRUE, write.forest = TRUE, probability = probability, ...)
    return(classifier)
} else {
	training.classes <- as.vector(x = training.classes)
    training.genes <- SetIfNull(x = training.genes, default = rownames(x = training))
    training.data <- as.data.frame(x = as.matrix(x = t(x = training[training.genes, ])))
    training.data$class <- factor(x = training.classes)
    if (verbose) {
        message("Training Classifier ...")
    }
    classifier <- ranger::ranger(data = training.data, dependent.variable.name = "class", 
        classification = TRUE, write.forest = TRUE, probability = probability, ...)
    return(classifier)
}
}