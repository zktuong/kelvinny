#' Random Forest Classifcation for single-cell expression data or just general matrix
#' 
#' @param classifier trained RF classifier.
#' @param test test expression matrix.
#' @param ... pass to predict.ranger.
#' @return Predicting using random forest classification.
#' @examples
#' classifier <- RFclassifier(seurat, training.classes = seurat@ident, importance = "impurity")
#' prediction <- RFpredictor(classifier, test@data)
#' @import ranger
#' @export
RFpredictor <- function(classifier, test, ...) {
if(class(test) == "seurat"){
    testData <- test
    features <- classifier$forest$independent.variable.names
    features.to.add <- tryCatch(
        setdiff(x = features, y = rownames(x = testData@data)), error = function(e) {
        setdiff(x = features, y = rownames(x = Seurat::GetAssayData(testData)))    
        }
        )
    data.to.add <- tryCatch(
        matrix(data = 0, nrow = length(x = features.to.add), ncol = ncol(x = testData@data)), error = function(e) {
        matrix(data = 0, nrow = length(x = features.to.add), ncol = ncol(x = Seurat::GetAssayData(testData)))    
        }
        )
    rownames(x = data.to.add) <- features.to.add
    testData.data <- tryCatch(
        rbind(testData@data, data.to.add), error = function(e) {
        rbind(Seurat::GetAssayData(testData), data.to.add)
        }
        )
    testData.data <- as.matrix(testData.data[features, ])
    message("Running Classifier ...")
    prediction <- predict(classifier, testData.data, ...)
    return(prediction)
} else {
    testData <- test
    features <- classifier$forest$independent.variable.names
    features.to.add <- setdiff(x = features, y = rownames(x = testData))
    data.to.add <- matrix(data = 0, nrow = length(x = features.to.add), ncol = ncol(x = testData))
    rownames(x = data.to.add) <- features.to.add
    testData <- rbind(testData, data.to.add)
    testData <- testData[features, ]
    testData <- as.matrix(x = t(testData))
    message("Running Classifier ...")
    prediction <- predict(classifier, testData, ...)
    return(prediction)
}
}