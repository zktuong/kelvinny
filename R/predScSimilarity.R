#' Predicts scRNA-seq data from trained multinomial cv.glmnet model
#' 
#' @param model trained model.
#' @param test SingleCellExperiment object for prediction
#' @param standardize a logical value specifying whether or not to standardize the test matrix
#' @param ... pass to predict. see ?predict.glmnet
#' @return Generates prediction.
#' @examples
#' pred <- predScSimilarity(model, test.sce)
#' @import glmnet
#' @import dplyr
#' @export
#'      
            
predScSimilarity <- function(model, test, standardize = TRUE, ...) {

    set.seed(42)
    preds <- list()
    trained.class <- names(model)

    standardizing <- function(X) {
        X <- X - mean(X)
        X <- X/sd(X)
        return(X)
    }

    newx <- as.matrix(t(assay(nathan.sce)))

    if (standardize == TRUE) {
        print("standardizing test data")
        newx <- t(apply(newx, 1, standardizing))
        }
    
    for(class in trained.class){
    message(sprintf("Predicting probabilities for %s", class))
    model.genes <- match(rownames(model[[class]]$glmnet.fit$beta), colnames(newx))

    preds[[class]] = predict(model[[class]], newx = newx[,model.genes], s = model[[class]]$lambda.1se, ...)
    colnames(preds[[class]]) <- class
    }
    return(preds)
}

    
