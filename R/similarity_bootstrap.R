#' Bootstraps training and prediction of scRNA-seq data using binomial/multinomial logistic regression with tunable penalties
#' 
#' @param trainingSet Seurat object, SummarizedExperiment object or expression matrix for training
#' @param trainingCellType The cell types/clusters in the training data set
#' @param testingSet Seurat object, SummarizedExperiment object or expression matrix for testing
#' @param nboots integet specifying number of bootstrap runs to run
#' @param simplify logical. collapse the predictions to a single data frame
#' @param nCores integer specifying number of cores for parallelization.
#' @param verbose logical. if TRUE, print the training/prediction steps. Otherwise, minimal messages will be returned on the screen
#' @param ... other functions pass to glmnet and predict. see ?glmnet and ?predict.glmnet
#' @return Generates prediction.
#' @examples
#' pred <- predScSimilarity(model, test.sce)
#' @import glmnet
#' @import parallel
#' @import SummarizedExperiment
#' @import Matrix
#' @import doMC
#' @export
#' 

similarity_bootstrap <- function(trainingSet, trainingCellType, testingSet, nboots = 50, simplify = FALSE, nCores = parallel::detectCores(), verbose = FALSE, ...){
    
    sdList <- function(pred.list) {
    n <- length(pred.list); 
    rc <- dim(pred.list[[1]]); 
    ar <- array(unlist(pred.list), c(rc, n), list(rownames(pred.list[[1]]), colnames(pred.list[[1]]))); 
    sdf <- apply(ar, c(1, 2), sd)
    return(sdf)}

    require(foreach)
    require(doMC)
    doMC::registerDoMC(cores = nCores)

    if(verbose){
        prediction <- foreach(i = 1:nboots) %dopar% {
            prediction <- foreach(i = 1:nboots) %dopar% {
            cat(paste0("Training bootstrap #", i))
            fit_model <- trainScSimilarity(train_data = trainingSet, train_cell_type = trainingCellType, test_data = testingSet, ...)
            cat(paste0("Predicting bootstrap #", i))
            pred <- predScSimilarity(model = fit_model, test = testingSet, ...)
            cat(green(paste0("Finished bootstrap #", i)))
            return(pred)
            }
        } else {
            print(paste0("Training bootstrap #", i))
            fit_model <- invisible(trainScSimilarity(train_data = trainingSet, train_cell_type = trainingCellType, test_data = testingSet, ...))
            cat(paste0("Predicting bootstrap #", i))
            pred <- invisible(predScSimilarity(model = fit_model, test = testingSet, ...))
            cat(paste0("Finished bootstrap #", i))
            return(pred)
            } 
        }
    
    if(length(prediction) > 1) {
        prediction_summary <- list()
        if(simplify) {
            cat("Averaging predictions")
            prediction_summary[["prediction"]] <- Reduce("+", prediction) / length(prediction)
            cat("Obtaining standard deviations")
            prediction_summary[["SD"]] <- sdList(prediction)
            return(prediction_summary)
        } else {
            cat("Returning all predictions")
            return(prediction)
        }
    } else {
        if(simplify) {
            warning("Cannot average predictions nor return SDs because only 1 run of training/prediction performed")
        }
        return(prediction)
    }
}
