#' Bootstraps training and prediction of scRNA-seq data using binomial/multinomial logistic regression with tunable penalties
#' 
#' @param trainingSet Seurat object, SummarizedExperiment object or expression matrix for training
#' @param trainingCellType The cell types/clusters in the training data set
#' @param testingSet Seurat object, SummarizedExperiment object or expression matrix for testing
#' @param nboots integet specifying number of bootstrap runs to run
#' @param trainingGenes Genes to use for training. If not provided, it will try to pick from all genes in the training dataset.
#' @param scale a logical value specifying whether or not to standardize the train matrix
#' @param nfold.cv integer specifying bin for cross validation. Use all samples if doing LOOCV.
#' @param a.parameter tunable regularization parameter. 0 = ridge (L2), 1 = LASSO (L1), in between = Elastic-net
#' @param lambda.min logical. Choose between lambda.min or lambda.1se
#' @param family.multinomial logical. Choose between family = 'binomial' or 'multinomial'.
#' @param nCores integer specifying number of cores for parallelization.
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

similarity_bootstrap <- function(trainingSet, trainingCellType, testingSet, nboots = 50, trainingGenes = NULL,
    scale = TRUE, nfold.cv = 10, a.parameter = 0.9, lambda.min = FALSE, family.multinomial = TRUE, 
    nCores = parallel::detectCores(), ...){
    require(pbmcapply)
    boots <- as.list(1:nboots)
    prediction <- mclapply(boots, function(x){
        fit <- trainScSimilarity(train_data = trainingSet, train_cell_type = trainingCellType, test_data = testingSet, train_genes = trainingGenes, standardize = scale, nfolds = nfold.cv, a = a.parameter, l.min = lambda.min, multinomial = family.multinomial, nParallel = nCores, ...)
        pred <- predScSimilarity(model = fit, test = testingSet, standardize = scale, l.min = lambda.min, ...)
        return(pred)
    }, mc.cores = nCores)
    return(prediction)
}
