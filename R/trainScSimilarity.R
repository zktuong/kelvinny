#' Trains scRNA-seq data to be used as multinomial glmnet::cv.glmnet model
#' 
#' @param train_data SingleCellExperiment object for training
#' @param train_cell_type The cell types/clusters in the training data set
#' @param train_genes Genes to use for training. If not provided, it will try to pick from all genes in the training dataset as per default glmnet.
#' @param standardize a logical value specifying whether or not to standardize the train matrix
#' @param nfolds integer specifying bin for cross validation. Use all samples if doing LOOCV.
#' @param nParallel integer specifying number of cores for parallelization.
#' @param ... pass to glmnet
#' @return Generates a trained model for predicting cell types for scRNAseq data
#' @examples
#' fit <- trainScSimilarity(data, cluster, genes)
#' @import glmnet
#' @import dplyr
#' @import SummarizedExperiment
#' @export
#'      
            
trainScSimilarity <- function(train_data, train_cell_type, train_genes = NULL, standardize = TRUE, nfolds = 10, nParallel = parallel::detectCores(), ...) 
{
    set.seed(42)
    fit <-  list()

    standardizing <- function(X) {
        X <- X - mean(X)
        X <- X/sd(X)
        return(X)
    }

    if(nParallel > 1)
    doMC::registerDoMC(cores = nParallel)

    if(is.null(train_genes)) {
        train_dat <- SummarizedExperiment::assay(train_data)
        print(paste0("No pre-defined genes provided. Submitting ", dim(train_dat)[1], " genes to glmnet for selecting predictors"))
        train_dat <- t(train_dat)
        Zero_col <- which(colSums(train_dat) == 0)
        duplicated_col <- which(duplicated(colnames(train_dat)) == TRUE)
        if (length(c(Zero_col, duplicated_col)) != 0) {
            print(paste0("removing ", length(c(Zero_col, duplicated_col)), " genes with no variance"))
            train_dat <- train_dat[, -c(Zero_col, duplicated_col)]
        }
        if (standardize == TRUE) {
            print("standardizing prediction/target dataset")
            train_dat <- t(apply(train_dat, 1, standardizing))
        }  
        labels = names(sort(table(as.character(train_cell_type))))
        for(label in labels){
            message(sprintf("Training model for %s", label))
            celltype = factor(train_cell_type == label)
            fit[[label]] = tryCatch(
                glmnet::cv.glmnet(train_dat, celltype, family = 'binomial', alpha = 0.99, nfolds = nfolds, type.measure = 'class', parallel = nParallel, ...), 
                error = function(e) {
                    tryCatch(
                        glmnet::cv.glmnet(train_dat, celltype,family = 'binomial', alpha = 0.99, nfolds = nfolds, type.measure = 'class', parallel = nParallel, lambda = exp(seq(-10, -3, length.out=100)), ...),
                        error = function(e) {
                            warning(sprintf("Could not train model for variable %s", label))
                            return(NULL)
                        }
                    )
                }
            )
          }
        print("extracting best gene features...")
        for (i in 1:length(fit)) {
            fit_out <- as.matrix(coef(fit[[i]], s = fit[[i]]$lambda.1se))
            fit_out <- as.data.frame(fit_out)
            fit_out <- fit_out[fit_out[,1] != 0, , drop = FALSE]
            cat(sep ="\n")
            print(paste0("Best genes for ", names(fit)[i]))
            print(fit_out)
            
        }
        return(fit)
    } else {
        all_genes <- elementMetadata(train_data)[, 1]
        train_dat <- SummarizedExperiment::assay(train_data[which(all_genes %in% train_genes)])
        print(paste0("using ", dim(train_dat)[1], " genes for training model"))
        train_dat <- t(train_dat)
        
        if (standardize == TRUE) {
            print("standardizing training data")
            train_dat <- t(apply(train_dat, 1, standardizing))
        }

        labels = names(sort(table(as.character(train_cell_type))))
        for(label in labels){
            message(sprintf("Training model for %s", label))
            celltype = factor(train_cell_type == label)
            fit[[label]] = tryCatch(
                glmnet::cv.glmnet(train_dat, celltype, family = 'binomial', alpha = 0.9, nfolds = nfolds, type.measure = 'class', parallel = nParallel, ...), 
                error = function(e) {
                    tryCatch(
                        glmnet::cv.glmnet(train_dat, celltype,family = 'binomial', alpha = 0.9, nfolds = nfolds, type.measure = 'class', parallel = nParallel, lambda = exp(seq(-10, -3, length.out=100)), ...),
                        error = function(e) {
                            warning(sprintf("Could not train model for variable %s", label))
                            return(NULL)
                        }
                    )
                }
            )
          }
        print("extracting best gene features...")
        for (i in 1:length(fit)) {
            fit_out <- as.matrix(coef(fit[[i]], s = fit[[i]]$lambda.1se))
            fit_out <- as.data.frame(fit_out)
            fit_out <- fit_out[fit_out[,1] != 0, , drop = FALSE]
            cat(sep ="\n")
            print(paste0("Best genes for ", names(fit)[i]))
            print(fit_out)
            
        }
        return(fit)
    }
}