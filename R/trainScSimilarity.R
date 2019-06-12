#' Trains scRNA-seq data to be used as binomial glmnet::cv.glmnet model
#' 
#' @param train_data seurat object, SummarizedExperiment object or expression matrix for training
#' @param train_cell_type The cell types/clusters in the training data set
#' @param train_genes Genes to use for training. If not provided, it will try to pick from all genes in the training dataset as per default glmnet.
#' @param standardize a logical value specifying whether or not to standardize the train matrix
#' @param nfolds integer specifying bin for cross validation. Use all samples if doing LOOCV.
#' @param alpha regularization parameter
#' @param nParallel integer specifying number of cores for parallelization.
#' @param ... pass to glmnet
#' @return Generates a trained model for predicting cell types for scRNAseq data
#' @examples
#' fit <- trainScSimilarity(data, cluster, genes)
#' @import glmnet
#' @import doMC
#' @export
#'      
            
trainScSimilarity <- function(train_data, train_cell_type, train_genes = NULL, standardize = TRUE, nfolds = 10, alpha = 0.99, nParallel = parallel::detectCores(),  ...) 
{
    fit <-  list()

    standardizing <- function(X) {
        X <- X - mean(X)
        X <- X/sd(X)
        return(X)
    }

    if(nParallel > 1) {
    doMC::registerDoMC(cores = nParallel)
    }
    
    if(is.null(train_genes)) {
        if(class(train_data) == "SummarizedExperiment"){
            train_dat <- SummarizedExperiment::assay(train_data)    
        } else if (class(train_data) == "Seurat"){
        train_dat <- tryCatch(
                train_data@data, error = function(e) {
                tryCatch(
                        GetAssayData(object = train_data), error = function(e) {
                warning(paste0("are you sure this is a seurat v3 object?"))
                return(NULL)
                })
            })
        } else {
            train_dat <- train_data
        }
        print(paste0("No pre-defined genes provided. Submitting ", dim(train_dat)[1], " genes to glmnet for selecting predictors"))
        print("converting to matrix")
        print(class(train_dat))
        train_dat <- as.matrix(train_dat)

        print("transposing matrix")
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
            print(paste0("Training model for ", label))
            celltype = factor(train_cell_type == label)
            getPopulationOffset = function(y) {
                if(!is.factor(y))
                    y = factor(y)
                if(length(levels(y)) != 2)
                    stop("y must be a two-level factor")
                    off = sum(y == levels(y)[2])/length(y)
                    off = log(off/(1-off))
            return(rep(off,length(y)))
        }
            fit[[label]] = tryCatch(
                glmnet::cv.glmnet(train_dat, celltype, offset = getPopulationOffset(celltype), family = 'binomial', alpha = alpha, nfolds = nfolds, type.measure = 'class', parallel = nParallel, ...), 
                error = function(e) {
                    tryCatch(
                        glmnet::cv.glmnet(train_dat, celltype, offset = getPopulationOffset(celltype), family = 'binomial', alpha = alpha, nfolds = nfolds, type.measure = 'class', parallel = nParallel, lambda = exp(seq(-10, -3, length.out=100)), ...),
                        error = function(e) {
                            warning(paste0("Could not train model for variable ", label))
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
        if(class(train_data) == "SummarizedExperiment"){
            all_genes <- elementMetadata(train_data)[, 1]
            train_dat <- SummarizedExperiment::assay(train_data[which(all_genes %in% train_genes)])
        } else if (class(train_data) == "Seurat"){
            all_genes <- tryCatch(
                all_genes <- rownames(train_data@data), error = function(e){
                all_genes <- tryCatch(all_genes <- rownames(GetAssayData(object = train_data)), error = function(e){
                warning(paste0("are you sure this is a seurat v3 object?"))
                return(NULL)    
                })})
            train_dat <- tryCatch(
            train_dat <- train_data@data[which(all_genes %in% train_genes), ], error = function(e){
            train_dat <- tryCatch(train_dat <- GetAssayData(object = train_data)[which(all_genes %in% train_genes), ], error = function(e){
            warning(paste0("are you sure this is a seurat v3 object?"))
            return(NULL)                    
            })})
        } else {
            all_genes <- rownames(train_data)
            train_dat <- train_data[which(all_genes %in% train_genes), ]
        }
        print(paste0("using ", dim(train_dat)[1], " genes for training model"))
        
        train_dat <- as.matrix(train_dat)
        train_dat <- t(train_dat)
        
        if (standardize == TRUE) {
            print("standardizing training data")
            train_dat <- t(apply(train_dat, 1, standardizing))
        }

        labels = names(sort(table(as.character(train_cell_type))))
        for(label in labels){
            print(paste0("Training model for ", label))
            celltype = factor(train_cell_type == label)
            fit[[label]] = tryCatch(
                glmnet::cv.glmnet(train_dat, celltype, family = 'binomial', alpha = alpha, nfolds = nfolds, type.measure = 'class', parallel = nParallel, ...), 
                error = function(e) {
                    tryCatch(
                        glmnet::cv.glmnet(train_dat, celltype,family = 'binomial', alpha = alpha, nfolds = nfolds, type.measure = 'class', parallel = nParallel, lambda = exp(seq(-10, -3, length.out=100)), ...),
                        error = function(e) {
                            warning(paste0("Could not train model for variable ", label))
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