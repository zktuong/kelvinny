#' Trains scRNA-seq data via tunable logistic regression model
#' 
#' @param train_data Seurat object, SummarizedExperiment object or expression matrix for training
#' @param train_cell_type The cell types/clusters in the training data set
#' @param test_data Seurat object, SummarizedExperiment object or expression matrix for testing later
#' @param train_genes Genes to use for training. If not provided, it will try to pick from all genes in the training dataset as per default glmnet.
#' @param standardize a logical value specifying whether or not to standardize the train matrix
#' @param nfolds integer specifying bin for cross validation. Use all samples if doing LOOCV.
#' @param a tunable regularization parameter. 0 = ridge (L2), 1 = LASSO (L1), in between = Elastic-net
#' @param l.min logical. Choose between lambda.min or lambda.1se
#' @param multinomial logical. Choose between family = 'binomial' or 'multinomial'.
#' @param nParallel integer specifying number of cores for parallelization.
#' @param ... other functions pass to glmnet
#' @return Generates a trained model for predicting cell types for scRNAseq data
#' @examples
#' fit <- trainScSimilarity(trainData, clusters, testData)
#' @import glmnet
#' @import doMC
#' @import SummarizedExperiment
#' @export
#'      
            
trainScSimilarity <- function(train_data, train_cell_type, test_data, train_genes = NULL, 
    standardize = TRUE, nfolds = 10, a = 0.9, l.min = FALSE, multinomial = TRUE, 
    nParallel = parallel::detectCores(), ...) {
    
    fit <- list()
    
    standardizeSparse <- function(A) {
        A@x <- A@x/rep.int(colSums(A), diff(A@p))
        return(A)
    }
    
    getPopulationOffset = function(y) {
        if (!is.factor(y)) 
            y = factor(y)
        if (length(levels(y)) != 2) 
            stop("y must be a two-level factor")
        off = sum(y == levels(y)[2])/length(y)
        off = log(off/(1 - off))
        return(rep(off, length(y)))
    }
    
    if (nParallel > 1) {
        doMC::registerDoMC(cores = nParallel)
    }
    
    if (is.null(train_genes)) {
        if (class(train_data) == "SummarizedExperiment") {
            train_dat <- SummarizedExperiment::assay(train_data)
        } else if (class(train_data) == "Seurat") {
            train_dat <- tryCatch(train_data@data, error = function(e) {
                tryCatch(GetAssayData(object = train_data), error = function(e) {
                  warning(paste0("are you sure this is a seurat v3 object?"))
                  return(NULL)
                })
            })
        } else {
            train_dat <- train_data
        }

        if (class(test_data) == "SummarizedExperiment") {
            test_dat <- SummarizedExperiment::assay(test_data)
        } else if (class(test_data) == "Seurat") {
            test_dat <- tryCatch(test_data@data, error = function(e) {
                tryCatch(GetAssayData(object = test_data), error = function(e) {
                  warning(paste0("are you sure this is a seurat v3 object?"))
                  return(NULL)
                })
            })
        } else {
            test_dat <- test_data
        }
        cat(paste0("No pre-defined genes provided. Filtering ", crayon::red(dim(train_dat)[1]), " genes for training"), sep = "\n")
        
        Zero_col <- which(colSums(train_dat) == 0)
        duplicated_col <- which(duplicated(colnames(train_dat)) == TRUE)
        if (length(c(Zero_col, duplicated_col)) != 0) {
            cat(paste0("Removing ", crayon::red(length(c(Zero_col, duplicated_col))), " genes with no variance"), sep = "\n")
            train_dat <- train_dat[, -c(Zero_col, duplicated_col)]
        }
        
        genes.intersect <- intersect(row.names(test_dat), row.names(train_dat))
        train_dat <- train_dat[which(row.names(train_dat) %in% genes.intersect), ]

        cat(paste0("Submitting ", crayon::red(length(genes.intersect)), " intersecting genes to glmnet for selecting predictors"), sep = "\n")

        cat("Transposing matrix", sep = "\n")
        if (class(train_dat) == "matrix") {
            train_dat <- Matrix::Matrix(train_dat, sparse = TRUE)
        }
        train_dat <- t(train_dat)
        
        
        if (standardize == TRUE) {
            cat("Standardizing training dataset", sep = "\n")
            train_dat <- standardizeSparse(train_dat)
        }
        
        labels = levels(train_cell_type)
        if (length(labels) == 0) {
            train_cell_type <- factor(train_cell_type)
            labels <- levels(train_cell_type)
        }
        
        if (multinomial == FALSE) {
            cat(paste0(crayon::magenta("Training model with family = "), crayon::yellow("binomial")), sep = "\n")
            for (label in labels) {
                cat(crayon::green(paste0("Training model for ", crayon::red(label))), sep = "\n")
                celltype = factor(train_cell_type == label)
                
                fit[[label]] = tryCatch(glmnet::cv.glmnet(train_dat, celltype, 
                  offset = getPopulationOffset(celltype), family = "binomial", 
                  alpha = a, nfolds = nfolds, type.measure = "class", parallel = nParallel, 
                  ...), error = function(e) {
                  tryCatch(glmnet::cv.glmnet(train_dat, celltype, offset = getPopulationOffset(celltype), 
                    family = "binomial", alpha = a, nfolds = nfolds, type.measure = "class", 
                    parallel = nParallel, lambda = exp(seq(-10, -3, length.out = 100)), 
                    ...), error = function(e) {
                    warning(paste0("Could not train model for variable ", 
                      label))
                    return(NULL)
                  })
                })
            }
            cat("Extracting best gene features...", sep = "\n")
            for (i in 1:length(fit)) {
                if (l.min) {
                  cat("Choosing the best model but with the caveat that may be too complex, may be slightly overfitted", sep = "\n")
                  fit_out <- as.matrix(coef(fit[[i]], s = fit[[i]]$lambda.min))
                } else {
                  cat("Choosing the simplest model that has comparable error to the best model given the uncertainty", sep = "\n")
                  fit_out <- as.matrix(coef(fit[[i]], s = fit[[i]]$lambda.1se))
                }
                fit_out <- as.data.frame(fit_out)
                fit_out <- fit_out[fit_out[, 1] != 0, , drop = FALSE]
                cat(sep = "\n")
                cat(paste0("Best genes for ", names(fit)[i]), sep = "\n")
                print(fit_out, sep = "\n")
            }
            return(fit)
        } else {
            cat(paste0(crayon::magenta("Training model with family = "), crayon::yellow("multinomial")), sep = "\n")
            fit <- tryCatch(glmnet::cv.glmnet(train_dat, train_cell_type, 
                family = "multinomial", alpha = a, nfolds = nfolds, type.measure = "class", 
                parallel = nParallel, ...), error = function(e) {
                tryCatch(glmnet::cv.glmnet(train_dat, train_cell_type, 
                  family = "multinomial", alpha = a, nfolds = nfolds, type.measure = "class", 
                  parallel = nParallel, lambda = exp(seq(-10, -3, length.out = 100)), 
                  ...), error = function(e) {
                  warning(paste0("Could not train model with family = multinomial. Try with multinomial = FALSE"))
                  return(NULL)
                })
            })
            cat("Extracting best gene features", sep = "\n")
            for (i in 1:length(fit)) {
                if (l.min) {
                  fit_out <- as.matrix(coef(fit, s = fit$lambda.min)[[i]])
                } else {
                  fit_out <- as.matrix(coef(fit, s = fit$lambda.1se)[[i]])
                }
                fit_out <- as.data.frame(fit_out)
                fit_out <- fit_out[fit_out[, 1] != 0, , drop = FALSE]
                cat(sep = "\n")
                cat(paste0("Best genes for ", levels(train_cell_type)[i]), sep = "\n")
                print(fit_out, sep = "\n")
            }
            return(fit)
        }
    } else {
        if (class(train_data) == "SummarizedExperiment") {
            all_genes <- elementMetadata(train_data)[, 1]
            train_dat <- SummarizedExperiment::assay(train_data[which(all_genes %in% 
                train_genes)])
        } else if (class(train_data) == "Seurat") {
            all_genes <- tryCatch(all_genes <- rownames(train_data@data), 
                error = function(e) {
                  all_genes <- tryCatch(all_genes <- rownames(GetAssayData(object = train_data)), 
                    error = function(e) {
                      warning(paste0("are you sure this is a seurat v3 object?"))
                      return(NULL)
                    })
                })
            train_dat <- tryCatch(train_dat <- train_data@data[which(all_genes %in% 
                train_genes), ], error = function(e) {
                train_dat <- tryCatch(train_dat <- GetAssayData(object = train_data)[which(all_genes %in% 
                  train_genes), ], error = function(e) {
                  warning(paste0("are you sure this is a seurat v3 object?"))
                  return(NULL)
                })
            })
        } else {
            all_genes <- rownames(train_data)
            train_dat <- train_data[which(all_genes %in% train_genes), ]
        }
        cat(paste0("using ", dim(train_dat)[1], " genes for training model"), sep = "\n")
        
        cat("Transposing matrix", sep = "\n")
        if (class(train_dat) == "matrix") {
            train_dat <- Matrix::Matrix(train_dat, sparse = TRUE)
        }
        train_dat <- t(train_dat)
        
        if (standardize == TRUE) {
            cat("Standardizing training dataset", sep = "\n")
            train_dat <- standardizeSparse(train_dat)
        }
        
        labels = levels(train_cell_type)
        if (length(labels) == 0) {
            train_cell_type <- factor(train_cell_type)
            labels <- levels(train_cell_type)
        }
        if (multinomial == FALSE) {
            cat(paste0(crayon::magenta("Training model with family = "), crayon::yellow("binomial")), sep = "\n")
            for (label in labels) {
                cat(crayon::green(paste0("Training model for ", crayon::red(label))), sep = "\n")
                celltype = factor(train_cell_type == label)
                fit[[label]] = tryCatch(glmnet::cv.glmnet(train_dat, celltype, 
                  family = "binomial", alpha = a, nfolds = nfolds, type.measure = "class", 
                  parallel = nParallel, ...), error = function(e) {
                  tryCatch(glmnet::cv.glmnet(train_dat, celltype, family = "binomial", 
                    alpha = a, nfolds = nfolds, type.measure = "class", 
                    parallel = nParallel, lambda = exp(seq(-10, -3, length.out = 100)), 
                    ...), error = function(e) {
                    warning(paste0("Could not train model for variable ", 
                      label))
                    return(NULL)
                  })
                })
            }
            cat("Extracting best gene features", sep = "\n")
            for (i in 1:length(fit)) {
                if (l.min) {
                  cat("Choosing the best model but with the caveat that may be too complex, may be slightly overfitted", sep = "\n")
                  fit_out <- as.matrix(coef(fit[[i]], s = fit[[i]]$lambda.min))
                } else {
                  cat("Choosing the simplest model that has comparable error to the best model given the uncertainty", sep = "\n")
                  fit_out <- as.matrix(coef(fit[[i]], s = fit[[i]]$lambda.1se))
                }
                fit_out <- as.data.frame(fit_out)
                fit_out <- fit_out[fit_out[, 1] != 0, , drop = FALSE]
                cat(sep = "\n")
                cat(paste0("Best genes for ", names(fit)[i]), sep = "\n")
                print(fit_out, sep = "\n")
            }
            return(fit)
        } else {
            cat(paste0(crayon::magenta("Training model with family = "), crayon::yellow("multinomial")), sep = "\n")
            fit <- tryCatch(glmnet::cv.glmnet(train_dat, train_cell_type, 
                family = "multinomial", alpha = a, nfolds = nfolds, type.measure = "class", 
                parallel = nParallel, ...), error = function(e) {
                tryCatch(glmnet::cv.glmnet(train_dat, train_cell_type, 
                  family = "multinomial", alpha = a, nfolds = nfolds, type.measure = "class", 
                  parallel = nParallel, lambda = exp(seq(-10, -3, length.out = 100)), 
                  ...), error = function(e) {
                  warning(paste0("Could not train model with family = multinomial. Try with multinomial = FALSE"))
                  return(NULL)
                })
            })
            cat("Extracting best gene features", sep = "\n")
            for (i in 1:length(fit)) {
                if (l.min) {
                  fit_out <- as.matrix(coef(fit, s = fit$lambda.min)[[i]])
                } else {
                  fit_out <- as.matrix(coef(fit, s = fit$lambda.1se)[[i]])
                }
                fit_out <- as.data.frame(fit_out)
                fit_out <- fit_out[fit_out[, 1] != 0, , drop = FALSE]
                cat(sep = "\n")
                cat(paste0("Best genes for ", levels(train_cell_type)[i]), sep = "\n")
                print(fit_out, sep = "\n")
            }
            return(fit)
        }
    }
}
