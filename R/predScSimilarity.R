#' Predicts expression data for class probability/similarity with a trained binomial/multinomial glmnet logistic regression model
#' 
#' @param model trained model.
#' @param test Seurat object, SummarizedExperiment object or expression matrix for prediction
#' @param standardize a logical value specifying whether or not to standardize the test matrix
#' @param l.min logical. FALSE will use lambda.1se for prediction
#' @param ... pass to predict function see ?predict.glmnet
#' @return Generates prediction.
#' @examples
#' pred <- predScSimilarity(model, test.sce)
#' @import glmnet
#' @import SummarizedExperiment
#' @export
#'      
            
predScSimilarity <- function(model, test, standardize = TRUE, l.min = FALSE, ...) {

    preds <- list()
    model.genes <- list()    
    
    standardizeSparse <- function(A) {
        A@x <- A@x/rep.int(colSums(A), diff(A@p))
        return(A)
    }

    if(class(test) == "SummarizedExperiment"){
        newx <- SummarizedExperiment::assay(test)
    } else if (class(test) == "Seurat"){
        newx <- tryCatch(test@data, error = function(e) {
            tryCatch(GetAssayData(object = test), error = function(e) {
                warning(paste0("are you sure this is a seurat v3 object?"))
                return(NULL)
    		})
        })    
	} else {
        newx <- test
    }

	if (class(newx) == "matrix") {
		newx <- Matrix::Matrix(newx, sparse = TRUE)
	}
	cat("Transposing test matrix")
	newx <- t(newx)    

	if (standardize == TRUE) {
		cat("Standardizing training dataset")
		newx <- standardizeSparse(newx)
	}

	if(class(model) == "list"){
		trained.class <- names(model)
		trained.class <- factor(trained.class, levels = trained.class)
    } else {
    	trained.class <- names(model$glmnet.fit$beta)
    	trained.class <- factor(trained.class, levels = trained.class)
    }

	if(class(model) == "list"){
    	if (l.min == FALSE) {
        	for(class in trained.class){
        	cat(red(paste0("Predicting probabilities for ", class)))
        	model.genes[[class]] <- match(rownames(model[[class]]$glmnet.fit$beta), colnames(newx))
        	model.genes[[class]] <- model.genes[[class]][!is.na(model.genes[[class]])]
        	preds[[class]] <- predict(model[[class]], newx = newx[,model.genes[[class]]], s = model[[class]]$lambda.1se, newoffset = rep(0, nrow(newx)), ...)
        	colnames(preds[[class]]) <- class        
        	} 
    	} else {
        	for(class in trained.class){
        	cat(green(paste0("Predicting probabilities for ", class)))
        	model.genes[[class]] <- match(rownames(model[[class]]$glmnet.fit$beta), colnames(newx))
			model.genes[[class]] <- model.genes[[class]][!is.na(model.genes[[class]])]
        	preds[[class]] = predict(model[[class]], newx = newx[,model.genes[[class]]], s = model[[class]]$lambda.min, newoffset = rep(0, nrow(newx)), ...)
        	colnames(preds[[class]]) <- class           
        	}
        }
    } else {
    	cat(green("Predicting probabilities"))
    	if (l.min == FALSE) {    		
    		model.genes <- match(rownames(model$glmnet.fit$beta[[1]]), colnames(newx))
    		model.genes <- model.genes[!is.na(model.genes)]
    		preds <- predict(model, newx = newx[,model.genes], s = model$lambda.1se, ...)
    	} else {
    		model.genes <- match(rownames(model$glmnet.fit$beta[[1]]), colnames(newx))
    		model.genes <- model.genes[!is.na(model.genes)]
    		preds <- predict(model, newx = newx[,model.genes], s = model$lambda.min, ...)
    	}

	}
    preds <- as.data.frame(preds)
    colnames(preds) <- gsub("[.]1$", "", colnames(preds))
    return(preds)
}


    
