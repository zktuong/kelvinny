#' Predicting new data using trained glmnet model
#' 
#' @param model trained model.
#' @param new_data new expression data.
#' @param type Type of prediction needed. see ?predict.glmnet
#' @param ... pass to predict. see ?predict.glmnet
#' @return Generates prediction.
#' @examples
#' # pred <- test_model_glmnet(model = model, new_data = newdat, type = "link")
#' @import glmnet
#' @export
test_model_glmnet <- function(model, new_data, type = "response", ...) {
    set.seed(42)
    # extract some of the features
    model.lambda_1se <- model$lambda.1se
    model.out <- as.data.frame(as.matrix(coef(model, model.lambda_1se)))
    model.out$name <- row.names(model.out)
    model.genes <- model.out$name
    model.genes <- model.genes[-1]

    # convert new data to matrix
    # match genes to model genes
    # subset new data according to model genes
    m <- match(model.genes, row.names(new_data))
    if(any(is.na(m))){
        unmatched <- which(is.na(m))
        subset_newdat <- new_data[m, ]
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print(paste('Following genes not found'))
        print(model.genes[unmatched])
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("Substituting the gene expression values = 0")
        row.names(subset_newdat) <- model.genes
        newx <- t(subset_newdat)
        s.newx  <- scale(newx)
        s.newx[,unmatched] <- rep(0)
    } else {
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("All genes matched! predicting now!")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        subset_newdat <- new_data[m, ]
        newx <- t(subset_newdat)
        s.newx  <- scale(newx)
    }
    
    # do prediction
    predict_samples <- predict(model, newx = s.newx, type = type, s = model.lambda_1se, ...)
    
    # save results
    return(predict_samples)
}
