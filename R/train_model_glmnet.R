#' Training a glmnet model for binomial prediction.
#' 
#' @param data expression data with genes in rows and samples in columns.
#' @param variable_colname Column name for binary class.
#' @param alpha Default = 1 for Lasso penalty. 0 = Ridge penalty. between 0 - 1 = Elastic net. See ?cv.glmnet.
#' @param type.measure See ?cv.glmnet.
#' @param standardize Default = FALSE. Expects data to be pre-scaled. if not, glmnet will do it for you.
#' @param type term for internal prediction type. 
#' @param cutOff value for internal prediction.
#' @param nfolds number of folds for cross validation. use ncol(data) for LOOCV.
#' @param ... pass to glmnet
#' @return Generates a glmnet model.
#' @examples
#' model <- train_model_glmnet(data, variable_colname = "disease", alpha = 0.5, cutOff = 0.5, nfolds = ncol(data)) # LOOCV Elastic net regularization
#' #           healthy disease
#   disease       1       9
#   healthy      17       0
# [1] "Prediction accuracy: 0.962962962962963"
# [1] "Genes used in the model"
#             Coefficients OddsRatios
# (Intercept)     -1.09922    0.33313
# IL1RN            0.13296    1.14220
# INHBA            0.11830    1.12558
# EGR2             0.06974    1.07223
# TGM2             0.35543    1.42680
# SPP1             0.18706    1.20570
# PLEK             0.19937    1.22064
# DUSP7            0.11981    1.12729
# SOCS3            0.37537    1.45553
# GDF15           -0.16623    0.84685
# CPEB4            0.14038    1.15071
# BCL2L1           0.26513    1.30359
# PDPN             0.18267    1.20042
# PTGS2            0.04174    1.04263
# IL1A             0.07565    1.07859
# EDN1            -0.21577    0.80592
# SLAMF1           0.17746    1.19419
# GALNT7           0.01012    1.01017
# GBP2             0.01272    1.01280
# DNAJC19         -0.16747    0.84580
# GSDMD            0.46654    1.59447
# MYO1F            0.28252    1.32647
# ARHGEF10L       -0.14357    0.86626
# APOL3            0.07749    1.08057
# BMF             -0.17314    0.84102
# TLR8             0.15920    1.17257
#' @import glmnet
#' @import doMC
#' @import parallel
#' @export
train_model_glmnet <- function(data, variable_colname , alpha = 1, type.measure = "mse", type = "response", standardize = FALSE, cutOff = 0.5, nfolds = 10, ...) {
    set.seed(42)
    require(doMC)
    registerDoMC(cores = parallel::detectCores())
    # get column index of predicted variable in dataset
    variableColNum.ori <- grep(variable_colname, names(data))    
    # train model
        fit <- glmnet::glmnet(as.matrix(data[, -variableColNum.ori]), data[, variableColNum.ori], alpha = alpha, family = "binomial", standardize = standardize, ...)
        cv.fit <- glmnet::cv.glmnet(as.matrix(data[, -variableColNum.ori]), data[, variableColNum.ori], alpha = alpha, family = "binomial", 
            type.measure = type.measure, standardize = standardize, nfolds = nfolds, parallel = TRUE, ...)
        # plot result
        plot(cv.fit)
        # min value of lambda lambda_min <- cv.fit$lambda.min best value of lambda
        lambda_1se <- cv.fit$lambda.1se
        # regression coefficients
        print(coef(fit, s = lambda_1se))
        cv.fit_out <- as.matrix(coef(fit, s = lambda_1se))
        cv.fit_out <- as.data.frame(cv.fit_out)
        cv.fit_out$name <- row.names(cv.fit_out)
        sub_cv.fit_out <- cv.fit_out[cv.fit_out$`1` != 0, ]
        # which classes do these probabilities refer to? What are 1 and 0?
        print(paste("Groups in data..."))
        print(contrasts(data[, variableColNum.ori]))
        print("~~~~~~~~~")
        # get test data
        newx <- as.matrix(data[, -variableColNum.ori])
        # predict class, type=”class”
        pred <- predict(fit, newx = newx, s = lambda_1se, type = type, ...)
        # translate probabilities to predictions
        predict <- rep(levels(data[, variableColNum.ori])[1], nrow(data))
        predict[pred > cutOff] <- levels(data[, variableColNum.ori])[2]        
        if (alpha == 1) {
            print("Results from LASSO classification...")
        } else if (alpha == 0.5) {
            print("Results from ElasticNet classification...")
        } else if (alpha == 0) {
            print("Results from Ridge classification...")
        }
        # confusion matrix
        predict[is.na(predict)] <- "unassigned"
        prediction <- data.frame(pred = predict, true = data[, variableColNum.ori])
        print(table(prediction$pred, prediction$true))
        print(paste0("Prediction accuracy: ", mean(prediction$pred == prediction$true, na.rm = FALSE)))
        print("Genes used in the model")
        # extract the coefficients
        Coefficients <- as.matrix(coef(fit, s = lambda_1se))
        # calculate the odds ratios
        OddsRatios <- as.matrix(exp(Coefficients))
        # store coefficients and odds ratios as a new data frame
        df <- data.frame(Coefficients = Coefficients, OddsRatio = OddsRatios)
        df <- round(df, 5)
        colnames(df) <- c("Coefficients", "OddsRatios")
        df <- df[cv.fit_out$`1` != 0, ]
        print(df)
        return(cv.fit)
}