setwd(system("echo $PROJ_DIR", intern = TRUE))
source("src/compute/utils.R")
setwd(system("echo $ZHU_DIR", intern = TRUE))
source("R/opt.R")
source("R/gen_data.R")
library(genlasso)
library(igraph)
library(splash)
# library(FGSG)
library(BigVAR)
library(glmnet)
# library(pracma)
library(Matrix)
library(caret)
setwd(system("echo $PROJ_DIR", intern = TRUE))

fit_splash.cv <- function(y, alpha, nlambdas, n_folds = 5, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]

    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    # Create cross-validation folds
    folds <- create_folds(y_train, nfolds = n_folds)

    # Fit the model on each fold and save the results
    errors <- matrix(NA, nrow = length(folds), ncol = nlambdas)
    for (i in 1:nfolds) {
        # Split the data into training and validation sets
        y_train_cv <- y_train[, folds$train[[i]]]
        y_val_cv <- y_train[, folds$test[[i]]]

        # Fit the model on the training set
        model_cv <- splash(t(y_train_cv), alphas = c(alpha), n_lambdas = nlambdas, ...)

        # Compute the prediction error on the validation set
        C_cv <- array(NA, dim = c(p, p, nlambda))
        A_cv <- array(NA, dim = c(p, p, nlambda))
        B_cv <- array(NA, dim = c(p, p, nlambda))
        y_pred_cv <- array(NA, dim = c(p, dim(y_val_cv)[2], nlambda))
        errors_cv <- array(NA, dim = c(nfolds, nlambdas))
        for (j in 1:nlambdas) {
            A_cv[, , j] <- model_cv$AB[, 1:p, j]
            B_cv[, , j] <- model_cv$AB[, (p + 1):(2 * p), j]
            C_cv[, , j] <- AB_to_C(A, B)
            y_pred_cv[, , j] <- predict_with_C(C[, , j], y_train_cv, y_val_cv)
            errors_cv[i, j] <- calc_msfe(y_pred_cv[, , j], y_val_cv)
        }
    }

    # Get the best lambda value
    best_lambda <- model$lambda[which.min(colMeans(errors_cv)), 1]

    # Fit the model on the entire training set
    t0 <- Sys.time()
    model <- splash(t(y_train), alphas = c(alpha), lambda = c(best_lambda), ...)
    t1 <- Sys.time()

    # Extract the model output
    A <- model$AB[, 1:p, ]
    B <- model$AB[, (p + 1):(2 * p), ]
    C <- AB_to_C(A, B)

    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = C,
        y_pred = predict_with_C(C, y_train, y_test),
        errors_cv = errors_cv,
        best_lambda = best_lambda
    )
    return(output_list)
}


create_folds <- function(y, nfolds = 5, test_size = 0.2) {
    # Create cross-validation folds
    time_indices <- 1:dim(y)[2]
    initialWindow <- dim(y)[2] - floor(dim(y)[2] * test_size)
    horizon <- dim(y)[2] * 0.2 / nfolds
    folds <- createTimeSlices(time_indices,
        initialWindow = InitialWindow,
        horizon = horizon,
        fixedWindow = TRUE,
        skip = horizon - 1
    )
    return(folds)
}
