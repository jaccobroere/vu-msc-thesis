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

fit_splash.cv <- function(y, alpha, nlambdas = 20, nfolds = 5, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]

    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    # Create cross-validation folds
    folds <- create_folds(y_train, nfolds = nfolds)

    # Fit the model on each fold and save the results
    errors <- matrix(NA, nrow = nfolds, ncol = nlambdas)
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
            C_cv[, , j] <- AB_to_C(A_cv[, , j], B_cv[, , j])
            y_pred_cv[, , j] <- predict_with_C(C_cv[, , j], y_train_cv, y_val_cv)
            errors_cv[i, j] <- calc_msfe(y_pred_cv[, , j], y_val_cv)
        }
    }

    # Get the best lambda value
    best_lambda <- model_cv$lambda[which.min(colMeans(errors_cv)), 1]

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
        best_lambda = best_lambda,
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
    return(output_list)
}

fit_fsplash.cv <- function(y, bandwidths, graph, Dtilde_inv, alpha, nlambdas = 20, nfolds = 5, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]
    m <- ecount(graph)
    k <- vcount(graph)

    # Parse bandwidths
    h0 <- bandwidths[1]
    h1 <- bandwidths[2]

    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    # Create cross-validation folds
    folds <- create_folds(y_train, nfolds = nfolds)

    # Fit the model on each fold and save the results
    errors <- matrix(NA, nrow = nfolds, ncol = nlambdas)
    for (i in 1:nfolds) {
        # Split the data into training and validation sets
        y_train_cv <- y_train[, folds$train[[i]]]
        y_val_cv <- y_train[, folds$test[[i]]]

        # Construct Vhat_d and sigma_hat from the training validation
        Sigma0 <- calc_Sigma_j(y_train_cv, 0)
        Sigma1 <- calc_Sigma_j(y_train_cv, 1)
        Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
        sigma_hat <- vec_sigma_h(Sigma1, h1)
        Vhat_d <- construct_Vhat_d(Vhat)

        # Transform Generalized Lasso problem into LASSO problem
        XD1 <- Vhat_d %*% Dtilde_inv # Same as Vhat_d * inv(Dtilde)
        X1 <- XD1[, 1:m]
        X2 <- XD1[, (m + 1):dim(XD1)[2]]
        X2_plus <- solve((t(X2) %*% X2), t(X2)) # Same as inv(t(X2) %*% X2) %*% t(X2)

        # Transform the input to LASSO objective
        IminP <- diag(nrow(X2)) - X2 %*% X2_plus # Calculate I - P directly, P = X2 %*% X2_plus
        ytilde <- as.vector(IminP %*% t(sigma_hat))
        Xtilde <- IminP %*% X1

        # Fit the model on the training set
        model_cv <- glmnet(XD1, ytilde, nlambda = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4, ...)

        # Compute the prediction error on the validation set
        C_cv <- array(NA, dim = c(p, p, nlambda))
        A_cv <- array(NA, dim = c(p, p, nlambda))
        B_cv <- array(NA, dim = c(p, p, nlambda))
        y_pred_cv <- array(NA, dim = c(p, dim(y_val_cv)[2], nlambda))
        errors_cv <- array(NA, dim = c(nfolds, nlambdas))
        for (j in 1:nlambdas) {
            theta1 <- as.vector(model$beta[, j])
            theta2 <- as.vector(X2_plus %*% (t(sigma_hat) - X1 %*% theta1))
            coef_cv <- Dtilde_inv %*% c(theta1, theta2)
            AB <- coef_to_AB(coef_cv, p)
            A_cv[, , j] <- AB$A
            B_cv[, , j] <- AB$B
            C_cv[, , j] <- AB_to_C(AB$A, AB$B)
            y_pred_cv[, , j] <- predict_with_C(C_cv[, , j], y_train_cv, y_val_cv)
            errors_cv[i, j] <- calc_msfe(y_pred_cv[, , j], y_val_cv)
        }
    }

    # Get the best lambda value
    best_lambda <- model_cv$lambda[which.min(colMeans(errors_cv))]

    # Construct Vhat_d and sigma_hat from the training validation
    Sigma0 <- calc_Sigma_j(y_train, 0)
    Sigma1 <- calc_Sigma_j(y_train, 1)
    Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
    sigma_hat <- vec_sigma_h(Sigma1, h1)
    Vhat_d <- construct_Vhat_d(Vhat)

    t0 <- Sys.time()
    XD1 <- Vhat_d %*% Dtilde_inv # Same as Vhat_d * inv(Dtilde)
    X1 <- XD1[, 1:m]
    X2 <- XD1[, (m + 1):dim(XD1)[2]]
    X2_plus <- solve((t(X2) %*% X2), t(X2)) # Same as inv(t(X2) %*% X2) %*% t(X2)

    # Transform the input to LASSO objective
    IminP <- diag(nrow(X2)) - X2 %*% X2_plus # Calculate I - P directly, P = X2 %*% X2_plus
    ytilde <- as.vector(IminP %*% t(sigma_hat))
    Xtilde <- IminP %*% X1
    # Fit the model on the entire training set

    # Fit the model on the training set
    model_cv <- glmnet(XD1, ytilde, lambda = best_lambda, alpha = 1, intercept = FALSE, ...)
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
        best_lambda = best_lambda,
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
    return(output_list)
}

fit_ssfsplash.cv <- function(y, bandwidths, graph, Dtilde_SSF_inv, alpha, nlambdas = 20, nfolds = 5, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]
    m <- ecount(graph)
    k <- vcount(graph)
    h <- floor(p / 4) # Bandwidth for the A and B matrix

    # Parse bandwidths
    h0 <- bandwidths[1]
    h1 <- bandwidths[2]

    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    # Create cross-validation folds
    folds <- create_folds(y_train, nfolds = nfolds)

    # Fit the model on each fold and save the results
    errors <- matrix(NA, nrow = nfolds, ncol = nlambdas)
    for (i in 1:nfolds) {
        # Split the data into training and validation sets
        y_train_cv <- y_train[, folds$train[[i]]]
        y_val_cv <- y_train[, folds$test[[i]]]

        # Construct Vhat_d and sigma_hat from the training validation
        Sigma0 <- calc_Sigma_j(y_train_cv, 0)
        Sigma1 <- calc_Sigma_j(y_train_cv, 1)
        Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
        sigma_hat <- vec_sigma_h(Sigma1, h1)
        Vhat_d <- construct_Vhat_d(Vhat)

        # Transform Generalized Lasso problem into LASSO problem
        multipliers <- c(
            (p - h):(p - 1),
            rev((p - h):(p - 1)),
            (p - h):(p - 1),
            (p),
            rev((p - h):(p - 1))
        )
        M_inv <- .sparseDiagonal(x = c(
            rep((1 / alpha), m),
            1 / ((1 - alpha) * sqrt(multipliers))
        ))
        Dtilde_SSF_inv_gamma <- Dtilde_SSF_inv %*% M_inv

        # Transform the input to LASSO objective (see Tibshirani and Taylor, 2011)
        XD1 <- Vhat_d %*% Dtilde_SSF_inv_gamma # Same as Vhat_d * inv(Dtilde)

        # Fit the model on the training set
        model_cv <- glmnet(XD1, as.vector(sigma_hat), nlambda = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4, ...)

        # Compute the prediction error on the validation set
        C_cv <- array(NA, dim = c(p, p, nlambda))
        A_cv <- array(NA, dim = c(p, p, nlambda))
        B_cv <- array(NA, dim = c(p, p, nlambda))
        y_pred_cv <- array(NA, dim = c(p, dim(y_val_cv)[2], nlambda))
        errors_cv <- array(NA, dim = c(nfolds, nlambdas))
        for (j in 1:nlambdas) {
            theta <- as.vector(model$beta[, j])
            coef_cv <- Dtilde_SSF_inv_gamma %*% theta
            AB <- coef_to_AB(coef_cv, p)
            A_cv[, , j] <- AB$A
            B_cv[, , j] <- AB$B
            C_cv[, , j] <- AB_to_C(AB$A, AB$B)
            y_pred_cv[, , j] <- predict_with_C(C_cv[, , j], y_train_cv, y_val_cv)
            errors_cv[i, j] <- calc_msfe(y_pred_cv[, , j], y_val_cv)
        }
    }

    # Get the best lambda value
    best_lambda <- model_cv$lambda[which.min(colMeans(errors_cv))]

    # Construct Vhat_d and sigma_hat from the training validation
    Sigma0 <- calc_Sigma_j(y_train, 0)
    Sigma1 <- calc_Sigma_j(y_train, 1)
    Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
    sigma_hat <- vec_sigma_h(Sigma1, h1)
    Vhat_d <- construct_Vhat_d(Vhat)

    t0 <- Sys.time()
    XD1 <- Vhat_d %*% Dtilde_inv # Same as Vhat_d * inv(Dtilde)
    X1 <- XD1[, 1:m]
    X2 <- XD1[, (m + 1):dim(XD1)[2]]
    X2_plus <- solve((t(X2) %*% X2), t(X2)) # Same as inv(t(X2) %*% X2) %*% t(X2)

    # Transform the input to LASSO objective
    IminP <- diag(nrow(X2)) - X2 %*% X2_plus # Calculate I - P directly, P = X2 %*% X2_plus
    ytilde <- as.vector(IminP %*% t(sigma_hat))
    Xtilde <- IminP %*% X1

    # Fit the model on the entire training set
    model_cv <- glmnet(XD1, ytilde, lambda = best_lambda, alpha = 1, intercept = FALSE, ...)
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
        best_lambda = best_lambda,
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
    return(output_list)
}

fit_pvar.cv <- function(y, nlambdas = 20, verbose = FALSE, ...) {
    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, (floor(dim(y)[2] / 5) * 4 + 1):dim(y)[2]]

    # Fit a single solution using PVAR(1) with the BigVAR package
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Create cross-validation folds
    folds <- create_folds(y_train, nfolds = nfolds)

    # Perform cross-validation for the selection of the penalty parameter
    model_setup <- constructModel(
        Y = t(y),
        p = 1,
        struct = "Basic",
        gran = c(50, nlambdas),
        loss = "L2",
        T1 = floor(dim(y_train)[2] / 5) * 4 + 1,
        T2 = dim(y_train)[2],
        window.size = flooor(dim(y_train)[2] / 5), # Use 5-fold cross-validation
        model.controls = list(
            intercept = FALSE,
            loss = "L2"
        )
    )
    model_cv <- cv.BigVAR(model_setup)

    # Fit PVAR(1) using the optimal value for lambda
    t0 <- Sys.time()
    model <- BigVAR.fit(
        Y = t(y_train), # Take transpose becasue BigVAR.fit expects a matrix with rows as observations and columns as variables
        p = 1,
        struct = "Basic",
        lambda = model_cv@OptimalLambda,
        intercept = FALSE,
        ...
    )
    t1 <- Sys.time()

    # Return the fitted model
    output_list <- list(
        model = model,
        cvmodel = model_cv,
        errors_cv = model_cv@InSampMSFE,
        C = model[, , 1][, -1], # Remove the first column (intercept)
        runtime = difftime(t1, t0, units = "secs")[[1]],
        y_pred = predict_with_C(model[, , 1][, -1], y_train, y_test)
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
