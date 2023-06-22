setwd(system("echo $PROJ_DIR", intern = TRUE))
source("src/compute/R/utils.R")
source("src/compute/R/vhat_sigmahat.R")
setwd(system("echo $ZHU_DIR", intern = TRUE))
source("R/opt.R")
source("R/gen_data.R")
library(genlasso)
library(igraph)
library(splash)
library(BigVAR)
library(glmnet)
library(Matrix)
library(Rcpp)
setwd(system("echo $PROJ_DIR", intern = TRUE))
Rcpp::sourceCpp("src/compute/R/Rcpp/cpp_utils.cpp")

fit_splash.cv <- function(y, alpha, nlambdas = 20, nfolds = 5, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]

    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    # Create cross-validation folds
    folds <- rolling_cv(y_train, nfolds = nfolds)

    # Fit the model on each fold and save the results
    C_cv <- array(NA, dim = c(p, p, nlambdas))
    A_cv <- array(NA, dim = c(p, p, nlambdas))
    B_cv <- array(NA, dim = c(p, p, nlambdas))
    y_pred_cv <- array(NA, dim = c(p, length(folds[[1]]$test), nlambdas))
    errors_cv <- array(NA, dim = c(nfolds, nlambdas))
    for (i in 1:nfolds) {
        # Split the data into training and validation sets
        y_train_cv <- y_train[, folds[[i]]$train]
        y_val_cv <- y_train[, folds[[i]]$test]

        # Fit the model on the training set
        model_cv <- splash(t(y_train_cv), alphas = c(alpha), n_lambdas = nlambdas, track_progress = FALSE, ...)

        # Compute the prediction error on the validation set
        for (j in 1:nlambdas) {
            A_cv[, , j] <- model_cv$AB[, 1:p, j]
            B_cv[, , j] <- model_cv$AB[, (p + 1):(2 * p), j]
            C_cv[, , j] <- AB_to_C(A_cv[, , j], B_cv[, , j])
            y_pred_cv[, , j] <- predict_with_C(C_cv[, , j], y_train_cv, y_val_cv)
            errors_cv[i, j] <- calc_msfe(y_pred_cv[, , j], y_val_cv)
        }
    }
    # Get the best lambda value
    best_idx <- which.min(colMeans(errors_cv))

    # Fit the model on the entire training set
    t0 <- Sys.time()
    model <- splash(t(y_train), alphas = c(alpha), n_lambdas = nlambdas, track_progress = FALSE, ...)
    t1 <- Sys.time()

    # Extract the model output
    A <- model$AB[, 1:p, best_idx]
    B <- model$AB[, (p + 1):(2 * p), best_idx]
    C <- AB_to_C(A, B)
    y_pred <- predict_with_C(C, y_train, y_test)

    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = C,
        y_pred = y_pred,
        errors_cv = errors_cv,
        msfe = calc_msfe(y_test, y_pred),
        best_lambda = model$lambda[best_idx, 1],
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
    return(output_list)
}

fit_fsplash.cv <- function(y, bandwidths, graph, Dtilde_inv, nlambdas = 20, nfolds = 5, ...) {
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
    folds <- rolling_cv(y_train, nfolds = nfolds)

    # Fit the model on each fold and save the results
    C_cv <- array(NA, dim = c(p, p, nlambdas))
    A_cv <- array(NA, dim = c(p, p, nlambdas))
    B_cv <- array(NA, dim = c(p, p, nlambdas))
    y_pred_cv <- array(NA, dim = c(p, length(folds[[1]]$test), nlambdas))
    errors_cv <- array(NA, dim = c(nfolds, nlambdas))
    for (i in 1:nfolds) {
        # Split the data into training and validation sets
        y_train_cv <- y_train[, folds[[i]]$train]
        y_val_cv <- y_train[, folds[[i]]$test]

        # Construct Vhat_d and sigma_hat from the training validation
        Sigma0 <- calc_Sigma_j(y_train_cv, 0)
        Sigma1 <- calc_Sigma_j(y_train_cv, 1)
        Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
        sigma_hat <- construct_sigma_hat(Sigma1, h1)
        Vhat_d <- construct_Vhat_d(Vhat)

        # Transform Generalized Lasso problem into LASSO problem, using Rcpp to speed up matrix multiplications
        res <- cpp_FSPLASH_precalc(Vhat_d, Dtilde_inv, sigma_hat, m)

        # Fit the model on the training set
        model_cv <- glmnet(res$Xtilde, res$ytilde, nlambda = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4)

        # Compute the prediction error on the validation set
        for (j in 1:dim(model_cv$beta)[2]) {
            theta1 <- as.vector(model_cv$beta[, j])
            theta2 <- as.vector(res$X2_plus %*% (sigma_hat - res$X1 %*% theta1))
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
    best_idx <- which.min(colMeans(errors_cv))

    t0 <- Sys.time()
    # Construct Vhat_d and sigma_hat from the training validation
    Sigma0 <- calc_Sigma_j(y_train, 0)
    Sigma1 <- calc_Sigma_j(y_train, 1)
    Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
    sigma_hat <- construct_sigma_hat(Sigma1, h1)
    Vhat_d <- construct_Vhat_d(Vhat)

    # Transform Generalized Lasso problem into LASSO problem
    res <- fastComputation(Vhat_d, Dtilde_inv, sigma_hat, m)

    # Fit the model on the training set
    # model <- glmnet(res$Xtilde, res$ytilde, nlambda = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4, ...)
    model <- glmnet(res$Xtilde, res$ytilde, nlambda = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4)
    theta1 <- as.vector(model$beta[, best_idx])
    theta2 <- as.vector(res$X2_plus %*% (sigma_hat - res$X1 %*% theta1))
    coef <- Dtilde_inv %*% c(theta1, theta2)
    t1 <- Sys.time()
    toc()
    # Extract the model output
    AB <- coef_to_AB(coef, p)
    A <- AB$A
    B <- AB$B
    C <- AB_to_C(A, B)
    y_pred <- predict_with_C(C, y_train, y_test)

    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = C,
        errors_cv = errors_cv,
        y_pred = y_pred,
        msfe = calc_msfe(y_test, y_pred),
        best_lambda = model$lambda[best_idx],
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
    return(output_list)
}

fit_ssfsplash.cv <- function(y, bandwidths, graph, Dtilde_SSF_inv, alpha, nlambdas = 20, nfolds = 5, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]
    m <- ecount(graph)
    k <- vcount(graph)
    h <- floor((p - 1) / 4) # Bandwidth for the A and B matrix

    # Parse bandwidths
    h0 <- bandwidths[1]
    h1 <- bandwidths[2]

    # Calculate the multiplier matrix to balance the l1 error and the fusion error in SSF model
    # It is only needed to calculate once for each value of alpha
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

    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    # Create cross-validation folds
    folds <- rolling_cv(y_train, nfolds = nfolds)

    # Fit the model on each fold and save the results
    C_cv <- array(NA, dim = c(p, p, nlambdas))
    A_cv <- array(NA, dim = c(p, p, nlambdas))
    B_cv <- array(NA, dim = c(p, p, nlambdas))
    y_pred_cv <- array(NA, dim = c(p, length(folds[[1]]$test), nlambdas))
    errors_cv <- array(NA, dim = c(nfolds, nlambdas))
    for (i in 1:nfolds) {
        # Split the data into training and validation sets
        y_train_cv <- y_train[, folds[[i]]$train]
        y_val_cv <- y_train[, folds[[i]]$test]

        # Construct Vhat_d and sigma_hat from the training validation
        Sigma0 <- calc_Sigma_j(y_train_cv, 0)
        Sigma1 <- calc_Sigma_j(y_train_cv, 1)
        Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
        sigma_hat <- construct_sigma_hat(Sigma1, h1)
        Vhat_d <- construct_Vhat_d(Vhat)

        # Transform Generalized Lasso problem into LASSO problem
        Dtilde_SSF_inv_gamma <- Dtilde_SSF_inv %*% M_inv

        # Transform the input to LASSO objective (see Tibshirani and Taylor, 2011)
        XD1 <- as.matrix(Vhat_d %*% Dtilde_SSF_inv_gamma) # Same as Vhat_d * inv(Dtilde)

        # Fit the model on the training set
        # model_cv <- glmnet(XD1, sigma_hat, nlambda = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4, ...)
        model_cv <- glmnet(as.matrix(XD1), sigma_hat, nlambda = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4)

        # Compute the prediction error on the validation set

        for (j in 1:dim(model_cv$beta)[2]) {
            theta <- as.vector(model_cv$beta[, j])
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
    best_idx <- which.min(colMeans(errors_cv))
    tic()
    t0 <- Sys.time()
    # Construct Vhat_d and sigma_hat from the training validation
    Sigma0 <- calc_Sigma_j(y_train, 0)
    Sigma1 <- calc_Sigma_j(y_train, 1)
    Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
    sigma_hat <- construct_sigma_hat(Sigma1, h1)
    Vhat_d <- construct_Vhat_d(Vhat)

    Dtilde_SSF_inv_gamma <- Dtilde_SSF_inv %*% M_inv

    # Transform the input to LASSO objective (see Tibshirani and Taylor, 2011)
    XD1 <- Vhat_d %*% Dtilde_SSF_inv_gamma # Same as Vhat_d * inv(Dtilde)

    # Fit the model on the entire training set
    model <- glmnet(as.matrix(XD1), sigma_hat, n_lambdas = 20, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-8)
    theta <- as.vector(model$beta[, best_idx])
    coef <- Dtilde_SSF_inv_gamma %*% theta
    t1 <- Sys.time()
    toc()
    # Extract the model output
    AB <- coef_to_AB(coef, p)
    A <- AB$A
    B <- AB$B
    C <- AB_to_C(A, B)
    y_pred <- predict_with_C(C, y_train, y_test)

    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = C,
        errors_cv = errors_cv,
        y_pred = y_pred,
        msfe = calc_msfe(y_test, y_pred),
        best_lambda = model$lambda[best_idx],
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
    return(output_list)
}

fit_sdfsplash.cv <- function(y, bandwidths, graph, Dtilde_SDF_inv, alpha, nlambdas = 20, nfolds = 5, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]
    m <- ecount(graph)
    k <- vcount(graph)
    h <- floor((p - 1) / 4) # Bandwidth for the A and B matrix

    # Parse bandwidths
    h0 <- bandwidths[1]
    h1 <- bandwidths[2]

    # Calculate the multiplier matrix to balance the l1 error and the fusion error in SSF model
    # It is only needed to calculate once for each value of alpha
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

    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    # Create cross-validation folds
    folds <- rolling_cv(y_train, nfolds = nfolds)

    # Fit the model on each fold and save the results
    C_cv <- array(NA, dim = c(p, p, nlambdas))
    A_cv <- array(NA, dim = c(p, p, nlambdas))
    B_cv <- array(NA, dim = c(p, p, nlambdas))
    y_pred_cv <- array(NA, dim = c(p, length(folds[[1]]$test), nlambdas))
    errors_cv <- array(NA, dim = c(nfolds, nlambdas))
    for (i in 1:nfolds) {
        # Split the data into training and validation sets
        y_train_cv <- y_train[, folds[[i]]$train]
        y_val_cv <- y_train[, folds[[i]]$test]

        # Construct Vhat_d and sigma_hat from the training validation
        Sigma0 <- calc_Sigma_j(y_train_cv, 0)
        Sigma1 <- calc_Sigma_j(y_train_cv, 1)
        Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
        sigma_hat <- construct_sigma_hat(Sigma1, h1)
        Vhat_d <- construct_Vhat_d(Vhat)

        # Transform Generalized Lasso problem into LASSO problem
        Dtilde_SDF_inv_gamma <- Dtilde_SDF_inv %*% M_inv

        # Transform the input to LASSO objective (see Tibshirani and Taylor, 2011)
        XD1 <- Vhat_d %*% Dtilde_SDF_inv_gamma # Same as Vhat_d * inv(Dtilde)

        # Fit the model on the training set
        model_cv <- glmnet(XD1, sigma_hat, nlambda = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4, ...)

        # Compute the prediction error on the validation set

        for (j in 1:dim(model_cv$beta)[2]) {
            theta <- as.vector(model_cv$beta[, j])
            coef_cv <- Dtilde_SDF_inv_gamma %*% theta
            AB <- coef_to_AB(coef_cv, p)
            A_cv[, , j] <- AB$A
            B_cv[, , j] <- AB$B
            C_cv[, , j] <- AB_to_C(AB$A, AB$B)
            y_pred_cv[, , j] <- predict_with_C(C_cv[, , j], y_train_cv, y_val_cv)
            errors_cv[i, j] <- calc_msfe(y_pred_cv[, , j], y_val_cv)
        }
    }

    # Get the best lambda value
    best_idx <- which.min(colMeans(errors_cv))

    t0 <- Sys.time()
    # Construct Vhat_d and sigma_hat from the training validation
    Sigma0 <- calc_Sigma_j(y_train, 0)
    Sigma1 <- calc_Sigma_j(y_train, 1)
    Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
    sigma_hat <- construct_sigma_hat(Sigma1, h1)
    Vhat_d <- construct_Vhat_d(Vhat)

    Dtilde_SDF_inv_gamma <- Dtilde_SDF_inv %*% M_inv

    # Transform the input to LASSO objective (see Tibshirani and Taylor, 2011)
    XD1 <- Vhat_d %*% Dtilde_SDF_inv_gamma # Same as Vhat_d * inv(Dtilde)

    # Fit the model on the entire training set
    model <- glmnet(XD1, sigma_hat, n_lambdas = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4, ...)
    theta <- as.vector(model$beta[, best_idx])
    coef <- Dtilde_SDF_inv_gamma %*% theta
    t1 <- Sys.time()

    # Extract the model output
    AB <- coef_to_AB(coef, p)
    A <- AB$A
    B <- AB$B
    C <- AB_to_C(A, B)
    y_pred <- predict_with_C(C, y_train, y_test)

    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = C,
        errors_cv = errors_cv,
        y_pred = y_pred,
        msfe = calc_msfe(y_test, y_pred),
        best_lambda = model$lambda[best_idx],
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
    return(output_list)
}


fit_pvar.cv <- function(y, nlambdas = 20, nfolds = 20, ...) {
    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, (floor(dim(y)[2] / 5) * 4 + 1):dim(y)[2]]

    # Fit a single solution using PVAR(1) with the BigVAR package
    # Retrieve the cross-sectional dimension of the problem
    p <- nrow(y)

    # Create cross-validation folds
    folds <- rolling_cv(y_train, nfolds = nfolds)
    t0 <- Sys.time()
    # Perform cross-validation for the selection of the penalty parameter
    model_setup <- constructModel(
        Y = t(y_train),
        p = 1,
        struct = "Basic",
        gran = c(50, nlambdas),
        loss = "L2",
        verbose = FALSE,
        T1 = dim(y_train)[2] - nfolds - 1,
        T2 = dim(y_train)[2] - 1,
        window.size = floor(dim(y_train)[2] / 5) * 3, # Use 5-fold cross-validation
        model.controls = list(
            intercept = FALSE,
            loss = "L2"
        )
    )
    model_cv <- cv.BigVAR(model_setup)
    t1 <- Sys.time()
    # Fit PVAR(1) using the optimal value for lambda
    model <- BigVAR.fit(
        Y = t(y_train), # Take transpose becasue BigVAR.fit expects a matrix with rows as observations and columns as variables
        p = 1,
        struct = "Basic",
        lambda = model_cv@OptimalLambda,
        intercept = FALSE,
        ...
    )

    C <- model[, , 1][, -1] # Remove the first column (intercept)
    y_pred <- predict_with_C(C, y_train, y_test)

    # Return the fitted model
    output_list <- list(
        model = model,
        cvmodel = model_cv,
        C = C,
        y_pred = y_pred,
        msfe = calc_msfe(y_test, y_pred),
        best_lambda = model_cv@OptimalLambda,
        runtime = difftime(t1, t0, units = "secs")[[1]] / nfolds
    )
    return(output_list)
}

fit_gfsplash.on_idx <- function(y, bandwidths, lam_idx, alpha, graph, nlambdas = 20, lambda.min.ratio = 1e-4, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]
    m <- ecount(graph)
    k <- vcount(graph)
    h <- floor((p - 1) / 4) # Bandwidth for the A and B matrix

    # Parse bandwidths
    h0 <- bandwidths[1]
    h1 <- bandwidths[2]

    # Construct the D matrix from the penalty graph
    D <- as(getDgSparse(graph = graph), "TsparseMatrix")

    # Reparameterize the lambda and alpha for the linreg_path_v2 function, see paper Zhu 2017 Augmented ADMM
    gamma <- alpha / (1 - alpha)

    # Split y into training and testing sets
    train_idx <- (floor(dim(y)[2] / 5) * 4)
    y_train <- y[, 1:train_idx]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    t0 <- Sys.time()
    Sigma0 <- calc_Sigma_j(y_train, 0)
    Sigma1 <- calc_Sigma_j(y_train, 1)
    Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
    sigma_hat <- construct_sigma_hat(Sigma1, h1)
    Vhat_d <- as.matrix(construct_Vhat_d(Vhat))

    # Calculate lambda_0 and construct lambda grid accordingly
    lambda_0 <- calc_lambda_0_gfsplash(sigma_hat, Vhat_d, graph, alpha) # This gives lambda_0 in the (lambda, alpha) parametrization
    lambda_grid_linreg <- rev(lambda_0 * 10^seq(log10(lambda.min.ratio), 0, length.out = nlambdas)) / (1 + gamma)
    lambda_selected <- lambda_grid_linreg[lam_idx]

    # Fit the model on the training set
    model <- linreg_path_v2(
        Y = as.vector(sigma_hat),
        X = Vhat_d,
        val = D@x,
        idx = D@i + 1,
        jdx = D@j + 1,
        lambda_graph = lambda_selected,
        gamma = gamma,
        p = dim(Vhat_d)[2],
        m = dim(D)[1],
        ...
    )
    t1 <- Sys.time()

    # Collect output
    coef <- as.vector(model$beta_path[, 1])
    AB <- coef_to_AB(coef, p)
    A <- AB$A
    B <- AB$B
    C <- AB_to_C(A, B)
    y_pred <- predict_with_C(C, y_train, y_test)

    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = C,
        y_pred = y_pred,
        best_lambda = lambda_selected,
        msfe = calc_msfe(y_test, y_pred),
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
}


fit_gfsplash.cv <- function(y, bandwidths, graph, alpha, nlambdas = 20, nfolds = 5, lambda.min.ratio = 1e-4, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]
    m <- ecount(graph)
    k <- vcount(graph)
    h <- floor((p - 1) / 4) # Bandwidth for the A and B matrix

    # Parse bandwidths
    h0 <- bandwidths[1]
    h1 <- bandwidths[2]

    # Construct the D matrix from the penalty graph
    D <- as(getDgSparse(graph = graph), "TsparseMatrix")

    # Reparameterize the lambda and alpha for the linreg_path_v2 function, see paper Zhu 2017 Augmented ADMM
    gamma <- alpha / (1 - alpha)

    # Split y into training and testing sets
    train_idx <- (floor(dim(y)[2] / 5) * 4)
    y_train <- y[, 1:train_idx]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    Sigma0 <- calc_Sigma_j(y_train, 0)
    Sigma1 <- calc_Sigma_j(y_train, 1)
    Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
    sigma_hat <- construct_sigma_hat(Sigma1, h1)
    Vhat_d <- as.matrix(construct_Vhat_d(Vhat))

    # Calculate lambda_0 and construct lambda grid accordingly
    lambda_0 <- calc_lambda_0_gfsplash(sigma_hat, Vhat_d, graph, alpha) # This gives lambda_0 in the (lambda, alpha) parametrization
    lambda_grid_linreg <- as.vector(rev(lambda_0 * 10^seq(log10(lambda.min.ratio), 0, length.out = nlambdas)) / (1 + gamma))

    # Create cross-validation folds
    folds <- rolling_cv(y_train, nfolds = nfolds)

    # Fit the model on each fold and save the results
    C_cv <- array(NA, dim = c(p, p, nlambdas))
    A_cv <- array(NA, dim = c(p, p, nlambdas))
    B_cv <- array(NA, dim = c(p, p, nlambdas))
    y_pred_cv <- array(NA, dim = c(p, length(folds[[1]]$test), nlambdas))
    errors_cv <- array(NA, dim = c(nfolds, nlambdas))
    for (i in 1:nfolds) {
        # Split the data into training and validation sets
        y_train_cv <- y_train[, folds[[i]]$train]
        y_val_cv <- y_train[, folds[[i]]$test]

        t0 <- Sys.time()
        # Construct Vhat_d and sigma_hat from the training validation
        Sigma0_cv <- calc_Sigma_j(y_train_cv, 0)
        Sigma1_cv <- calc_Sigma_j(y_train_cv, 1)
        Vhat_cv <- construct_Vhat(Sigma0_cv, Sigma1_cv, h0, h1)
        sigma_hat_cv <- construct_sigma_hat(Sigma1_cv, h1)
        Vhat_d_cv <- construct_Vhat_d(Vhat_cv)

        # Fit the model on the training set
        model_cv <- linreg_path_v2(
            Y = as.vector(sigma_hat_cv),
            X = as.matrix(Vhat_d_cv),
            val = D@x,
            idx = D@i + 1,
            jdx = D@j + 1,
            lambda_graph = lambda_grid_linreg,
            gamma = gamma,
            p = dim(Vhat_d)[2],
            m = dim(D)[1],
            ...
        )
        t1 <- Sys.time()
        # Compute the prediction error on the validation set
        for (j in 1:dim(model_cv$beta)[2]) {
            coef_cv <- as.vector(model_cv$beta_path[, j])
            AB <- coef_to_AB(coef_cv, p)
            A_cv[, , j] <- AB$A
            B_cv[, , j] <- AB$B
            C_cv[, , j] <- AB_to_C(AB$A, AB$B)
            y_pred_cv[, , j] <- predict_with_C(C_cv[, , j], y_train_cv, y_val_cv)
            errors_cv[i, j] <- calc_msfe(y_pred_cv[, , j], y_val_cv)
        }
    }

    # Get the best lambda value
    best_idx <- which.min(colMeans(errors_cv))
    best_lambda <- lambda_grid_linreg[best_idx]
    # Fit the model on the training set
    model <- linreg_path_v2(
        Y = as.vector(sigma_hat),
        X = as.matrix(Vhat_d),
        val = D@x,
        idx = D@i + 1,
        jdx = D@j + 1,
        lambda_graph = best_lambda,
        gamma = gamma,
        p = dim(Vhat_d)[2],
        m = dim(D)[1],
        ...
    )
    coef <- as.vector(model$beta_path[, 1])

    # Extract the model output
    AB <- coef_to_AB(coef, p)
    A <- AB$A
    B <- AB$B
    C <- AB_to_C(A, B)
    y_pred <- predict_with_C(C, y_train, y_test)

    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = C,
        errors_cv = errors_cv,
        y_pred = y_pred,
        msfe = calc_msfe(y_test, y_pred),
        best_lambda = best_lambda,
        lambda_grid = lambda_grid_linreg,
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
    return(output_list)
}
