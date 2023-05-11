fit_fgsg_gsplash <- function(sigma_hat, Vhat_d, graph, lambda1, lambda2, ...) {
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Create edge vector
    edge_vector <- as.vector(t(as_edgelist(graph)))

    # Fit FGSG GFLASSO implementation
    t0 <- Sys.time()
    model <- FGSG::gflasso(y = sigma_hat, A = Vhat_d, tp = edge_vector, s1 = lambda2, s2 = lambda1, ...) # Notice that \lambda1 and \lambda2 are swapped here
    runtime <- difftime(Sys.time(), t0, units = "secs")[[1]]

    # Print how long it took to run formatted with a message
    message(paste0("FGSG took ", round(runtime, 2), " seconds to run."))

    # Transform back to A, B
    coef <- model$weight
    AB <- coef_to_AB(coef, p)
    A <- AB$A
    B <- AB$B

    # Return the fitted model
    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = AB_to_C(A, B),
        runtime = runtime
    )
}

get_lam_best <- function(df, model) {
    return(df[df$model == model, "lambda"])
}

save_fitting_results <- function(model, prefix, fit_dir) {
    fwrite(data.table(model$A), file = file.path(fit_dir, paste0(prefix, "_estimate_A.csv")))
    fwrite(data.table(model$B), file = file.path(fit_dir, paste0(prefix, "_estimate_B.csv")))
    fwrite(data.table(model$C), file = file.path(fit_dir, paste0(prefix, "_estimate_C.csv")))
    fwrite(data.table(model$yhat), file = file.path(fit_dir, paste0(prefix, "_yhat.csv")))
}

run_lambda_finder_gfsplash <- function(y, sigma_hat, Vhat_d, C_true, graph, alpha, path) {
    train_idx <- (floor(dim(y)[2] / 5) * 4)
    y_train <- y[, 1:train_idx]
    y_test <- y[, (train_idx + 1):ncol(y)]
    # Calculate lambda_0 for the GSPLASH
    lam0 <- calc_lambda_0_gfsplash(sigma_hat, Vhat_d, graph, alpha = alpha)
    # Generate grid of values for lambda
    grid_lam <- gen_lambda_grid(lam0, length.out = 20)
    # Call the function for each specification and create file if it does not exist yet
    df_lam <- create_lambda_df(grid_lam, path)
    # Fit the models and save the prediction results in the data.frame
    for (i in 1:length(grid_lam)) {
        lam <- grid_lam[i]
        model <- fit_admm_gsplash(sigma_hat, Vhat_d, graph, lam, alpha = alpha)
        # Calculate predictions and error metric
        y_pred <- predict_with_C(model$C, y)
        y_true_pred <- predict_with_C(C_true, y)
        error_metric <- calc_msfe(y_test, y_pred)
        df_lam[1, colnames(df_lam)[i]] <- error_metric
    }

    # Append the results to the table
    write.table(df_lam, path, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    return(df_lam)
}

run_lambda_finder_fsplash <- function(y, sigma_hat, Vhat_d, C_ture, graph, Dtilde_inv, path) {
    train_idx <- (floor(dim(y)[2] / 5) * 4)
    y_train <- y[, 1:train_idx]
    y_test <- y[, (train_idx + 1):ncol(y)]
    # Calculate lambda_0 for the GSPLASH
    df_lam <- create_lambda_df(grid_lam, path)
    # Fit the models and save the prediction results in the data.frame
    model <- fit_fsplash(sigma_hat, Vhat_d, graph, Dtilde_inv, lambda = NULL)
    # Calculate predictions and error metric
    for (i in 1:length(grid_lam)) {
        lam <- model$model$lambda[i]
        y_pred <- predict_with_C(res$C[, , i], y)
        error_metric <- calc_msfe(y_test, y_pred)
        df_lam[1, colnames(df_lam)[i]] <- error_metric
    }

    # Append the results to the table
    write.table(df_lam, path, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
}

run_lambda_finder_splash <- function(y, alpha, C_true, path, lambda_min_mult = 1e-4) {
    train_idx <- (floor(dim(y)[2] / 5) * 4)
    y_train <- y[, 1:train_idx]
    y_test <- y[, (train_idx + 1):ncol(y)]
    # Split the training set off
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    # Run the SPLASH model with the given lambda_grid
    # model <- splash::splash(t(y_train), banded_covs = c(TRUE, TRUE), B = 500, n_lambdas = 20, alphas = c(alpha), lambda_min_mult = lambda_min_mult)
    lambda_grid <- gen_lambda_grid(lambda_min_mult, length.out = 20)
    model <- splash::splash(t(y_train), banded_covs = c(TRUE, TRUE), B = 500, lambdas = lambda_grid, alphas = c(alpha))
    p <- dim(model$AB)[1]
    # Create placeholder dataframe
    df_lam <- create_lambda_df(c(model$lambdas), path)
    # Generate and save predictions for each lambda value
    for (i in 1:dim(model$AB)[3]) {
        A <- model$AB[, 1:p, i]
        B <- model$AB[, (p + 1):(2 * p), i]
        C <- AB_to_C(A, B)
        # Calculate predictions and error metric
        y_pred <- predict_with_C(C, y)
        y_true_pred <- predict_with_C(C_true, y)
        error_metric <- calc_msfe(y_test, y_pred)
        df_lam[1, colnames(df_lam)[i]] <- error_metric
    }

    # Append the results to the table
    write.table(df_lam, path, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    return(df_lam)
}

# Create dataframes with lambda values as column names
create_lambda_df <- function(lambda_grid, filename) {
    df <- data.frame(matrix(ncol = length(lambda_grid), nrow = 0))
    colnames(df) <- c(1:length(lambda_grid))

    # If the file does not exist, create it
    if (!file.exists(filename)) {
        # Save the dataframe as a .csv file
        write.csv(df, file = filename, row.names = FALSE, quote = FALSE)
    }
    return(df)
}
