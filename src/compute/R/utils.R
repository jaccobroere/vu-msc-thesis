library(Matrix)

coef_to_AB <- function(coef, p) {
    bandwidth <- (p - 1) %/% 4
    AB <- matrix(0, nrow = p, ncol = 2 * p)
    cnt <- 1
    for (i in 1:p) {
        for (j in 1:(2 * p)) {
            if (j <= p) { # Belongs to A
                if (abs(i - j) <= bandwidth && i != j) {
                    AB[i, j] <- coef[cnt]
                    cnt <- cnt + 1
                }
            } else { # Belongs to B
                if (abs(i - abs(j - p)) <= bandwidth) {
                    AB[i, j] <- coef[cnt]
                    cnt <- cnt + 1
                }
            }
        }
    }
    output_list <- list(A = AB[, 1:p], B = AB[, (p + 1):(2 * p)])
    return(output_list)
}

AB_to_coef <- function(AB, p) {
    bandwidth <- (p - 1) %/% 4
    coef <- rep(0, p * (4 * bandwidth + 1) - 2 * (bandwidth^2 + bandwidth))
    cnt <- 1
    for (i in 1:p) {
        for (j in 1:(2 * p)) {
            if (j <= p) { # Belongs to A
                if (abs(i - j) <= bandwidth && i != j) {
                    coef[cnt] <- AB[i, j]
                    cnt <- cnt + 1
                }
            } else { # Belongs to B
                if (abs(i - abs(j - p)) <= bandwidth) {
                    coef[cnt] <- AB[i, j]
                    cnt <- cnt + 1
                }
            }
        }
    }
    return(coef)
}

AB_to_C <- function(A, B) {
    C <- solve(diag(nrow(A)) - A) %*% B
    return(C)
}

calc_rmsfe <- function(y_true, y_pred, y_true_hat) {
    # Calculate the relative mean squared forecast error
    sum_hat <- sum((y_true - y_pred)^2)
    sum_true <- sum((y_true - y_true_hat)^2)
    return(sum_hat / sum_true)
}

calc_msfe <- function(y_true, y_pred) {
    # Calculate the mean squared forecast error
    return(mean((y_true - y_pred)^2))
}

predict_with_C <- function(C_hat, y_train, y_test) {
    if (ncol(y_test) == 1) {
        return(C_hat %*% y_train[, ncol(y_train)])
    }
    predictions <- matrix(0, nrow = nrow(y_train), ncol = ncol(y_test))

    # Predictions
    predictions[, 1] <- C_hat %*% y_train[, ncol(y_train)] # First prediction is based on last element of y_train
    for (i in 2:ncol(y_test)) {
        predictions[, i] <- C_hat %*% y_test[, i - 1]
    }

    return(predictions)
}

# Calcualte lambda_0 for the generalized lasso cases
calc_lambda_0_gfsplash <- function(sigma_hat, Vhat_d, graph, alpha, ...) {
    # Check validity of alpha value
    if (alpha < 0 | alpha >= 1) stop("alpha must be between in [0, 1), for alpha=1, use a standard lasso solver instead")

    # Fit genlasso for only the first two iterations (1 iteration does not work)
    model <- fusedlasso(sigma_hat, Vhat_d, graph = graph, gamma = alpha / (1 - alpha), maxsteps = 2, ...)

    # Retreive the lambda from the first iteration
    lambda_0 <- coef(model)$lambda[1]

    # Transform lambda_0 back to (lambda, alpha) parametrization
    lambda_res <- lambda_0 * (1 + alpha / (1 - alpha))
    return(lambda_res)
}

calc_Dtilde <- function(graph) {
    null_vecs <- null_space_graph(graph)
    Dtilde <- rbind(as.matrix(getDg(graph)), t(null_vecs))
    return(Dtilde)
}

calc_Dtilde_sparse <- function(graph) {
    null_vecs <- null_space_graph(graph)
    Dtilde <- rbind(as.matrix(getDg(graph)), t(null_vecs))
    Dtilde <- Matrix(Dtilde, sparse = TRUE)
    return(Dtilde)
}

null_space_graph <- function(graph) {
    null_vecs <- zeros(vcount(graph), (vcount(graph) - ecount(graph)))
    conn <- components(graph)
    mem <- conn$membership
    for (i in 1:length(mem)) {
        null_vecs[i, mem[i]] <- 1
    }
    return(null_vecs)
}

calc_EE <- function(M_true, M_hat, type = "2") {
    # type="F" for Frobenius, type="1" for max absolute column norm,  uses Matrix package
    return(norm(M_true - M_hat, type = "2"))
}

rolling_cv <- function(y, nfolds, test_size = 0.2) {
    # Get the number of rows in the matrix
    ncols <- ncol(y)

    # Initialize a list to store the training and validation sets
    cv_folds <- list()

    # Calculate the window sizes
    train_window_size <- ncols - (floor(ncols * test_size))
    test_window_size <- floor(ncols * test_size / nfolds)

    idx <- 1
    for (i in 1:nfolds) {
        # The training set is the current window
        training_set <- idx:(idx + train_window_size - 1)

        # The validation set is the row immediately after the window
        test_set <- (idx + train_window_size):(idx + train_window_size + test_window_size - 1)

        # Add the training and validation sets to the list
        cv_folds[[i]] <- list("train" = training_set, "test" = test_set)

        # Update the index
        idx <- idx + test_window_size
    }

    # Return the list of training and validation sets
    return(cv_folds)
}
