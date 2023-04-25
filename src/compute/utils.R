coef_to_AB <- function(coef, p) {
    bandwidth <- p %/% 4
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
    bandwidth <- p %/% 4
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

predict_with_C <- function(C, y) {
    train_idx <- 1:floor(dim(y)[2] / 5) * 4
    y_train <- y[, 1:train_idx]
    y_test <- y[, train_idx:ncol(y)] # Notice leaving out + 1 here to have a predictions for the first element of y_test
    predictions <- matrix(0, nrow = nrow(y), ncol = ncol(y_test))

    # Predictions
    predictions[, 1] <- C %*% y_train[, ncol(y_train)] # First prediction is based on last element of y_train
    for (i in 2:ncol(y_test)) {
        predictions[, i] <- C %*% y_test[, i - 1]
    }
    return(predictions)
}

# Calcualte lambda_0 for the generalized lasso cases
calc_lambda_0 <- function(sigma_hat, Vhat_d, graph, alpha, ...) {
    # Check validity of alpha value
    if (alpha < 0 | alpha >= 1) stop("alpha must be between 0 and 1")

    # Fit genlasso for only the first two iterations (1 iteration does not work)
    model <- fusedlasso(sigma_hat, Vhat_d, graph = graph, gamma = alpha / (1 - alpha), maxsteps = 2, ...)

    # Retreive the lambda from the first iteration
    lambda_0 <- coef(model)$lambda[1]

    # Transform lambda_0 back to (lambda, alpha) parametrization
    lambda_res <- lambda_0 * (1 + alpha / (1 - alpha))
    return(lambda_res)
}
