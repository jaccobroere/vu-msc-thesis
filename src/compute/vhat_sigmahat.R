library(Matrix)

calc_Sigma_j <- function(y, j) {
    p <- nrow(y)
    T <- ncol(y)
    sigma_j <- matrix(0, p, p)
    ymean <- rowMeans(y)
    for (t in 1:(T - j)) {
        sigma_j <- sigma_j + (y[, t] - ymean) %*% t(y[, t + j] - ymean)
    }
    return(sigma_j / T)
}

band_matrix <- function(A, h) {
    N <- nrow(A)
    T <- ncol(A)
    B <- matrix(0, N, T)
    for (i in 1:N) {
        for (j in max((i - h), 1):min((i + h), T)) {
            B[i, j] <- A[i, j]
        }
    }
    return(B)
}

constr_Vhat <- function(Sigma0, Sigma1, h0 = 0, h1 = 0) {
    if (h0 == 0 || h1 == 0) {
        return(cbind(t(Sigma1), Sigma0))
    }
    Sigma0 <- band_matrix(Sigma0, h0)
    Sigma1 <- band_matrix(Sigma1, h1)
    return(cbind(t(Sigma1), Sigma0))
}

vec_sigma_h <- function(Sigma1, h1 = 0) {
    if (h1 == 0) {
        return(as.vector(t(Sigma1)))
    }
    res <- band_matrix(Sigma1, h1)
    return(as.vector(t(res)))
}

active_cols <- function(p, bandwidth = 0) {
    if (bandwidth == 0) {
        bandwidth <- floor(p / 4)
    }

    active_set <- matrix(FALSE, p, p * 2)

    for (i in 1:p) {
        full_a <- max(1, i - bandwidth):min(i + bandwidth, p)
        full_b <- (p + max(1, (i - bandwidth))):(p + min(i + bandwidth, p))
        selection <- c(full_a, full_b)

        active_set[i, selection] <- TRUE
        active_set[i, i] <- FALSE
    }
    return(active_set)
}

constr_Vhat_d <- function(V, bandwidth = 0) {
    p <- nrow(V)
    active <- active_cols(p, bandwidth)
    Vhat_d <- Matrix(0, p^2, sum(active), sparse = TRUE)

    row_index <- col_index <- idx <- 1
    for (idx in 1:p) {
        M <- V[, active[idx, ]]
        Vhat_d[row_index:(row_index + p - 1), col_index:(col_index + dim(M)[2] - 1)] <- M
        row_index <- row_index + p
        col_index <- col_index + dim(M)[2]
    }
    return(Vhat_d)
}



library(caret)

dum <- c(1:500)

cv <- createTimeSlices(dum, initialWindow = 400, fixedWindow = TRUE, horizon = 20, skip = 19)

cv$train[[3]]
cv$test[[5]]
