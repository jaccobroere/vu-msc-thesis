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

construct_Vhat <- function(Sigma0, Sigma1, h0 = 0, h1 = 0) {
    if (h0 == 0 || h1 == 0) {
        return(cbind(t(Sigma1), Sigma0))
    }
    Sigma0 <- band_matrix(Sigma0, h0)
    Sigma1 <- band_matrix(Sigma1, h1)
    return(cbind(t(Sigma1), Sigma0))
}

construct_sigma_hat <- function(Sigma1, h1 = 0) {
    if (h1 == 0) {
        return(as.vector(t(Sigma1)))
    }
    res <- band_matrix(Sigma1, h1)
    return(as.vector(t(res)))
}

active_cols <- function(p, bandwidth = 0) {
    if (bandwidth == 0) {
        bandwidth <- floor((p - 1) / 4)
    }

    active_set <- matrix(FALSE, p, p * 2)

    for (i in 1:p) {
        full_a <- max(1, i - bandwidth):min(i + bandwidth, p)
        full_b <- (p + max(1, (i - bandwidth))):(p + min(i + bandwidth, p))
        selection <- c(full_a, full_b)

        active_set[i, selection] <- TRUE

        # Diagonal elements of A are set to FALSE
        active_set[i, i] <- FALSE
    }
    return(active_set)
}

active_cols_alt <- function(p, m) {
    active_set_A <- matrix(FALSE, p, p * 2)
    active_set_B <- matrix(FALSE, p, p * 2)

    for (i in 1:p) {
        for (j in 1:p) {
            # Vertical neighbours
            if (abs(i - j) == m) {
                active_set_A[i, j] <- TRUE
                active_set_B[i, j] <- TRUE
            }

            # Conditions
            not_side_edge <- !((i %% m == 0 && (j - 1) %% m == 0) || ((i - 1) %% m == 0 && j %% m == 0))
            not_top_edge <- (i > m) && (j > m)
            not_bottom_edge <- (i <= (p - m)) && (j <= (m * (m - 1)))

            # Horizontal neighbours
            if (abs(i - j) == 1 && not_side_edge) {
                active_set_A[i, j] <- TRUE
                active_set_B[i, j] <- TRUE
            }

            # Bottom right diagonal neighbour
            if (not_side_edge && not_bottom_edge && abs(i - j) == (m + 1)) {
                active_set_A[i, j] <- TRUE
                active_set_B[i, j] <- TRUE
            }

            # Top right diagonal neighbour
            if (not_side_edge && not_top_edge && abs(i - j) == (m - 1)) {
                active_set_A[i, j] <- TRUE
                active_set_B[i, j] <- TRUE
            }

            # Bottom left diagonal neighbour
            if (not_bottom_edge && not_side_edge && abs(i - j) == (m - 1)) {
                active_set_A[i, j] <- TRUE
                active_set_B[i, j] <- TRUE
            }

            # Top left diagonal neighbour
            if (not_top_edge && not_side_edge && abs(i - j) == (m + 1)) {
                active_set_A[i, j] <- TRUE
                active_set_B[i, j] <- TRUE
            }

            # Zero diagonal for A
            if (i == j) {
                active_set_A[i, j] <- FALSE
            }
        }
    }
    active_set <- cbind(active_set_A, active_set_B)

    return(active_set)
}

construct_Vhat_d <- function(V, bandwidth = 0) {
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
