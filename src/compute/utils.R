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
