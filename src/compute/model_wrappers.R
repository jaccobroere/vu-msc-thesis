setwd(system("echo $PROJ_DIR", intern = TRUE))
source("src/compute/utils.R")
setwd(system("echo $ZHU_DIR", intern = TRUE))
source("R/opt.R")
source("R/gen_data.R")
library(genlasso)
library(igraph)
library(splash)
library(FGSG)
library(BigVAR)
library(glmnet)
library(pracma)
library(Rlinsolve)
library(Matrix)
setwd(system("echo $PROJ_DIR", intern = TRUE))

fit_gfsplash <- function(sigma_hat, Vhat_d, graph = NULL, D = NULL, lambda, alpha, verbose = FALSE, scale = TRUE, ...) {
    if (scale) {
        # Scale the data
        sd_y <- sd(sigma_hat)
        sigma_hat <- scale(as.vector(sigma_hat), center = FALSE)
        sd_x <- apply(Vhat_d, 2, sd)
        Vhat_d <- scale(Vhat_d)
    }
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Reparameterize the lambda and alpha
    gamma <- alpha / (1 - alpha)
    lambda_admm <- lambda / (1 + gamma)

    # Convert graph to sparse matrix
    t0 <- Sys.time()
    D <- as(getDgSparse(graph = graph), "TsparseMatrix")
    idx <- D@i + 1
    jdx <- D@j + 1
    val <- D@x
    runtimeD <- difftime(Sys.time(), t0, units = "secs")[[1]]

    # Fit linear regression model using ADMM
    t0 <- Sys.time()
    model <- linreg_path_v2(
        Y = as.vector(sigma_hat),
        X = Vhat_d,
        val = val,
        idx = idx,
        jdx = jdx,
        lambda_graph = lambda_admm,
        gamma = gamma,
        p = dim(Vhat_d)[2],
        m = dim(D)[1],
        ...
    )
    runtimeM <- difftime(Sys.time(), t0, units = "secs")[[1]]

    t0 <- Sys.time()
    coef <- as.vector(model$beta_path[, 1])
    if (scale) {
        coef <- coef / sd_x * sd_y
    }
    AB <- coef_to_AB(coef, p)
    A <- AB$A
    B <- AB$B
    runtimeC <- difftime(Sys.time(), t0, units = "secs")[[1]]

    # Print how long it took to run formatted with a message
    if (verbose) {
        message(paste0("ADMM took ", round(runtimeM, 2), " seconds to run for the model."))
    }

    # Return the fitted model
    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = AB_to_C(A, B),
        runtimeD = runtimeD,
        runtimeM = runtimeM,
        runtimeC = runtimeC
    )
    return(output_list)
}

fit_splash <- function(y, lambda, alpha, verbose = FALSE, ...) {
    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(dim(y)[1])

    # Fit SPLASH from Reuvers and Wijler (2021)
    t0 <- Sys.time()
    model <- splash(t(y_train), alphas = c(alpha), lambdas = c(lambda), ...) # Take transpose because splash() accepts T x p matrix
    runtime <- difftime(Sys.time(), t0, units = "secs")[[1]]

    # Print how long it took to run formatted with a message
    if (verbose) {
        message(paste0("SPLASH took ", round(runtime, 2), " seconds to run."))
    }

    A <- model$AB[, 1:p, ]
    B <- model$AB[, (p + 1):(2 * p), ]
    # Return the fitted model
    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = AB_to_C(A, B),
        runtime = runtime
    )
    return(output_list)
}

fit_pvar <- function(y, verbose = FALSE, ...) {
    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, (floor(dim(y)[2] / 5) * 4 + 1):dim(y)[2]]

    # Fit a single solution using PVAR(1) with the BigVAR package
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Compute lambda grid
    lambda_grid <- gen_lambda_grid(20)

    # Perform cross-validation for the selection of the penalty parameter
    t0 <- Sys.time()
    model <- constructModel(
        Y = t(y),
        p = 1,
        struct = "Basic",
        gran = lambda_grid,
        loss = "L2",
        T1 = floor(dim(y_train)[2] / 5) * 4 + 1,
        T2 = dim(y_train)[2],
        ownlambdas = TRUE,
        model.controls = list(
            intercept = FALSE,
            loss = "L2"
        )
    )
    cvmodel <- cv.BigVAR(model)
    runtimeCV <- difftime(Sys.time(), t0, units = "secs")[[1]]

    # Fit PVAR(1) using the optimal value for lambda
    t0 <- Sys.time()
    model <- BigVAR.fit(
        Y = t(y_train), # Take transpose becasue BigVAR.fit expects a matrix with rows as observations and columns as variables
        p = 1,
        struct = "Basic",
        lambda = cvmodel@OptimalLambda,
        intercept = FALSE,
        ...
    )
    runtime <- difftime(Sys.time(), t0, units = "secs")[[1]]

    # Print how long it took to run formatted with a message
    if (verbose) {
        message(paste0("PVAR took ", round(runtime, 2), " seconds to run for the model."))
        message(paste0("PVAR took ", round(runtimeCV, 2), " seconds to run for cross-validation."))
    }

    # Return the fitted model
    output_list <- list(
        model = model,
        cvmodel = cvmodel,
        C = model[, , 1][, -1], # Remove the first column (intercept)
        runtime = runtime,
        runtimeCV = runtimeCV
    )
    return(output_list)
}

fit_fsplash <- function(sigma_hat, Vhat_d, graph, Dtilde_inv, lambda = NULL, nlambda = 20) {
    m <- ecount(graph)
    k <- vcount(graph)
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    t0 <- Sys.time()
    # Use linear system solvers for faster computation of the change of variables (see Tibshirani and Taylor, 2011)
    XD1 <- Vhat_d %*% Dtilde_inv # Same as Vhat_d * inv(Dtilde)
    X1 <- XD1[, 1:m]
    X2 <- XD1[, (m + 1):dim(XD1)[2]]
    X2_plus <- solve((t(X2) %*% X2), t(X2)) # Same as inv(t(X2) %*% X2) %*% t(X2)

    # Transform the input to LASSO objective
    P <- X2 %*% X2_plus
    ytilde <- as.vector((diag(nrow(P)) - P) %*% t(sigma_hat))
    Xtilde <- (diag(nrow(P)) - P) %*% X1
    runtimeXtilde <- difftime(Sys.time(), t0, units = "secs")[[1]]

    t0 <- Sys.time()
    if (is.null(lambda)) {
        model <- glmnet(Xtilde, ytilde, nlambda = nlambda, alpha = 1, intercept = FALSE, standardize = TRUE)
        # Loop over fitted solution for all lambdas
        C <- array(NA, dim = c(p, p, nlambda))
        A <- array(NA, dim = c(p, p, nlambda))
        B <- array(NA, dim = c(p, p, nlambda))
        coef <- array(NA, dim = c(k, nlambda))
        for (i in 1:nlambda) {
            theta1 <- as.vector(model$beta[, i])
            theta2 <- as.vector(X2_plus %*% (t(sigma_hat) - X1 %*% theta1))
            c <- Dtilde_inv %*% c(theta1, theta2)
            AB <- coef_to_AB(c, p)
            A[, , i] <- AB$A
            B[, , i] <- AB$B
            C[, , i] <- AB_to_C(AB$A, AB$B)
            coef[, i] <- as.vector(c)
        }
    } else {
        # model <- glmnet(Xtilde, ytilde, lambda = lambda, alpha = 1, intercept = FALSE, standardize = FALSE)
        model <- glmnet(Xtilde, ytilde, lambda = lambda, alpha = 1, intercept = FALSE, standardize = TRUE)
        theta1 <- as.vector(model$beta)
        theta2 <- as.vector(X2_plus %*% (t(sigma_hat) - X1 %*% theta1))
        coef <- Dtilde_inv %*% c(theta1, theta2)
        AB <- coef_to_AB(coef, p)
        A <- AB$A
        B <- AB$B
        C <- AB_to_C(A, B)
    }

    runtimeM <- difftime(Sys.time(), t0, units = "secs")[[1]]

    # Return the fitted model
    output_list <- list(
        model = model,
        coef = coef,
        A = A,
        B = B,
        C = C,
        runtimeM = runtimeM,
        runtimeXtilde = runtimeXtilde
    )
}

fit_ssfsplash <- function(sigma_hat, Vhat_d, Dtilde_SSF_inv, lambda, alpha) {
    # Get dimensions of the problem
    # m <- ecount(graph)
    # k <- vcount(graph)
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Reparametrize to gamma
    gamma <- alpha / (1 - alpha)
    lambda_star <- lambda * (1 - alpha)

    # Calculate the scaled version of Dtilde_SSF_inv, the multiplier matrix M, has to be inverted
    # but is a diagonal matrix, so we now the explicit form of the inverse and calculate it directly
    # M_inv <- .sparseDiagonal(x = c(rep(1, m), rep((1 / gamma), k - m)))
    Dtilde_SSF_inv_gamma <- Dtilde_SSF_inv # %*% M_inv

    # Transform the input to LASSO objective (see Tibshirani and Taylor, 2011)
    XD1 <- Vhat_d %*% Dtilde_SSF_inv_gamma # Same as Vhat_d * inv(Dtilde)

    t0 <- Sys.time()
    if (is.null(lambda)) {
        model <- glmnet(XD1, as.vector(sigma_hat), nlambda = nlambda, alpha = 1, intercept = FALSE, standardize = TRUE)
        # Loop over fitted solution for all lambdas
        C <- array(NA, dim = c(p, p, nlambda))
        A <- array(NA, dim = c(p, p, nlambda))
        B <- array(NA, dim = c(p, p, nlambda))
        coef <- array(NA, dim = c(k, nlambda))
        for (i in 1:nlambda) {
            theta <- as.vector(model$beta[, i])
            c <- Dtilde_SSF_inv_gamma %*% theta
            AB <- coef_to_AB(c, p)
            A[, , i] <- AB$A
            B[, , i] <- AB$B
            C[, , i] <- AB_to_C(AB$A, AB$B)
            coef[, i] <- as.vector(c)
        }
    } else {
        # model <- glmnet(XD1, as.vector(sigma_hat), lambda = lambda, alpha = 1, intercept = FALSE, standardize = FALSE)
        model <- glmnet(XD1, as.vector(sigma_hat), lambda = lambda, alpha = 1, intercept = FALSE, standardize = TRUE)
        theta <- as.vector(model$beta)
        coef <- Dtilde_SSF_inv_gamma %*% theta
        AB <- coef_to_AB(coef, p)
        A <- AB$A
        B <- AB$B
        C <- AB_to_C(A, B)
    }

    runtimeM <- difftime(Sys.time(), t0, units = "secs")[[1]]
    # Return the fitted model
    output_list <- list(
        model = model,
        coef = coef,
        A = A,
        B = B,
        C = C,
        runtimeM = runtimeM
    )
}

fit_fgsg_gsplash <- function(sigma_hat, Vhat_d, graph, lambda, alpha, verbose = FALSE, ...) {
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Reparamaterize lambda
    lambda2 <- lambda * alpha # l1 penalty
    lambda1 <- lambda * (1 - alpha) # Fusion penalty

    # Create edge vector
    edge_vector <- as.vector(t(as_edgelist(graph)))

    # Fit FGSG GFLASSO implementation
    t0 <- Sys.time()
    model <- FGSG::gflasso(y = sigma_hat, A = Vhat_d, tp = edge_vector, s1 = lambda1, s2 = lambda2, ...) # Notice that \lambda1 and \lambda2 are swapped here
    runtime <- difftime(Sys.time(), t0, units = "secs")[[1]]

    # Print how long it took to run formatted with a message
    if (verbose) {
        message(paste0("FGSG took ", round(runtime, 2), " seconds to run."))
    }

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
