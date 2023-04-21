setwd(system("echo $ROOT_DIR", intern = TRUE))
source("src/compute/utils.R")
setwd(system("echo $ZHU_DIR", intern = TRUE))
source("R/opt.R")
source("R/gen_data.R")
library(genlasso)
library(igraph)
library(splash)

admm_gsplash <- function(sigma_hat, Vhat_d, graph, lambda, gamma, ...) {
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Convert graph to sparse matrix
    t0 <- Sys.time()
    D <- as(getDgSparse(graph = graph), "TsparseMatrix")
    idx <- D@i + 1
    jdx <- D@j + 1
    val <- D@x
    runtimeD <- Sys.time() - t0

    # Fit linear regression model using ADMM
    t0 <- Sys.time()
    model <- linreg_path_v2(
        Y = sigma_hat,
        X = Vhat_d,
        val = val,
        idx = idx,
        jdx = jdx,
        lambda_graph = lambda,
        gamma = gamma,
        p = dim(Vhat_d)[2],
        m = dim(D)[1],
        ...
    )
    runtimeM <- Sys.time() - t0

    t0 <- Sys.time()
    A <- coef_to_AB(model$beta_path[, 1], p)$A
    B <- coef_to_AB(model$beta_path[, 1], p)$B
    runtimeC <- Sys.time() - t0

    # Print how long it took to run formatted with a message
    message(paste0("ADMM took ", round(runtimeD, 2), " seconds to run for the calculation of D."))
    message(paste0("ADMM took ", round(runtimeM, 2), " seconds to run for the model."))
    message(paste0("ADMM took ", round(runtimeC, 2), " seconds to run for conversion to A, B."))

    # Return the fitted model
    output_list <- list(
        model = model,
        A = A,
        B = B,
        runtimeD = runtimeD,
        runtimeM = runtimeM,
        runtimeC = runtimeC
    )
    return(output_list)
}

regular_splash <- function(y, ...) {
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(dim(y)[1])

    # Fit SPLASH from Reuvers and Wijler (2021)
    t0 <- Sys.time()
    model <- splash::splash(t(y), ...) # Take transpose because splash() accepts T x p matrix
    runtime <- Sys.time() - t0

    # Print how long it took to run formatted with a message
    message(paste0("SPLASH took ", round(runtime, 2), " seconds to run."))

    # Return the fitted model
    output_list <- list(
        model = model,
        A = model$AB[, 1:p, ],
        B = model$AB[, (p + 1):(2 * p), ],
        runtime = runtime
    )
    return(output_list)
}
