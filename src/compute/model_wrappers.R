setwd("/Users/jacco/Documents/repos/vu-msc-thesis/admm_src_zhu")
source("R/opt.R")
source("R/gen_data.R")
library(genlasso)
library(igraph)
library(splash)

admm_gsplash <- function(sigma_hat, Vhat_d, graph, lambda, gamma, ...) {
    # Convert graph to sparse matrix
    D <- as(getDgSparse(graph = graph), "TsparseMatrix")
    idx <- D@i + 1
    jdx <- D@j + 1
    val <- D@x

    # Fit linear regression model using ADMM
    t0 <- Sys.time()
    admmmodel <- linreg_path_v2(
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
    runtime <- Sys.time() - t0

    # Print how long it took to run formatted with a message
    message(paste0("ADMM took ", round(runtime, 2), " seconds to run."))

    # Return the fitted model
    return(admmmodel)
}

regular_splash <- function(y, ...) {
    # Fit SPLASH from Reuvers and Wijler (2021)
    t0 <- Sys.time()
    splashmodel <- splash::splash(y, ...)
    runtime <- Sys.time() - t0

    # Print how long it took to run formatted with a message
    message(paste0("SPLASH took ", round(runtime, 2), " seconds to run."))

    # Return the fitted model
    return(splashmodel)
}
