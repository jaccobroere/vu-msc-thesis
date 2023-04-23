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
setwd(system("echo $PROJ_DIR", intern = TRUE))

fit_admm_gsplash <- function(sigma_hat, Vhat_d, graph, lambda1, labmda2, ...) {
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Translate lambdas to gamma
    lambda <- lambda1
    gamma <- lambda2 / lambda1

    # Convert graph to sparse matrix
    t0 <- Sys.time()
    D <- as(getDgSparse(graph = graph), "TsparseMatrix")
    idx <- D@i + 1
    jdx <- D@j + 1
    val <- D@x
    runtimeD <- difftime(Sys.time() - t0, units = "secs")[[1]]

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
    AB <- coef_to_AB(model$beta_path[, 1], p)
    A <- AB$A
    B <- AB$B
    runtimeC <- difftime(Sys.time(), t0, units = "secs")[[1]]

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

fit_regular_splash <- function(y, ...) {
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(dim(y)[1])

    # Fit SPLASH from Reuvers and Wijler (2021)
    t0 <- Sys.time()
    model <- splash::splash(t(y), ...) # Take transpose because splash() accepts T x p matrix
    runtime <- difftime(Sys.time() - t0, units = "secs")[[1]]

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

fit_fgsg_gsplash <- function(sigma_hat, Vhat_d, graph, lambda1, lambda2, ...) {
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Create edge vector
    edge_vector <- as.vector(t(as_edgelist(graph)))

    # Fit FGSG GFLASSO implementation
    t0 <- Sys.time()
    model <- FGSG::gflasso(y = sigma_hat, A = Vhat_d, tp = edge_vector, s1 = lambda2, s2 = lambda1, ...) # Notice that \lambda1 and \lambda2 are swapped here
    runtime <- difftime(Sys.time() - t0, units = "secs")[[1]]

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
        runtime = runtime
    )
}

fit_pvar_bigvar <- function(y, lambda,  ...) {
    # Split y into training and testing sets
    y_train = y[, 1:floor(dim(y)[2] / 5) * 4]
    y_test = y[, (floor(dim(y)[2] / 5) * 4 + 1):dim(y)[2]]

    # Fit a single solution using PVAR(1) with the BigVAR package
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Perform cross-validation for the selection of the penalty parameter
    t0 <- Sys.time()
    model <- constructModel(
        Y = t(y_train),
        p = 1,
        struct = "Basic",
        gran = c(100, 10, 1, 0.1, 0.01),
        loss = "L1",
        T1 = floor(dim(y_train)[2] / 5) * 4 + 1,
        ownlambdas = TRUE,
        model.controls = list(
            intercept = FALSE,
            loss = "L1"
        )
    )
    cvmodel <- cv.BigVAR(model)
    runtimeCV <- difftime(Sys.time() - t0, units = "secs")[[1]]

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
    runtime <- difftime(Sys.time() - t0, units = "secs")[[1]]

    # Print how long it took to run formatted with a message
    message(paste0("PVAR took ", round(runtime, 2), " seconds to run for the model."))

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

