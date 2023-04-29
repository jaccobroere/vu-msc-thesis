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

fit_admm_gsplash <- function(sigma_hat, Vhat_d, graph, lambda, alpha, verbose = FALSE, ...) {
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
        Y = sigma_hat,
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
    AB <- coef_to_AB(model$beta_path[, 1], p)
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

fit_regular_splash <- function(y, lambda, alpha, verbose = FALSE, ...) {
    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(dim(y)[1])

    # Fit SPLASH from Reuvers and Wijler (2021)
    t0 <- Sys.time()
    model <- splash::splash(t(y_train), alphas = c(alpha), lambdas = c(lambda), ...) # Take transpose because splash() accepts T x p matrix
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

fit_pvar_bigvar <- function(y, verbose = FALSE, ...) {
    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, (floor(dim(y)[2] / 5) * 4 + 1):dim(y)[2]]

    # Fit a single solution using PVAR(1) with the BigVAR package
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Compute lambda grid
    lambda_grid <- gen_lambda_grid(100)

    # Perform cross-validation for the selection of the penalty parameter
    t0 <- Sys.time()
    model <- constructModel(
        Y = t(y),
        p = 1,
        struct = "Basic",
        gran = lambda_grid,
        loss = "L2",
        T1 = floor(dim(y_train)[2] / 5) * 4 + 1,
        T2 = dim(y_train)[2] + 1,
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
    model <- FGSG::gflasso(y = sigma_hat, A = Vhat_d, tp = edge_vector, s1 = lambda2, s2 = lambda1, ...) # Notice that \lambda1 and \lambda2 are swapped here
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
