wrap_ssfsplash <- function() {
    tic()
    model <- fit_ssfsplash(sigma_hat, Vhat_d, Dtilde_SSF_inv, lambda = 0.05, alpha = 0)
    toc()
    print(calc_rmsfe(y_test, predict_with_C(model$C, y), predict_with_C(C_true, y)))
    print(calc_EE(A_true, model$A, type = "2"))
    print(calc_EE(B_true, model$B, type = "2"))
    return(model)
}

wrap_gfsplash <- function() {
    tic()
    model <- fit_gfsplash(sigma_hat, Vhat_d, reg_gr, alpha = 0, lambda = 0.05, scale = TRUE)
    toc()
    print(calc_rmsfe(y_test, predict_with_C(model$C, y), predict_with_C(C_true, y)))
    print(calc_EE(A_true, model$A, type = "2"))
    print(calc_EE(B_true, model$B, type = "2"))
    return(model)
}

wrap_fsplash <- function() {
    tic()
    model <- fit_fsplash(sigma_hat, Vhat_d, reg_gr, Dtilde_inv, lambda = 0.05)
    toc()
    print(calc_rmsfe(y_test, predict_with_C(model$C, y), predict_with_C(C_true, y)))
    print(calc_EE(A_true, model$A, type = "2"))
    print(calc_EE(B_true, model$B, type = "2"))
    return(model)
}

wrap_splash <- function(lambda = 0.1) {
    tic()
    model <- fit_splash(y, lambda = lambda, alpha = 0)
    toc()
    print(calc_rmsfe(y_test, predict_with_C(model$C, y), predict_with_C(C_true, y)))
    print(calc_EE(A_true, model$A, type = "2"))
    print(calc_EE(B_true, model$B, type = "2"))
    return(model)
}

wrap_pvar <- function() {
    tic()
    model <- fit_pvar(y)
    toc()
    print(calc_rmsfe(y_test, predict_with_C(model$C, y), predict_with_C(C_true, y)))
    print(calc_EE(C_true, model$C, type = "2"))
    return(model)
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
