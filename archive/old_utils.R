fit_fgsg_gsplash <- function(sigma_hat, Vhat_d, graph, lambda1, lambda2, ...) {
    # Retrieve the cross-sectional dimension of the problem
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    # Create edge vector
    edge_vector <- as.vector(t(as_edgelist(graph)))

    # Fit FGSG GFLASSO implementation
    t0 <- Sys.time()
    model <- FGSG::gflasso(y = sigma_hat, A = Vhat_d, tp = edge_vector, s1 = lambda2, s2 = lambda1, ...) # Notice that \lambda1 and \lambda2 are swapped here
    runtime <- difftime(Sys.time(), t0, units = "secs")[[1]]

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
        C = AB_to_C(A, B),
        runtime = runtime
    )
}
