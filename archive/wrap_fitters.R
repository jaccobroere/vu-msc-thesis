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
