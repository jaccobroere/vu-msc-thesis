PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/R/utils.R")
source("src/compute/R/model_cv_wrappers.R")
library(data.table)
library(tictoc)
library(matrixStats)
library(Matrix)
library(tictoc)
library(splash)
setwd(PROJ_DIR)
################################################################################
# PATHING AND DATA LOADING
################################################################################
read_data_design <- function(sim_design_id) {
    # Set up directories
    data_dir <<- file.path(PROJ_DIR, "data/simulation/runtime", sim_design_id)
    fit_dir <<- file.path(PROJ_DIR, "out/simulation/", sim_design_id)

    # Parse paths
    path_reg_graph <<- file.path(data_dir, "reg_graph.graphml")
    path_sym_graph <<- file.path(data_dir, "sym_graph.graphml")
    path_y <<- file.path(data_dir, "y.csv")
    path_A <<- file.path(data_dir, "A.csv")
    path_B <<- file.path(data_dir, "B.csv")
    path_Dtilde_inv <<- file.path(data_dir, "Dtilde_inv.mtx")
    path_Dtilde <<- file.path(data_dir, "Dtilde.mtx")
    path_Dtilde_SSF_inv <<- file.path(data_dir, "Dtilde_SSF_inv.mtx")
    path_Dtilde_SSF <<- file.path(data_dir, "Dtilde_SSF.mtx")
    path_bandwidth <<- file.path(data_dir, "bootstrap_bandwidths.csv")

    # Load the data
    reg_gr <<- read_graph(path_reg_graph, format = "graphml")
    sym_gr <<- read_graph(path_sym_graph, format = "graphml")
    Dtilde_inv <<- readMM(path_Dtilde_inv)
    Dtilde <<- readMM(path_Dtilde)
    Dtilde_SSF_inv <<- readMM(path_Dtilde_SSF_inv)
    Dtilde_SSF <<- readMM(path_Dtilde_SSF)
    bandwidths <<- as.data.frame(fread(path_bandwidth, header = T, skip = 0))
    # Load the true values
    A_true <<- as.matrix(fread(path_A, header = T, skip = 0))
    B_true <<- as.matrix(fread(path_B, header = T, skip = 0))
    C_true <<- AB_to_C(A_true, B_true)
    y <<- as.matrix(fread(path_y, header = T, skip = 0))
    # Split the data into train and test (Test set is only one value as we are MC simulating one-step ahead prediction error)
    train_idx <<- (floor(dim(y)[2] / 5) * 4)
    y_train <<- y[, 1:train_idx]
    y_test <<- y[, (train_idx + 1):dim(y)[2]]
}

fit_cv_and_time <- function(res_tensor, design, i) {
    # Fit FSPLASH model and measure runtime
    t0 <- Sys.time()
    Dtilde_inv <- Matrix::solve(Dtilde, diag(dim(Dtilde)[1]), sparse = TRUE)
    model_fsplash <- fit_fsplash.cv(y, bandwidths, reg_gr, Dtilde_inv, nlambdas = 10, nfolds = 1)
    t1 <- Sys.time()
    runtime_fsplash <- difftime(t1, t0, units = "secs")[[1]]

    # Fit SSFSPLASH model and measure runtime
    t0 <- Sys.time()
    Dtilde_SSF_inv <- Matrix::solve(Dtilde_SSF, diag(dim(Dtilde_SSF)[1]), sparse = TRUE)
    model_ssfsplash <- fit_ssfsplash.cv(y, bandwidths, reg_gr, Dtilde_SSF_inv, alpha = 0.5, nlambdas = 10, nfolds = 1)
    t1 <- Sys.time()
    runtime_ssfsplash <- difftime(t1, t0, units = "secs")[[1]]

    # Fit SPLASH models and measure runtime
    t0 <- Sys.time()
    model_splash_a0 <- fit_splash.cv(y, alpha = 0, nlambdas = 10, nfolds = 1)
    t1 <- Sys.time()
    runtime_splash_a0 <- difftime(t1, t0, units = "secs")[[1]]

    # Fit GF-SPLASH models and measure runtime
    t0 <- Sys.time()
    model_gfsplash_sym_a0 <- fit_gfsplash.cv(y, bandwidths, graph = sym_gr, alpha = 0, nlambdas = 10, nfolds = 1)
    t1 <- Sys.time()
    runtime_gfsplash_sym_a0 <- difftime(t1, t0, units = "secs")[[1]]

    # PVAR
    t0 <- Sys.time()
    model_pvar <- fit_pvar.cv(y, nlambdas = 10, nfolds = 1)
    t1 <- Sys.time()
    runtime_pvar <- difftime(t1, t0, units = "secs")[[1]]

    # Create a data frame with the runtimes
    res_tensor[design, , i] <- c(runtime_fsplash, runtime_ssfsplash, runtime_splash_a0, runtime_gfsplash_sym_a0, runtime_pvar)

    # Print the runtimes
    print(paste("FSPLASH runtime:", runtime_fsplash))
    print(paste("SSFSPLASH runtime:", runtime_ssfsplash))
    print(paste("SPLASH a=0 runtime:", runtime_splash_a0))
    print(paste("GF-SPLASH sym a=0 runtime:", runtime_gfsplash_sym_a0))
    print(paste("PVAR runtime:", runtime_pvar))

    return(res_tensor)
}

# Data directory
root_data_dir <- file.path(PROJ_DIR, "data/simulation/runtime")

# Initialize the result tensor
designs <- sort(list.dirs(root_data_dir, recursive = FALSE, full.names = FALSE))
models <- c("FSPLASH", "SSFSPLASH", "SPLASH_a0", "GF_SPLASH_sym_a0", "PVAR")
K <- 10
iterations <- seq(1:K)
res_tensor <- array(NA, dim = c(length(designs), length(models), K), dimnames = list(designs, models, iterations))

for (i in iterations) {
    for (design in designs) {
        tic()
        print(paste("Design:", design, "Iteration:", i))
        read_data_design(design)
        res_tensor <- fit_cv_and_time(res_tensor, design, i)
        toc()
    }
}

mean_res <- apply(res_tensor, c(1, 2), mean)
std_res <- apply(res_tensor, c(1, 2), sd)

# Save the results
# Save mean_res as CSV
dir.create(file.path("out", "simulation", "runtime"), recursive = TRUE)
write.csv(mean_res, file.path("out", "simulation", "runtime", "runtime_mean.csv"), row.names = TRUE)

# Save std_res as CSV
write.csv(std_res, file.path("out", "simulation", "runtime", "runtime_std.csv"), row.names = TRUE)


# FASTER MATRIX MULTIPLICATION
fit_fsplash.cv <- function(y, bandwidths, graph, Dtilde_inv, nlambdas = 20, nfolds = 5, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]
    m <- ecount(graph)
    k <- vcount(graph)

    # Parse bandwidths
    h0 <- bandwidths[1]
    h1 <- bandwidths[2]

    # Split y into training and testing sets
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]

    # Create cross-validation folds
    folds <- rolling_cv(y_train, nfolds = nfolds)

    # Fit the model on each fold and save the results
    C_cv <- array(NA, dim = c(p, p, nlambdas))
    A_cv <- array(NA, dim = c(p, p, nlambdas))
    B_cv <- array(NA, dim = c(p, p, nlambdas))
    y_pred_cv <- array(NA, dim = c(p, length(folds[[1]]$test), nlambdas))
    errors_cv <- array(NA, dim = c(nfolds, nlambdas))
    for (i in 1:nfolds) {
        # Split the data into training and validation sets
        y_train_cv <- y_train[, folds[[i]]$train]
        y_val_cv <- y_train[, folds[[i]]$test]

        # Construct Vhat_d and sigma_hat from the training validation
        Sigma0 <- calc_Sigma_j(y_train_cv, 0)
        Sigma1 <- calc_Sigma_j(y_train_cv, 1)
        Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
        sigma_hat <- construct_sigma_hat(Sigma1, h1)
        Vhat_d <- construct_Vhat_d(Vhat)

        # Transform Generalized Lasso problem into LASSO problem
        XD1 <- Vhat_d %*% Dtilde_inv # Same as Vhat_d * inv(Dtilde)
        X1 <- XD1[, 1:m]
        X2 <- XD1[, (m + 1):dim(XD1)[2]]
        X2_plus <- solve((t(X2) %*% X2), t(X2)) # Same as inv(t(X2) %*% X2) %*% t(X2)

        # Transform the input to LASSO objective
        IminP <- as.matrix(diag(nrow(X2)) - as.matrix(X2) %*% as.matrix(X2_plus)) # Calculate I - P directly, P = X2 %*% X2_plus
        ytilde <- as.vector(IminP %*% sigma_hat)
        Xtilde <- IminP %*% as.matrix(X1)
        # Fit the model on the training set
        model_cv <- glmnet(Xtilde, ytilde, nlambda = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4, ...)

        # Compute the prediction error on the validation set
        for (j in 1:dim(model_cv$beta)[2]) {
            theta1 <- as.vector(model_cv$beta[, j])
            theta2 <- as.vector(X2_plus %*% (sigma_hat - X1 %*% theta1))
            coef_cv <- Dtilde_inv %*% c(theta1, theta2)
            AB <- coef_to_AB(coef_cv, p)
            A_cv[, , j] <- AB$A
            B_cv[, , j] <- AB$B
            C_cv[, , j] <- AB_to_C(AB$A, AB$B)
            y_pred_cv[, , j] <- predict_with_C(C_cv[, , j], y_train_cv, y_val_cv)
            errors_cv[i, j] <- calc_msfe(y_pred_cv[, , j], y_val_cv)
        }
    }

    # Get the best lambda value
    best_idx <- which.min(colMeans(errors_cv))

    t0 <- Sys.time()
    # Construct Vhat_d and sigma_hat from the training validation
    Sigma0 <- calc_Sigma_j(y_train, 0)
    Sigma1 <- calc_Sigma_j(y_train, 1)
    Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
    sigma_hat <- construct_sigma_hat(Sigma1, h1)
    Vhat_d <- construct_Vhat_d(Vhat)

    XD1 <- Vhat_d %*% Dtilde_inv # Same as Vhat_d * inv(Dtilde)
    X1 <- XD1[, 1:m]
    X2 <- XD1[, (m + 1):dim(XD1)[2]]
    X2_plus <- solve((t(X2) %*% X2), t(X2)) # Same as inv(t(X2) %*% X2) %*% t(X2)

    # Transform the input to LASSO objective
    IminP <- as.matrix(diag(nrow(X2)) - X2 %*% X2_plus) # Calculate I - P directly, P = X2 %*% X2_plus
    ytilde <- as.vector(IminP %*% sigma_hat)
    Xtilde <- IminP %*% as.matrix(X1)

    # Fit the model on the training set
    model <- glmnet(Xtilde, ytilde, nlambda = nlambdas, alpha = 1, intercept = FALSE, lambda.min.ratio = 1e-4, ...)
    theta1 <- as.vector(model$beta[, best_idx])
    theta2 <- as.vector(X2_plus %*% (sigma_hat - X1 %*% theta1))
    coef <- Dtilde_inv %*% c(theta1, theta2)
    t1 <- Sys.time()

    # Extract the model output
    AB <- coef_to_AB(coef, p)
    A <- AB$A
    B <- AB$B
    C <- AB_to_C(A, B)
    y_pred <- predict_with_C(C, y_train, y_test)

    output_list <- list(
        model = model,
        A = A,
        B = B,
        C = C,
        errors_cv = errors_cv,
        y_pred = y_pred,
        msfe = calc_msfe(y_test, y_pred),
        best_lambda = model$lambda[best_idx],
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
    return(output_list)
}



read_data_design("sim_runtime_p56_T1000")
graph <- reg_gr
nlambdas <- 20
alpha <- 0.5
nfolds <- 1

# Read problem dimensionality
p <- dim(y)[1]
m <- ecount(graph)
k <- vcount(graph)

# Parse bandwidths
h0 <- bandwidths[1]
h1 <- bandwidths[2]

# Split y into training and testing sets
y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]
Sigma0 <- calc_Sigma_j(y_train, 0)
Sigma1 <- calc_Sigma_j(y_train, 1)
Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
sigma_hat <- construct_sigma_hat(Sigma1, h1)
Vhat_d <- construct_Vhat_d(Vhat)

XD1 <- Vhat_d %*% Dtilde_inv # Same as Vhat_d * inv(Dtilde)
X1 <- XD1[, 1:m]
X2 <- XD1[, (m + 1):dim(XD1)[2]]
X2_plus <- solve((t(X2) %*% X2), t(X2)) # Same as inv(t(X2) %*% X2) %*% t(X2)

# Transform the input to LASSO objective
IminP <- as.matrix(diag(nrow(X2)) - X2 %*% X2_plus) # Calculate I - P directly, P = X2 %*% X2_plus
ytilde <- as.vector(IminP %*% sigma_hat)
Xtilde <- IminP %*% as.matrix(X1)

library(microbenchmark)

microbenchmark(IminP <- as.matrix(diag(nrow(X2)) - X2 %*% X2_plus), times = 3L)
microbenchmark(IminP <- as.matrix(diag(nrow(X2)) - fastMatMult(as.matrix(X2), as.matrix(X2_plus))), times = 3L)

microbenchmark(Xtilde <- IminP %*% as.matrix(X1), times = 3L)
microbenchmark(Xtilde <- fastMatMult(as.matrix(IminP), as.matrix(X1)), times = 3L)


microbenchmark(fastMatMult(as.matrix(Vhat_d), as.matrix(Dtilde_inv)), times = 3L)

microbenchmark(
    MMs(IminP, X1),
    times = 3L
)
MM(IminP, X1)

library(Rcpp)
# include Rcpp and Armadillo libraries
cppFunction(depends = "RcppArmadillo", code = "

// Function for normal matrix multiplication
arma::mat MM(const arma::mat & A, const arma::mat & B) {
  arma::mat C = A * B;
  return C;
}

// Function for sparse matrix multiplication
arma::sp_mat MMs(const arma::sp_mat & A, const arma::sp_mat & B) {
  arma::sp_mat C = A * B;
  return C;
}

")
