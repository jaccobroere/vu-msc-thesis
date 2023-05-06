# This script fits the models using the best lambdas for the GF-SPLASH, SGF-SPLASH, and SPLASH


PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/utils.R")
source("src/compute/model_wrappers.R")
library(data.table)
library(tictoc)
library(matrixStats)
library(Matrix)
setwd(PROJ_DIR)

# Read CLI arguments
args <- commandArgs(trailingOnly = TRUE)
sim_design_id <- ifelse(length(args) < 1, "designB_T500_p49", args[1])
uuidtag <- ifelse(length(args) < 2, "ED04541C-0149-45C6-8F59-6D5A21BFB0C2", args[2])

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "mc", uuidtag)
fit_dir <- file.path(PROJ_DIR, "out/simulation/fit", sim_design_id, uuidtag)
lambdas_dir <- file.path(PROJ_DIR, "out/simulation/lambdas", sim_design_id)

# Parse paths
path_sigma_hat <- file.path(data_dir, "sigma_hat.csv")
path_Vhat_d <- file.path(data_dir, "Vhat_d.mtx")
path_reg_graph <- file.path(data_dir, "reg_graph.graphml")
path_sym_graph <- file.path(data_dir, "sym_graph.graphml")
path_y <- file.path(data_dir, "y.csv")
path_A <- file.path(data_dir, "A.csv")
path_B <- file.path(data_dir, "B.csv")

# Load the data
sigma_hat <- t(fread(path_sigma_hat, header = T, skip = 0))
Vhat_d <- as.matrix(readMM(path_Vhat_d))
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")
# Load the true values
A_true <- as.matrix(fread(path_A, header = T, skip = 0))
B_true <- as.matrix(fread(path_B, header = T, skip = 0))
y <- as.matrix(fread(path_y, header = T, skip = 0))

# Read best_lambdas
best_lam_df <- as.data.frame(fread(file.path(lambdas_dir, "best_lambdas.csv"), header = T, skip = 0))

# Fit the a single solution using (Augmented) ADMM of GSPLASH
model_gsplash_a0 <- fit_gfsplash(sigma_hat, Vhat_d, reg_gr, get_lam_best(best_lam_df, "best_lam_reg_a0"), alpha = 0, standard_ADMM = TRUE)
model_gsplash_a05 <- fit_gfsplash(sigma_hat, Vhat_d, reg_gr, get_lam_best(best_lam_df, "best_lam_reg_a05"), alpha = 0.5, standard_ADMM = TRUE)
# Fit a single solution of symmetric_GSPLASH
model_sym_gsplash_a0 <- fit_gfsplash(sigma_hat, Vhat_d, sym_gr, get_lam_best(best_lam_df, "best_lam_sym_a0"), alpha = 0, standard_ADMM = TRUE)
model_sym_gsplash_a05 <- fit_gfsplash(sigma_hat, Vhat_d, sym_gr, get_lam_best(best_lam_df, "best_lam_sym_a05"), alpha = 0.5, standard_ADMM = TRUE)
# Fit the a single solution using SPLASH
model_splash_a0 <- fit_splash(y, get_lam_best(best_lam_df, "best_lam_spl_a0"), alpha = 0)
model_splash_a05 <- fit_splash(y, get_lam_best(best_lam_df, "best_lam_spl_a05"), alpha = 0.5)
# Fit a single solution using PVAR(1) with the BigVAR package, also runs a grid search
model_pvar <- fit_pvar(y)

# Model fast fusion
model_fast_fusion <- fit_fsplash(sigma_hat, Vhat_d, reg_gr, lambda = get_lam_best(best_lam_df, "best_lam_reg_a0"))

# Compute and save the predictions
C_true <- AB_to_C(A_true, B_true)
model_gsplash_a0$yhat <- predict_with_C(model_gsplash_a0$C, y)
model_gsplash_a05$yhat <- predict_with_C(model_gsplash_a05$C, y)
model_splash_a0$yhat <- predict_with_C(model_splash_a0$C, y)
model_splash_a05$yhat <- predict_with_C(model_splash_a05$C, y)
model_sym_gsplash_a0$yhat <- predict_with_C(model_sym_gsplash_a0$C, y)
model_sym_gsplash_a05$yhat <- predict_with_C(model_sym_gsplash_a05$C, y)
model_pvar$yhat <- predict_with_C(model_pvar$C, y)
y_hat_true <- predict_with_C(C_true, y)
model_fast_fusion$yhat <- predict_with_C(model_fast_fusion$C, y)

# Save the results
save_fitting_results(model_gsplash_a0, "reg_a0", fit_dir)
save_fitting_results(model_gsplash_a05, "reg_a05", fit_dir)
save_fitting_results(model_splash_a0, "spl_a0", fit_dir)
save_fitting_results(model_splash_a05, "spl_a05", fit_dir)
save_fitting_results(model_sym_gsplash_a0, "sym_a0", fit_dir)
save_fitting_results(model_sym_gsplash_a05, "sym_a05", fit_dir)
save_fitting_results(model_pvar, "pvar", fit_dir)

# Save true values as well ,y_hat_true are the next time step predictions based on the true C matrix
fwrite(data.table(y_hat_true), file = file.path(fit_dir, "y_hat_true.csv"))
fwrite(data.table(A_true), file = file.path(fit_dir, "A_true.csv"))
fwrite(data.table(B_true), file = file.path(fit_dir, "B_true.csv"))


model_fast_fusion$runtimeM
model_fast_fusion$runtimeXtilde

train_idx <- (floor(dim(y)[2] / 5) * 4)
y_train <- y[, 1:train_idx]
y_test <- y[, (train_idx + 1):ncol(y)]

calc_msfe(y_test, model_gsplash_a0$yhat)
calc_msfe(y_test, model_gsplash_a05$yhat)
calc_msfe(y_test, model_splash_a0$yhat)
calc_msfe(y_test, model_splash_a05$yhat)
calc_msfe(y_test, model_sym_gsplash_a0$yhat)
calc_msfe(y_test, model_sym_gsplash_a05$yhat)
calc_msfe(y_test, model_pvar$yhat)
calc_msfe(y_test, y_hat_true)
calc_msfe(y_test, model_fast_fusion$yhat)


model_gsplash_a0$model$beta[1:10]
model_fast_fusion$coef[1:10]



### DEVELOPING THE FAST FUSION METHOD

class(Vhat_d)

Vhat_d_sparse <- Matrix(Vhat_d, sparse = TRUE)
class(Vhat_d_sparse)

library(Rlinsolve)
library(Matrix)

t(Vhat_d_sparse)
Dtilde_sparse <- calc_Dtilde_sparse(reg_gr)

# Trying solvers
res1 <- lsolve.bicg(t(Dtilde_sparse), t(Vhat_d_sparse))

res2 <- lsolve.bicgstab(t(Dtilde_sparse), t(Vhat_d_sparse))



fit_faster_fusion <- function(sigma_hat, Vhat_d, graph, lambda) {
    Vhat_d <- Matrix(Vhat_d, sparse = TRUE)
    t0 <- Sys.time()
    Dtilde <- calc_Dtilde_sparse(graph)
    runtimeDprime <- difftime(Sys.time(), t0, units = "secs")[[1]]

    m <- ecount(graph)
    p <- as.integer(sqrt(dim(Vhat_d)[1]))

    t0 <- Sys.time()
    # Use linear system solvers for faster computation of the change of variables (see Tibshirani and Taylor, 2011)
    XD1 <- t(lsolve.sor(t(Dtilde), t(Vhat_d))$x) # Same as Vhat_d * inv(Dtilde)
    X1 <- XD1[, 1:m]
    X2 <- XD1[, (m + 1):dim(XD1)[2]]
    X2_plus <- lsolve.sor((t(X2) %*% X2), t(X2))$x # Same as inv(t(X2) %*% X2) %*% t(X2)

    # Transform the input to LASSO objective
    P <- X2 %*% X2_plus
    ytilde <- as.vector((diag(nrow(P)) - P) %*% t(sigma_hat))
    Xtilde <- (diag(nrow(P)) - P) %*% X1
    runtimeXtilde_fast <- difftime(Sys.time(), t0, units = "secs")[[1]]


    t0 <- Sys.time()
    # Fit the LASSO model and back-transform the coefficients
    model <- glmnet(Xtilde, ytilde, lambda = lambda, alpha = 1, intercept = FALSE, standardize = FALSE)
    theta1 <- as.vector(model$beta)
    theta2 <- as.vector(X2_plus %*% (t(sigma_hat) - X1 %*% theta1))
    coef <- lsolve.sor(Dtilde, c(theta1, theta2))$x
    runtimeM <- difftime(Sys.time(), t0, units = "secs")[[1]]

    AB <- coef_to_AB(coef, p)
    A <- AB$A
    B <- AB$B

    # Return the fitted model
    output_list <- list(
        model = model,
        coef = coef,
        A = A,
        B = B,
        C = AB_to_C(A, B),
        runtimeM = runtimeM,
        runtimeDprime = runtimeDprime,
        runtimeXtilde = runtimeXtilde,
        runtimeXtilde_fast = runtimeXtilde_fast
    )
}
