# This script fits the models using the best lambdas for the GF-SPLASH, SGF-SPLASH, and SPLASH
PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/utils.R")
source("src/compute/model_wrappers.R")
source("src/compute/model_cv_wrappers.R")
library(data.table)
library(tictoc)
library(matrixStats)
library(Matrix)
library(tictoc)
setwd(PROJ_DIR)

# Read CLI arguments
args <- commandArgs(trailingOnly = TRUE)
sim_design_id <- ifelse(length(args) < 1, "designB_T1000_p25", args[1])
uuidtag <- ifelse(length(args) < 2, "D67DFC71-FF34-4731-B009-6C9668E5DA4E", args[2])

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "mc", uuidtag)
fit_dir <- file.path(PROJ_DIR, "out/simulation/fit", sim_design_id, uuidtag)
lambdas_dir <- file.path(PROJ_DIR, "out/simulation/lambdas", sim_design_id)
sim_id_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "mc")

# Parse paths
path_sigma_hat <- file.path(data_dir, "sigma_hat.csv")
path_Vhat_d <- file.path(data_dir, "Vhat_d.mtx")
path_reg_graph <- file.path(data_dir, "reg_graph.graphml")
path_sym_graph <- file.path(data_dir, "sym_graph.graphml")
path_y <- file.path(data_dir, "y.csv")
path_A <- file.path(data_dir, "A.csv")
path_B <- file.path(data_dir, "B.csv")
path_Dtilde_inv <- file.path(sim_id_dir, "Dtilde_inv.mtx")
path_Dtilde <- file.path(sim_id_dir, "Dtilde.mtx")
path_Dtilde_SSF_inv <- file.path(sim_id_dir, "Dtilde_SSF_inv.mtx")
path_Dtilde_SSF <- file.path(sim_id_dir, "Dtilde_SSF.mtx")
path_bandwidth <- file.path(data_dir, "bootstrap_bandwidths.csv")

# Load the data
sigma_hat <- t(fread(path_sigma_hat, header = T, skip = 0))
Vhat_d <- as.matrix(readMM(path_Vhat_d))
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")
Dtilde_inv <- readMM(path_Dtilde_inv)
Dtilde <- readMM(path_Dtilde)
Dtilde_SSF_inv <- readMM(path_Dtilde_SSF_inv)
Dtilde_SSF <- readMM(path_Dtilde_SSF)
bandwidths <- as.data.frame(fread(path_bandwidth, header = T, skip = 0))
# Load the true values
A_true <- as.matrix(fread(path_A, header = T, skip = 0))
B_true <- as.matrix(fread(path_B, header = T, skip = 0))
C_true <- AB_to_C(A_true, B_true)
y <- as.matrix(fread(path_y, header = T, skip = 0))


################################################################################
# MODEL FITTING
################################################################################
# Read best_lambdas
best_lam_df <- as.data.frame(fread(file.path(lambdas_dir, "best_lambdas.csv"), header = T, skip = 0))

# Fit the a single solution using (Augmented) ADMM of GSPLASH
model_gsplash_a05 <- fit_gfsplash(sigma_hat, Vhat_d, reg_gr, get_lam_best(best_lam_df, "best_lam_reg_a05"), alpha = 0.5, standard_ADMM = TRUE)

# Fit a single solution of symmetric_GSPLASH
model_sym_gsplash_a0 <- fit_gfsplash(sigma_hat, Vhat_d, sym_gr, get_lam_best(best_lam_df, "best_lam_sym_a0"), alpha = 0, standard_ADMM = TRUE)
model_sym_gsplash_a05 <- fit_gfsplash(sigma_hat, Vhat_d, sym_gr, get_lam_best(best_lam_df, "best_lam_sym_a05"), alpha = 0.5, standard_ADMM = TRUE)

# Fit the a single solution using SPLASH
model_splash_a0 <- fit_splash(y, get_lam_best(best_lam_df, "best_lam_spl_a0"), alpha = 0)
model_splash_a05 <- fit_splash(y, get_lam_best(best_lam_df, "best_lam_spl_a05"), alpha = 0.5)

# Fit a single solution for F-SPLASH and SSF-SPLASH
model_fsplash <- fit_fsplash(sigma_hat, Vhat_d, reg_gr, lambda = get_lam_best(best_lam_df, "best_lam_fsp"))
model_ssfsplash <- fit_ssfsplash(sigma_hat, Vhat_d, Dtilde_SSF_inv, lambda = get_lam_best(best_lam_df, "best_lam_ssf"), alpha = 0.5)

# Fit a single solution using PVAR(1) with the BigVAR package, also runs a grid search
model_pvar <- fit_pvar(y)

# Compute and save the predictions
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


A_true[1:10, 1:10]
model_splash_a05$model$lambda

model_fsplash <- fit_fsplash(sigma_hat, Vhat_d, reg_gr, Dtilde_inv, lambda = c(0.1, 0.2))

model_fsplash$model$lambda


mmm <- glmnet(Vhat_d, as.vector(sigma_hat), alpha = 1)
mmm$lambda

model_pvar$cvmodel


# Testing cv functions
source("src/compute/model_cv_wrappers.R")
cv_splash_a0 <- fit_splash.cv(y, alpha = 0, nfolds = 5, nlambdas = 20)
cv_fsplash <- fit_fsplash.cv(y, bandwidths, reg_gr, Dtilde_inv, nfolds = 5, nlambdas = 20)
cv_ssfsplash <- fit_ssfsplash.cv(y, bandwidths, reg_gr, Dtilde_SSF_inv, alpha = 0.5, nfolds = 5, nlambdas = 20)
cv_pvar <- fit_pvar.cv(y, nfolds = 5, nlambdas = 20)
cv_splash_a0$msfe
cv_fsplash$msfe
cv_ssfsplash$msfe
cv_pvar$msfe

cv_pvar
dim(cv_pvar$errors_cv)

alpha <- 0
nfolds <- 5
nlambdas <- 20
