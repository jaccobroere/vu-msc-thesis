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

################################################################################
# PATHING AND DATA LOADING
################################################################################
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
# Split the data into train and test (Test set is only one value as we are MC simulating one-step ahead prediction error)
train_idx <- train_idx <- (floor(dim(y)[2] / 5) * 4)
y_train <- y[, 1:train_idx]
y_test <- y[, (train_idx + 1):dim(y)[2]]

################################################################################
# Load the best lambdas for GF-SPLASH models
################################################################################
grid_reg_a05 <- as.data.frame(fread(file.path(lambdas_dir, "grid_reg_a05.csv"), header = T, skip = 0))
grid_sym_a0 <- as.data.frame(fread(file.path(lambdas_dir, "grid_sym_a0.csv"), header = T, skip = 0))
grid_sym_a05 <- as.data.frame(fread(file.path(lambdas_dir, "grid_sym_a05.csv"), header = T, skip = 0))

get_best_lam_idx <- function(grid) {
    # Get the best lambda index
    best_lam_idx <- which.min(colMeans(grid))
    return(best_lam_idx)
}

lam_idx_reg_a05 <- get_best_lam_idx(grid_reg_a05)
lam_idx_sym_a0 <- get_best_lam_idx(grid_sym_a0)
lam_idx_sym_a05 <- get_best_lam_idx(grid_sym_a05)

################################################################################
# MODEL FITTING
################################################################################
# Fit F-SPLASH and SSF-SPLASH models (Uses CV at each iteration)
model_fsplash <- fit_fsplash.cv(y, bandwidths, reg_gr, Dtilde_inv, nlambdas = 20, nfolds = 5)
model_ssfsplash <- fit_ssfsplash.cv(y, bandwidths, reg_gr, Dtilde_SSF_inv, alpha = 0.5, nlambdas = 20, nfolds = 5)

# Fit SPLASH model (Uses CV at each iteration)
model_splash_a0 <- fit_splash.cv(y, alpha = 0, nlambdas = 20, nfolds = 5)
model_splash_a05 <- fit_splash.cv(y, alpha = 0.5, nlambdas = 20, nfolds = 5)

# Fit GF-SPLASH models (Preliminary CV was employed to find optimal lambda_index)
model_gfsplash_reg_a05 <- fit_gfsplash.on_idx(y, bandwidths, lam_idx_reg_a05, alpha = 0.5, graph = reg_gr, nlambdas = 20, nfolds = 5)
model_gfsplash_sym_a0 <- fit_gfsplash.on_idx(y, bandwidths, lam_idx_sym_a0, alpha = 0, graph = sym_gr, nlambdas = 20, nfolds = 5)
model_gfsplash_sym_a05 <- fit_gfsplash.on_idx(y, bandwidths, lam_idx_sym_a05, alpha = 0.5, graph = sym_gr, nlambdas = 20, nfolds = 5)

# Fit PVAR (Uses CV at each iteration)
model_pvar <- fit_pvar.cv(y, nlambdas = 20, nfolds = 5)

# Compute the predictions stemming from the true value of C
y_hat_true <- predict_with_C(C_true, y_train, y_test)

################################################################################
# SAVING RESULTS
################################################################################
save_fitting_results <- function(model, prefix, fit_dir) {
    rmsfe <- calc_rmsfe(y_test[, 1], model$y_pred[, 1], y_hat_true[, 1]) # One step ahead RMSFE
    fwrite(data.table(model$A), file = file.path(fit_dir, paste0(prefix, "_estimate_A.csv")))
    fwrite(data.table(model$B), file = file.path(fit_dir, paste0(prefix, "_estimate_B.csv")))
    fwrite(data.table(model$C), file = file.path(fit_dir, paste0(prefix, "_estimate_C.csv")))
    fwrite(data.table(model$y_pred), file = file.path(fit_dir, paste0(prefix, "_y_pred.csv")))
    fwrite(data.table(rmsfe), file = file.path(fit_dir, paste0(prefix, "_rmsfe.csv")))
}

save_fitting_results(model_gfsplash_reg_a05, "reg_a05", fit_dir)
save_fitting_results(model_gfsplash_sym_a0, "sym_a0", fit_dir)
save_fitting_results(model_gfsplash_sym_a05, "sym_a05", fit_dir)
save_fitting_results(model_fsplash, "fsplash", fit_dir)
save_fitting_results(model_ssfsplash, "ssfsplash", fit_dir)
save_fitting_results(model_splash_a0, "splash_a0", fit_dir)
save_fitting_results(model_splash_a05, "splash_a05", fit_dir)
save_fitting_results(model_pvar, "pvar", fit_dir)

# Save true values as well ,y_hat_true are the next time step predictions based on the true C matrix
fwrite(data.table(y_hat_true), file = file.path(fit_dir, "y_hat_true.csv"))
fwrite(data.table(A_true), file = file.path(fit_dir, "A_true.csv"))
fwrite(data.table(B_true), file = file.path(fit_dir, "B_true.csv"))
