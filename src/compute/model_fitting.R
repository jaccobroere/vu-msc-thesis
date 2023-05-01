PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/utils.R")
source("src/compute/model_wrappers.R")
library(data.table)
library(tictoc)
library(matrixStats)

# Read CLI arguments
args <- commandArgs(trailingOnly = TRUE)
# sim_design_id <- "designB_T500_p25"
sim_design_id <- args[1]
uuidtag <- args[2]

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation/", sim_design_id, uuidtag)
fit_dir <- file.path(PROJ_DIR, "out/simulation/fit/", sim_design_id, uuidtag)

# Parse paths
path_sigma_hat <- paste0(data_dir, sim_design_id, "_sigma_hat.csv")
path_Vhat_d <- paste0(data_dir, sim_design_id, "_Vhat_d.csv")
path_reg_graph <- paste0(data_dir, sim_design_id, "_graph.graphml")
path_sym_graph <- paste0(data_dir, sim_design_id, "_sym_graph.graphml")
path_y <- paste0(data_dir, sim_design_id, "_y.csv")
path_A <- paste0(data_dir, sim_design_id, "_A.csv")
path_B <- paste0(data_dir, sim_design_id, "_B.csv")

# Load the data
sigma_hat <- t(fread(path_sigma_hat, header = T, skip = 0))
Vhat_d <- as.matrix(fread(path_Vhat_d, header = T, skip = 0))
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")
# Load the true values
A <- as.matrix(fread(path_A, header = T, skip = 0))
B <- as.matrix(fread(path_B, header = T, skip = 0))
y <- as.matrix(fread(path_y, header = T, skip = 0))

# Fit the a single solution using (Augmented) ADMM of GSPLASH
model_gsplash_a0 <- fit_admm_gsplash(sigma_hat, Vhat_d, reg_gr, lambda, alpha = 0, standard_ADMM = TRUE)
model_gsplash_a05 <- fit_admm_gsplash(sigma_hat, Vhat_d, reg_gr, lambda, alpha = 0.5, standard_ADMM = TRUE)
# Fit a single solution of symmetric_GSPLASH
model_sym_gsplash_a0 <- fit_admm_gsplash(sigma_hat, Vhat_d, sym_gr, lambda, alpha = 0, standard_ADMM = TRUE)
model_sym_gsplash_a05 <- fit_admm_gsplash(sigma_hat, Vhat_d, sym_gr, lambda, alpha = 0.5, standard_ADMM = TRUE)
# Fit the a single solution using SPLASH
model_splash_a0 <- fit_regular_splash(y, lambda, alpha = 0)
model_splash_a05 <- fit_regular_splash(y, lambda, alpha = 0.5)
# Fit a single solution using PVAR(1) with the BigVAR package, also runs a grid search
model_pvar <- fit_pvar_bigvar(y)

# Compute and save the predictions
C <- AB_to_C(A, B)
yhat_gsplash <- predict_with_C(model_gsplash$C, y)
yhat_splash <- predict_with_C(model_splash$C, y)
yhat_sym_gsplash <- predict_with_C(model_sym_gsplash$C, y)
yhat_pvar <- predict_with_C(model_pvar$C, y)
y_hat_true <- predict_with_C(C, y)

# Save the results
fwrite(data.table(model_gsplash$A), file = paste0(coef_dir, sim_design_id, "_gsplash_estimate_A.csv"))
fwrite(data.table(model_gsplash$B), file = paste0(coef_dir, sim_design_id, "_gsplash_estimate_B.csv"))
fwrite(data.table(model_splash$A), file = paste0(coef_dir, sim_design_id, "_splash_estimate_A.csv"))
fwrite(data.table(model_splash$B), file = paste0(coef_dir, sim_design_id, "_splash_estimate_B.csv"))
fwrite(data.table(model_sym_gsplash$A), file = paste0(coef_dir, sim_design_id, "_sym_gsplash_estimate_A.csv"))
fwrite(data.table(model_sym_gsplash$B), file = paste0(coef_dir, sim_design_id, "_sym_gsplash_estimate_B.csv"))
fwrite(data.table(model_pvar$A), file = paste0(coef_dir, sim_design_id, "_pvar_estimate_C.csv"))
fwrite(data.table(yhat_gsplash), file = paste0(coef_dir, sim_design_id, "_gsplash_estimate_yhat.csv"))
fwrite(data.table(yhat_splash), file = paste0(coef_dir, sim_design_id, "_splash_estimate_yhat.csv"))
fwrite(data.table(yhat_sym_gsplash), file = paste0(coef_dir, sim_design_id, "_sym_gsplash_estimate_yhat.csv"))
fwrite(data.table(yhat_pvar), file = paste0(coef_dir, sim_design_id, "_pvar_estimate_yhat.csv"))

# train_idx <- (floor(dim(y)[2] / 5) * 4)
# y_train <- y[, 1:train_idx]
# y_test <- y[, train_idx:ncol(y)]

# calc_rmsfe(y_test, yhat_gsplash, y_hat_true)
# calc_rmsfe(y_test, yhat_splash, y_hat_true)
# calc_rmsfe(y_test, yhat_sym_gsplash, y_hat_true)
# calc_rmsfe(y_test, yhat_pvar, y_hat_true)

# calc_msfe(y_test, yhat_gsplash)
# calc_msfe(y_test, yhat_splash)
# calc_msfe(y_test, yhat_sym_gsplash)
# calc_msfe(y_test, yhat_pvar)
