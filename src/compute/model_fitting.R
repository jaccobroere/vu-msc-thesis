PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/utils.R")
source("src/compute/model_wrappers.R")
library(data.table)
library(tictoc)
library(matrixStats)

# Set up directories
data_dir <- paste0(PROJ_DIR, "/data/simulation/")
coef_dir <- paste0(PROJ_DIR, "/out/simulation/coef/")

# Read CLI arguments
args <- commandArgs(trailingOnly = TRUE)
sim_design_id <- "designB_T500_p25"
# sim_design_id <- args[1]

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

# Print the dimensions of the data
message(cat("The dimension of y: ", dim(y)[1], dim(y)[2]))


#
alpha <- 0.5
lambda <- 0.01

# Fit the a single solution using (Augmented) ADMM of GSPLASH
model_gsplash <- fit_admm_gsplash(sigma_hat, Vhat_d, reg_gr, lambda, alpha, standard_ADMM = TRUE)

# Fit a single solution of symmetric_GSPLASH
model_sym_gsplash <- fit_admm_gsplash(sigma_hat, Vhat_d, sym_gr, lambda, alpha, standard_ADMM = TRUE)

# Fit the a single solution using SPLASH
model_splash <- fit_regular_splash(y, lambda, alpha)

# Fit a single solution using PVAR(1) with the BigVAR package
model_pvar <- fit_pvar_bigvar(y, lambda)

# Compute and save the predictions
yhat_gsplash <- predict_with_C(model_gsplash$C, y)
yhat_splash <- predict_with_C(model_splash$C, y)
yhat_sym_gsplash <- predict_with_C(model_sym_gsplash$C, y)
yhat_pvar <- predict_with_C(model_pvar$C, y)

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


spl <- splash::splash(t(y), banded_covs = c(TRUE, TRUE), B = 500, n_lambdas = 5, alphas = c(0.5), lambda_min_mult = 1e-4)


run_lambda_finder_gfsplash(sigma_hat, Vhat_d, reg_gr, 0.5, "gsplash_runner.csv")
run_lambda_finder_splash(y, 0.5, "splash_runner.csv")


spl$AB
