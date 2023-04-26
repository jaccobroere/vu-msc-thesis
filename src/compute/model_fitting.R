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
sim_design_id <- "designA_T500_p50"
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
sigma_hat <- t(fread(path1, header = T, skip = 0))
Vhat_d <- as.matrix(fread(path2, header = T, skip = 0))
gr <- read_graph(path3, format = "graphml")
sym_gr <- read_graph(path4, format = "graphml")

# Load the true values
A <- as.matrix(fread(path5, header = T, skip = 0))
B <- as.matrix(fread(path6, header = T, skip = 0))
y <- as.matrix(fread(path7, header = T, skip = 0))

# Print the dimensions of the data
message(cat("The dimension of y: ", dim(y)[1], dim(y)[2]))

# Fit the a single solution using (Augmented) ADMM of GSPLASH
model_gsplash <- fit_admm_gsplash(sigma_hat, Vhat_d, gr, lambda1, lambda2, standard_ADMM = TRUE)

# Fit a single solution of symmetric_GSPLASH
model_sym_gsplash <- fit_admm_gsplash(sigma_hat, Vhat_d, sym_gr, lambda1, lambda2, standard_ADMM = TRUE)

# Fit the a single solution using SPLASH
model_splash <- fit_regular_splash(y, banded_covs = c(TRUE, TRUE), B = 500, alphas = c(0.5), lambdas = c(lambda_splash))

# Fit a single solution using PVAR(1) with the BigVAR package
model_pvar <- fit_pvar_bigvar(y, lambda_pvar)

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
