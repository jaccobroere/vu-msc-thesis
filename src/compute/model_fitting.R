PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/utils.R")
source("src/compute/model_wrappers.R")
library(data.table)
library(tictoc)
library(matrixStats)

# Read CLI arguments
args <- commandArgs(trailingOnly = TRUE)

# Set up directories
data_dir <- paste0(PROJ_DIR, "/data/simulation/")
out_dir <- paste0(PROJ_DIR, "/out/simulation/coef/")

path_prefix <- "designB_T1000_p100"
# Read arguments from command line input
# path_prefix <- args[1]

# Parse paths
path1 <- paste0(data_dir, path_prefix, "_sigma_hat.csv")
path2 <- paste0(data_dir, path_prefix, "_Vhat_d.csv")
path3 <- paste0(data_dir, path_prefix, "_graph.graphml")
path4 <- paste0(data_dir, path_prefix, "_sym_graph.graphml")
path5 <- paste0(data_dir, path_prefix, "_A.csv")
path6 <- paste0(data_dir, path_prefix, "_B.csv")
path7 <- paste0(data_dir, path_prefix, "_y.csv")

# Load the data
sigma_hat <- t(fread(path1, header = T, skip = 0))
Vhat_d <- as.matrix(fread(path2, header = T, skip = 0))
gr <- read_graph(path3, format = "graphml")
sym_gr <- read_graph(path4, format = "graphml")

# Load the true values
A <- as.matrix(fread(path4, header = T, skip = 0))
B <- as.matrix(fread(path5, header = T, skip = 0))
y <- as.matrix(fread(path6, header = T, skip = 0))

# Print the dimensions of the data
message(cat("The dimension of y: ", dim(y)[1], dim(y)[2]))

# Set the regularization parameter
lambda1 <- 0.05 # Fused penalty
lambda2 <- 0.01 # Lasso penalty

lambda_splash <- 0.01
lambda_pvar <- 1

# Fit the a single solution using (Augmented) ADMM of GSPLASH
gsplash <- fit_admm_gsplash(sigma_hat, Vhat_d, gr, lambda1, lambda2, standard_ADMM = TRUE)

# Fit a single solution of symmetric_GSPLASH
sym_gsplash <- fit_sym_gsplash(sigma_hat, Vhat_d, gr, lambda1, lambda2, standard_ADMM = TRUE)

# Fit the a single solution using SPLASH
splash <- fit_regular_splash(y, banded_covs = c(TRUE, TRUE), B = 500, alphas = c(0.5), lambdas = c(lambda_splash))

# Fit a single solution using PVAR(1) with the BigVAR package
pvar <- fit_pvar_bigvar(y_train, lambda_pvar)

# Save the results
fwrite(data.table(gsplash$A), file = paste0(out_dir, path_prefix, "_admm_gsplash_estimate_A.csv"))
fwrite(data.table(gsplash$B), file = paste0(out_dir, path_prefix, "_admm_gsplash_estimate_B.csv"))
fwrite(data.table(splash$A), file = paste0(out_dir, path_prefix, "_splash_estimate_A.csv"))
fwrite(data.table(splash$B), file = paste0(out_dir, path_prefix, "_splash_estimate_B.csv"))


norm(gsplash$A - A, "2")
norm(splash$A - A, "2")

norm(gsplash$B - B, "2")
norm(splash$B - B, "2")
