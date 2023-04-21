PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/utils.R")
source("src/compute/model_wrappers.R")
library(data.table)
library(tictoc)

# Read CLI arguments
args <- commandArgs(trailingOnly = TRUE)

# Set up directories
data_dir <- paste0(PROJ_DIR, "/data/simulation/")
out_dir <- paste0(PROJ_DIR, "/out/simulation/coef/")

path_prefix <- "designA_T500_p25"
# Read arguments from command line input
# path_prefix <- args[1]

# Parse paths
path1 <- paste0(data_dir, path_prefix, "_sigma_hat.csv")
path2 <- paste0(data_dir, path_prefix, "_Vhat_d.csv")
path3 <- paste0(data_dir, path_prefix, "_graph.graphml")
path4 <- paste0(data_dir, path_prefix, "_A.csv")
path5 <- paste0(data_dir, path_prefix, "_B.csv")
path6 <- paste0(data_dir, path_prefix, "_y.csv")

# Load the data
sigma_hat <- t(fread(path1, header = T, skip = 0))
Vhat_d <- as.matrix(fread(path2, header = T, skip = 0))
gr <- read_graph(path3, format = "graphml")
edge_vector <- as.vector(t(as_edgelist(gr)))

# Load the true values
A <- as.matrix(fread(path4, header = T, skip = 0))
B <- as.matrix(fread(path5, header = T, skip = 0))
y <- as.matrix(fread(path6, header = T, skip = 0))

# Print the dimensions of the data
message(cat("The dimension of y: ", dim(y)[1], dim(y)[2]))

# Set the regularization parameter
lambda1 <- 0.01
lambda2 <- 0

lambda_splash <- 0.05

# Fit the a single solution using (Augmented) ADMM
gsplash <- admm_gsplash(sigma_hat, Vhat_d, gr, lambda1, lambda2, standard_ADMM = TRUE)

# Fit a single solution using FGSG
gsplash_2 <- fgsg_gsplash(sigma_hat, Vhat_d, gr, lambda1, lambda2)

# Fit the a single solution using SPLASH
splash <- regular_splash(y, banded_covs = c(FALSE, FALSE), B = 500, alphas = c(0.5), lambdas = c(lambda_splash))

# Fit a single solution using PVAR(1)
# pvar <- pvar1(y, lambda)

# Fit a single solution using GMWK TODO

# Fit a single solution using GMWK(k_0) TODO


# Save the results
fwrite(data.table(gsplash$A), file = paste0(out_dir, path_prefix, "_admm_gsplash_estimate_A.csv"))
fwrite(data.table(gsplash$B), file = paste0(out_dir, path_prefix, "_admm_gsplash_estimate_B.csv"))
fwrite(data.table(splash$A), file = paste0(out_dir, path_prefix, "_splash_estimate_A.csv"))
fwrite(data.table(splash$B), file = paste0(out_dir, path_prefix, "_splash_estimate_B.csv"))


gsplash$A
gsplash$B

data.frame(gsplash_2$A)
