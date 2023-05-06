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
sim_design_id <- ifelse(length(args) < 1, "designB_T500_p9", args[1])
uuidtag <- ifelse(length(args) < 2, "5c9c2f78-3338-4b4a-8121-2adc01e8990d", args[2])

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id)
fit_dir <- file.path(PROJ_DIR, "out/simulation/fit", sim_design_id)
lambdas_dir <- file.path(PROJ_DIR, "out/simulation/lambdas", sim_design_id)

# Parse paths
path_sigma_hat <- file.path(data_dir, "sigma_hat.csv")
path_Vhat_d <- file.path(data_dir, "Vhat_d.mtx")
path_reg_graph <- file.path(data_dir, "reg_graph.graphml")
path_sym_graph <- file.path(data_dir, "sym_graph.graphml")
path_y <- file.path(data_dir, "y.csv")
path_A <- file.path(data_dir, "A.csv")
path_B <- file.path(data_dir, "B.csv")
path_H <- file.path(data_dir, "H.csv")
path_z <- file.path(data_dir, "z.csv")

# Load the data
sigma_hat <- t(fread(path_sigma_hat, header = T, skip = 0))
Vhat_d <- as.matrix(readMM(path_Vhat_d))
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")
# Load the true values
A_true <- as.matrix(fread(path_A, header = T, skip = 0))
B_true <- as.matrix(fread(path_B, header = T, skip = 0))
H <- as.matrix(fread(path_H, header = T, skip = 0))
z <- as.matrix(fread(path_z, header = T, skip = 0))
y <- as.matrix(fread(path_y, header = T, skip = 0))

# Set the regularization parameter
lambda <- 0.005050447
alpha <- 0

# Fit the a single solution using (Augmented) ADMM
gsplash <- fit_admm_gsplash(sigma_hat, Vhat_d, reg_gr, lambda, alpha)
flasso <- fusedlasso(sigma_hat, Vhat_d, D = D, gamma = 0)

flasso$beta[, 50]
flasso$lambda[50]
as.vector(gsplash$model$beta_path)

# Fit the a single solution using SPLASH
splash <- regular_splash(y, banded_covs = c(FALSE, FALSE), B = 500, alphas = c(0.5), lambdas = c(lambda))

# Save the results
fwrite(data.table(gsplash$A), file = paste0(out_dir, file_prefix, "_admm_gsplash_estimate_A.csv"))
fwrite(data.table(gsplash$B), file = paste0(out_dir, file_prefix, "_admm_gsplash_estimate_B.csv"))
fwrite(data.table(splash$A), file = paste0(out_dir, file_prefix, "_splash_estimate_A.csv"))
fwrite(data.table(splash$B), file = paste0(out_dir, file_prefix, "_splash_estimate_B.csv"))


# Trying regular splash with transformed data
library(pracma)
n <- nrow(Vhat_d)
p <- ncol(Vhat_d)
D <- getDg(reg_gr)
Phi <- rbind(Vhat_d, D)
QR <- qr(Phi)
Q <- as.matrix(qr.Q(QR))
R <- as.matrix(qr.R(QR))
Q1 <- Q[1:n, ]
Q2 <- Q[(n + 1):nrow(Q), ]

pQ2 <- t(Q2) %*% solve(Q2 %*% t(Q2))
H <- Q1 %*% pQ2
z <- as.vector((diag(n) - Q1 %*% (diag(p) - pQ2 %*% Q2) %*% t(Q1)) %*% t(sigma_hat))


path <- glmnet(H, z, lambda = lambda, alpha = 1, intercept = FALSE, standardize = FALSE)
theta <- path$beta

beta <- solve(R) %*% (diag(p) - pQ2 %*% Q2) %*% t(Q1) %*% t(sigma_hat) + solve(R) %*% pQ2 %*% beta
as.vector(beta)

model <- fit_fast_fusion(sigma_hat, Vhat_d, reg_gr, lambda)


model
model$coef
