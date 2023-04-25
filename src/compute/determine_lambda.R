setwd(system("echo $PROJ_DIR", intern = TRUE))
source("src/compute/utils.R")
setwd(system("echo $ZHU_DIR", intern = TRUE))
source("R/opt.R")
source("R/gen_data.R")
library(genlasso)
library(igraph)
library(splash)
library(FGSG)
library(BigVAR)
setwd(system("echo $PROJ_DIR", intern = TRUE))

# Set up directories
data_dir <- paste0(PROJ_DIR, "/data/simulation/")
out_dir <- paste0(PROJ_DIR, "/out/simulation/coef/")

# Calculate lambda_0 for the GSPLASH
path_sigma_hat <- paste0(data_dir, sim_id_dir, path_prefix, "_sigma_hat.csv")
path_Vhat_d <- paste0(data_dir, sim_id_dir, path_prefix, "_Vhat_d.csv")
path_reg_graph <- paste0(data_dir, sim_id_dir, path_prefix, "_graph.graphml")
path_sym_graph <- paste0(data_dir, sim_id_dir, path_prefix, "_sym_graph.graphml")

# Load the data
sigma_hat <- t(fread(path_sigma_hat, header = T, skip = 0))
Vhat_d <- as.matrix(fread(path_Vhat_d, header = T, skip = 0))
gr <- read_graph(path_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")

# Calculate lambda_0 for the GSPLASH
lam0_reg_a0 <- calc_lambda_0(sigma_hat, Vhat_d, gr, alpha = 0)
lam0_reg_a05 <- calc_lambda_0(sigma_hat, Vhat_d, gr, alpha = 0.5)
lam0_sym_a0 <- calc_lambda_0(sigma_hat, Vhat_d, sym_gr, alpha = 0)
lam0_sym_a05 <- calc_lambda_0(sigma_hat, Vhat_d, sym_gr, alpha = 0.5)

# Generate grid of values for lambda
