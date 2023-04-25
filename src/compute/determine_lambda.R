PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/utils.R")
source("src/compute/model_wrappers.R")
setwd(system("echo $ZHU_DIR", intern = TRUE))
source("R/opt.R")
source("R/gen_data.R")
library(genlasso)
library(data.table)
library(igraph)
library(splash)
library(FGSG)
library(BigVAR)
setwd(PROJ_DIR)

# Read arguments from CLI
args <- commandArgs(trailingOnly = TRUE)
path_prefix <- "designA_T500_p25" # path_prefix <- args[1]
sim_id <- "1343" # sim_id_dir <- args[2]

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation/")
out_dir <- file.path(PROJ_DIR, "out/")
lambdas_dir <- file.path(PROJ_DIR, "out/simulation/lambdas/")

# Create directory for the simulation ID
sim_id_dir <- file.path(lambdas_dir, sim_id)
dir.create(sim_id_dir)

# Calculate lambda_0 for the GSPLASH
path_sigma_hat <- paste0(data_dir, path_prefix, "_sigma_hat.csv")
path_Vhat_d <- paste0(data_dir, path_prefix, "_Vhat_d.csv")
path_reg_graph <- paste0(data_dir, path_prefix, "_graph.graphml")
path_sym_graph <- paste0(data_dir, path_prefix, "_sym_graph.graphml")
path_y <- paste0(data_dir, path_prefix, "_y.csv")

# Load the data
sigma_hat <- t(fread(path_sigma_hat, header = T, skip = 0))
Vhat_d <- as.matrix(fread(path_Vhat_d, header = T, skip = 0))
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")
y <- as.matrix(fread(path_y, header = T, skip = 0))

# Calculate lambda_0 for the GSPLASH
# lam0_reg_a0 <- calc_lambda_0(sigma_hat, Vhat_d, reg_gr, alpha = 0)
# lam0_reg_a05 <- calc_lambda_0(sigma_hat, Vhat_d, reg_gr, alpha = 0.5)
# lam0_sym_a0 <- calc_lambda_0(sigma_hat, Vhat_d, sym_gr, alpha = 0)
# lam0_sym_a05 <- calc_lambda_0(sigma_hat, Vhat_d, sym_gr, alpha = 0.5)

# # Generate grid of values for lambda
# grid_lam_reg_a0 <- gen_lambda_grid(lam0_reg_a0)
# grid_lam_reg_a05 <- gen_lambda_grid(lam0_reg_a05)
# grid_lam_sym_a0 <- gen_lambda_grid(lam0_sym_a0)
# grid_lam_sym_a05 <- gen_lambda_grid(lam0_sym_a05)

# # Call the function for each specification and create file if it does not exist yet
# df_lam_reg_a0 <- create_lambda_df(grid_lam_reg_a0, file.path(sim_id_dir, "reg_a0.csv"))
# df_lam_reg_a05 <- create_lambda_df(grid_lam_reg_a05, file.path(sim_id_dir, "reg_a05.csv"))
# df_lam_sym_a0 <- create_lambda_df(grid_lam_sym_a0, file.path(sim_id_dir, "sym_a0.csv"))
# df_lam_sym_a05 <- create_lambda_df(grid_lam_sym_a05, file.path(sim_id_dir, "sym_a05.csv"))


res_reg_a0 <- run_lambda_finder(sigma_hat, Vhat_d, reg_gr, alpha = 0, path = file.path(sim_id_dir, "reg_a0.csv"))
res_reg_a05 <- run_lambda_finder(sigma_hat, Vhat_d, reg_gr, alpha = 0.5, path = file.path(sim_id_dir, "reg_a05.csv"))
res_sym_a0 <- run_lambda_finder(sigma_hat, Vhat_d, sym_gr, alpha = 0, path = file.path(sim_id_dir, "sym_a0.csv"))
res_sym_a05 <- run_lambda_finder(sigma_hat, Vhat_d, sym_gr, alpha = 0.5, path = file.path(sim_id_dir, "sym_a05.csv"))
