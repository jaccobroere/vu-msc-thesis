# Load all the necessary packages
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
sim_design_id <- args[1]
uuidtag <- args[2]

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "detlam", uuidtag)
lambdas_dir <- file.path(PROJ_DIR, "out/simulation/lambdas", sim_design_id, uuidtag)

# Calculate lambda_0 for the GSPLASH
path_sigma_hat <- file.path(data_dir, "sigma_hat.csv")
path_Vhat_d <- file.path(data_dir, "Vhat_d.mtx")
path_reg_graph <- file.path(data_dir, "reg_graph.graphml")
path_sym_graph <- file.path(data_dir, "sym_graph.graphml")
path_y <- file.path(data_dir, "y.csv")
path_A <- file.path(data_dir, "A.csv")
path_B <- file.path(data_dir, "B.csv")


# Load the data
sigma_hat <- t(fread(path_sigma_hat, header = T, skip = 0))
Vhat_d <- as.matrix(readMM(path_Vhat_d))
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")
y <- as.matrix(fread(path_y, header = T, skip = 0))
A <- as.matrix(fread(path_A, header = T, skip = 0))
B <- as.matrix(fread(path_B, header = T, skip = 0))

# Calculate the C_true matrix for the RMSFE calcualtion
C_true <- AB_to_C(A, B)

# Calculate the RMSFE for each of the lambdas and add them to the .csv file to later aggregate
res_reg_a0 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, reg_gr, alpha = 0, path = file.path(lambdas_dir, "reg_a0.csv"))
res_reg_a05 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, reg_gr, alpha = 0.5, path = file.path(lambdas_dir, "reg_a05.csv"))
res_sym_a0 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, sym_gr, alpha = 0, path = file.path(lambdas_dir, "sym_a0.csv"))
res_sym_a05 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, sym_gr, alpha = 0.5, path = file.path(lambdas_dir, "sym_a05.csv"))
res_spl_a0 <- run_lambda_finder_splash(y, alpha = 0, C_true, path = file.path(lambdas_dir, "spl_a0.csv"))
res_spl_a05 <- run_lambda_finder_splash(y, alpha = 0.5, C_true, path = file.path(lambdas_dir, "spl_a05.csv"))
