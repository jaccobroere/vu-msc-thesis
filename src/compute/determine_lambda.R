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

# Define the functions to run the lambda finder
run_lambda_finder_fsplash <- function(y, sigma_hat, Vhat_d, C_true, graph, Dtilde_inv, path) {
    train_idx <- (floor(dim(y)[2] / 5) * 4)
    y_train <- y[, 1:train_idx]
    y_test <- y[, (train_idx + 1):ncol(y)]
    # Fit the models and save the prediction results in the data.frame
    model <- fit_fsplash(sigma_hat, Vhat_d, graph, Dtilde_inv, lambda = NULL, nlambda = 20)
    df_lam <- create_lambda_df(model$model$lambda, path)
    # Calculate predictions and error metric
    for (i in 1:length(grid_lam)) {
        lam <- model$model$lambda[i]
        y_hat <- predict_with_C(res$C[, , i], y)
        error_metric <- calc_msfe(y_test, y_hat)
        df_lam[1, colnames(df_lam)[i]] <- error_metric
    }

    # Append the results to the table
    # write.table(df_lam, path, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
}

# Read CLI arguments
args <- commandArgs(trailingOnly = TRUE)
sim_design_id <- ifelse(length(args) < 1, "designB_T500_p25", args[1])
uuidtag <- ifelse(length(args) < 2, "503FD119-142F-4B9B-AD46-CA0A417B03E6", args[2])

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "detlam", uuidtag)
data_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "mc", uuidtag)
lambdas_dir <- file.path(PROJ_DIR, "out/simulation/lambdas", sim_design_id, uuidtag)
sim_id_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "mc")

# Parse paths
path_sigma_hat <- file.path(data_dir, "sigma_hat.csv")
path_Vhat_d <- file.path(data_dir, "Vhat_d.mtx")
path_reg_graph <- file.path(data_dir, "reg_graph.graphml")
path_sym_graph <- file.path(data_dir, "sym_graph.graphml")
path_y <- file.path(data_dir, "y.csv")
path_A <- file.path(data_dir, "A.csv")
path_B <- file.path(data_dir, "B.csv")
path_Dtilde_inv <- file.path(sim_id_dir, "Dtilde_inv.mtx")
path_Dtilde <- file.path(sim_id_dir, "Dtilde.mtx")

# Load the data
sigma_hat <- t(fread(path_sigma_hat, header = T, skip = 0))
Vhat_d <- as.matrix(readMM(path_Vhat_d))
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")
Dtilde_inv <- readMM(path_Dtilde_inv)
Dtilde <- readMM(path_Dtilde)
A_true <- as.matrix(fread(path_A, header = T, skip = 0))
B_true <- as.matrix(fread(path_B, header = T, skip = 0))
y <- as.matrix(fread(path_y, header = T, skip = 0))

# Calculate the C_true matrix for the RMSFE calcualtion
C_true <- AB_to_C(A, B)

# Calculate the RMSFE for each of the lambdas and add them to the .csv file to later aggregate
res_reg_a0 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, reg_gr, alpha = 0, path = file.path(lambdas_dir, "reg_a0.csv"))
res_reg_a05 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, reg_gr, alpha = 0.5, path = file.path(lambdas_dir, "reg_a05.csv"))
res_sym_a0 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, sym_gr, alpha = 0, path = file.path(lambdas_dir, "sym_a0.csv"))
res_sym_a05 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, sym_gr, alpha = 0.5, path = file.path(lambdas_dir, "sym_a05.csv"))
res_spl_a0 <- run_lambda_finder_splash(y, alpha = 0, C_true, path = file.path(lambdas_dir, "spl_a0.csv"))
res_spl_a05 <- run_lambda_finder_splash(y, alpha = 0.5, C_true, path = file.path(lambdas_dir, "spl_a05.csv"))



runner <- run_lambda_finder_fsplash(y, sigma_hat, Vhat_d, C_true, reg_gr, Dtilde_inv, path = file.path(lambdas_dir, "fsplash.csv"))


dim(Vhat_d)
dim(Dtilde_inv)
