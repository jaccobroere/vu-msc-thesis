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
run_lambda_finder_fsplash(y, sigma_hat, Vhat_d, C_ture, graph, Dtilde_inv, path) {
    train_idx <- (floor(dim(y)[2] / 5) * 4)
    y_train <- y[, 1:train_idx]
    y_test <- y[, (train_idx + 1):ncol(y)] 
    # Calculate lambda_0 for the GSPLASH
    df_lam <- create_lambda_df(grid_lam, path)
    # Fit the models and save the prediction results in the data.frame
    model <- fit_fsplash(sigma_hat, Vhat_d, graph, Dtilde_inv, lambda=NULL)
    # Calculate predictions and error metric
    for (i in 1:length(grid_lam)) {
        lam <- model$model$lambda[i]
        y_hat <- predict_with_C(res$C[, , i], y)
        error_metric <- calc_msfe(y_test, y_hat)
        df_lam[1, i] <- error_metric
    }

    # Append the results to the table
    write.table(df_lam, path, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
}

run_lambda_finder_gfsplash <- function(y, sigma_hat, Vhat_d, C_true, graph, alpha, path) {
    train_idx <- (floor(dim(y)[2] / 5) * 4)
    y_train <- y[, 1:train_idx]
    y_test <- y[, (train_idx + 1):ncol(y)]
    # Calculate lambda_0 for the GSPLASH
    lam0 <- calc_lambda_0_gfsplash(sigma_hat, Vhat_d, graph, alpha = alpha)
    # Generate grid of values for lambda
    grid_lam <- gen_lambda_grid(lam0, length.out = 20)
    # Call the function for each specification and create file if it does not exist yet
    df_lam <- create_lambda_df(grid_lam, path)
    # Fit the models and save the prediction results in the data.frame
    for (i in 1:length(grid_lam)) {
        lam <- grid_lam[i]
        model <- fit_admm_gsplash(sigma_hat, Vhat_d, graph, lam, alpha = alpha)
        # Calculate predictions and error metric
        y_hat <- predict_with_C(model$C, y)
        y_true_pred <- predict_with_C(C_true, y)
        error_metric <- calc_msfe(y_test, y_hat)
        df_lam[1, colnames(df_lam)[i]] <- error_metric
    }

    # Append the results to the table
    write.table(df_lam, path, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    return(df_lam)
}

run_lambda_finder_splash <- function(y, alpha, C_true, path, lambda_min_mult = 1e-4) {
    train_idx <- (floor(dim(y)[2] / 5) * 4)
    y_train <- y[, 1:train_idx]
    y_test <- y[, (train_idx + 1):ncol(y)]
    # Split the training set off
    y_train <- y[, 1:(floor(dim(y)[2] / 5) * 4)]
    # Run the SPLASH model with the given lambda_grid
    # model <- splash::splash(t(y_train), banded_covs = c(TRUE, TRUE), B = 500, n_lambdas = 20, alphas = c(alpha), lambda_min_mult = lambda_min_mult)
    lambda_grid <- gen_lambda_grid(lambda_min_mult, length.out = 20)
    model <- splash::splash(t(y_train), banded_covs = c(TRUE, TRUE), B = 500, lambdas = lambda_grid, alphas = c(alpha))
    p <- dim(model$AB)[1]
    # Create placeholder dataframe
    df_lam <- create_lambda_df(c(model$lambdas), path)
    # Generate and save predictions for each lambda value
    for (i in 1:dim(model$AB)[3]) {
        A <- model$AB[, 1:p, i]
        B <- model$AB[, (p + 1):(2 * p), i]
        C <- AB_to_C(A, B)
        # Calculate predictions and error metric
        y_hat <- predict_with_C(C, y)
        y_true_pred <- predict_with_C(C_true, y)
        error_metric <- calc_msfe(y_test, y_hat)
        df_lam[1, colnames(df_lam)[i]] <- error_metric
    }

    # Append the results to the table
    write.table(df_lam, path, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    return(df_lam)
}


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
