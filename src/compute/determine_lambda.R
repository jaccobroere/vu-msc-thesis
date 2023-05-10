# Load all the necessary packages
PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/utils.R")
source("src/compute/model_wrappers.R")
source("src/compute/model_cv_wrappers.R")
setwd(system("echo $ZHU_DIR", intern = TRUE))
source("R/opt.R")
source("R/gen_data.R")
library(genlasso)
library(data.table)
library(igraph)
library(splash)
library(BigVAR)
setwd(PROJ_DIR)

# Read CLI arguments
args <- commandArgs(trailingOnly = TRUE)
sim_design_id <- ifelse(length(args) < 1, "designB_T500_p25", args[1])
uuidtag <- ifelse(length(args) < 2, "503FD119-142F-4B9B-AD46-CA0A417B03E6", args[2])

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "detlam", uuidtag)
lambdas_dir <- file.path(PROJ_DIR, "out/simulation/lambdas", sim_design_id, uuidtag)
sim_id_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id)

# Parse paths
path_reg_graph <- file.path(data_dir, "reg_graph.graphml")
path_sym_graph <- file.path(data_dir, "sym_graph.graphml")
path_y <- file.path(data_dir, "y.csv")
path_A <- file.path(data_dir, "A.csv")
path_B <- file.path(data_dir, "B.csv")
path_Dtilde_inv <- file.path(sim_id_dir, "Dtilde_inv.mtx")
path_Dtilde <- file.path(sim_id_dir, "Dtilde.mtx")
path_bandwidth <- file.path(sim_id_dir, "bootstrap_bandwidths.csv")

# Load the data
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")
Dtilde_inv <- readMM(path_Dtilde_inv)
Dtilde <- readMM(path_Dtilde)
A_true <- as.matrix(fread(path_A, header = T, skip = 0))
B_true <- as.matrix(fread(path_B, header = T, skip = 0))
C_true <- AB_to_C(A, B)
y <- as.matrix(fread(path_y, header = T, skip = 0))

# Calculate the RMSFE for each of the lambdas and add them to the .csv file to later aggregate
res_reg_a05 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, reg_gr, alpha = 0.5, path = file.path(lambdas_dir, "reg_a05.csv"))
res_sym_a0 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, sym_gr, alpha = 0, path = file.path(lambdas_dir, "sym_a0.csv"))
res_sym_a05 <- run_lambda_finder_gfsplash(y, sigma_hat, Vhat_d, C_true, sym_gr, alpha = 0.5, path = file.path(lambdas_dir, "sym_a05.csv"))


#
fit_gfsplash.lambda_grid <- function(y, bandwidths, graph, alpha, path, nlambdas = 20, lambda.min.ratio = 1e-4) {
    # Call the function for each specification and create file if it does not exist yet
    df_lam <- data.frame(matrix(ncol = nlambdas, nrow = 0))
    colnames(df_lam) <- paste0("ilam", 1:nlambdas)
    if (!file.exists(path)) {
        write.table(df_lam, path, row.names = FALSE, col.names = TRUE, append = FALSE, sep = ",")
    }
    # Get the errors per lambda index from the model
    model_grid <- fit_gfsplash.cv.iter(y, bandwidths, graph, alpha = alpha, nlambdas = nlambdas, lambda.min.ratio = lambda.min.ratio)
    errors <- data.frame(model_grid$errors)
    # Append the results to the table
    write.table(df_lam, path, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    return(df_lam)
}


grid_reg_a05 <- fit_gfsplash.lambda_grid(y, bandwidths, reg_gr, alpha = 0.5, path = "tryout_df_lam.csv")

graph <- reg_gr
alpha <- 0.5
nlambdas <- 20
lambda.min.ratio <- 1e-4


log10(lambda.min.ratio)


model_grid$beta_path[, 10]

best_idx
best_lambda
best_msfe

dim(y)
train_idx
