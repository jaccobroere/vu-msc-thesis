PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/R/utils.R")
source("src/compute/R/model_cv_wrappers.R")
library(data.table)
library(tictoc)
library(matrixStats)
library(Matrix)
library(tictoc)
library(splash)
setwd(PROJ_DIR)
################################################################################
# PATHING AND DATA LOADING
################################################################################
read_data_design <- function(sim_design_id) {
    # Set up directories
    data_dir <<- file.path(PROJ_DIR, "data/simulation/runtime", sim_design_id)
    fit_dir <<- file.path(PROJ_DIR, "out/simulation/", sim_design_id)

    # Parse paths
    path_reg_graph <<- file.path(data_dir, "reg_graph.graphml")
    path_sym_graph <<- file.path(data_dir, "sym_graph.graphml")
    path_y <<- file.path(data_dir, "y.csv")
    path_A <<- file.path(data_dir, "A.csv")
    path_B <<- file.path(data_dir, "B.csv")
    path_Dtilde_inv <<- file.path(data_dir, "Dtilde_inv.mtx")
    path_Dtilde <<- file.path(data_dir, "Dtilde.mtx")
    path_Dtilde_SSF_inv <<- file.path(data_dir, "Dtilde_SSF_inv.mtx")
    path_Dtilde_SSF <<- file.path(data_dir, "Dtilde_SSF.mtx")
    path_bandwidth <<- file.path(data_dir, "bootstrap_bandwidths.csv")

    # Load the data
    reg_gr <<- read_graph(path_reg_graph, format = "graphml")
    sym_gr <<- read_graph(path_sym_graph, format = "graphml")
    Dtilde_inv <<- readMM(path_Dtilde_inv)
    Dtilde <<- readMM(path_Dtilde)
    Dtilde_SSF_inv <<- readMM(path_Dtilde_SSF_inv)
    Dtilde_SSF <<- readMM(path_Dtilde_SSF)
    bandwidths <<- as.data.frame(fread(path_bandwidth, header = T, skip = 0))
    # Load the true values
    A_true <<- as.matrix(fread(path_A, header = T, skip = 0))
    B_true <<- as.matrix(fread(path_B, header = T, skip = 0))
    C_true <<- AB_to_C(A_true, B_true)
    y <<- as.matrix(fread(path_y, header = T, skip = 0))
    # Split the data into train and test (Test set is only one value as we are MC simulating one-step ahead prediction error)
    train_idx <<- (floor(dim(y)[2] / 5) * 4)
    y_train <<- y[, 1:train_idx]
    y_test <<- y[, (train_idx + 1):dim(y)[2]]
}

fit_cv_and_time <- function(res_tensor, design, i) {
    # Fit FSPLASH model and measure runtime
    t0 <- Sys.time()
    Dtilde_inv <- Matrix::solve(Dtilde, diag(dim(Dtilde)[1]), sparse = TRUE)
    model_fsplash <- fit_fsplash.cv(y, bandwidths, reg_gr, Dtilde_inv, nlambdas = 10, nfolds = 1)
    t1 <- Sys.time()
    runtime_fsplash <- difftime(t1, t0, units = "secs")[[1]]

    # Fit SSFSPLASH model and measure runtime
    t0 <- Sys.time()
    Dtilde_SSF_inv <- Matrix::solve(Dtilde_SSF, diag(dim(Dtilde_SSF)[1]), sparse = TRUE)
    model_ssfsplash <- fit_ssfsplash.cv(y, bandwidths, reg_gr, Dtilde_SSF_inv, alpha = 0.5, nlambdas = 10, nfolds = 1)
    t1 <- Sys.time()
    runtime_ssfsplash <- difftime(t1, t0, units = "secs")[[1]]

    # Fit SPLASH models and measure runtime
    t0 <- Sys.time()
    model_splash_a0 <- fit_splash.cv(y, alpha = 0, nlambdas = 10, nfolds = 1)
    t1 <- Sys.time()
    runtime_splash_a0 <- difftime(t1, t0, units = "secs")[[1]]

    # Fit GF-SPLASH models and measure runtime
    t0 <- Sys.time()
    model_gfsplash_sym_a0 <- fit_gfsplash.cv(y, bandwidths, graph = sym_gr, alpha = 0, nlambdas = 10, nfolds = 1)
    t1 <- Sys.time()
    runtime_gfsplash_sym_a0 <- difftime(t1, t0, units = "secs")[[1]]

    # PVAR
    t0 <- Sys.time()
    model_pvar <- fit_pvar.cv(y, nlambdas = 20, nfolds = 1)
    t1 <- Sys.time()
    runtime_pvar <- difftime(t1, t0, units = "secs")[[1]]

    # Create a data frame with the runtimes
    res_tensor[design, , i] <- c(runtime_fsplash, runtime_ssfsplash, runtime_splash_a0, runtime_gfsplash_sym_a0, runtime_pvar)

    # Print the runtimes
    print(paste("FSPLASH runtime:", runtime_fsplash))
    print(paste("SSFSPLASH runtime:", runtime_ssfsplash))
    print(paste("SPLASH a=0 runtime:", runtime_splash_a0))
    print(paste("GF-SPLASH sym a=0 runtime:", runtime_gfsplash_sym_a0))
    print(paste("PVAR runtime:", runtime_pvar))

    return(res_tensor)
}

# Data directory
root_data_dir <- file.path(PROJ_DIR, "data/simulation/runtime")

# Initialize the result tensor
designs <- list.dirs(root_data_dir, recursive = FALSE, full.names = FALSE)
models <- c("FSPLASH", "SSFSPLASH", "SPLASH_a0", "GF_SPLASH_sym_a0", "PVAR")
K <- 10
iterations <- seq(1:K)
res_tensor <- array(NA, dim = c(length(designs), length(models), K), dimnames = list(designs, models, iterations))

for (i in iterations) {
    for (design in designs) {
        print(paste("Design:", design, "Iteration:", i, "\n"), sep = "\n")
        read_data_design(design)
        res_tensor <- fit_cv_and_time(res_tensor, design, i)
    }
}

mean_res <- apply(res_tensor, c(1, 2), mean)
std_res <- apply(res_tensor, c(1, 2), sd)

# Save the results
# Save mean_res as CSV
dir.create(file.path("out", "simulation", "runtime"), recursive = TRUE)
write.csv(mean_res, file.path("out", "simulation", "runtime", "runtime_mean.csv"), row.names = TRUE)

# Save std_res as CSV
write.csv(std_res, file.path("out", "simulation", "runtime", "runtime_std.csv"), row.names = TRUE)
