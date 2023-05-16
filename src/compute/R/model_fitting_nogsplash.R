# This script fits the models using the best lambdas for the GF-SPLASH, SGF-SPLASH, and SPLASH
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
# Read CLI arguments
args <- commandArgs(trailingOnly = TRUE)
sim_design_id <- ifelse(length(args) < 1, "designB_T500_p16", args[1])
uuidtag <- ifelse(length(args) < 2, "f8dac225-066f-4fc4-95e8-b3a838dcc96a", args[2])

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "mc", uuidtag)
fit_dir <- file.path(PROJ_DIR, "out/simulation/fit", sim_design_id, uuidtag)
lambdas_dir <- file.path(PROJ_DIR, "out/simulation/lambdas", sim_design_id)
sim_id_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "mc")

# Parse paths
path_reg_graph <- file.path(sim_id_dir, "reg_graph.graphml")
path_sym_graph <- file.path(sim_id_dir, "sym_graph.graphml")
path_y <- file.path(data_dir, "y.csv")
path_A <- file.path(data_dir, "A.csv")
path_B <- file.path(data_dir, "B.csv")
path_Dtilde_inv <- file.path(sim_id_dir, "Dtilde_inv.mtx")
path_Dtilde <- file.path(sim_id_dir, "Dtilde.mtx")
path_Dtilde_SSF_inv <- file.path(sim_id_dir, "Dtilde_SSF_inv.mtx")
path_Dtilde_SSF <- file.path(sim_id_dir, "Dtilde_SSF.mtx")
path_bandwidth <- file.path(sim_id_dir, "bootstrap_bandwidths.csv")

# Load the data
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")
Dtilde_inv <- readMM(path_Dtilde_inv)
Dtilde <- readMM(path_Dtilde)
Dtilde_SSF_inv <- readMM(path_Dtilde_SSF_inv)
Dtilde_SSF <- readMM(path_Dtilde_SSF)
bandwidths <- as.data.frame(fread(path_bandwidth, header = T, skip = 0))
# Load the true values
A_true <- as.matrix(fread(path_A, header = T, skip = 0))
B_true <- as.matrix(fread(path_B, header = T, skip = 0))
C_true <- AB_to_C(A_true, B_true)
y <- as.matrix(fread(path_y, header = T, skip = 0))
# Split the data into train and test (Test set is only one value as we are MC simulating one-step ahead prediction error)
train_idx <- (floor(dim(y)[2] / 5) * 4)
y_train <- y[, 1:train_idx]
y_test <- y[, (train_idx + 1):dim(y)[2]]

################################################################################
# MODEL FITTING
################################################################################
# Fit F-SPLASH and SSF-SPLASH models (Uses CV at each iteration)
model_fsplash <- fit_fsplash.cv(y, bandwidths, reg_gr, Dtilde_inv, nlambdas = 20, nfolds = 5)
model_ssfsplash <- fit_ssfsplash.cv(y, bandwidths, reg_gr, Dtilde_SSF_inv, alpha = 0.5, nlambdas = 20, nfolds = 5)

# Fit SPLASH model (Uses CV at each iteration)
model_splash_a0 <- fit_splash.cv(y, alpha = 0, nlambdas = 20, nfolds = 5)
model_splash_a05 <- fit_splash.cv(y, alpha = 0.5, nlambdas = 20, nfolds = 5)

# Fit PVAR (Uses CV at each iteration)
model_pvar <- fit_pvar.cv(y, nlambdas = 20, nfolds = 5)

# Compute the predictions stemming from the true value of C
y_hat_true <- predict_with_C(C_true, y_train, y_test)

################################################################################
# SAVING RESULTS
################################################################################
save_fitting_results <- function(model, prefix, fit_dir, save_AB = TRUE) {
    rmsfe <- calc_rmsfe(y_test[, 1], model$y_pred[, 1], y_hat_true[, 1]) # One step ahead RMSFE
    if (save_AB) {
        fwrite(data.table(model$A), file = file.path(fit_dir, paste0(prefix, "__estimate_A.csv")))
        fwrite(data.table(model$B), file = file.path(fit_dir, paste0(prefix, "__estimate_B.csv")))
    }
    fwrite(data.table(model$C), file = file.path(fit_dir, paste0(prefix, "__estimate_C.csv")))
    fwrite(data.table(model$y_pred), file = file.path(fit_dir, paste0(prefix, "__y_pred.csv")))
    fwrite(data.table(rmsfe), file = file.path(fit_dir, paste0(prefix, "__rmsfe.csv")))
}

save_fitting_results(model_fsplash, "fsplash", fit_dir)
save_fitting_results(model_ssfsplash, "ssfsplash", fit_dir)
save_fitting_results(model_splash_a0, "splash_a0", fit_dir)
save_fitting_results(model_splash_a05, "splash_a05", fit_dir)
save_fitting_results(model_pvar, "pvar", fit_dir, save_AB = FALSE)

# Save true values as well ,y_hat_true are the next time step predictions based on the true C matrix
fwrite(data.table(y_hat_true), file = file.path(fit_dir, "y_hat_true.csv"))
fwrite(data.table(A_true), file = file.path(fit_dir, "A_true.csv"))
fwrite(data.table(B_true), file = file.path(fit_dir, "B_true.csv"))
fwrite(data.table(y), file = file.path(fit_dir, "y_true.csv"))
