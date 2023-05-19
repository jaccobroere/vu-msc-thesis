# Load all the necessary packages
PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/R/vhat_sigmahat.R")
source("src/compute/R/utils.R")
setwd(system("echo $ZHU_DIR", intern = TRUE))
source("R/opt.R")
source("R/gen_data.R")
library(data.table)
library(igraph)
library(genlasso)
setwd(PROJ_DIR)

# Read CLI arguments
args <- commandArgs(trailingOnly = TRUE)
sim_design_id <- ifelse(length(args) < 1, "designB_T500_p36", args[1])
uuidtag <- ifelse(length(args) < 2, "A650C52A-1CB2-40C9-A358-D55204F8E54A", args[2])

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "detlam", uuidtag)
lambdas_dir <- file.path(PROJ_DIR, "out/simulation/lambdas", sim_design_id)
sim_id_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id, "detlam")

# Parse paths
path_reg_graph <- file.path(sim_id_dir, "reg_graph.graphml")
path_sym_graph <- file.path(sim_id_dir, "sym_graph.graphml")
path_y <- file.path(data_dir, "y.csv")
path_bandwidth <- file.path(sim_id_dir, "bootstrap_bandwidths.csv")

# Load the data
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")
y <- as.matrix(fread(path_y, header = T, skip = 0))
bandwidths <- as.data.frame(fread(path_bandwidth, header = T, skip = 0))

# Define functions to perform the search
fit_gfsplash.cv.iter <- function(y, bandwidths, graph, alpha, nlambdas = 20, lambda.min.ratio = 1e-4, ...) {
    # Read problem dimensionality
    p <- dim(y)[1]
    m <- ecount(graph)
    k <- vcount(graph)
    h <- floor(p / 4) # Bandwidth for the A and B matrix

    # Parse bandwidths
    h0 <- bandwidths[1]
    h1 <- bandwidths[2]

    # Reparameterize the lambda and alpha for the linreg_path_v2 function, see paper Zhu 2017 Augmented ADMM
    gamma <- alpha / (1 - alpha)

    # Split y into training and testing sets
    train_idx <- (floor(dim(y)[2] / 5) * 4)
    y_train <- y[, 1:train_idx]
    # y_test <- y[, ((floor(dim(y)[2] / 5) * 4) + 1):dim(y)[2]]
    y_test <- y[, (train_idx + 1):(train_idx + 20)]

    # Construct the D matrix from the penalty graph
    D <- as(getDgSparse(graph = graph), "TsparseMatrix")

    # Fit the model on each fold and save the results
    C <- array(NA, dim = c(p, p, nlambdas))
    A <- array(NA, dim = c(p, p, nlambdas))
    B <- array(NA, dim = c(p, p, nlambdas))
    y_pred <- array(NA, dim = c(p, dim(y_test)[2], nlambdas))
    errors <- rep(NA, nlambdas)

    # Construct Vhat_d and sigma_hat from the training validation
    t0 <- Sys.time()
    Sigma0 <- calc_Sigma_j(y_train, 0)
    Sigma1 <- calc_Sigma_j(y_train, 1)
    Vhat <- construct_Vhat(Sigma0, Sigma1, h0, h1)
    sigma_hat <- construct_sigma_hat(Sigma1, h1)
    Vhat_d <- as.matrix(construct_Vhat_d(Vhat))

    # Calculate lambda_0 and construct lambda grid accordingly
    lambda_0 <- calc_lambda_0_gfsplash(sigma_hat, Vhat_d, graph, alpha) # This gives lambda_0 in the (lambda, alpha) parametrization
    lambda_grid_linreg <- rev(lambda_0 * 10^seq(log10(lambda.min.ratio), 0, length.out = nlambdas)) / (1 + gamma)

    # Fit the model on the training set
    model_grid <- linreg_path_v2(
        Y = as.vector(sigma_hat),
        X = Vhat_d,
        val = D@x,
        idx = D@i + 1,
        jdx = D@j + 1,
        lambda_graph = as.vector(lambda_grid_linreg),
        gamma = gamma,
        p = dim(Vhat_d)[2],
        m = dim(D)[1],
        ...
    )
    t1 <- Sys.time()

    # Compute the prediction error on the validation set
    for (j in 1:nlambdas) {
        coef <- as.vector(model_grid$beta_path[, j])
        AB <- coef_to_AB(coef, p)
        A[, , j] <- AB$A
        B[, , j] <- AB$B
        C[, , j] <- AB_to_C(AB$A, AB$B)
        y_pred[, , j] <- predict_with_C(C[, , j], y_train, y_test)
        errors[j] <- calc_msfe(y_pred[, , j], y_test)
    }

    # Get the best lambda value
    best_idx <- which.min(errors)
    best_msfe <- errors[best_idx]
    best_lambda <- lambda_grid_linreg[best_idx] * (1 + gamma)

    output_list <- list(
        model = model_grid,
        A = A,
        B = B,
        C = C,
        errors = errors,
        y_pred = y_pred,
        best_msfe = best_msfe,
        best_lambda = best_lambda,
        best_idx = best_idx,
        runtime = difftime(t1, t0, units = "secs")[[1]]
    )
    return(output_list)
}

# Function to determine the optimal lambda for a given specification of GF-SPLASH
fit_gfsplash.lambda_grid <- function(y, bandwidths, graph, alpha, path, nlambdas = 20, lambda.min.ratio = 1e-4) {
    # Call the function for each specification and create file if it does not exist yet
    df_lam <- data.frame(matrix(0, ncol = nlambdas, nrow = 0))
    colnames(df_lam) <- paste0("ilam", 1:nlambdas)
    if (!file.exists(path)) {
        write.table(df_lam, path, row.names = FALSE, col.names = TRUE, append = FALSE, sep = ",")
    }
    # Get the errors per lambda index from the model
    model_grid <- fit_gfsplash.cv.iter(y, bandwidths, graph, alpha = alpha, nlambdas = nlambdas, lambda.min.ratio = lambda.min.ratio, standard_ADMM = TRUE, max_num_iter = 1e3)
    out <- matrix(data = model_grid$errors, nrow = 1, ncol = nlambdas)
    # Append the results to the table
    write.table(out, path, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    return(df_lam)
}

# Fit the grids and save the results
grid_gfsplash_a05 <- fit_gfsplash.lambda_grid(y, bandwidths, reg_gr, alpha = 0.5, path = file.path(lambdas_dir, "grid_gfsplash_a05.csv"))
grid_gfsplash_sym_a0 <- fit_gfsplash.lambda_grid(y, bandwidths, sym_gr, alpha = 0, path = file.path(lambdas_dir, "grid_gfsplash_sym_a0.csv"))
grid_gfsplash_sym_a05 <- fit_gfsplash.lambda_grid(y, bandwidths, sym_gr, alpha = 0.5, path = file.path(lambdas_dir, "grid_gfsplash_sym_a05.csv"))
