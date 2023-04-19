source("src/compute/utils.R")
source("src/compute/model_wrappers.R")
library(data.table)
library(FGSG)
library(tictoc)

# Set up directories
data_dir <- "/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/"
path_prefix <- "designB_T1000_p100"

# Parse paths
path1 <- paste0(data_dir, path_prefix, "_sigma_hat.csv")
path2 <- paste0(data_dir, path_prefix, "_Vhat_d.csv")
path3 <- paste0(data_dir, path_prefix, "_graph.graphml")
path4 <- paste0(data_dir, path_prefix, "_A.csv")
path5 <- paste0(data_dir, path_prefix, "_B.csv")
path6 <- paste0(data_dir, path_prefix, "_y.csv")

# Load the data
sigma_hat <- t(fread(path1, header = T, skip = 0))
Vhat_d <- as.matrix(fread(path2, header = T, skip = 0))
gr <- read_graph(path3, format = "graphml")

# Load the true values
A <- as.matrix(fread(path4, header = T, skip = 0))
B <- as.matrix(fread(path5, header = T, skip = 0))
y <- as.matrix(fread(path6, header = T, skip = 0))

# Convert the graph to an edge vector
edge_vector <- as.vector(t(as_edgelist(gr)))

# Print the dimensions of the data
message(cat("The dimension of y: ", dim(y)[1], dim(y)[2]))

# Fit a single solution of the GFLASSO
lambda <- 0.086
tic()
smodel <- gflasso(y = sigma_hat, A = Vhat_d, tp = edge_vector, s1 = 0, s2 = lambda)
toc()
length(smodel$weight)

# Fit the entire solution path of the GFLASSO
# tic()
# fmodel <- genlasso::fusedlasso(y = sigma_hat, X = Vhat_d, graph = gr, gamma = 1, verbose = FALSE)
# toc()

# Fit the a single solution using (Augmented) ADMM
m1 <- admm_gsplash(sigma_hat, Vhat_d, gr, lambda, 1, standard_ADMM = TRUE)
m2 <- admm_gsplash(sigma_hat, Vhat_d, gr, lambda, 1, standard_ADMM = FALSE)

m3 <- regular_splash(y, banded_covs = c(TRUE, TRUE), B = 500, alphas = c(0.5), lambdas = c(lambda))
