# Import libraries and set working directory for Zhu's code
# setwd("/Users/jacco/Documents/repos/vu-msc-thesis/admm_src_zhu")
source("src/compute/model_wrappers.R")

library(splash)
library(genlasso)
library(data.table)
library(igraph)
library(tictoc)
library(FGSG)

# Set up directories
data_dir <- "/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/"
path_prefix <- "designB_T500_p25"

# Parse paths
path1 <- paste0(data_dir, path_prefix, "_sigma_hat.csv")
path2 <- paste0(data_dir, path_prefix, "_Vhat_d.csv")
path3 <- paste0(data_dir, path_prefix, "_graph.graphml")

# Load the data
sigma_hat <- t(fread(path1, header = T, skip = 0))
Vhat_d <- as.matrix(fread(path2, header = T, skip = 0))
gr <- read_graph(path3, format = "graphml")
edge_vector <- as.vector(t(as_edgelist(gr)))

# Print dimensions of the data
print(dim(sigma_hat))
print(dim(Vhat_d))
print(vcount(gr))
print(length(edge_vector))

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
