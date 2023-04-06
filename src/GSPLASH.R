# Install necessary libraries
# install.packages(c("genlasso", "data.table", "arrow", "tictoc", "igraph"))

# Import necessary libraries
library(genlasso)
library(data.table)
library(arrow)
library(tictoc)
library(igraph)

# Set working directory
# setwd("/Users/jacco/Documents/repos/vu-msc-thesis/")

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
path_prefix <- args[1]
gamma <- as.numeric(args[2])
path1 <- paste0("out/", path_prefix, "_sigma_hat.csv")
path2 <- paste0("out/", path_prefix, "_Vhat_d.csv")
path3 <- paste0("out/", path_prefix, "_graph.csv")
# Read data
sigma_hat <- t(fread(path1, header = T, skip = 0))
Vhat_d <- as.matrix(fread(path2, header = T, skip = 0)
gr <- read_graph(path3, format = "graphml")

# Fit fused lasso model for 'genlasso'
tic()
model <- genlasso::fusedlasso(y = sigma_hat, X = Vhat_d, graph = gr, gamma = 0, verbose = TRUE)
toc()

# TODO: Save results