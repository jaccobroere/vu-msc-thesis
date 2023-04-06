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

# Read data
sigma_hat <- t(fread("out/v1_sigma_hat.csv", header = T, skip = 0))
Vhat_d <- as.matrix(fread("out/v1_Vhat_d.csv", header = T, skip = 0))
D <- as.matrix(fread("out/v1_D.csv", header = T, skip = 0))

# Fit fused lasso model for 'genlasso'
tic()
model <- genlasso::fusedlasso1d(y = sigma_hat, X = Vhat_d, verbose = TRUE)
toc()


# Defining igraph object that will be used to create the graph connecting elements on the same non-zero diagonals of the matrix
gr <- read_graph("out/gsplash_graph_p7_h2.graphml", format = "graphml")

plot(gr)
