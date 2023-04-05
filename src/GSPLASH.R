# Import necessary libraries
library(genlasso)
library(data.table)
library(arrow)
library(tictoc)
library(igraph)

# Set working directory
setwd("/Users/jacco/Documents/repos/vu-msc-thesis/")

# Read data
sigma_hat <- t(fread("out/v1_sigma_hat.csv", header = T, skip = 0))
Vhat_d <- as.matrix(fread("out/v1_Vhat_d.csv", header = T, skip = 0))
D <- as.matrix(fread("out/v1_D.csv", header = T, skip = 0))

# Fit fused lasso model for 'genlasso'
tic()
model <- genlasso::fusedlasso1d(y = sigma_hat, X = Vhat_d, verbose = TRUE)
toc()


# Defining igraph object that will be used to create the graph connecting elements on the same non-zero diagonals of the matrix
edges <- c(1, 2, 1, 3, 1, 5, 2, 4, 2, 5, 3, 6, 3, 7, 3, 8, 6, 7, 6, 8)
gr <- graph(edges = edges, directed = FALSE)

plot(gr)
