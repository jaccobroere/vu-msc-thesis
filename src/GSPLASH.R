# Import necessary libraries
library(genlasso)
library(data.table)
library(arrow)
library(tictoc)

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
