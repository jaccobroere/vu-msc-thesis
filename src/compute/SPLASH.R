# Install ncecessary libraries
# install.packages("splash_1.0.tar.gz", repos = NULL, type = "source")
# install.packages(c("data.table", "tictoc"))

# Import necessary libraries
library(splash)
library(data.table)
library(tictoc)

# Set working directory
# setwd("/Users/jacco/Documents/repos/vu-msc-thesis/")

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
path_prefix <- args[1]
# alpha <- as.numeric(args[2])

y <- t(fread("data/simulation/designA_T500_p100_y.csv", header = T, skip = 0))

# Fit SPLASH
tic()
model <- splash(y, alpha = 0.5, n_lambdas = 1, banded_covs=c(TRUE, TRUE), B = 500)
toc()

# Save R environment
save.image(file = paste0("out/", path_prefix, "_splash_env.RData"))
