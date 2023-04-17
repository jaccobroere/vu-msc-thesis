setwd("/Users/jacco/Documents/repos/vu-msc-thesis/Zhu")
source("R/opt.R")
source("R/gen_data.R")
# Install the packages if necessary
# install_github("monogenea/gflasso")
# install.packages(c("genlasso","igraph","tictoc","FGSG","data.table"))
# library(splash)
# library(gflasso)
library(genlasso)
library(data.table)
library(igraph)
library(tictoc)
library(FGSG)

# Set up directories
data_dir <- "/Users/jacco/Documents/repos/vu-msc-thesis/out/"
path_prefix <- "exp_p10"

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
tic()
fmodel <- genlasso::fusedlasso(y = sigma_hat, X = Vhat_d, graph = gr, gamma = 1, verbose = TRUE)
toc()

# Fit the a single solution using (Augmented) ADMM
D <- as(getDgSparse(graph = gr), "TsparseMatrix")
idx <- D@i + 1
jdx <- D@j + 1
val <- D@x

tic()
admmmodel <- linreg_path_v2(Y = sigma_hat, X = Vhat_d, val = val, idx = idx, jdx = jdx, lambda_graph = lambda, gamma = 1, p = dim(Vhat_d)[2], m = dim(D)[1], standard_ADMM = TRUE)
toc()