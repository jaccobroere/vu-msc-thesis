PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
library(Matrix)
library(igraph)
library(genlasso)

# Set up directories
sim_design_id <- "designB_T1000_p9"
data_dir <- file.path(PROJ_DIR, "data/simulation", sim_design_id)
save_dir <- file.path(PROJ_DIR, "out/figures/penalty_matrices")

# Parse paths
path_reg_graph <- file.path(data_dir, "reg_graph.graphml")
path_sym_graph <- file.path(data_dir, "sym_graph.graphml")

# Load the data
reg_gr <- read_graph(path_reg_graph, format = "graphml")
sym_gr <- read_graph(path_sym_graph, format = "graphml")

# Compute the Dmatrices for GF-SPLASH
alpha <- 0.3
D <- as(getDg(graph = reg_gr), "TsparseMatrix")
D_alpha <- rbind(D * -1, diag(ncol(D)))
writeMM(D_alpha, file.path(save_dir, "D_alpha_GFSPLASH.mtx"))

D_sym <- as(getDg(graph = sym_gr), "TsparseMatrix")
D_alpha_sym <- rbind(D_sym * -1, diag(ncol(D_sym)))
writeMM(D_alpha_sym, file.path(save_dir, "D_alpha_sym_GFSPLASH.mtx"))
