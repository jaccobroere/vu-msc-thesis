PROJ_DIR <- system("echo $PROJ_DIR", intern = TRUE)
setwd(PROJ_DIR)
source("src/compute/utils.R")
library(data.table)

# Read arguments from CLI
args <- commandArgs(trailingOnly = TRUE)
path_prefix <- "designA_T500_p25" # path_prefix <- args[1]
sim_id <- "1343" # sim_id_dir <- args[2]

# Set up directories
data_dir <- file.path(PROJ_DIR, "data/simulation/")
out_dir <- file.path(PROJ_DIR, "out/")
lambdas_dir <- file.path(PROJ_DIR, "out/simulation/lambdas/")
sim_id_dir <- file.path(lambdas_dir, sim_id)

# Read in the lambda dataframes
grid_lam_reg_a0 <- fread(file.path(sim_id_dir, "reg_a0.csv"), header = T, sep = ",", skip = 0)
grid_lam_reg_a05 <- fread(file.path(sim_id_dir, "reg_a05.csv"), sep = ",", skip = 0)
grid_lam_sym_a0 <- fread(file.path(sim_id_dir, "sym_a0.csv"), sep = ",", skip = 0)
grid_lam_sym_a05 <- fread(file.path(sim_id_dir, "sym_a05.csv"), sep = ",", skip = 0)

# Calculate the average for each column
avg_lam_reg_a0 <- colMeans(grid_lam_reg_a0)
avg_lam_reg_a05 <- colMeans(grid_lam_reg_a05)
avg_lam_sym_a0 <- colMeans(grid_lam_sym_a0)
avg_lam_sym_a05 <- colMeans(grid_lam_sym_a05)

# Calculate what lambda was the best
best_lam_reg_a0 <- avg_lam_reg_a0[which.min(avg_lam_reg_a0)]
best_lam_reg_a05 <- avg_lam_reg_a05[which.min(avg_lam_reg_a05)]
best_lam_sym_a0 <- avg_lam_sym_a0[which.min(avg_lam_sym_a0)]
best_lam_sym_a05 <- avg_lam_sym_a05[which.min(avg_lam_sym_a05)]

# Save the grid of best lambdas selected lambda values
best_lambdas <- data.frame(
    variable = c("best_lam_reg_a0", "best_lam_reg_a05", "best_lam_sym_a0", "best_lam_sym_a05"),
    lambda = c(best_lam_reg_a0, best_lam_reg_a05, best_lam_sym_a0, best_lam_sym_a05)
)
write.table(best_lambdas, paste0(out_dir, sim_id_dir, path_prefix, "_best_lambdas.csv"), row.names = FALSE, append = TRUE, col.names = FALSE)
