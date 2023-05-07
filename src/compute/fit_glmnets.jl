using Pkg
Pkg.activate(joinpath(ENV["PROJ_DIR"], "juliaenv"))

using GLMNet
using SparseArrays, MatrixMarket
using CSV, DataFrames

PROJ_DIR = ENV["PROJ_DIR"]
include(joinpath(PROJ_DIR, "src", "compute", "transform_bootstrap_graph.jl"))

# Read CLI arguments
args = ARGS
sim_design_id = length(args) < 1 ? "designB_T500_p49" : args[1]
uuidtag = length(args) < 2 ? "ED04541C-0149-45C6-8F59-6D5A21BFB0C2" : args[2]

# Set up directories
data_dir = joinpath(PROJ_DIR, "data/simulation", sim_design_id, "mc", uuidtag)
fit_dir = joinpath(PROJ_DIR, "out/simulation/fit", sim_design_id, uuidtag)
lambdas_dir = joinpath(PROJ_DIR, "out/simulation/lambdas", sim_design_id)
sim_id_dir = joinpath(PROJ_DIR, "out/simulation", sim_design_id)

# Parse paths
path_sigma_hat = joinpath(data_dir, "sigma_hat.csv")
path_Vhat_d = joinpath(data_dir, "Vhat_d.mtx")
path_reg_graph = joinpath(data_dir, "reg_graph.graphml")
path_sym_graph = joinpath(data_dir, "sym_graph.graphml")
path_y = joinpath(data_dir, "y.csv")
path_A = joinpath(data_dir, "A.csv")
path_B = joinpath(data_dir, "B.csv")
path_Dtilde_inv = joinpath(sim_id_dir, "Dtilde_inv.mtx")
path_Dtilde = joinpath(sim_id_dir, "Dtilde.mtx")

# Read in the CSV or matrix files
sigma_hat = read_data(path_sigma_hat)
Vhat_d = mmread(path_Vhat_d)
y = read_data(path_y)
A = read_data(path_A)
B = read_data(path_B)
Dtilde_inv = mmread(path_Dtilde_inv)
Dtilde = mmread(path_Dtilde)

# 