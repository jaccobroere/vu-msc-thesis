using Pkg
Pkg.activate(joinpath(ENV["PROJ_DIR"], "juliaenv"))

using GLMNet
using SparseArrays, MatrixMarket
using CSV, DataFrames
using RCall
Pkg.add("RCall")


PROJ_DIR = ENV["PROJ_DIR"]
include(joinpath(PROJ_DIR, "src", "compute", "transform_bootstrap_graph.jl"))

# Read CLI arguments
args = ARGS
sim_design_id = length(args) < 1 ? "designB_T500_p49" : args[1]
uuidtag = length(args) < 2 ? "" : args[2]

# Set up directories
data_dir = joinpath(PROJ_DIR, "data/simulation", sim_design_id, "mc", uuidtag)
fit_dir = joinpath(PROJ_DIR, "out/simulation/fit", sim_design_id, uuidtag)
lambdas_dir = joinpath(PROJ_DIR, "out/simulation/lambdas", sim_design_id)
sim_id_dir = joinpath(PROJ_DIR, "data/simulation", sim_design_id)

# Read in the CSV or matrix files
sigma_hat = read_data(joinpath(data_dir, "sigma_hat.csv"))
Vhat_d = mmread(joinpath(data_dir, "Vhat_d.mtx"))
y = read_data(joinpath(data_dir, "y.csv"))
A = read_data(joinpath(data_dir, "A.csv"))
B = read_data(joinpath(data_dir, "B.csv"))

Dtilde_inv = mmread(joinpath(sim_id_dir, "Dtilde_inv.mtx"))
Dtilde = mmread(joinpath(sim_id_dir, "Dtilde.mtx"))
Dtilde_SSF = mmread(joinpath(sim_id_dir, "Dtilde_SSF.mtx"))
Dtilde_SSF_inv = mmread(joinpath(sim_id_dir, "Dtilde_SSF_inv.mtx"))

graph = create_gsplash_graph(size(y, 1))

# Prepare the model input
m = ne(graph)
k = nv(graph)
p = Integer(sqrt(size(Vhat_d)[1]))

# Use linear system solvers for faster computation of the change of variables (see Tibshirani and Taylor, 2011)
XD1 = Vhat_d * Dtilde_inv # Same as Vhat_d * inv(Dtilde)
X1 = XD1[:, 1:m]
X2 = XD1[:, (m+1):size(XD1)[2]]
X2_plus = pinv(Matrix(X2))
# Transform the input to LASSO objective
P = X2 * X2_plus
ytilde = vec((I - P) * sigma_hat)
Xtilde = (I - P) * X1

@time model = glmnet(Xtilde, ytilde, alpha=1, intercept=false, standardize=false)
theta1 = model.betas[:, 50]
theta2 = vec(X2_plus * (sigma_hat - X1 * theta1))
coef = Dtilde_inv * [theta1; theta2]


using RCall

R"""
source('src/compute/utils.R')
"""
using BenchmarkTools
@benchmark rcopy(R"coef_to_AB($coef, $p)$A")


AB = coef_to_AB(coef, p)
A = AB.A
B = AB.B
C = AB_to_C(A, B)
