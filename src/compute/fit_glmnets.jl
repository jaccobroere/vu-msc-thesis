using Pkg
Pkg.activate(joinpath(ENV["PROJ_DIR"], "juliaenv"))

using GLMNet
using SparseArrays, MatrixMarket
using CSV, DataFrames
using RCall
using CrossValidation

R"""
source('src/compute/utils.R')
"""


PROJ_DIR = ENV["PROJ_DIR"]
include(joinpath(PROJ_DIR, "src", "compute", "transform_bootstrap_graph.jl"))

# Read CLI arguments
args = ARGS
sim_design_id = length(args) < 1 ? "designB_T500_p81" : args[1]
uuidtag = length(args) < 2 ? "29D1E659-10D3-42F9-BF51-B83B26709C63" : args[2]

# Set up directories
data_dir = joinpath(PROJ_DIR, "data/simulation", sim_design_id, "mc", uuidtag)
fit_dir = joinpath(PROJ_DIR, "out/simulation/fit", sim_design_id, uuidtag)
lambdas_dir = joinpath(PROJ_DIR, "out/simulation/lambdas", sim_design_id)
sim_id_dir = joinpath(PROJ_DIR, "data/simulation", sim_design_id, "mc")

# Read in the CSV or matrix files
sigma_hat = read_data(joinpath(data_dir, "sigma_hat.csv"))
Vhat_d = mmread(joinpath(data_dir, "Vhat_d.mtx"))
h0, h1 = read_data(joinpath(data_dir, "bootstrap_bandwidths.csv")) |> Matrix{Int}
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
start_time = time()
XD1 = Vhat_d * Dtilde_inv # Same as Vhat_d * inv(Dtilde)
X1 = XD1[:, 1:m]
X2 = XD1[:, (m+1):size(XD1)[2]]
X2_plus = pinv(Matrix(X2))
# Transform the input to LASSO objective
P = X2 * X2_plus
ytilde = vec((I - P) * sigma_hat)
Xtilde = (I - P) * X1

@time model = glmnet(Xtilde, ytilde, alpha=1, intercept=false, standardize=false, nlambda=20)
theta1 = model.betas[:, 50]
theta2 = vec(X2_plus * (sigma_hat - X1 * theta1))
coef = Dtilde_inv * [theta1; theta2]
elapsed_time = (time() - start_time)

# Convert to matrices A and B
A = rcopy(R"coef_to_AB($coef, $p)$A")
B = rcopy(R"coef_to_AB($coef, $p)$B")
C = rcopy(R"AB_to_C($A, $B)")

function cv_fsplash(y, Dtilde_inv, graph, h0, h1, n_fold, n_lambda)
    train_size = div(size(y, 2), 5) * 4
    y_train, y_test = y[:, 1:train_size], y[:, (train_size+1):size(y, 2)]
    cv = SlidingWindow(y_train, floor(Int, train_size * 0.8), floor(Int, train_size * 0.2 / n_fold))

    results = zeros(n_fold, n_lambda)
    for (fold_idx, (y_trainval, y_val)) in enumerate(cv)
        # Calcualte the new Vhat_d
        Vhat_d, sigma_hat = calc_Vhat_d_sigma_hat_nb(y_trainval, h0, h1)
        # Prepare the model input
        m, k = ne(graph), nv(graph)
        p = Integer(sqrt(size(Vhat_d)[1]))

        # Use linear system solvers for faster computation of the change of variables (see Tibshirani and Taylor, 2011)
        XD1 = Vhat_d * Dtilde_inv # Same as Vhat_d * inv(Dtilde)
        X1 = XD1[:, 1:m]
        X2 = XD1[:, (m+1):size(XD1)[2]]
        X2_plus = pinv(Matrix(X2))
        # Transform the input to LASSO objective
        IminusP = I - X2 * X2_plus
        ytilde = vec(IminusP * sigma_hat)
        Xtilde = IminusP * X1

        model = glmnet(Xtilde, ytilde, alpha=1, intercept=false, standardize=false, nlambda=n_lambda)

        for (lam_idx, col) ∈ enumerate(eachcol(model.betas))
            theta1 = col
            theta2 = vec(X2_plus * (sigma_hat - X1 * theta1))
            coef = Dtilde_inv * [theta1; theta2]
            A = rcopy(R"coef_to_AB($coef, $p)$A")
            B = rcopy(R"coef_to_AB($coef, $p)$B")
            C = rcopy(R"AB_to_C($A, $B)")
            # Prepare input to R functions
            y_in = hcat(y_trainval, y_val)
            train_idx = size(y_trainval, 2)
            # Calculate validation error and save in results matrix
            y_hat = rcopy(R"predict_with_C($C, $y_in, train_idx=$train_idx)")
            msfe = rcopy(R"calc_msfe($y_val, $y_hat)")
            results[fold_idx, lam_idx] = msfe
        end
    end
    return results
end


# Cross-validation
@time cv_fsplash(y, Dtilde_inv, graph, h0, h1, 5, 20)

div(size(y, 2), 5)

div(train_size, 5)


SlidingWindow(y_train,)