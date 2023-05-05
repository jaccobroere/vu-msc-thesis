using Graphs: SimpleGraphs, connected_components!
using LinearAlgebra
using Base: vect, create_expr_cache
using SparseArrays
using GLMNet
using MatrixMarket
using CSV, Tables
using IterativeSolvers

include(joinpath("construct_graph.jl"))
include(joinpath("transform_bootstrap_graph.jl"))
include(joinpath("utils.jl"))

function null_space_graph_sparse(graph::SimpleGraph{Int64})::SparseMatrixCSC{Int64}
    null_vecs = spzeros(Int64, nv(graph), nv(graph) - ne(graph))
    conn = connected_components(graph)
    for (j, vec) in enumerate(conn)
        for i in vec
            null_vecs[i, j] = 1
        end
    end
    return null_vecs
end

function null_space_graph(graph::SimpleGraph{Int64})::Matrix{Int64}
    null_vecs = zeros(Int64, nv(graph), nv(graph) - ne(graph))
    conn = connected_components(graph)
    for (j, vec) in enumerate(conn)
        for i in vec
            null_vecs[i, j] = 1
        end
    end
    return null_vecs
end

function calc_Dtilde(graph::SimpleGraph)::Matrix{Int64}
    D_prime = Matrix(incidence_matrix(graph, oriented=true)) # Incidence matrix needs to be transposed before obtaining D^(G)
    null = null_space_graph(graph)
    return vcat(D_prime', null') # Thus transpose here
end

function calc_Dtilde_sparse(graph::SimpleGraph)::SparseMatrixCSC
    D_prime = incidence_matrix(graph, oriented=true) # Incidence matrix needs to be transposed before obtaining D^(G)
    null = null_space_graph(graph)
    return vcat(D_prime', null') # Thus transpose here
end

function calc_inv_LU(A)
    F = lu(A)
    A_inv = inv(F.U) * inv(F.L) * F.P
    return A_inv
end

function calc_inv_LU_sparse(A::SparseMatrixCSC)::SparseMatrixCSC
    F = lu(A)
    A_inv = inv(F.U) * inv(F.L) * F.P
    return A_inv
end

###### PLAYGROUND #######
y = read_data(joinpath("/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/designB_T500_p9/mc/6A446516-84EF-49C9-9435-23E00DB22756", "y.csv"))
Vhat_d = mmread(joinpath("/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/designB_T500_p9/mc/6A446516-84EF-49C9-9435-23E00DB22756", "Vhat_d.mtx"))
sigma_hat = read_data(joinpath("/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/designB_T500_p9/mc/6A446516-84EF-49C9-9435-23E00DB22756", "sigma_hat.csv"))


# Vhat_d, sigma_hat = calc_Vhat_d_sigma_hat(y)

G = create_gsplash_graph(9)
# D = Matrix(incidence_matrix(G, oriented=true))'
D = Matrix(incidence_matrix(G, oriented=true))'
Dtilde = calc_Dtilde(G)

null_vecs_sparse = null_space_graph_sparse(G)
D_sparse = incidence_matrix(G, oriented=true)'
Dtilde_sparse = calc_Dtilde_sparse(G)

@time inv(Dtilde)
@time calc_inv_LU(Dtilde)
# @btime calc_inv_LU_sparse(Dtilde_sparse)

F = lu(Dtilde_sparse)
inv(Matrix(F.U)) * inv(Matrix(F.L)) * F.P
Dtilde_sparse


# @btime (Dtilde' \ Vhat_d')'
# @btime Vhat_d * inv(Dtilde)

function model_fast_fusion(sigma_hat::Matrix{Float64}, Vhat_d::SparseMatrixCSC{Float64}, graph::SimpleGraph{Int64})
    # Calculate D_tilde by extending it with orthogonal rows to a square matrix
    Dtilde = calc_Dtilde(graph)
    m, p = ne(graph), nv(graph)

    # Use linear system solvers for faster computation of the change of variables (see Tibshirani and Taylor, 2011)
    XD1 = (Dtilde' \ Vhat_d')' # Same as Vhat_d * inv(Dtilde)
    X1, X2 = XD1[:, 1:m], XD1[:, (m+1):end]
    X2_plus = (X2' * X2) \ X2' # Same as inv(X2' * X2) * X2'

    # Transform the input to LASSO objective
    P = X2 * X2_plus
    ytilde = vec((I - P) * sigma_hat)
    Xtilde = (I - P) * X1

    # Solve LASSO
    path = glmnet(Xtilde, ytilde, intercept=false, lambda=[0.615848211066027], alpha=1, standardize=false)

    # Transform back to original variables
    theta1 = vec(path.betas)
    theta2 = X2_plus * (sigma_hat - X1 * theta1)
    theta = vcat(theta1, theta2)
    coef = Dtilde \ theta # Same as inv(Dtilde) * theta

    return coef
end

@time model_fast_fusion(sigma_hat, Vhat_d, G)

# Duan approachn, p, m = size(Vhat_d)[1], size(Vhat_d)[2], size(D)[1]
using IterativeSolvers

@time gmres(Dtilde', Vhat_d', verbose=false, log=false)
@time Dtilde' \ Vhat_d'