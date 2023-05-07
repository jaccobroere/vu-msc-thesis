using SparseArrays
using Graphs: SimpleGraphs, connected_components!
using LinearAlgebra
using Base: vect, create_expr_cache
using SparseArrays
using GLMNet
using MatrixMarket
using CSV, Tables
using IterativeSolvers

PROJ_DIR = ENV["PROJ_DIR"]
include(joinpath(PROJ_DIR, "src", "compute", "construct_graph.jl"))
include(joinpath(PROJ_DIR, "src", "compute", "utils.jl"))
include(joinpath(PROJ_DIR, "src", "compute", "transform_bootstrap_graph.jl"))

function null_space_graph(graph::SimpleGraph{Int64})::Matrix{Int64}
    null_vecs = zeros(Float64, nv(graph), nv(graph) - ne(graph))
    conn = connected_components(graph)
    for (j, vec) in enumerate(conn)
        for i in vec
            null_vecs[i, j] = 1.0
        end
    end
    return null_vecs
end

function calc_Dtilde(graph::SimpleGraph)::Matrix{Int64}
    D_prime = Matrix(incidence_matrix(graph, oriented=true)) # Incidence matrix needs to be transposed before obtaining D^(G)
    null = null_space_graph(graph)
    return vcat(D_prime', null') # Thus transpose here
end

function calc_Dtilde_sparse(graph::SimpleGraph)::SparseMatrixCSC{Float64}
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


p = 5
G = create_gsplash_graph(p)
D = Matrix(incidence_matrix(G, oriented=true))'
k, m = nv(G), ne(G)

