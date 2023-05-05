using Graphs: SimpleGraphs, connected_components!
using LinearAlgebra
using Base: vect, create_expr_cache
using SparseArrays
using GLMNet
using MatrixMarket
using CSV, Tables

include(joinpath("construct_graph.jl"))
include(joinpath("transform_bootstrap_graph.jl"))
include(joinpath("utils.jl"))

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

###### PLAYGROUND #######
y = read_data(joinpath("/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/designB_T500_p9", "y.csv"))
Vhat_d = mmread(joinpath("/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/designB_T500_p9", "Vhat_d.mtx"))
sigma_hat = read_data(joinpath("/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/designB_T500_p9", "sigma_hat.csv"))


# Vhat_d, sigma_hat = calc_Vhat_d_sigma_hat(y)

G = create_gsplash_graph(9)
D = Matrix(incidence_matrix(G, oriented=true))'
Dtilde = calc_Dtilde(G)
null_vecs = null_space_graph(G)

D

for i in 1:9
    for j in 1:60
        println(dot(vec(D[j, :]), vec(null_vecs[:, i])))
    end
end

# Duan approach
n, p, m = size(Vhat_d)[1], size(Vhat_d)[2], size(D)[1]

Φ = vcat(Vhat_d, D)
QR = qr(Φ)
Q = Matrix(QR.Q)
R = Matrix(QR.R)
Q1, Q2 = Q[1:n, :], Q[n+1:end, :]

pQ2 = Q2' * inv(Q2 * Q2')

pinv(Q2)

pQ2
H = Q1 * pQ2
z = vec((I - Q1 * (I - pQ2 * Q2) * Q1') * sigma_hat)


using GLMNet, MLBase, Lasso
path = glmnet(H, z, intercept=false, lambda=[0.005050447], alpha=1, standardize=false)
theta = path.betas
beta = R^(-1) * (I - pQ2 * Q2) * Q1' * sigma_hat + R^(-1) * pQ2 * beta


## Tibshirani approach
n, p, m = size(Vhat_d)[1], size(Vhat_d)[2], size(D)[1]
graph = create_gsplash_graph(9)
Dtilde = calc_Dtilde(G)
XD1 = Vhat_d * inv(Dtilde)
X1, X2 = XD1[:, 1:m], XD1[:, (m+1):end]
P = X2 * inv(X2' * X2) * X2'

ytilde = vec((I - P) * sigma_hat)
Xtilde = (I - P) * X1
path_tib = glmnet(Xtilde, ytilde, intercept=false, lambda=[0.005050447], alpha=1, standardize=false)

theta1 = path_tib.betas
theta2 = inv(X2' * X2) * X2' * (sigma_hat - X1 * theta1)
beta = inv(Dtilde) * vcat(theta1, theta2)

## Write data
CSV.write(joinpath("/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/designB_T500_p9", "H.csv"), Tables.table(H))
CSV.write(joinpath("/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/designB_T500_p9", "z.csv"), Tables.table(z))

connected_components(G)