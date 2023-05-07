# isage example: julia construct_V_sigma.jl "data\y.csv" "version1"
# Output: version1_Vhat_d.csv, version1_sigma_hat.csv, version1_Vhat_d.feather, version1_sigma_hat.feather
using Pkg
Pkg.activate(joinpath(ENV["PROJ_DIR"], "juliaenv"))
using Distributed
Threads.nthreads()

using LoopVectorization
using Tables
using Base: setindex
using LinearAlgebra
using Random
using CSV
using StatsBase
using Statistics
using Distributions
using SparseArrays, MatrixMarket
using DataFrames

# Set random seed
# Random.seed!(2023)

path = dirname(abspath(@__FILE__))
include(joinpath(path, "construct_graph.jl"))
include(joinpath(path, "utils.jl"))

"""
Read data from a csv file and return a matrix.

## Arguments
- `path::String`: path to the csv file.

## Returns
- `Matrix{Float64}`: matrix of shape (N, T) where N is the number of variables and T is the number of time periods.
"""
function read_data(path::String)::Matrix{Float64}
    return CSV.read(path, DataFrame; types=Float64) |> Matrix
end


"""
Calculate covariance matrix of lag j, this assumes that the time index is over the columns.
Thus, y is a N x T matrix, where N is the number of variables and T is the number of time periods.

## Arguments
- `y::Matrix{Float64}`: matrix of shape (N, T) where N is the number of variables and T is the number of time periods.
- `j::int`: number of lags.
- `bias::bool=false`: whether to use bias correction in covariance calculation.

## Returns
- `Matrix{Float64}`: covariance matrix of lag j with shape (N, N).
"""
function calc_Σj(y::Matrix{Float64}, j::Int)::Matrix{Float64}
    p, T = size(y)
    sigma_j = zeros(p, p)
    ymean = mean(y, dims=2)
    for t in 1:(T-j)
        @inbounds sigma_j += (y[:, t] - ymean) * (y[:, t+j] - ymean)'
    end
    return sigma_j / T
end


"""
Return a banded matrix with bandwidth h on both sides of the main diagonal.

## Arguments
- `A::Matrix{Float64}`: The input matrix.
- `h::int`: The half-bandwidth to use.

## Returns
- `Matrix{Float64}`: A banded matrix with bandwidth h on both sides of the main diagonal.
"""
function band_matrix(A::Matrix{Float64}, h::Int)::Matrix{Float64}
    N, T = size(A)
    B = zeros(N, T)
    for i in 1:N
        for j in max((i - h), 1):min((i + h), T)
            @inbounds B[i, j] = A[i, j]
        end
    end
    return B
end

"""
Construct the V matrix containing the columns of C' = [A B]', that are nonzero.

## Arguments
- `Σ0::Matrix{Float64}`: The covariance matrix of y.
- `Σ1::Matrix{Float64}`: The covariance matrix of lag 1.
- `h::int=0`: The bandwidth to use, if 0, then no bandwidth is used.
- `prebanded::bool=false`: Whether Σ0 and Σ1 are already banded. Default is false.

## Returns
- `Matrix{Float64}`: The V matrix with dimensions N x (N+1), containing the columns of C' = [A B]', that are nonzero.
"""
function constr_Vhat(Σ0::Matrix{Float64}, Σ1::Matrix{Float64}, h0::Int=0, h1::Int=0)::Matrix{Float64}
    if h0 == 0 || h1 == 0
        return [Σ1' Σ0]
    end

    Σ0 = band_matrix(Σ0, h0)
    Σ1 = band_matrix(Σ1, h1)
    return [Σ1' Σ0]
end

"""
Construct the vectorized autocovariance of lag 1 from the banded autocovariance matrix.

## Arguments
- `Σ1::Matrix{Float64}`: The banded autocovariance matrix of lag 1.
- `h::Int=0`: The half-bandwidth to use. Default is size(Σ1, 1)/4.

## Returns
- `Vector{Float64}`: The vectorized autocovariance of lag 1.
"""
function vec_sigma_h(Σ1::Matrix{Float64}, h1::Int=0)::Vector{Float64}
    if h1 == 0
        return vec(Σ1')
    end
    Σ1 = band_matrix(Σ1, h1)
    return vec(Σ1')
end


"""
Calculate the active column indices for the matrix V_h^(d).

## Arguments
- `p::Int`: The dimension of the matrix V.
- `bandwidth::Int=0`: The bandwidth to use. Default is p/4.

## Returns
- `Vector{Vector{Bool}}`: A p-length vector containing Boolean vectors representing the active column indices for each row.
"""
function active_cols(p::Int, bandwidth::Int=0)::Vector{Vector{Bool}}
    if bandwidth == 0
        bandwidth = div(p, 4)
    end

    active_set = [zeros(Bool, p * 2) for _ in 1:p]

    for i in 1:p
        lower_a = collect(max(1, i - bandwidth):max(0, i - 1))
        upper_a = collect((i+1):min(i + bandwidth, p))
        full_b = collect((p+max(1, (i - bandwidth))):(p+min(i + bandwidth, p)))
        selection = vcat(lower_a, upper_a, full_b)

        setindex!(active_set[i], ones(Bool, length(selection)), selection)
    end
    return active_set
end

"""
Construct a diagonal block matrix from the active columns of V.

## Arguments
- `V::Matrix{Float64}`: The V matrix with dimensions N x (N+1), containing the columns of C' = [A B]', that are nonzero.
- `bandwidth::Int=0`: The half-bandwidth to use. Default is size(V, 1)/4.

## Returns
- `SparseMatrixCSC{Float64}`: The resulting diagonal block matrix with dimensions p^2 x K, where K is the number of nonzero columns in V.
"""
function constr_Vhat_d(V::Matrix{Float64}, bandwidth::Int=0)::SparseMatrixCSC{Float64}
    p = size(V, 1)
    active = active_cols(p, bandwidth)
    res = [V[:, active[i]] for i in 1:p]
    Vhat_d = spzeros(p^2, sum(sum(active)))

    row_index = col_index = 1
    for (M, act) in zip(res, active)
        @inbounds setindex!(Vhat_d, M, collect(row_index:(row_index+p-1)), collect(col_index:(col_index+sum(act)-1)))
        row_index += p
        col_index += sum(act)
    end
    return Vhat_d
end

"""
Estimate the covariance matrix of lag j using a bootstrap method (Guo et al. 2016).

## Arguments
- `y::Matrix{Float64}`: The data matrix with dimensions N x T.
- `j::Int`: The lag of the covariance matrix to estimate.

## Returns
- `Matrix{Float64}`: The estimated covariance matrix of lag j with dimensions N x N.
"""
function bootstrap_estimator_Σj_exp(y::Matrix{Float64}, j::Int)::Matrix{Float64}
    N, T = size(y)
    Σj = zeros(N, N)
    y_mean = mean(y, dims=2)
    for t in 1:(T-j)
        u_t = rand(Exponential(1))
        @inbounds Σj += u_t * (y[:, t] - y_mean) * (y[:, t+j] - y_mean)'
    end
    return Σj / T
end

function bootstrap_estimator_Σj_norm(y::Matrix{Float64}, j::Int)::Matrix{Float64}
    N, T = size(y)
    Σj = zeros(N, N)
    y_mean = mean(y, dims=2)
    for t in 1:(T-j)
        u_t = diagm(rand(Normal(0, 1) + 1, N))
        @inbounds Σj += u_t * (y[:, t] - y_mean) * (y[:, t+j] - y_mean)'
    end
    return Σj / T
end


"""
Estimate the bandwidth for banded autocovariance estimation using a bootstrap method (Guo et al. 2016).

## Arguments
- `y::Matrix{Float64}`: The data matrix with dimensions N x T.
- `q::Int=500`: The number of bootstrap samples to use. Default is 500.

## Returns
- `Tuple{Int, Int}`: The estimated bandwidths (h_Σ0, h_Σ1).
"""
function bootstrap_estimator_R(y::Matrix{Float64}, q::Int=500)::Tuple{Int,Int}
    N, T = size(y)
    Σ0 = calc_Σj(y, 0)
    Σ1 = calc_Σj(y, 1)
    R0 = zeros(Float64, q, N - 1)
    R1 = zeros(Float64, q, N - 1)
    Threads.@threads for i in 1:q
        bootstrap_Σ0 = bootstrap_estimator_Σj_norm(y, 0)
        bootstrap_Σ1 = bootstrap_estimator_Σj_norm(y, 1)
        for h in 1:(N-1)
            @inbounds R0[i, h] += norm((band_matrix(bootstrap_Σ0, h) - Σ0), 1) / q
            @inbounds R1[i, h] += norm((band_matrix(bootstrap_Σ1, h) - Σ1), 1) / q
        end
    end
    return (argmin(vec(sum(R0, dims=1))), argmin(vec(sum(R1, dims=1))))
end

# """
# Constuct Vhat_d without bootstrapping, meant for use in the cross-validation steps
# """
# function calc_Vhat_d_nb(y::Matrix{Float64})::SparseMatrixCSC{Float64}
#     N, T = size(y)
#     Σ0 = calc_Σj(y, 0)
#     Σ1 = calc_Σj(y, 1)
#     Vhat = constr_Vhat(Σ0, Σ1)
#     return constr_Vhat_d(Vhat)
# end

function main(sim_design_id, uuidtag)
    if uuidtag !== nothing
        path = joinpath("data", "simulation", sim_design_id, uuidtag)
        if !isdir(path)
            mkpath(path)
        end
    else
        path = joinpath("data", "simulation", sim_design_id)
    end

    # Read data 
    y = read_data(joinpath(path, "y.csv"))
    p = size(y, 1)
    h = div(p, 4)

    # Subset the first 80% of the data
    y = y[:, 1:div(size(y, 2), 5)*4]

    # Bootstrap the bandwidth
    h0, h1 = bootstrap_estimator_R(y, 500)

    # Do calculations
    Σ1 = calc_Σj(y, 1)
    Σ0 = calc_Σj(y, 0)

    Vhat = constr_Vhat(Σ0, Σ1, h0, h1)
    sigma_hat = vec_sigma_h(Σ1, h1)
    Vhat_d = constr_Vhat_d(Vhat) # Bandwitdh is set to floor(p/4) by default

    # Construct underlying graph 
    regular_graph = create_gsplash_graph(size(y, 1), symmetric=false) # Bandwitdh is set to floor(p/4) by default
    symmetric_graph = create_gsplash_graph(size(y, 1), symmetric=true)

    # Create and invert Dtilde if it does not exist for this dimension (only need to be calculated once)
    path_sim = joinpath("data", "simulation", sim_design_id)
    if !isfile(joinpath(path_sim, "Dtilde.mtx"))
        # Calculate the Dtilde matrix and its inverse, for F-SPLASH and SSF-SPLASH, respectively
        Dtilde = calc_Dtilde_sparse(regular_graph)
        Dtilde_inv = inv_Dtilde_sparse(regular_graph)
        Dtilde_SSF = calc_Dtilde_SSF_sparse(regular_graph, h)
        Dtilde_SSF_inv = inv_Dtilde_SSF_sparse(regular_graph, h)
        # Save the matrices in sparse matrix format
        mmwrite(joinpath("data", "simulation", sim_design_id, "Dtilde.mtx"), Dtilde)
        mmwrite(joinpath("data", "simulation", sim_design_id, "Dtilde_inv.mtx"), Dtilde_inv)
        mmwrite(joinpath("data", "simulation", sim_design_id, "Dtilde_SSF.mtx"), Dtilde_SSF)
        mmwrite(joinpath("data", "simulation", sim_design_id, "Dtilde_SSF_inv.mtx"), Dtilde_SSF_inv)
    end

    # Write output
    mmwrite(joinpath(path, "Vhat_d.mtx"), Vhat_d)
    CSV.write(joinpath(path, "sigma_hat.csv"), Tables.table(sigma_hat))
    save_graph_as_gml(regular_graph, joinpath(path, "reg_graph.graphml"))
    save_graph_as_gml(symmetric_graph, joinpath(path, "sym_graph.graphml"))

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    sim_design_id = ARGS[1]
    uuidtag = length(ARGS) >= 2 ? ARGS[2] : nothing
    main(sim_design_id, uuidtag)
end

# ## TESTING
# y = read_data(joinpath("/Users/jacco/Documents/repos/vu-msc-thesis/data/simulation/designB_T500_p49/mc/ED04541C-0149-45C6-8F59-6D5A21BFB0C2", "y.csv"))


# y_train = y[:, 1:div(size(y, 2), 5)*4]
# bootstrap_estimator_R(y_train, 999)

# Σ = calc_Σj(y_train, 1)

# band_matrix(Σ, 23)



# # Random 5x5 matrix
# A = rand(5, 5)
# band_matrix(A, 2)
