# isage example: julia construct_V_sigma.jl "data\y.csv" "version1"
# Output: version1_Vhat_d.csv, version1_sigma_hat.csv, version1_Vhat_d.feather, version1_sigma_hat.feather
using Pkg
Pkg.activate(joinpath(ENV["JULIA_DIR"]), io=devnull)
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

# function bootstrap_estimator_Σj_norm(y::Matrix{Float64}, j::Int)::Matrix{Float64}
#     N, T = size(y)
#     Σj = zeros(N, N)
#     y_mean = mean(y, dims=2)
#     for t in 1:(T-j)
#         u_t = diagm(rand(Normal(0, 1) + 1, N))
#         @inbounds Σj += u_t * (y[:, t] - y_mean) * (y[:, t+j] - y_mean)'
#     end
#     return Σj / T
# end

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
        bootstrap_Σ0 = bootstrap_estimator_Σj_exp(y, 0)
        bootstrap_Σ1 = bootstrap_estimator_Σj_exp(y, 1)
        for h in 1:(N-1)
            @inbounds R0[i, h] += opnorm((band_matrix(bootstrap_Σ0, h) - Σ0), 1) / q
            @inbounds R1[i, h] += opnorm((band_matrix(bootstrap_Σ1, h) - Σ1), 1) / q
        end
    end
    return (argmin(vec(sum(R0, dims=1))), argmin(vec(sum(R1, dims=1))))
end

function main(sim_design_id, uuidtag)
    if uuidtag !== nothing
        path = joinpath("data", "simulation", sim_design_id, uuidtag)
        if !isdir(path)
            mkpath(path)
        end
    else
        path = joinpath("data", "simulation", sim_design_id)
    end

    # Create and invert Dtilde if it does not exist for this dimension (only need to be calculated once)
    path_sim = joinpath("data", "simulation", sim_design_id)
    if !isfile(joinpath(path_sim, "Dtilde.mtx"))
        # Read data 
        y_read = read_data(joinpath(path, "y.csv"))
        p = size(y_read, 1)
        h = div((p - 1), 4)
        # Subset the first 80% of the data
        y = y_read[:, 1:div(size(y_read, 2), 5)*4]
        # Bootstrap the bandwidth
        h0, h1 = bootstrap_estimator_R(y, 2000)
        # Construct underlying graph 
        regular_graph = create_gsplash_graph(size(y, 1), symmetric=false) # Bandwitdh is set to floor(p/4) by default
        symmetric_graph = create_gsplash_graph(size(y, 1), symmetric=true)
        # Calculate the Dtilde matrix and its inverse, for F-SPLASH and SSF-SPLASH, respectively
        Dtilde = calc_Dtilde_sparse(regular_graph)
        Dtilde_inv = inv_Dtilde_sparse(regular_graph)
        Dtilde_SSF = calc_Dtilde_SSF_sparse(regular_graph, h)
        Dtilde_SSF_inv = inv_Dtilde_SSF_sparse(regular_graph, h)
        
        # Save the matrices in sparse matrix format
        mmwrite(joinpath(path_sim, "Dtilde.mtx"), Dtilde)
        mmwrite(joinpath(path_sim, "Dtilde_inv.mtx"), Dtilde_inv)
        mmwrite(joinpath(path_sim, "Dtilde_SSF.mtx"), Dtilde_SSF)
        mmwrite(joinpath(path_sim, "Dtilde_SSF_inv.mtx"), Dtilde_SSF_inv)
        mmwrite(joinpath(path_sim, "Dtilde_SDF.mtx"), Dtilde_SDF)
        mmwrite(joinpath(path_sim, "Dtilde_SDF_inv.mtx"), Dtilde_SDF_inv)

        # Save the graphs
        save_graph_as_gml(regular_graph, joinpath(path_sim, "reg_graph.graphml"))
        save_graph_as_gml(symmetric_graph, joinpath(path_sim, "sym_graph.graphml"))
        # Save the bandwidths
        save_bandwiths_bootstrap(joinpath(path_sim, "bootstrap_bandwidths.csv"), h0, h1)
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line argument
    sim_design_id = ARGS[1]
    uuidtag = length(ARGS) >= 2 ? ARGS[2] : nothing

    # Run main function
    main(sim_design_id, uuidtag)
end