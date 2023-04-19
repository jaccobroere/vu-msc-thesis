# isage example: julia construct_V_sigma.jl "data\y.csv" "version1"
# Output: version1_Vhat_d.csv, version1_sigma_hat.csv, version1_Vhat_d.feather, version1_sigma_hat.feather
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
using SparseArrays
using DataFrames

# Set random seed
Random.seed!(2023)

path = dirname(abspath(@__FILE__))
include(joinpath(path, "construct_graph.jl"))

"""
Read data from a csv file and return a matrix.

## Arguments
- `path::String`: path to the csv file.

## Returns
- `Matrix{Float64}`: matrix of shape (N, T) where N is the number of variables and T is the number of time periods.
"""
function read_data(path::String)::Matrix{Float64}
    return Matrix(CSV.read(path, DataFrame))
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
function calc_Σj(y::Matrix{Float64}, j::Int, bias::Bool=false)::Matrix{Float64}
    # Calculate covariance matrix of lag j
    return cov(y[:, 1:end-j], y[:, j+1:end], dims=2, corrected=bias)
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
- `h::int=0`: The half-bandwidth to use. Default is size(Σ1, 1)/4.
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
function bootstrap_estimator_Σj(y::Matrix{Float64}, j::Int)::Matrix{Float64}
    N, T = size(y)
    Σj = zeros(N, N)
    for t in 1:(T-j)
        u_t = rand(Exponential(1))
        @inbounds Σj += u_t * y[:, t] * y[:, t+j]'
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
    R0 = zeros(Float64, q, div(N, 4))
    R1 = zeros(Float64, q, div(N, 4))
    Threads.@threads for i in 1:q
        bootstrap_Σ0 = bootstrap_estimator_Σj(y, 0)
        bootstrap_Σ1 = bootstrap_estimator_Σj(y, 1)
        for h in 1:div(N, 4)
            @inbounds R0[i, h] += norm((band_matrix(bootstrap_Σ0, h) - Σ0), 1)
            @inbounds R1[i, h] += norm((band_matrix(bootstrap_Σ1, h) - Σ1), 1)
        end
    end
    return (argmin(vec(sum(R0, dims=1))), argmin(vec(sum(R1, dims=1))))
end

function main(prefix)
    # Read data 
    y = read_data(joinpath("data", "simulation", "$(prefix)_y.csv"))

    # Bootstrap the bandwidth
    h0, h1 = bootstrap_estimator_R(y, 500)

    # Do calculations
    Σ1 = calc_Σj(y, 1)
    Σ0 = calc_Σj(y, 0)

    Vhat = constr_Vhat(Σ0, Σ1, h0, h1)
    sigma_hat = vec_sigma_h(Σ1, h1)
    Vhat_d = constr_Vhat_d(Vhat) # Bandwitdh is set to floor(p/4) by default

    # Construct underlying graph 
    graph = create_gsplash_graph(size(y, 1)) # Bandwitdh is set to floor(p/4) by default

    # Write output
    CSV.write(joinpath("data", "simulation", "$(prefix)_Vhat_d.csv"), Tables.table(Vhat_d))
    CSV.write(joinpath("data", "simulation", "$(prefix)_sigma_hat.csv"), Tables.table(sigma_hat))
    save_graph_as_gml(graph, joinpath("data", "simulation", "$(prefix)_graph.graphml"))

    return nothing
end

main(ARGS[1])

## TESTING
y = read_data(joinpath("data", "simulation", "designA_T500_p100_y.csv"))

bootstrap_estimator_R(y, 500)