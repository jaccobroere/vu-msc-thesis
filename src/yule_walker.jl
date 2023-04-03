using Base: setindex
using BlockDiagonals
using BenchmarkTools
using LinearAlgebra
using Distributions
using Random
using CSV
using DataFrames
using StatsBase
using ConfParser
using Statistics
using SparseArrays

# Set random seed
Random.seed!(2023)

"""
Read data from a csv file and return a matrix.
"""
function read_data(path::String)::Matrix{Float64}
    return Matrix(CSV.read(path, DataFrame))
end


"""
Calculate covariance matrix of lag 0, this assumes that the time index is over the columns.
Thus, y is a N x T matrix, where N is the number of variables and T is the number of time periods.
"""
function calc_Σ0(y::Matrix{Float64}, bias::Bool=false)::Matrix{Float64}
    """Calculate covariance matrix of lag 0, this assumes that the time index is over the columns.
    Thus, y is a N x T matrix, where N is the number of variables and T is the number of time periods.
    """
    # Calculate contemporaneous covariance matrix
    return cov(y, dims=2, corrected=bias)
end

"""
Calculate covariance matrix of lag 1, this assumes that the time index is over the columns.
Thus, y is a N x T matrix, where N is the number of variables and T is the number of time periods.
"""
function calc_Σ1(y::Matrix{Float64}, bias::Bool=false)::Matrix{Float64}
    # Calculate covariance matrix of lag 1
    return cov(y[:, 1:end-1], y[:, 2:end], dims=2, corrected=bias)
end


"""
Return a banded matrix with bandwidth h on both sides of the main diagonal.
"""
function band_matrix(A::Matrix{Float64}, h::Int)::Matrix{Float64}
    N, T = size(A)
    B = zeros(N, T)
    for i in 1:N
        for j in max((i - h), 1):min((i + h), T)
            B[i, j] = A[i, j]
        end
    end
    return B
end


"""
Construct the V matrix containing the columns of C' = [A B]', that are nonzero
"""
function constr_Vhat(Σ0::Matrix{Float64}, Σ1::Matrix{Float64}; prebanded::Bool=false)::Matrix{Float64}
    if !prebanded
        if h == 0
            h = floor(Int, size(Σ1, 1) / 4)
        end
        Σ1 = band_matrix(Σ1, h)
        Σ0 = band_matrix(Σ0, h)
    end
    return [Σ1' Σ0]
end

"""
Construct the vectorized autocovariance of lag 1 from the banded autocovariance matrix.
"""
function vec_sigma_h(Σ1::Matrix{Float64}; prebanded::Bool=false, h::Int=0)::Vector{Float64}
    if !prebanded
        if h == 0
            h = floor(Int, size(Σ1, 1) / 4)
        end
        Σ1 = band_matrix(Σ1, h)
    end
    return vec(Σ1')
end

"""
This function calculates the active column indices for the matrix V_h^(d).
As a function of the dimesion p and the bandwidth h.
"""
function active_cols(p::Int, h::Int=0)::Vector{Vector{Bool}}
    if h == 0
        h = floor(Int, p / 4)
    end

    active_set = [zeros(Bool, p * 2) for _ in 1:p]

    for i in 1:p
        lower_a = collect(max(1, i - h):max(0, i - 1))
        upper_a = collect(min(i + 1, p):min(i + h, p - 1))
        full_b = collect((p+max(1, (i - h))):(p+min(i + h, p - 1)))
        selection = vcat(lower_a, upper_a, full_b)
        setindex!(active_set[i], ones(Bool, length(selection)), selection)
    end
    return active_set
end

"""
This function takes the active columns of V and stacks them into a diagonal block matrix
"""
function constr_Vhat_d(V::Matrix{Float64})::SparseMatrixCSC{Float64}
    p = size(V, 1)
    active = active_cols(p)
    res = [V[:, active[i]] for i in 1:p]
    Vhat_d = spzeros(p^2, sum(sum(active)))

    row_index = col_index = 1
    for (M, act) in zip(res, active)
        setindex!(Vhat_d, M, collect(row_index:(row_index+p-1)), collect(col_index:(col_index+sum(act)-1)))
        row_index += p
        col_index += sum(act)
    end
    return Vhat_d
end

function main(path::String)::nothing
    # Read data 
    y = read_data(path)

    Σ1 = calc_Σ1(y)
    Σ0 = calc_Σ0(y)

    Vhat = constr_V(s0, s1)
    sigma_hat = vec_sigma_h(s1)
    Vhat_d = constr_Vhat_d(Vhat)

    CSV.write("$(ARGS[2])_Vhat_d.csv", Vhat_d)
    CSV.write("$(ARGS[2])_sigma_hat.csv", sigma_vec)
end

main(ARGS[1])