using BenchmarkTools
using LinearAlgebra
using Distributions
using Random
using CSV
using DataFrames
using StatsBase
using ConfParser
using Statistics

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
function constr_V(Σ0::Matrix{Float64}, Σ1::Matrix{Float64})::Matrix{Float64}
    return [Σ1' Σ0]
end

"""
Construct the vectorized autocovariance of lag 1 from the banded autocovariance matrix.
"""
function vec_sigma_h(Σ1::Matrix{Float64}, h::Int; prebanded::Bool=false)::Vector{Float64}
    if !prebanded
        Σ1 = band_matrix(Σ1, h)
    end
    return vec(Σ1')
end

"""
This function calculates the active column indices for the matrix V_h^(d).
As a function of the dimesion p and the bandwidth h.
"""
function active_cols(p::Int, h::Int)::Vector{Int}
    
    return [i for i in 1:p if (i - h) > 0 && (i + h) <= p]
end

"""


# Playground
y = read_data("out/sim_y.csv")
A = read_data("out/sim_A.csv")

calc_Σ0(y)

s1 = calc_Σ1(y)
s0 = calc_Σ0(y)

N, T = size(A)

band_matrix(s1, 2)

vec(A)

constr_V(s0, s1)

vec_sigma_h(s1, 2, prebanded=false)