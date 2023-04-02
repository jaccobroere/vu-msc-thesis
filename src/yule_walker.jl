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

function read_data(path::String)::Matrix{Float64}
    return Matrix(CSV.read(path, DataFrame))
end

function band_matrix!(A::Matrix{Float64}, h::Int)::Matrix{Float64}
    """Return a banded matrix with bandwidth h on both sides of the main diagonal."""
    N, T = size(A)
    for i in 1:N
        for j in 1:T
            if abs(i - j) > h
                A[i, j] = 0.0
            end
        end
    end
    return A
end

y = read_data("out/sim_y.csv")
A = read_data("out/sim_A.csv")

calc_Σ0(y)

s1 = calc_Σ1(y)

N, T = size(A)

band_matrix!(s1, 2)

s1