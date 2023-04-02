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

function calc_Σ0(y::Matrix{Float64}, bias::Bool=false)::Matrix{Float64}
    # Calculate contemporaneous covariance matrix
    Σ = zeros(size(y, 1), size(y, 1))
    mean_time = mean(y, dims=2)
    for i in eachcol(y)
        Σ += (i .- mean_time) * (i .- mean_time)'
    end

    # Return biased estimate if bias is true
    if bias
        return Σ / size(y, 2)
    end

    # Otherwise return unbiased estimate
    return Σ / (size(y, 2) - 1)
end


function calc_Σ1(y::Int)::Matrix{Float64}
    """
    Calculate covariance matrix of lag 1, this assumes that the time index is over the columns.
    Thus, y is a N x T matrix, where N is the number of variables and T is the number of time periods."""
    # Calculate covariance matrix of lag 1
    return cov(y[:, 1:end-1], y[:, 2:end], dims=2)
end

function read_data(path::String)::Matrix{Float64}
    return Matrix(CSV.read(path, DataFrame))
end


y = read_data("out/sim_y.csv")


cov(y)


mean(y)

mean(y, dims=2)

calc_Σ0(y)

size(y)[1]

cov(y, dims=2)

calc_Σ0(y) - cov(y, dims=2)

cov(y[:, 1:end-1], y[:, 2:end], dims=2)

# Select all but the last column
y[:, 1:end-1]