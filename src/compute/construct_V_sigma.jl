# isage example: julia construct_V_sigma.jl "data\y.csv" "version1"
# Output: version1_Vhat_d.csv, version1_sigma_hat.csv, version1_Vhat_d.feather, version1_sigma_hat.feather

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

"""
Read data from a csv file and return a matrix.

Args:
    path (str): path to the csv file.

Returns:
    Matrix{Float64}: matrix of shape (N, T) where N is the number of variables and T is the number of time periods.
"""
function read_data(path::String)::Matrix{Float64}
    return Matrix(CSV.read(path, DataFrame))
end


"""
Calculate covariance matrix of lag j, this assumes that the time index is over the columns.
Thus, y is a N x T matrix, where N is the number of variables and T is the number of time periods.

Args:
    y (Matrix{Float64}): matrix of shape (N, T) where N is the number of variables and T is the number of time periods.
    j (int): number of lags.
    bias (bool): whether to use bias correction in covariance calculation.

Returns:
    Matrix{Float64}: covariance matrix of lag j with shape (N, N).
"""
function calc_Σj(y::Matrix{Float64}, j::Int, bias::Bool=false)::Matrix{Float64}
    # Calculate covariance matrix of lag j
    return cov(y[:, 1:end-j], y[:, j+1:end], dims=2, corrected=bias)
end

"""
Return a banded matrix with bandwidth h on both sides of the main diagonal.

Parameters
----------
A : Matrix{Float64}
    The input matrix.
h : int
    The half-bandwidth to use.

Returns
-------
Matrix{Float64}
    A banded matrix with bandwidth h on both sides of the main diagonal.
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

Parameters
----------
Σ0 : Matrix{Float64}
    The covariance matrix of y.
Σ1 : Matrix{Float64}
    The covariance matrix of lag 1.
h : int, optional
    The half-bandwidth to use. Default is size(Σ1, 1)/4.
prebanded : bool, optional
    Whether Σ0 and Σ1 are already banded. Default is false.

Returns
-------
Matrix{Float64}
    The V matrix with dimensions N x (N+1), containing the columns of C' = [A B]', that are nonzero.
"""
function constr_Vhat(Σ0::Matrix{Float64}, Σ1::Matrix{Float64}; h::Int=0, prebanded::Bool=false)::Matrix{Float64}
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

Parameters
----------
Σ1 : Matrix{Float64}
    The banded autocovariance matrix of lag 1.
prebanded : bool, optional
    Whether Σ1 is already banded. Default is false.
h : int, optional
    The half-bandwidth to use. Default is size(Σ1, 1)/4.

Returns
-------
Vector{Float64}
    The vectorized autocovariance of lag 1.
"""
function vec_sigma_h(Σ1::Matrix{Float64}; prebanded::Bool=false, h::Int=0)::Vector{Float64}
    if !prebanded
        if h == 0
            h = div(size(Σ1, 1), 4)
        end
        Σ1 = band_matrix(Σ1, h)
    end
    return vec(Σ1')
end


"""
Calculate the active column indices for the matrix V_h^(d).

Parameters
----------
p : int
    The dimension of the matrix V.
h : int, optional
    The bandwidth to use. Default is p/4.

Returns
-------
Vector{Vector{Bool}}
    A p-length vector containing Boolean vectors representing the active column indices for each row.
"""
function active_cols(p::Int, h::Int=0)::Vector{Vector{Bool}}
    if h == 0
        h = div(p, 4)
    end

    active_set = [zeros(Bool, p * 2) for _ in 1:p]

    for i in 1:p
        lower_a = collect(max(1, i - h):max(0, i - 1))
        upper_a = collect((i+1):min(i + h, p))
        full_b = collect((p+max(1, (i - h))):(p+min(i + h, p)))
        selection = vcat(lower_a, upper_a, full_b)

        setindex!(active_set[i], ones(Bool, length(selection)), selection)
    end
    return active_set
end

"""
Construct a diagonal block matrix from the active columns of V.

Parameters
----------
V : Matrix{Float64}
    The V matrix with dimensions N x (N+1), containing the columns of C' = [A B]', that are nonzero.

Returns
-------
SparseMatrixCSC{Float64}
    The resulting diagonal block matrix with dimensions p^2 x K, where K is the number of nonzero columns in V.
"""
function constr_Vhat_d(V::Matrix{Float64})::SparseMatrixCSC{Float64}
    p = size(V, 1)
    active = active_cols(p)
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

Parameters
----------
y : Matrix{Float64}
    The data matrix with dimensions N x T.
j : int
    The lag of the covariance matrix to estimate.

Returns
-------
Matrix{Float64}
    The estimated covariance matrix of lag j with dimensions N x N.
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

Parameters
----------
y : Matrix{Float64}
    The data matrix with dimensions N x T.
j : int
    The lag of the covariance matrix to use in the estimation.
q : int, optional
    The number of bootstrap samples to use. Default is 500.

Returns
-------
int
    The estimated bandwidth.
"""
function bootstrap_estimator_Rj(y::Matrix{Float64}, j::Int, q::Int=500)::Int
    N, T = size(y)
    Σj = calc_Σj(y, j)
    Rj = zeros(Float64, div(N, 4))
    for i in 1:q
        bootstrap_Σj = bootstrap_estimator_Σj(y, j)
        for h in 1:div(N, 4)
            Rj[h] += norm((band_matrix(bootstrap_Σj, h) - Σj), 1)
        end
    end
    return argmin(Rj)
end


function main(prefix)
    # Read data 
    y = read_data(joinpath("out", "$(prefix)_y.csv"))

    # Do calculations
    Σ1 = calc_Σj(y, 1)
    Σ0 = calc_Σj(y, 0)

    Vhat = constr_Vhat(Σ0, Σ1)
    sigma_hat = vec_sigma_h(Σ1)
    Vhat_d = constr_Vhat_d(Vhat)

    # Write output
    CSV.write(joinpath("data", "simulation", "$(prefix)_Vhat_d.csv"), Tables.table(Vhat_d))
    CSV.write(joinpath("data", "simulation", "$(prefix)_sigma_hat.csv"), Tables.table(sigma_hat))
    return nothing
end

main(ARGS[1])



