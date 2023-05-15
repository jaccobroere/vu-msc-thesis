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