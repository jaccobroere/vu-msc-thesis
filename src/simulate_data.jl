using BenchmarkTools
using LinearAlgebra
using Distributions
using Random
using CSV
using DataFrames

# Set random seed
Random.seed!(2023)


function D_fused(p::Int)::Matrix{Int8}
    D = -Matrix{Int8}(I, p, p)
    for i in 1:p-1
        D[i, i+1] = 1
    end
    return D[1:p-1, :]
end

# Generates a banded matrix with random values between -2 and 2, with bandwidth h on both sides of the main diagonal
function generate_banded_matrix(p::Int, h::Int)::Matrix{Float64}
    matrix = zeros(p, p)
    for i in 1:p
        for j in 1:p
            if abs(i - j) <= h
                matrix[i, j] = rand(Uniform(-2, 2))
            end
        end
    end
    spectral_norm = maximum(eigvals(matrix' * matrix))
    return 0.5 * matrix / sqrt(spectral_norm)
end

function generate_random_vector(p::Int)::Array{Float64,1}
    return rand(Normal(0, 1), p)
end

function generate_errors_over_time(T::Int, p::Int)::Array{Array{Float64,1},1}
    errors = [generate_random_vector(p) for _ in 1:T]
    return errors
end

function simulate_svar(A::Matrix{Float64}, B::Matrix{Float64}, errors::Array{Array{Float64,1},1})::Matrix{Float64}
    p = size(A, 1)
    T = length(errors)
    y = zeros(p, T)
    C = inv(I - A) * B
    print(dim(C))
    for t in 1:T
        if t == 1
            y[:, t] = errors[t]
            continue
        end

        y[:, t] = C * y[:, t-1] + errors[t]
    end

    return y
end

function run_simulation(p::Int, T::Int, h_A::Int, h_B::Int, path_prefix::String="sim", write::Bool=true)::Matrix{Float64}
    # Define the constants
    p = 100
    T = 500
    h_A = 3
    h_B = 3

    A = generate_banded_matrix(p, h_A)
    B = generate_banded_matrix(p, h_B)
    errors = generate_errors_over_time(T, p)
    y = simulate_svar(A, B, errors)

    if write == true
        CSV.write(joinpath("out", path_prefix * "_" * "A.csv"), DataFrame(A, :auto))
        CSV.write(joinpath("out", path_prefix * "_" * "B.csv"), DataFrame(B, :auto))
        CSV.write(joinpath("out", path_prefix * "_" * "y.csv"), DataFrame(y, :auto))
    end

    return y
end

run_simulation(100, 500, 3, 3, "sim", true)

e = generate_errors_over_time(500, 100)

sigma_e = cov(e)

normm = norm(sigma_e)

function calc_sigma_e(e::Vector{Vector{Float64}})::Matrix{Float64}
    sigma_e = zeros(size(e[1], 1), size(e[1], 1))
    emean = mean(e)
    for i in 1:length(e)
        sigma_e += (e[i] .- emean) * (e[i] .- emean)'
    end
    return sigma_e / (length(e) - 1)
end

@benchmark calc_sigma_e(e)

@benchmark cov(e)

calc_sigma_e(e) == cov(e)


# Check similarity between cov and calc_sigma_e
calc_sigma_e(e) .== cov(e)

calc_sigma_e(e)

cov(e)

mean(e)