"""
This simulation design replicates the simulation design A from Reuvers and Wijler (2021)
"""

using LinearAlgebra
using Distributions
using Random
using CSV
using DataFrames

# Include simulation utils
include(joinpath(dirname(abspath(@__FILE__)), "simulation_utils.jl"))

# Generates a banded matrix with random values between -2 and 2, with bandwidth h on both sides of the main diagonal
function design_A_generate_A(p::Int, h::Int)::Matrix{Float64}
    matrix = zeros(p, p)
    for i in 1:p
        for j in 1:p
            if abs(i - j) == h
                matrix[i, j] = rand(Uniform(-2, 2))
            end
            if abs(i - j) < h && i != j
                u = rand(Bernoulli(0.4))
                matrix[i, j] = u * 0 + (1 - u) * rand(Normal(0, 1))
            end
        end
    end
    eta = rand(Uniform(0.4, 0.8))
    spectral_norm = maximum(eigvals(matrix' * matrix))
    return eta * matrix / sqrt(spectral_norm)
end

function design_A_generate_B(p::Int, h::Int)::Matrix{Float64}
    matrix = zeros(p, p)
    for i in 1:p
        for j in 1:p
            if abs(i - j) == h
                matrix[i, j] = rand(Uniform(-2, 2))
            end
            if abs(i - j) < h
                u = rand(Bernoulli(0.4))
                matrix[i, j] = u * 0 + (1 - u) * rand(Normal(0, 1))
            end
        end
    end
    eta = rand(Uniform(0.4, 0.8))
    spectral_norm = maximum(eigvals(matrix' * matrix))
    return eta * matrix / sqrt(spectral_norm)
end

function run_simulation(p::Int, T::Int, h_A::Int, h_B::Int, path_prefix::String="sim", write::Bool=true, train_size::Float64=0.8)::Matrix{Float64}
    T = Int(T / train_size)
    A = design_A_generate_A(p, h_A)
    B = design_A_generate_B(p, h_B)
    errors = generate_errors_over_time(T, p)
    y = simulate_svar(A, B, errors)

    if write == true
        write_simulation_output(A, B, y, path_prefix)
    end

    return y
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line argument        
    p = parse(Int, ARGS[1])
    T = parse(Int, ARGS[2])
    h_A = parse(Int, ARGS[3])
    h_B = parse(Int, ARGS[4])
    path_prefix = ARGS[5]

    # Run simulation
    run_simulation(p, T, h_A, h_B, path_prefix, true)
end


