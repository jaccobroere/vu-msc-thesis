"""
This simulation design replicates the simulation design B from Reuvers and Wijler (2021).
In this simulation design:
- A has first horizontal and first vertical interactions between neighbors
- B is a diagonal matrix
"""

using BenchmarkTools
using LinearAlgebra
using Distributions
using Random
using CSV
using DataFrames

# Include simulation utils
include(joinpath(dirname(abspath(@__FILE__)), "simulation_utils.jl"))

function design_B_generate_A(m::Int)::Matrix{Float64}
    p = m^2
    A = zeros(p, p)
    for i in 1:p
        for j in 1:p
            if abs(i - j) == 1 || abs(i - j) == m
                setindex!(A, 0.2, i, j)
            end
        end
    end
    return A
end

function design_B_generate_B(m::Int)::Matrix{Float64}
    p = m^2
    B = zeros(p, p)
    for i in 1:p
        setindex!(B, 0.2, i, i)
    end
    return B
end

function run_simulation(m::Int, T::Int, path_prefix::String="sim", write::Bool=true, train_size::Float64=0.8)::Nothing
    T = Int(T / train_size)
    A = design_B_generate_A(m)
    B = design_B_generate_B(m)
    errors = generate_errors_over_time(T, p)
    y = simulate_svar(A, B, errors)

    if write == true
        write_simulation_output(A, B, y, path_prefix)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line argument
    p = parse(Int, ARGS[1])
    m = Int(sqrt(p))
    T = parse(Int, ARGS[2])
    path_prefix = ARGS[5]

    # Run simulation
    run_simulation(m, T, path_prefix, true)
end



