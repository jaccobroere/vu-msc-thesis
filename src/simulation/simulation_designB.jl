"""
This simulation design replicates the simulation design B from Reuvers and Wijler (2021).
In this simulation design:
- A has first horizontal and first vertical interactions between neighbors
- B is a diagonal matrix
"""

using LinearAlgebra
using Distributions
using Random
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

function run_simulation(p::Int, T::Int, file_prefix::String="sim", write::Bool=true, train_size::Float64=0.8, h_A::Int=3, h_B::Int=3)::Matrix{Float64}
    T = Int(T / train_size)
    A = design_A_generate_A(p, h_A)
    B = design_A_generate_B(p, h_B)
    errors = generate_errors_over_time(T, p)
    y = simulate_svar(A, B, errors)

    if write == true
        write_simulation_output(A, B, y, file_prefix)
    end

    return y
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line argument
    p = parse(Int, ARGS[1])
    m = Int(sqrt(p))
    T = parse(Int, ARGS[2])
    file_prefix = ARGS[5]
    uuidtag = ARGS[6]

    print(uuidtag)

    # Run simulation
    run_simulation(m, T, file_prefix, true)
end



