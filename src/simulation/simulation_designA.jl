#=
This simulation design replicates the simulation design A from Reuvers and Wijler (2021)
=#
using Pkg
Pkg.activate(joinpath(ENV["JULIA_DIR"]), io=devnull)
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

function draw_until_valid(p::Int, h_A::Int, h_B::Int)::Tuple{Matrix{Float64},Matrix{Float64}}
    # First try to draw A and B
    A = design_A_generate_A(p, h_A)
    B = design_A_generate_B(p, h_B)

    #  If the spectral radius of A or B is larger than 1, draw again
    stability = norm(inv(I - A) * B, 2) <= 0.95
    while !(stability)
        A = design_A_generate_A(size(A, 1), h_A)
        B = design_A_generate_B(size(B, 1), h_B)
        stability = norm(inv(I - A) * B, 2) <= 0.95
    end
    return A, B
end

function run_simulation(p::Int, T::Int, sim_design_id::String="sim", write::Bool=true, uuidtag::Union{String,Nothing}=nothing, train_size::Float64=0.8, h_A::Int=3, h_B::Int=3)::Matrix{Float64}
    T = Int(T / train_size)
    A, B = draw_until_valid(p, h_A, h_B)
    errors = generate_errors_over_time(T, p)
    y = simulate_svar(A, B, errors)

    if write == true
        write_simulation_output(A, B, y, sim_design_id, uuidtag)
    end

    return y
end


if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line argument
    sim_design_id = ARGS[1]
    T, p = parse_sim_design_id(sim_design_id)
    uuidtag = length(ARGS) >= 2 ? ARGS[2] : nothing

    # Run simulation
    run_simulation(p, T, sim_design_id, true, uuidtag)
end


