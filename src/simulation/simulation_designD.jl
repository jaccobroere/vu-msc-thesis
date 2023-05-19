using Distributions: continuous_distributions
#=
This simulation design replicates the simulation design B from Reuvers and Wijler (2021).
In this simulation design:
- A has first horizontal and first vertical interactions between neighbors
- B is a diagonal matrix
=#
using Pkg
Pkg.activate(joinpath(ENV["JULIA_DIR"]), io=devnull)
using LinearAlgebra
using Distributions
using Random
using DataFrames

# Include simulation utils
include(joinpath(dirname(abspath(@__FILE__)), "simulation_utils.jl"))
include(joinpath(dirname(abspath(@__FILE__)), "simulation_designC.jl"))

function design_D_generate_A(m::Int)::Matrix{Float64}
    # Generate spatial matrix with vertical, horizontal, and diagonal neighbour interactions
    A = generate_spatial_neighbour_matrix(m, true)
    for i in axes(A, 1)
        for j in axes(A, 2)
            if i < j
                A[i, j] = A[i, j] * 1.5
            end
        end
    end
    return A
end

function design_D_generate_B(m::Int)::Matrix{Float64}
    # Generate diagonal matrix
    B = zeros(m^2, m^2)
    for i in axes(B, 1)
        B[i, i] = 0.2
    end
    return B
end

function run_simulation(p::Int, m::Int, T::Int, sim_design_id::String="sim", write::Bool=true, uuidtag::Union{String,Nothing}=nothing, train_size::Float64=0.8, h_A::Int=3, h_B::Int=3)::Matrix{Float64}
    T = Int(T / train_size)
    A = design_B_generate_A(m)
    B = design_B_generate_B(m)
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
    m = Int(sqrt(p))
    uuidtag = length(ARGS) >= 2 ? ARGS[2] : nothing

    # Run simulation
    run_simulation(p, m, T, sim_design_id, true, uuidtag)
end