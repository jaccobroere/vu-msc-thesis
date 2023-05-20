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

function generate_spatial_neighbour_matrix(m::Int, zero_diagonal::Bool)::Matrix{Float64}
    p = m^2
    A = zeros(p, p)
    for i in 1:p
        for j in 1:p
            # Skip the diagonal for matrix A
            if i == j && zero_diagonal == true
                continue
            end

            # Vertical neighbours
            if abs(i - j) == m
                setindex!(A, 0.2, i, j)
            end

            # Horizontal neighbours
            not_side_edge = !((i % m == 0) && ((j - 1) % m == 0) || ((j % m == 0) && ((i - 1) % m == 0)))
            if abs(i - j) == 1 && not_side_edge
                setindex!(A, 0.2, i, j)
            end

            # Diagonal neighbours
            # First create reusable conditions 
            not_top_edge = (i > m) && (j > m)
            not_bottom_edge = (i <= (p - m)) && (j <= (p - m))
            diag_value = 0.15
            # Bottom right diagonal neighbour
            if (not_side_edge && not_bottom_edge)
                if abs(i - j) == (m + 1)
                    setindex!(A, diag_value, i, j)
                end
            end
            # Top right diagonal neighbour
            if (not_side_edge && not_top_edge)
                if abs(i - j) == (m - 1)
                    setindex!(A, diag_value, i, j)
                end
            end
            # Bottom left diagonal neighbour
            if (not_bottom_edge && not_side_edge)
                if abs(i - j) == (m - 1)
                    setindex!(A, diag_value, i, j)
                end
            end
            # Top left diagonal neighbour
            if (not_top_edge && not_side_edge)
                if abs(i - j) == (m + 1)
                    setindex!(A, diag_value, i, j)
                end
            end
        end
    end
    return A
end

function design_D_generate_A(m::Int)::Matrix{Float64}
    # Generate spatial matrix with vertical, horizontal, and diagonal neighbour interactions
    A = generate_spatial_neighbour_matrix(m, true)
    for i in axes(A, 1)
        for j in axes(A, 2)
            if i < j
                A[i, j] = A[i, j] / 3
            end
        end
    end
    return A
end

function design_D_generate_B(m::Int)::Matrix{Float64}
    B_vals = Dict(
        4 => 0.3,
        5 => 0.25,
        6 => 0.22,
    )
    p = m^2
    # Generate diagonal matrix
    B = zeros(p, p)
    for i in axes(B, 1)
        B[i, i] = B_vals[m]
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