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

function design_B_generate_A(width::Int, height::Int)::Matrix{Float64}
    p = width * height
    A = zeros(p, p)
    for i in 1:p
        for j in 1:p
            if abs(i - j) == width
                setindex!(A, 0.2, i, j)
            end

            if abs(i - j) == 1 && (((i % width == 0) && ((j - 1) % width == 0)) || ((j % width == 0) && ((i - 1) % width == 0)))
                continue
            end

            if abs(i - j) == 1
                setindex!(A, 0.2, i, j)
            end
        end
    end
    return A
end

function design_B_generate_B(width::Int, height::Int, dimension_dict::Union{Nothing,Dict}=nothing)::Matrix{Float64}
    if isnothing(dimension_dict)
        dimension_dict = Dict(
            4 => 0.28,
            5 => 0.25,
            6 => 0.23
        )
    end

    p = width * height
    B = zeros(p, p)
    for i in 1:p
        setindex!(B, dimension_dict[(width, height)], i, i)
    end
    return B
end

function run_simulation(width::Int, height::Int, T::Int, dimension_dict::Union{Nothing,Dict}=nothing, write::Bool=true, train_size::Float64=0.8, h_A::Int=3, h_B::Int=3)::Matrix{Float64}
    T = Int(T / train_size)
    A = design_B_generate_A(width, height)
    B = design_B_generate_B(width, height, dimension_dict)
    errors = generate_errors_over_time(T, width * height)
    y = simulate_svar(A, B, errors)

    sim_design_id = "sim_runtime_p$(width * height)_T$(Int(T * train_size))"

    if write == true
        path = joinpath("data", "simulation", "runtime", sim_design_id)

        if !isdir(path)
            mkpath(path)
        end

        CSV.write(joinpath(path, "A.csv"), DataFrame(A, :auto))
        CSV.write(joinpath(path, "B.csv"), DataFrame(B, :auto))
        CSV.write(joinpath(path, "y.csv"), DataFrame(y, :auto))
    end

    run(`julia src/compute/jl/precalculations_and_write.jl runtime/$sim_design_id`)
    return y
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line argument
    uuidtag = length(ARGS) >= 2 ? ARGS[2] : nothing

    # Dimension dictionary
    dimension_dict = Dict(
        (4, 4) => 0.28,
        (4, 5) => 0.28,
        (5, 5) => 0.25,
        (5, 6) => 0.25,
        (6, 6) => 0.23,
        (6, 7) => 0.23,
        (7, 7) => 0.22,
        (7, 8) => 0.22,
        (6, 10) => 0.21,
        (8, 8) => 0.21,
        (8, 9) => 0.20,
        (9, 9) => 0.20,
        (9, 10) => 0.19,
        (10, 10) => 0.19
    )

    T = 1000

    # Run simulation for every dimension in dimension_dict
    for (width, height) in keys(dimension_dict)
        run_simulation(width, height, T, dimension_dict, true)
    end
end