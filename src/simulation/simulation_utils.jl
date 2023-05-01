using LinearAlgebra
using Random
using DataFrames
using CSV

"""
generate_random_vector(p)

Return a random vector of length p with values drawn from a standard normal distribution.

# Arguments
- `p::Int`: The length of the vector to generate.

# Returns
- An array of `p` random floating-point numbers drawn from a standard normal distribution.
"""
function generate_random_vector(p::Int)::Array{Float64,1}
    return rand(Normal(0, 1), p)
end

"""
generate_errors_over_time(T, p)

Return an array of error vectors of length `p` for each time period `T`.

# Arguments
- `T::Int`: The number of time periods to simulate.
- `p::Int`: The length of the error vectors.

# Returns
- An array of length `T`, where each element is an array of `p` random floating-point numbers drawn from a standard normal distribution.
"""
function generate_errors_over_time(T::Int, p::Int)::Array{Array{Float64,1},1}
    errors = [generate_random_vector(p) for _ in 1:T]
    return errors
end

"""
simulate_svar(A, B, errors)

Simulate a Structural Vector Autoregression (SVAR) model.

# Arguments
- `A::Matrix{Float64}`: A matrix representing the structural coefficients of the model.
- `B::Matrix{Float64}`: A matrix representing the impact of the errors on the variables.
- `errors::Array{Array{Float64,1},1}`: An array of error vectors for each time period.

# Returns
- A matrix where each row represents the values of a variable at each time period.
"""
function simulate_svar(A::Matrix{Float64}, B::Matrix{Float64}, errors::Array{Array{Float64,1},1})::Matrix{Float64}
    p = size(A, 1)
    T = length(errors)
    y = zeros(p, T)
    C = inv(I - A) * B
    for t in 1:T
        if t == 1
            y[:, t] = errors[t]
            continue
        end

        y[:, t] = C * y[:, t-1] + errors[t]
    end

    return y
end

"""
Write simulation outputs to CSV files.

# Arguments
- `A::Matrix{Float64}`: Matrix A from the SVAR model.
- `B::Matrix{Float64}`: Matrix B from the SVAR model.
- `y::Matrix{Float64}`: Matrix of simulated data from the SVAR model.
- `file_prefix::String`: Prefix to use for the file names.
- `uuidtag::Union{String, Nothing}`: UUID tag to use in the directory name (optional).

# Returns
- None
"""
function write_simulation_output(A::Matrix{Float64}, B::Matrix{Float64}, y::Matrix{Float64}, sim_design_id::String, uuidtag::Union{String,Nothing}=nothing)::Nothing
    if uuidtag !== nothing
        path = joinpath("data", "simulation", sim_design_id, uuidtag)
        if !isdir(path)
            mkpath(path)
        end
    else
        path = joinpath("data", "simulation", sim_design_id)
    end

    CSV.write(joinpath(path, "A.csv"), DataFrame(A, :auto))
    CSV.write(joinpath(path, "B.csv"), DataFrame(B, :auto))
    CSV.write(joinpath(path, "y.csv"), DataFrame(y, :auto))
    return nothing
end

