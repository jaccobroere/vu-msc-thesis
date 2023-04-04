using CSV
using Tables
using LinearAlgebra
using DataFrames


# Read data
y = Matrix(CSV.read("out/sim_y.csv", DataFrame))
Vhat_d = Matrix(CSV.read("out/v2_Vhat_d.csv", DataFrame))
sigma_hat = Matrix(CSV.read("out/v2_sigma_hat.csv", DataFrame))

"""
Calculates the amount of active elemements, i.e. the number of columns of Vhat_d
"""
function total_active_elements(p::Int, h::Int)::Int
    return -2h^2 - 2h + p * (4h + 1)
end

total_active_elements(25, 6)

# Test the active elements selection logic
h = 2
p = 5

tot = 0
for i in 1:p
    lower_a = collect(max(1, i - h):max(0, i - 1))
    upper_a = collect((i+1):min(i + h, p))
    full_a = vcat(lower_a, upper_a)
    full_b = collect((p+max(1, (i - h))):(p+min(i + h, p)))
    selection = vcat(lower_a, upper_a, full_b)
    print("Length selection $(length(selection)) length full_a $(length(full_a)) length full_b $(length(full_b))\n")
    tot += length(selection)
end

print(tot)