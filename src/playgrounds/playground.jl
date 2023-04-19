using Base: vect
using CSV
using Tables
using LinearAlgebra
using DataFrames

"""
Calculates the amount of active elemements, i.e. the number of columns of Vhat_d
"""
function total_active_elements(p::Int, h::Int)::Int
    return p * (4h + 1) - 2(h^2 + h)
end

"""
Calculates the total number of elements in the vector c = vec(C')
"""
function tot_per_ci(p::Int, h::Int)::Int
    tot = 0
    for i in 1:p
        tot += 2h + 1 + 2 * min(abs(i - 1), h, abs(p - i))
    end
    return tot
end

"""
Calculates the number of nonzero elements per equation, that is the number of elements
in the subvector c_i of c = vec(C')
"""
function nonzero_elements_per_equation(i::Int, p::Int, h::Int)::Int
    return 2h + 1 + 2 * min(abs(i - 1), h, abs(p - i))
end


# total_active_elements(25, 6)

# Test the active elements selection logic
# h = 2
# p = 5

# tot = 0
# for i in 1:p
#     lower_a = collect(max(1, i - h):max(0, i - 1))
#     upper_a = collect((i+1):min(i + h, p))
#     full_a = vcat(lower_a, upper_a)
#     full_b = collect((p+max(1, (i - h))):(p+min(i + h, p)))
#     selection = vcat(lower_a, upper_a, full_b)
#     print("Length selection $(length(selection)) length full_a $(length(full_a)) length full_b $(length(full_b))\n")
#     tot += length(selection)
# end

# print(tot)

# tot_per_ci(25, 6)
# Read data
# y = Matrix(CSV.read("out/sim_y.csv", DataFrame))
# Vhat_d = Matrix(CSV.read("out/v2_Vhat_d.csv", DataFrame))
# sigma_hat = Matrix(CSV.read("out/v2_sigma_hat.csv", DataFrame))

using CSV, DataFrames, Tables

coef = Matrix(CSV.read("/Users/jacco/Documents/repos/vu-msc-thesis/out/simulation/coef/designB_T500_p25_admm_estimate.csv", DataFrame))
coef = vec(coef)

function coef_to_AB(coef::Vector{Float64}, p::Int)::Tuple{Matrix{Float64},Matrix{Float64}}
    rcoef = reverse(coef)
    bandwidth = div(p, 4)
    AB = zeros(p, 2p)
    cnt = 1
    for i in 1:p
        for j in 1:2p
            if j <= p # Belongs to A
                if abs(i - j) <= bandwidth && i != j
                    AB[i, j] = pop!(rcoef)
                end
            else # Belongs to B
                if abs(i - abs(j - p)) <= bandwidth
                    AB[i, j] = pop!(rcoef)
                end
            end
        end
    end
    return AB[:, 1:p], AB[:, (p+1):end]
end

A, B = coef_to_AB(coef, 25)

A