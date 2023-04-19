function D_fusedlasso(p::Int)::Matrix{Float64}
    D = zeros(p, p)
    for i in 1:p
        D[i, i] = 1
        if i < p
            D[i, i+1] = -1
        end
    end
    return D
end

function calc_sigma_e(e::Vector{Vector{Float64}})::Matrix{Float64}
    sigma_e = zeros(size(e[1], 1), size(e[1], 1))
    emean = mean(e)
    for i in 1:length(e)
        sigma_e += (e[i] .- emean) * (e[i] .- emean)'
    end
    return sigma_e / (length(e) - 1)
end