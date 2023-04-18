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