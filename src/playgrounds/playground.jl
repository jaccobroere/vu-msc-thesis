using Graphs: SimpleGraphs, connected_components!
using LinearAlgebra
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
Calculates the total number of edges in a non-symmetric GF-SPLASH graph,
namely the total number of nonzero elements minus (4 * bandwith + 1)
"""
function total_number_edges(p::Int, h::Int)::Int
    return p * (4h + 1) - 2h^2 - 6h - 1
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



function transform_genlasso(graph::SimpleGraph{Int64})::Matrix{Float64}
    D = Matrix(incidence_matrix(graph, oriented=true))
    Dtilde = vcat(D, nullspace(D)')

    return inv(Dtilde)
end


function qr_transform_genlasso(graph::SimpleGraph{Int64})
    D = Matrix(incidence_matrix(graph, oriented=true))
    p, m = size(D) # m < p in pure fusion case
    QR = qr(D)

    return Matrix(QR.Q)
end

function null_space_graph(graph::SimpleGraph{Int64})::Matrix{Int64}
    null_vecs = zeros(Int64, nv(graph), nv(graph) - ne(graph))
    conn = connected_components(graph)
    for (j, vec) in enumerate(conn)
        for i in vec
            null_vecs[i, j] = 1
        end
    end
    return null_vecs
end

function construct_Dtilde(graph::SimpleGraph)::Matrix{Int64}
    null = null_space_graph(graph)
    return hcat(D, null)
end

function generate_diagonal_collection(p::Int, bandwidth::Int)::Dict{Tuple{Int,Char},Vector{Int}}
    diagonal_collection = Dict{Tuple{Int,Char},Vector{Int}}()
    idx = 1
    for i = 1:p
        for j = 1:2p
            if j <= p # Belongs to A
                if abs(i - j) <= bandwidth && i != j
                    try
                        push!(diagonal_collection[(i - j, 'A')], idx)
                    catch
                        diagonal_collection[(i - j, 'A')] = [idx]
                    end
                    idx += 1
                    continue
                end
            else # Belongs to B
                if abs(i - abs(j - p)) <= bandwidth
                    try
                        push!(diagonal_collection[(i - j, 'B')], idx)
                    catch
                        diagonal_collection[(i - j, 'B')] = [idx]
                    end
                    idx += 1
                    continue
                end
            end
        end
    end
    return diagonal_collection
end

function entire_diagonal_penalties(D::Matrix{Int64}, p::Int, h::Int)::Matrix{Int}
    n_diagonals = 4h + 1
    diagonal_collection = generate_diagonal_collection(p, h)
    penalties = spzeros(Int, n_diagonals, size(D, 2))
    for (i, (key, value)) in enumerate(diagonal_collection)
        setindex!(penalties, ones(Int, length(value)), i, value)
    end

    return penalties
end

function calc_Dtilde_SDF(D::Matrix{Int}, p::Int, h::Int)::Matrix{Int}
    return vcat(D, entire_diagonal_penalties(D, p, h))
end


function Dtilde_per_p(p::Int)
    h = div(p - 1, 4)
    graph = create_gsplash_graph(p)
    D = Matrix(incidence_matrix(graph, oriented=true)')
    n_diagonals = 4h + 1

    diagonal_collection = generate_diagonal_collection(p, h)
    EX = entire_diagonal_penalties(D, p, h)
    Dtilde = calc_Dtilde_SDF(D, p, h)

    try
        res = inv(Dtilde)
        print("SUCCESS")
        return res
    catch
        println("p = $p is not invertible")
        return nothing
    end
end


for p in 95:100
    Dtilde = Dtilde_per_p(p)
end