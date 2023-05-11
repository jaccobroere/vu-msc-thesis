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


include(joinpath("..", "compute", "construct_graph.jl"))


graph = create_gsplash_graph(5)
sym_graph = create_gsplash_graph(9, symmetric=true)

using Graphs, Colors, Plots, LinearAlgebra

D = incidence_matrix(graph, oriented=true)
D = Matrix(D)




null_vecs

null_vecs' * D


s

A = nullspace(D)

Dtilde = vcat(D, A')
inv(Dtilde)

function compute_null_space(graph::SimpleGraph)
    conn = connected_components(graph)
    p = nv(graph)



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

Q = qr_transform_genlasso(graph)




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

graph = create_gsplash_graph(5)
null = null_space_graph(graph)

Dtilde = construct_Dtilde(D, null)


