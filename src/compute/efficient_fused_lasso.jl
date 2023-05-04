using Graphs: SimpleGraphs, connected_components!
using LinearAlgebra
using Base: vect
using SparseArrays

include(joinpath("construct_graph.jl"))

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

