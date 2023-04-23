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

"""
Create an edge between two vertices in the graph.

Args:
i (int): The index of the first vertex.
p (int): The dimension of the square matrix.
h (int): The bandwidth of the matrix.

Returns:
Union{Tuple{Int, Int}, Nothing}: A tuple of two integers representing the indices of the two vertices connected by the edge. If an edge cannot be created between the vertices, returns nothing.
"""
function create_edge(i::Int, p::Int, h::Int)::Union{Tuple{Int,Int},Nothing}
    if is_edge_case(i, p, h)
        return nothing
    end

    row_num = row_number_of_element(i, p, h)
    n_elements_this_row = nonzero_elements_per_equation(row_num, p, h)
    n_elements_next_row = nonzero_elements_per_equation(row_num + 1, p, h)
    delta_elements = n_elements_next_row - n_elements_this_row

    if belongs_to_A(i, p, h)
        if delta_elements > 0
            correction = -1
        elseif delta_elements == 0
            correction = 0
        else
            correction = 2
        end
    else  # Belongs to B
        if delta_elements < 0
            correction = 1
        else
            correction = 0
        end
    end

    return (i, i + n_elements_this_row + delta_elements + correction)
end

"""
Create a G-SPLASH graph for a given system of equations.

Args:
p (int): The dimension of the square matrix.
h (int): The bandwidth of the matrix.

Returns:
SimpleGraphs.SimpleGraph{Int}: A G-SPLASH graph representing the system of equations.
"""
function create_gsplash_graph(p::Int, h::Int=0)::SimpleGraph{Int}
    if h == 0
        h = div(p, 4)
    end

    # Define the maximum numbers of vertices
    M::Int = total_number_nonzero_elements(p, h)
    graph = SimpleGraph(M)

    # Add edges
    for i::Int in 1:M
        edge = create_edge(i, p, h)
        if edge !== nothing
            add_edge!(graph, edge)
        end
    end

    return graph
end