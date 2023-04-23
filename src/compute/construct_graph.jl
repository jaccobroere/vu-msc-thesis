using Graphs, GraphIO, EzXML

path = dirname(abspath(@__FILE__))
include(joinpath(path, "utils.jl"))

"""
generate_idx_mapping(p::Int, bandwidth::Int) -> Dict

Generate a dictionary that maps the indices of two matrices A and B to the 
indices of a vector c. The mapping is determined by the size of the matrices 
(p) and the maximum distance allowed between the indices (bandwidth).

Arguments:
- `p::Int`: The size of the matrices A and B.
- `bandwidth::Int`: The maximum distance allowed between the indices.

Returns:
- A dictionary that maps a tuple of two indices (i, j) to an integer index in c.
"""
function generate_idx_mapping(p::Int, bandwidth::Int)::Dict{Tuple{Int,Int},Int}
    mapping_dict = Dict{Tuple{Int,Int},Int}()
    idx = 1
    for i = 1:p
        for j = 1:2p
            if j <= p # Belongs to A
                if abs(i - j) <= bandwidth && i != j
                    mapping_dict[(i, j)] = idx
                    idx += 1
                    continue
                end
            else # Belongs to B
                if abs(i - abs(j - p)) <= bandwidth
                    mapping_dict[(i, j)] = idx
                    idx += 1
                    continue
                end
            end
        end
    end
    return mapping_dict
end


"""
create_diagonal_edge(i, j, mapping_dict) -> Tuple or nothing

Create a diagonal edge in a graph based on the index mapping between two matrices.

Arguments:
- `i`: The row index in the matrices.
- `j`: The column index in the matrices.
- `mapping_dict`: A dictionary that maps a tuple of two indices (i, j) to an integer index in c.

Returns:
- A tuple representing a diagonal edge in the graph, or nothing if the edge cannot be created.
"""
function create_diagonal_edge(i, j, mapping_dict)
    current = get(mapping_dict, (i, j), -1)
    neighbor = get(mapping_dict, (i + 1, j + 1), -1)

    if current == -1 || neighbor == -1
        return nothing
    end

    return (current, neighbor)
end


"""
create_bidirectional_edge(i, j, mapping_dict) -> Tuple or nothing

Create a bidirectional edge in a graph based on the index mapping between two matrices.

Arguments:
- `i`: The row index in the matrices.
- `j`: The column index in the matrices.
- `mapping_dict`: A dictionary that maps a tuple of two indices (i, j) to an integer index in c.

Returns:
- A tuple representing a bidirectional edge in the graph, or nothing if the edge cannot be created.
"""
function create_bidirectional_edge(i, j, mapping_dict)
    current = get(mapping_dict, (i, j), -1)
    neighbor = get(mapping_dict, (j, i), -1)

    if current == -1 || neighbor == -1 || i < j
        return nothing
    end

    return (current, neighbor)
end




"""
Create a G-SPLASH graph for a given system of equations.

Args:
p (int): The dimension of the square matrix.
h (int): The bandwidth of the matrix.

Returns:
SimpleGraphs.SimpleGraph{Int}: A G-SPLASH graph representing the system of equations.
"""
function create_gsplash_graph2(p::Int, h::Int=0, bidirectional::Bool=false)::SimpleGraph{Int}
    if h == 0
        h = div(p, 4)
    end

    # Define the maximum numbers of vertices
    M::Int = total_number_nonzero_elements(p, h)
    graph = SimpleGraph(M)

    # Generate the index mapping
    mapping_dict = generate_idx_mapping(p, h)

    # Add edges
    for key in keys(mapping_dict)
        i, j = key

        # If bidirectional, add bidirectional edges as well
        if bidirectional
            edge = create_bidirectional_edge(i, j, mapping_dict)
            if edge !== nothing
                add_edge!(graph, edge)
            end
        end

        # Add diagonal edges, always 
        edge = create_diagonal_edge(i, j, mapping_dict)
        if edge !== nothing
            add_edge!(graph, edge)
        end
    end

    return graph
end

"""
Prints out a list of edges for the given graph.

Args:
g (AbstractGraph): The graph to print the edges of.
"""
function print_edges(g::AbstractGraph)
    for edge in edges(g)
        println("Edge $(edge.src) -> $(edge.dst)")
    end
end

"""
Prints the vertices of the given graph.

Args:
g (AbstractGraph{T}): The graph to print the vertices of.
"""
function print_vertices(g::AbstractGraph{T}) where {T}
    for v in vertices(g)
        println("Vertex $v")
    end
end


"""
Save a G-SPLASH graph to a file in GraphML format.

Args:
graph (SimpleGraphs.SimpleGraph): The G-SPLASH graph to save.
path (str, optional): The path to save the graph to. Defaults to "graph.graphml".
"""
function save_graph_as_gml(graph, path::AbstractString="graph.graphml")
    savegraph(path, graph, GraphIO.GraphMLFormat())
end

if abspath(PROGRAM_FILE) == @__FILE__
    p = parse(Int, ARGS[1])
    path_prefix = ARGS[2]
    h = div(p, 4)
    graph = create_gsplash_graph(p, h)
    save_graph_as_gml(graph, path_prefix)
end