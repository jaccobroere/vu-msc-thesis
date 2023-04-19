using Graphs, GraphIO, EzXML

path = dirname(abspath(@__FILE__))
include(joinpath(path, "utils.jl"))

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