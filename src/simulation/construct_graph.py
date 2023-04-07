import sys

import igraph as ig
from julia import Main

Main.eval('include("src/simulation/utils.jl")')


def create_edge(i, p, h):
    """
    Create an edge between two vertices in the graph.

    Parameters
    ----------
    i : int
        The index of the first vertex.
    p : int
        The dimension of the square matrix.
    h : int
        The bandwidth of the matrix.

    Returns
    -------
    tuple
        A tuple of two integers representing the indices of the two vertices connected by the edge.
    """
    if Main.is_edge_case(i, p, h):
        return None

    row_num = Main.row_number_of_element(i, p, h)
    n_elements_this_row = Main.nonzero_elements_per_equation(row_num, p, h)
    delta_elements = (
        Main.nonzero_elements_per_equation(row_num + 1, p, h) - n_elements_this_row
    )

    if Main.belongs_to_A(i, p, h):
        if delta_elements > 0:
            edge = (i, i + n_elements_this_row + delta_elements - 1)
        elif delta_elements == 0:
            edge = (i, i + n_elements_this_row + delta_elements)
        else:  # delta_elements < 0
            edge = (i, i + n_elements_this_row + delta_elements + 2)
    else:  # Belongs to B
        if delta_elements >= 0:
            edge = (i, i + n_elements_this_row + delta_elements)
        else:  # delta_elements < 0
            edge = (i, i + n_elements_this_row + delta_elements + 1)

    return edge[0] - 1, edge[1] - 1


def create_gsplash_graph(p, h):
    """
    Create a G-SPLASH graph for a given system of equations.

    Parameters
    ----------
    p : int
        The dimension of the square matrix.
    h : int
        The bandwidth of the matrix.

    Returns
    -------
    ig.Graph
        A G-SPLASH graph representing the system of equations.
    """
    # Define the maximum numbers of vertices
    M = Main.total_number_nonzero_elements(p, h)
    graph = ig.Graph()

    # Add vertices
    graph.add_vertices([i for i in range(1, M + 1)])

    # Add edges
    edges = [create_edge(i, p, h) for i in range(1, M + 1)]
    edges = [edge for edge in edges if edge is not None]
    graph.add_edges(edges)

    return graph


def save_graph_as_gml(graph, path_prefix="out/"):
    """
    Save a G-SPLASH graph to a file in GraphML format.

    Parameters
    ----------
    graph : ig.Graph
        The G-SPLASH graph to save.
    p : int
        The dimension of the square matrix.
    h : int
        The bandwidth of the matrix.
    """
    graph.write(f"out/{path_prefix}_graph.graphml")


if __name__ == "__main__":
    p = int(sys.argv[1])
    path_prefix = sys.argv[2]
    h = p // 4
    graph = create_gsplash_graph(p, h)
    save_graph_as_gml(graph, path_prefix=path_prefix)
