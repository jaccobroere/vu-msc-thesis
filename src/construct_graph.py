# Import relevant libraries
import sys

import igraph as ig
from julia import Main  # Julia API

Main.eval('include("src/utils.jl")')  # Import utils.jl


def create_edge(i: int, p: int, h: int) -> tuple:
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


def save_graph_as_gml(graph, p, h):
    graph.write(f"out/gsplash_graph_p{p}_h{h}.graphml")


if __name__ == "__main__":
    p = int(sys.argv[1])
    h = p // 4
    graph = create_gsplash_graph(p, h)
    save_graph_as_gml(graph, p, h)
