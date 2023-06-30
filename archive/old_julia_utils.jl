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

function fast_inv_Dtilde(p::Int, h::Int, graph::SimpleGraph)::SparseMatrixCSC
    # h == 0 ? h = div(p, 4) : h = h # Set h if not provided
    # Calculate the Dtilde matrix
    Dtilde = calc_Dtilde_sparse(graph)
    # Calculate the index of the first element before first skip
    iF = n_elements_to_first_skip(p, h) - 1

    # Index the blocks of Dtilde
    F = UpperTriangular(Dtilde[1:iF, 1:iF])
    G = Dtilde[1:iF, (iF+1):end]
    H = Dtilde[(iF+1):end, 1:iF]
    J = Dtilde[(iF+1):end, (iF+1):end]

    # Calculate the inverse of the blocks F and J
    F_inv = inv(lu(F))
    J_inv = inv(lu(J))

    # Calculate the inverse of Dtilde based on block inverses
    M_11 = inv(lu((F - G * J_inv * H)))
    M_22 = inv(lu((J - H * F_inv * G)))

    D_inv = [M_11 (M_11*-G*J_inv)
        (M_22*-H*F_inv) M_22]

    return D_inv
end

# TODO: Implement this using fast sparse solvers 
function model_fast_fusion(sigma_hat::Matrix{Float64}, Vhat_d::SparseMatrixCSC{Float64}, graph::SimpleGraph{Int64})
    # Calculate D_tilde by extending it with orthogonal rows to a square matrix
    Dtilde = calc_Dtilde(graph)
    m, p = ne(graph), nv(graph)

    # Use linear system solvers for faster computation of the change of variables (see Tibshirani and Taylor, 2011)
    XD1 = (Dtilde' \ Vhat_d')' # Same as Vhat_d * inv(Dtilde)
    X1, X2 = XD1[:, 1:m], XD1[:, (m+1):end]
    X2_plus = (X2' * X2) \ X2' # Same as inv(X2' * X2) * X2'

    # Transform the input to LASSO objective
    P = X2 * X2_plus
    ytilde = vec((I - P) * sigma_hat)
    Xtilde = (I - P) * X1

    # Solve LASSO
    path = glmnet(Xtilde, ytilde, intercept=false, lambda=[0.615848211066027], alpha=1, standardize=false)

    # Transform back to original variables
    theta1 = vec(path.betas)
    theta2 = X2_plus * (sigma_hat - X1 * theta1)
    theta = vcat(theta1, theta2)
    coef = Dtilde \ theta # Same as inv(Dtilde) * theta

    return coef
end

function entire_diagonal_penalties(D::SparseMatrixCSC, p::Int, h::Int)::Matrix{Int}
    n_diagonals = 4h + 1
    diagonal_collection = generate_diagonal_collection(p, h)
    penalties = spzeros(Int, n_diagonals, size(D, 2))
    sorted_keys = sort(collect(keys(diagonal_collection)))
    for (i, key) in enumerate(sorted_keys)
        value = diagonal_collection[key]
        setindex!(penalties, ones(Int, length(value)), i, value)
    end

    return penalties
end

function calc_Dtilde_SDF_sparse(graph::SimpleGraph, p::Int, h::Int=0)::SparseMatrixCSC{Float64}
    D = sparse(incidence_matrix(graph, oriented=true)') # Incidence matrix needs to be transposed before obtaining D^(G)
    EX = entire_diagonal_penalties(D, p, h)
    return vcat(D, EX)
end

function inv_Dtilde_SDF_sparse(graph::SimpleGraph, p::Int, h::Int=0)::SparseMatrixCSC{Float64}
    Dtilde = calc_Dtilde_SDF_sparse(graph, p, h)
    return sparse(inv(lu(Dtilde)))
end

