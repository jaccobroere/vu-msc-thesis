"""
Calculates the amount of active elemements, i.e. the number of columns of Vhat_d, or the number of elements in vec(C')
"""
function total_number_nonzero_elements(p::Int, h::Int)::Int
    return p * (4h + 1) - 2(h^2 + h)
end

"""
Calculates the total number of elements in the vector c = vec(C')
"""
function total_number_nonzero_elements_check(p::Int, h::Int)::Int
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
function nonzero_elements_per_equation(j::Int, p::Int, h::Int)::Int
    return 2h + 1 + 2 * min(abs(j - 1), h, abs(p - j))
end

"""
This function calculates the number of elemens per equation, that are in the matrix A 
"""
function nonzero_elements_per_equation_A(i::Int, p::Int, h::Int)::Int
    return h + min(abs(i - 1), h, abs(p - i))
end

"""
This function calculates the number of elemens per equation, that are in the matrix B 
"""
function nonzero_elements_per_equation_B(i::Int, p::Int, h::Int)::Int
    return h + 1 + min(abs(i - 1), h, abs(p - i))
end

"""
Checks whether an element belongs to A (true) or B (false)
"""
function belongs_to_A(i::Int, p::Int, h::Int)::Bool
    j = row_number_of_element(i, p, h)
    lower_bound = sum(nonzero_elements_per_equation.(1:(j-1), p, h))
    upper_bound = lower_bound + nonzero_elements_per_equation_A(j, p, h)
    return (lower_bound < i <= upper_bound)
end

"""
Determines the row number the element is on
"""
function row_number_of_element(i::Int, p::Int, h::Int)::Int
    if i > total_number_nonzero_elements(p, h) || i < 1
        error("Index out of bounds")
    end

    cur = 0
    row = 0
    while i > cur
        row += 1
        cur += nonzero_elements_per_equation(row, p, h)
    end
    return row
end


"""
Determines whether the element is on the bottom row of A or Base
"""
function is_on_bottom_row(i::Int, p::Int, h::Int)::Bool
    j = row_number_of_element(i, p, h)
    return (j == p)
end

"""
Determines whether the element is on the last column of A or B
"""
function is_on_last_column(i::Int, p::Int, h::Int)::Bool
    if i > total_number_nonzero_elements(p, h) || i < 1
        error("Index out of bounds")
    end

    j = row_number_of_element(i, p, h)
    edge_A = sum(nonzero_elements_per_equation.(1:(j-1), p, h)) + nonzero_elements_per_equation_A(j, p, h)
    edge_B = sum(nonzero_elements_per_equation.(1:j, p, h))

    if (h + j) >= p && (i == edge_A || i == edge_B)
        return true
    else
        return false
    end
end

"""
Determines whether the element is an edge case to be excluded
"""
function is_edge_case(i::Int, p::Int, h::Int)::Bool
    if is_on_bottom_row(i, p, h) || is_on_last_column(i, p, h)
        return true
    else
        return false
    end
end

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
            # mapping_dict[(i, j)] = -1
        end
    end
    return mapping_dict
end

function create_diagonal_edge(i, j, mapping_dict)
    current = get(mapping_dict, (i, j), -1)
    neighbor = get(mapping_dict, (i + 1, j + 1), -1)

    if current == -1 || neighbor == -1
        return nothing
    end

    return (current, neighbor)
end

function create_bidirectional_edge(i, j, mapping_dict)
    current = get(mapping_dict, (i, j), -1)
    neighbor = get(mapping_dict, (j, i), -1)

    if current == -1 || neighbor == -1 || i < j
        return nothing
    end

    return (current, neighbor)
end
