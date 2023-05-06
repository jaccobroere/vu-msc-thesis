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

"""
This function calculates the number of the first element that has 
no fusion partner in the Fused Lasso penalty associated with it. That is, it is
the number of elements in the first (p - h) rows that are put into c, minus 2h + 1,
which are the number of elements that come after this element in the (p -h)th row.
"""
function n_elements_to_first_skip(p::Int, h::Int)::Int
    return (-5h^2 - 2h + 4h * p + p) - (2h + 1)
end



p = 5
h = 2



tot = 0
for i in 1:p
    tot += nonzero_elements_per_equation(i, p, h)
    # Print i along with the tot
    println("$i: $tot")
    # println(nonzero_elements_per_equation(i, p, h))
end



function sum_elements_row(p::Int, h::Int)::Int
    tot = 0
    for i in 1:(p-h)
        tot += nonzero_elements_per_equation(i, p, h)
    end
    return tot
end

sum_elements_row(p, h)