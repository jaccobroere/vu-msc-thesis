using CSV
using Tables



function main()
    p = 8749
    D = D_fusedlasso(p)

    CSV.write("out/D_fusedlasso_$(p).csv", Tables.table(D))
    return nothing
end

main()