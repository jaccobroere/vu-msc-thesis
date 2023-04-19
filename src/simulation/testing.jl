using Test

# include(joinpath(dirname(abspath(@__FILE__)), "simulation_designA.jl"))
include(joinpath(dirname(abspath(@__FILE__)), "simulation_designB.jl"))
# include(joinpath(dirname(abspath(@__FILE__)), "simulation_designC.jl"))

generate_A(5)