path = dirname(abspath(@__FILE__))
include(joinpath(path, "utils.jl"))
include(joinpath(path, "construct_graph.jl"))
include(joinpath(path, "precalculations_and_write.jl.jl"))
include(joinpath("../simulation/", "simulation_utils.jl"))
include(joinpath("../simulation/", "simulation_designA.jl"))
include(joinpath("../simulation/", "simulation_designB.jl"))

