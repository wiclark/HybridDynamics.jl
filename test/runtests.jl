import HybridDynamics as HD
using Test

# Tolerance for testing the output of solution structs
soltol = 0.5

@testset "HybridDynamics.jl" begin
    include("general.jl")
    include("linear.jl")
    include("affine.jl")
    include("mechanical.jl")
    include("nonholonomic.jl")
    include("stochastic.jl")
    include("filippov.jl")

    include("solvers.jl")
end