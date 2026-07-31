import HybridDynamics as HD
using Test

# Tolerance for testing the output of solution structs
soltol = 0.01

@testset "HybridDynamics.jl" begin
    include("general.jl")
    include("solvers.jl")
end
