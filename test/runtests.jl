import HybridDynamics as HD
using Test

# Tolerance for testing the output of solution structs
soltol = 1e-1

@testset "HybridDynamics.jl" begin
    include("general.jl")
    include("solvers.jl")
end
