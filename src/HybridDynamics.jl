module HybridDynamics

# Write your package code here.

export forward_euler, greet_your_package_name
include("functions.jl")
include("LagrangianDynamics.jl")
include("ODE_solvers.jl")
end
