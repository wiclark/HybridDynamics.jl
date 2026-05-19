module HybridDynamics

using LinearAlgebra
using ForwardDiff

export forward_euler, greet_your_package_name
include("LagHamDynamics.jl")
include("ODE_solvers.jl")
include("HybridSystemsDS.jl")
include("functions.jl")

export HybridLinearSystem, HybridAffineSystem
export HybridLinearProblem, HybridAffineProblem
export solve_hybrid_system, solve_hybrid_system_exp
export beating_and_blocking_sets, is_trivially_blocking, basis_beating_and_blocking_sets
export forward_euler_step, rk_23_step, modified_midpoint_step, modified_trap_step

end
