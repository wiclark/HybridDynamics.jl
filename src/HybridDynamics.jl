module HybridDynamics

using LinearAlgebra
using ForwardDiff

export forward_euler, greet_your_package_name
include("LagHamDynamics.jl")
include("ODE_solvers.jl")
include("HybridSystemsDS.jl")
include("functions.jl")

export HybridLinearSystem, HybridAffineSystem, 
       HybridLinearProblem, HybridAffineProblem, 
       HybridLinearSolution, HybridAffineSolution,
       CreateSystem, CreateProblem, CreateSolution,
       solve_hybrid_system_exp,
       solve, solveloop, ForwardEuler, ModifiedTrap, ModifiedMidpoint, BisectionLocator

end
