module HybridDynamics

using LinearAlgebra
using ForwardDiff

#The template for any solver method we use. I think this makes it nicer?
abstract type AbstractODESolver end

#Euler method tag and the others. I wont comment each as its pretty self explanatory
struct ForwardEuler <: AbstractODESolver end
struct ModifiedTrap <: AbstractODESolver end
struct ModifiedMidpoint <: AbstractODESolver end

abstract type AbstractEventLocator end
struct LinearLocator <: AbstractEventLocator end
struct BisectionLocator <: AbstractEventLocator end

include("ODE_solvers.jl")

include("Systems/Lag_Ham.jl")
include("Systems/Linear_Affine.jl")

export HybridLinearSystem, HybridAffineSystem, 
       HybridLinearProblem, HybridAffineProblem, 
       HybridLinearSolution, HybridAffineSolution,
       CreateSystem, CreateProblem, CreateSolution,
       solve_hybrid_system_exp,
       solve, solveloop, ForwardEuler, ModifiedTrap, ModifiedMidpoint, BisectionLocator, LinearLocator,
       forward_euler

end
