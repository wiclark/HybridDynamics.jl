module HybridDynamics

using LinearAlgebra
using ForwardDiff

include("Definitions.jl")
include("Systems/Linear_Affine.jl")
include("Systems/Lag_Ham.jl")
include("Systems/General.jl")
include("ODE_solvers.jl")

#Basic Def Structs
export AbstractHybridSystem, AbstractHybridProblem

# Lagrangian system struct
export LagSys

#Linear/Affine System/Problem structs
export LinearSystem, AffineSystem
#Problem struct
export Problem

#General System/Problem structs
export GeneralSystem, GeneralProblem

#ODE Step solvers - Interpolation
export solve, ForwardEuler, ModifiedTrap, ModifiedMidpoint, ExponentialSolver

#ODE Step solvers - Extrapolation 
export AdamsBashforth2, AdamsBashforth3, RichardsonExtrapolation

#EventDetection locators
export LinearLocator, BisectionLocator

#Linear/Affine additives
export beating_and_blocking_sets, is_trivially_blocking

end


"""
INTEGRATION GUIDE

This method of doing things uses Multiple Dispatch to seperate the physics from the solvers. 
To integrate a new system or solver, follow this guide that hopefully explains everything. 
Note: if you want more details on the specific things going on, go to the comments I have in each section as this will just be the overview

1. ADDING A NEW SYSTEM TYPE
To add a new system (e.g., 'MySystem.jl'), define a struct that subtypes 'AbstractHybridSystem' in a new file within '/src/Systems/'.

Create the corresponding problem and solution structs for the new system. Then create createfunctions that construct the problem/solution structs from inputs. 
    Define a init_solution function that initializes the solution struct with the initial conditions.
        Construct the solve function that accomodates your system. Along with any secondary functions you need to define that may be special to your system. 

2) ADDING A NEW SOLVER
If you need to add a new integration method (e.g., RK4, etc):

WHERE: /src/ODE_solvers.jl
HOW: 
    1) Define a tag: 'struct MyNewSolver <: AbstractODESolver end'
    2) Implement the engine:
        function take_step(::MyNewSolver, sys, f, xₖ, tₖ, Δt, tol)
            1) Compute x_predict (the integration math)
            2) Compute h_now = guard(sys, xₖ)
            3) Compute h_next = guard(sys, x_predict)
            4) eventrigger = (signbit(h_now) != signbit(h_next))
            5) dt_next = Δt * 1.2 (or some other logic for adaptive stepping)
            6) Return (x_predict, eventtriggered, h_now, dt_next)
        

3) ADDING A NEW EVENT LOCATOR (INTER/EXTRAPOLATION METHODS)
If you need a specialized root finding stragedy for event, or higher/lower order methods:

WHERE: /src/ODE_solvers.jl
HOW: 
    1) Define a tag: 'struct MyNewLocator <: AbstractEventLocator end'
    2) Implement the search:
    function locate_event(::MyNewLocator, sys, f, xₖ, tₖ, Δt, h_now, tol)
        1) Define search interval: [0, Δt]
        2) Perform root-finding: (e.g., Bisection or Newton's) until range < tol
        3) Calculate t_star = tₖ + τ_star (where τ_star is the found crossing time)
        4) Calculate x_star = (the state at t_star via interpolation, extrapolation or integration)
        5) Return (t_star, x_star)
    end


VARIABLE DICTIONARY:
* sys:  The physical system object (e.g., LinearSystem)
* f:    The vector field function '(x,t) -> dx/dt
* xₖ:   The state vector (current position in state space)
* tₖ:   The current time in the simulation
* Δt:   The time step (duration of current step)
* tol:  The tolerance (error threshold for adaptive steps or root-finding)
* h_now:    The guard calue at the start of the step (sign indicates position relative to guard surface)
* h_next:   The guard value at the end of the predicted step
* eventtriggered: A Bool flag indicating if the guard was crossed. 
* t_star:   The precise time of impact/event
* x_star:   The precise state vector at the moment of impact  
"""