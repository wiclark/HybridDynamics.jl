module HybridDynamics

using LinearAlgebra
using ForwardDiff

include("Definitions.jl")

include("Systems/Linear_Affine.jl")
include("Systems/Lag_Ham.jl")
include("Systems/General.jl")
include("Systems/Filippov.jl")

include("ODE_solvers.jl")

#Basic Def Structs
export AbstractHybridSystem, AbstractHybridProblem

#Problem struct
export prob

########
# System types
####

# Filippov
export FilippovSys
# Lagrangian system struct
export LagSys
#Linear/Affine System/Problem structs
export LinearSystem, AffineSystem
#General System/Problem structs
export GeneralSystem, GeneralProblem

########
#ODE Step solvers - Interpolation with Fixed Step size
export solve, ForwardEuler, ModifiedTrap, ModifiedMidpoint, ExponentialSolver, RichardsonExtrapolation
#With Adaptive step size
export RK23, RK45

#ODE Step solvers - Extrapolation 
export AdamsBashforth2, AdamsBashforth3

#EventDetection locators
export LinearLocator, BisectionLocator, QuadraticLocator

#Linear/Affine additives
export beating_and_blocking_sets, is_trivially_blocking

end


"""
INTEGRATION GUIDE

This method of doing things uses Multiple Dispatch to seperate the physics from the solvers. 
To integrate a new system or solver, follow this guide that hopefully explains everything. 
Note: if you want more details on the specific things going on, go to the comments I have in each section as this will just be the overview

For examples on how to run anything try the Demo files within /test/...!

1. ADDING A NEW SYSTEM TYPE
To add a new system (e.g., 'MySystem.jl'), define a struct that subtypes 'AbstractHybridSystem' in a new file within '/src/Systems/'.

    • Step 1: Define the core system parameters and functions (vector field 'f', guard function 'h', and discrete map 'Δ', or whatever may be special to your system).
    • Step 2: Create the corresponding problem ('prob') and solution ('sol') structs for the new system type.
    • Step 3: Implement an 'init_solution' function that intializes whatever your custom solution struct with the initial condition requires.
    • Step 4: Construct the outer 'solve' function signature to accomodate your new system type if it requires certain wacky behavior. 

2) ADDING A NEW SOLVER
If you need to add a new integration method (e.g., RK4, etc):

WHERE: /src/ODE_solvers.jl
HOW: 
    1) Define a tag: 'struct MyNewSolver <: AbstractODESolver end'
    2) Implement the engine with the 'take_step' function using the quadratic guard crossing layout a mock is below.
        function take_step(::MyNewSolver, sys, f , xₖ, tₖ, Δt, tol)
            # 1) compute x_predict (the full step integration)
            x_predict = my_solver_step(f, xₖ, Δt, tₖ)

            # 2) Compute x_mid (the exact midpoint state at Δt / 2.0)
            x_mid = my_solver_step(f, xₖ, Δt / 2.0, tₖ)
            #NOTE: you may need to use another solver to do this step, particularly if you are using a Linear Multistep Method. 

            # 3) Evaluate the guard function at start, middle, and predicted end
            h_now = guard(sys, xₖ)
            h_mid = guard(sys, x_mid)
            h_next = guard(sys, x_predict)

            # 4) Evaluate event using quadratic interpolation method
            eventtriggered = crossed_guard(h_now, h_mid, h_next; tol=tol)

            # 5) Determine step size progression (adaptive or fixed) 
            #This will be improved in the future hopefully
            dt_next = Δt * 1.2

            # 6) Return the standard step package
            return x_predict, eventtriggered, h_now, dt_next
        end
        

3) ADDING A NEW EVENT LOCATOR (INTER/EXTRAPOLATION METHODS)
If you need a specialized root finding stragedy for event, or higher/lower order methods:

WHERE: /src/ODE_solvers.jl
HOW:
    1) Define a type tag: struct MyNewLocator <: AbstractEventLocator end
    2) Implement the engine with the 'locate_event' function to pinpoint root boundaries
        function locate_event(::MyNewLocator, sys, solver, f, xₖ, tₖ, Δt, h_now, tol, sol)
            # 1) Your Custom Root-Finding Strategy
            # Run your custom algorithm to isolate the fractional step size τ_star ∈ [0, Δt].

            # To evaluate the guard at any intermediate guess (τ) during your loop, query the solver engine like this:
            #   x_guess, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, τ, tol, sol)
            #   h_guess = guard(sys, x_guess)
    
    τ_star = # ... Your logic computes the precise fractional step here

    # 2) Calculate precise crossing variables at the isolated root
    t_star = tₖ + τ_star
    x_star, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, τ_star, tol, sol)
    
    # 3) Return the pinpointed impact conditions
    return t_star, x_star
end
## VARIABLE DICTIONARY

| Variable | Description |
| `sys` | The physical system object subtyping `AbstractHybridSystem`
| `f` | The continuous vector field function `(x, t) -> dx/dt`
| `xₖ` | The state vector at the start of the current step
| `x_mid` | The calculated state vector exactly halfway through the step time 
| `x_predict` | The predicted state vector at the end of the full step
| `tₖ` | The current simulation time
| `Δt` | The current time step duration
| `tol` | The numerical error/root-finding tolerance threshold
| `h_now` | Guard evaluation at the start of the step
| `h_mid` | Guard evaluation at the exact midpoint of the step (used for quadratic event detection)
| `h_next` | Guard evaluation at the predicted end of the step
| `eventtriggered` | Boolean indicating a confirmed linear crossing or quadratic dip
| `t_star` | The pinpointed time of guard impact
| `x_star` | The state vector on the guard surface at `t_star`
| `sol` | The complete solution object containing trajectory history
"""