module HybridDynamics

using LinearAlgebra
using ForwardDiff
# Needed for Stochastic Systems
using Random
using Distributions

include("ODE_Solvers/Definitions.jl")

include("Hermite_Interp.jl")
include("Analysis/Pathology.jl")

include("Systems/Linear_Affine.jl")
include("Systems/Mechanical.jl")
include("Systems/Nonholonomic.jl")
include("Systems/General.jl")
include("Systems/Filippov.jl")
include("Systems/Stochastic.jl")

include("ODE_Solvers/Linear_Multistep_Methods/Adaptive_LMM_Steps.jl")
include("ODE_Solvers/Linear_Multistep_Methods/Fixed_LMM_Steps.jl")
include("ODE_Solvers/Runge_Kutta_Methods/Fixed_RK_Steps.jl")
include("ODE_Solvers/Runge_Kutta_Methods/Adaptive_RK_Steps.jl")
include("ODE_Solvers/Stochastic_Methods/Fixed_SDE_Steps.jl")
include("ODE_Solvers/Bonus_Methods.jl")
include("ODE_Solvers/Event_Locators.jl")

include("Analysis/VariationalEquation.jl")
include("Analysis/Utilities.jl")




#Basic Def Structs
export AbstractHybridSystem, AbstractHybridProblem, AbstractHybridSolution

#Problem struct
export prob

########
# System types
####
export StochasticSystem
# Filippov
export FilippovSystem
# Mechanical system struct
export MechanicalSystem
# Nonholonomic system struct
export NonholonomicSystem
#Linear/Affine System/Problem structs
export LinearSystem, AffineSystem
#General System/Problem structs
export GeneralSystem

########
#ODE Step solvers - Interpolation with Fixed Step size
export solve, ForwardEuler, ModifiedTrap, ModifiedMidpoint, ExponentialSolver, RichardsonExtrapolationm, ImplicitEuler, RK4
#With Adaptive step size
export RK23, RK45

# Stochastic solver
export EulerMaruyama

#Fixed LMM
export AdamsBashforth2, AdamsBashforth3, BDF2
#Adaptive LMM
export AdaptiveABM2, AdaptiveABM3

#Extra Solvers
export MagnusLeapfrog

#EventDetection locators
export LinearLocator, BisectionLocator, QuadraticLocator, NewtonLocator

#Linear/Affine additives
export beating_and_blocking_sets, is_trivially_blocking

#Variational Equation
export variational_vector_field, compute_pushforward, apply_variational_jump

#Plotting Help
export split_jumps
end


"""
INTEGRATION GUIDE

### 1. ADDING A NEW SYSTEM TYPE
Create a new file in `/src/Systems/` (e.g., `MySystem.jl`).
* **Step 1:** Define your `struct` as a subtype of `AbstractHybridSystem`.
* **Step 2:** Define system-specific functions (e.g., drift `f`, diffusion `g`, constraint matrix `A`).
* **Step 3:** Create a constructor that handles optional parameters and uses `ForwardDiff` 
  to generate the guard normal if not provided.
* **Step 4:** Define a custom `Sol` struct (subtype `AbstractHybridSolution`) with an 
  initialization method.
* **Step 5:** Implement a `take_step_mysystem!` function to centralize logic 
  (e.g., sliding modes, Zeno points) by wrapping standard solver calls.
* **Step 6"** Implement your solve dispatch.

### 2. ADDING A NEW SOLVER
Numerical solvers are agnostic of the physical system type.
* **WHERE:** `/src/ODE_solvers.jl/...`
* **HOW:**
    1. Define a tag: `struct MyNewSolver <: AbstractODESolver end`.
    2. Implement `take_step(solver::MyNewSolver, prob, f, x, t, Δt, tol, sol; kwargs...)`.
    3. Return `x_predict`, `eventtrigger` (Bool), root-finding metrics, and `dt_next`.

    NOTE: Currently we use subsets of these solver types (Runge Kutta and Linear Multistep Methods)
    If yours fits within on of those then add it there. Otherwise add it as you see fit. 

### 3. ADDING A NEW EVENT LOCATOR
Locators handle root-finding strategies after a guard crossing is detected.
* **WHERE:** `/src/Event_Locators.jl`
* **HOW:**
    1. Define a tag: `struct MyNewLocator <: AbstractEventLocator end`.
    2. Implement `locate_event(::MyNewLocator, prob, solver, f, xₖ, tₖ, Δt, h_val, tol, sol, stepper)`.
    3. Isolate the fractional time τ* ∈ [0, Δt] where the impact occurs, 
       periodically querying the solver via `take_step` to probe the guard state.

### VARIABLE DICTIONARY

Below are some of the most common variables we use.

| Variable | Description |
| `sys` | The physical system object (subtype of `AbstractHybridSystem`) |
| `f` / `f_λ` | The continuous vector field (optionally including multipliers) |
| `xₖ` | The state vector at the start of the current step |
| `tₖ` | The current simulation time |
| `Δt` | The current time step duration |
| `tol` | Numerical tolerance for integration and event detection |
| `h` / `∇h` | Guard function and its gradient (surface normal) |
| `x_predict` | The predicted state at the end of the step |
| `t_star` | The pinpointed time of guard impact |
| `x_star` | The interpolated state vector exactly on the guard surface |
| `sol` | The `AbstractHybridSolution` object storing history/events |
"""