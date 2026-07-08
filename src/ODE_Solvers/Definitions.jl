#---------------------------
#ABSTRACT TYPES
#These are the empty categories that hold nothing themselves. By restricting our solver function args to these types, we can make sure the solver can accept ANY system or solver that is in these families.

#Parent Abstract Type
abstract type AbstractODESolver end
abstract type RK <: AbstractODESolver end
abstract type LMM <: AbstractODESolver end
abstract type STO <: AbstractODESolver end
abstract type ME <: AbstractODESolver end
abstract type NH <: AbstractODESolver end

#Parent Category for any root-finding algorithm used to pinpoint the exact time of a guard crossing.
abstract type AbstractEventLocator end

#Parent Category for any physical system that features both continuous flow and discrete jumps 
abstract type AbstractHybridSystem end
abstract type AbstractHybridProblem end

#Parent Category for Solution types.
abstract type AbstractHybridSolution end

struct prob{F <: AbstractHybridSystem, I  <: AbstractArray{Float64}, T <: Tuple{Float64, Float64}} <: AbstractHybridProblem
    sys::F
    init::I
    tspan::T
end

"""
    solve(prob; kwargs...)

Solve a hybrid dynamical system.

The available keyword arguments depend on the type of system stored in
`prob.sys`. See the method-specific documentation for details.

Supported systems include:

- `GeneralSystem`
- `MechanicalSystem`
- `FilippovSys`
- ...

Examples
--------
solve(prob, RK4())
solve(prob, RK4(); dense_output=true)
"""
function solve end