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

#LOCATOR TAGS Temp Location while I figure out where they can go in the include load order
#Tag to use Linear Interpolation that uses quadratic interpolation if linear missed the event. 
struct LinearQuadratic <: AbstractEventLocator end

#Tage to use linear interpolation with Hermite interpolation if linear missed the event.
struct LinearHermite <: AbstractEventLocator end

#Tag to use JUST linear interpolation
struct Linear <: AbstractEventLocator end


struct prob{S <: AbstractHybridSystem, I  <: AbstractArray{Float64}, T <: Tuple{Float64, Float64}} <: AbstractHybridProblem
    sys::S
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
- `LinearSystem`
- `AffineSystem`
- `MechanicalSystem`
- `NonholonomicSystem`
- `StochasticSystem`
- `FilippovSystem`

Examples
--------
    solve(prob, RK4())
"""
function solve end