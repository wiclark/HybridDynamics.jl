
#This file is to establish the foundation of everything we use with Multiple Dispatch. It doesnt need to be here exactly but I think it will be a good place to keep it organized. 

#---------------------------
#ABSTRACT TYPES
#These are the empty categories that hold nothing themselves. By restricting our solver function args to these types, we can make sure the solver can accept ANY system or solver that is in these families.

#Parent Abstract Type
abstract type AbstractODESolver end

#Parent Category for any root-finding algorithm used to pinpoint the exact time of a guard crossing.
abstract type AbstractEventLocator end

#Parent Category for any physical system that features both continuous flow and discrete jumps (unsure how Filippov will integrate but I hope its not too bad)
abstract type AbstractHybridSystem end

abstract type AbstractHybridProblem end

#Parent Category for Solution types.
abstract type AbstractHybridSolution end

# General problem
struct prob{F <: AbstractHybridSystem, I  <: AbstractVector{Float64}, T <: Tuple{Float64, Float64}} <: AbstractHybridProblem
    sys::F
    init::I
    tspan::T
end


