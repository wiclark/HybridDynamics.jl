
#This file is to establish the foundation of everything we use with Multiple Dispatch. It doesnt need to be here exactly but I think it will be a good place to keep it organized. 

#---------------------------
#ABSTRACT TYPES
#These are the empty categories that hold nothing themselves. By restricting our solver function args to these types, we can make sure the solver can accept ANY system or solver that is in these families.

#Parent category for any continuous ODE integration method
abstract type AbstractODESolver end

#Parent Category for any root-finding algorithm used to pinpoint the exact time of a guard crossing.
abstract type AbstractEventLocator end

#Parent Category for any physical system that features both continuous flow and discrete jumps (unsure how Filippov will integrate but I hope its not too bad)
abstract type AbstractHybridSystem end

abstract type AbstractHybridProblem end

#Parent Category for Solution types.
abstract type AbstractHybridSolution end

#-----------------
#SOLVE TAGS
#These are empty structs. Their purpose is to act as a tag for the Julia compiler. When a user passes 'ModifiedMidpoint()' to the solver, the compiler uses multiple dispatch to route to the 'take_step' function written for that tag. 

#Tag for standard Forward Euler integration 
struct ForwardEuler <: AbstractODESolver end

#Tag for the Modifed Trapezoidal method 
struct ModifiedTrap <: AbstractODESolver end

#Tag for the Modifed Midpoint method 
struct ModifiedMidpoint <: AbstractODESolver end

#Tag for exact integration of Linear systems via the matrix Exponential 
struct ExponentialSolver <: AbstractODESolver end

#Extrapolation struct tags

#Richardson Extrapolation tag
struct RichardsonExtrapolation <: AbstractODESolver end

struct AdamsBashforth2 <: AbstractODESolver end

struct AdamsBashforth3 <: AbstractODESolver end

#-----------------------
#LOCATOR TAGS
#These are similar to the solver tags. But these define how the solver finds the exact crossing time once an event is detected. 

#Tag to use Linear Interpolation. Very fast but can be innacurate for higher order methods. 
struct LinearLocator <: AbstractEventLocator end

#Tag to use a bisection method serach. Can be very accurate but also very slow with complex systems
struct BisectionLocator <: AbstractEventLocator end

#Tag for quadratic event locator
struct QuadraticLocator <: AbstractEventLocator end

# General problem
struct prob{F <: AbstractHybridSystem, I  <: AbstractVector{Float64}, T <: Tuple{Float64, Float64}} <: AbstractHybridProblem
    sys::F
    init::I
    tspan::T
end

#Internal for now
#Check Zeno function for now:
function check_zeno(jump_interval, last_jump_interval, zeno_count, t_star, tol, zeno_ratio, max_zeno_jumps)
    status =:continue
    if  jump_interval < last_jump_interval * zeno_ratio #contraction checking
        zeno_count += 1
        if zeno_count >= max_zeno_jumps
            @warn "Zeno Behavior Detected: System reached max zeno jumps ($max_zeno_jumps)"
            status =:terminate
        elseif jump_interval <= tol && zeno_count > 7
            @info "Zeno Accumulation Point Reached at t = $t_star. Terminating."
            status=:terminate
        end
    else
        zeno_count = 0
    end
    return zeno_count, status
end
