
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

######
### WC: I would break down these structures more.
######
#=
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
=#
struct RK <: AbstractODESolver end
struct ForwardEuler <: RK end
struct ModifiedTrap <: RK end
struct ModifiedMidpoint <: RK end
struct RichardsonExtrapolation <: RK end # I think so, at least

struct ExponentialSolver <: AbstractODESolver end

struct LMM <: AbstractODESolver end
struct AdamsBashforth2 <: LMM end
struct AdamsBashforth3 <: LMM end

######
### WC: Would it make more sense to define these structures in the ODE_solvers? This way you don't have to edit multiple files to add a new solver.
######

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


######
### WC: Why are these functions in this file? Shouldn't these be all located in General.jl?
######
#Internal for now
#Check Zeno function for now:
function check_zeno(jump_interval, last_jump_interval, zeno_count, t_star, tol, zeno_ratio, max_zeno_jumps)
    status =:continue

    is_contracting = jump_interval < last_jump_interval * zeno_ratio
    is_already_at_accumulation = (jump_interval <= tol && last_jump_interval <= tol)

    if  is_contracting || is_already_at_accumulation
        zeno_count += 1
        if zeno_count >= max_zeno_jumps
            @warn "Zeno Behavior Detected: System reached max zeno jumps ($max_zeno_jumps)"
            status =:terminate
        elseif jump_interval <= tol && zeno_count > 3
            @info "Zeno Accumulation Point Reached at t = $t_star. Terminating."
            status=:terminate
        end
    else
        zeno_count = 0
    end
    return zeno_count, status
end

function check_system_pathology(
    jump_interval, last_intervals, 
    state_step,
    zeno_count, instant_jump_count,
    t_star, tol, zeno_ratio, max_zeno_jumps,
    beating_warn_threshold, max_instant_jumps,
    min_zeno_confirmations = 3)

    status = :continue

    #Time min checks
    t_at_floor = jump_interval <= tol

    #time contraction check. Is it smaller than the max of recent history? 
    t_contracting = !isempty(last_intervals) && jump_interval < maximum(last_intervals) * zeno_ratio

    #State min check
    x_at_floor = state_step <= tol

    #Case 1: Zeno Termination
    #Time is effectivly flatlined and discrete jump is gone. If we have the history of tight jumps, we call it Zeno
    if t_at_floor && x_at_floor && zeno_count >= min_zeno_confirmations
        @info "Zeno Accumulation Point Reached at t = $t_star. Terminating."
        return zeno_count, instant_jump_count, :terminate
    end

    #Case 2:  Active Zeno tracking
    if t_contracting || (t_at_floor && !x_at_floor)
        if t_contracting
            zeno_count += 1
        end
        instant_jump_count = 0

        if zeno_count >= max_zeno_jumps
            @warn "Zeno Behaviour Detected: System reached max zeno jumps ($max_zeno_jumps)."
            return zeno_count, instant_jump_count, :terminate
        end
    else 
        if !t_at_floor
            zeno_count = 0 
        end
    end

    #Case 3: Beating and blocking detection
    if t_at_floor && !x_at_floor
        instant_jump_count += 1
        if instant_jump_count == beating_warn_threshold
            @info "Beating Detected: System has undergone $beating_warn_threshold instant jumps at t = $t_star."
        elseif instant_jump_count >= max_instant_jumps
            @warn "Blocking Detected: System Trapped on Guard ($max_instant_jumps instant jumps). Terminating."
            return zeno_count, instant_jump_count, :terminate
        end
    else
        instant_jump_count = 0
    end
    return zeno_count, instant_jump_count, status
end
