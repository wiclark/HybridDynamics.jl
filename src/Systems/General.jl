# A general hybrid dynamical system
struct GeneralSystem <: AbstractHybridSystem
    f::Function     #Continuous Dynamics: (x,t) -> dx/dt
    h::Function     #Guard Surface: x -> real
    Δ::Function     #Reset map: x -> x⁺
    direction::Int  #Direction for the guard
end

#Constructor for above
function GeneralSystem(f, h, Δ; direction::Int=0)
    return GeneralSystem(f, h, Δ, direction)
end
# GeneralSystem(f,h,Δ,direction::Int=0) = GeneralSystem(f,h,Δ,direction)

struct GeneralSol{T, X, DX} <: AbstractHybridSolution
    t::T                        #Time points of sim
    x::X                        #State traj: T is generic to support varying state types
    dx::DX 
    jump_times::Vector{Float64} #Explicit storage of timestamps where resets occurred
    jump_indices::Vector{Int}   #Map of jump_times to indices in the x and t vectors
end

#Internal - to initialize the solution struct
function GeneralSol(prob::prob{F, I, T}) where {F<:GeneralSystem, I, T}
    return GeneralSol([prob.tspan[1]], 
        [prob.init], 
        Vector{Vector{Float64}}(),
        Float64[],
        Int[])
end
#=
function init_solution(prob::prob{F, I, T}) where {F<:GeneralSystem, I, T}
    return GeneralSolution([prob.tspan[1]], [prob.init], Vector{Vector{Float64}}(), Float64[], Int[])
end
=#
#Internal
function guard(sys::GeneralSystem, x::AbstractArray)
    x_phys = (x isa AbstractMatrix) ? x[:, 1] : x
    val = sys.h(x_phys)
    # Return the most negative value (the most violated constraint)
    # If all are positive, it returns the smallest positive value.
    return val isa AbstractVector ? minimum(val) : val
end
#Internal
function apply_reset(sys::GeneralSystem, x::AbstractArray)
    return sys.Δ(x)
end

"""
    solve(prob::prob{F, I, T}, solver::AbstractODESolver=RK45(); kwargs...) where {F<:GeneralSystem, I, T}

ARGUMENT KEY

## Required:
'prob': The problem definition containing the system dynamics 'sys', initial state, and time span.
'solver': The numerical integration method used for continous steps. Defaults to RK45().

## Optional:

### Simulation and Step Controls:
* 'dt_initial' (Float64, default '0.01'): The starting time step for the continuous solver.
* 'dt_min' (Float64, default '1e-6'): The absolute minimum allowable time step. If the solver tries to go below this, the simulation terminates.
* 'max_iter' (Int, default '10^6'): The maximum number of continuous integration steps allowed before forcing a timeout.
'tol' (Float64, default '1e-6'): The baseline numerical tolerance used across the solver. Acts as the foundational unit for multipliers below.

### Event Handling
* 'event_method' (AbstractEventLocator, default 'LinearLocator()'): The algorithm used to pinpoint the exact time and state of a guard crossing.
* 'stepper' (AbstractODESolver. default 'RK4()'): The secondary ODE solver used internally by the event detection locator to pinpoint the impact state.
* 'guard_direction' (Int, default 'default_guard_direction(prob.sys)'): 0 -> detects crossings in either direction. 1 -> detects increasing crossings. -1 -> detects decreating crossings.

### Pathology Tuning
* 'zeno_ratio' (Float64, default '.90'): The ratio threshold for Zeno detection. If consecutive jump intervals contract by this ratio (or faster) it triggers a Zeno classification. 
* 'max_zeno_jumps' (Int, default '3'): The maximum number of consecutuve Zeno contractions allowed before terminating the simulation. 
* 'max_instant_jumps' (Int, default '5'): The maximum number of instantaneous jumps allowed before classifying the system as blocked and terminating. 
* 'max_buffer_size' (Int, default '5'): The number of previous jump intervals stored in memory to evaluate Zeno contractions.

### Fine-Tuning Multipliers (scaled against 'tol')
* 'min_zeno_history' (Int, default '2'): The minimum number of recorded jumps required before the solver will attempt to calculate a Zeno contraction ratio. 
* 'zeno_floor_mult' (Float64, default '2.0'): Defines the numerical floor ('tol * zeno_floor_mult'). If a jump interval falls below this, it maintains a Zeno state to prevent
machine precision drops into beating blocks.
* 'zeno_time_threshold' (Float64, default '1e-2'): The absolute maximum duration a jump interval can be while still being eligible for Zeno contraction classification.
* 'zeno_reset_mult' (Float64, default '100.0'): If a jump interval exceeds 'tol * zeno_reset_mult', the system is deemed safe and the Zeno counter is decremented.
* 'beating_tol_mult' (Float64, default '1.0'): Defines the time window ('tol * beating_tol_mult'). Jumps occurring within this window are classifed as instantaneous jumps or 'beating'
* 'adaptive_tol_mult' (Float64, defualt '100.0'): The boundary distance multiplier ('tol * adaptive_tol_mult'). When the system enters this distance from the guard, it shrinks step size to increase resolution.
* 'sliding_tol_mult' (Float64, default '10.0'): If the post-impact guard value is within 'tol * sliding_tol_mult', the solver enters sliding mode to suppress immediate erroneous events. This helps us avoid strange chattering.  

"""
function solve(prob::prob{F, I, T}, solver::AbstractODESolver=RK4();
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial=0.01, dt_min = 1e-6, max_iter = 10^6,
               tol = 1e-6,
               zeno_ratio = 0.90, max_zeno_jumps = 10,
               stepper::AbstractODESolver=RK4(),
               max_buffer_size=10,
               max_instant_jumps = 5,
               guard_direction = prob.sys.direction,
               ) where {F<:GeneralSystem, I, T}

    sys = prob.sys
    f = sys.f
    sol = GeneralSol(prob)
    t_start, t_end = prob.tspan
    Δt = dt_initial
    iter = 0

    # Pathology trackers
    instant_jump_count = 0
    zeno_count = 0
    last_jump_time = t_start
    last_intervals = Float64[]


    while sol.t[end] < t_end
        # Halt if we hit the iteration limit 
        iter += 1
        if iter > max_iter
            @info "Maximum iterations $max_iter reached."
            break
        end

        # Terminate if time is below machine precision
        if t_end - sol.t[end] <= eps(t_end)
            @info "Time step below minimum threshold $dt_min. Terminating."
            break
        end

        # The current state/time
        xₖ = sol.x[end]
        tₖ = sol.t[end]

        # Truncate time step if we overshoot
        Δt = (tₖ + Δt > t_end) ? (t_end - tₖ) : Δt
        x_predict, eventtriggered, t_root, dt_used, dt_next = take_step(solver, prob, f, xₖ, tₖ, Δt, tol, sol; guard_direction=guard_direction)

        if eventtriggered
            # An event has been found, time to locate it
            t_star, x_star = locate_event(event_method, prob, solver, f, xₖ, tₖ, Δt, guard(sys, xₖ), tol, sol, stepper)
            # Record the impact information
            if abs(guard(sys, x_star)) > 1e-3
                @warn "Event location is not located on the guard."
            end
            x⁺ = apply_reset(sys, x_star)

            # Perform a naive pathological check
            if (guard_direction == 0) && (abs(guard(sys, x⁺)>tol))
                @warn "An instantaneous reset is detected. Possible beating/blocking/Zeno. Terminating."
                break
            end
            # Pre-reset
            push!(sol.jump_times, t_star)
            push!(sol.t, t_star)
            push!(sol.x, x_star)
            # Post-reset
            push!(sol.t, t_star)
            push!(sol.x, x⁺)
            # println(x⁺)
            if dense_out
                push!(sol.dx, f(x_star, t_star))
                push!(sol.dx, f(x⁺, t_star))
            end
            
            # Check for fast switchings
            if (length(sol.jump_times) > 2) && (sol.jump_times[end]-sol.jump_times[end-1]<1e-3)
                @warn "Fast switching is detected. Possible beating/blocking/Zeno. Terminating"
                break
            end

        else
            push!(sol.t, tₖ+dt_used)
            push!(sol.x, x_predict)
            if dense_out
                push!(sol.dx, f(x_predict, tₖ+dt_used))
            end
        end
        Δt = dt_next

        #=
        #Discrete event logic
        if eventtriggered
            
            t_star = t_root
            τ = t_star - tₖ

            x_star, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ, tol, sol; guard_direction=guard_direction)

            jump_interval = t_star - last_jump_time

            status = check_general_pathology(
                jump_interval,
                last_intervals, t_star;
                tol=tol, 
                zeno_ratio=zeno_ratio,
                max_zeno_jumps=max_zeno_jumps, 
                max_buffer_size=max_buffer_size)

            #Terminate sim if pathological behavior is confirmed
            if status == :terminate
                break
            end
            last_jump_time = t_star

            #apply state jump (reset map)
            x⁺ = apply_reset(sys, x_star)
            h⁺ = guard(sys, x⁺)

            #Record impact point and post reset state
            push!(sol.t, t_star, t_star)
            push!(sol.x, x_star, x⁺)

            if hasproperty(sol, :jump_times)
                push!(sol.jump_times, t_star)
                push!(sol.jump_indices, length(sol.x))
            end
            #Reset step size to default after event
            Δt = dt_initial
        else
            #Normal cont updated
            if dense_out
                #Added phys_x for Var Eq
                phys_x = (xₖ isa AbstractMatrix) ? xₖ[:, 1] : xₖ
                push!(sol.dx, f(phys_x, tₖ))
            end
            push!(sol.t, tₖ + dt_used)
            push!(sol.x, x_predict)
            Δt = dt_next
        end =#
    end
    return sol
end