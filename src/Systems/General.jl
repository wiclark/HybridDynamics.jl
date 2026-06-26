
struct GeneralSystem <: AbstractHybridSystem
    f::Function     #Continuous Dynamics: (x,t) -> dx/dt
    h::Function     #Guard Surface: x-> real
    Δ::Function     #Reset map: x-> x⁺
    direction::Int  #Direction for the guard
end

#Constructor for above
GeneralSystem(f,h,Δ,direction::Int=0) = GeneralSystem(f,h,Δ,direction)

struct GeneralSolution{X, DX} <: AbstractHybridSolution
    t::Vector{Float64}          #Time points of sim
    x::Vector{X}                #State traj: T is generic to support varying state types
    dx::DX 
    jump_times::Vector{Float64} #Explicit storage of timestamps where resets occurred
    jump_indices::Vector{Int}   #Map of jump_times to indices in the x and t vectors
end

#Internal
function init_solution(prob::prob{F, I, T}) where {F<:GeneralSystem, I, T}
    return GeneralSolution([prob.tspan[1]], [prob.init], Vector{Vector{Float64}}(), Float64[], Int[])
end
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
function solve(prob::prob{F, I, T}, solver::AbstractODESolver=RK45();
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial=0.01, dt_min = 1e-6, max_iter = 10^6,
               tol = 1e-6,
               zeno_ratio = 0.90, max_zeno_jumps = 5,
               stepper::AbstractODESolver=RK4(),
               max_buffer_size=10,
               max_instant_jumps = 5,
               guard_direction = prob.sys.direction,
               ) where {F<:GeneralSystem, I, T}

    sys = prob.sys
    f = sys.f
    sol = init_solution(prob)
    t_start, t_end = prob.tspan
    Δt = dt_initial
    iter = 0
    sol.x = []

    # Pathology trackers
    instant_jump_count = 0
    zeno_count = 0
    last_jump_time = t_start
    last_intervals = Float64[]

    #State flag: if true, the system is constrained on the guard. This is to avoid weird sliding and it does some work for falling through the guard too!
    #I will add, for beating and blocking we will most likely need to fine tune the tolerances. I think this will mainly be done when we really test things
    #Perhaps I will also find a smarter way to do this but for now we are good. Liam can mess with it
    in_sliding_mode = false

    while sol.t[end] < t_end
        iter += 1
        if iter > max_iter
            @info "Maximum iterations $max_iter reached."
            break
        end

        # Terminate if time is below machine precision
        if t_end - sol.t[end] < dt_min
            @info "Time step below minimum threshold $dt_min. Terminating."
            break
        end

        # Truncate time step if we overshoot
        dt_step = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

        dt_step = min(Δt, t_end - sol.t[end])

        # Continuous integration
        xₖ = sol.x[end]
        tₖ = sol.t[end]
        hₖ = guard(sys, xₖ)

        x_predict, eventtriggered, t_root, dt_used, dt_next = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol; guard_direction=guard_direction)
   
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

            in_sliding_mode = false

            #Record impact point and post reset state
            push!(sol.t, t_star, t_star)
            push!(sol.x, x_star, x⁺)

            if hasproperty(sol, :jump_times)
                push!(sol.jump_times, t_star)
                push!(sol.jump_indices, length(sol.x))
            end
            #Reset step size to default after event
            Δt = dt_initial
            in_sliding_mode = false
        else
            #Normal cont updated
            if dense_out
                phys_x = (xₖ isa AbstractMatrix) ? xₖ[:, 1] : xₖ
                push!(sol.dx, f(phys_x, tₖ))
            end
            push!(sol.t, tₖ + dt_used)
            push!(sol.x, x_predict)
            Δt = dt_next
            in_sliding_mode = false
        end
    end
    return sol
end