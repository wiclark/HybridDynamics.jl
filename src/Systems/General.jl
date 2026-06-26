
struct GeneralSystem <: AbstractHybridSystem
    f::Function     #Continuous Dynamics: (x,t) -> dx/dt
    h::Function     #Guard Surface: x-> real
    Δ::Function     #Reset map: x-> x⁺
end

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
    check_system_pathology(jump_interval, last_intervals, zeno_count, instant_jump_count, 
                           t_star, tol, zeno_ratio, max_zeno_jumps, max_instant_jumps, max_buffer_size)

Evaluates the discrete event history of a hybrid system to detect and prevent pathological behaviors 
such as Zeno accumulation (infinite jumps in finite time), Beating (finite jumps instantaneously) and Blocking (infinite instantaneous jumps).

How it works:
1) History update: Adds the current jump interval to a rolling buffer ('last_intervals').
2) Zeno Detection: Checks if the last three intervals are contracting by a factor of 'zeno_ratio'.
   If a contraction is detected it increments 'zeno_count'. Once 'zeno_count' reaches 'max_zeno_jumps',
   we confirm Zeno behavior and terminate the simulation (Accumulation point logic will be added at a later date.)
3) Blocking Detection: If the system is not Zeno but the jump interval is effectively zero ( ≤ 'tol'), it is considered an instant jump. 
   if 'instant_jump_count' exceeds 'max_instant_jumps' the system is considered trapped and the simulation terminates.
4) Beating: Similar to Blocking check we check instant jumps but allow the system to potentially escape. 
5) Reset: IF continuous interval is large and no contraction is detected we proceed as normal and counters are reset. 

# Arguments
- `jump_interval::Float64`: The time elapsed since the last discrete jump event.
- `last_intervals::Vector{Float64}`: A rolling buffer storing the history of recent jump intervals.
- `zeno_count::Int`: The current number of consecutive Zeno contractions detected.
- `instant_jump_count::Int`: The current number of consecutive instant jumps detected.
- `t_star::Float64`: The timestamp of the current event (used for logging and warnings).
- `tol::Float64`: The numerical tolerance used to define an "instantaneous" jump.
- `zeno_ratio::Float64`: The multiplier used to define a contraction (e.g., `0.90`).
- `max_zeno_jumps::Int`: The number of consecutive contractions required to confirm a Zeno state and terminate.
- `max_instant_jumps::Int`: The maximum allowed instantaneous jumps before the system is declared blocked.
- `max_buffer_size::Int`: The maximum number of historical jump intervals to store in `last_intervals`.

# Returns
- `zeno_count::Int`: The updated counter for consecutive Zeno events.
- `instant_jump_count::Int`: The updated counter for consecutive instantaneous jumps.
- `status::Symbol`: Returns `:continue` to proceed, or `:terminate` if a pathology limit is reached.

"""

#NEW ZENO PLEASE WORK!!!!!!!!!!!!!!!!!!!
function check_system_pathology(
    jump_interval, last_intervals, 
    zeno_count, instant_jump_count,
    t_star, tol, zeno_ratio, max_zeno_jumps, max_instant_jumps,
    max_buffer_size;
    #Tunable zeno values - mainly for Linear/Affine
    zeno_floor_mult = 2.0,
    zeno_time_threshold = 1e-2,
    zeno_reset_mult = 100.0,
    beating_tol_mult = 1.0,
    min_zeno_history = 2)

    #Update history FIRST so Zeno can be evaluated
    push!(last_intervals, jump_interval)
    if length(last_intervals) > max_buffer_size
        popfirst!(last_intervals)
    end

    #Zeno Check before anything else
    is_contracting = length(last_intervals) >= min_zeno_history && 
                     (last_intervals[end] <= last_intervals[end-1] * zeno_ratio)

    # If we are already in a Zeno and hit the numerical floor, 
    # maintain the Zeno classification instead of dropping to Blocking.THIS HAPPENED SO MANY TIMES
    hit_zeno_floor = (jump_interval <= tol * zeno_floor_mult)

    in_zeno_state = (is_contracting && jump_interval < zeno_time_threshold) || hit_zeno_floor

    if in_zeno_state
        zeno_count += 1
        instant_jump_count = 0 # Explicitly bypass and reset the trap if it happens
        @info "Zeno contraction detected. count: $zeno_count"
        
        if zeno_count >= max_zeno_jumps
            @warn "Max Zeno jumps reached at t = $t_star. Terminating."
            return zeno_count, instant_jump_count, :terminate
        end
        
        return zeno_count, instant_jump_count, :continue
    elseif jump_interval > tol * zeno_reset_mult
        # Only reset if the interval genuinely grows or stabilizes outside Zeno
        zeno_count = max(0, zeno_count - 1)
    end

    #Beating and Blocking Check (Only evaluated if NOT Zeno)
    if jump_interval <= tol * beating_tol_mult
        instant_jump_count += 1
        
        if instant_jump_count >= max_instant_jumps
            @warn "Blocking Detected at t = $t_star (Exceeded max instant jumps). Terminating."
            return zeno_count, instant_jump_count, :terminate
        end
        
        @info "Beating event $instant_jump_count at t = $t_star"
        return zeno_count, instant_jump_count, :continue
    end

    #Continuous movement
    instant_jump_count = 0
    return zeno_count, instant_jump_count, :continue
end


#Function calcs the continuous time derivative for the augmented state. Runs the usual dx and dΦ.
function variational_vector_field(f, U::AbstractMatrix, t)
    # Unpack state
    x = U[:, 1]
    Φ = U[:, 2:end]

    # Base dynamics: x' = f(x, p, t)
    dx = f(x, t)

    # Variational dynamics: Φ' = A(t)Φ
    A = ForwardDiff.jacobian(y -> f(y, t), x)
    dΦ = A * Φ

    # Return augmented derivative
    return hcat(dx, dΦ)
end

#Calcs the Δ_*^f at a boundary. We look how how much the vf mismatches before and after the jump (f⁺/f⁻), how much it scales things (DΔ⁻), and how the trajectory hits boundary (dh⁻).
#Then we create a matrix that takes all of this.
function compute_pushforward(f, Δ, h_guard, x⁻, t)
    n = length(x⁻)
    Id = I(n)

    # Eval field at boundaries (using p for parameters)
    f⁻ = f(x⁻, t)
    x⁺ = Δ(x⁻, t)
    f⁺ = f(x⁺, t)

    # Compute grads and jacob via ForwardDiff
    dh⁻ = ForwardDiff.gradient(h_guard, x⁻)
    DΔ⁻ = ForwardDiff.jacobian(y -> Δ(y, t), x⁻)

    # Check dh(x) * f(x) = 0
    denom = dot(dh⁻, f⁻)
    if abs(denom) < 1e-6
        @warn "Non-transversal crossing detected: Trajectory is tangent to guard surface."
    end

    # Outer prods
    term1 = Id - (f⁻ * dh⁻') ./ denom
    term2 = (f⁺ * dh⁻') ./ denom

    # Full pushforward Δᶠ_*
    Δ_star_f = DΔ⁻ * term1 + term2

    return Δ_star_f
end

#When our solver hits h(x)=0 this is called. We get the pushforward matrix and multiply by Φ⁻ to get new Φ⁺. 
function apply_variational_jump(U::AbstractMatrix, f, Δ, h_guard, t)
    x⁻ = U[:, 1]
    Φ⁻ = U[:, 2:end]

    # Compute the pushforward before state updates
    Δ_star_f = compute_pushforward(f, Δ, h_guard, x⁻, t)

    # Apply disc jump to base state x⁺ = Δ(x⁻)
    x⁺ = Δ(x⁻, t)

    # Apply pf mapping to fund matrix: Φ⁺ = Δ_*^f * Φ⁻
    Φ⁺ = Δ_star_f * Φ⁻

    # Update state vector in-place
    U[:, 1] .= x⁺
    U[:, 2:end] .= Φ⁺
    
    return U
end

#Plotting Lines Expelled!
function split_jumps(sol::AbstractHybridSolution) 
    states = sol.x
    
    # Preallocate an array of exactly whatever type the state is
    T = typeof(states[1])
    data_list = T[]
    t_list = Float64[]
    
    # Create a NaN container of the exact same size and shape as the state
    nan_state = fill!(similar(states[1]), NaN)
    
    for i in 1:length(states)
        push!(data_list, states[i])
        push!(t_list, sol.t[i])
        
        # Insert the NaN break if a jump occurred (zero-time transition)
        if i < length(sol.t) && sol.t[i] == sol.t[i+1]
            push!(data_list, nan_state)
            push!(t_list, NaN)
        end
    end
    
    return t_list, data_list
end

#External

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
* 'guard_direction' (Symbol, default 'default_guard_direction(prob.sys)'): Dictates which zero-crossing direction triggers an event (e.g., positive-to-negative).

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
               zeno_ratio = 0.90, max_zeno_jumps = 3,
               stepper::AbstractODESolver=RK4(),
               max_buffer_size=5,
               max_instant_jumps = 5,
               guard_direction::Symbol = default_guard_direction(prob.sys),
               #Tolerance tuning
               min_zeno_history = 2,
               zeno_floor_mult = 2.0,
               zeno_time_threshold = 1e-2,
               zeno_reset_mult = 100.0,
               beating_tol_mult = 1.0,
               adaptive_tol_mult = 100.0,
               adaptive_dt_mult = 10.0,
               sliding_tol_mult = 10.0) where {F<:GeneralSystem, I, T}
    
    sys = prob.sys
    f = sys.f
    sol = init_solution(prob)
    t_start, t_end = prob.tspan
    Δt = dt_initial 
    iter = 0

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

        #Adaptive step sizing: slow down when near guard to increase resolution
        if abs(guard(sys, sol.x[end])) < tol * adaptive_tol_mult
            dt_step = min(dt_step, dt_min * adaptive_dt_mult)
        end

        # Continuous integration
        xₖ = sol.x[end]
        tₖ = sol.t[end]

        #Exit sliding mod if the trajectory moves back to safe place
        if in_sliding_mode && guard(sys, xₖ) > 0
            in_sliding_mode = false
        end

        # Attempt continuous step
        x_predict, eventtriggered, h_now, dt_used, dt_next = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol; guard_direction=guard_direction)

        #Logic: if we are sliding, we ignore further events to avoid breaking everything
        if in_sliding_mode
            eventtriggered = false
        end

        #Check for crossing: either explict event flag or a sign flip in the guard function
        if !eventtriggered && !in_sliding_mode && (guard(sys, xₖ) * guard(sys, x_predict) < 0 || guard(sys, x_predict) <= 0)
            eventtriggered = true
        end

        #Discrete event logic
        if eventtriggered
            # Pinpoint exact time and state event happened
            t_star, x_star = locate_event(event_method, prob, solver, f, xₖ, tₖ, dt_used, h_now, tol, sol, stepper)

            #Pathology Checks
            jump_interval = t_star - last_jump_time
            zeno_count, instant_jump_count, status = check_system_pathology(
                jump_interval, last_intervals, 
                zeno_count, instant_jump_count,
                t_star, tol, zeno_ratio, max_zeno_jumps, max_instant_jumps,
                max_buffer_size;
                min_zeno_history=min_zeno_history,
                zeno_floor_mult=zeno_floor_mult,
                zeno_time_threshold=zeno_time_threshold,
                zeno_reset_mult=zeno_reset_mult,
                beating_tol_mult=beating_tol_mult
            )

            #Terminate sim if pathological behavior is confirmed
            if status == :terminate
                break
            end
            last_jump_time = t_star

            #Trigger sliding mode after impact
            in_sliding_mode = true
            if abs(h_now) < tol * sliding_tol_mult
                in_sliding_mode = true
            end

            #apply state jump (reset map)
            x⁺ = apply_reset(sys, x_star)
           
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
                push!(sol.dx, f(xₖ, tₖ)) # hey, it works (usually)
            end
            push!(sol.t, tₖ + dt_used)
            push!(sol.x, x_predict)
            Δt = dt_next
        end
    end
    return sol
end