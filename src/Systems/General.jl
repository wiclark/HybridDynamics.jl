
struct GeneralSystem <: AbstractHybridSystem
    f::Function     #Continuous Dynamics: (x,t) -> dx/dt
    h::Function     #Guard Surface: x-> real
    Δ::Function     #Reset map: x-> x⁺
end

struct GeneralSolution{T} <: AbstractHybridSolution
    t::Vector{Float64}
    x::Vector{T}
    jump_times::Vector{Float64}
    jump_indices::Vector{Int}
end

#Internal
function init_solution(prob::prob{F, I, T}) where {F<:GeneralSystem, I, T}
    return GeneralSolution{I}([prob.tspan[1]], [prob.init], Float64[], Int[])
end
#Internal
function guard(sys::GeneralSystem, x::AbstractArray)
    # If it's a matrix (like augmented state U), extract the physical state
    x_phys = (x isa AbstractMatrix) ? x[:, 1] : x
    val = sys.h(x_phys)
    return val isa AbstractVector ? minimum(abs.(val)) : val
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
    max_buffer_size)
    #There is room for a lot of tolerance stuff here. I will work on that at some point - DS

    # 1. Update history FIRST so Zeno can be evaluated
    push!(last_intervals, jump_interval)
    if length(last_intervals) > max_buffer_size
        popfirst!(last_intervals)
    end

    # 2. Zeno Check before anything else
    is_contracting = length(last_intervals) >= 3 && 
                     (last_intervals[end] < last_intervals[end-1] * zeno_ratio) &&
                     (last_intervals[end-1] < last_intervals[end-2] * zeno_ratio)

    # If we are already in a Zeno and hit the numerical floor, 
    # maintain the Zeno classification instead of dropping to Blocking.THIS HAPPENED SO MANY TIMES
    hit_zeno_floor = (zeno_count > 0) && (jump_interval <= tol)

    if (is_contracting && jump_interval < 1e-2) || hit_zeno_floor
        zeno_count += 1
        instant_jump_count = 0 # Explicitly bypass and reset the blocking trap
        @info "Zeno contraction detected. count: $zeno_count"
        
        if zeno_count >= max_zeno_jumps
            @warn "Zeno Accumulation Point Reached at t = $t_star. Terminating."
            return zeno_count, instant_jump_count, :terminate
        end
        
        return zeno_count, instant_jump_count, :continue
    else
        # Only reset if the interval genuinely grows or stabilizes outside Zeno
        zeno_count = 0 
    end

    # 3. Beating and Blocking Check (Only evaluated if NOT Zeno)
    if jump_interval <= tol
        instant_jump_count += 1
        
        if instant_jump_count >= max_instant_jumps
            @warn "Blocking Detected at t = $t_star (Exceeded max instant jumps). Terminating."
            return zeno_count, instant_jump_count, :terminate
        end
        
        @info "Beating event $instant_jump_count at t = $t_star"
        return zeno_count, instant_jump_count, :continue
    end

    # 4. Continuous movement
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
    t_list = sol.t
    states = sol.x
    
    n_total = length(states[1])

    data_list = Vector{Float64}[]
    
    for i in 1:length(states)
        push!(data_list, vec(states[i])) 
        if i < length(sol.t) && sol.t[i] == sol.t[i+1]
            push!(data_list, fill(NaN, n_total))
        end
    end

    data = reduce(hcat, data_list)'
    
    return t_list, data
end

#External
function solve(prob::prob{F, I, T}, solver::AbstractODESolver=RK45(); 
               event_method::AbstractEventLocator=LinearLocator(), 
               dt_initial=.01, dt_min = 1e-6, max_iter = 10^6, 
               tol = 1e-6, 
               zeno_ratio = 0.90, max_zeno_jumps = 15,
               stepper::AbstractODESolver=ModifiedTrap(),
               max_buffer_size=5,
               beating_warn_threshold=3,
               max_instant_jumps = 5,
               guard_direction::Symbol = default_guard_direction(prob.sys)) where {F<:GeneralSystem, I, T}
    
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

        # Continuous integration
        xₖ = sol.x[end]
        tₖ = sol.t[end]

        # Attempt continuous step
        x_predict, eventtriggered, h_now, dt_used, dt_next = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol; guard_direction=guard_direction)

        # Discrete event logic
        if eventtriggered
            # Pinpoint exact time and state event happened
            t_star, x_star = locate_event(event_method, prob, solver, f, xₖ, tₖ, dt_used, h_now, tol, sol, stepper)

            #Pathology Checks
            jump_interval = t_star - last_jump_time
            zeno_count, instant_jump_count, status = check_system_pathology(
                jump_interval, last_intervals, 
                zeno_count, instant_jump_count,
                t_star, tol, zeno_ratio, max_zeno_jumps, max_instant_jumps,
                max_buffer_size
            )

            if status == :terminate
                break
            end
            last_jump_time = t_star

            x⁺ = apply_reset(sys, x_star)
           
            push!(sol.t, t_star, t_star)
            push!(sol.x, x_star, x⁺)

            if hasproperty(sol, :jump_times)
                push!(sol.jump_times, t_star)
                push!(sol.jump_indices, length(sol.x))
            end

            Δt = dt_initial

        else 
            push!(sol.t, tₖ + dt_used)
            push!(sol.x, x_predict)
            Δt = dt_next
        end
    end
    return sol
end