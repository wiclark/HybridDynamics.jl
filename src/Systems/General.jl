
struct GeneralSystem <: AbstractHybridSystem
    f::Function     #Continuous Dynamics: (x,t) -> dx/dt
    h::Function     #Guard Surface: x-> real
    Δ::Function     #Reset map: x-> x⁺
end

struct GeneralSolution <: AbstractHybridSolution
    t::Vector{Float64}
    x::Vector{Vector{Float64}}
    jump_times::Vector{Float64}
    jump_indices::Vector{Int}
end

#Internal
function init_solution(prob::prob{GeneralSystem, I, T}) where {I, T}
    return GeneralSolution([prob.tspan[1]], [prob.init], Float64[], Int[])
end
#Internal
function guard(sys::GeneralSystem, x::AbstractVector)
    return sys.h(x)
end
#Internal
function apply_reset(sys::GeneralSystem, x::AbstractVector)
    return sys.Δ(x)
end

#Internal
######
### WC: Should there be somesort of for loop in this to count up?
######
function check_beating_blocking(jump_interval, instant_jump_count, t_star, tol, beating_warn_threshold, max_instant_jumps)
    status = :continue
    if jump_interval <= tol
        instant_jump_count += 1
        if instant_jump_count == beating_warn_threshold
            @info "Beating Detected: System has undergone $beating_warn_threshold instant jumps at t = $t_star"
        elseif instant_jump_count >= max_instant_jumps
            @warn "Blocking Detected: System trapped on guard ($max_instant_jumps instant jumps). Terminating."
            status =:terminate
        end
    else
        instant_jump_count = 0
    end
    return instant_jump_count, status
end

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
    t_contracting = !isempty(last_intervals) && (jump_interval > tol) && (jump_interval < maximum(last_intervals) * zeno_ratio)

    #State min check
    x_at_floor = state_step <= tol

    #Case 1: Pure Zeno Termination
    #Time is effectivly flatlined and discrete jump is gone. If we have the history of tight jumps, we call it Zeno
    if t_at_floor && zeno_count >= min_zeno_confirmations
        @info "Zeno Accumulation Point Reached at t = $t_star. Terminating."
        return zeno_count, instant_jump_count, :terminate
    end

    #Case 2: Active Zeno tracking
    if t_contracting
        zeno_count += 1
        instant_jump_count = 0
        if zeno_count >= max_zeno_jumps
            @warn "Zeno Behavior Detected: System reach max zeno jumps ($max_zeno_jumps)."
            return zeno_count, instant_jump_count, :terminate
        end
    else 
        #if we take normal step that isnt shrinking, reset zeno 
        if !t_at_floor
            zeno_count = 0
        end
    end

    #Case 3: Beating and blocking detection
    #if we are trapped in time, but it is NOT end of Zeno sequence
    if t_at_floor && zeno_count < min_zeno_confirmations
        instant_jump_count += 1
        if instant_jump_count == beating_warn_threshold
            @info "Beating Detected: System has undergone $beating_warn_threshold instant jumps at t = $t_star."
        elseif instant_jump_count >= max_instant_jumps
            @warn "Blocking Detected: System trapped on guard ($max_instant_jumps instant jumps). Terminating."
            return zeno_count, instant_jump_count, :terminate
        end
    else 
        instant_jump_count = 0
    end
    return zeno_count, instant_jump_count, status
end

#External
function solve(prob::prob{F, I, T}, solver::AbstractODESolver=RK45(); 
               event_method::AbstractEventLocator=LinearLocator(), dt_initial=.01, dt_min = 1e-6, 
               max_iter = 10^6, tol = 1e-6, beating_warn_threshold = 3, max_instant_jumps = 100, 
               zeno_ratio = .99, max_zeno_jumps = 100, history_window = 3, min_zeno_confirmations = 3) where {F<:GeneralSystem, I, T}
    
    sys = prob.sys
    f = sys.f
    sol = init_solution(prob)
    t_start, t_end = prob.tspan
    Δt = dt_initial 
    iter = 0

    # Pathology trackers
    last_jump_time = -Inf
    last_intervals = Float64[] # Replaces last_jump_interval with a rolling history window
    instant_jump_count = 0
    zeno_count = 0

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
        x_predict, eventtriggered, h_now, dt_next = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol)

        # Discrete event logic
        if eventtriggered
            # Pinpoint exact time and state event happened
            t_star, x_star = locate_event(event_method, prob, solver, f, xₖ, tₖ, dt_step, h_now, tol, sol)

            # Apply reset early to calculate the spatial step size
            x⁺ = sys.Δ(x_star)
            state_step = norm(x⁺ - x_star)

            # Calculate temporal step size
            jump_interval = t_star - last_jump_time

            # Unified Pathology Engine
            if last_jump_time == -Inf
                zeno_count = 0
                instant_jump_count = 0
            else 
                zeno_count, instant_jump_count, pathology_status = check_system_pathology(
                    jump_interval, last_intervals, state_step,
                    zeno_count, instant_jump_count, t_star, tol,
                    zeno_ratio, max_zeno_jumps, beating_warn_threshold,
                    max_instant_jumps, min_zeno_confirmations
                )
                
                if pathology_status == :terminate
                    break
                end

                # Update tracking histories
                #Put here to prevent pushin Inf into window on first jump.
                push!(last_intervals, jump_interval)
                if length(last_intervals) > history_window
                    popfirst!(last_intervals)
                end
            end
            last_jump_time = t_star

            # Track data: pre and post jump states
            push!(sol.t, t_star, t_star)
            push!(sol.x, x_star, x⁺)

            # Record for event analysis
            if hasproperty(sol, :jump_times)
                push!(sol.jump_times, t_star)
                push!(sol.jump_indices, length(sol.t))
            end

            # Shrink main step size to avoid overshooting
            Δt = dt_min

        else
            t_next = tₖ + dt_step
            push!(sol.t, t_next)
            push!(sol.x, x_predict)
            Δt = dt_initial
        end
    end
    return sol
end