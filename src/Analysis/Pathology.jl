#NEW ZENO PLEASE WORK!!!!!!!!!!!!!!!!!!!
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
#External

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

######
### WC: This function is far too complicated.
### (Additionally, it does not work. It incorrectly identifies the NL bouncing ball as Zeno quite prematurely.)
######
function check_general_pathology(
    jump_interval,
    last_intervals,
    t_star;
    tol = 1e-6,
    zeno_ratio = 0.90,
    max_zeno_jumps = 8, 
    max_buffer_size = 10
    )
    #Adds current time between events to buffer
    push!(last_intervals, jump_interval)

    #keeps buffer at fixed size to ensure we only analyze revent stuff to avoid memeory overload. 
    if length(last_intervals) > max_buffer_size
        popfirst!(last_intervals)
    end

    #This rarely triggers but I think its good to have just in case of numerical evil
    if jump_interval < tol
        @warn "Micro jumps detected at t = $t_star. Terminating."
        return :terminate
    end

    len = length(last_intervals)
    if len >= max_zeno_jumps
        #We only perform analysis once we have enough data
        recent = last_intervals[end-max_zeno_jumps+1:end]
        n = length(recent)
        
        #Calc mean duration of recent jumps
        μ = sum(recent) / n

        #Part of the vairiance calc
        sum_sq_diff = 0.0
        for val in recent
            sum_sq_diff += (val - μ)^2
        end
        #standard deviation, representing how much jump intervals flucutate. The logic here seems to work but definitely double check my work - DS
        σ = sqrt(sum_sq_diff / n)

        #This detects beating or blocking. mu < 1e-3 for small jumps (I will add tols here later)
        #σ/μ coeff of variation less than 5%.
        #This means intervals are not shrinking, they get stuck at a fixed tiny values. hence beating/blocking. 
        if μ < 1e-3 && (μ > 0 ? σ/μ < 0.05 : σ < 0.05)
            @warn "Exceptional Pathology Detected at t = $t_star. Beating and Blocking suspected. Terminating."
            return :terminate
        end

        #Calc ratios of consecutive jump intervals. 
        #In Zeno this is < 1 (so its contracting)
        sum_ratios = 0.0
        for i in 2:n
            sum_ratios += recent[i] / recent[i-1]
        end
        #Avg rate of decary of intervals
        avg_ratio = sum_ratios / (n - 1)
        
        #IF avg decay is faster than the threshold we confirm Zeno. 
        if avg_ratio < zeno_ratio
            @warn "Zeno Detected at t = $t_star. Terminating."
            return :terminate
        end
    end
    
    return :continue
end