
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
    val = sys.h(x)
    return val isa AbstractVector ? minimum(abs.(val)) : val
end
#Internal
function apply_reset(sys::GeneralSystem, x::AbstractVector)
    return sys.Δ(x)
end

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

function variational_vector_field(f, u, t, n)
    # Unpack state
    x = u[1:n]
    Φ = reshape(u[n+1:end], n, n)

    # Base dynamics: x' = f(x, p, t)
    dx = f(x, t)

    # Variational dynamics: Φ' = A(t)Φ
    A = ForwardDiff.jacobian(y -> f(y, t), x)
    dΦ = A * Φ

    # Return augmented derivative
    return vcat(dx, vec(dΦ))
end

function compute_pushforward(f, Δ, h_guard, x⁻, t)
    n = length(x⁻)
    Id = I(n)

    # Eval field at boundaries (using p for parameters)
    f⁻ = f(x⁻, t)
    x⁺ = Δ(x⁻, t)
    f⁺ = f(x⁺, t)

    # Compute grads and jacob via ForwardDiff
    dh⁻ = ForwardDiff.gradient(h_guard, x⁻)
    DΔ⁻ = ForwardDiff.jacobian(Δ, x⁻)

    # Check dh(x) * f(x) = 0
    denom = dot(dh⁻, f⁻)
    if abs(denom) < 1e-12
        @warn "Non-transversal crossing detected: Trajectory is tangent to guard surface."
    end

    # Outer prods
    term1 = Id - (f⁻ * dh⁻') ./ denom
    term2 = (f⁺ * dh⁻') ./ denom

    # Full pushforward Δᶠ_*
    Δ_star_f = DΔ⁻ * term1 + term2

    return Δ_star_f
end

function apply_variational_jump(u, n, f, Δ, h_guard, t)
    x⁻ = u[1:n]
    Φ⁻ = reshape(u[n+1:end], n, n)

    # Compute the pushforward before state updates
    Δ_star_f = compute_pushforward(f, Δ, h_guard, x⁻, t)

    # Apply disc jump to base state x⁺ = Δ(x⁻)
    x⁺ = Δ(x⁻, t)

    # Apply pf mapping to fund matrix: Φ⁺ = Δ_*^f * Φ⁻
    Φ⁺ = Δ_star_f * Φ⁻

    # Update state vector in-place
    u[1:n] .= x⁺
    u[n+1:end] .= vec(Φ⁺)
    
    return u
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
               max_instant_jumps = 5) where {F<:GeneralSystem, I, T}
    
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
        x_predict, eventtriggered, h_now, dt_used, dt_next = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol)

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

            # Apply reset early to calculate the spatial step size
            x⁺ = sys.Δ(x_star)
           
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