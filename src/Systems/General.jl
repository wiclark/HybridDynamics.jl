
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

function init_solution(prob::prob{GeneralSystem, I, T}) where {I, T}
    return GeneralSolution([prob.tspan[1]], [prob.init], Float64[], Int[])
end

function guard(sys::GeneralSystem, x::AbstractVector)
    return sys.h(x)
end
function apply_reset(sys::GeneralSystem, x::AbstractVector)
    return sys.Δ(x)
end

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

function solve(prob::prob{F, I, T}, solver::AbstractODESolver=ModifiedMidpoint(); event_method::AbstractEventLocator=QuadraticLocator(), dt_initial=.01, dt_min = 1e-6, max_iter = 10^6, tol = 1e-6, beating_warn_threshold = 3, max_instant_jumps = 100, zeno_ratio = .99, max_zeno_jumps = 100) where {F<:GeneralSystem, I, T}
    sys = prob.sys
    f = sys.f
    sol = init_solution(prob)
    t_start, t_end = prob.tspan
    Δt = dt_initial 
    iter = 0

    #pathology trackers
    last_jump_time = -Inf
    last_jump_interval = Inf
    instant_jump_count = 0
    zeno_count = 0

    while sol.t[end] < t_end 
        iter += 1
        if iter > max_iter 
            @info "Maximum iterations $max_iter reached."
            break
        end

        #terminaite if time is below machine precision
        if t_end - sol.t[end] < dt_min
            @info "Time step below minimum threshold $dt_min. Terminating."
            break
        end

        #truncate time step if we overshoot
        dt_step = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

        #continuous integration
        xₖ = sol.x[end]
        tₖ = sol.t[end]

        #attempt cont step
        x_predict, eventtriggered, h_now, dt_next = take_step(solver, sys, f, xₖ, tₖ, dt_step, tol, sol)

        #discrete event logic
        if eventtriggered
            #pinpoint exact time and state even happened
            t_star, x_star = locate_event(event_method, sys, solver, f, xₖ, tₖ, dt_step, h_now, tol, sol)

            #Pathology Checks
            jump_interval = t_star - last_jump_time

            if last_jump_time == Inf
                zeno_count = 0
            else 
                zeno_count, zeno_status = check_zeno(jump_interval, last_jump_interval, zeno_count, t_star, tol, zeno_ratio, max_zeno_jumps)
                if zeno_status == :terminate
                    break
                end
            end

            if zeno_count == 0 && last_jump_time != -Inf
                instant_jump_count, beat_status = check_beating_blocking(jump_interval, instant_jump_count, t_star, tol, beating_warn_threshold, max_instant_jumps)
                if beat_status == :terminate
                    break
                end
            else 
                instant_jump_count = 0
            end

            last_jump_interval = jump_interval
            last_jump_time = t_star

            #apply reset
            x⁺ = sys.Δ(x_star)

            #track data: pre and post jump states
            push!(sol.t, t_star, t_star)
            push!(sol.x, x_star, x⁺)

            #record for event analysis
            if hasproperty(sol, :jump_times)
                push!(sol.jump_times, t_star)
                push!(sol.jump_indices, length(sol.t))
            end

            #shrink main step size to avoid overshooting
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