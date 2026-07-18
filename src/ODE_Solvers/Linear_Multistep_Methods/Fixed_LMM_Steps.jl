abstract type FixedLMM <: LMM end
#comment
"""
Fixed second order Adams Bashforth multistep method. RK4 is used to caluclate starting values.
"""
struct AdamsBashforth2 <: FixedLMM end
"""
Fixed third order Adams Bashforth multistep method. RK4 is used to caluclate starting values.
"""
struct AdamsBashforth3 <: FixedLMM end
"""
Fixed second order backwards differentiation multistep method. RK4 is used to caluclate starting values.
"""
struct BDF2 <: FixedLMM end

lmm_order(::AdamsBashforth2) = 2
lmm_order(::AdamsBashforth3) = 3
lmm_order(::BDF2) = 2

function take_step(solver::FixedLMM, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver = RK4(); check=true, guard_direction=default_guard_direction(prob.sys))

    if !check
        return take_step(stepper, prob, f, xₖ, tₖ, Δt, tol, sol, stepper; check=false, guard_direction=guard_direction)
    end

    sys = prob.sys
    k = lmm_order(solver)

    # Determine how many continuous steps we have since the last jump
    history_len = isempty(sol.event_indices) ? length(sol.x) : (length(sol.x) - sol.event_indices[end])

    if history_len <= k
        return take_step(stepper, prob, f, xₖ, tₖ, Δt, tol, sol, stepper; check=check, guard_direction=guard_direction)
    end

    # Multistep phase: We do have history so we extract prev states. 
    x_history = sol.x[end - k + 1 : end - 1]
    t_history = sol.t[end - k + 1 : end - 1]

    time_diffs = diff(vcat(t_history, tₖ))

    if any(time_diffs .<= 1e-12)
        return take_step(stepper, prob, f, xₖ, tₖ, Δt, tol, sol, stepper; check=check, guard_direction=guard_direction)
    end

    # Pass history arrays forward
    x_predict = compute_lmm_step(solver, f, xₖ, tₖ, Δt, x_history, t_history)
    
    if check 
        h_now = guard(sys, xₖ)
        h_next = guard(sys, x_predict)

        if history_len > 1
            idx = length(sol.x) - 1
            t_prev = sol.t[idx]
            h_prev = guard(sys, sol.x[idx])
        else
            t_prev = tₖ - Δt
            h_prev = h_now
        end

        eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt; tol = tol, direction=guard_direction)
    
        return x_predict, eventtrigger, t_root, Δt, Δt
    else
        # Fallback (Structurally unreachable due to top delegation)
        return x_predict, false, NaN, Δt, Δt
    end
end

function compute_lmm_step(::AdamsBashforth2, f, xₖ, tₖ, Δt, x_history, t_history)
    x_prev = x_history[end]
    t_prev = t_history[end]

    t0 = Δt
    t1 = tₖ - t_prev

    fₖ = f(xₖ, tₖ)
    f_minus1 = f(x_prev, t_prev)

    β0 = t0*(2t1+t0)/(2t1)
    β1 = -t0^2/(2t1)

    return xₖ .+ β0 .* fₖ .+ β1 .* f_minus1
end

function compute_lmm_step(::AdamsBashforth3, f, xₖ, tₖ, Δt, x_history, t_history)
    # x_history[end] is x_{k-1}, x_history[end-1] is x_{k-2}
    x_prev1 = x_history[end]
    t_prev1 = t_history[end]
    
    x_prev2 = x_history[end-1]
    t_prev2 = t_history[end-1]
    
    fₖ = f(xₖ, tₖ)
    f_prev1 = f(x_prev1, t_prev1)
    f_prev2 = f(x_prev2, t_prev2)
    
    return xₖ .+ Δt .* ( (23/12) .* fₖ .- (16/12) .* f_prev1 .+ (5/12) .* f_prev2)
end

function compute_lmm_step(::BDF2, f, xₖ, tₖ, Δt, x_history, t_history)
    #retrieve state at previous time step.
    x_prev = x_history[end]
    t_new = tₖ + Δt

    #doing some algebra to isolate the implicit part of BDF2 and "c" is the right side
    c = (4.0 / 3.0) .* xₖ .- (1.0 / 3.0) .* x_prev
    #coeff scaling step size 
    α = 2.0 / 3.0

    #Gen initial guess using Exp Euler to help the Newton Root finder
    x_guess = xₖ .+ Δt .* f(xₖ, tₖ)
    #solve implicit system G(z) = 0 using Newtons
    return implicit_newton_solve(f, x_guess, c, α, Δt, t_new)
end