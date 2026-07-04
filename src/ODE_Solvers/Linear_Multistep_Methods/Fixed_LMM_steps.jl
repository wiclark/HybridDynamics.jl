abstract type FixedLMM <: LMM end

struct AdamsBashforth2 <: FixedLMM end
struct AdamsBashforth3 <: FixedLMM end
struct BDF2 <: FixedLMM end

lmm_order(::AdamsBashforth2) = 2
lmm_order(::AdamsBashforth3) = 3
lmm_order(::BDF2) = 2

function take_step(solver::FixedLMM, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver = RK4();  guard_direction=default_guard_direction(prob.sys))
    sys = prob.sys
    k = lmm_order(solver)

    #Determine how many continuous steps we have since the last jump
    history_len = isempty(sol.event_indices) ? length(sol.x) : (length(sol.x) - sol.event_indices[end])

    h_now = guard(sys, xₖ)

    if history_len < k
        # Startup phase: Use single step predictor
        x_predict = compute_step(stepper, f, xₖ, Δt, tₖ)
        h_next = guard(sys, x_predict)

        if history_len > 1
            idx = length(sol.x) - 1
            t_prev = sol.t[idx]
            h_prev = guard(sys, sol.x[idx])
        else 
            t_prev = tₖ - Δt
            h_prev = h_now
        end

        eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt; tol=tol, direction=guard_direction)

        if eventtrigger
            if (t_root - tₖ) < (1e-4 * Δt) # Add a small buffer
                eventtrigger = false
                t_root = tₖ + Δt # Reset t_root to end of step
            end
        end
        return x_predict, eventtrigger, t_root, Δt, Δt
    else
        #Multistep phase: We do have history so we extract prev states. 
        x_history = sol.x[end - k + 1 : end - 1]
        t_history = sol.t[end - k + 1 : end - 1]

        #pass history arrays forward
        x_predict = compute_lmm_step(solver, f, xₖ, tₖ, Δt, x_history, t_history)
        
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
        if eventtrigger
            if (t_root - tₖ) < (1e-4 * Δt) # Add a small buffer
                eventtrigger = false
                t_root = tₖ + Δt # Reset t_root to end of step
            end
        end

        return x_predict, eventtrigger, t_root, Δt, Δt
    end
end

function compute_lmm_step(::AdamsBashforth2, f, xₖ, tₖ, Δt, x_history, t_history)
    #x_history[end] is x_{k-1}
    x_prev = x_history[end]
    t_prev = t_history[end]

    #Calc previous step size
    dt_previous = tₖ - t_prev

    fₖ = f(xₖ, tₖ)
    f_prev_val = f(x_prev, t_prev)

    #Variable step AB2 Formula
    α = Δt / dt_previous
    return xₖ .+ Δt .* ((1.0 + .5 * α) .* fₖ .- (.5 * α) .* f_prev_val)
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

function compute_lmm_step(::BDF2, f, zₖ, tₖ, Δt, x_history, t_history)
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