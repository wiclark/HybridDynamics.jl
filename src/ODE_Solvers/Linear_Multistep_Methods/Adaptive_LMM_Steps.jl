abstract type AdaptiveLMM <: LMM end

struct AdaptiveABM2 <: AdaptiveLMM end
struct AdaptiveABM3 <: AdaptiveLMM end

lmm_order(::AdaptiveABM2) = 2
lmm_order(::AdaptiveABM3) = 3

"""
    take_step(solver::AdaptiveLMM, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::RK=RK4())

Executes single variable-step Linear Multstep Method update for hybrid systems
How adaptive LMMs work:
Adaptivity requires us to know how much error we haeve when we make the current step so we can shrink or grow the step size h. LMMs do this using a Predictor/Corrector system. 
Predictor: Uses histoy points to project a rough guess for next state
Corrector: Takes that guess from predictor and refines it using the current dynamics, acting as a stabilizer.
Because predictor and corrector have known, and linked, error bounds, the difference between their outputs gives us an estimate of the Local Truncation Error (LTE).

If LTE is too high Δt is shrunk, and history is interpolated or reset, and we try to step again. If LTE is below our tolerance the step is accepted and solver calcs a slightly larger Δt for next step. 

How it works:
1) History Validation: The solver checks the length of continuous history since the last disc jump. If the history buffer is smaller than the order of the LMM 'k', 
the step is delegated to the single-step stepper (default RK4).

2)Predictor/Corrector and Error Est: Once history is established, the solver extracts the past `k` states and calls `compute_lmm_step`. This helper function executes the explicit prediction, the implicit correction, and extracts the isolated LTE using Milne's device.

3) Adaptivity Loop: Operates within 'while true' rejection loop. If LTE exceeds 'tol' the step size 'Δt' shrinks, and the step is recalculated. If accepted, it computes the optimal Δt_next for the subsequent step

WHY I DID THINGS:
We use Milne's device for error est because it is easy to compute. Since we already are performing an explicit prediction and an implicit correction, the difference serves as a solid estimate.
"""
#user can specify if they want RK4 here
function take_step(solver::AdaptiveLMM, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::RK=RK4(); check=true, guard_direction=default_guard_direction(prob.sys))
    
    # 1. Probing Override: If checking is disabled, LMM history assumptions are violated. Route to RK stepper.
    if !check
        return take_step(stepper, prob, f, xₖ, tₖ, Δt, tol, sol, stepper; check=false, guard_direction=guard_direction)
    end

    sys = prob.sys
    k = lmm_order(solver)
    tf = prob.tspan[2]

    Δt = minimum([Δt, tf - tₖ])

    # Determine how many cont steps we have since last jump
    history_len = isempty(sol.event_indices) ? length(sol.x) : (length(sol.x) - sol.event_indices[end])

    # Startup phase: IF we dont have history we use RK stepper
    if history_len <= k
        return take_step(stepper, prob, f, xₖ, tₖ, Δt, tol, sol, stepper; check=check, guard_direction=guard_direction)
    end

    # Extract history for LMM
    x_history = sol.x[end - k + 1 : end - 1]
    t_history = sol.t[end - k + 1 : end - 1]

    time_diffs = diff(vcat(t_history, tₖ))

    if any(time_diffs .<= 1e-12)
        return take_step(stepper, prob, f, xₖ, tₖ, Δt, tol, sol, stepper; check=check, guard_direction=guard_direction)
    end

    # Adaptive step loop
    while true
        # Compute step and retrieve LTE 
        x_next, x_predict, LTE = compute_lmm_step(solver, f, xₖ, tₖ, Δt, x_history, t_history)

        h_now = guard(sys, xₖ)
        h_corr = guard(sys, x_next)

        # Calc proposed next step size using helper from beginning 
        dt_next = 0.9*Δt*(tol/LTE)^(1/(k+1))

        # Prevent violent changes in step size.
        growth = 1.25
        shrink = 0.8

        dt_next = min(dt_next, growth*Δt)
        dt_next = max(dt_next, shrink*Δt)
        if LTE < tol
            if check
                h_now = guard(sys, xₖ)
                h_next = guard(sys, x_next)

                # Guard eval looking back to prev step for the quad check 
                if history_len > 1
                    idx = length(sol.x) - 1
                    t_prev = sol.t[idx]
                    h_prev = guard(sys, sol.x[idx])
                else
                    t_prev = tₖ - Δt
                    h_prev = h_now
                end
                
                eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt; tol=tol, direction=guard_direction)
                return x_next, eventtrigger, t_root, Δt, dt_next
            else 
                # Fallback (Structurally unreachable due to top delegation)
                return x_next, false, NaN, Δt, dt_next
            end
        else
            # Step rejected: shrink and try again 
            Δt = max(dt_next, 0.8*Δt)
            if Δt < 1e-6
                @warn "LMM Step size has decreased below 1e-6"
                # Force break to avoid inf loops
                return x_next, false, NaN, Δt, dt_next
            end
        end
    end
end

function compute_lmm_step(::AdaptiveABM2, f, xₖ, tₖ, Δt, x_history, t_history)
    x₋₁ = x_history[end]
    t₋₁ = t_history[end]

    fₖ  = f(xₖ, tₖ)
    f₋₁ = f(x₋₁, t₋₁)

    h = Δt
    h1 = tₖ - t₋₁

    # Variable-step AB2 predictor
    df1 = (fₖ .- f₋₁) ./ h1
    x_predict = xₖ .+ h .* fₖ .+ (h^2 / 2) .* df1

    f_predict = f(x_predict, tₖ + h)

    # Trapezoidal corrector (AM2) - inherently step-size independent
    x_correct = xₖ .+ (h / 2) .* (f_predict .+ fₖ)
    f_correct = f(x_correct, tₖ + h)

    x_correct2 = xₖ .+ (h / 2) .* (f_correct .+ fₖ)

    LTE = norm(x_correct2 .- x_predict) / max(norm(x_correct2), 1.0)

    return x_correct2, x_predict, LTE
end

function compute_lmm_step(::AdaptiveABM3, f, xₖ, tₖ, Δt, x_history, t_history)
    x_prev1 = x_history[end]
    t_prev1 = t_history[end]
    x_prev2 = x_history[end-1]
    t_prev2 = t_history[end-1]

    fₖ      = f(xₖ, tₖ)
    f_prev1 = f(x_prev1, t_prev1)
    f_prev2 = f(x_prev2, t_prev2)

    h = Δt
    h1 = tₖ - t_prev1
    h2 = t_prev1 - t_prev2

    # Divided differences for history
    df1 = (fₖ .- f_prev1) ./ h1
    df2 = (f_prev1 .- f_prev2) ./ h2
    D2_pred = (df1 .- df2) ./ (h1 + h2)

    # Variable-step AB3 predictor
    x_predict = xₖ .+ h .* fₖ .+ (h^2 / 2) .* df1 .+ (h^2 * (h/3 + h1/2)) .* D2_pred
    f_predict = f(x_predict, tₖ + h)

    # Variable-step AM3 corrector
    df_p = (f_predict .- fₖ) ./ h
    D2_corr = (df_p .- df1) ./ (h + h1)
    
    x_correct = xₖ .+ (h / 2) .* (f_predict .+ fₖ) .- (h^3 / 6) .* D2_corr
    f_correct = f(x_correct, tₖ + h)

    # second correction
    df_c = (f_correct .- fₖ) ./ h
    D2_corr2 = (df_c .- df1) ./ (h + h1)
    x_correct2 = xₖ .+ (h / 2) .* (f_correct .+ fₖ) .- (h^3 / 6) .* D2_corr2

    LTE = norm(x_correct2 .- x_predict) / max(norm(x_correct2), 1.0)

    return x_correct2, x_predict, LTE
end