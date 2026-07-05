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
                @show Δt dt_next LTE
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

    # history
    xₖ₋₁ = x_history[end]
    tₖ₋₁ = t_history[end]

    # slopes
    fₖ   = f(xₖ, tₖ)
    fₖ₋₁ = f(xₖ₋₁, tₖ₋₁)

    x_predict = xₖ .+ Δt .* ( (3/2).*fₖ .- (1/2).*fₖ₋₁ )

    f_predict = f(x_predict, tₖ + Δt)

    x_correct = xₖ .+ (Δt/2) .* (f_predict .+ fₖ)

    f_correct = f(x_correct, tₖ + Δt)

    x_correct2 = xₖ .+ (Δt/2) .* (f_correct .+ fₖ)

    LTE = norm(x_correct2 .- x_predict) /
          max(norm(x_correct2), 1.0)

    return x_correct2, x_predict, LTE
end

function compute_lmm_step(::AdaptiveABM3, f, xₖ, tₖ, Δt, x_history, t_history)

    x_prev1 = x_history[end]
    t_prev1 = t_history[end]

    x_prev2 = x_history[end-1]
    t_prev2 = t_history[end-1]

    fₖ      = f(xₖ,tₖ)
    f_prev1 = f(x_prev1,t_prev1)
    f_prev2 = f(x_prev2,t_prev2)

    x_predict = xₖ .+ Δt .* ((23/12).*fₖ .- (16/12).*f_prev1 .+  (5/12).*f_prev2)

    f_predict = f(x_predict,tₖ+Δt)

    x_correct = xₖ .+ Δt .* ((5/12).*f_predict.+ (8/12).*fₖ.- (1/12).*f_prev1)

    f_correct = f(x_correct,tₖ+Δt)

    x_correct2 = xₖ .+ Δt .* ((5/12).*f_correct .+ (8/12).*fₖ .- (1/12).*f_prev1)

    LTE = norm(x_correct2 .- x_predict) / max(norm(x_correct2),1.0)

    return x_correct2, x_predict, LTE
end