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
    if history_len < k
        return take_step(stepper, prob, f, xₖ, tₖ, Δt, tol, sol, stepper; check=check, guard_direction=guard_direction)
    end

    # Extract history for LMM
    x_history = sol.x[end - k + 1 : end - 1]
    t_history = sol.t[end - k + 1 : end - 1]

    # Adaptive step loop
    while true
        # Compute step and retrieve LTE 
        x_predict, LTE = compute_lmm_step(solver, f, xₖ, tₖ, Δt, x_history, t_history)

        # Calc proposed next step size using helper from beginning 
        dt_next = updated_step(LTE, tol, Δt, k)

        if LTE < tol
            if check
                h_now = guard(sys, xₖ)
                h_next = guard(sys, x_predict)

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

                # Buffer
                if eventtrigger
                    if (t_root - tₖ) < (1e-4 * Δt)
                        eventtrigger = false
                        t_root = tₖ + Δt
                    end
                end

                return x_predict, eventtrigger, t_root, Δt, dt_next
            else 
                # Fallback (Structurally unreachable due to top delegation)
                return x_predict, false, NaN, Δt, dt_next
            end
        else
            # Step rejected: shrink and try again 
            Δt = dt_next
            if Δt < 1e-6
                @warn "LMM Step size has decreased below 1e-6"
                # Force break to avoid inf loops
                return x_predict, false, NaN, Δt, dt_next
            end
        end
    end
end

function compute_lmm_step(::AdaptiveABM2, f, xₖ, tₖ, Δt, x_history, t_history)
    #Extract history
    x_prev = x_history[end]
    t_prev = t_history[end]

    #Calc previous step size
    dt_prev = tₖ - t_prev
    α = Δt / dt_prev

    fₖ = f(xₖ, tₖ)
    f_prev = f(x_prev, t_prev)

    #Predictor for variable step AB2
    x_predict = xₖ .+ Δt .* ((1.0 + 0.5 * α) .* fₖ .- (0.5 * α) .* f_prev)

    #Eval vector field at pred state
    f_predict = f(x_predict, tₖ + Δt)

    #Corrector: Implicit AM2 via predicted vf
    x_correct = xₖ .+ Δt .* (0.5 .* f_predict .+ 0.5 .* fₖ)

    #Error est: Difference between corrector and predictor give local truncation error
    LTE = norm(x_correct .- x_predict)

    return x_correct, LTE
end

function compute_lmm_step(::AdaptiveABM3, f, xₖ, tₖ, Δt, x_history, t_history)
    #extract history
    x_prev1 = x_history[end]
    t_prev1 = t_history[end]

    x_prev2 = x_history[end - 1]
    t_prev2 = t_history[end - 1]

    #Eval vf at past known coords
    fₖ      = f(xₖ, tₖ)
    f_prev1 = f(x_prev1, t_prev1)
    f_prev2 = f(x_prev2, t_prev2)

    #Calc non uniform time grid ints
    hk = Δt                     #Current proposed step size t_{k+1} - t_k
    hk1 = tₖ - t_prev1           #Previous step size duration t_k - t_{k-1}
    hk2 = t_prev1 - t_prev2      #Two steps back duraction t_{k-1} - t_{k-2}

    #Predictor Step: Fully variable-step explicit AB3
    β₀ = hk * (hk^2 / 3.0 + hk * (2.0 * hk1 + hk2) / 2.0 + hk1 * (hk1 + hk2)) / (hk1 * (hk1 + hk2))
    β₁ = -hk^2 * (2.0 * hk + 3.0 * hk1 + 3.0 * hk2) / (6.0 * hk1 * hk2)
    β₂ = hk^2 * (2.0 * hk + 3.0 * hk1) / (6.0 * hk2 * (hk1 + hk2))

    x_predict = xₖ .+ (β₀ .* fₖ .+ β₁ .* f_prev1 .+ β₂ .* f_prev2)

    #Eval vf at predicted state
    f_predict = f(x_predict, tₖ + hk)

    #Corrector step: Adams Moulton 3
    c_β_p1 = hk * (hk / 3.0 + hk1 / 2.0) / (hk + hk1)
    c_β₀   = hk * (hk + 3.0 * hk1) / (6.0 * hk1)
    c_β_m1 = -hk^3 / (6.0 * hk1 * (hk + hk1))
    
    x_correct = xₖ .+ (c_β_p1 .* f_predict .+ c_β₀ .* fₖ .+ c_β_m1 .* f_prev1)

    #ERROR EST
    LTE = norm(x_correct .- x_predict)

    return x_correct, LTE
end