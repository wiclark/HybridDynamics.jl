abstract type AdaptiveLMM <: LMM end

"""
Adaptive second order Adams Bashforth Moulton multistep method. RK4 is used to caluclate starting values.
"""
struct AdaptiveABM2 <: AdaptiveLMM end
"""
Adaptive third order Adams Bashforth Moulton multistep method. RK4 is used to caluclate starting values.
"""
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
function take_step(solver::AdaptiveLMM, prob::AbstractHybridProblem, f, Df, xₖ, tₖ, Δt, tol, sol, stepper::RK=RK4(); check=true, guard_direction=default_guard_direction(prob.sys))
    
    # Probing Override: If checking is disabled, LMM history assumptions are violated. Route to RK stepper.
    if !check
        return take_step(stepper, prob, f, Df, xₖ, tₖ, Δt, tol, sol, stepper; check=false, guard_direction=guard_direction)
    end

    sys = prob.sys
    k = lmm_order(solver)
    tf = prob.tspan[2]

    Δt = minimum([Δt, tf - tₖ])

    # Determine how many cont steps we have since last jump
    history_len = isempty(sol.event_indices) ? length(sol.x) : (length(sol.x) - sol.event_indices[end])

    # Startup phase: IF we dont have history we use RK stepper
    if history_len <= k
        if length(sol.t) >= 2
            h_startup = sol.t[end] - sol.t[end - 1]

            if h_startup > 1e-12
                Δt = min(Δt, h_startup)
            end
        end
        return take_step(stepper, prob, f, Df, xₖ, tₖ, Δt, tol, sol, stepper; check=check, guard_direction=guard_direction)
    end

    # Extract history for LMM
    x_history = sol.x[end - k + 1 : end - 1]
    t_history = sol.t[end - k + 1 : end - 1]

    time_diffs = diff(vcat(t_history, tₖ))

    if any(time_diffs .<= 1e-12)
        return take_step(stepper, prob, f, Df, xₖ, tₖ, Δt, tol, sol, stepper; check=check, guard_direction=guard_direction)
    end

    h_last = time_diffs[end]

    if Δt > 2.0 * h_last
        Δt = 2.0 * h_last
    end

    # Adaptive step loop
    while true
        # Compute step and retrieve LTE 
        x_next, x_predict, LTE = compute_lmm_step(solver, f, Df, xₖ, tₖ, Δt, x_history, t_history)

        # Safeguard against numerical explosions with fancy guards or weird history details.  
        if any(isnan, x_next) || any(isinf, x_next) || any(isnan, x_predict) || any(isinf, x_predict)
            return take_step(stepper, prob, f, Df, xₖ, tₖ, Δt * 0.5, tol, sol, stepper; check=false, guard_direction=guard_direction)
        end

        # Guard boundary rejection logic. 
        h_now = guard(sys, xₖ)
        h_predict = guard(sys, x_predict)
        h_end = guard(sys, x_next)

        if !isnothing(h_now) && !isnothing(h_predict) && !isnothing(h_end)
            x_probe1 = (3.0 / 4.0) .* xₖ .+ (1.0 / 4.0) .* x_next
            x_probe2 = (1.0 / 2.0) .* xₖ .+ (1.0 / 2.0) .* x_next
            x_probe3 = (1.0 / 4.0) .* xₖ .+ (3.0 / 4.0) .* x_next

            h_probe1 = guard(sys, x_probe1)
            h_probe2 = guard(sys, x_probe2)
            h_probe3 = guard(sys, x_probe3)

            samples = (h_now, h_probe1, h_probe2, h_probe3, h_end)

            crossing_count = 0

            for i in 1:(length(samples) - 1)
                if samples[i] * samples[i + 1] < 0
                    crossing_count += 1
                end
            end

            multiple_crossings = crossing_count > 1

            # Predictor/end-point checks
            stage_crossed = (h_now * h_predict < 0)
            end_missed    = (h_now * h_end > 0)
            crossed       = (h_now * h_end < 0)

            # Check if system is a Filippov system using type name inspection
            # All this does is swap logic if we have a Filippov system to keep track of sliding mode entry and exit. 
            # Mechanical is the root of the issue as the logic each system needs doesnt fully match so this is what we have (IT SUCKS I KNOW - DS)
            is_filippov_sys = occursin("Filippov", string(typeof(sys)))
            should_enforce_penetration = !is_filippov_sys

            if multiple_crossings ||
               (stage_crossed && end_missed) ||
               (h_predict * h_end < 0) ||
               (should_enforce_penetration && crossed && abs(h_end) > tol)

                Δt = Δt / 2.0
                continue
            end
        end

        # Prevent violent changes in step size.
        growth = 2.0
        shrink = 0.25

        # Calc proposed next step size using helper from beginning 
        if LTE == 0
            dt_next = growth * Δt
        else
            dt_next = 0.9 * Δt * (tol / LTE)^(1 / (k + 1))
        end

        dt_next = min(dt_next, growth * Δt)
        dt_next = max(dt_next, shrink * Δt)

        if LTE < tol
            if check
                # Evaluate looking back to previous accepted step
                # for the event locator.
                if length(sol.x) >= 2
                    idx = max(1, length(sol.x) - 1)

                    t_prev = sol.t[idx]
                    h_prev = guard(sys, sol.x[idx])

                    # Prevent degenerate quadratic interpolation if things get small. 
                    if (tₖ - t_prev) > 10 * Δt
                        t_prev = tₖ - Δt
                        # Linear backward extrapolation gives the locator a stable local bracket.
                        h_prev = h_now - (h_end - h_now)
                    end
                else
                    t_prev = tₖ - Δt
                    h_prev = h_now
                end

                eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_end, t_prev, tₖ, tₖ + Δt; tol=tol, direction=guard_direction)
                return x_next, eventtrigger, t_root, Δt, dt_next
            else
                return x_next, false, NaN, Δt, dt_next
            end
        else
            # Step rejected due to LTE error: shrink and retry
            Δt = dt_next

            if Δt < 1e-12
                @warn "LMM Step size has decreased below 1e-12"
                return x_next, false, NaN, Δt, dt_next
            end
        end
    end
end

function compute_lmm_step(::AdaptiveABM2, f, Df, xₖ, tₖ, Δt, x_history, t_history)
    x_prev = x_history[end]
    t_prev = t_history[end]

    fₖ = f(xₖ, tₖ)
    f_prev = f(x_prev, t_prev)

    ab_times = [t_prev, tₖ]
    # Pass tₖ and Δt to integrate from tₖ to tₖ + Δt
    ab_coeff = variable_adams_coefficients(ab_times, tₖ, Δt)

    # AB Predictor
    x_predict = xₖ .+ Δt .* (ab_coeff[2] .* fₖ .+ ab_coeff[1] .* f_prev)

    # AM corrector
    f_predict = f(x_predict, tₖ + Δt)

    am_times = [tₖ, tₖ + Δt]
    # Pass tₖ and Δt to integrate from tₖ to tₖ + Δt
    am_coeff = variable_adams_coefficients(am_times, tₖ, Δt)
    
    x_correct = xₖ .+ Δt .* (am_coeff[2] .* f_predict .+ am_coeff[1] .* fₖ)

    # Milne LTE estimate
    scale = max(norm(x_correct), 1.0)
    LTE = (1 / 6) * norm(x_correct .- x_predict) / scale

    return x_correct, x_predict, LTE
end

function compute_lmm_step(::AdaptiveABM3, f, Df, xₖ, tₖ, Δt, x_history, t_history)
    x_prev1 = x_history[end]
    t_prev1 = t_history[end]

    x_prev2 = x_history[end-1]
    t_prev2 = t_history[end-1]

    fₖ = f(xₖ, tₖ)
    f_prev1 = f(x_prev1, t_prev1)
    f_prev2 = f(x_prev2, t_prev2)

    # AB predictor
    ab_times = [t_prev2, t_prev1, tₖ]
    ab_coeff = variable_adams_coefficients(ab_times, tₖ, Δt)

    x_predict = xₖ .+ Δt .* (ab_coeff[3] .* fₖ .+ ab_coeff[2] .* f_prev1 .+ ab_coeff[1] .* f_prev2)

    # AM corrector
    f_predict = f(x_predict, tₖ+Δt)

    am_times = [t_prev1, tₖ, tₖ + Δt]
    am_coeff = variable_adams_coefficients(am_times, tₖ, Δt)

    x_correct = xₖ .+ Δt .* (am_coeff[3] .* f_predict .+ am_coeff[2] .* fₖ .+ am_coeff[1] .* f_prev1)

    # Milne LTE estimate
    scale = max(norm(x_correct), 1.0)
    LTE = (1 / 10) * norm(x_correct .- x_predict) / scale

    return x_correct, x_predict, LTE
end

# Compute Adams coeffs for arbitrary time spacing and target step size
function variable_adams_coefficients(times, t_n, Δt)
    n = length(times)
    
    # Normalize the time points based on the starting time t_n and step size Δt
    ξ = (times .- t_n) ./ Δt

    coeffs = zeros(n)

    for i in 1:n
        # Build Lagrange basis polynomial
        poly = [1.0]
        denom = 1.0

        for j in 1:n
            if j != i
                denom *= (ξ[i] - ξ[j])

                newpoly = zeros(length(poly)+1)

                for k in 1:length(poly)
                    newpoly[k] += -ξ[j]*poly[k]
                    newpoly[k+1] += poly[k]
                end

                poly = newpoly
            end
        end

        poly ./= denom

        # Integrate from 0 to 1 (which corresponds to t_n to t_n + Δt)
        integral = 0.0
        for k in 1:length(poly)
            integral += poly[k]/k
        end

        coeffs[i] = integral
    end

    return coeffs
end