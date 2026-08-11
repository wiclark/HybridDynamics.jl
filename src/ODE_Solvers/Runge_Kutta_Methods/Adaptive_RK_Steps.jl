abstract type AdaptiveRK <: RK end

#Adaptive solvers
"""
Adaptive fourth-fifth order Runge-Kutta method.
"""
struct RK45 <: AdaptiveRK end
"""
Adaptive second-third order Runge-Kutta method.
"""
struct RK23 <: AdaptiveRK end

#Adaptive step methods
#Helper function to take the step via multiple dispatch
compute_step(::RK23, f, x, Δt, t, tf, sys, tol; adaptive=true) = rk_23_step(f, x, Δt, t, tf, sys, tol)
compute_step(::RK45, f, x, Δt, t, tf, sys, tol; adaptive=true) = rk_45_step(f, x, Δt, t, tf, sys, tol)

function take_step(solver::AdaptiveRK, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver=ModifiedMidpoint();
        check=true, guard_direction=default_guard_direction(prob.sys))
    sys = prob.sys
    tf = prob.tspan[2] #terminal time

    #Take adaptive step (passes tf to prevent overshooting)
    x_predict, Δt, dt_next = compute_step(solver, f, xₖ, Δt, tₖ, tf, sys, tol)

    if check
        # Evaluate guards
        h_now  = guard(sys, xₖ)
        h_next = guard(sys, x_predict)

        idx = max(1, length(sol.x) - 1)
        t_prev = sol.t[idx]
        x_prev = sol.x[idx]
        h_prev = guard(sys, x_prev)
        # Use cross guard check
        eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt; tol=tol, direction=guard_direction)

        if eventtrigger
            if (t_root - tₖ) < (1e-4 * Δt)
                eventtrigger = false
                t_root = tₖ + Δt
            end
        end

        return x_predict, eventtrigger, t_root, Δt, dt_next
    else
        return x_predict, false, NaN, Δt, dt_next
    end
end

# A helper function to determine the adapted step size
function updated_step(LTE::AbstractFloat, tol::AbstractFloat, Δt::AbstractFloat, n::Int)
    # The safety parameters
    facmax = 3.
    facmin = 1/3
    fac    = 0.9
    # The predicted multiplier
    ε = abs( tol / LTE ) ^ (1/n)
    # The updated step
    return Δt * minimum( [ facmax, maximum( [ facmin, fac*ε ] ) ] )
end

# Runge-Kutta 23 (Bogacki-Shampine 3(2))
function rk_23_step(f::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat, tf::AbstractFloat, sys, tol)
    dt_step = minimum([Δt, tf - t])

    x2 = xₖ; x3 = xₖ
    while true
        h_now = guard(sys, xₖ)

        # Stage 1
        k1 = f(xₖ, t)

        # Stage 2
        x2 = xₖ + dt_step * (1/2) * k1
        k2 = f(x2, t + dt_step * 1/2)
        h2 = guard(sys, x2)

        # Stage 3
        x3 = xₖ + dt_step * (3/4) * k2
        k3 = f(x3, t + dt_step * 3/4)
        h3 = guard(sys, x3)

        # Compute 3rd-order state y3
        y3 = xₖ + dt_step * ((2/9) * k1 + (1/3) * k2 + (4/9) * k3)

        # Evaluate k4 at y3
        k4 = f(y3, t + dt_step)

        # Direct error calculation: err = y3 - y2
        err = dt_step * ((-5/72) * k1 + (1/12) * k2 + (1/9) * k3 - (1/8) * k4)

        # Propagate 2nd-order solution y2 = y3 - err 
        x_predict = y3 - err

        scale = max(norm(x_predict), 1.0)
        LTE = norm(err) / scale

        dt_next = updated_step(LTE, tol, dt_step, 3)

        h_end = guard(sys, x_predict)

        stage_crossed = (h_now * h2 < 0) || (h_now * h3 < 0)
        end_missed = (h_now * h_end > 0)
        crossed = (h_now * h_end < 0)

        # Step rejection
        if (stage_crossed && end_missed) || (crossed && abs(h_end) > max(abs(h_now), 10 * tol))
            dt_step = dt_step / 2.0
            continue
        end

        if LTE < tol
            return x_predict, dt_step, dt_next
        else
            dt_step = dt_next
        end

        if dt_step < 1e-12
            @warn "Step size has decreased below 1e-12"
            return x_predict, dt_step, dt_next
        end
    end
end

# Runge-Kutta 45 (Dormand-Prince 5(4))
function rk_45_step(f::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat, tf::AbstractFloat, sys, tol)
    dt_step = minimum([Δt, tf - t])

    x2 = xₖ; x3 = xₖ; x4 = xₖ; x5 = xₖ; x6 = xₖ
    while true
        h_now = guard(sys, xₖ)

        # Stage 1
        k1 = f(xₖ, t)

        # Stage 2
        x2 = xₖ + dt_step * (1/5) * k1
        k2 = f(x2, t + dt_step * 1/5)
        h2 = guard(sys, x2)

        # Stage 3
        x3 = xₖ + dt_step * (3/40 * k1 + 9/40 * k2)
        k3 = f(x3, t + dt_step * 3/10)
        h3 = guard(sys, x3)

        # Stage 4
        x4 = xₖ + dt_step * (44/45 * k1 - 56/15 * k2 + 32/9 * k3)
        k4 = f(x4, t + dt_step * 4/5)
        h4 = guard(sys, x4)

        # Stage 5
        x5 = xₖ + dt_step * (19372/6561 * k1 - 25360/2187 * k2 + 64448/6561 * k3 - 212/729 * k4)
        k5 = f(x5, t + dt_step * 8/9)
        h5 = guard(sys, x5)

        # Stage 6
        x6 = xₖ + dt_step * (9017/3168 * k1 - 355/33 * k2 + 46732/5247 * k3 + 49/176 * k4 - 5103/18656 * k5)
        k6 = f(x6, t + dt_step)
        h6 = guard(sys, x6)

        # Compute 5th-order state y5
        y5 = xₖ + dt_step * ((35/384) * k1 + (500/1113) * k3 + (125/192) * k4 - (2187/6784) * k5 + (11/84) * k6)

        # Evaluate k7 at y5
        k7 = f(y5, t + dt_step)

        # Direct error calculation: err = y5 - y4
        err = dt_step * ((71/57600) * k1 - (71/16695) * k3 + (71/1920) * k4 - (17253/339200) * k5 + (22/525) * k6 - (1/40) * k7)

        # Propagate 4th-order solution y4 = y5 - err
        x_predict = y5 - err

        scale = max(norm(x_predict), 1.0)
        LTE = norm(err) / scale

        dt_next = updated_step(LTE, tol, dt_step, 5)

        h_end = guard(sys, x_predict)

        stage_crossed = (h_now * h2 < 0) || (h_now * h3 < 0) || (h_now * h4 < 0) || (h_now * h5 < 0) || (h_now * h6 < 0)
        end_missed = (h_now * h_end > 0)
        crossed = (h_now * h_end < 0)

        # Step rejection
        if (stage_crossed && end_missed) || (crossed && abs(h_end) > max(abs(h_now), 10 * tol))
            dt_step = dt_step / 2.0
            continue
        end
        
        if LTE < tol
            return x_predict, dt_step, dt_next
        else
            dt_step = dt_next
        end

        if dt_step < 1e-12
            @warn "Step size has decreased below 1e-12"
            return x_predict, dt_step, dt_next
        end
    end
end
