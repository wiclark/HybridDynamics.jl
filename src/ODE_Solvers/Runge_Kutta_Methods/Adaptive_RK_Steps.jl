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

# Runge-Kutta 23
function rk_23_step(f::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat, tf::AbstractFloat, sys, tol)
    dt_step = minimum([Δt, tf-t])

    while true
        h_now = guard(sys, xₖ)
        
        k1 = f(xₖ, t)
        
        x2 = xₖ .+ dt_step .* (1/2) .* k1
        k2 = f(x2, t + dt_step/2)
        h2 = guard(sys, x2)
        
        x3 = xₖ .+ dt_step .* (3/4) .* k2
        k3 = f(x3, t + dt_step*3/4)
        h3 = guard(sys, x3)
        
        # 3rd-order advance
        x_predict = xₖ .+ dt_step .* (2/9 .* k1 .+ 1/3 .* k2 .+ 4/9 .* k3)
        h_end = guard(sys, x_predict)
        
        k4 = f(x_predict, t + dt_step) # FSAL setup
        
        # 2nd-order error estimate
        x_err = xₖ .+ dt_step .* (7/24 .* k1 .+ 1/4 .* k2 .+ 1/3 .* k3 .+ 1/8 .* k4)
        
        LTE = norm(x_predict .- x_err)
        
        # Error is bounded by the 2nd-order estimate, so it scales as O(h^3)
        dt_next = updated_step(LTE, tol, dt_step, 3) 
        
        stage_crossed = (h_now * h2 < 0) || (h_now * h3 < 0)
        end_missed = (h_now * h_end > 0)
        
        if stage_crossed && end_missed
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

# Runge-Kutta 45
function rk_45_step(f::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat, tf::AbstractFloat, sys, tol)
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    dt_step = minimum([Δt, tf-t])

    #Initialization to make them exist
    x2 = xₖ; x3 = xₖ; x4 = xₖ; x5 = xₖ; x6 = xₖ

    # Loop through to find an acceptable step
    while true
        h_now = guard(sys, xₖ)
        # Compute the two predictions and their difference
        k1 = f(xₖ, t)

        x2 = xₖ + dt_step*1/5*k1
        k2 = f(x2, t+dt_step*1/5)
        h2 = guard(sys, x2)

        x3 = xₖ + dt_step*(3/40*k1 + 9/40*k2)
        k3 = f(x3, t+dt_step*3/10)
        h3 = guard(sys, x3)

        x4 = xₖ + dt_step*(44/45*k1 - 56/15*k2 + 32/9*k3)
        k4 = f(x4, t+dt_step*4/5)
        h4 = guard(sys, x4)

        x5 = xₖ + dt_step*(19372/6561*k1 - 25360/2187*k2 + 64448/6561*k3 - 212/729*k4)
        k5 = f(x5,t+dt_step*8/9)
        h5 = guard(sys, x5)

        x6 = xₖ + dt_step*(9017/3168*k1 - 355/33*k2 + 46732/5247*k3 + 49/176*k4 - 5105/18656*k5)
        k6 = f(x6, t+dt_step)
        h6 = guard(sys, x6)

        x_predict = xₖ + dt_step*(35/384*k1 + 0*k2 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6)

        k7 = f(x_predict, t+dt_step)

        x1_5 = xₖ + dt_step*(5179/57600*k1 + 0*k2 + 7571/16695*k3 + 393/640*k4 - 92097/339200*k5 + 187/2100*k6 + 1/40*k7)

        LTE = norm(x_predict - x1_5)
        # Reject or accept?
        dt_next = updated_step(LTE, tol, dt_step, 5)
        
        h_end = guard(sys, x_predict)

        stage_crossed = (h_now * h2 < 0) || (h_now * h3 < 0) || (h_now * h4 < 0) || (h_now * h5 < 0) || (h_now * h6 < 0)
        end_missed = (h_now * h_end > 0)

        if stage_crossed && end_missed
            dt_step = dt_step / 2.0 #force smaller step
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
