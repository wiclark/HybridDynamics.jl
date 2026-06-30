abstract type AdaptiveRK <: RK end

#Adaptive solvers
struct RK45 <: AdaptiveRK end
struct RK23 <: AdaptiveRK end

#Adaptive step methods
#Helper function to take the step via multiple dispatch
compute_step(::RK23, f, x, Δt, t, tf, sys, tol; adaptive=true) = rk_23_step(f, x, Δt, t, tf, sys, tol; adaptive=adaptive)
compute_step(::RK45, f, x, Δt, t, tf, sys, tol; adaptive=true) = rk_45_step(f, x, Δt, t, tf, sys, tol; adaptive=adaptive)

function take_step(solver::AdaptiveRK, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver=ModifiedMidpoint();
        check=true, guard_direction=default_guard_direction(prob.sys))
    sys = prob.sys
    tf = prob.tspan[2] #terminal time

    #Take adaptive step (passes tf to prevent overshooting)
    ######
    ### WC: Why was adaptive=false here?
    ###
    ### DS: Not supposed to be, good catch. 
    ######
    x_predict, dt_used, dt_next = compute_step(solver, f, xₖ, Δt, tₖ, tf, sys, tol)
    # Get midpoint for quad guard check
    # For a half-step the local ceiling is just midpoint of time 
    # Disable adaptive loop to get the true midpoint. Without it it would reset the loop and ruin or previous adaptivity
    x_mid, _, _ = compute_step(solver, f, xₖ, dt_used / 2.0, tₖ, tf, sys, tol; adaptive=false)

    if check
        # Evaluate guards
        h_now  = guard(sys, xₖ)
        h_mid  = guard(sys, x_mid)
        h_next = guard(sys, x_predict)
        # Use cross guard check
        eventtrigger, t_root, _ = crossed_guard(sys, h_now, h_mid, h_next, tₖ, tₖ + dt_used / 2.0, tₖ + dt_used; tol=tol, direction=guard_direction)

        return x_predict, eventtrigger, t_root, dt_used, dt_next
    else
        return x_predict, false, NaN, Δt, Δt
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
function rk_23_step(f::Function, z::AbstractArray, Δt::AbstractFloat, t::AbstractFloat, tf::AbstractFloat, sys, tol; adaptive=true)
    # As this is an adaptive step solver, Δt is the step size from the previous step
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    dt = minimum([Δt, tf-t])

    #Intialization to make sure they exist
    z2 = z; z3 = z

    # Loop through to find an acceptable step
    while true
        h_now = guard(sys, z)
        # Compute the two predictions and their difference
        k1 = f(z, t)

        z2 = z + dt*k1
        k2 = f(z2, t+dt)
        h2 = guard(sys, z2)

        z3 = z + dt/4*(k1+k2)
        k3 = f(z3, t+dt/2)
        h3 = guard(sys, z3)

        z1_3 = z + dt*(1/6*k1+1/6*k2+2/3*k3)
        z1_2 = z + dt*(1/2*k1+1/2*k2)

        #used to find true midpoint for quad guard check. 
        #This prevents solver from spinning up a whole new adaptive loop so x_mid is truly in the mid which is what we want for quad matrix. 
        if !adaptive
            return z1_2, dt, dt
        end

        LTE = norm(z1_2 - z1_3)
        # Reject or accept?
        dt_new = updated_step(LTE, tol, dt, 3)
        
        h_end = guard(sys, z1_2)

        #did any intermediate stage cross guard?
        stage_crossed = (h_now * h2 < 0) || (h_now * h3 < 0)
        #Did final state completely miss crossing?
        end_missed = (h_now * h_end > 0)

        #if we crossed inside step but missed at end, force rejection
        if stage_crossed && end_missed
            dt = dt / 2.0 #force smaller step
            continue
        end

        if LTE < tol
            return z1_2, dt, dt_new
        else
            # We reject and repeat the loop with an updated step
            dt = dt_new
        end
        if dt < 1e-12
            @warn "Step size has decreased below 1e-12"
            return z1_2, dt, dt_new #force break to avoid looping forever
        end
    end
end

# Runge-Kutta 45
function rk_45_step(f::Function, z::AbstractArray, Δt::AbstractFloat, t::AbstractFloat, tf::AbstractFloat, sys, tol; adaptive=true)
    # As this is an adaptive step solver, h is the step size from the pervious step
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    dt = minimum([Δt, tf-t])

    #Initialization to make them exist
    z2 = z; z3 = z; z4 = z; z5 = z; z6 = z

    # Loop through to find an acceptable step
    while true
        h_now = guard(sys, z)
        # Compute the two predictions and their difference
        k1 = f(z, t)

        z2 = z + dt*1/5*k1
        k2 = f(z2, t+dt*1/5)
        h2 = guard(sys, z2)

        z3 = z + dt*(3/40*k1 + 9/40*k2)
        k3 = f(z3, t+dt*3/10)
        h3 = guard(sys, z3)

        z4 = z + dt*(44/45*k1 - 56/15*k2 + 32/9*k3)
        k4 = f(z4, t+dt*4/5)
        h4 = guard(sys, z4)

        z5 = z + dt*(19372/6561*k1 - 25360/2187*k2 + 64448/6561*k3 - 212/729*k4)
        k5 = f(z5,t+dt*8/9)
        h5 = guard(sys, z5)

        z6 = z + dt*(9017/3168*k1 - 355/33*k2 + 46732/5247*k3 + 49/176*k4 - 5105/18656*k5)
        k6 = f(z6, t+dt)
        h6 = guard(sys, z6)

        k7 = f(z + dt*(35/384*k1 + 0*k2 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6), t+dt)

        # The two updates
        z1_4 = z + dt*k7
        z1_5 = z + dt*(5179/57600*k1 + 0*k2 + 7571/16695*k3 + 393/640*k4 - 92097/339200*k5 + 187/2100*k6 + 1/40*k7)

        #used to find true midpoint for quad guard check. 
        #This prevents solver from spinning up a whole new adaptive loop so x_mid is truly in the mid which is what we want for quad matrix. 
        if !adaptive
            return z1_4, dt, dt
        end

        LTE = norm(z1_4 - z1_5)
        # Reject or accept?
        dt_new = updated_step(LTE, tol, dt, 5)
        
        h_end = guard(sys, z1_4)

        stage_crossed = (h_now * h2 < 0) || (h_now * h3 < 0) || (h_now * h4 < 0) || (h_now * h5 < 0) || (h_now * h6 < 0)
        end_missed = (h_now * h_end > 0)

        if stage_crossed && end_missed
            dt = dt / 2.0 #force smaller step
            continue
        end
        

        if LTE < tol
            return z1_4, dt, dt_new
        else
            dt = dt_new
        end
        if dt < 1e-12
            @warn "Step size has decreased below 1e-12"
            return z1_4, dt, dt_new
        end
    end
end
