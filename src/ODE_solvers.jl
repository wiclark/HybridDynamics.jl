# A collection of miscellaneous ODE integrators
# For all that follows:
#  1) f::Function is the vector field
#  2) z::Vector is the current state
#  3) h::Float is the step size (if fixed)
#  4) t::Float is the current time

#Parent Abstract Type
abstract type AbstractODESolver end

#Runge Kutta / Single Step Family Solvers
abstract type RK <: AbstractODESolver end

struct ForwardEuler <: RK end
struct ModifiedTrap <: RK end
struct ModifiedMidpoint <: RK end
struct RichardsonExtrapolation <: RK end

#Exponential Solver
struct ExponentialSolver <: AbstractODESolver end

#Linear Multistep Method Family Solvers
abstract type LMM <: AbstractODESolver end

struct AdamsBashforth2 <: LMM end
struct AdamsBashforth3 <: LMM end


## Single step, fully explicit methods

# Forward Euler
function forward_euler_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat)
    return z .+ h*f(z, t)
end

#NEW FORWARD EULER 
function forward_euler(f::Function, u0, tspan::Tuple{Float64, Float64}; dt::Float64 = .01)
    t_start, t_end = tspan

    t = Float64[t_start]
    u = typeof(u0)[u0]

    current_t = t_start
    current_u = u0

    while current_t < t_end
        #Ensure last step doesnt overshoot t_end
        step_size = min(dt, t_end - current_t)

        #Calling step function
        current_u = forward_euler_step(f, current_u, step_size, current_t)
        current_t += step_size

        push!(t, current_t)
        push!(u, current_u)
    end
    return t, u
end

# Modified (fully explicit) Trapezoid Rule
function modified_trap_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat)
    z_guess = z .+ h*f(z,t)
    return z .+ 1/2*h*( f(z, t) + f(z_guess, t+h) )
end

# Modified (fully explicit) Midpoint Rule
function modified_midpoint_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat)
    z_guess = z .+ h/2*f(z, t)
    return z .+ h*f(z_guess, t+h/2)
end

## Adaptive Runge-Kutta methods

# A helper function to determine the adapted step size
function updated_step(LTE::AbstractFloat, tol::AbstractFloat, h::AbstractFloat, n::Int)
    # The safety parameters
    facmax = 3.
    facmin = 1/3
    fac    = 0.9
    # The predicted multiplier
    ε = abs( tol / LTE ) ^ (1/n)
    # The updated step
    return h * minimum( [ facmax, maximum( [ facmin, fac*ε ] ) ] )
end

# Runge-Kutta 23
function rk_23_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat, tf::AbstractFloat)
    # As this is an adaptive step solver, h is the step size from the previous step
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    h = minimum([h, tf-t])
    # Set the tolerence
    tol = 1e-4

    # Loop through to find an acceptable step
    while true
        # Compute the two predictions and their difference
        k1 = f(z, t)
        k2 = f(z+h*k1, t+h)
        k3 = f(z+h/4*(k1+k2), t+h/2)
        z1_3 = z + h*(1/6*k1+1/6*k2+2/3*k3)
        z1_2 = z + h*(1/2*k1+1/2*k2)
        LTE = norm(z1_2 - z1_3)
        # Reject or accept?
        h_new = updated_step(LTE, tol, h, 3)
        if LTE < tol
            return z1_2, h, h_new
        else
            # We reject and repeat the loop with an updated step
            h = h_new
        end
        if h < 1e-12
            @warn "Step size has decreased below 1e-12"
        end
    end
end

# Runge-Kutta 45
function rk_45_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat, tf::AbstractFloat)
    # As this is an adaptive step solver, h is the step size from the pervious step
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    h = minimum([h, tf-t])
    # Set the tolerence
    tol = 1e-4

    # Loop through to find an acceptable step
    while true
        # Compute the two predictions and their difference
        k1 = f(z, t)
        k2 = f(z+h*1/5*k1, t+h*1/5)
        k3 = f(z+h*(3/40*k1+9/40*k2), t+h*3/10)
        k4 = f(z+h*(44/45*k1-56/15*k2+32/9*k3), t+h*4/5)
        k5 = f(z+h*(19372/6561*k1-25360/2187*k2+64448/6561*k3-212/729*k4),t+h*8/9)
        k6 = f(z+h*(9017/3168*k1-355/33*k2+46732/5247*k3+49/176*k4-5105/18656*k5), t+h)
        k7 = f(z+h*(35/384*k1+0*k2+500/1113*k3+125/192*k4-2187/6784*k5+11/84*k6),t+h)
        # The two updates
        z1_4 = z + h*k7
        z1_5 = z + h*(5179/57600*k1 + 0*k2 + 7571/16695*k3 + 393/640*k4 - 92097/339200*k5 + 187/2100*k6 + 1/40*k7)
        LTE = norm(z1_4 - z1_5)
        # Reject or accept?
        h_new = updated_step(LTE, tol, h, 5)
        if LTE < tol
            return z1_4, h, h_new
        else
            h = h_new
        end
        if h < 1e-12
            @warn "Step size has decreased below 1e-12"
        end
    end
end

#Richardson Extrapolation based on modifed midpoint from Wiki
function richardson_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat)
    #Take one full step size h
    z1 = modified_midpoint_step(f,z,h,t)

    #Take two smaller steps of h/2
    h_half = h / 2
    z_half = modified_midpoint_step(f,z,h_half,t)
    z2 = modified_midpoint_step(f, z_half,h_half,t+h_half)

    #Extrapolate to get rid of lower order error
    return (4 .* z2 .- z1) ./ 3.0
end

#====================================#
#TAKE STEP SOLVERS CONSOLIDATION

#Helper function to take the step using multiple dispatch.
compute_step(::ForwardEuler, f, x, Δt, t) = forward_euler_step(f, x, Δt, t)
compute_step(::ModifiedTrap, f, x, Δt, t) = modified_trap_step(f, x, Δt, t)
compute_step(::ModifiedMidpoint, f, x, Δt, t) = modified_midpoint_step(f, x, Δt, t)
compute_step(::RichardsonExtrapolation, f, x, Δt, t) = richardson_step(f, x, Δt, t)

#Single take_step for RK methods
#Note sol is not used, we do this to make using the function easier. We would need an if/else statement everytime we use this function without it
function take_step(solver::Union{RK, ExponentialSolver}, sys, f, xₖ, tₖ, Δt, tol, sol) 
    #Handle Exp vs RK math
    if solver isa ExponentialSolver
        flowmap   = LinearFlow(sys.A)
        x_predict = flow(flowmap, Δt, xₖ)
        x_mid     = flow(flowmap, Δt / 2.0, xₖ)
    else 
        x_predict = compute_step(solver, f, xₖ, Δt, tₖ)
        x_mid     = compute_step(solver, f, xₖ, Δt / 2.0, tₖ)
    end

    #Evaluate Guards
    h_now  = guard(sys, xₖ)
    h_mid  = guard(sys, x_mid)
    h_next = guard(sys, x_predict)

    #Use cross guard check
    eventtrigger, dt_next, _ = crossed_guard_will(h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol)

    return x_predict, eventtrigger, h_now, dt_next
end


#Helper function to tell engine how many history steps are needed for LMM
lmm_order(::AdamsBashforth2) = 2
lmm_order(::AdamsBashforth3) = 3

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
    f_prev1_val = f(x_prev1, t_prev1)
    f_prev2_val = f(x_prev2, t_prev2)
    
    # Currently using fixed-step coefficients. 
    # To upgrade to full variable-step AB3 later, you would calculate the 
    # ratios between (tₖ - t_prev1) and (t_prev1 - t_prev2) to dynamically adjust these weights.
    #This will be done eventually but for now this is fine. 
    return xₖ .+ Δt .* ( (23/12) .* fₖ .- (16/12) .* f_prev1_val .+ (5/12) .* f_prev2_val )
end

#One take_step for LMM
function take_step(solver::LMM, sys, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver = ModifiedTrap())
    k = lmm_order(solver)

    #Determine how many continuous steps we have since the last jump
    history_len = isempty(sol.jump_indices) ? length(sol.x) : (length(sol.x) - sol.jump_indices[end] + 1)

    h_now = guard(sys, xₖ)

    if history_len < k
        #Startup phase: Use single step predictor
        #using midpoint here to feed the quadratic guard check
        x_predict = compute_step(stepper, f, xₖ, Δt, tₖ)
        x_mid     = compute_step(stepper, f, xₖ, Δt / 2.0, tₖ)

        h_mid  = guard(sys, x_mid)
        h_next = guard(sys, x_predict)

        eventtrigger, dt_next, _ = crossed_guard_will(h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol)
        return x_predict, eventtrigger, h_now, dt_next
    else
        #Multistep phase: We do have rich enough history. Extract past states
        #If sol.x[end] is xₖ, then sol.x[end-1] is x_{k-1}
        x_history = sol.x[end - k + 1 : end - 1]
        t_history = sol.t[end - k + 1 : end - 1]

        #pass history arrays forward
        x_predict = compute_lmm_step(solver, f, xₖ, tₖ, Δt, x_history, t_history)
        
        h_next = guard(sys, x_predict)

        #Clark Fix: No calc of midpoint. 
        #Look back to the prev step evaluation to build quad guard
        t_prev = sol.t[end-1]
        x_prev = sol.x[end-1]
        h_prev = guard(sys, x_prev)

        eventtrigger, dt_next, _ = crossed_guard_will(h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt; tol = tol)

        return x_predict, eventtrigger, h_now, dt_next
    end

    #Guard Eval
    h_now  = guard(sys, xₖ)
    h_mid  = guard(sys, x_mid)
    h_next = guard(sys, x_predict)

    eventtrigger, dt_next, _ = crossed_guard_will(h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol)

    return x_predict, eventtrigger, h_now, dt_next
end



#Event detection utility. 
#If a guard surface was crossed and during the ODE step. We check for a sign change between start and end of the step.

######
### WC: Your parabolic version is incorrect. You do not necessarily know that the points have a uniform Δt (e.g., RK45). 
###     See my version below:
######

### WC: My version: Returns true/false along with the predicted impact time
###     I am padding the returns by NaN if I don't care about that result.
function crossed_guard_will(h_prev, h_now, h_next, t_prev, t_now, t_next; tol=1e-6)
    # I'm not sure why you have this, so I'll keep it anyway
    if isnothing(h_now)
        return false, t_next-t_now, NaN
    end

    # Performing the linear version. Why is the 'or' part required? You have the absolute value of now, but not for the next?
    if (h_now * h_next < 0)
        return true, t_now-h_now*(t_next-t_now)/(h_next-h_now), NaN
    end

    # Performing the quadratic version. If the discriminant is positive, there are roots.
    # Recall, that by the IVT, the linear test guarantees a crossing. The quadratic test does not guarantee one. This triggering should be treated as a warning.
    a, b, c = [t_prev^2 t_prev 1;t_now^2 t_now 1;t_next^2 t_next 1] \ [h_prev, h_now, h_next]
    if b^2-4*a*c > 0
        # We will return the smallest root. This will correspond to the '=' solution
        proposed_root = (-b-sqrt(b^2-4*a*c))/(2*a)
        # Another good point to text would be the critical point. If there is a crossing, this point is predicted to have crossed the most.
        # We will record this number. As a crossing should be the most extreme, between t_now and this critical point, the linear test should work.
        critical_point = -b/(2*a)
        if t_now < proposed_root < t_next
            return true, proposed_root, critical_point
        end
    end
end

#Locator Dispatches
#Isolates the root finding mathematics inside each one. This gets rid of global helpers so when we add new locators its really easy

#IF a locator is called while using an LMM, we temporaily swap to a RK method . 
#Since an event implies impending jump (which wipes the LMM history)
# using an RK method to pinpoint the impact should be sound.
#We will add the stepper to the arg list in each event locator even if not used to allow us to do this. (Should be better for the user)
function locate_event(locator, sys, solver::LMM, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::AbstractODESolver = ModifiedTrap())
    # Route through to the RK version, passing the stepper forward
    return locate_event(locator, sys, stepper, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper)
end

#Bisection Method (Iterative)
######
### WC: This requires that the LMMs can have variable step size as τ_m is going to vary and not be equal to Δt.
### Either fix the LMMs or only allow RK methods to work here.
######

# By restricting these to ::RK, we ensure that LMMs can only reach 
# these methods if they have been successfully intercepted and "converted" 
# to an RK solver.
function locate_event(::BisectionLocator, sys, solver::RK, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::AbstractODESolver = ModifiedTrap())
    τ_l, τ_r = 0.0, Δt
    h_l = h_now
    
    for _ in 1:100 # max_iter
        if (τ_r - τ_l) < tol break end

        #test midpoint
        τ_m = (τ_l + τ_r) / 2.0

        x_m, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, τ_m, tol, sol)
        h_m = guard(sys, x_m)
    
        if signbit(h_l) != signbit(h_m)
            τ_r = τ_m
        else
            τ_l = τ_m
            h_l = h_m
        end
    end
    
    t_star = tₖ + τ_l
    x_star, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, τ_l, tol, sol)
    return t_star, x_star
end
#Linear Interpolation
######
### WC: Why is there no for loop here? This should generate iterations on approximations to the crossing.
######
function locate_event(::LinearLocator, sys, solver::RK, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::AbstractODESolver = ModifiedTrap())
    x_predict, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, Δt, tol, sol)
    h_next = guard(sys, x_predict) 
    θ = -h_now / (h_next - h_now)
    t_star = tₖ + θ * Δt
    x_star, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, θ * Δt, tol, sol)
    return t_star, x_star
end

######
### WC: I am not convinced that quadratic needs to be implemented here. It is much more important to be utilized in 'crossed_guard'
### Newton's method may actually be a fun mathod to implement, but I suspect that linear will suffice.
######
function locate_event(::QuadraticLocator, sys, solver::RK, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::AbstractODESolver = ModifiedTrap())
    #Get three points 
    h₀ = h_now

    #middle point. We just take a half step instead of going one before the start point. I think itll be more stable. 
    x₁, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, Δt / 2.0, tol, sol)
    h₁ = guard(sys, x₁)

    #endpoint
    x₂, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, Δt, tol, sol)
    h₂ = guard(sys, x₂)

    #compute parabola coeffs
    c = h₀
    b = (-3.0 * h₀ + 4.0 * h₁ - h₂) / Δt
    a = 2.0 * (h₀ - 2.0 * h₁ + h₂) / (Δt ^ 2)

    #root finding with Fallbacks
    if abs(a) < tol
        #curve is basically zero, fallback to linear interp
        θ = -h₀ / (h₂ - h₀)
        τ_star = θ * Δt
    else 
        discriminant = b^2 - 4.0 * a * c

        if discriminant < 0
            #fallback to linear interp
            θ = -h₀ / (h₂ - h₀)
            τ_star = θ * Δt
        else
            #stable quad root finding
            q = -.5 * (b + sign(b) * sqrt(discriminant))
            root_1 = q / a
            root_2 = c / q

            valid_1 = 0.0 <= root_1 < Δt
            valid_2 = 0.0 <= root_2 < Δt

            #select right root
            if valid_1 && valid_2
                #if parabola crosses twice we pick first one
                τ_star = min(root_1, root_2)
            elseif valid_1
                τ_star = root_1
            elseif valid_2
                τ_star = root_2
            else
                #roots drifted to narnia? 
                θ = -h₀ / (h₂ - h₀)
                τ_star = θ * Δt
            end
        end
    end
    t_star = tₖ + τ_star
    x_star, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, τ_star, tol, sol)

    return t_star, x_star
end

