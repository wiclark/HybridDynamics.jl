# A collection of miscellaneous ODE integrators
# For all that follows:
#  1) f::Function is the vector field
#  2) z::Vector is the current state
#  3) h::Float is the step size (if fixed)
#  4) t::Float is the current time

#Runge Kutta / Single Step Family Solvers
abstract type RK <: AbstractODESolver end

struct ForwardEuler <: RK end
struct ModifiedTrap <: RK end
struct ModifiedMidpoint <: RK end
struct RichardsonExtrapolation <: RK end

#Adaptive Solvers
struct RK45 <: RK end
struct RK23 <: RK end
 
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
function rk_23_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat, tf::AbstractFloat, sys, tol; adaptive=true)
    # As this is an adaptive step solver, h is the step size from the previous step
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    h = minimum([h, tf-t])

    #Intialization to make sure they exist
    z2 = z; z3 = z

    # Loop through to find an acceptable step
    while true
        h_now = guard(sys, z)
        # Compute the two predictions and their difference
        k1 = f(z, t)

        z2 = z + h*k1
        k2 = f(z2, t+h)
        h2    = guard(sys, z2)

        z3 = z + h/4*(k1+k2)
        k3 = f(z+h/4*(k1+k2), t+h/2)
        h3    = guard(sys, z3)

        z1_3 = z + h*(1/6*k1+1/6*k2+2/3*k3)
        z1_2 = z + h*(1/2*k1+1/2*k2)

        #used to find true midpoint for quad guard check. 
        #This prevents solver from spinning up a whole new adaptive loop so x_mid is truly in the mid which is what we want for quad matrix. 
        if !adaptive
            return z1_2, h, h 
        end

        LTE = norm(z1_2 - z1_3)
        # Reject or accept?
        h_new = updated_step(LTE, tol, h, 3)
        
        h_end = guard(sys, z1_2)

        #did any intermediate stage cross guard?
        stage_crossed = (h_now * h2 < 0) || (h_now * h3 < 0)
        #Did final state completely miss crossing?
        end_missed = (h_now * h_end > 0)

        #if we crossed inside step but missed at end, force rejection
        if stage_crossed && end_missed
            h = h / 2.0 #force smaller step
            continue
        end

        if LTE < tol
            return z1_2, h, h_new
        else
            # We reject and repeat the loop with an updated step
            h = h_new
        end
        if h < 1e-12
            @warn "Step size has decreased below 1e-12"
            return z1_2, h, h_new #force break to avoid looping forever
        end
    end
end

# Runge-Kutta 45
function rk_45_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat, tf::AbstractFloat, sys, tol; adaptive=true)
    # As this is an adaptive step solver, h is the step size from the pervious step
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    h = minimum([h, tf-t])

    #Initialization to make them exist
    z2 = z; z3 = z; z4 = z; z5 = z; z6 = z
    # Loop through to find an acceptable step
    while true
        h_now = guard(sys, z)
        # Compute the two predictions and their difference
        k1 = f(z, t)

        z2 = z + h*1/5*k1
        k2 = f(z+h*1/5*k1, t+h*1/5)
        h2 = guard(sys, z2)

        z3 = z + h*(3/40*k1 + 9/40*k2)
        k3 = f(z+h*(3/40*k1+9/40*k2), t+h*3/10)
        h3 = guard(sys, z3)

        z4 = z + h*(44/45*k1 - 56/15*k2 + 32/9*k3)
        k4 = f(z+h*(44/45*k1-56/15*k2+32/9*k3), t+h*4/5)
        h4 = guard(sys, z4)

        z5 = z + h*(19372/6561*k1 - 25360/2187*k2 + 64448/6561*k3 - 212/729*k4)
        k5 = f(z+h*(19372/6561*k1-25360/2187*k2+64448/6561*k3-212/729*k4),t+h*8/9)
        h5 = guard(sys, z5)

        z6 = z + h*(9017/3168*k1 - 355/33*k2 + 46732/5247*k3 + 49/176*k4 - 5105/18656*k5)
        k6 = f(z+h*(9017/3168*k1-355/33*k2+46732/5247*k3+49/176*k4-5105/18656*k5), t+h)
        h6 = guard(sys, z6)

        k7 = f(z + h*(35/384*k1 + 0*k2 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6), t+h)

        # The two updates
        z1_4 = z + h*k7
        z1_5 = z + h*(5179/57600*k1 + 0*k2 + 7571/16695*k3 + 393/640*k4 - 92097/339200*k5 + 187/2100*k6 + 1/40*k7)

        #used to find true midpoint for quad guard check. 
        #This prevents solver from spinning up a whole new adaptive loop so x_mid is truly in the mid which is what we want for quad matrix. 
        if !adaptive
            return z1_4, h, h
        end

        LTE = norm(z1_4 - z1_5)
        # Reject or accept?
        h_new = updated_step(LTE, tol, h, 5)

        h_end = guard(sys, z1_4)

        stage_crossed = (h_now * h2 < 0) || (h_now * h3 < 0) || (h_now * h4 < 0) || (h_now * h5 < 0) || (h_now * h6 < 0)
        end_missed = (h_now * h_end > 0)

        if stage_crossed && end_missed
            h = h / 2.0 #force smaller step
            continue
        end

        if LTE < tol
            return z1_4, h, h_new
        else
            h = h_new
        end
        if h < 1e-12
            @warn "Step size has decreased below 1e-12"
            return z1_4, h, h_new
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

######
### WC: At some point, it would be good to add a stiff solver to the roster.
######

#Single take_step for RK methods
#Fixed step methods 
const FixedRK = Union{ForwardEuler, ModifiedTrap, ModifiedMidpoint, RichardsonExtrapolation}
#Helper function to take the step using multiple dispatch.
compute_step(::ForwardEuler, f, x, Δt, t) = forward_euler_step(f, x, Δt, t)
compute_step(::ModifiedTrap, f, x, Δt, t) = modified_trap_step(f, x, Δt, t)
compute_step(::ModifiedMidpoint, f, x, Δt, t) = modified_midpoint_step(f, x, Δt, t)
compute_step(::RichardsonExtrapolation, f, x, Δt, t) = richardson_step(f, x, Δt, t)

######
### WC: Add an RK4 method (not adaptive) and set this as a default for the LMM initializations.
######

#Adaptive step methods
const AdaptiveRK = Union{RK23, RK45}
#Helper function to take the step via multiple dispatch
compute_step(::RK23, f, x, Δt, t, tf, sys, tol; adaptive=true) = rk_23_step(f, x, Δt, t, tf, sys, tol; adaptive=adaptive)
compute_step(::RK45, f, x, Δt, t, tf, sys, tol; adaptive=true) = rk_45_step(f, x, Δt, t, tf, sys, tol; adaptive=adaptive)

#Note sol is not used, we do this to make using the function easier. We would need an if/else statement everytime we use this function without it
function take_step(solver::FixedRK, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver=ModifiedMidpoint()) 
    sys = prob.sys
    x_predict = compute_step(solver, f, xₖ, Δt, tₖ)
    x_mid     = compute_step(solver, f, xₖ, Δt / 2.0, tₖ)

    #Evaluate Guards
    h_now  = guard(sys, xₖ)
    h_mid  = guard(sys, x_mid)
    h_next = guard(sys, x_predict)

    #Use cross guard check
    eventtrigger, dt_next, _ = crossed_guard(h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol)

    return x_predict, eventtrigger, h_now, Δt, dt_next
end

function take_step(solver::AdaptiveRK, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol,stepper::AbstractODESolver=ModifiedMidpoint())
    sys = prob.sys
    tf = prob.tspan[2] #terminal time

    #Take adaptive step (passes tf to prevent overshooting)
    #Disable adaptive loop to get the true midpoint. Without it it would reset the loop and ruin or previous adaptivity
    x_predict, dt_used, dt_next = compute_step(solver, f, xₖ, Δt, tₖ, tf, sys, tol; adaptive=false)

    #Get midpoint for quad guard check
    #For a half-step the local ceiling is just midpoint of time 
    x_mid, _, _ = compute_step(solver, f, xₖ, dt_used / 2.0, tₖ, tf, sys, tol)

    #eval guards
    h_now  = guard(sys, xₖ)
    h_mid  = guard(sys, x_mid)
    h_next = guard(sys, x_predict)

    eventtrigger, _, _ = crossed_guard(h_now, h_mid, h_next, tₖ, tₖ + dt_used / 2.0, tₖ + dt_used; tol=tol)

    return x_predict, eventtrigger, h_now, dt_used, dt_next
end

#===========================#
######
### WC: Have you had any success in learning about adaptive step size LMMs?
######
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
function take_step(solver::LMM, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver = ModifiedTrap())
    sys = prob.sys
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

        eventtrigger, dt_next, _ = crossed_guard(h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol)
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

        eventtrigger, dt_next, _ = crossed_guard(h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt; tol = tol)

        return x_predict, eventtrigger, h_now, dt_next
    end

    #Guard Eval
    h_now  = guard(sys, xₖ)
    h_mid  = guard(sys, x_mid)
    h_next = guard(sys, x_predict)

    eventtrigger, dt_next, _ = crossed_guard(h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol)

    return x_predict, eventtrigger, h_now, Δt, dt_next
end


#====================================#
#Event detection utility. 
#If a guard surface was crossed and during the ODE step. We check for a sign change between start and end of the step.
function crossed_guard(h_prev, h_now, h_next, t_prev, t_now, t_next; tol=1e-6)

    # Linear Crossing check 
    #if the sign changes between now and next a root must exist. 
    #updated to account for Logical error
    if (h_prev * h_now < 0)
        #Linear interp
        t_root = t_prev - h_prev * (t_now - t_prev) / (h_now - h_prev)
        return true, t_root, NaN
    elseif (h_now * h_next < 0)
        t_root = t_now - h_now * (t_next - t_now) / (h_next - h_now)
        return true, t_root, NaN 
    end


    # Performing the quadratic version. If the discriminant is positive, there are roots.
    # Recall, that by the IVT, the linear test guarantees a crossing. The quadratic test does not guarantee one. This triggering should be treated as a warning.

    try 
        #matrix of three points we use per WC
        A = [t_prev^2 t_prev 1; t_now^2 t_now 1; t_next^2 t_next 1]
        a, b, c = A \ [h_prev, h_now, h_next]

        #ensure a is valid (a != 0)
        if abs(a) > 1e-10
            discriminant = b^2 - 4*a*c

            #if disc is positive, there are roots
            if discriminant > 0 
                #calc roots
                sqrt_d = sqrt(discriminant)
                r1 = (-b - sqrt_d) / (2*a)
                r2 = (-b + sqrt_d) / (2*a)

                #Fixed logic error where r1 not in tspan but r2 is
                valid_roots = Float64[]
                if t_prev <= r1 <= t_next push!(valid_roots, r1) end
                if t_prev <= r2 <= t_next push!(valid_roots, r2) end

                if !isempty(valid_roots)
                    proposed_root = minimum(valid_roots) #smallest valid root
                    critical_point = -b / (2*a)
                    return true, proposed_root, critical_point
                end
            end
        end
    catch 
        #if matrix is singular or calc fails we ignore quad warning
        return false, NaN, NaN
    end
    return false, NaN, NaN
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
######
### WC: Should it be 'stepper::RK'?
######
function locate_event(::BisectionLocator, prob, solver::RK, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::AbstractODESolver = ModifiedTrap())
    sys = prob.sys
    τ_l, τ_r = 0.0, Δt
    h_l = h_now
    
    for _ in 1:100 # max_iter
        if (τ_r - τ_l) < tol break end

        #test midpoint
        τ_m = (τ_l + τ_r) / 2.0

        x_m, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_m, tol, sol, stepper)
        h_m = guard(sys, x_m)
    
        if signbit(h_l) != signbit(h_m)
            τ_r = τ_m
        else
            τ_l = τ_m
            h_l = h_m
        end
    end
    
    t_star = tₖ + τ_l
    x_star, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_l, tol, sol, stepper)
    return t_star, x_star
end
#Linear Interpolation
######
### WC: Why is there no for loop here? This should generate iterations on approximations to the crossing.
######
function locate_event(::LinearLocator, prob, solver::RK, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::AbstractODESolver = ModifiedTrap())
    sys = prob.sys
    x_predict, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, Δt, tol, sol, stepper)
    h_next = guard(sys, x_predict) 
    θ = -h_now / (h_next - h_now)
    t_star = tₖ + θ * Δt
    x_star, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, θ * Δt, tol, sol, stepper)
    return t_star, x_star
end

######
### WC: I am not convinced that quadratic needs to be implemented here. It is much more important to be utilized in 'crossed_guard'
### Newton's method may actually be a fun mathod to implement, but I suspect that linear will suffice.
######
function locate_event(::QuadraticLocator, prob, solver::RK, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::AbstractODESolver = ModifiedTrap())
    sys = prob.sys

    #Get three points 
    h₀ = h_now

    #middle point. We just take a half step instead of going one before the start point. I think itll be more stable. 
    x₁, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, Δt / 2.0, tol, sol, stepper)
    h₁ = guard(sys, x₁)

    #endpoint
    x₂, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, Δt, tol, sol, stepper)
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
    x_star, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_star, tol, sol, stepper)

    return t_star, x_star
end


######
# I don't know where to put this
# Cubic Hermite interpolant for dense output stuff

struct HermiteInterp
    t0::Float64             # times
    t1::Float64
    x0::Vector{Float64}     # states
    x1::Vector{Float64}
    f0::Vector{Float64}     # derivatives
    f1::Vector{Float64}
end
# CK: This still needs to be combined and integrated into solve dispatches

# A function that computes the third-order Hermite interpolation
# This only works for scalar x_data
function hermite_interpolation(t_data::Vector, x_data::Vector, f::Function, t::AbstractFloat)
    # The first step is to make sure that t∈t_data
    if t < t_data[1] || t > t_data[end]
        @warn "Time is out of bounds"
        return NaN
    end
    
    # Next, determine the interval t lives in
    idx_first = searchsortedfirst(t_data, t) - 1
    idx_second = idx_first + 1
    # Gather all of the useful information
    t₁, t₂ = t_data[idx_first], t_data[idx_second]
    x₁, x₂ = x_data[idx_first], x_data[idx_second]
    f₁, f₂ = f(x₁), f(x₂)
    Δt = t₂ - t₁
    # Determine the coefficients
    Aⁱ = [2/Δt^3 -2/Δt^3  1/Δt^2 1/Δt^2;
         -3/Δt^2  3/Δt^2 -2/Δt  -1/Δt;
          0       0       1      0;
          1       0       0      0]
    bⁱ = [x₁, x₂, f₁, f₂]
    α, β, γ, δ = Aⁱ * bⁱ
    ts = t - t₁
    return α*ts^3 + β*ts^2 + γ*ts + δ
end