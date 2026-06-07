# A collection of miscellaneous ODE integrators
# For all that follows:
#  1) f::Function is the vector field
#  2) z::Vector is the current state
#  3) h::Float is the step size (if fixed)
#  4) t::Float is the current time

## Single step, fully explicit methods

# Forward Euler
function forward_euler_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat)
    return z .+ h*f(z, t)
end

# NOT SINGLE STEP: Forward Euler method. We solve an ODE defined by $f(u,t)$ starting at u0 and over tspan with step size dt. 
function forward_euler(f::Function,u0,tspan::Tuple{Float64,Float64}; dt::Float64 = 0.01)
    t_start, t_end = tspan

    #Create a time vector:
    t = collect(t_start:dt:t_end) #gets the range of values for time and puts thm into an array (vector). So we can index. 
    num_steps = length(t)

    #intialize the solution array to match the type of initial condition 
    u = Vector{typeof(u0)}(undef, num_steps) #typeof so we can keep things straight. undef is supposedly faster than zeros() but I didnt fact check that. Makes sense though as we skip some values
    u[1] = u0

    for i in 1:(num_steps-1)
        #Eulers update: u_next = u_now + dt * slope
        u[i+1] = u[i] + dt* f(u[i], t[i])
    end
    return t,u
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

#Extrapolation Step
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

function adamsbashforth2_step(f::Function, xₖ::Vector, tₖ::AbstractFloat, Δt::AbstractFloat, x_prev::Vector, t_prev::AbstractFloat)
    #calc previous time step size
    dt_previous = tₖ - t_prev

    #eval current and past derivatives
    fₖ = f(xₖ, tₖ)
    f_prev = f(x_prev, t_prev)

    #AB2 formula 
    α = Δt / dt_previous
    return xₖ .+ Δt .* ((1.0 + .5 * α) .* fₖ .- (0.5 * α) .* f_prev)
end

function adamsbashforth3_step(f::Function, xₖ::Vector, tₖ::AbstractFloat, Δt::AbstractFloat, x_prev1::Vector, t_prev1::AbstractFloat, x_prev2::Vector, t_prev2::AbstractFloat)
    fₖ = f(xₖ, tₖ)
    f_prev1 = f(x_prev1, t_prev1)
    f_prev2 = f(x_prev2, t_prev2)

    return xₖ .+ Δt .* ( (23/12) .* fₖ .- (16/12) .* f_prev1 .+ (5/12) .* f_prev2 )
end
    


#Event detection utility. 
#If a guard surface was crossed and during the ODE step. We check for a sign change between start and end of the step. 
#Made as a function so we can use it to all solvers with the same logic. We will need another idea for event detection when signs dont change but I havent gotten that far
function crossed_guard(h_now, h_next; tol=1e-6)
    # Handles the case when no guard is provided
    if isnothing(h_now)
        return false
    else
    #returns true if sign flipped, or if we start on the guard and push through
    return (h_now * h_next < 0) || (abs(h_now) <= tol && h_next < -tol)
    end
end

#solver steps we can have. Should be easy to implement more by just adding on. 
#only for FEuler
function take_step(::ForwardEuler, sys, f, xₖ, tₖ, Δt, tol)
    x_predict = forward_euler_step(f, xₖ, Δt, tₖ) #next state based on the linear/affine dynamics

    h_now = guard(sys, xₖ) #evaluates guard function at current position to see how far we are from it
    h_next = guard(sys, x_predict) #Evalutes guard function at predicted next position to check if we moved through it
    eventtrigger = crossed_guard(h_now, h_next; tol=tol) #compares the above to check if we crossed guard. 

    dt_next = Δt * 1.2 

    return x_predict, eventtrigger, h_now, dt_next
end

function take_step(::ModifiedTrap, sys, f, xₖ, tₖ, Δt, tol)
    x_predict = modified_trap_step(f, xₖ, Δt, tₖ) #calc next state using modifed trap method

    #eval guard condition at the start and predicted positions
    h_now = guard(sys, xₖ) 
    h_next = guard(sys, x_predict)

    #check for sign change (I hope to make htis more dignified later)
    eventtrigger = crossed_guard(h_now, h_next; tol=tol)

    dt_next = Δt * 1.2

    return x_predict, eventtrigger, h_now, dt_next
end

function take_step(::ExponentialSolver, sys, f, xₖ, tₖ, Δt, tol)
    flowmap = LinearFlow(sys.A)
    x_predict = flow(flowmap, Δt, xₖ)
    h_now = guard(sys, xₖ)
    h_next = guard(sys, x_predict)
    eventtrigger = crossed_guard(h_now, h_next; tol=tol)
    dt_next = Δt * 1.2
    return x_predict, eventtrigger, h_now, dt_next
end

#Extrapolation 
function take_step(::RichardsonExtrapolation, sys, f, xₖ, tₖ, Δt, tol)
    #calc with richardson extrapolation
    x_predict = richardson_step(f, xₖ, Δt, tₖ)

    #eval guard conditions 
    h_now = guard(sys, xₖ)
    h_next = guard(sys, x_predict)

    #check event 
    eventtrigger = crossed_guard(h_now, h_next; tol=tol)

    dt_next = Δt * 1.2

    return x_predict, eventtrigger, h_now, dt_next
end

#LINEAR MULTISTEP METHODS STUFF
#This take_step is used to differentiate between the LMM and the other types. Passing sol is what does it. Then we use the LMM only when we have this type of arguments
function take_step(solver::AbstractODESolver, sys, f, xₖ, tₖ, Δt, tol, sol)
    return take_step(solver, sys, f, xₖ, tₖ, Δt, tol)
end

#AB2
function take_step(::AdamsBashforth2, sys, f, xₖ, tₖ, Δt, tol, sol)
    #check how many continuous steps have occurred since last jump
    #sol.jump_indices[end] tells us where the current trajectory started
    history_len = isempty(sol.jump_indices) ? length(sol.x) : (length(sol.x) - sol.jump_indices[end] + 1)

    if history_len < 2
        #post jump reset phase:
        #we dont have enough points yet. we use Forward Euler to guess for now
        x_predict = modified_trap_step(f, xₖ, Δt, tₖ)
    else
        #Multistep phase: extract current and previous data points
        x_prev = sol.x[end-1]
        t_prev = sol.t[end-1]
        x_predict = adamsbashforth2_step(f, xₖ, tₖ, Δt, x_prev, t_prev)
    end

    #usual guard and event logic
    h_now = guard(sys, xₖ)
    h_next = guard(sys, x_predict)
    eventtrigger = crossed_guard(h_now, h_next; tol=tol)
    dt_next = Δt * 1.2
    return x_predict, eventtrigger, h_now, dt_next
end

function take_step(::AdamsBashforth3, sys, f, xₖ, tₖ, Δt, tol, sol)
    #Determine how many cont steps we have since last jump
    history_len = isempty(sol.jump_indices) ? length(sol.x) : (length(sol.x) - sol.jump_indices[end] + 1)

    if history_len < 3
        #startup phase: less than 3 points of cont history we need. Use Richardson extra to build
        x_predict = richardson_step(f, xₖ, Δt, tₖ)
    else 
        #multistep phase: get history time
        x_prev1 = sol.x[end-1]
        t_prev1 = sol.t[end-1]

        x_prev2 = sol.x[end-2]
        t_prev2 = sol.t[end-2]

        x_predict = adamsbashforth3_step(f, xₖ, tₖ, Δt, x_prev1, t_prev1, x_prev2, t_prev2)
    end
    #guard stuff
    h_now = guard(sys,xₖ)
    h_next = guard(sys, x_predict)
    eventtrigger = crossed_guard(h_now, h_next; tol=tol)
    dt_next = Δt * 1.2

    return x_predict, eventtrigger, h_now, dt_next
end


#Locator Dispatches
#Isolates the root finding mathematics inside each one. This gets rid of global helpers so when we add new locators its really easy

#Bisection Method (Iterative)
function locate_event(::BisectionLocator, sys, solver, f, xₖ, tₖ, Δt, h_now, tol, sol)
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
function locate_event(::LinearLocator, sys, solver, f, xₖ, tₖ, Δt, h_now, tol, sol)
    x_predict, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, Δt, tol, sol)
    h_next = guard(sys, x_predict) 
    θ = -h_now / (h_next - h_now)
    t_star = tₖ + θ * Δt
    x_star, _, _, _ = take_step(solver, sys, f, xₖ, tₖ, θ * Δt, tol, sol)
    return t_star, x_star
end

function locate_event(::QuadraticLocator, sys, solver, f, xₖ, tₖ, Δt, h_now, tol, sol)
    #Get three points 
    h₀ = h_now

    #middle point 
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

