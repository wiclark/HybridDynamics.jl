# A collection of miscellaneous ODE integrators
# For all that follows:
#  1) f::Function is the vector field
#  2) z::Vector is the current state
#  3) h::Float is the step size (if fixed)
#  4) t::Float is the current time

## Single step, fully explicit methods

######
### WC: what all is being used here? 'forward_euler' does not call 'forward_euler_step'. Fat should be trimmed.
######

# Forward Euler
function forward_euler_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat)
    return z .+ h*f(z, t)
end

######
### WC: Using t_start:dt:t_end will cause issues. Suppose, for example, that we have A = 0.0:0.1:2.12. What is A[end]?
######

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

######
### WC: These two methods have a fixed value of Δt. There is no need to require t_prev (or t_prev1 and t_prev2)
######

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

######
### WC: Your parabolic version is incorrect. You do not necessarily know that the points have a uniform Δt (e.g., RK45). 
###     See my version below:
######

function crossed_guard(h_now, h_mid, h_next; tol=1e-6)\
    if isnothing(h_now)
        return false
    end

    #Standardendpoint check returns true if sign flipped or if we start on guard and push through
    if (h_now * h_next < 0) || (abs(h_now) <= tol && h_next < -tol)
        return true
    end
    
    #Quad step. Fit P(τ) = aτ^2 + bτ + c for τ∈[0,1] using start mid and end points
    c = h_now
    b = 4 * h_mid - 3 * h_now - h_next
    a = 2 * h_next + 2 * h_now - 4 * h_mid
    
    #if a is positive, parabola is concave up
    #this should be only shape that dips below zero to come back with positive endpoints
    if a > tol
        τ_vertex = -b / (2 * a)

        #check if vertex (lowest point) is in current step int
        if 0.0 < τ_vertex < 1.0
            h_min = a * τ_vertex^2 + b * τ_vertex + c

            #if dips below guard, trigger event
            if h_min <= tol
                return true
            end
        end
    end
    return false
end


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
    a, b, c = [t_prev^2 t_pref 1;t_now^2 t_now 1;t_next^2 t_next 1] \ [h_prev, h_now, h_next]
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


#solver steps we can have. Should be easy to implement more by just adding on. 
#Look at FEuler to see the format of everything else. 

function take_step(::ForwardEuler, sys, f, xₖ, tₖ, Δt, tol)
    x_predict = forward_euler_step(f, xₖ, Δt, tₖ) #next state based on the linear/affine dynamics

    ### WC: no, don't do this. the midpoint as the average is a garbage approximation
    ###     If anything, take a half Euler step. (Better would be to take in the pervious state.)
    ###     Upon reflection, you really should take the previous state as an additional input to this function (you don't know that Δt is constant)

    #Calc midpoint for quad
    #x_mid = (xₖ .+ x_predict) ./ 2.0
    x_mid = forward_euler_step(f, xₖ, Δt/2, tₖ)

    h_now = guard(sys, xₖ) #evaluates guard function at current position to see how far we are from it
    h_mid = guard(sys, x_mid) #eval guard at midpoint between start and end point for quad check
    h_next = guard(sys, x_predict) #Evalutes guard function at predicted next position to check if we moved through it
    #eventtrigger = crossed_guard(h_now, h_mid, h_next; tol=tol) #compares the above to check if we crossed guard. 
    eventtrigger, dt_next, _ = crossed_guard(h_now, h_mid, h_next; tol=tol)

    #dt_next = Δt * 1.2 

    return x_predict, eventtrigger, h_now, dt_next
end

######
### WC: This seems too simliar to the previous version to warrant another dispatch
######
function take_step(::ModifiedTrap, sys, f, xₖ, tₖ, Δt, tol)
    x_predict = modified_trap_step(f, xₖ, Δt, tₖ) #calc next state using modifed trap method

    x_mid = modified_trap_step(f, xₖ, Δt / 2.0, tₖ)

    #eval guard condition at the start and predicted positions
    h_now = guard(sys, xₖ) 
    h_mid = guard(sys, x_mid)
    h_next = guard(sys, x_predict)

    #check for sign change (I hope to make htis more dignified later)
    eventtrigger = crossed_guard(h_now, h_mid, h_next; tol=tol)

    dt_next = Δt * 1.2

    return x_predict, eventtrigger, h_now, dt_next
end

######
### WC: I feel like there should be a smarter way to implement this case. One that actually uses the matrix exponential.
######
function take_step(::ExponentialSolver, sys, f, xₖ, tₖ, Δt, tol)
    flowmap = LinearFlow(sys.A)
    x_predict = flow(flowmap, Δt, xₖ)
    x_mid = flow(flowmap, Δt / 2.0, xₖ)
    h_now = guard(sys, xₖ)
    h_mid = guard(sys, x_mid)
    h_next = guard(sys, x_predict)
    eventtrigger = crossed_guard(h_now, h_mid, h_next; tol=tol)
    dt_next = Δt * 1.2
    return x_predict, eventtrigger, h_now, dt_next
end


######
### WC: Same comment as above. Why another dispatch?
######
#Extrapolation 
function take_step(::RichardsonExtrapolation, sys, f, xₖ, tₖ, Δt, tol)
    #calc with richardson extrapolation
    x_predict = richardson_step(f, xₖ, Δt, tₖ)

    x_mid = richardson_step(f, xₖ, Δt / 2.0, tₖ)

    #eval guard conditions 
    h_now = guard(sys, xₖ)
    h_mid = guard(sys, x_mid)
    h_next = guard(sys, x_predict)

    #check event 
    eventtrigger = crossed_guard(h_now, h_mid, h_next; tol=tol)

    dt_next = Δt * 1.2

    return x_predict, eventtrigger, h_now, dt_next
end

#LINEAR MULTISTEP METHODS STUFF
#This take_step is used to differentiate between the LMM and the other types. Passing sol is what does it. Then we use the LMM only when we have this type of arguments
function take_step(solver::AbstractODESolver, sys, f, xₖ, tₖ, Δt, tol, sol)
    return take_step(solver, sys, f, xₖ, tₖ, Δt, tol)
end

######
### WC: This actually feels quite different from the pervious dispatches. Is it because these are LMMs and you are tracking the history?
###     In that case, you should only have to have two dispatches for 'take_step'. One for Runge-Kutta methods and one for LMMs. Possibly a third/fourth for adaptive step versions
######
#AB2
function take_step(::AdamsBashforth2, sys, f, xₖ, tₖ, Δt, tol, sol)
    #check how many continuous steps have occurred since last jump
    #sol.jump_indices[end] tells us where the current trajectory started
    history_len = isempty(sol.jump_indices) ? length(sol.x) : (length(sol.x) - sol.jump_indices[end] + 1)

    if history_len < 2
        #post jump reset phase:
        #we dont have enough points yet. we use Forward Euler to guess for now
        ######
        ### WC: Don't you mean modified trapezoid?
        ######
        x_predict = modified_trap_step(f, xₖ, Δt, tₖ)
        x_mid = modified_trap_step(f, xₖ, Δt / 2.0, tₖ)
    else
        #Multistep phase: extract current and previous data points
        x_prev = sol.x[end-1]
        t_prev = sol.t[end-1]
        x_predict = adamsbashforth2_step(f, xₖ, tₖ, Δt, x_prev, t_prev)
        x_mid = modified_trap_step(f, xₖ, Δt / 2.0, tₖ)
    end

    #usual guard and event logic
    h_now = guard(sys, xₖ)
    h_mid = guard(sys, x_mid)
    h_next = guard(sys, x_predict)
    ######
    ### WC: This should be updated as above.
    ######
    eventtrigger = crossed_guard(h_now, h_mid, h_next; tol=tol)
    dt_next = Δt * 1.2
    return x_predict, eventtrigger, h_now, dt_next
end

######
### WC: This shouldn't really require an additional dispatch. You should be write one for an abstract LMM (of arbitrary step length)
### Actually, do you actually need xₖ, tₖ as inputs. Aren't these already stored in 'sys'?
### I suppose that you would have to supply a predicting algorithm
######
function take_step(::AdamsBashforth3, sys, f, xₖ, tₖ, Δt, tol, sol)
    #Determine how many cont steps we have since last jump
    history_len = isempty(sol.jump_indices) ? length(sol.x) : (length(sol.x) - sol.jump_indices[end] + 1)

    if history_len < 3
        #startup phase: less than 3 points of cont history we need. Use Richardson extra to build
        x_predict = richardson_step(f, xₖ, Δt, tₖ)
        x_mid = richardson_step(f, xₖ, Δt / 2.0, tₖ)
    else 
        #multistep phase: get history time
        x_prev1 = sol.x[end-1]
        t_prev1 = sol.t[end-1]

        x_prev2 = sol.x[end-2]
        t_prev2 = sol.t[end-2]

        x_predict = adamsbashforth3_step(f, xₖ, tₖ, Δt, x_prev1, t_prev1, x_prev2, t_prev2)

        x_mid = richardson_step(f, xₖ, Δt / 2.0, tₖ)
    end
    #guard stuff
    h_now = guard(sys,xₖ)
    h_mid = guard(sys, x_mid)
    h_next = guard(sys, x_predict)
    eventtrigger = crossed_guard(h_now, h_mid, h_next; tol=tol)
    dt_next = Δt * 1.2

    return x_predict, eventtrigger, h_now, dt_next
end

######
### WC: Observe. Also do something with take_step(::RK,...)
######
function take_step(::LMM, sys, f, xₖ, tₖ, Δt, tol, sol, prediction_stepper, k)
    # Determine how many cont steps we have since last jump
    history_len = isempty(sol.jump_indices) ? length(sol.x) : (length(sol.x) - sol.jump_indices[end] + 1)

    # Do we have a sufficiently rich history?
    if history_len < k
        x_predict = prediction_stepper(f, xₖ, Δt, tₖ)
        x_mid = prediction_stepper(f, xₖ, Δt / 2.0, tₖ)
    else
        # This is the part where you need a little creativity
        x_prev = sol.x[end-k:end]
        t_prev = sol.t[end-k:end]
        x_predict = LLM(f, xₖ, tₖ, x_prev, t_prev)
        # I still don't like this. You should be comparing to the previous evaluation, not the next half one
        x_mid = prediction_stepper(f, xₖ, Δt / 2.0, tₖ)
    end
    # Guard stuff
    h_now = guard(sys,xₖ)
    h_mid = guard(sys, x_mid)
    h_next = guard(sys, x_predict)
    eventtrigger, dt_next, _ = crossed_guard(h_now, h_mid, h_next; tol=tol)

    return x_predict, event_trigger, h_now, dt_next
end


#Locator Dispatches
#Isolates the root finding mathematics inside each one. This gets rid of global helpers so when we add new locators its really easy

#Bisection Method (Iterative)
######
### WC: This requires that the LMMs can have variable step size as τ_m is going to vary and not be equal to Δt.
### Either fix the LMMs or only allow RK methods to work here.
######
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
######
### WC: Why is there no for loop here? This should generate iterations on approximations to the crossing.
######
function locate_event(::LinearLocator, sys, solver, f, xₖ, tₖ, Δt, h_now, tol, sol)
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
function locate_event(::QuadraticLocator, sys, solver, f, xₖ, tₖ, Δt, h_now, tol, sol)
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

