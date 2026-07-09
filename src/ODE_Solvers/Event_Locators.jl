#====================================#
#Event detection utility. 
#Assigns a default crossing direction to a specific system to reduce issues (this may be made better for the front end later but for now this works)
default_guard_direction(sys::MechanicalSystem) = sys.direction
default_guard_direction(sys::NonholonomicSystem) = sys.direction
#Gen system lack specific constraints so we monitor crossings in both directions
default_guard_direction(sys::GeneralSystem) = sys.direction
default_guard_direction(sys::LinearSystem) = sys.direction
default_guard_direction(sys::AffineSystem) = sys.direction
default_guard_direction(sys::StochasticSystem) = sys.direction
default_guard_direction(sys) = 0

#Wrapper function to interface between the system state and core logic
function crossed_guard(sys, h_prev, h_now, h_next, t_prev, t_now, t_next; 
                       tol=1e-6, direction=default_guard_direction(sys))   
    # Calls evaluator with the directionality
    return evaluate_crossing(h_prev, h_now, h_next, t_prev, t_now, t_next, direction; tol=tol)
end

#Core engine: determines if/when the guard function 'h' changes sign. 
function evaluate_crossing(h_prev, h_now, h_next, t_prev, t_now, t_next, direction::Int; tol=1e-6)
    #Helper to validate a linear sign change based on the required direction
    valid_linear(h1, h2) = 
        (direction == 0 && h1 * h2 < 0) ||
        (direction == -1 && h1 > 0 && h2 < 0) ||
        (direction == 1 && h1 < 0 && h2 > 0)

    #First check: Simple linear crossing detection. Between previous and current
    if valid_linear(h_prev, h_now)
        #Linear interpolation to find the root
        t_root = t_prev - h_prev * (t_now - t_prev) / (h_now - h_prev)
        # println([h_prev, h_now, h_next])
        return true, t_root, NaN
    #Second check: Linear beween current and next
    elseif valid_linear(h_now, h_next)
        t_root = t_now - h_now * (t_next - t_now) / (h_next - h_now)
        return true, t_root, NaN 
    end

    #Quad version
    try
        #Set up linear system to solve for parabola coeffs  
        A = [t_prev^2 t_prev 1; t_now^2 t_now 1; t_next^2 t_next 1]
        a, b, c  = A \ [h_prev, h_now, h_next]

        #Ensure it is actually a parabola then check disc
        if abs(a) > tol
            discriminant = b^2 - 4*a*c
            if discriminant > 0
                sqrt_d = sqrt(discriminant)
                r1 = (-b - sqrt_d) / (2*a)
                r2 = (-b + sqrt_d) / (2*a)

                valid_roots = Float64[]
                eps_t = 1e-9 #Buffer to ensure root isnt some weird artifact (it does seem necessary)

                #Derivative of parabola eval the slope of root. 
                h_prime(t) = 2*a*t + b

                #Helper: Does quad curve cross in the correct direction?
                valid_quad(r) = 
                    (direction == 0) ||
                    (direction == -1 && h_prime(r) < 0) ||
                    (direction == 1 && h_prime(r) > 0)

                #Check bounds and direction
                if (t_prev + eps_t) <= r1 <= t_next && valid_quad(r1) push!(valid_roots, r1) end 
                if (t_prev + eps_t) <= r2 <= t_next && valid_quad(r2) push!(valid_roots, r2) end

                #If root exists, return earliest and the parabolas critical point
                if !isempty(valid_roots)
                    proposed_root = minimum(valid_roots)
                    critical_point = -b / (2*a)
                    return true, proposed_root, critical_point
                end
            end
        end
    catch
        #If linear system is singular or math fails, return failure to cross. 
        return false, NaN, NaN
    end
    return false, NaN, NaN
end

#Locator Dispatches
#Isolates the root finding mathematics inside each one. This gets rid of global helpers so when we add new locators its really easy

#-----------------------
#LOCATOR TAGS
#These are similar to the solver tags. But these define how the solver finds the exact crossing time once an event is detected. 

#Tag to use Linear Interpolation. Very fast but can be innacurate for higher order methods. 
struct LinearLocator <: AbstractEventLocator end

#Tag to use a bisection method serach. Can be very accurate but also very slow with complex systems
struct BisectionLocator <: AbstractEventLocator end

#Tag for quadratic event locator
struct QuadraticLocator <: AbstractEventLocator end

#Tag to use Newtons method for event locators
struct NewtonLocator <: AbstractEventLocator end

#Bisection Method (Iterative)

function locate_event(::BisectionLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys
    τ_l, τ_r = 0.0, Δt
    h_l = h_now
    
    for _ in 1:100 # max_iter
        if (τ_r - τ_l) < tol break end

        #test midpoint
        τ_m = (τ_l + τ_r) / 2.0

        x_m, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_m, tol, sol, stepper; check=false)
        h_m = guard(sys, x_m)
    
        if signbit(h_l) != signbit(h_m)
            τ_r = τ_m
        else
            τ_l = τ_m
            h_l = h_m
        end
    end
    
    t_star = tₖ + τ_l
    x_star, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_l, tol, sol, stepper; check=false)
    return t_star, x_star
end
#Linear Interpolation

function locate_event(::LinearLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    # Extract System
    sys = prob.sys
    # Extract left boundary to 0 and right to Δt
    τ_l, τ_r = 0.0, Δt
    # Set guard function value at the left boundary to the current value.
    h_l = h_now

    # Check if state already satisfies the event condition, if so we return the state and time.
    if abs(h_now) < tol 
        return tₖ, xₖ
    end

    # Get right side of boundary
    x_r, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_r, tol, sol, stepper; check=false)
    # Eval guard at this new right side state
    h_r = guard(sys, x_r)
    
    # If we ever get a step not bracketing a root we exit to avoid iterating garbage. 
    if signbit(h_l) == signbit(h_r)
        # Return the time and state at the end of the full step since no event took place
        return tₖ + Δt, x_r
    end

    # Initialize our best guess for the event time offset and state to the right boundary
    τ_star = Δt
    x_star = x_r

    for _ in 1:100
        # Linear Interpolation
        # Calc time offset where event occurs
        τ_m = τ_r - h_r * (τ_r - τ_l) / (h_r - h_l)

        # Step ODE solver forward by time offset τ_m
        x_m, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_m, tol, sol, stepper; check=false)
        # Eval guard at this new state
        h_m = guard(sys, x_m)

        # Check if this interpolated state satifies the event tolerance 
        if abs(h_m) < tol
            # If so we save current as our final answer and break out of the loop. 
            τ_star = τ_m
            x_star = x_m
            break
        end

        # keep root bracketed
        # If the sign at the left boudnary is diff form the sign at the new point, root is in left half 
        if signbit(h_l) != signbit(h_m)
            # Update right boundary time to the new time
            τ_r = τ_m
            # Update the right boundary guard value
            h_r = h_m
        else 
            # Otherwise root is in the right half, so update left boundary time
            τ_l = τ_m
            h_l = h_m
        end
    end
    # Calc absolute time of the event by adding base time to the found offset
    t_star = tₖ + τ_star
    # return abs time and its corresponding state. 
    return t_star, x_star
end

function locate_event(::QuadraticLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys

    # Get three points 
    # Start of interval guard value. 
    h₀ = h_now

    # middle point. We just take a half step instead of going one before the start point. I think itll be more stable. 
    x₁, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, Δt / 2.0, tol, sol, stepper; check=false)
    h₁ = guard(sys, x₁)

    # endpoint
    x₂, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, Δt, tol, sol, stepper; check=false)
    h₂ = guard(sys, x₂)

    # Must actually bracket a root
    # Check if the start and end points have same sign (meaning they likely dont have a root between them)
    if signbit(h₀) == signbit(h₂)
        # Compute linear interp guess as failsafe
        θ = -h₀ / (h₂ - h₀)
        # Ensure guessed fraction doesnt go outside the time interval 
        τ_star = clamp(θ * Δt, 0.0, Δt)

        # Step solver to time we got above. 
        x_star, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_star, tol, sol, stepper; check=false)

        # Return time and state bypassing the quadratic logic.
        return tₖ + τ_star, x_star
    end

    # compute parabola coeffs
    # y int is the starting value
    c = h₀
    # Compute Linear coeffs based on 3 points
    b = (-3.0 * h₀ + 4.0 * h₁ - h₂) / Δt
    # Compute quad coeffs based on 3 points
    a = 2.0 * (h₀ - 2.0 * h₁ + h₂) / (Δt ^ 2)

    # Check if quadratic term if effectively zero. If so we fall back to linear interpolation
    if abs(a) < tol
        θ = -h₀ / (h₂ - h₀)
        τ_star = clamp(θ * Δt, 0.0, Δt)
    else
        # Calc discriminant of quad formula
        disc = b^2 - 4.0*a*c
        # If disc is negative, parabola doesnt cross the axis with no real roots.
        if disc < 0
            # Fallback to linear interpolation if no real roots exist.
            θ = -h₀ / (h₂ - h₀)
            τ_star = clamp(θ * Δt, 0.0, Δt)
        else
            # Use quadratic formula (this is the stable version)
            q = -0.5 * (b + sign(b)*sqrt(disc))

            # First potential root
            root1 = q/a
            # Second potential root
            root2 = c/q

            # Verify if first root is a finite number and lies within the time interval
            valid1 = isfinite(root1) && (0.0 <= root1 <= Δt)
            # See above but for second root
            valid2 = isfinite(root2) && (0.0 <= root2 <= Δt)

            # If both roots are valid inside the interval...
            if valid1 && valid2
                # Choose earliest event time as the root. 
                τ_star = min(root1, root2)
            # If only the first root is valid....
            elseif valid1
                # Assign this first root as the root
                τ_star = root1
            # If only second is valid...
            elseif valid2
                # Assign this second root as the root
                τ_star = root2
            # If neither root is valid...
            else 
                #Fall back to linear interp. 
                θ = -h₀ / (h₂ - h₀)
                τ_star = clamp(θ * Δt, 0.0, Δt)
            end
        end
    end

    # Verify the quadratic prediction. We dont just trust it with our heart of hearts. 
    x_test, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_star, tol, sol, stepper; check=false)
    h_test = guard(sys, x_test)

    # If prediction is worse than midpoint we abandon it and linear interp. Just a failsafe.
    if abs(h_test) > abs(h₁)
        θ = -h₀ / (h₂ - h₀)
        τ_star = clamp(θ * Δt, 0.0, Δt)

        x_test, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_star, tol, sol, stepper; check=false)
        h_test = guard(sys, x_test)

    end

    # If the prediction STILL doesnt meet the tolerance requirements, attempt one more linear step.
    if abs(h_test) > tol
        # Reset left boundary to start of step
        τ_l = 0.0
        # Reset right boundary to end of the step
        τ_r = Δt
        # Reset left boundary guard value
        h_l = h₀
        # Reset right boundary guard value
        h_r = h₂
        
        # Check if the root is between the start and the test point 
        if signbit(h_l) != signbit(h_test)
            # Shrink right boundary to test point
            τ_r = τ_star
            h_r = h_test
        else
            # Otherwise shrink left boundary to test point
            τ_l = τ_star
            h_l = h_test
        end
        # Perform one final linear interp on new bounds 
        τ_star = τ_r - h_r*(τ_r - τ_l)/(h_r - h_l)
        # ensure final pred is within clamped bounds
        τ_star = clamp(τ_star, 0.0, Δt)
        # Step solver to this final refined time
        x_test, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_star, tol, sol, stepper; check=false)
    end
    # Calc the absolute time of the verified event
    t_star = tₖ + τ_star
    # Return abs time and its corresponding state. 
    return t_star, x_test
end

function locate_event(::NewtonLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys

    τ_prev = 0.0
    h_prev = h_now

    τ_curr = Δt
    x_curr, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_curr, tol, sol, stepper; check=false)
    h_curr = guard(sys, x_curr)

    for _ in 1:100
        if abs(h_curr) < tol || abs(τ_curr - τ_prev) < tol
            break
        end

        #Estimate derivative 
        dh_dτ = (h_curr - h_prev) / (τ_curr - τ_prev)

        if abs(dh_dτ) < tol
            @warn "Derivative less than tolerance: Newton method unstable. Terminating."
            break
        end

        #Newton step
        τ_next = τ_curr - h_curr / dh_dτ

        #Bound next step to prevent over shooting 
        τ_next = clamp(τ_next, 0.0, Δt)

        #update for next iteration
        τ_prev = τ_curr
        h_prev = h_curr

        τ_curr = τ_next
        x_curr, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_curr, tol, sol, stepper; check=false)
        h_curr = guard(sys, x_curr)
    end
    t_star = tₖ + τ_curr
    return t_star, x_curr
end
