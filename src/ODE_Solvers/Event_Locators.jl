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
#HERMITE VERSION
function crossed_guard(sys, h_now, h_next, hp_now, hp_next, t_now, t_next; 
                       tol=1e-6, direction=default_guard_direction(sys))   
    # new call for the hermite locator
    return evaluate_crossing(h_now, h_next, hp_now, hp_next, t_now, t_next, direction; tol=tol)
end

#OLD VERSION
#=
function crossed_guard(sys, h_prev, h_now, h_next, t_prev, t_now, t_next;
                       tol=1e-6, direction=default_guard_direction(sys))  
    # Calls evaluator with the directionality
    return evaluate_crossing(h_prev, h_now, h_next, t_prev, t_now, t_next, direction; tol=tol)

end 
=#

#Calcs dh/dt = ∇h(x) * dx
function guard_derivatives(sys, x, dx, h_val; ε=1e-7)
    return (guard(sys, x .+ ε .* dx) - h_val) / ε
end

#HERMITE VERSION
function evaluate_crossing(h_now, h_next, hp_now, hp_next, t_now, t_next, direction::Int; tol=1e-6)
    Δt = t_next - t_now
    #maybe useful?
    if Δt <= 0
        return false, NaN, NaN
    end
    # The first check is simply the linear crossing - this guarantees a crossing by the MVT
    valid_linear(h1, h2) = 
        (direction == 0 && h1 * h2 < 0) ||
        (direction == -1 && h1 > 0 && h2 < 0) ||
        (direction == 1 && h1 < 0 && h2 > 0)

    #I rewrote the formulas here but they are the same. Just happened as I was bugfixing - DS
    α = (3 * (h_next - h_now) - (2 * hp_now + hp_next) * Δt) / (Δt^2)
    β = (2 * (h_now - h_next) + (hp_now + hp_next) * Δt) / (Δt^3)

    # The original code used linear secant lines to approximate the root.
    # Because we already have the Hermite coefficients, we can solve the cubic polynomial - DS
    function find_cubic_root(τ_a, τ_b)
        p(τ) = h_now + hp_now * τ + α * τ^2 + β * τ^3
        dp(τ) = hp_now + 2 * α * τ + 3 * β * τ^2
        
        h_a, h_b = p(τ_a), p(τ_b)
        τ = clamp(τ_a + (τ_b - τ_a) * (-h_a / (h_b - h_a)), τ_a, τ_b) # Initial secant guess
        for _ in 1:6
            der = dp(τ)
            abs(der) < 1e-14 && break
            τ = clamp(τ - p(τ) / der, τ_a, τ_b)
        end
        return τ
    end

    # The next function attempts to find an event based off of Hermite interpolation
    function valid_hermite()
        # If the system trajectory is purely quadratic (constant acceleration), β = 0.
        # The original code calculated (6*β) in the denominator, which produced NaN and missed events. - DS
        if abs(β) < 1e-12
            if abs(α) > 1e-12
                R1 = -hp_now / (2 * α) # Quadratic vertex formula
                if 0 < R1 < Δt
                    y1 = h_now + hp_now * R1 + α * R1^2
                    if valid_linear(h_now, y1)
                        return true, find_cubic_root(0.0, R1)
                    end
                end
            end
            return false, NaN
        end

        D = 4 * α^2 - 12 * β * hp_now
        if D < 0 # Negative discriminant <- No real roots of the derivative
            return false, NaN
        else
            # Is there a (or two) turning point(s) within the interval?
            r1 = (-2 * α + sqrt(D)) / (6 * β)
            r2 = (-2 * α - sqrt(D)) / (6 * β)
            R1 = minimum([r1, r2])
            R2 = maximum([r1, r2])
            if 0 < R1 < Δt
                # Find the value at the first turning point
                y1 = h_now + hp_now * R1 + α * R1^2 + β * R1^3
                if valid_linear(h_now, y1)
                    return true, find_cubic_root(0.0, R1)
                end
            end
            if 0 < R2 < Δt
                # Find the value at the second turning point
                y2 = h_now + hp_now * R2 + α * R2^2 + β * R2^3
                if valid_linear(h_now, y2)
                    start_τ = (0 < R1 < Δt) ? R1 : 0.0
                    return true, find_cubic_root(start_τ, R2)
                end
            end 
            return false, NaN
        end 
    end
    # Do we suspect a crossing?
    if valid_linear(h_now, h_next)
        #I replaced the linear secant root finder with the cubic one - DS
        τ_root = find_cubic_root(0.0, Δt)
        return true, t_now + τ_root, NaN
    else
        is_valid, offset = valid_hermite()
        if is_valid
            return true, t_now + offset, NaN
        else
            return false, NaN, NaN
        end
    end
end

#OLD VERSION MAYBE GET RID OF ONCE GOOD ON HERMITE? 
#=
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
        if abs(a) > eps(Float64)
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
=#
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

#Tage to use Hermite locator
struct HermiteLocator <: AbstractEventLocator end

#Linear Interpolation

function locate_event(::LinearLocator, prob, solver::AbstractODESolver, f, Df, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
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
    x_r, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_r, tol, sol, stepper; check=false)
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
        x_m, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_m, tol, sol, stepper; check=false)
        # Eval guard at this new state
        h_m = guard(sys, x_m)

        # Check if this interpolated state satifies the event tolerance 
        if abs(h_m) < tol
            # If so we save current as our final answer and break out of the loop. 
            τ_star = τ_m
            x_star = x_m
            break
        end
        # Safeguard?
        τ_star = τ_m
        x_star = x_m

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

#Hermite Locator
function locate_event(::HermiteLocator, prob, solver::AbstractODESolver, f, Df, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys
    
    # Boundary setup
    τ_l, τ_r = 0.0, Δt
    x_l = xₖ
    h_l = h_now
    dx_l = f(x_l, tₖ)
    hp_l = guard_derivatives(sys, x_l, dx_l, h_l)

    if abs(h_now) < tol
        return tₖ, xₖ
    end

    x_r, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_r, tol, sol, stepper; check=false)
    h_r = guard(sys, x_r)
    dx_r = f(x_r, tₖ + Δt)
    hp_r = guard_derivatives(sys, x_r, dx_r, h_r)

    # Return step if no root is bracketed
    if signbit(h_l) == signbit(h_r)
        return tₖ + Δt, x_r
    end

    τ_m = 0.5 * Δt
    x_m = x_r

    for _ in 1:100
        dτ = τ_r - τ_l
        if dτ < 1e-12
            break
        end

        # Hermite Coeffs rewritten
        α = (3 * (h_r - h_l) - (2 * hp_l + hp_r) * dτ) / (dτ^2)
        β = (2 * (h_l - h_r) + (hp_l + hp_r) * dτ) / (dτ^3)

        p(τ) = h_l + hp_l * τ + α * τ^2 + β * τ^3
        dp(τ) = hp_l + 2 * α * τ + 3 * β * τ^2

        # Initial secant guess
        τ_step = clamp(dτ * (-h_l / (h_r - h_l)), 0.05 * dτ, 0.95 * dτ)

        # Newton refinement
        for _ in 1:5
            val = p(τ_step)
            der = dp(τ_step)
            abs(der) < 1e-14 && break
            τ_step = clamp(τ_step - val / der, 0.0, dτ)
        end

        τ_m = τ_l + τ_step

        # Step ODE solver to intermediate candidate point
        x_m, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_m, tol, sol, stepper; check=false)
        h_m = guard(sys, x_m)

        if abs(h_m) < tol
            return tₖ + τ_m, x_m
        end

        # Update bracket and boundary derivatives
        dx_m = f(x_m, tₖ + τ_m)
        hp_m = guard_derivatives(sys, x_m, dx_m, h_m)

        if signbit(h_l) != signbit(h_m)
            τ_r, x_r, h_r, hp_r = τ_m, x_m, h_m, hp_m
        else
            τ_l, x_l, h_l, hp_l = τ_m, x_m, h_m, hp_m
        end
    end

    return tₖ + τ_m, x_m
end

#Bisection Method (Iterative)
function locate_event(::BisectionLocator, prob, solver::AbstractODESolver, f, Df, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys
    τ_l, τ_r = 0.0, Δt
    h_l = h_now
    
    for _ in 1:100 # max_iter
        if (τ_r - τ_l) < tol break end

        #test midpoint
        τ_m = (τ_l + τ_r) / 2.0

        x_m, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_m, tol, sol, stepper; check=false)
        h_m = guard(sys, x_m)
    
        if signbit(h_l) != signbit(h_m)
            τ_r = τ_m
        else
            τ_l = τ_m
            h_l = h_m
        end
    end
    
    t_star = tₖ + τ_l
    x_star, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_l, tol, sol, stepper; check=false)
    return t_star, x_star
end

function locate_event(::QuadraticLocator, prob, solver::AbstractODESolver, f, Df, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys

    # Get three points 
    # Start of interval guard value. 
    h₀ = h_now

    # middle point. We just take a half step instead of going one before the start point. I think itll be more stable. 
    x₁, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, Δt / 2.0, tol, sol, stepper; check=false)
    h₁ = guard(sys, x₁)

    # endpoint
    x₂, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, Δt, tol, sol, stepper; check=false)
    h₂ = guard(sys, x₂)

    # Must actually bracket a root
    # Check if the start and end points have same sign (meaning they likely dont have a root between them)
    if signbit(h₀) == signbit(h₂)
        # Compute linear interp guess as failsafe
        θ = -h₀ / (h₂ - h₀)
        # Ensure guessed fraction doesnt go outside the time interval 
        τ_star = clamp(θ * Δt, 0.0, Δt)

        # Step solver to time we got above. 
        x_star, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_star, tol, sol, stepper; check=false)

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
    x_test, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_star, tol, sol, stepper; check=false)
    h_test = guard(sys, x_test)

    # If prediction is worse than midpoint we abandon it and linear interp. Just a failsafe.
    if abs(h_test) > abs(h₁)
        θ = -h₀ / (h₂ - h₀)
        τ_star = clamp(θ * Δt, 0.0, Δt)

        x_test, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_star, tol, sol, stepper; check=false)
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
        x_test, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_star, tol, sol, stepper; check=false)
    end
    # Calc the absolute time of the verified event
    t_star = tₖ + τ_star
    # Return abs time and its corresponding state. 
    return t_star, x_test
end

function locate_event(::NewtonLocator, prob, solver::AbstractODESolver, f, Df, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys

    τ_prev = 0.0
    h_prev = h_now

    τ_curr = Δt
    x_curr, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_curr, tol, sol, stepper; check=false)
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
        x_curr, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_curr, tol, sol, stepper; check=false)
        h_curr = guard(sys, x_curr)
    end
    t_star = tₖ + τ_curr
    return t_star, x_curr
end
