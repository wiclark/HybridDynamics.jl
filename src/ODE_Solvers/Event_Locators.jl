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

# HERMITE VERSION
function crossed_guard(event_method::LinearHermite, sys, f, sol, xₖ, tₖ, x_next, Δt;
                       tol=1e-6, direction=default_guard_direction(sys), guard_derivative=nothing)
    h_now  = guard(sys, xₖ)
    h_next = guard(sys, x_next)

    hp_now  = guard_derivatives(sys, xₖ,     f(xₖ, tₖ),          h_now; guard_derivative=guard_derivative)
    hp_next = guard_derivatives(sys, x_next, f(x_next, tₖ + Δt), h_next; guard_derivative=guard_derivative)

    return evaluate_crossing(event_method, h_now, h_next, hp_now, hp_next, tₖ, tₖ + Δt, direction; tol=tol)
end

# LINEAR/QUADRATIC VERSION
function crossed_guard(event_method::LinearQuadratic, sys, f, sol, xₖ, tₖ, x_next, Δt;
                       tol=1e-6, direction=default_guard_direction(sys), guard_derivative=guard_derivative)
    h_now  = guard(sys, xₖ)
    h_next = guard(sys, x_next)

    # Needs the previous accepted point from the solution history
    idx    = max(1, length(sol.x) - 1)
    t_prev = sol.t[idx]
    h_prev = guard(sys, sol.x[idx])

    return evaluate_crossing(event_method, h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt, direction; tol=tol)
end

#JUST LINEAR VERSION
function crossed_guard(event_method::Linear, sys, f, sol, xₖ, tₖ, x_next, Δt; tol=1e-6, direction=default_guard_direction(sys), guard_derivative=guard_derivative)
    h_now = guard(sys, xₖ)
    h_next = guard(sys, x_next)

    return evaluate_crossing(event_method, h_now, h_next, tₖ, tₖ + Δt, direction; tol=tol)
end 

#Calcs dh/dt = ∇h(x) * dx
function guard_derivatives(sys, x, dx, h_val; guard_derivative=nothing)
    dhdx = if guard_derivative == nothing
        ForwardDiff.gradient(z -> guard(sys, z), x)
    else
        guard_derivative(x)
    end

    return dot(dhdx, dx)
end

#LINEAR VERSION
function evaluate_crossing(::Linear, h_now, h_next, t_now, t_next, direction::Int; tol=1e-6)
    Δt = t_next - t_now

    if Δt <= 0
        return false, NaN, NaN
    end

    valid_linear(h1, h2) = 
        (direction == 0 && h1 * h2 < 0) ||
        (direction == -1 && h1 > 0 && h2 < 0) ||
        (direction == 1 && h1 < 0 && h2 > 0)

        if valid_linear(h_now, h_next)
            τ_root = -h_now * Δt / (h_next - h_now)
            return true, t_now + τ_root, NaN
        end
    return false, NaN, NaN
end

#HERMITE VERSION
function evaluate_crossing(::LinearHermite, h_now, h_next, hp_now, hp_next, t_now, t_next, direction::Int; tol=1e-6)
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
        denom = h_b - h_a
        L = τ_b - τ_a
        
        if abs(denom) < 1e-14
            τ = τ_a + 0.5 * L
        else
            τ_secant = τ_a - h_a * (L / denom)
            τ = clamp(τ_secant, τ_a + 0.05 * L, τ_b - 0.05 * L)
        end

        for _ in 1:6
            der = dp(τ)
            abs(der) < 1e-14 && break
            τ = clamp(τ - p(τ) / der, τ_a, τ_b)
        end
        return τ
    end

    # The next function attempts to find an event based off of Hermite interpolation
    function valid_hermite()
        extrema = Float64[]
        if abs(β) < 1e-12
            if abs(α) > 1e-12
                R1 = -hp_now / (2 * α)
                if 0 < R1 < Δt; push!(extrema, R1); end
            end
        else
            D = 4 * α^2 - 12 * β * hp_now
            if D >= 0
                r1 = (-2 * α + sqrt(D)) / (6 * β)
                r2 = (-2 * α - sqrt(D)) / (6 * β)
                R1, R2 = min(r1, r2), max(r1, r2)
                if 0 < R1 < Δt; push!(extrema, R1); end
                if 0 < R2 < Δt; push!(extrema, R2); end
            end
        end
        # build sub intervals
        pts = [0.0; extrema; Δt]
        p(τ) = h_now + hp_now * τ + α * τ^2 + β * τ^3
        vals = [p(t) for t in pts]

        # check each segment for linear crossing
        for i in 1:(length(pts)-1)
            if valid_linear(vals[i], vals[i+1])
                return true, find_cubic_root(pts[i], pts[i+1])
            end
        end
        return false, NaN
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

#QUADRATIC VERSION 
function evaluate_crossing(::LinearQuadratic, h_prev, h_now, h_next, t_prev, t_now, t_next, direction::Int; tol=1e-6)    #Helper to validate a linear sign change based on the required direction
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

#Locator Dispatches
function locate_event(::LinearQuadratic, prob, solver::AbstractODESolver, f, Df, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4(); event_max_iters::Int = 100, force_iters::Bool = false, guard_derivative=nothing)
    # Extract System
    sys = prob.sys
    # Extract left boundary to 0 and right to Δt
    τ_l, τ_r = 0.0, Δt
    # Set guard function value at the left boundary to the current value.
    h_l = h_now

    # Check if state already satisfies the event condition, if so we return the state and time.
    if abs(h_now) < tol 
        return tₖ, xₖ, 0
    end

    # Get right side of boundary
    x_r, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_r, tol, sol, stepper; check=false)
    # Eval guard at this new right side state
    h_r = guard(sys, x_r)
    
    # If we ever get a step not bracketing a root we exit to avoid iterating garbage. 
    if signbit(h_l) == signbit(h_r)
        # Return the time and state at the end of the full step since no event took place
        return tₖ + Δt, x_r, 0
    end

    # Initialize our best guess for the event time offset and state to the right boundary
    τ_star = Δt
    x_star = x_r
    iterations = 0

    for i in 1:event_max_iters
        iterations = i

        if (h_r - h_l) == 0
            # avoiding division by zero
            break
        end

        # Linear Interpolation
        # Calc time offset where event occurs
        τ_m = τ_r - h_r * (τ_r - τ_l) / (h_r - h_l)

        # Step ODE solver forward by time offset τ_m
        x_m, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_m, tol, sol, stepper; check=false)
        # Eval guard at this new state
        h_m = guard(sys, x_m)

        # Check if this interpolated state satifies the event tolerance 
        if !force_iters && abs(h_m) < tol
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
    return t_star, x_star, iterations
end

#Linear Locator
function locate_event(::Linear, prob, solver::AbstractODESolver, f, Df, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK=RK4(); event_max_iters::Int=100, force_iters::Bool=false, guard_derivative=nothing)
    sys = prob.sys

    τ_l, τ_r = 0.0, Δt
    h_l = h_now

    if abs(h_now) < tol
        return tₖ, xₖ, 0
    end

    x_r, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_r, tol, sol, stepper; check=false)
    h_r = guard(sys, x_r)

    if signbit(h_l) == signbit(h_r)
        return tₖ + Δt, x_r, 0
    end

    τ_star = Δt
    x_star = x_r
    iterations = 0

    for i in 1:event_max_iters
        iterations = i

        if (h_r - h_l) == 0
            break
        end

        τ_m = τ_r - h_r * (τ_r - τ_l) / (h_r - h_l)
        x_m, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_m, tol, sol, stepper; check=false)
        h_m = guard(sys, x_m)

        if !force_iters && abs(h_m) < tol
            τ_star = τ_m
            x_star = x_m 
            break
        end

        τ_star = τ_m
        x_star = x_m

        if signbit(h_l) != signbit(h_m)
            τ_r = τ_m
            h_r = h_m
        else
            τ_l = τ_m
            h_l = h_m
        end
    end
    t_star = tₖ + τ_star
    return t_star, x_star, iterations
end

#Hermite Locator
function locate_event(::LinearHermite, prob, solver::AbstractODESolver, f, Df, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4(); event_max_iters::Int = 100, force_iters::Bool = false, guard_derivative=nothing)
    sys = prob.sys
    
    # Boundary setup
    τ_l, τ_r = 0.0, Δt
    x_l = xₖ
    h_l = h_now
    dx_l = f(x_l, tₖ)
    hp_l = guard_derivatives(sys, x_l, dx_l, h_l; guard_derivative=guard_derivative)

    if abs(h_now) < tol
        return tₖ, xₖ, 0
    end

    x_r, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_r, tol, sol, stepper; check=false)
    h_r = guard(sys, x_r)
    dx_r = f(x_r, tₖ + Δt)
    hp_r = guard_derivatives(sys, x_r, dx_r, h_r; guard_derivative=guard_derivative)

    τ_m = 0.5 * Δt
    x_m = x_r
    iterations = 0

    for i in 1:event_max_iters
        iterations = i
        dτ = τ_r - τ_l

        if !force_iters && dτ < 1e-12
            #division by zero
            break
        end

        if dτ < 1e-12
            break
        end

        # Hermite Coeffs rewritten
        α = (3 * (h_r - h_l) - (2 * hp_l + hp_r) * dτ) / (dτ^2)
        β = (2 * (h_l - h_r) + (hp_l + hp_r) * dτ) / (dτ^3)

        p(τ) = h_l + hp_l * τ + α * τ^2 + β * τ^3
        dp(τ) = hp_l + 2 * α * τ + 3 * β * τ^2

        # Initial secant guess
        denom = h_r - h_l
        τ_step = abs(denom) < 1e-14 ? 0.5 * dτ : clamp(dτ * (-h_l / denom), 0.05 * dτ, 0.95 * dτ)

        # Newton refinement
        for _ in 1:5
            val = p(τ_step)
            der = dp(τ_step)
            abs(der) < 1e-14 && break
            τ_step = clamp(τ_step - val / der, 0.05 * dτ, 0.95 * dτ)
        end

        τ_m = τ_l + τ_step

        # Step ODE solver to intermediate candidate point
        x_m, _, _, _, _ = take_step(solver, prob, f, Df, xₖ, tₖ, τ_m, tol, sol, stepper; check=false)
        h_m = guard(sys, x_m)

        if !force_iters && abs(h_m) < tol
            return tₖ + τ_m, x_m, iterations
        end

        # Update bracket and boundary derivatives
        dx_m = f(x_m, tₖ + τ_m)
        hp_m = guard_derivatives(sys, x_m, dx_m, h_m; guard_derivative=guard_derivative)

        if signbit(h_l) != signbit(h_m)
            τ_r, x_r, h_r, hp_r = τ_m, x_m, h_m, hp_m
        else
            τ_l, x_l, h_l, hp_l = τ_m, x_m, h_m, hp_m
        end
    end

    return tₖ + τ_m, x_m, iterations
end