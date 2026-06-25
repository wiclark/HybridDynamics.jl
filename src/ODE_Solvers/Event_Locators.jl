#====================================#
#Event detection utility. 
#Assigns a default crossing direction to a specific system to reduce issues (this may be made better for the front end later but for now this works)
default_guard_direction(sys::MechanicalSystem) = :falling
default_guard_direction(sys::NonholonomicSystem) = :falling
#Gen system lack specific constraints so we monitor crossings in both directions
default_guard_direction(sys::GeneralSystem) = :both
default_guard_direction(sys) = :both

#Wrapper function to interface between the system state and core logic
function crossed_guard(sys, h_prev, h_now, h_next, t_prev, t_now, t_next; 
                       tol=1e-6, direction=default_guard_direction(sys))   
    # Calls evaluator with the directionality
    return evaluate_crossing(h_prev, h_now, h_next, t_prev, t_now, t_next, direction; tol=tol)
end

#Core engine: determines if/when the guard function 'h' changes sign. 
function evaluate_crossing(h_prev, h_now, h_next, t_prev, t_now, t_next, direction::Symbol; tol=1e-6)
    #Helper to validate a linear sign change based on the required direction
    valid_linear(h1, h2) = 
        (direction == :both && h1 * h2 < 0) ||
        (direction == :falling && h1 > 0 && h2 < 0) ||
        (direction == :rising && h1 < 0 && h2 > 0)

    #First check: Simple linear crossing detection. Between previous and current
    if valid_linear(h_prev, h_now)
        #Linear interpolation to find the root
        t_root = t_prev - h_prev * (t_now - t_prev) / (h_now - h_prev)
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
        if abs(a) > 1e-6
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
                    (direction == :both) ||
                    (direction == :falling && h_prime(r) < 0) ||
                    (direction == :rising && h_prime(r) > 0)

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

function locate_event(::LinearLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys
    τ_l, τ_r = 0.0, Δt
    h_l = h_now

    if abs(h_now) < tol || h_now < 0
        return tₖ, xₖ
    end

    #Get right side of boundary
    x_r, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_r, tol, sol, stepper)
    h_r = guard(sys, x_r)

    τ_star = Δt
    x_star = x_r

    for _ in 1:100
        if abs(τ_r - τ_l) < tol || abs(h_r) < tol
            τ_star = τ_r
            x_star = x_r
            break
        end

        #Linear Interp
        τ_m = τ_r - h_r * (τ_r - τ_l) / (h_r - h_l)

        x_m, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_m, tol, sol, stepper)
        h_m = guard(sys, x_m)

        #keep root bracketed
        if signbit(h_l) != signbit(h_m)
            τ_r = τ_m
            h_r = h_m
        else 
            τ_l = τ_m
            h_l = h_m
        end
    end
    t_star = tₖ + τ_star
    return t_star, x_star
end

function locate_event(::QuadraticLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
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

function locate_event(::NewtonLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys

    τ_prev = 0.0
    h_prev = h_now

    τ_curr = Δt
    x_curr, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_curr, tol, sol, stepper)
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
        x_curr, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_curr, tol, sol, stepper)
        h_curr = guard(sys, x_curr)
    end
    t_star = tₖ + τ_curr
    return t_star, x_curr
end
