
struct FilippovSys{F, G, H, N} <: AbstractHybridSystem
    F::F    # Function one, H(x) > 0
    G::G    # Function two, H(x) < 0
    H::H    # Guard
    N::N    # Normal to the guard, ∇H
end

# Default to auto diff to find ∇H
function FilippovSys(F, G, H; N= (x-> ForwardDiff.gradient(H,x)))
    return FilippovSys(F, G, H, N)
end
struct FilippovSol{T, X, S}
    T::T    # Time data
    X::X    # Position data
    S::S    # When sliding perhaps?
end

######
### WC: What kind of data is S? Intervals? Start-end points? A collection of ordered pairs?
# CK: Idk yet, I haven't actually put anything there
######

# Constructor for solution struct
function FilippovSol(T, X; S = NaN)
    return FilippovSol(T, X, S)
end

# Initialize solution struct
function initsol(prob::prob{FilippovSys})
    return FilippovSol([prob.tspan[1]], [prob.init])
end

# INTERNAL
# Returns a vector field at state `x` for a Filippov system
######
### WC: This function confused me for a bit. It is taking in a location x and returning the vector field *function* where x is located. This should be made clearer (via comments).
######
function filippov_vector_field(sys, x;
        Ftol=1e-7,
        atol=1e-7)

    F, G, H, N = sys

    h = H(x)

# Away from the guard
    if h > Ftol
        return F
    elseif h < -Ftol
        return G
    end

# Near the guard

    a(x) = dot(N, F(x))
    b(x) = dot(N, G(x))

    # Attracting sliding
    if a < -atol && b > atol

        λ(x) = a(x) / (a(x) - b(x))
        return x -> (1 - λ(x)) * F(x) + λ(x) * G(x)
    end

    # Repelling sliding (non-unique, but still gotta go somewhere)
    if a > atol && b < -atol

        λ(x) = a(x) / (a(x) - b(x))
        return x -> (1 - λ(x)) * F(x) + λ(x) * G(x)
    end

    # Direct crossings
    if a > atol && b > atol
        return F
    end
    if a < -atol && b < -atol
        return G
    end

# Otherwise
    # fallback: choose least transverse field
    if abs(a) < abs(b)
        return G
    else
        return F
    end
end

# EXTERNAL
# Filippov-specific solve
function solve(prob::prob{<:FilippovSys}, solver; dt_initial = 0.01, max_iter = 10^6, tol = 1e-6, kwargs...)
    
    sys = prob.sys
    F = sys.F
    G = sys.G
    H = sys.H
    sol = initsol(prob)

    t_start, t_end = prob.tspan     # Extract start and end times for bounds

    Δt = dt_initial                 #Initialize current time step with user input
    iter = 0                        #Start iteration counter

    # Run sim until end of specified time span
    while sol.T[end] < t_end 
        
    # Safties
        # Stop if we hit the iteration limit to avoid memory doomsday
        iter += 1
        if iter > max_iter 
            @warn "Maximum Iteration Count ($max_iter) exceeded."
            break
        end

        # Terminate if the remaining time is below machine precision
        if t_end - sol.T[end] <= eps(t_end)
            break
        end

        #Truncate time step if we overshoot the final sim time
        dt_step = (sol.T[end] + Δt > t_end) ? (t_end - sol.T[end]) : Δt

    # Actually solve now

        xₖ = sol.X[end] #Retrieve current state at start of step
        tₖ = sol.T[end] #Retrieve current time at start of step

        # Choose vector field for current step based on current position relative to the guard
        vf = Filippov_vec_field(sys, x)
        x_predict, event_triggered, h_now = take_step(solver, prob, vf, xₖ, tₖ, dt_step, tol, sol)

        t_next = tₖ + dt_step
        tₖ += dt_step
        xₖ = x_predict
        push!(sol.T, tₖ)
        push!(sol.X, xₖ)

        dt = min(dt_next, dt_initial)
        
        ######
        ### WC: You have a handful of variables that are never used. You can use _ 
        ######

    end

    return sol
end