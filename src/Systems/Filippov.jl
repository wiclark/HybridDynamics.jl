
struct FilippovSys{F, G, H, N}
    F::F    # Function one, H(x) > 0
    G::G    # Function two, H(x) < 0
    H::H    # Guard
    N::N    # Normal to the guard, ΔH
end

######
### WC: The normal is ∇H, not ΔH. 
######

# Default to auto diff to find ΔH
function FilippovSys(F, G, H; N=ForwardDiff.gradient(H,x))
    return FilippovSys(F, G, H, N)
end
struct FilippovSol{T, X, S}
    T::T    # Time data
    X::X    # Position data
    S::S    # When sliding perhaps?
end

######
### WC: What kind of data is S? Intervals? Start-end points? A collection of ordered pairs?
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
### WC: These tolerences are very small... You would want these tolerences to be commensurate with the tolerences from the ode solvers.
######
######
### WC: This function confused me for a bit. It is taking in a location x and returning the vector field *function* where x is located. This should be made clearer (via comments).
######
function filippov_vector_field(sys, x;
        Ftol=1e-10,
        atol=1e-10)

    F, G, H, N = sys

    h = H(x)

# Away from the guard
    if h > Ftol
        return F
    elseif h < -Ftol
        return G
    end

# Near the guard
    Fx = F(x)
    Gx = G(x)

    a = dot(N, Fx)
    b = dot(N, Gx)

    # Attracting sliding
    ######
    ### WC: Will this warning appear for every time step that sliding occurs?
    ### This could get very obnoxious. I would only have this trigger when sliding begins.
    ######
    if a < -atol && b > atol

        λ = a / (a - b)
        @warn "Attracting sliding mode" a b λ
        return x -> (1 - λ) * F(x) + λ * G(x)
    end
    ######
    ### WC: You're returning a function of x, not just its evaluation? As a function, this is incorrect as λ is also a function of x.
    ######

    # Repelling sliding (non-unique, but still gotta go somewhere)
    if a > atol && b < -atol

        λ = a / (a - b)
        @warn "Repelling sliding mode" a b λ
        return x -> (1 - λ) * F(x) + λ * G(x)
    end

    # Direct crossings
    if a > atol && b > atol
        return F
    end
    if a < -atol && b < -atol
        return G
    end

# Otherwise
    @warn "Degenerate something something tangency" a b h

    # fallback: choose least transverse field
    if abs(a) < abs(b)
        return G
    else
        return F
    end
end

# EXTERNAL
# Filippov-specific solve
######
### WC: We've gone over tolerences being this small....
######
function solve(prob::prob{<:FilippovSys}, solver; dt_initial = 0.01, max_iter = 10^6, tol = 1e-12, kwargs...)
    
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