
struct FilippovSys{F, G, H}
    F::F    # Function one, H(x) > 0
    G::G    # Function two, H(x) < 0
    H::H    # Guard
end

struct FilippovSol{T, X, S}
    T::T    # Time data
    X::X    # Position data
    S::S    # When sliding perhaps?
end

# Constructor for solution struct
function FilippovSol(T, X; S = NaN)
    return FilippovSol(T, X, S)
end

# Initialize solution struct
function initsol(prob::prob{FilippovSys})
    return FilippovSol([prob.tspan[1]], [prob.init])
end

# Filippov-specific solve
function solve(prob::prob{FilippovSys}, solver; tol=1e-10, max_iter=10^6, kwargs...) 

    sys = prob.sys
    F, G, H = sys   # Unpack system

    sol = initsol(prob)

    t_start, t_end = prob.tspan     # Extract start and end times for bounds

    Δt = dt_initial                 #Initialize current time step with user input
    iter = 0                        #Start iteration counter

    # Run sim until end of specified time span
    while sol.t[end] < t_end 
        
    # Safties
        # Stop if we hit the iteration limit to avoid memory doomsday
        iter += 1
        if iter > max_iter 
            @warn "Maximum Iteration Count ($max_iter) exceeded."
            break
        end

        # Terminate if the remaining time is below machine precision
        if t_end - sol.t[end] <= eps(t_end)
            break
        end

        #Truncate time step if we overshoot the final sim time
        Δt_step = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

    # Actually solve now

        xₖ = sol.x[end] #Retrieve current state at start of step
        tₖ = sol.t[end] #Retrieve current time at start of step





        x_predict, eventtrigger, h_now = take_step(solver, sys, something, xₖ, tₖ, Δt_step, tol)
        t_next = tₖ + Δt_step

    # Store solution 
        push!(sol.x, x_predict)
        push!(sol.t, t_next)
        Δt = dt_initial
        
    end

    return sol
end