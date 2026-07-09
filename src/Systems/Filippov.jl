# A Filippov (discontinuous) dynamical system
struct FilippovSystem{F, G, H, N} <: AbstractHybridSystem
    F::F    # Function one, H(x) > 0
    G::G    # Function two, H(x) < 0
    h::H    # Guard
    N::N    # Normal to the guard, ∇H
end

# Default to auto diff to find ∇H
"""
    FilippovSystem(F, G, H; N= (x-> ForwardDiff.gradient(H,x)))

Construct a Filippov system of the form:


"""
function FilippovSystem(F, G, H; N= (x-> ForwardDiff.gradient(H,x)))
    return FilippovSystem(F, G, H, N)
end

struct FilippovSol{T, X, DX, S} <: AbstractHybridSolution
    t::T    # Time data
    x::X    # Position data
    dx::DX  # f(x) Derivative at each state x - only filled when dense_out = true
    s::S    # Time indices while sliding (this still needs added)
end

function FilippovSol(prob::prob{F, I, T}) where {F<:FilippovSystem, I, T}
    return FilippovSol([prob.tspan[1]],
        [prob.init],
        Vector{Vector{Float64}}(),
        Float64[])
end

# INTERNAL
# Returns a vector field function at state `x` for a Filippov system
######
### WC: I do not believe that there is any reason for 'atol'
######
function filippov_vector_field(sys, x;
        Ftol=1e-4,
        atol=0.0)

    F = sys.F
    G = sys.G
    H = sys.h
    N = sys.N

    h = H(x)

# Away from the guard
    if h > Ftol
        return F, false
    elseif h < -Ftol
        return G, false
    end

# Near the guard

    a(y) = dot(N(y), F(y))
    b(y) = dot(N(y), G(y))

    if (a(x) < -atol && b(x) > atol) || (a(x) > atol && b(x) < -atol)
        λ(y) = a(y)/(a(y)-b(y))
        return y -> (1-λ(y))*F(y) + λ(y)*G(y), true


    elseif a(x) > atol && b(x) > atol
        return F, false

    elseif a(x) < -atol && b(x) < -atol
        return G, false
    else
        println([a(x), b(x)])
        error("Failed vector field determination")
    end
end

function guard(sys::FilippovSystem, x)
    if isnothing(sys.h)
        return nothing
    else
        return sys.h(x)
    end
end

# EXTERNAL
# Filippov-specific solve
"""
    solve(prob::prob{S,I,T}, solver::AbstractODESolver=RK45();
    dense_out = true,
    dt_initial = 0.001, max_iter = 10^6, tol = 1e-6, 
    stepper::AbstractODESolver=RK4(), event_method::AbstractEventLocator=LinearLocator(),
    kwargs...) where {S<:FilippovSystem, I, T}

Solve a Filippov system.

"""
function solve(prob::prob{S,I,T}, solver::AbstractODESolver=RK45();
    dense_out = true,
    dt_initial = 0.001, max_iter = 10^6, tol = 1e-6, 
    stepper::AbstractODESolver=RK4(), event_method::AbstractEventLocator=LinearLocator(),
    kwargs...) where {S<:FilippovSystem, I, T}
    
    sys = prob.sys
    sol = FilippovSol(prob)

    _, t_end = prob.tspan     # Extract start and end times for bounds

    Δt = dt_initial                 #Initialize current time step with user input
    iter = 0                        #Start iteration counter
    sliding_prev = false

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
        dt_step = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

    # Actually solve now

        xₖ = sol.x[end] #Retrieve current state at start of step
        tₖ = sol.t[end] #Retrieve current time at start of step

        # Choose vector field for current step based on current position relative to the guard
        vf_fun, sliding_now = filippov_vector_field(sys, xₖ)

        ## This doesn't quite work due to tolerances in the above function
        if sliding_now && !sliding_prev
            @warn "Sliding mode entered at t = $(tₖ)"
        end

        sliding_prev = sliding_now
        if sliding_now
            push!(sol.s, tₖ)
        end

        vf(x,t) = vf_fun(x)
        x_predict, _, _, dt_used, dt_next = take_step(solver, prob, vf, xₖ, tₖ, dt_step, tol, sol)

        tₖ += dt_used
        xₖ = x_predict
        if dense_out
            push!(sol.dx, vf(xₖ, tₖ))
        end
        push!(sol.t, tₖ)
        push!(sol.x, xₖ)

        Δt = dt_next

    end

    return sol
end