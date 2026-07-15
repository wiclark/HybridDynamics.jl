# General stochastic system
struct StochasticSystem{F,G,H,N,D} <: AbstractHybridSystem
    f::F          # Drift dynamics
    g::G          # Diffusion term
    h::H          # Guard/event function
    Δ::D          # Reset map
    normal::N     # Normal to the guard, ∇G
    direction::Int # Direction of impacts on the guard
end

# EXTERNAL
# Make the the guard, reset map, and coefficient of restitution optional
"""
    StochasticSystem(f, g, h, Δ; 
    normal = nothing, direction=0)

Construct a stochastic system of the form:

"""
function StochasticSystem(f, g, h, Δ; 
    normal = nothing, direction=0)

    if isnothing(h) && !isnothing(normal)
        error("Normal to guard was provided, but a guard was not")
    end

    if !isnothing(h) && isnothing(normal)
        normal = q -> ForwardDiff.gradient(h, q)
    end

    # Wrap g to ensure its output is strictly a matrix
    g_safe = (args...) -> begin
        res = g(args...)
        if res isa Number
            return [res;;] # Convert scalar to 1x1 matrix
        elseif res isa AbstractVector
            return reshape(res, length(res), 1) # Convert 1D Vector to Nx1 Matrix
        else
            return res # Already a matrix, pass it through
        end
    end

    return StochasticSystem(
        f, g_safe, h, Δ, normal, direction
    )
end

# General solution struct for mechanical systems
struct StochasticSol{T, X, DX, I, E, EI} <: AbstractHybridSolution
    t::T        # Time data
    x::X        # x, the state
    dx::DX      # f(x) Derivative at each state x - only filled when dense_out = true
    prob::I     # Remember the problem - to aid interpolation
    event_times::E    # Times where an event has occurred
    event_indices::EI #Indices where an event has occurred
end

# Function to initialize solution struct
function StochasticSol(prob)
    return StochasticSol([prob.tspan[1]],
        [prob.init],
        Vector{Vector{Float64}}(),      
        prob,
        Float64[],
        Int[])
end

#####################################################
function take_step_stochastic!(solver, prob::prob{S, I, T}, Δt, 
    tol, sol, event_method;
    guard_direction=default_guard_direction(prob.sys)) where {S<:StochasticSystem, I, T}

    @assert event_method isa LinearLocator "Stochastic systems only support LinearLocator."

    # Extract out the state
    xₖ, tₖ = sol.x[end], sol.t[end]
    # Extract out the problem details
    sys = prob.sys
    h = sys.h
    Δ = sys.Δ

    f, g = sys.f, sys.g

    # This is going to be quite naive
    x_predict, _, _, _, _ = take_step(solver, prob, f, g, xₖ, tₖ, Δt, tol, sol; check=false)

    # Record information
    push!(sol.x, x_predict)
    push!(sol.t, tₖ + Δt)
    # Is there an impact?
    valid_linear(h1, h2) = 
        (guard_direction == 0 && h1 * h2 < 0) ||
        (guard_direction == -1 && h1 > 0 && h2 < 0) ||
        (guard_direction == 1 && h1 < 0 && h2 > 0)
    if valid_linear(h(xₖ), h(x_predict))
        x⁺ = Δ(x_predict)
        push!(sol.x, x⁺)
        push!(sol.t, tₖ + Δt)
        push!(sol.event_times, tₖ + Δt)
        push!(sol.event_indices, length(sol.t))
    end
    return sol.x[end], Δt
end
#####################################################

"""
    solve(prob::prob{S, I, T},
               solver::AbstractODESolver=EulerMaruyama();
               dense_out = false,
               dt_initial = 0.01, max_iter = 10^5, 
               tol = 1e-6, guard_direction=default_guard_direction(prob.sys),
               kwargs...) where {S<:StochasticSystem, I, T}

Solve a stochastic hybrid system.

"""
function solve(prob::prob{S, I, T},
               solver::AbstractODESolver=EulerMaruyama();
               event_method=LinearLocator(),
               dense_out = false,
               dt_initial = 0.01, max_iter = 10^5, 
               tol = 1e-6, guard_direction=default_guard_direction(prob.sys),
               kwargs...) where {S<:StochasticSystem, I, T}

    # Initialize solution struct
    sol = StochasticSol(prob)

     _, t_end = prob.tspan           # Extract the terminal time of the problem
    Δt = dt_initial
    iter = 0

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
        Δt = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

    # Actually solve now

        take_step_stochastic!(solver, prob, Δt, tol, sol, event_method; guard_direction = guard_direction)

    end

    return sol

end