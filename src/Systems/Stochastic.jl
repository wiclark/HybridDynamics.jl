# General mechanical system
struct StochasticSystem{F,H,N,D} <: AbstractHybridSystem
    f::F          # Continuous dynamics
    h::H          # Guard/event function
    Δ::D          # Reset map
    normal::N     # Normal to the guard, ∇G
end

# EXTERNAL
# Make the the guard, reset map, and coefficient of restitution optional
"""
Stochastic System
 - f
"""
function StochasticSystem(f, h, Δ; normal = nothing)

    if isnothing(h) && !isnothing(normal)
        error("Normal to guard was provided, but a guard was not")
    end

    if !isnothing(h) && isnothing(normal)
        normal = q -> ForwardDiff.gradient(h, q)
    end

    return StochasticSystem{typeof(f), typeof(h), typeof(Δ), typeof(normal)}(
        f, h, Δ, normal
    )
end

# General solution struct for mechanical systems
struct StochasticSol{T, X, DX, I, E, Z}
    t::T        # Time data
    x::X        # x, the state
    dx::DX      # f(x) Derivative at each state x - only filled when dense_out = true
    prob::I     # Remember the problem - to aid interpolation
    event::E    # Times where an event has occurred 
    zeno::Z     # Times of Zeno points
end

# Function to initialize solution struct
function StochasticSol(prob)
    return StochasticSol([prob.tspan[1]],
        [prob.init],
        Vector{Vector{Float64}}(),      
        prob,
        Float64[],
        Float64[])
end



function solve(prob::prob{S, I, T};
               solver::AbstractODESolver=RK4(),
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial = 0.01, max_iter = 10^5, 
               tol = 1e-6, ztol = 1e-4,
               kwargs...) where {S<:StochasticSystem, I, T}
    
    sys = prob.sys
    h = sys.h
    ∇h = sys.normal
    Δ = sys.Δ

    # Initialize solution struct
    sol = StochasticSol(prob)



end