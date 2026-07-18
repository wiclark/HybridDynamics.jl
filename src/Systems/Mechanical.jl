
# General mechanical system
struct MechanicalSystem{M,V,G,N,R,E} <: AbstractHybridSystem
    M::M          # Mass matrix function M(q)
    V::V          # Potential energy function V(q)
    guard::G      # guard/event function
    normal::N     # Normal to the guard, ∇G
    reset::R      # reset map
    e::E          # coefficient of restitution
    direction::Int # Directional Support for Guard
end

# EXTERNAL
# Make the the guard, reset map, and coefficient of restitution optional; default to fully elastic specular reflection
"""
    MechanicalSystem(M, V; guard=nothing, normal=nothing, reset=specular_refl,
                     e=1.0, direction=-1)

Create a mechanical hybrid system of the form:
```math
\\begin{cases}
M(q) \\ddot{q} + \\nabla V(q) = 0, \\quad & h(q) \\neq 0 \\
\\text{reset}, & h(q) = 0
\\end{cases}
```

"""
function MechanicalSystem(M, V;
                guard = nothing,
                normal = nothing,
                reset = (x, Mfun, dh, sys::MechanicalSystem) -> specular_refl(x, Mfun, dh, sys),
                e = 1.0,
                direction = -1)

    if isnothing(guard) && !isnothing(normal)
        error("Normal to guard was provided, but a guard was not")
    end

    if !isnothing(guard) && isnothing(normal)
        normal = q -> ForwardDiff.gradient(guard, q)
    end

    return MechanicalSystem{typeof(M), typeof(V), typeof(guard), typeof(normal), typeof(reset), typeof(e)}(
        M, V, guard, normal, reset, e, direction
    )
end

# General solution struct for mechanical systems
"""

Solution object for mechanical systems.
"""
struct MechanicalSol{T, X, DX, I, E, EI, Z} <: AbstractHybridSolution
    t::T        # Time data
    x::X        # x = (q,p), the state and momentum
    dx::DX      # f(x) Derivative at each state x - only filled when dense_out = true
    prob::I     # Remember the problem - to aid interpolation
    event_times::E    # Times where an event has occurred 
    event_indices::EI #Indices of where an event occurred
    zeno::Z     # Times of Zeno points
end

# Function to initialize solution struct
function MechanicalSol(prob)
    return MechanicalSol([prob.tspan[1]],
        [prob.init],
        Vector{Vector{Float64}}(),      
        prob,
        Float64[],
        Int[],
        Float64[])
end

# INTERNAL
# Default reset map: specular reflection with coefficient of restitution. 
# Ames, Aaron & Zheng, Haiyang & Gregg, Robert & Sastry, Shankar. (2006). 
# Is there life after Zeno? Taking executions past the breaking (Zeno) point. 2006. 6 pp.. 10.1109/ACC.2006.1656623
function specular_refl(x, M, dh, sys)

    e = sys.e
    n = length(x) ÷ 2

    q = x[1:n]          # positions
    Mq = M(q)           # Mass matrix
    p = x[n+1:end]      # velocities

    # Constraint normal (row -> column) # <- It just is a column vector
    normal = dh(q)

    # Denominator
    # println(typeof(normal))
    # println("q = ", q)
    # println("p before = ", p)
    # println("normal = ", normal)

    Minv_n = Mq \ normal

    denom = dot(normal, Minv_n)

    # Full equation
    pnew = p - (1+e) * dot(p, Minv_n) / denom * normal

    # impulse = (1+e) * dot(p, Minv_n)/denom
    # println("impulse = ", impulse)
    # println("expected dp = ", -impulse .* normal)
    # println("new p = ", pnew)

    return vcat(q, pnew)
end

# We can solve for the multiplier directly (this does require differentiation)
function find_multiplier(prob::prob{S, I, T}) where {S<:MechanicalSystem, I, T}
    sys = prob.sys
    ∇h = sys.normal
    M(q) = sys.M(q)
    V(q) = sys.V(q)
    # Create vector field for ODE solving
    H(q,p) = 1/2*dot(p, M(q) \ p) + V(q)
    # Create the (unforced) vector fields
    q_dot(q, p) = ForwardDiff.gradient(p -> H(q,p),p)
    p_dot(q, p) = ForwardDiff.gradient(q -> -H(q,p), q)
    # Work on the multiplier
    # The function below needs to be treated as a vector
    ρ(q) = length(q)>1 ? ∇h(q)' * inv(M(q)) : [∇h(q)' * inv(M(q))]
    dρ(q) = ForwardDiff.jacobian(q->ρ(q), q)
    λ(q,p) = 1/dot(∇h(q), M(q)\∇h(q)) * ( -dot(ρ(q), p_dot(q,p)) - dot(dρ(q)*q_dot(q,p), p) )
    return λ
end
###############################################

# Internal
function guard(sys::MechanicalSystem, x::AbstractArray)
    x_phys = (x isa AbstractMatrix) ? x[:, 1] : x
    val = sys.guard(x_phys)
    # Returns the most negative value <-- Why is this all required here?
    return val isa AbstractVector ? minimum(val) : val
end

###############################################
function take_step_mechanical!(solver, prob::prob{S, I, T}, f_λ, Δt,
    tol, ztol, sol; stepper::AbstractODESolver=ModifiedMidpoint(), dense_out = true, event_method::AbstractEventLocator=LinearLocator(), 
    guard_direction=default_guard_direction(prob.sys)) where {S<:MechanicalSystem, I, T}
    # Extract out the state
    xₖ, tₖ = sol.x[end], sol.t[end]
    n = length(xₖ) ÷ 2
    qₖ, pₖ = xₖ[1:n], xₖ[n+1:end]
    # Extract out the problem details
    sys = prob.sys
    h = sys.guard
    ∇h = sys.normal
    M(q) = sys.M(q)
    V(q) = sys.V(q)
    Δ = sys.reset

    # Check the guard direction
    if guard_direction ≠ -1
        @warn "Invalid guard direction for mechanical systems"
    end
    # Check to make sure we're on the right side of the guard
    if h(qₖ) < ztol
        # If we're moving in, we should be moving out and count it as a reset
        if dot(∇h(qₖ), M(qₖ) \ pₖ) < 0
            x_new = Δ(xₖ, M, ∇h, sys)
            push!(sol.t, tₖ)
            push!(sol.x, x_new)
            if dense_out
                push!(sol.dx, f_λ(x_new[1:n], x_new[n+1:end], 0.0))
            end
            # Update the state
            qₖ, pₖ = x_new[1:n], x_new[n+1:end]
        end
    end
    # Determine whether or not we are on the sliding/post-Zeno regime
    if h(qₖ) < ztol  &&  abs(dot(∇h(qₖ), M(qₖ) \ pₖ)) < ztol
        # To what degree does λ preserve the holonomic constraint?
        function guard_error(λ)
            F(z, t) = f_λ(z[1:n], z[n+1:end], λ)
            x_predict, _, _, dt_used, dt_next = take_step(solver, prob, F, vcat(qₖ,pₖ), tₖ, Δt, tol, sol; check=false)
            q_next, p_next = x_predict[1:n], x_predict[n+1:end]
            # (Tangent) constraint violation
            return dot(∇h(q_next), M(q_next) \ p_next)
        end
        # If ∇h(q)̇q>0 with λ=0, then we are escaping the guard (inwards) and are escaping the sliding mode
        if guard_error(0.0) > 0
            if solver isa LMM
                push!(sol.event_indices, length(sol.x))
            end
            F(z, t) = f_λ(z[1:n], z[n+1:end], 0.0)
            x_predict, _, _, dt_used, dt_next = take_step(solver, prob, F, vcat(qₖ,pₖ), tₖ, Δt, tol, sol; check=false)
            # Record the derivative
            if dense_out
                push!(sol.dx, F(x_predict, tₖ+dt_used))
            end
        else
            if solver isa LMM
                push!(sol.event_indices, length(sol.x))
            end
            # We have the expression for the multiplier
            λ_constrained = find_multiplier(prob)
            F2(z, t) = f_λ(z[1:n], z[n+1:end], λ_constrained(z[1:n], z[n+1:end]))
            x_predict, _, _, dt_used, dt_next = take_step(solver, prob, F2, vcat(qₖ,pₖ), tₖ, Δt, tol, sol; check=false)
            # Record the derivative
            if dense_out
                push!(sol.dx, F2(x_predict, tₖ+dt_used))
            end
        end
        # Collect our results
        push!(sol.x, x_predict)
        push!(sol.t, tₖ+dt_used)
        push!(sol.zeno, tₖ)
        return x_predict, dt_used, dt_next, true
    else # No Zeno stuff is present
        f(z, t) = f_λ(z[1:n], z[n+1:end], 0.0)
        x_predict, eventtrigger, t_root, dt_used, dt_next = take_step(solver, prob, f, vcat(qₖ, pₖ), tₖ, Δt, tol, sol)
        # Was there an impact?
        if eventtrigger
            t_star, x_star = locate_event(event_method, prob, solver, f, vcat(qₖ, pₖ), tₖ, Δt, guard(sys, xₖ), tol, sol, stepper)
            x_predict = Δ(x_star, M, ∇h, sys)

            push!(sol.event_times, t_star)

            if solver isa LMM
                push!(sol.event_indices, length(sol.x))
            end
            push!(sol.t, t_star)
            push!(sol.x, x_predict)

            if dense_out
                push!(sol.dx, f(x_predict, t_star))
            end

        else
            push!(sol.x, x_predict)
            push!(sol.t, tₖ+dt_used)
            if dense_out
                push!(sol.dx, f(x_predict, tₖ+dt_used))
            end
        end
        return x_predict, dt_used, dt_next, false
    end
end

###############################################
"""
    solve(prob; kwargs...)

Solve a mechanical hybrid system.
"""
function solve(prob::prob{S, I, T},
               solver::AbstractODESolver=RK4();
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial = 0.01, max_iter = 10^6, 
               tol = 1e-6, ztol = 1e-3,
               guard_direction = default_guard_direction(prob.sys),
               kwargs...) where {S<:MechanicalSystem, I, T}
    
    sys = prob.sys
    ∇h = sys.normal
    M(q) = sys.M(q)
    V(q) = sys.V(q)

    # Initialize solution struct
    sol = MechanicalSol(prob)

    # Create vector field for ODE solving
    H(q,p) = 1/2*dot(p, M(q) \ p) + V(q)
    # λ is the unknown multiplier to enforce the constraint ∇h(q)̇q = ∇h(q)M(q)\p = 0
    q_dot(q, p) = ForwardDiff.gradient(p -> H(q,p),p)
    p_dot(q, p, λ) = ForwardDiff.gradient(q -> -H(q,p), q) .+ λ*∇h(q)
    # Combining these vector fields together
    f_λ(q, p, λ) = [q_dot(q, p); p_dot(q, p, λ)]

    _, t_end = prob.tspan           # Extract the terminal time of the problem

    Δt = dt_initial                 # Initialize current time step with user input
    iter = 0                        # Start iteration counter

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

        # Perform the step
        _, _, Δt, _ = take_step_mechanical!(solver, prob, f_λ, Δt, tol, ztol, sol; 
                        dense_out = dense_out, event_method=event_method, guard_direction = guard_direction)
    end

    return sol
end