
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
Mechanical System
 - M(q): mass matrix
 - V(q): potential energy
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
struct MechanicalSol{T, X, DX, I, E, Z} <: AbstractHybridSolution
    t::T        # Time data
    x::X        # x = (q,p), the state and momentum
    dx::DX      # f(x) Derivative at each state x - only filled when dense_out = true
    prob::I     # Remember the problem - to aid interpolation
    event::E    # Times where an event has occurred 
    zeno::Z     # Times of Zeno points
end

# Function to initialize solution struct
function MechanicalSol(prob)
    return MechanicalSol([prob.tspan[1]],
        [prob.init],
        Vector{Vector{Float64}}(),      
        prob,
        Float64[],
        Float64[])
end

# INTERNAL
# Default reset map: specular reflection with coefficient of restitution. 
# Ames, Aaron & Zheng, Haiyang & Gregg, Robert & Sastry, Shankar. (2006). Is there life after Zeno? Taking executions past the breaking (Zeno) point. 2006. 6 pp.. 10.1109/ACC.2006.1656623
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
    denom = dot(normal, Mq \ normal)

    # Full equation
    pnew = p - (1+e) * dot(p, Mq\normal) / denom * normal

    return vcat(q, pnew)
end

######
### WC: Incorporate this
######
# function locate_event_mechanical(f, z, Δt, h)
function crossed_guard_mechanical(h_now, h_next, t_now, t_next)
    # Two point crossing condition
    if h_now > 0 && h_next < 0
        t_root = t_now - h_now * (t_next-t_now) / (h_next-h_now)
        return true, t_root
    else
        return false, NaN
    end
end
function locate_event_mechanical(solver, prob, f, z, tₖ, Δt, tol, sol, h)
    # We will implement a linear finder as that is what triggered 'crossed_guard'
    t₀, t₁ = 0.0, Δt
    # The event function as a function of time
    # E(t) = h(solver(f,z,t))
    E(t) = h(take_step(solver, prob, f, z, tₖ, t, tol, sol; check=false)[1])

    # Loop until convergence
    for _ ∈ 1:100
        t₂ = (t₀*E(t₁) - t₁*E(t₀)) / (E(t₁) - E(t₀))
        # Keep the two points we like
        if E(t₀)*E(t₂) > 0
            t_old, t_new = t₁, t₂
        else
            t_old, t_new = t₀, t₂
        end
        t₀, t₁ = t_old, t_new
        # Have we converged? The tolerence is currently hard coded.
        if abs(E(t₁)) < 1e-3
            break
        end
    end
    # What if we have failed to converge?
    if abs(E(t₁)) > 1e-3
        @warn "Failed to converge to an event time"
    end
    # The predicted impact time
    return t₁
end


function solve(prob::prob{S, I, T};
               solver::AbstractODESolver=RK4(),
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial = 0.01, max_iter = 10^6, 
               tol = 1e-6, ztol = 1e-3,
               guard_direction = default_guard_direction(prob.sys),
               kwargs...) where {S<:MechanicalSystem, I, T}
    
    sys = prob.sys
    h = sys.guard
    ∇h = sys.normal
    Δ = sys.reset
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
    f(x,t) = f_λ(x[1:div(length(x), 2)], x[(div(length(x),2)+1):end], 0.0)

    t_start, t_end = prob.tspan     # Extract start and end times for bounds

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

    # Actually solve now

        tₖ = sol.t[end] # Current time at start of step
        xₖ = sol.x[end]
        qₖ = xₖ[1:div(length(xₖ), 2)]
        pₖ = xₖ[(div(length(xₖ), 2) + 1):end]

        # Recall that we want h(z)≈0 and -ε<dh(q)̇q<0
        # First, is this a (post) Zeno state?
        if h(qₖ) < ztol  &&  abs(dot(∇h(qₖ), M(qₖ) \ pₖ))<ztol #-ztol < dot(∇h(qₖ), M(qₖ) \ pₖ) < 0
            # Does λ preserve the constraint?
            function guard_error(λ)
                F(z, t) = f_λ(z[1:div(length(z), 2)],z[(div(length(z), 2) + 1):end], λ)
                x_next, _, _ = take_step(solver, prob, F, xₖ, tₖ, Δt, tol, sol; check=false)
                q_next, p_next = x_next[1:div(length(x_next), 2)], x_next[(div(length(x_next), 2) + 1):end]
                
                # Constraint?
                return dot(∇h(q_next), M(q_next) \ p_next)
            end

            # If ∇h(q)̇q>0 (with λ=0), then we are moving to the interior of the state-space and the constraint need not be applied
            if guard_error(0) > 0
                F(z, t) = f_λ(z[1:div(length(z), 2)],z[(div(length(z), 2) + 1):end], 0.0)
                x_next, _, _ = take_step(solver, prob, F, xₖ, tₖ, Δt, tol, sol; check=false)
            else
                # This is the fun part; we need to actually solve for λ
                # We will solve for λ via the method of false position. This needs two initial guesses for λ
                # λ₀ = 0, because why not.
                # λ₁ = the answer predicted by applying symplectic Euler. 
                λ₀ = 0.0
                λ₁ = dot(∇h(qₖ), M(qₖ)\(-p_dot(qₖ,pₖ,0)-pₖ)) / dot(∇h(qₖ), M(qₖ)\∇h(qₖ))
                # Repeat until we have an answer
                for _ ∈ 1:100
                    λ₂ = (λ₀*guard_error(λ₁) - λ₁*guard_error(λ₀)) / (guard_error(λ₁) - guard_error(λ₀))
                    # Determine which two to keep
                    if guard_error(λ₀) * guard_error(λ₁) > 0
                        λ_old, λ_new = λ₁, λ₂
                    else
                        λ_old, λ_new = λ₀, λ₂
                    end
                    λ₀, λ₁ = λ_old, λ_new
                    # Have we converged? The tolerence is currently hard coded.
                    if abs(guard_error(λ₁)) < 1e-3
                        break
                    end
                end
                # We have our multiplier!
                f_constrained(z, t) = f_λ(z[1:div(length(z), 2)],z[(div(length(z), 2) + 1):end], λ₁)
                # x_next = solver(F, xₖ, Δt)
                x_next, _, _ = take_step(solver, prob, f_constrained, xₖ, tₖ, Δt, tol, sol; check=false)
            end            
            Δt_found = Δt

            # No sliding occurs, 'normal' hybrid situation
            else
                # Propose a step assuming no impacts
                # z_proposed = solver(f, xₖ, Δt)
                z_proposed, _, _ = take_step(solver, prob, f, xₖ, tₖ, Δt, tol, sol; check=false)

                # Is there a crossing detected?
                if crossed_guard_mechanical(h(qₖ), h(z_proposed[1:div(length(z_proposed), 2)]), 0.0, Δt)[1]
                    #Δt_found = locate_event_mechanical(f, xₖ, Δt, h)
                    Δt_found = locate_event_mechanical(solver, prob, f, xₖ, tₖ, Δt, tol, sol, h)
                    # z_impact = solver(f, xₖ, Δt_found)
                    z_impact, _, _ = take_step(solver, prob, f, xₖ, tₖ, Δt_found, tol, sol; check=false)
                    x_next = Δ(z_impact, M, ∇h, sys)
                else
                    x_next = z_proposed
                    Δt_found = Δt
                end
            end

        # Record
        if dense_out
            push!(sol.dx, f(xₖ, tₖ)) # hey, it works (usually)
        end
        push!(sol.t, Δt_found+tₖ)
        push!(sol.x, x_next)
    end

    return sol
end