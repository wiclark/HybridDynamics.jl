# Computes a Nonholonomic system subject to mechanical impacts. 
# See "The Bouncing Penny and Nonholonomic Impacts" by Clark and Bloch

# The current version implements linear constraints - affine are forthcoming

# Define the Nonholonomic structure
struct NonholonomicSystem{M, V, A, G, N, R, E} <: AbstractHybridSystem
    M::M      # Mass matrix function M(q)
    V::V      # Potential energy function V(q)
    A::A      # Constraint matrix function A(q)
    guard::G  # Guard/event function (only of q)
    normal::N # Normal to the guard, ∇G
    reset::R  # The reset map
    e::E      # Coefficient of restitution
    direction::Int #Directional support 
end

# Create the system
function NonholonomicSystem(M, V;
                        A = nothing,
                        guard = nothing,
                        normal = nothing,
                        reset = (x, Mfun, Afun, ∇h, sys::NonholonomicSystem) -> specular_refl(x, Mfun, Afun, ∇h, sys),
                        e = 1,
                        direction = -1)
    # Standard error checking
    if isnothing(guard) && !isnothing(normal)
        error("Normal to guard was provided, but the guard was not.")
    end
    # Compute the gradient through automatic differentiation
    if !isnothing(guard) && isnothing(normal)
        normal = q -> ForwardDiff.gradient(guard, q)
    end
    # Create the structure
    return NonholonomicSystem{typeof(M), typeof(V), typeof(A), typeof(guard), typeof(normal), typeof(reset), typeof(e)}(
        M, V, A, guard, normal, reset, e, direction
    )
end

# General solution struct for nonholonomic systems
struct NonholonomicSol{T, X, DX, I, E, EI, Z}
    t::T     # Time data
    x::X     # x = (q,p), the state and momentum
    dx::DX   # f(x) derivatve at each state x - only filled out when dense_out = true
    prob::I  # Remember the problem - to aid interpolation
    event_times::E # Times where an event occurs
    event_indices::EI #Indices of event times
    zeno::Z  # Times of Zeno points
end

# Function to initialize the solution struct
function NonholonomicSol(prob)
    return NonholonomicSol(
        [prob.tspan[1]],
        [prob.init],
        Vector{Vector{Float64}}(),
        prob,
        Float64[],
        Int[],
        Float64[])
end

# Determine the nonholonomic reset map
function specular_refl(x, M, A, ∇h, sys)
    # We are working in the momentum coordinates
    e = sys.e
    n = length(x) ÷ 2
    q, p = x[1:n], x[n+1:end]
    # The mass, constraint matrix, and normal vector
    Mq, Aq, ∇hq = M(q), A(q), ∇h(q)
    sub_matrix = Aq * (Mq \ Aq')
    constraint_vector = vec(∇hq' * (Mq \ Aq'))
    # The denominator
    denom = dot(∇hq, Mq \ ∇hq) - dot(constraint_vector, sub_matrix \ constraint_vector)
    # The 'external' multiplier
    ε = (1+e) * (dot(constraint_vector, sub_matrix \ (Aq * (Mq\p))) - dot(p, Mq \ ∇hq)) / denom
    # The 'intermal' multipliers
    λ = -ε * inv(Aq * inv(Mq) * Aq') * Aq * inv(Mq) * ∇hq
    # The new momentum
    p_new = p + ε*∇hq + Aq'*λ
    return vcat(q, p_new)
end

function find_multiplier(prob::prob{S, I, T}; sliding=false) where {S<:NonholonomicSystem, I, T}
    sys = prob.sys
    ∇h = sys.normal
    M(q) = sys.M(q)
    V(q) = sys.V(q)
    # Is the contraint matrix sliding or not?
    A(q) = sliding ? vcat(sys.A(q), ∇h(q)') : sys.A(q)
    # Create the vector field for ODE solving
    H(q,p) = 1/2*dot(p, M(q) \ p) + V(q)
    q_dot(q, p) = ForwardDiff.gradient(p -> H(q,p),p)
    p_dot(q, p) = ForwardDiff.gradient(q -> -H(q,p), q)
    # Work out the value of the k (or k+1) multipliers
    # Differentiate the constraint and solve for the multiplier
    A_dot(q, p) = ForwardDiff.derivative(ε -> A(q .+ ε*inv(M(q))*p), 0.0)
    λ(q, p) = inv(A(q)*inv(M(q))*(A(q)')) * (-A_dot(q,p)*inv(M(q))*p + A(q)*inv(M(q))*p_dot(q,p))
    # The version below works, but is slower.
    #=
    ρ(q, p) = A(q) * inv(M(q)) *  p
    dρₓ(q,p) = ForwardDiff.jacobian(q -> ρ(q,p), q) * q_dot(q, p)
    dρₚ(q,p) = ForwardDiff.jacobian(p -> ρ(q,p), p) * p_dot(q, p)
    λ(q, p) = -inv(A(q) * inv(M(q)) * A(q)') * (dρₓ(q, p) + dρₚ(q, p))
    =#
    return λ
end

# Internal
function guard(sys::NonholonomicSystem, x::AbstractArray)
    x_phys = (x isa AbstractMatrix) ? x[:,1] : x
    val = sys.guard(x_phys)
    return val isa AbstractVector ? minimum(val) : val
end

# One-sided guard crossing detection
function crossed_guard_nonholonomic(h_now, h_next, t_now, t_next)
    # Two point condition
    if h_now > 0 && h_next < 0
        t_root = t_now - h_now * (t_next-t_now) / (h_next-h_now)
        return true, t_root
    else
        return false, NaN
    end
end

##############################################################

##############################################################

function take_step_nonholonomic!(solver, prob::prob{S, I, T}, f_λ, Δt,
    tol, ztol, sol; stepper::AbstractODESolver=ModifiedMidpoint(), dense_out=true, event_method=AbstractEventLocator=LinearLocator(),
    guard_direction=default_guard_direction(prob.sys)) where {S<:NonholonomicSystem, I, T}
    # Extract out the state
    xₖ, tₖ = sol.x[end], sol.t[end]
    n = length(xₖ) ÷ 2
    qₖ, pₖ = xₖ[1:n], xₖ[n+1:end]
    # Extract out the problem details
    sys = prob.sys
    h, ∇h = sys.guard, sys.normal
    M(q) = sys.M(q)
    V(q) = sys.V(q)
    A(q) = sys.A(q)
    Δ = sys.reset

    # Check the guard direction
    if guard_direction ≠ -1
        @warn "Invalid guard direction for nonholonomic systems"
    end
    
    # The multipliers for both the free and sliding modes
    λ_free = find_multiplier(prob)
    λ_dh   = find_multiplier(prob; sliding=true)

    # Check to make sure we're on the correct side of the guard
    if h(qₖ) < ztol
        # If we're moving in, we should be moving out and count that as a reset
        if dot(∇h(qₖ), M(qₖ) \ pₖ) < 0
            x_new = Δ(xₖ, M, A, ∇h, sys)
            push!(sol.t, tₖ)
            push!(sol.x, x_new)
            if dense_out
                push!(sol.dx, f_λ(x_new[1:n], x_new[n+1:end], λ_free(x_new[1:n], x_new[n+1:end]), 0.0))
            end
            # Update the state
            qₖ, pₖ = x_new[1:n], x_new[n+1:end]
        end
    end

    # Determine whether or not we are on the sliding/post-Zeno regime
    if h(qₖ) < ztol  &&  abs(dot(∇h(qₖ), M(qₖ) \ pₖ)) < ztol
        # To what degree does λ preserve the holonomic constraint?
        function guard_error_nh(Λ)
            F(z, t) = f_λ(z[1:n], z[n+1:end], λ_free(x_new[1:n], x_new[n+1:end]), Λ)
            x_predict, _, _, dt_used, dt_next = take_step(solver, prob, F, vcat(qₖ, pₖ), tₖ, Δt, tol, sol; check=false)
            q_next, p_next = x_predict[1:n], x_predict[n+1:end]
            # (Tangent) constraint violation
            return dot(∇h(q_next), M(q_next) \ p_next)
        end
        # If ∇h(q)̇q>0 with λ=0, then we are escaping the guard (inwards) and are escaping the sliding mode
        if guard_error_nh(0.0) > 0
            F(z, t) = f_λ(z[1:n], z[n+1:end], λ_free(z[1:n], z[n+1:end]), 0.0)
            x_predict, _, _, dt_used, dt_next = take_step(solver, prob, F, vcat(qₖ, pₖ), tₖ, Δt, tol, sol; check=false)
            # Record the derivative
            if dense_out
                push!(sol.dx, F(x_predict, tₖ+dt_used))
            end
        else
            # We have the augmented multiplier
            F2(z, t) = f_λ(z[1:n], z[n+1:end], λ_free(z[1:n], z[n+1:end]), λ_dh(z[1:n], z[n+1:end])[end])
            x_predict, _, _, dt_used, dt_next = take_step(solver, prob, F2, vcat(qₖ, pₖ), tₖ, Δt, tol, sol; check=false)
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
        f(z, t) = f_λ(z[1:n], z[n+1:end], λ_free(z[1:n], z[n+1:end]), 0.0)
        x_predict, eventtrigger, t_root, dt_used, dt_next = take_step(solver, prob, f, vcat(qₖ, pₖ), tₖ, Δt, tol, sol)
        # Was there an impact?
        if eventtrigger
            t_star, x_star = locate_event(event_method, prob, solver, f, vcat(qₖ, pₖ), tₖ, Δt, guard(sys, xₖ), tol, sol, stepper)
            x⁺ = Δ(x_star, M, A, ∇h, sys)
            push!(sol.t, t_star)
            push!(sol.x, x_star)

            push!(sol.event_times, t_star)
            push!(sol.event_indices, length(sol.t))

            push!(sol.t, t_star)
            push!(sol.x, x⁺)

            if dense_out
                push!(sol.dx, f(x_star, t_star))
                push!(sol.dx, f(x⁺, t_star))
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

function solve(prob::prob{S, I, T},
               solver::AbstractODESolver=RK4();
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial = 0.01, max_iter = 10^6, 
               tol = 1e-6, ztol = 1e-3,
               guard_direction = default_guard_direction(prob.sys),
               kwargs...) where {S<:NonholonomicSystem, I, T}
    
    sys = prob.sys
    ∇h = sys.normal
    M(q) = sys.M(q)
    V(q) = sys.V(q)
    A(q) = sys.A(q)

    # Initialize solution struct
    sol = NonholonomicSol(prob)

    # Create vector field for ODE solving
    H(q,p) = 1/2*dot(p, M(q) \ p) + V(q)
    # λ is the unknown multiplier to enforce the constraint ∇h(q)̇q = ∇h(q)M(q)\p = 0
    q_dot(q, p) = ForwardDiff.gradient(p -> H(q,p),p)
    p_dot(q, p, λ, β) = ForwardDiff.gradient(q -> -H(q,p), q) .+ A(q)' * λ .+ β*∇h(q)
    # Combining these vector fields together
    f_λ(q, p, λ, β) = vec([q_dot(q, p); p_dot(q, p, λ, β)])

    _, t_end = prob.tspan           # Extract the terminal time of the problem

    Δt = dt_initial                 # Initialize current time step with user input
    iter = 0                        # Start iteration counter

    # Make sure the initial conditions are compatable with the constraints
    x₀ = sol.x[end]
    n = length(x₀) ÷ 2
    q₀, p₀ = x₀[1:n], x₀[n+1:end]
    
    # Orthogonal projection to distribution
    B_Δ(q) = A(q) * inv(M(q))
    # This ↓ functions don't work...
    # P_Δ(q) = I - B_Δ(q)' * inv(B_Δ(q) * B_Δ(q)') * B_Δ(q) 
    if norm( A(q₀) * (M(q₀) \ p₀) ) > tol
        @warn "Initial conditions do not satisfy the constraint. Projecting to distribution."
        p₀ = p₀ - B_Δ(q₀)' * inv(B_Δ(q₀) * B_Δ(q₀)') * B_Δ(q₀) * p₀
        sol.x[end] = vcat(q₀, p₀)
    end

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
        _, _, Δt, _ = take_step_nonholonomic!(solver, prob, f_λ, Δt, tol, ztol, sol; 
                        dense_out = dense_out, event_method=event_method, guard_direction = guard_direction)
    end

    return sol
end