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
end

# Create the system
function NonholonomicSystem(M, V;
                        A = nothing,
                        guard = nothing,
                        normal = nothing,
                        reset = (x, Mfun, Afun, ∇h, sys::NonholonomicSystem) -> specular_refl(x, Mfun, Afun, ∇h, sys),
                        e = 1)
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
        M, V, A, guard, normal, reset, e
    )
end

# General solution struct for nonholonomic systems
struct NonholonomicSol{T, X, DX, I, E, Z}
    t::T     # Time data
    x::X     # x = (q,p), the state and momentum
    dx::DX   # f(x) derivatve at each state x - only filled out when dense_out = true
    prob::I  # Remember the problem - to aid interpolation
    event::E # Times where an event occurs
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
    ε = -(1+e) * dot(p, Mq \ normal) / denom
    # The 'intermal' multipliers
    λ = -ε * sub_matrix \ constraint_vector
    # The new momentum
    p_new = p + ε*∇hq + Aq'*λ
    return vcat(q, p_new)
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

# Determine the nonholonomic integrator step
# This will essentially be another dispatch of 'take_step'
function take_nh_step(solver, prob, xₖ, tₖ, dt_step, tol, sol; Ah = nothing)
    # Extract out the pertinent information
    sys = prob.sys
    M(q) = sys.M(q)
    V(q) = sys.V(q)
    # We have the option to overwrite the contraint matrix
    if isnothing(Ah)
        A(q) = sys.A(q)
    else
        A(q) = Ah(q)
    end

    # Extract out the state and momentum
    n = length(xₖ) ÷ 2
    qₖ, pₖ = xₖ[1:n], xₖ[n+1:end]
    # Create the Hamiltonian
    H(q, p) = 1/2*dot(p, M(q) \ p) + V(q)
    # Create the vector field
    q_dot(q, p) = ForwardDiff.gradient(p -> H(q,p), p)
    p_dot(q, p, λ) = ForwardDiff.gradient(q -> -H(q,p), q) .+ A(q)'*λ
    # We need to solve for λ such that the constraint is conserved, A(q)v=0
    # A first-order analysis predicts that λ should be
    λ₀ = (A(qₖ) * (M(qₖ)\A(qₖ))) \ (A(qₖ) * (M(qₖ) \ (ForwardDiff.gradient(q -> V(q), qₖ) - 1/dt_step*pₖ)))
    # Determine the constraint error as a function of λ
    function constraint_error(λ)
        f_λ(q, p, λ) = [q_dot(q,p); p_dot(q,p,λ)]
        f(x, t) = f_λ(x[1:div(length(x), 2)], x[(div(length(x),2)+1):end], λ)
        x_next, _, _ = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol; check=false)
        q_next, p_next = x_next[1:n], x[n+1:end]
        return A(q_next) * (M(q_next) \ p_next)
    end
    # Let's be brave and solve this via Newton's method
    jacobian_constraint_error(λ) = ForwardDiff.jacobian(constraint_error, λ)
    for _ ∈ 1:10
        λ₀ = λ₀ - jacobian_constraint_error(λ₀) \ constraint_error(λ₀)
        if norm(constraint_error(λ₀)) < tol
            break
        end
    end
    if norm(constraint_error(λ₀)) > tol
        @warn "Unable to solve for nonholonomic constraint multiplier"
    end
    # We have the multiplier. Let's use that to obtain the next step
    f(x, t) = f_λ(x[1:div(length(x), 2)], x[(div(length(x),2)+1):end], λ₀)
    x_next, _, _ = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol; check=false)
    return x_next
end

# Find the time of a crossing
function locate_event_nonholonomic(solver, prob, xₖ, tₖ, Δt, tol, sol, h)
    # We will implement a linear finder that triggered 'crossed_guard_nonholonomic'
    t₀, t₁ = 0.0, Δt
    # The event function as a function of time
    E(t) = h(take_nh_step(solver, prob, xₖ, tₖ, t, tol, sol))

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

# The solve function
function solve(prob::prob{S, I, T};
               solver::AbstractODESolver=RK4(),
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial = 0.01, max_iter = 10^5,
               tol = 1e-6, ztol = 1e-4,
               kwargs...) where {S<:NonholonomicSystem, I, T}
    
    # Unpack pertinent information
    sys = prob.sys
    h = sys.guard
    ∇h = sys.normal
    Δ = sys.reset
    M(q) = sys.M(q)
    V(q) = sys.V(q)
    A(q) = sys.A(q)

    # Initialize the solution struct
    sol = NonholonomicSol(prob)

    # Time bounds
    t_start, t_end = prob.tspan
    Δt = dt_initial
    iter = 0

    # Make sure the initial conditions are compatable with the constraints
    x₀ = sol.x[end]
    n = length(x₀) ÷ 2
    q₀, p₀ = x₀[1:n], x₀[n+1,end]
    # Orthogonal projection to distribution
    B_Δ(q) = A(q) * inv(M(q))
    P_Δ(q) = I - B_Δ(q)' * ( (B_Δ(q) * B_Δ(q)') \ B_Δ(q) )
    if norm( A(q₀) * (M(q₀) \ p₀) ) > tol
        @warn "Initial conditions do not satisfy the constraint. Projecting to distribution."
        p₀ = P_Δ(q₀) * p₀
        sol.x[end] = vcat(q₀, p₀)
    end

    # Run until end of specified time span
    while sol.t[end] < t_end
        # Safties
        iter += 1
        if iter > max_iter
            @warn "Maximum Iteration Count ($max_iter) exceeded."
            break
        end
        # Terminate if the remaining time is below machine precision
        if t_end - sol.t[end] <= eps(t_end)
            break
        end

        # Truncate the time step if we overshoot the final time
        dt_step = (sol.t[end]+Δt > t_end) ? (t_end-sol.t[end]) : Δt

        # Perform the step
        tₖ, xₖ = sol.t[end], sol.x[end]
        qₖ, pₖ = xₖ[1:div(length(xₖ), 2)], xₖ[div(length(xₖ),2)+1:end]
        
        # We want to allow for Zeno/post Zeno solutions
        if (h(qₖ) < ztol && -ztol < dot(∇h(qₖ), M(qₖ) \ pₖ) < 0) || h(qₖ) < -ztol
            # The sliding mode step, is the (holonomic) constraint satisfied?
            # We can simply append the constraint ∇h⋅v = 0 to the matrix A(q)
            A_h(q) = vcat(A(q), ∇h(q))
            # And take the nonholonomic step with this new matrix
            # But! we only want to apply this constraint unilaterally
            # Not applying this constraint first
            x_next = take_nh_step(solver, prob, xₖ, tₖ, dt_step, tol, sol)
            q_next, p_next = x_next[1:div(length(x_next), 2)], x_next[(div(length(x_next), 2) + 1):end]
            if dot(∇h(q_next), M(q_next) \ p_next) < 0 # <- So we actually do need a reaction force
                x_next = take_nh_step(solver, prob, xₖ, tₖ, dt_step, tol, sol; Ah = A_h)
            end
        else
            # The 'usual' hybrid step
            x_proposed = take_nh_step(solver, prob, xₖ, tₖ, dt_step, tol, sol)
            # Is there a crossing detected?
            if crossed_guard_nonholonomic(h(qₖ), h(x_proposed[1:div(length(x_proposed), 2)]), 0.0, Δt)[1]
                Δt_found = locate_event_nonholonomic(solver, prob, xₖ, tₖ, dt_step, tol, sol, h)
                x_impact = take_nh_step(solver, prob, xₖ, tₖ, Δt_found, tol, sol)
                x_post = Δ(x_impact, M, A, ∇h, sys)
                x_next = take_nh_step(solver, prob, x_post, tₖ, dt_step-Δt_found, tol, sol)
            else
                x_next = x_proposed
            end
        end
        # Record
        push!(sol.t, dt_step+tₖ)
        push!(sol.x, x_next)
    end

    return sol
end