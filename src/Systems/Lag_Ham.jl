
# INTERNAL
# Default reset map: specular reflection with coefficient of restitution. 
# Ames, Aaron & Zheng, Haiyang & Gregg, Robert & Sastry, Shankar. (2006). Is there life after Zeno? Taking executions past the breaking (Zeno) point. 2006. 6 pp.. 10.1109/ACC.2006.1656623

function specular_refl(x, M, dh, sys)

    e = sys.e
    n = length(x) ÷ 2

    q = x[1:n]          # positions
    v = x[n+1:end]      # velocities

    Mq = M(q)           # Mass matrix

    # Constraint normal (row -> column)
    normal = vec(dh(q))

    # Denominator
    denom = normal' * (Mq \ normal)

    # Full equation
    P = I - (1 + e) * ((Mq \ (normal * normal')) / denom)
    vnew = P * v

    return vcat(q, vnew)
end

# General solution struct for Lagrangian and Hamiltonian systems
struct LHSol{T, X, DX, I, T_e1, T_z}
    t::T        # Time data
    x::X        # Position data
    dx::DX      # f(x) Derivative at each state x - only filled when dense_out = true
    prob::I     # Remember the problem - to aid interpolation
    T_e1::T_e1  # Time of first event
    T_z::T_z    # Time of Zeno pts
end
######
### WC: There can be multiple Zeno events. We'll discuss this more in the future.
### As it is currently written, you halt when Zeno occurs, correct?
# CK: Currently, yes
######

# Does this need to be a general constructor or can I go straight to initializing the solution state?
function LHSol(prob)
    return LHSol([prob.tspan[1]], [prob.init], Vector{Vector{Float64}}(), prob, Float64[], Float64[])
end

################################
################################
## Lagrangian dynamics

# General Lagrangian system
struct LagSys{M,V,G,N,R,E} <: AbstractHybridSystem
    M::M          # Mass matrix function M(q)
    V::V          # Potential energy function V(q)
    guard::G      # guard/event function
    normal::N     # Normal to the guard, ΔG ### <- ∇G
    reset::R      # reset map
    e::E          # coefficient of restitution
end

# EXTERNAL
# Make the the guard, reset map, and coefficient of restitution optional; default to fully elastic specular reflection
"""
Lagrangian System
 - M(q): mass matrix
 - V(q): potential energy
"""
function LagSys(M, V;
                guard = nothing,
                normal = nothing,
                reset = (x, Mfun, dh, sys::LagSys) -> specular_refl(x, Mfun, dh, sys),
                e = 1.0)

    if isnothing(guard) && !isnothing(normal)
        error("Normal to guard was provided, but a guard was not")
    end

    if !isnothing(guard) && isnothing(normal)
        normal = q -> ForwardDiff.gradient(guard, q)
    end

    return LagSys{typeof(M), typeof(V), typeof(guard), typeof(normal), typeof(reset), typeof(e)}(
        M, V, guard, normal, reset, e
    )
end

# Equations of motion from Euler-Lagrange expressed via M(q) and V(q)
function lagrangian_force(M, V, q::AbstractVector, qdot::AbstractVector)
    ForwardDiff.gradient(q_ -> 0.5 * dot(qdot, M(q_) * qdot) - V(q_), q)
end

function lagrangian_momentum(M, q::AbstractVector, qdot::AbstractVector)
    M(q) * qdot
end

function lagrangian_mass_matrix(M, q::AbstractVector)
    M(q)
end

function lagrangian_coriolis(M, q::AbstractVector, qdot::AbstractVector)
    jac = ForwardDiff.jacobian(q_ -> M(q_) * qdot, q)
    return jac * qdot
end

function lagrangian_acceleration(M, V, q::AbstractVector, qdot::AbstractVector)
    Mq = lagrangian_mass_matrix(M, q)
    C = lagrangian_coriolis(M, q, qdot)
    F = lagrangian_force(M, V, q, qdot)
    return Mq \ (F - C)
end

function vec_field(sys::LagSys, x::AbstractVector, t)
    n = length(x) ÷ 2
    q = x[1:n]
    qdot = x[n+1:end]
    qddot = lagrangian_acceleration(sys.M, sys.V, q, qdot)
    return vcat(qdot, qddot)
end


################################
################################

# IGNORE FOR NOW
# Hamiltonian dynamics

    # Define a general Hamiltonian problem
    struct HamSys{H,G,R,E}
        H::H          # Hamiltonian
        guard::G          # guard/event function
        reset::R      # reset map
        e::E          # coefficient of restitution 
    end

# Make the the guard, reset map, and coefficient of restitution optional; default to fully elastic spectral reflection
function HamSys(H; guard=nothing, reset = (x,e) -> specularl_refl(x,e), e = 1.0)
    HamSys(H, guard, reset, e)
end

    # Find the vector field from Hamilton's equations
    function vec_field(sys::HamSys, x::AbstractVector, t)

        n = length(x) ÷ 2

        # Split up state vector
        q = x[1:n]
        p = x[n+1:end]

        # Gradient of H with respect to state vector x = [q; p]
        gradH = hamiltonian_gradient(sys, x)

        dqdt = gradH[n+1:end]      # ∂H/∂p
        dpdt = -gradH[1:n]         # -∂H/∂q

        # Recombine derivatives into single vecotr field
        return vcat(dqdt, dpdt)
    end

# Find the gradient of Hamiltonians that can be differentiated using ForwardDiff
function hamiltonian_gradient(HamSys, x)
    ForwardDiff.gradient(HamSys.H, x)
end


################################
################################
## Solver

# INTERNAL
# Evaluates guard function if present. Required for event locating.
function guard(sys::Union{LagSys, HamSys}, x)
    if isnothing(sys.guard)
        return nothing
    else
        return sys.guard(x)
    end
end

# EXTERNAL
# Solver specifically for Lagrangian and Hamiltonian systems
######
### WC: Why don't you transfrom this problem type into a general hybrid system and pass to that solver? What is gained by writing another dispatch?
######
function solve(prob::prob{S, I, T};
               solver::AbstractODESolver=RK45(),
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial = 0.01, max_iter = 10^5, tol = 1e-6, kwargs...) where {S<:Union{LagSys, HamSys}, I, T}
    
    sys = prob.sys
    G = sys.guard
    dh = sys.normal
    reset = sys.reset

    # Create vector field for ODE solving
    f(x,t) = vec_field(sys, x, t)

    # Mass matrix function for event resets
    M(x) = sys.M(x)

    # Initialize solution struct
    sol = LHSol(prob)

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

        # Stagnation error
        if length(sol.t) > 6
            Δt = sol.t[end] - sol.t[end-5]
            Δx = norm(sol.x[end] - sol.x[end-5])

            if Δt < tol && Δx < tol
                error("Stagnation detected: no meaningful time/state progression over 5 steps")
            end
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

        #ATTEMPT CONTINUOUS STEP
        #Dispatch calls the specific math for the chosen solver. Returns pre state and boolean flag for if guard was crossed. 
        x_predict, event_triggered, h_now = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol)
        t_next = tₖ + dt_step

   if event_triggered

        #ZENO DETECTION
        veloapprox = (xₖ .- sol.x[end - 1]) / dt_step
        if abs(G(xₖ)) < tol && abs(dot(dh(xₖ), veloapprox))
            
            # Halt for now, the following few lines approximate continuing after zeno
            error("Zeno condition detected", xₖ)

            # v = veloapprox

            # # Compute normal
            # n = dh(x₀)

            # # Project velocity into tangent space
            # vₜ = v - (dot(n, v) / dot(n, n)) * n

            # # Store solution and continue
            # # NEED TO ADD PUSHES 
            # push!(sol.dx, NaN) # idk if this is how I want this to work but there shouldn't be interpolation between discrete jumps

        # Not Zeno, normal reset map
        else
            #Pinpoint the exact impact time and state using the chosen locator strategy 
            t_star, x_star = locate_event(event_method, prob, solver, f, xₖ, tₖ, Δt, h_now, tol, sol)            

            # Reset map - defaults to spectral reflection
            x⁺ = reset(xₖ, M, dh, sys)

            # Store solution
            push!(sol.t, t_star, t_star)
            push!(sol.x, x_star, x⁺)
            if dense_out
                push!(sol.dx, similar(sol.x[end]) .= NaN)  # Don't interpolate across discrete jumps
            end

            # Lock in post-impact state to continue the loop
            xₖ = x⁺ 
            tₖ = t_star
        
            # Shrink step size for next step to avoid overshooting and missing possible events
            dt_step = dt_step * 0.5
        end

    #IF NOT EVENT go to next step Log it then Loop.
        else

            # this is technically redundant, but it would take a lot of rewriting to get f out of take_step
            if dense_out
                push!(sol.dx, f(xₖ, tₖ))
            end

            tₖ += dt_step
            xₖ = x_predict
            push!(sol.t, tₖ)
            push!(sol.x, xₖ)

            dt_step = min(dt_step, dt_initial)
        end
    end

    return sol
end
