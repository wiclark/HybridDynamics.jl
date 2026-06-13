
# General mechanical system
struct MechanicalSystem{M,V,G,N,R,E} <: AbstractHybridSystem
    M::M          # Mass matrix function M(q)
    V::V          # Potential energy function V(q)
    guard::G      # guard/event function
    normal::N     # Normal to the guard, ∇G
    reset::R      # reset map
    e::E          # coefficient of restitution
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
                reset = (x, Mfun, dh, sys::LagSys) -> specular_refl(x, Mfun, dh, sys),
                e = 1.0)

    if isnothing(guard) && !isnothing(normal)
        error("Normal to guard was provided, but a guard was not")
    end

    if !isnothing(guard) && isnothing(normal)
        normal = q -> ForwardDiff.gradient(guard, q)
    end

    return MechanicalSystem{typeof(M), typeof(V), typeof(guard), typeof(normal), typeof(reset), typeof(e)}(
        M, V, guard, normal, reset, e
    )
end

# General solution struct for mechanical systems
struct MechanicalSol{T, Q, V, P, DX, I, E, Z}
    t::T        # Time data
    q::Q        # Position data
    v::V        # Velocity
    p::P        # Momentum
    dx::DX      # f(x) Derivative at each state x - only filled when dense_out = true
    prob::I     # Remember the problem - to aid interpolation
    event::E    # Times where an event has occurred 
    zeno::Z     # Times of Zeno points
end

# Function to initialize solution struct
function MechanicalSol(prob)
    return MechanicalSol(prob.tspan[1],
        [prob.init],
        Vector{Vector{Float64}}(),      # Velocity
        Vector{Vector{Float64}}(),      # Momentum
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

function solve(prob::prob{S, I, T};
               solver::AbstractODESolver=RK4(),
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial = 0.01, max_iter = 10^5, 
               tol = 1e-6, ztol = 1e-4,
               kwargs...) where {S<:MechanicalSystem, I, T}
    
    sys = prob.sys
    G = sys.guard
    dh = sys.normal
    reset = sys.reset
    M(q) = sys.M(q)

    # Initialize solution struct
    sol = MechanicalSol(prob)

    # Create vector field for ODE solving
    H(q,p) = 1/2*dot(p, M(q) \ p) + V(q)
    # Requires ForwardDiff (I'll add the option to user define these later -CK)
    # λ is the unknown multiplier to enforce the constraint ∇h(q)̇q = ∇h(q)M(q)\p = 0
    q_dot(q, p) = ForwardDiff.gradient(p -> H(q,p),p)
    p_dot(q, p, λ) = ForwardDiff.gradient(q -> -H(q,p), q) + λ*∇h(q)
    # Combining these vector fields together
    f_λ(q, p, λ) = [q_dot(q[1:n], q[n+1:end]); p_dot(z[1:n], z[n+1:end], λ)]
    f(x,t) = f_λ(q, p, λ) # Maybe? I think x needs to contain q and p in spirit before the previous line

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

        tₖ = sol.t[end] # Current time at start of step
        qₖ = sol.q[end] # Current state at start of step
        pₖ = sol.p[end] # Current momentum at start of step

        # Recall that we want h(z)≈0 and -ε<dh(q)̇q<0
        # First, is this a (post) Zeno state?
        if h(qₖ) < ztol  &&  -ztol < dot(∇h(q), M(q) \ p) < 0
            # Does λ preserve the constraint?
            function guard_error(λ)
                F(z, t) = f_λ(z, λ)
                #### Make sure that the syntax below works ####
                q_next = solver(F, z, Δt) 
                q_next, p_next = q_next[1:n], q_next[n+1:end]
                # Constraint?
                return dot(∇h(q_next), M(q_next) \ p_next)
            end

            # If ∇h(q)̇q>0 (with λ=0), then we are moving to the interior of the state-space and the constraint need not be applied
            if guard_error(0) > 0
                F(z, t) = f_λ(z, 0.0)
                q_next = solver(F, z, Δt) #### <--- This should still check for impacts!!!!
            else
                # This is the fun part; we need to actually solve for λ
                # We will solve for λ via the method of false position. This needs two initial guesses for λ
                # λ₀ = 0, because why not.
                # λ₁ = the answer predicted by applying symplectic Euler. 
                λ₀ = 0.0
                λ₁ = dot(∇h(q), M(q)\(-p_dot(q,p,0)-p)) / dot(∇h(q), M(q)\∇h(q))
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
                f_constrained(z, t) = f_λ(z, λ₁)
                q_next = solver(F, z, Δt)
            end            
            Δt_found = Δt

            # No sliding occurs, 'normal' hybrid situation
            else
                # Propose a step assuming no impacts
                z_proposed = solver(f, z, Δt)

                # Is there a crossing detected?
                if crossed_guard(h(z[1:n]), h(z_proposed[1:n]), 0.0, Δt)[1]
                    Δt_found = locate_event(f, z, Δt, h)
                    z_impact = solver(f, z, Δt_found)
                    q_next = Δ(z_impact)
                else
                    q_next = q_proposed
                    Δt_found = Δt
                end
            end

        # Record
        push!(sol.t, Δt_found+tₖ)
        push!(sol.q, q_next)
    end

    return sol
end


##############
##############
### Original WC code
# As we are only interested in Lagrangian/Hamiltonian systems of mechanical type, it would be completely redundant to write up both cases
# I will nominally write everything in Hamiltonian form, but the two versions can be translated back and forth

# The data for the problem will be the mass matrix M(q) and the potential energy function V(q). 
# Everything will be assumed to be time-invariant (the ODE will be autonomous). Down the line, I'd like to incorporate time-dependent systems.

# The function written here will perform a single step. It will have to be iterated to make the full solver.

# The step will consist of two parts:
#  1) Check if we're sliding along the guard (Zeno / post-Zeno)
#  2) Apply the usual dynamics with impacts

# This requires a 'solver' of the form solver(f::Vector field, z::State, Δt::Time step)

# As things below are written, h(q) is a function of q, not z = [q,p]

"""
A function that computes the dynamics for sliding along the guard (to be used post-Zeno)
    This function will allow the particle to leave the surface.
    It is assumed that Q = {h(q)≥0} and thus an impact occurs when h(q)=0 and dh(̇q)<0
        Inputs:
            M::Matrix function (mass matrix)
            V::Potential function (scalar valued)
            h::Guard function (scalar valued)
            ∇h::Its gradient (technically, its differential, but as a vector)
            solver::The solver type, e.g. RK4
            x::The current state
            Δt::The proposed time step
            *There is no time input as the system is assumed to be autonomous*
        Outputs:
            x_out::The computed output state
            Δt_out::The actual time step
"""
function sliding_dynamics_step(M, V, h, ∇h, solver, x, Δt)
    # Creating the Hamiltonian and vector field
    # It would most probably be better to precompute this vector field and pass it to this function.
    n = length(x) ÷ 2
    q = x[1:n]
    p = x[n+1:end]
    H(q,p) = 1/2*dot(p, M(q) \ p) + V(q)
    # The dynamics - requires ForwardDiff
    # λ is the unknown multiplier to enforce the constraint ∇h(q)̇q = ∇h(q)M(q)\p = 0
    q_dot(q, p) = ForwardDiff.gradient(p -> H(q,p),p)
    p_dot(q, p, λ) = ForwardDiff.gradient(q -> -H(q,p), q) + λ*∇h(q)
    # Combining these vector fields together
    f_λ(z, λ) = [q_dot(z[1:n], z[n+1:end]); p_dot(z[1:n], z[n+1:end], λ)]
    # Does λ preserve the constraint?
    function guard_error(λ)
        F(z, t) = f_λ(z, λ)
        #### Make sure that the syntax below works ####
        q_next = solver(F, z, Δt) 
        q_next, p_next = q_next[1:n], q_next[n+1:end]
        # Constraint?
        return dot(∇h(q_next), M(q_next) \ p_next)
    end

    # If ∇h(q)̇q>0 (with λ=0), then we are moving to the interior of the state-space and the constraint need not be applied
    if guard_error(0) > 0
        F(z, t) = f_λ(z, 0.0)
        q_next = solver(F, z, Δt) #### <--- This should still check for impacts!!!!
    else
        # This is the fun part; we need to actually solve for λ
        # We will solve for λ via the method of false position. This needs two initial guesses for λ
        # λ₀ = 0, because why not.
        # λ₁ = the answer predicted by applying symplectic Euler. 
        λ₀ = 0.0
        λ₁ = dot(∇h(q), M(q)\(-p_dot(q,p,0)-p)) / dot(∇h(q), M(q)\∇h(q))
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
        f_constrained(z, t) = f_λ(z, λ₁)
        q_next = solver(F, z, Δt)
    end
    return q_next
end

"""
A function which determines whether or not the guard is crossed. 
    This is going to be lazy and only utilize the linear locater.
    Moreover, we only want to record one-sided impacts. That is, h(q)=0 and dh(q)̇q<0
        Inputs:
            h_now::The current value of h(q_now)
            h_next::The next predicted value of h(q_next)
            t_now::The current time
            t_next::The next time
        Outputs:
            A boolean::Whether or not there is a hit
            A scalar::The predicted hit time (NaN if false)
"""
function crossed_guard(h_now, h_next, t_now, t_next)
    # Two point crossing condition
    if h_now > 0 && h_next < 0
        t_root = t_now - h_now * (t_next-t_now) / (h_next-h_now)
        return true, t_root
    else
        return false, NaN
    end
end

"""
A function that locates the actual time of impact if one is predicted by 'crossed_guard'
    This method utilized a line search along various step sizes. 
    *This requires methods with variable step size.*
        Inputs:
            f::The vector field
            z::The current state
            Δt::The proposed step size
            h::The event function
        Outputs:
            The predicted time t^* ∈ (0.0, Δt)
"""
function locate_event(f, z, Δt, h)
    # We will implement a linear finder as that is what triggered 'crossed_guard'
    t₀, t₁ = 0.0, Δt
    # The event function as a function of time
    E(t) = h(solver(f,z,t))
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

"""
A function that computed a single step with a Zeno test
    This step will either compute a sliding & constrained mode, a continuous mode, or an impact mode
    Currently this is lacking as it assumes that only a single impact occurs. 
    The solver will halt at a premature time if an impact occurs.
        Inputs: 
            M::Mass matrix
            V::Potential energy
            f::Vector field
            z::Current state
            Δt::Proposed time step
            h::Event function
            ∇h::Gradient of event function
            solver::Choise of stepper
            Δ::The impact map
        Outputs:
            q_next::The computed state
            Δt_found::The actually used time step
"""
function allowing_zeno_step(M, V, f, z, Δt, h, ∇h, solver, Δ)
    # First off, is this a (post) Zeno state?
    # Recall that we want h(z)≈0 and -ε<dh(q)̇q<0
    n = length(x) ÷ 2
    q = x[1:n]
    p = x[n+1:end]
    if h(q)<1e-3 && -1e-3<dot(∇h(q),M(q) \ p)<0
        q_next =  sliding_dynamics_step(M, V, h, ∇h, solver, z, Δt)
        Δt_found = Δt
    else # <- No sliding occurs, 'normal' hybrid situation
        # Propose a step assuming no impacts
        z_proposed = solver(f, z, Δt)
        # Is there a crossing detected?
        if crossed_guard(h(z[1:n]), h(z_proposed[1:n]), 0.0, Δt)[1]
            Δt_found = locate_event(f, z, Δt, h)
            z_impact = solver(f, z, Δt_found)
            q_next = Δ(z_impact)
        else
            q_next = z_proposed
            Δt_found = Δt
        end
    end
    # Return the state and time step
    return q_next, Δt_found
end

"""
A function that iterates over 'allowing_zeno_step' to compute the full trajectory
    A variable step solver is required
        Inputs:
            M::Mass matrix
            V::Potential energy
            f::Vector field
            z0::Initial condition
            t_f::Final time (intitial time is t0=0.0)
            h::Event function
            ∇h::Gradient of event function
            solver::Proposed ODE stepper
            Δt::Proposed time step
            Δ::Reset map
        Outputs:
"""
function solve_allowing_zeno(M, V, f, z0, t_f, h, ∇h, solver, Δt, Δ)
    # Initialize the data
    solt = [0.0]
    solz = [z0]
    # Loop over the steps
    while solt[end] < tf - 1e-3
        # How far do we want to go?
        dt = minimum([Δt, t_f-solt[end]])
        # Update the state and actually see how far we go
        z_now = solz[end]
        q_next, Δt_found = allowing_zeno_step(M, V, f, z_now, dt, h, ∇h, solver, Δ)
        # Record
        push!(solt, Δt_found+solt[end])
        push!(solz, q_next)
    end
    return solt, solz
end