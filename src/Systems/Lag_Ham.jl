
# INTERNAL
# Default reset map: spectral reflection with coefficient of restitution. See "Is There Life After Zeno? paper
######
### WC: You should include a complete reference.
######
function specular_refl(x, M, dh; e=1.0)

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
######
### WC: This requires 'using' LinearAlgebra. Have you settled on 'using' rather than 'import'?
######

# General solution struct for Lagrangian and Hamiltonian systems
struct LHSol{T, X, T_e1, T_z}
    T::T        # Time data
    X::X        # Position data
    T_e1::T_e1  # Time of first event
    T_z::T_z    # Time of Zeno
end
######
### WC: There can be multiple Zeno events. We'll discuss this more in the future.
### As it is currently written, you halt when Zeno occurs, correct?
######

# Does this need to be a general constructor or can I go straight to initializing the solution state?
function LHSol(prob)
    return LHSol([prob.tspan[1]], [prob.init], Float64[], Float64[])
end

################################
################################
## Lagrangian dynamics

######
### WC: These Lagrangian systems are not of the form we discussed. 
### LagSys(M,V,h,e) should suffice
### The reset map does not need to be given, so I am confused about how to are constructing everything here.
######

# General Lagrangian system
struct LagSys{L,G,N,R,E}
    L::L          # Lagrangian
    guard::G      # guard/event function
    normal::N     # Normal to the guard, ΔG ### <- ∇G
    reset::R      # reset map
    e::E          # coefficient of restitution
end

# EXTERNAL
# Make the the guard, reset map, and coefficient of restitution optional; default to fully elastic specular reflection
"""
Lagrangian System
 - L

"""
function LagSys(L;
            guard = nothing,
            normal = nothing,
            reset = (x, M, dh, e) -> specular_refl(x, M, dh; e),
            e = 1.0)

    if isnothing(guard) &&  !isnothing(normal)
        error("Normal to guard was provided, but a guard was not")
    end

    # Default to auto diff
    if !isnothing(guard) && isnothing(normal)
        normal= q -> ForwardDiff.gradient(guard, q)
    end

    return LagSys(L, guard, normal, reset, e)
end

# INTERNAL
# Equations of motion from Euler-Lagrange equations using ForwardDiff

    # dL/dq is the force
    function lagrangian_force(L::Function, q::AbstractVector, qdot::AbstractVector)
        ForwardDiff.gradient(q -> L(q, qdot), q)
    end

    # dL/dqdot is the momentum
    function lagrangian_momentum(L::Function, q::AbstractVector, qdot::AbstractVector)
        ForwardDiff.gradient(qdot -> L(q, qdot), qdot)
    end

    # Mass matrix
    function lagrangian_mass_matrix(L::Function, q::AbstractVector, qdot::AbstractVector)
        ForwardDiff.hessian(qdot -> L(q, qdot), qdot)
    end

    # Coriolis matrix
    function lagrangian_coriolis(L::Function, q::AbstractVector, qdot::AbstractVector)
        jac = ForwardDiff.jacobian(q -> ForwardDiff.gradient(qdot -> L(q, qdot), qdot), q)
        jac * qdot
    end

    # Inertial acceleration
    function lagrangian_acceleration(L::Function, q::AbstractVector, qdot::AbstractVector)
        M = lagrangian_mass_matrix(L, q, qdot)
        C = lagrangian_coriolis(L, q, qdot)
        F = lagrangian_force(L, q, qdot)
        M \ (F - C)
    end

    # INTERNAL
    # Find the complete vector field from a Lagrangian using automatic differentiation
    function vec_field(sys::LagSys, x::AbstractVector, t)

        L = sys.L

        # Integer division
        n = length(x) ÷ 2

        q = Vector(x[1:n])
        qdot = Vector(x[n+1:end])

        # Wrap L(x,t) into internal L(q,qdot), so I don't have to rewrite everything internally for now
        Lsplit(q, qdot) = L(vcat(q, qdot), t)

        # Solve for acceleration
        qddot = lagrangian_acceleration(Lsplit, q, qdot)

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

 # Some way to find the gradient without ForwardDiff (untested)
# function hamiltonian_gradient(Hsys, x)
    
#     H = Hsys.H
#     h = Hsys.h

#     # initialize place to store the output
#     grad = similar(x)

#     # This feels like a bad way of doing this but I think it's correct
#     for i in eachindex(x)

#         xp = copy(x)
#         xm = copy(x)

#         # perturb the point forward and backward
#         xp[i] += h
#         xm[i] -= h

#         # Central finite diff approximation
#         grad[i] = (H(xp) - H(xm)) / (2h)
#     end

#     return grad
# end

# Could add another dispatch method to better deal with interpolating data style hams


################################
################################
## Solver

# INTERNAL
# Zero on guard. I guess this works?
######
### WC: This is just evaluation of the 'guard function'. The guard is the zero level set of this function.
######
function guard(sys::Union{LagSys, HamSys}, x)
    if isnothing(sys.guard)
        return nothing
    else
        return sys.guard(x)
    end
end

# EXTERNAL
# solver specifically for Lagrangian and Hamiltonian systems
######
### WC: Again with the tolerences....
### Also, there is no 'M'?
######
######
### WC: Why don't you transfrom this problem type into a general hybrid system and pass to that solver? What is gained by writing another dispatch?
######
function solve(prob::prob{S, I, T}, solver; event_method::AbstractEventLocator=BisectionLocator(), dt_initial = 0.01, max_iter = 10^6, tol = 1e-12, kwargs...) where {S<:Union{LagSys, HamSys}, I, T}
    
    sys = prob.sys
    G = sys.guard
    dh = sys.normal

    # Create vector field for ODE solving
    f(x,t) = vec_field(sys, x, t)
    # Need mass matrix
    M = 
    #Define reset map    
    reset(x, M, dh; e=1.0) = sys.R(x, M, dh; e=1.0)

    sol = LHSol(prob)     # See line 57

    t_start, t_end = prob.tspan     # Extract start and end times for bounds

    Δt = dt_initial                 #Initialize current time step with user input
    iter = 0                        #Start iteration counter

    # Run sim until end of specified time span
    while sol.T[end] < t_end 
        
    # Safties
        # Stop if we hit the iteration limit to avoid memory doomsday
        iter += 1
        if iter > max_iter 
            @warn "Maximum Iteration Count ($max_iter) exceeded."
            break
        end

        # Terminate if the remaining time is below machine precision
        if t_end - sol.T[end] <= eps(t_end)
            break
        end

        #Truncate time step if we overshoot the final sim time
        dt_step = (sol.T[end] + Δt > t_end) ? (t_end - sol.T[end]) : Δt

    # Actually solve now

        xₖ = sol.X[end] #Retrieve current state at start of step
        tₖ = sol.T[end] #Retrieve current time at start of step

        #ATTEMPT CONTINUOUS STEP
        #Dispatch calls the specific math for the chosen solver. Returns pre state and boolean flag for if guard was crossed. 
        x_predict, event_triggered, h_now = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol)
        t_next = tₖ + dt_step

   if event_triggered

        #ZENO DETECTION
        veloapprox = (xₖ .- sol.X[end - 1]) / dt_step
        if abs(G(xₖ)) < tol && abs(dot(dh(xₖ), veloapprox))
            
            @warn "Zeno condition detected" xₖ

            # # Reinitialize on guard
            # x₀ = project_to_guard(xₖ, G, dh)

            v = veloapprox
            # Compute normal
            n = dh(x₀)

            # Project velocity into tangent space
            vₜ = v - (dot(n, v) / dot(n, n)) * n

            # # Replace current state with constrained state
            # xₖ = x₀

            # Continue with reduced timestep since weird stuff is happening
            dt_step = dt_step * 0.5

        # Not Zeno, normal reset map
        else
            #Pinpoint the exact impact time and state using the chosen locator strategy 
            t_star, x_star = locate_event(event_method, sys, solver, f, xₖ, tₖ, Δt, h_now, tol, sol)            

            # Reset map - defaults to specular reflection
            x⁺ = reset(xₖ, M, dh; e=1.0)

            #PLOTTING ARCH
            #We explicitly push both the pre-impact state x_star and post-impact state x⁺ to the same timestamp t_star.
            #This allows our post-processing functions to give NaN values between the points preventing lines between plots when we do that.
            #Also old code, can be gotten rid of/ altered?
            push!(sol.T, t_star, t_star)
            push!(sol.X, x_star, x⁺)

            #Lock in post-impact state to continue the loop
            xₖ = x⁺ 
            tₖ = t_star
        
            #Shrink step size for next step to avoid overshooting and missing possible events. 
            dt = dt_min
        end

    #IF NOT EVENT go to next step Log it then Loop.
        else
            tₖ += dt_step
            xₖ = x_predict
            push!(sol.T, tₖ)
            push!(sol.X, xₖ)

            dt = min(dt_step, dt_initial)
        end
    end

    return sol
end
