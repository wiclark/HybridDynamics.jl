# I want to understand hybrid Lagrangian systems (Hamiltonian will follow)

using ForwardDiff

struct HybridLagrangianSystem
    L::Function
    h::Function
    e::Function
end

# I want the following:
# 1. Generate trajectories (you can call forward diff)
# 2. Implement an intelligent way to perform event detection
# 3. Determine whether or not the system is Zeno and state when/where that occurs.


# # Overarching type, might not be necessary
# abstract type AbstractDynamics end


# Allow options for how the solver arrives at the equations of motion (EOM)
abstract type Backend end
struct AutoForwardDiff <: Backend end   # Automatic differentiation using ForwardDiff
struct AutoFiniteDiff  <: Backend end   # Finite differences
struct ManualEOM       <: Backend end   # Manually provide sumthin?

# General problem structure
struct Prob{F, I, T}
    sys::F
    init::I
    tspan::T
end

################################
################################
"""
Lagrangian dynamics
"""

# General Lagrangian system
struct LagSys{L,G,R,E,B}
    L::L          # Lagrangian
    guard::G      # guard/event function
    reset::R      # reset map
    e::E          # coefficient of restitution
    B::B          # backend use to arrive at EOM
end

# Make the the guard, reset map, and coefficient of restitution optional; default to fully elastic spectral reflection
function LagSys(L;
            guard=nothing,
            reset = (x,e) -> spectral_refl(x,e),
            e = 1.0,
            B = AutoForwardDiff())

    LagSys(L, guard, reset, e, B)
end

# Solve a Lagrangian problem
function solve(prob::Prob{<:LagSys}, solver; kwargs...) # solver specifically for LagProb struct, kwargs for step size or tolerances (optional / solver dependent)
    
    fun = prob.sys
    F(x, t) = lagrangian_vec_field(fun.L, x, t)

    sys = prob.sys

    # If no guard is provided, return the normal continuous solve
    if isnothing(sys.guard)
        return solver(F, prob.init, prob.tspan; kwargs...)
    end

    # If a guard is provided, return the hybrid specific solve function
    return hybrid_solve(F, prob, solver; kwargs...)
end

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

    # Find the complete vector field from a Lagrangian using automatic differentiation
    function lagrangian_vec_field(L, x::AbstractVector, t)

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

# *Still need to make this go somewhere* Option to input M, C, and F directly
struct ManualLag{M, C, F}
    M::M    # Mass matrix, could be noninvertable
    C::C    # Coriolis matrix
    F::F    # Force
end


### Hybrid Lagrangian systems:

# Solve dispatch specific to hybrid systems
function hybrid_solve(prob, solver; kwargs...)

    sys = prob.sys

    x = prob.init
    t0, tf = prob.tspan

    

    return sol
end
"""
This was my first thought. I want the event detection to be within the solver though (I think)
    # function solve(prob::LagProb{<:HybridLagSys}, solver; kwargs...)

    #     sys = prob.sys

    #     x = prob.init
    #     t0, tf = prob.tspan

    #     trajectory = []

    #     while t0 < tf

    #         # Continuous dynamics
    #         F(x,t) = lagrangian_vec_field(sys.L, x, t)

    #         # Integrate until event
    #         sol = solver(F, x, (t0, tf); event = sys.h, kwargs...)

    #         push!(trajectory, sol)

    #         # No event occurred
    #         if terminal(sol)
    #             break
    #         end

    #         # Apply reset map after impact
    #         x = sys.reset(sol.x[end])

    #         # Restart after event time
    #         t0 = sol.t[end]
    #     end

    #     return trajectory
    # end
"""

# Default reset map: spectral reflection with coefficient of restitution
function spectral_refl(x, e)

    n = length(x) ÷ 2

    q = x[1:n]          # positions
    qdot = x[n+1:end]   # velocities

    # positions unchanged, velocities reflected




    return vcat(q, qdot)
end


################################
################################
"""
Hamiltonian dynamics
   q̇ =  ∂H/∂p
   -ṗ = ∂H/∂q
"""

# Define a general Hamiltonian problem
struct HamSys{H,G,R,E}
    H::H          # Hamiltonian
    guard::G          # guard/event function
    reset::R      # reset map
    e::E          # coefficient of restitution 
end

# Make the the guard, reset map, and coefficient of restitution optional; default to fully elastic spectral reflection
function HamSys(H; guard=nothing, reset = (x,e) -> spectral_refl(x,e), e = 1.0)
    HamSys(H, guard, reset, e)
end

# Solve a Hamiltonian problem
function solve(prob::Prob{HamSys}, solver; kwargs...)
    
    # Extract the Hamiltonian
    system = prob.sys
    H = system.H

    # Create the vector field from the hamiltonian
    F(x, t) = hamiltonian_vec_field(H, x, t)

    # Solve for trajectories along the above vector field
    return solver(F, prob.init, prob.tspan; kwargs...)
end

# Find the vector field from Hamilton's equations
function hamiltonian_vec_field(Hsys, x::AbstractVector, t)

    n = length(x) ÷ 2

    # Split up state vector
    q = x[1:n]
    p = x[n+1:end]

    # Gradient of H with respect to state vector x = [q; p], this is where multiple dispatch takes care of different kind of hams
    gradH = hamiltonian_gradient(Hsys, x)

    dqdt = gradH[n+1:end]      # ∂H/∂p
    dpdt = -gradH[1:n]         # -∂H/∂q

    # Recombine derivatives into single vecotr field
    return vcat(dqdt, dpdt)
end

# Find the gradient of Hamiltonians that can be differentiated using ForwardDiff
function hamiltonian_gradient(Hsys, x)
    ForwardDiff.gradient(Hsys.H, x)
end

# # Some way to find the gradient without ForwardDiff (untested)
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