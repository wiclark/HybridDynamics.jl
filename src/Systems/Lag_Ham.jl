# I want to understand hybrid Lagrangian systems (Hamiltonian will follow)

# using ForwardDiff

struct HybridLagrangianSystem
    L::Function
    h::Function
    e::Function
end

# I want the following:
# 1. Generate trajectories (you can call forward diff)
# 2. Implement an intelligent way to perform event detection
# 3. Determine whether or not the system is Zeno and state when/where that occurs.

# Allow options for how the solver arrives at the equations of motion (EOM)
abstract type Backend end
struct AutoForwardDiff <: Backend end   # Automatic differentiation using ForwardDiff
struct AutoFiniteDiff  <: Backend end   # Finite differences
struct ManualEOM       <: Backend end   # Manually provide sumthin?


# Default reset map: spectral reflection with coefficient of restitution
""" 
See "Is There Life After Zeno? paper
"""
function spectral_refl(x, M, dh; e=1.0)

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

################################
################################
## Lagrangian dynamics

# General Lagrangian system
struct LagSys{L,G,R,E,B}
    L::L          # Lagrangian
    guard::G      # guard/event function
    reset::R      # reset map
    e::E          # coefficient of restitution
    B::B          # backend use to arrive at EOM
end

# Make the the guard, reset map, and coefficient of restitution optional; default to fully elastic spectral reflection
"""
explaining the thing and the kwargs
 - L

"""
function LagSys(L;
            guard=nothing,
            reset = (x, M, dh, e) -> spectral_refl(x, M, dh; e),
            e = 1.0,
            B = AutoForwardDiff())

    return LagSys(L, guard, reset, e, B)
end

# Solution object for Lagrangian problems
struct LagSol{T, X, T_e1, T_z}
    T::T        # Time data
    X::X        # Position data
    T_e1::T_e1  # Time of first event
    T_z::T_z    # Time of Zeno
end

# Does this need to be a general constructor or can I go straight to initializing the solution state?
function LagSol

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
    function lagrangian_vec_field(L::Function, x::AbstractVector, t)

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

# Option to input M, C, and F directly
    struct ManualLag{M, C, F}
        M::M    # Mass matrix, could be noninvertable
        C::C    # Coriolis matrix
        F::F    # Force
    end

    function lagrangian_vec_field(L::ManualLag, x::AbstractVector, t)
        # Unpack matrices
        M, C, F = L

        # Calculate the vector field


        return
    end

# Solve a Lagrangian problem
function solve(prob::prob{LagSys}, solver; kwargs...) # solver specifically for LagProb struct, kwargs for step size or tolerances (optional / solver dependent)
    
    sys = prob.sys
    f = lagrangian_vec_field(sys.L, x, t)
    sol = initsol(prob)     # See line 85

    



    return 
end

################################
################################
## Hamiltonian dynamics

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
function solve(prob::prob{HamSys}, solver; kwargs...)
    
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