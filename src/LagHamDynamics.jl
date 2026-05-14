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

########
"""
Lagrangian dynamics
"""

### Solving for trajectories from normal Lagrangians:

struct LagProb{F, I, T}
    sys::F
    init::I
    tspan::T
end

function solve(prob::LagProb, solver; kwargs...) # solver specifically for LagProb struct, kwargs for step size or tolerances (optional / solver dependent)
    F(x, t) = lagrangian_vec_field(prob.sys, x, t)
    return solver(F, prob.init, prob.tspan; kwargs...)
end

### Trajectories from hybrid Lagrangian systems:

# struct HybridLagSys{L, h; e = 1.0}  # is this how kwargs work with structs?
#     L::L
#     h::h
#     e::e
# end

# struct LagProb{F::HybridLagSys, I, T}
#     sys::F
#     init::I
#     tspan::T
# end

# function solve(prob::HybridLagSys)

# end

### Find the equations of motion using Euler-Lagrange equations

# dL/dq is the force
function lagrangian_force(L::Function, q::AbstractVector, qdot::AbstractVector)
    ForwardDiff.gradient(q -> L(q, qdot), q)
end

# dL/dqdot is the momentum
function lagrangian_momentum(L::Function, q::AbstractVector, qdot::AbstractVector)
    ForwardDiff.gradient(qdot -> L(q, qdot), qdot)
end

# find mass matrix
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


############
"""
Hamiltonian dynamics
   q̇ =  ∂H/∂p
   -ṗ = ∂H/∂q
"""

struct HamProb{H, I, T}
    sys::H   # Hamiltonian
    init::I   # Initial condition
    tspan::T   # Time span
end

# Solver dispatch specifically for HamProb struct
function solve(prob::HamProb, solver; kwargs...)
    
    # Create the vector field from the hamiltonian
    F(x, t) = hamiltonian_vec_field(prob.sys, x, t)

    # Solve for trajectories along the above vector field
    return solver(F, prob.init, prob.tspan; kwargs...)
end

function hamiltonian_vec_field(Hsys::AbstractHamiltonian, x::AbstractVector, t)

    n = length(x) ÷ 2

    # Split up state vector
    q = x[1:n]
    p = x[n+1:end]

    # Gradient of H with respect to state vector x = [q; p]
    gradH = hamiltonian_gradient(Hsys, x)

    dqdt = gradH[n+1:end]      # ∂H/∂p
    dpdt = -gradH[1:n]         # -∂H/∂q

    # Recombine derivatives into single vecotr field
    return vcat(dqdt, dpdt)
end


### Allow for different representations of the Hamiltonian
abstract type AbstractHamiltonian end

# Auto differentiation Ham struct; for closed form Hamiltonian functions
struct ADHamiltonian{F} <: AbstractHamiltonian
    H::F
end

# dispatch to find the gradient of Hamiltonians that can be differentiated using ForwardDiff
function hamiltonian_gradient(Hsys::ADHamiltonian, x)
    ForwardDiff.gradient(Hsys.H, x)
end

# Finite differences struct for Hamiltonians that can't be cleanly differentiated
struct FDHamiltonian{F,T} <: AbstractHamiltonian
    H::F
    h::T   # finite difference step size
end

# Some way to find the gradient without ForwardDiff (untested)
function hamiltonian_gradient(Hsys::FDHamiltonian, x)
    
    H = Hsys.H
    h = Hsys.h

    # initialize place to store the output
    grad = similar(x)

    # This feels like a bad way of doing this but I think it's correct
    for i in eachindex(x)

        xp = copy(x)
        xm = copy(x)

        # perturb the point forward and backward
        xp[i] += h
        xm[i] -= h

        # Central finite diff approximation
        grad[i] = (H(xp) - H(xm)) / (2h)
    end

    return grad
end


# Could add another dispatch method to better deal with interpolating data style hams