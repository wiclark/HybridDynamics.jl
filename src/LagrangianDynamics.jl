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

######

### Solving for trajectories from normal Lagrangians:

struct LagProb{F, I, T}
    sys::F
    init::I
    tspan::T
end

function solve(prob::LagProb, solver; kwargs...) # solver specifically for LagProb struct, kwargs for step size or tolerances (optional / solver dependent)
    F(x) = lagrangian_vec_field(prob.sys, x)
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

function lagrangian_vec_field(L::Function, x::AbstractVector)
    n = length(x) ÷ 2       #integrer division
    q = Vector(x[1:n])
    qdot = Vector(x[n+1:end])
    qddot = lagrangian_acceleration(L, q, qdot)
    vcat(qdot, qddot)
end

