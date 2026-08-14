#Non-adaptive solvers
abstract type FixedRK <: RK end

"""
Forward (explicit) Euler method.
"""
struct ForwardEuler <: FixedRK end
"""
Backward (implicit) Euler method.
"""
struct BackwardEuler <: FixedRK end
"""
Modified trapezoidal method.
"""
struct ModifiedTrap <: FixedRK end
"""
Modified midpoint method.
"""
struct ModifiedMidpoint <: FixedRK end
"""
Richardson Extrapolation.
"""
struct RichardsonExtrapolation <: FixedRK end
"""
Fourth-order Runge-Kutta.
"""
struct RK4 <: FixedRK end
"""
Implicit trapezoid (Lobatto IIIA, order 2)
"""
struct ImplicitTrap <: FixedRK end
"""
Implicit RK (Radau IIA, order 3)
"""
struct RadauIIA <: FixedRK end


#Single take_step for RK methods
#Fixed step methods 
#Helper function to take the step using multiple dispatch.
compute_step(::ForwardEuler, f, Df, xₖ, Δt, t) = forward_euler_step(f, xₖ, Δt, t)
compute_step(::ModifiedTrap, f, Df, xₖ, Δt, t) = modified_trap_step(f, xₖ, Δt, t)
compute_step(::ModifiedMidpoint, f, Df, xₖ, Δt, t) = modified_midpoint_step(f, xₖ, Δt, t)
compute_step(::RichardsonExtrapolation, f, Df, xₖ, Δt, t) = richardson_step(f, xₖ, Δt, t)
compute_step(::RK4, f, Df, xₖ, Δt, t) = rk_4_step(f, xₖ, Δt, t)

# Implicit methods - these require the Jacobian
compute_step(::BackwardEuler, f, Df, xₖ, Δt, t) = backward_euler_step(f, Df, xₖ, Δt, t)
compute_step(::ImplicitTrap, f, Df, xₖ, Δt, t) = implicit_trap_step(f, Df, xₖ, Δt, t)
compute_step(::RadauIIA, f, Df, xₖ, Δt, t) = radauiia_step(f, Df, xₖ, Δt, t)
# In case DF is omitted
compute_step(solver::BackwardEuler, f, xₖ, Δt, t) = backward_euler_step(f, nothing, xₖ, Δt, t)
compute_step(solver::ImplicitTrap, f, xₖ, Δt, t)  = implicit_trap_step(f, nothing, xₖ, Δt, t)
compute_step(solver::RadauIIA, f, xₖ, Δt, t)     = radauiia_step(f, nothing, xₖ, Δt, t)

#Note sol is not used, we do this to make using the function easier. We would need an if/else statement everytime we use this function without it
function take_step(solver::FixedRK, prob::AbstractHybridProblem, f, Df, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver=ModifiedMidpoint(); 
        check=true, guard_direction=default_guard_direction(prob.sys)) 
        
    sys = prob.sys
    x_predict = compute_step(solver, f, Df, xₖ, Δt, tₖ)
    
    if check
        #Evaluate Guards
        h_now  = guard(sys, xₖ)
        h_next = guard(sys, x_predict)

        idx = max(1, length(sol.x) - 1)
        t_prev = sol.t[idx]
        x_prev = sol.x[idx]
        h_prev = guard(sys, x_prev)

        #Use cross guard check
        eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt; tol=tol, direction=guard_direction)

        if eventtrigger
            if (t_root - tₖ) < (1e-4 * Δt)
                eventtrigger = false
                t_root = tₖ + Δt
            end
        end

        return x_predict, eventtrigger, t_root, Δt, Δt
    else
        return x_predict, false, NaN, Δt, Δt
    end
end

##################
# Explicit Solvers
##################

# Forward Euler
function forward_euler_step(f::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat)
    return xₖ .+ Δt*f(xₖ, t)
end

# Modified (fully explicit) Trapezoid Rule
function modified_trap_step(f::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat)
    x_guess = xₖ .+ Δt*f(xₖ, t)
    return xₖ .+ 1/2*Δt*(f(xₖ, t) + f(x_guess, t+Δt))
end

# Modified (fully explicit) Midpoint Rule
function modified_midpoint_step(f::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat)
    x_guess = xₖ .+ Δt/2*f(xₖ, t)
    return xₖ .+ Δt*f(x_guess, t+Δt/2)
end

#Classic Runge-Kutta 4
function rk_4_step(f::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat)
    k1 = f(xₖ,t)
    k2 = f(xₖ + Δt/2 * k1, t + Δt/2)
    k3 = f(xₖ + Δt/2 * k2, t + Δt/2)
    k4 = f(xₖ + Δt * k3, t + Δt)

    #Shouldnt need to check guard in these inner stages here as no adpative step size.

    return xₖ + Δt/6 * (k1 + 2*k2 + 2*k3 + k4)
end

#Richardson Extrapolation based on modifed midpoint from Wiki
function richardson_step(f::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat)
    #Take one full step size h
    x1 = modified_midpoint_step(f,xₖ,Δt,t)

    #Take two smaller steps of h/2
    Δt_half = Δt / 2
    x_half = modified_midpoint_step(f, xₖ, Δt_half, t)
    x2 = modified_midpoint_step(f, x_half, Δt_half, t+Δt_half)

    #Extrapolate to get rid of lower order error
    return (4 .* x2 .- x1) ./ 3.0
end

##################
# Implicit Solvers
##################

# Helper for jacobian evaluations
@inline function eval_jacobian(f, Df, x, t)
	if Df !== nothing
		return Df(x,t)
	else 
		return ForwardDiff.jacobian(y -> f(y, t), x)
	end
end

# Backward Euler when Df is omitted. 
backward_euler_step(f, xₖ, Δt, t; kwargs...) = backward_euler_step(f, nothing, xₖ, Δt, t; kwargs...)

# Backward (implicit) Euler
function backward_euler_step(f::Function, Df::Union{Function, Nothing}, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat; maxIter::Int=10, tol::AbstractFloat=1e-9)
    y = copy(xₖ)
    t_next = t + Δt

    for _ in 1:maxIter
        G = y .- xₖ .- Δt .* f(y, t_next)
        
        # Absolute check to avoid Jacobians if we are at equilibrium
        if norm(G) < 1e-14
            return y
        end

        J = eval_jacobian(f, Df, y, t_next)
        DG = I - Δt .* J
        update = DG \ G
        y = y .- update
        
        # Exit when the Newton step stops meaningfully changing the state
        # this is different from the original as it would cause us to fail to take any step in some cases.
        if norm(update) < tol
            return y
        end
    end

    # Final residual check
    G_final = y .- xₖ .- Δt .* f(y, t_next)
    if norm(G_final) > max(tol * 10, 1e-6)
        @warn "Backward Euler step did not converge."
    end
    return y
end

# Implicit Trapezoid when Df is omitted.
implicit_trap_step(f, xₖ, Δt, t; kwargs...) = implicit_trap_step(f, nothing, xₖ, Δt, t; kwargs...)

# The implicit trapezoid method (Lobatto IIIA, order 2)
function implicit_trap_step(f::Function, Df::Union{Function, Nothing}, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat; maxIter::Int=10, tol::AbstractFloat=1e-9)
    y = copy(xₖ)
    t_next = t + Δt
    fₖ = f(xₖ, t)

    for _ in 1:maxIter
        G = y .- xₖ .- 1/2 .* Δt .* (fₖ .+ f(y, t_next))

        if norm(G) < 1e-14
            return y
        end

        J = eval_jacobian(f, Df, y, t_next)
        DG = I - 1/2 .* Δt .* J
        update = DG \ G
        y = y .- update

        # Exit when the Newton step stops meaningfully changing the state
        # this is different from the original as it would cause us to fail to take any step in some cases.
        if norm(update) < tol
            return y
        end
    end

    G_final = y .- xₖ .- 1/2 .* Δt .* (fₖ .+ f(y, t_next))
    if norm(G_final) > max(tol * 10, 1e-6)
        @warn "Implicit Trapezoid step did not converge."
    end
    return y
end

# Radau IIA when DF is omitted.
radauiia_step(f, xₖ, Δt, t; kwargs...) = radauiia_step(f, nothing, xₖ, Δt, t; kwargs...)

# The Radau IIA method (2 stages, order 3)
function radauiia_step(f, Df::Union{Function, Nothing}, xₖ, Δt, t; maxIter::Int=10, tol=1e-9)
    n = length(xₖ)
    t1 = t + Δt / 3
    t2 = t + Δt

    # Initialize stage derivatives k1, k2 with current derivative guess
    f_k = f(xₖ, t)
    k = vcat(f_k, f_k)

    for _ in 1:maxIter
        k1 = k[1:n]
        k2 = k[n+1:end]

        z1 = xₖ .+ Δt .* ((5/12) .* k1 .- (1/12) .* k2)
        z2 = xₖ .+ Δt .* ((3/4) .* k1 .+ (1/4) .* k2)

        G1 = k1 .- f(z1, t1)
        G2 = k2 .- f(z2, t2)
        G = vcat(G1, G2)

        if norm(G) < 1e-14
            return xₖ .+ Δt .* ((3/4) .* k1 .+ (1/4) .* k2)
        end

        J1 = eval_jacobian(f, Df, z1, t1)
        J2 = eval_jacobian(f, Df, z2, t2)

        In = I(n)
        DG11 = In .- (5/12 * Δt) .* J1
        DG12 = (1/12 * Δt) .* J1
        DG21 = -(3/4 * Δt) .* J2
        DG22 = In .- (1/4 * Δt) .* J2

        DG = [DG11 DG12; DG21 DG22]
        update = DG \ G
        k = k .- update

        # Exit when the Newton step stops meaningfully changing the state
        # this is different from the original as it would cause us to fail to take any step in some cases.
        if norm(update) < tol
            break
        end
    end

    k1 = k[1:n]
    k2 = k[n+1:end]
    z1 = xₖ .+ Δt .* ((5/12) .* k1 .- (1/12) .* k2)
    z2 = xₖ .+ Δt .* ((3/4) .* k1 .+ (1/4) .* k2)
    G_final = vcat(k1 .- f(z1, t1), k2 .- f(z2, t2))

    if norm(G_final) > max(tol * 10, 1e-6)
        @warn "Implicit Radau IIA step did not converge."
    end

    return xₖ .+ Δt .* ((3/4) .* k1 .+ (1/4) .* k2)
end