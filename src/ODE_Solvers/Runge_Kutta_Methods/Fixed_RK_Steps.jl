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
compute_step(::ForwardEuler, f, xₖ, Δt, t) = forward_euler_step(f, xₖ, Δt, t)
compute_step(::ModifiedTrap, f, xₖ, Δt, t) = modified_trap_step(f, xₖ, Δt, t)
compute_step(::ModifiedMidpoint, f, xₖ, Δt, t) = modified_midpoint_step(f, xₖ, Δt, t)
compute_step(::RichardsonExtrapolation, f, xₖ, Δt, t) = richardson_step(f, xₖ, Δt, t)
compute_step(::RK4, f, xₖ, Δt, t) = rk_4_step(f, xₖ, Δt, t)
# Implicit methods - these require the Jacobian
compute_step(::BackwardEuler, f, Df, xₖ, Δt, t) = backward_euler_step(f, Df, xₖ, Δt, t)
compute_step(::ImplicitTrap, f, Df, xₖ, Δt, t) = implicit_trap_step(f, Df, xₖ, Δt, t)
compute_step(::RadauIIA, f, Df, xₖ, Δt, t) = radauiia_step(f, Df, xₖ, Δt, t)

#Note sol is not used, we do this to make using the function easier. We would need an if/else statement everytime we use this function without it
function take_step(solver::FixedRK, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver=ModifiedMidpoint(); 
        check=true, guard_direction=default_guard_direction(prob.sys)) 
        
    sys = prob.sys
    x_predict = compute_step(solver, f, xₖ, Δt, tₖ)
    
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

# Implicit Methods - If Df is not given, it is reasonable to use ForwardDiff
# Backward (implicit) Euler
function backward_euler_step(f::Function, Df::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat; maxIter::Int=10, tol::AbstractFloat=1e-3)
	# Set up the function and its derivative
	G(y) = y - xₖ - Δt*f(t+Δt, y)
	DG(y) = I - Δt*Df(t+Δt, y)
	# Loop
	y0 = xₖ
	for i ∈ 1:maxIter
		y0 = y0 - DG(y0) \ G(y0)
	end
	if norm(G(y0)) > tol
		@warn "Implicit Euler step did not converge"
	end
	return y0
end

# The implicit trapezoid method (Lobatto IIIA, order 2)
function implicit_trapezoid_step(f::Function, Df::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat; maxIter::Int=10, tol::AbstractFloat=1e-3)
	# Set up the function and its derivative
	G(y) = y - xₖ - 1/2*Δt*(f(t, xₖ) + f(t+Δt, y))
	DG(y) = I - 1/2*Δt*Df(t+Δt, y)
	# Loop
	y0 = xₖ
	for i ∈ 1:maxIter
		y0 = y0 - DG(y0) \ G(y0)
	end
	if norm(G(y0)) > tol
		@warn "Implicit trapezoid step did not converge"
	end
	return y0
end

# The Radau IIA method (2 stages, order 3)
function implicit_RadauIIA_step(f::Function, Df::Function, xₖ::AbstractArray, Δt::AbstractFloat, t::AbstractFloat; maxIter::Int=10, tol::AbstractFloat=1e-3)
	n = length(xₖ)
	# Set up the function
	G1(k1, k2) = k1 .- f(t+1/3*Δt, xₖ+Δt*(5/12*k1-1/12*k2))
	G2(k1, k2) = k2 .- f(t+Δt, xₖ+Δt*(3/4*k1+1/4*k2))
	G(k) = vcat( G1(k[1:n], k[n+1:end]), G2(k[1:n], k[n+1:end]) )
	# and its derivative
	DG11(k1, k2) = I - 5/12*Δt*Df(t+1/3*Δt, xₖ+Δt*(5/12*k1-1/12*k2))
	DG12(k1, k2) = 1/12*Δt*Df(t+1/3*Δt, xₖ+Δt*(5/12*k1-1/12*k2))
	DG21(k1, k2) = -3/4*Δt*Df(t+Δt, xₖ+Δt*(3/4*k1+1/4*k2))
	DG22(k1, k2) = I - 1/4*Δt*Df(t+Δt, xₖ+Δt*(3/4*k1+1/4*k2))
	DGM(k1, k2) = [DG11(k1, k2) DG12(k1, k2);DG21(k1, k2) DG22(k1, k2)]
	DG(k) = DGM(k[1:n], k[n+1:end])
	# Loop
	guess = f(t, xₖ)
	y0 = vcat(guess, guess)
	for i ∈ 1:maxIter
		y0 = y0 - DG(y0) \ G(y0)
	end
	if norm(G(y0)) > tol
		@warn "Implicit Radau step did not converge"
	end
	K1, K2 = y0[1:n], y0[n+1:end]
	return xₖ + 1/4*Δt*(3*K1+K2)
end