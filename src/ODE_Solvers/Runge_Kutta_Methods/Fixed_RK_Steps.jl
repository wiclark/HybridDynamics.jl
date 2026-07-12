#Non-adaptive solvers
abstract type FixedRK <: RK end
struct ForwardEuler <: FixedRK end
struct ModifiedTrap <: FixedRK end
struct ModifiedMidpoint <: FixedRK end
struct RichardsonExtrapolation <: FixedRK end
struct RK4 <: FixedRK end
struct BackwardEuler <: FixedRK end

#Single take_step for RK methods
#Fixed step methods 
#Helper function to take the step using multiple dispatch.
compute_step(::ForwardEuler, f, xₖ, Δt, t) = forward_euler_step(f, xₖ, Δt, t)
compute_step(::ModifiedTrap, f, xₖ, Δt, t) = modified_trap_step(f, xₖ, Δt, t)
compute_step(::ModifiedMidpoint, f, xₖ, Δt, t) = modified_midpoint_step(f, xₖ, Δt, t)
compute_step(::RichardsonExtrapolation, f, xₖ, Δt, t) = richardson_step(f, xₖ, Δt, t)
compute_step(::RK4, f, xₖ, Δt, t) = rk_4_step(f, xₖ, Δt, t)

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
