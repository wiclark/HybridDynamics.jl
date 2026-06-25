#Non-adaptive solvers
abstract type FixedRK <: RK end
struct ForwardEuler <: FixedRK end
struct ModifiedTrap <: FixedRK end
struct ModifiedMidpoint <: FixedRK end
struct RichardsonExtrapolation <: FixedRK end
struct RK4 <: FixedRK end
struct ImplicitEuler <: FixedRK end

#Single take_step for RK methods
#Fixed step methods 
#Helper function to take the step using multiple dispatch.
compute_step(::ForwardEuler, f, x, Δt, t) = forward_euler_step(f, x, Δt, t)
compute_step(::ModifiedTrap, f, x, Δt, t) = modified_trap_step(f, x, Δt, t)
compute_step(::ModifiedMidpoint, f, x, Δt, t) = modified_midpoint_step(f, x, Δt, t)
compute_step(::RichardsonExtrapolation, f, x, Δt, t) = richardson_step(f, x, Δt, t)
compute_step(::RK4, f, x, Δt, t) = rk_4_step(f, x, Δt, t)

#Note sol is not used, we do this to make using the function easier. We would need an if/else statement everytime we use this function without it
function take_step(solver::FixedRK, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver=ModifiedMidpoint(); check=true, guard_direction=default_guard_direction(prob.sys)) 
    sys = prob.sys
    x_predict = compute_step(solver, f, xₖ, Δt, tₖ)
    x_mid     = compute_step(solver, f, xₖ, Δt / 2.0, tₖ)

    if check
        #Evaluate Guards
        h_now  = guard(sys, xₖ)
        h_mid  = guard(sys, x_mid)
        h_next = guard(sys, x_predict)

        #Use cross guard check
        eventtrigger, t_root, _ = crossed_guard(sys, h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol, direction=guard_direction)

        return x_predict, eventtrigger, t_root, Δt, Δt
    else
        return x_predict, NaN, NaN, NaN, NaN
    end
end

# Forward Euler
function forward_euler_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    return z .+ h*f(z, t)
end

# Modified (fully explicit) Trapezoid Rule
function modified_trap_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    z_guess = z .+ h*f(z,t)
    return z .+ 1/2*h*( f(z, t) + f(z_guess, t+h) )
end

# Modified (fully explicit) Midpoint Rule
function modified_midpoint_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    z_guess = z .+ h/2*f(z, t)
    return z .+ h*f(z_guess, t+h/2)
end

#Classic Runge-Kutta 4
function rk_4_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    k1 = f(z,t)
    k2 = f(z + h/2 * k1, t + h/2)
    k3 = f(z + h/2 * k2, t + h/2)
    k4 = f(z + h * k3, t + h)

    #Shouldnt need to check guard in these inner stages here as no adpative step size.

    return z + h/6 * (k1 + 2*k2 + 2*k3 + k4)
end

#Richardson Extrapolation based on modifed midpoint from Wiki
function richardson_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    #Take one full step size h
    z1 = modified_midpoint_step(f,z,h,t)

    #Take two smaller steps of h/2
    h_half = h / 2
    z_half = modified_midpoint_step(f,z,h_half,t)
    z2 = modified_midpoint_step(f, z_half,h_half,t+h_half)

    #Extrapolate to get rid of lower order error
    return (4 .* z2 .- z1) ./ 3.0
end
