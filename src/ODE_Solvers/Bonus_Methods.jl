#Exponential Solver
struct ExponentialSolver <: AbstractODESolver end

struct MagnusLeapfrog <: FixedRK end

#Supposedly useful as when doing variational equation stuff the magnus expansion preserves Lie group structure (I dont know enough about Lie groups to tell you what that means)
#so we maintain the properties we want like volume and determinants which is very good for variational equations and eventual Lyapunov Exponents. 
function take_step(solver::MagnusLeapfrog, prob::AbstractHybridProblem, f, Uₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver=ModifiedMidpoint(), guard_direction=default_guard_direction(prob.sys))
    sys = prob.sys

    #Compute augmented matrix predictions
    U_predict = compute_step(solver, f, Uₖ, Δt, tₖ)
    U_mid     = compute_step(solver, f, Uₖ, Δt / 2.0, tₖ)

    #Guard checks only care about phys state
    xₖ = Uₖ[:, 1]
    x_mid = U_mid[:, 1]
    x_predict = U_predict[:, 1]

    #Eval guards
    h_now = guard(sys, xₖ)
    h_mid = guard(sys, x_mid)
    h_next = guard(sys, x_predict)

    #Use check
    eventtrigger, t_root, _ = crossed_guard(sys, h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol, direction=guard_direction)

    return U_predict, eventtrigger, t_root, Δt, Δt
end

#Newton Raphson solver to find roots for implicit integration steps. 
#Solves for z in G(z) = z - c - α*h*f(z, t_new) = 0
function implicit_newton_solve(f::Function, z_guess::Vector, c::Vector, α::AbstractFloat, h::AbstractFloat, t_new::AbstractFloat; max_iter=50, tol=1e-6)
    z_curr = copy(z_guess)

    for _ in 1:max_iter
        #Eval vector field at curent guess for z_{n+1}
        val = f(z_curr, t_new)
        #construct residual function G(z) = z - c - αhf(z). When G(z) = 0 we have a sol
        G = z_curr .- c .- (α*h) .* val

        #Convergence check
        if norm(G) < tol
            return z_curr
        end

        #Calc Jacobian of vector field at our current guess
        J_f = ForwardDiff.jacobian(y -> f(y, t_new), z_curr)
        #Calc Jacobian of the residual G(z). 
        J_G = I - (α*h) * J_f

        #Update based on Newtons method. 
        z_curr = z_curr .- J_G \ G
    end
    #Just in case
    @warn "Newton Raphson solver failed to converge at t = $t_new. Consider a smaller step size."
    #return NaN to signal the sim that this step is invalid. 
    return fill(NaN, size(z_curr))
end


compute_step(::ImplicitEuler, f, x, Δt, t) = implicit_euler_step(f, x, Δt, t)
function implicit_euler_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    t_new = t + h
    #Initial Guess via exp Euler
    z_guess = forward_euler_step(f, z, h, t)

    #Imp euler 
    return implicit_newton_solve(f, z_guess, z, 1.0, h, t_new)
end

#Extra solvers
compute_step(::MagnusLeapfrog, f, U::AbstractMatrix, Δt, t) = magnus_leapfrog_step(f, U, Δt, t)
#MagnusLeapfrog step
function magnus_leapfrog_step(f::Function, U::AbstractMatrix, h::AbstractFloat, t::AbstractFloat)
    #Extract state and fund matrix
    xₖ = U[:, 1]
    Φₖ = U[:, 2:end]

    #Midpoint Eval for base trajectory 
    x_mid = xₖ .+ (h / 2.0) .* f(xₖ, t)
    t_mid = t + h / 2.0

    #Advance base state fully
    x_next = xₖ .+ h .* f(x_mid, t_mid)

    #Compute Jacobian exactly at midpoint
    A_mid = ForwardDiff.jacobian(y -> f(y, t_mid), x_mid)

    #Lie-Group step for fund matrix 
    Φ_next = exp(h * A_mid) * Φₖ

    return hcat(x_next, Φ_next)
end