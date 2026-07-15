#Exponential Solver
struct ExponentialSolver <: AbstractODESolver end

struct MagnusLeapfrog <: FixedRK end

#Supposedly useful as when doing variational equation stuff the magnus expansion preserves Lie group structure (I dont know enough about Lie groups to tell you what that means)
#so we maintain the properties we want like volume and determinants which is very good for variational equations and eventual Lyapunov Exponents. 
function take_step(solver::MagnusLeapfrog, prob::AbstractHybridProblem, f, Uₖ, tₖ, Δt, tol, sol, stepper=nothing; guard_direction=default_guard_direction(prob.sys))    
    sys = prob.sys

    #Compute augmented matrix predictions
    U_predict = compute_step(solver, f, Uₖ, Δt, tₖ)

    #Guard checks only care about phys state
    xₖ = Uₖ[:, 1]
    x_predict = U_predict[:, 1]

    h_now = guard(sys, xₖ)
    h_next = guard(sys, x_predict)

    idx = max(1, length(sol.x) - 1)
    t_prev = sol.t[idx]
    U_prev = sol.x[idx]
    x_prev = U_prev[:, 1]
    h_prev = guard(sys, x_prev)

    #Use check
    eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt; tol=tol, direction=guard_direction)

    if eventtrigger
        if (t_root - tₖ) < (1e-4 * Δt)
            eventtrigger = false
            t_root = tₖ + Δt
        end
    end

    return U_predict, eventtrigger, t_root, Δt, Δt
end

#Newton Raphson solver to find roots for implicit integration steps. 
# Solves for x in G(x) = x - c - α*h*f(x, t_new) = 0
function implicit_newton_solve(f::Function, x_guess::Vector, c::Vector, α::AbstractFloat, h::AbstractFloat, t_new::AbstractFloat; max_iter=50, tol=1e-12)
    x_curr = copy(x_guess)
    for _ in 1:max_iter
        #Eval vector field at curent guess for z_{n+1}
        val = f(x_curr, t_new)
        #construct residual function G(z) = z - c - αhf(z). When G(z) = 0 we have a sol
        G = x_curr .- c .- (α*h) .* val

        #Convergence check
        if norm(G) < tol
            return x_curr
        end

        #Calc Jacobian of vector field at our current guess
        J_f = ForwardDiff.jacobian(y -> f(y, t_new), x_curr)
        #Calc Jacobian of the residual G(z). 
        J_G = I - (α*h) * J_f

        #Update based on Newtons method. 
        x_curr = x_curr .- J_G \ G
    end
    #Just in case
    @warn "Newton Raphson solver failed to converge at t = $t_new. Consider a smaller step size."
    #return NaN to signal the sim that this step is invalid. 
    return fill(NaN, size(x_curr))
end


compute_step(::BackwardEuler, f, xₖ, Δt, t) = implicit_euler_step(f, xₖ, Δt, t)
function implicit_euler_step(f::Function, xₖ::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    t_new = t + h
    #Initial Guess via exp Euler
    x_guess = forward_euler_step(f, xₖ, h, t)

    #Imp euler 
    return implicit_newton_solve(f, x_guess, xₖ, 1.0, h, t_new)
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