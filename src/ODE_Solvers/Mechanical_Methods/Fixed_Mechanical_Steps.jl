# This is not where I want it yet. Ignore this file.

abstract type FixHam <: ME end

"""
First order, fixed step, implicit symplectic method.
"""
struct SymplecticEuler <: FixHam end

compute_step(::SymplecticEuler, f, g, qₖ, pₖ, Δt, t) = symplectic_euler_step(f, g, qₖ, pₖ, Δt, t)

function take_step(solver::FixHam, prob::MechanicalSystem, f, g, qₖ, pₖ, tₖ, Δt, tol, sol; check=true, guard_direction = default_guard_direction(prob.sys))
    sys = prob.sys
    #Compute predicted position and momentum for the next time step
    q_predict, p_predict = compute_step(solver, f, g, qₖ, pₖ, Δt, tₖ)
    # Check whether or not there is an impact. Notice that h is only a function of q, rather than (q,p)
    if check
        # Evaluate guards
        h_now  = guard(sys, qₖ)
        h_next = guard(sys, q_predict)

        idx = max(1, length(sol.q) - 1)
        t_prev = sol.t[idx]
        q_prev = sol.q[idx]
        h_prev = guard(sys, q_prev)
        # Use cross guard check
        eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_next, t_prev, tₖ, tₖ+Δt; tol=tol, direction=guard_direction)

        if eventtrigger
            if (t_root - tₖ) < (1e-4 * Δt)
                eventtrigger = false
                t_root = tₖ + Δt
            end
        end

        return q_predict, p_predict, eventtrigger, t_root, Δt, Δt
    else
        return q_predict, p_predict, false, NaN, Δt, Δt
    end
end

# The actual steppers
function symplectic_euler_step(f, g, q, p, Δt, t)
    # This is an implicit method
    # We want to solve for the root of
    root(q_find) = q_find - q - Δt*f(q_find, p, t)
    root_deriv(q_find) = ForwardDiff.jacobian(root, q_find)
    q₀ = q
    for _ ∈ 1:10
        q₀ = q₀ - root_deriv(q₀) \ root(q₀)
        if norm(root(q₀)) < 1e-6
            break
        end
    end
    # p is found explicitly
    p₀ = p - Δt*g(q₀, p, t)
    return q₀, p₀
end