# This is not where I want it yet. Ignore this file.

abstract type FixHam <: ME end

struct SymplecticEuler <: FixHam end

compute_step(::SymplecticEuler, f, g, q, p, Δt, t) = symplectic_euler_step(f, g, q, p, Δt, t)

function take_step(solver::FixHam, prob::MechanicalSystem, f, g, qₖ, pₖ, tₖ, Δt, sol; check-true, guard_direction = default_guard_direction(prob.sys))
    sys = prob.sys
    q_predict, p_predict = compute_step(solver, f, g, qₖ, pₖ, Δt, tₖ)
    q_mid, _             = compute_step(solver, f, g, qₖ, pₖ, Δt/2.0, tₖ)
    # Check whether or not there is an impact. Notice that h is only a function of q, rather than (q,p)
    if check
        # Evaluate guards
        h_now  = guard(sys, qₖ)
        h_mid  = guard(sys, q_mid)
        h_next = guard(sys, q_predict)
        # Use cross guard check
        eventtrigger, t_root, _ = crossed_guard(sys, h_now, h_mid, h_next, tₖ, tₖ+Δt/2.0, tₖ+Δt; tol=tol, direction=guard_direction)
        return q_predict, p_predict, eventtrigger, t_root, Δt, Δt
    else
        return q_predict, p_predict, NaN, NaN, NaN
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