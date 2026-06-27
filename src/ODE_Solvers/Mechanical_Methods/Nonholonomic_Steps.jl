# Modify various ODE solvers to handle nonholonomic constraints - a special type of DAE

struct NH_RK4 <: NH end

compute_step(::NH_RK4, f, A, zₖ, Δt, t) = nonholonomic_rk4(f, A, zₖ, Δt, t)

function take_step(solver::NH, prob::NonholonomicSystem, f, g, A, qₖ, pₖ, Δt, t, sol; check=true, guard_direction=default_guard_direction(prob.sys))
    sys = prob.sys
    z_predict = compute_step(solver, f, A, zₖ, Δt, t)
    z_mid     = compute_step(solver, f, A, zₖ, Δt/2.0, t)

    if check
        # Evaluate guards (only along the position)
        h_now = guard(sys, z[1:div(length(z), 2)])
        h_mid = guard(sys, z_mid[1:div(length(z), 2)])
        h_next = guard(sys, z_predict[1:div(length(z), 2)])

        # Check
        eventtrigger, t_root, _ = crossed_guard(sys, h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol, direction=guard_direction)
        return z_predict, eventtrigger, t_root, Δt, Δt
    else
        return z_predict, NaN, NaN, NaN, NaN
    end
end

# Compute the reaction force
function find_reaction(z, M, A)
    q, p = z[1:div(length(z), 2)], z[(div(length(z),2)+1):end]
    A_now, M_now = A(q), M(q)
    # Solving the linear system
    LHS = A_now * (M_now \ p)
    RHS = -A_now * (M_now \ A_now')
    λ = RHS \ LHS
    F = A_now' * λ
    p_forced = p + F
    return vcat(q, p_forced) - z
end

function nonholonomic_rk4(f::Function, A::Function, z::AbstractArray, Δt::AbstractFloat, t::AbstractFloat)
    k1 = find_reaction(z+Δt*f(z, t), M, A)
    k2 = find_reaction(z+Δt*f(z+Δt/2*k1, t+Δt/2), M, A)
    k3 = find_reaction(z+Δt*f(z+Δt/2*k2, t+Δt/2), M, A)
    k4 = find_reaction(z+Δt*f(z+Δt*k3, t+Δt), M, A)
    
    z_guessed = z + Δt/6 * (k1 + 2*k2 + 2*k3 + k4)

    return z_guessed + find_reaction(z_guessed, M, A)
end