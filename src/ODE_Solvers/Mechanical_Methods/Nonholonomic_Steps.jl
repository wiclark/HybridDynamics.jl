# This file is not where I want it. Ignore the contents.
# Modify various ODE solvers to handle nonholonomic constraints - a special type of DAE

"""
RK4 specific to nonholonomic systems.
"""
struct NH_RK4 <: NH end

compute_step(::NH_RK4, f, M, A, zₖ, Δt, t) = nonholonomic_rk4(f, M, A, zₖ, Δt, t)

function take_step(solver::NH, prob::NonholonomicSystem, f, g, M, A, qₖ, pₖ, Δt, tₖ, tol, sol; check=true, guard_direction=default_guard_direction(prob.sys))
    sys = prob.sys
    
    #Explicitly combine q and p into the state zₖ
    zₖ = vcat(qₖ, pₖ)
    
    z_predict = compute_step(solver, f, M, A, zₖ, Δt, tₖ)

    if check
        # Evaluate guards (only along the position portion of the state vector)
        h_now = guard(sys, zₖ[1:div(length(zₖ), 2)])
        
        # Evaluate guard for the predicted position
        h_next = guard(sys, z_predict[1:div(length(z_predict), 2)])

        # extract the previous state from history
        idx = max(1, length(sol.x) - 1)
        t_prev = sol.t[idx]
        z_prev = sol.x[idx]
        
        # Evaluate guard for the previous position
        h_prev = guard(sys, z_prev[1:div(length(z_prev), 2)])
        # Check
        eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt; tol=tol, direction=guard_direction)

        if eventtrigger 
            if (t_root - tₖ) < (1e-4 * Δt)
                eventtrigger = false
                t_root = tₖ + Δt
            end
        end
        return z_predict, eventtrigger, t_root, Δt, Δt
    else
        return z_predict, false, NaN, Δt, Δt
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

function nonholonomic_rk4(f::Function, M::Function, A::Function, z::AbstractArray, Δt::AbstractFloat, t::AbstractFloat)
    k1 = find_reaction(z + Δt*f(z, t), M, A)
    k2 = find_reaction(z + Δt*f(z + Δt/2*k1, t + Δt/2), M, A)
    k3 = find_reaction(z + Δt*f(z + Δt/2*k2, t + Δt/2), M, A)
    k4 = find_reaction(z + Δt*f(z + Δt*k3, t + Δt), M, A)
    z_guessed = z + Δt/6 * (k1 + 2*k2 + 2*k3 + k4)
    return z_guessed + find_reaction(z_guessed, M, A)
end