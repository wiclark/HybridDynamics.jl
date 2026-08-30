## A function that computes the tangent dynamics (variational equation) as well as the Lyapunov spectrum
# This only works for GeneralSystems. I would like to do Filippov and Mechanical next.
# Bugs: The current sol struct doesn't record the locations of resets, only their times. Moreover, if we evaluate 
#       the solution at the impact time, it provides the state immedeately post reset, not pre. 
#       This results in having to go slightly back in time, which is clumsy and awkward.
# Notes: I think I like the archecture of how this is laid out. All differentiation requires ForwardDiff.
#        There is currently no way to manually input the Jacobians.

## The solver will be a midpoint exponential solver. (I am ignoring Magnus currently.)
function exp_step(A, Φ, Δt, t)
    return exp(A(t+Δt/2)*Δt) * Φ
end

## I want to compute the tangent dynamics as well as the Lyapunov exponents
# Inputs:
#   sol::GeneralSol, This is written only for general systems. 
#                    The variational equation requires integrating along a trajectory
#   sys::GeneralSystem, To obtain the variational equation, we need to differentiate data
#                       from the original system. (This uses ForwardDiff.)
# Outputs:
#   Φ, the state-transition matrix at the terminal time (with I as the initial condition)
#   pre_mature, if switchings become too fast, integration will halt and this will be true
"""
    tangent_dynamics(sol::GeneralSol, sys::GeneralSystem; 
            Df = nothing, res_deriv = nothing, 
            guard_deriv = nothing, t_s = 1000)

Compute the state transition matrix across a trajectory from a General system.

"""
function tangent_dynamics(sol::GeneralSol, sys::GeneralSystem; Df = nothing, Dg = nothing, res_deriv = nothing, guard_deriv = nothing, t_s = 1000)
    # Extract out the event times
    T_events = sol.event_times

    # Create the tangent vector field
    A_int(t) =  isnothing(Df) ? ForwardDiff.jacobian(y->sys.f(y,t), sol(t)) : Df(sol(t))
    
    Φ = Matrix(I(size(A_int(sol.t[end]))[1]))
    t_past = sol.t[1]

    # Create the augmented differential matrix
    Δ_star(x) = ForwardDiff.jacobian(y->sys.Δ(y), x)
    adjust_vector(x,t) = sys.f(sys.Δ(x),t) - Δ_star(x)*sys.f(x,t)
    adjust_covector(x) = ForwardDiff.gradient(sys.h, x)'
    Δ_augmented_f(x,t) = Δ_star(x) .+ 1/(adjust_covector(x) * sys.f(x,t)) * (adjust_vector(x,t)*adjust_covector(x))

    # Loop through events using their index to get pre/post states
    for k in 1:length(T_events)
        t_e = T_events[k]

        times_segments = LinRange(t_past, t_e, t_s)
        dt = times_segments[2] - times_segments[1]
        for i ∈ 1:t_s -1
            Φ = exp_step(A_int, Φ, dt, times_segments[i])
        end

        # Exact states at the event boundary
        idx = sol.event_indices[k]
        x⁻ = sol.x[idx-1]   # pre-impact state

        # Eval augmented diff matrices at the event
        Δ_augmented = Δ_augmented_f(x⁻, t_e)

        Φ = Δ_augmented * Φ
        t_past = t_e
    end

    # Determine the tail
    t_current = isempty(T_events) ? sol.t[1] : T_events[end]
    t_final = sol.t[end]
    times_segment = LinRange(t_current, t_final, t_s)
    dt = times_segment[2] - times_segment[1]
    for i ∈ 1:t_s -1
        Φ = exp_step(A_int, Φ, dt, times_segment[i])
    end
    return Φ
end

## Tangent map for mechanical systems. Zeno does not currently work
"""
    tangent_dynamics(sol::MechanicalSol, sys::MechanicalSystem; 
            Df = nothing, res_deriv = nothing, 
            guard_deriv = nothing, t_s = 1000)

Compute the state transition matrix across a trajectory from a Mechanical system.

"""
function tangent_dynamics(sol::MechanicalSol, sys::MechanicalSystem; Df = nothing, Dg = nothing, res_deriv = nothing, guard_deriv = nothing, t_s = 1000)
    # Extract out the event times
    T_events = sol.event_times

    # Before we compute the variations across this trajectory, we check for Zeno
    if length(sol.zeno) > 0
        Φ = fill(NaN, size(A_int(sol.t[end])))
    end

    n = Int(length(sol(sol.t[1])) / 2)

    # We need the vector field as well as its Jacobian
    M(q) = sys.M(q)
    V(q) = sys.V(q)
    H(q,p) = 1/2*dot(p, M(q) \ p) + V(q)
    q_dot(q, p) = ForwardDiff.gradient(p -> H(q,p),p)
    p_dot(q, p) = ForwardDiff.gradient(q -> -H(q,p), q)
    f(y) = [q_dot(y[1:n],y[n+1:end]); p_dot(y[1:n],y[n+1:end])]

    # Create the tangent vector field
    A_int(t) =  isnothing(Df) ? ForwardDiff.jacobian(y->f(y), sol(t)) : Df(sol(t))
    
    Φ = Matrix(I(size(A_int(sol.t[end]))[1]))
    t_past = sol.t[1]

    # Create the augmented differential matrix
    Δ(x) = sys.reset(x, sys.M, sys.normal, sys)
    Δ_star(x) = ForwardDiff.jacobian(y->Δ(y), x)
    adjust_vector(x) = f(Δ(x)) - Δ_star(x)*f(x)
    adjust_covector(x) = ForwardDiff.gradient(sys.guard, x)'
    Δ_augmented_f(x) = Δ_star(x) .+ 1/(adjust_covector(x) * f(x)) * (adjust_vector(x)*adjust_covector(x))

    # Loop through events using their index to get pre/post states
    for k in 1:length(T_events)
        t_e = T_events[k]

        times_segments = LinRange(t_past, t_e, t_s)
        dt = times_segments[2] - times_segments[1]
        for i ∈ 1:t_s -1
            Φ = exp_step(A_int, Φ, dt, times_segments[i])
        end

        # Exact states at the event boundary
        idx = sol.event_indices[k]
        x⁻ = sol.x[idx]   # pre-impact state (indexing is different from General)

        # Eval augmented diff matrices at the event
        Δ_augmented = Δ_augmented_f(x⁻)
        Φ = Δ_augmented * Φ
        t_past = t_e
    end

    # Determine the tail
    t_current = isempty(T_events) ? sol.t[1] : T_events[end]
    t_final = sol.t[end]
    times_segment = LinRange(t_current, t_final, t_s)
    dt = times_segment[2] - times_segment[1]
    for i ∈ 1:t_s -1
        Φ = exp_step(A_int, Φ, dt, times_segment[i])
    end
    # The STM should always have determinate 1 (if elastic)
    if abs(det(Φ)-1) > 0.1 && sys.e == 1.0
        @warn "The STM does not have determinate 1."
    end
    return Φ
end

## Tangent dynamics for a Filippov system
"""
    tangent_dynamics(sol::FilippovSol, sys::FilippovSystem; 
            Df = nothing, Dg = nothing, 
            guard_deriv = nothing, t_s = 1000)

Compute the state transition matrix across a trajectory from a General system.

"""
function tangent_dynamics(sol::FilippovSol, sys::FilippovSystem; Df = nothing, Dg = nothing, res_deriv = nothing, guard_deriv = nothing, t_s = 1000)
    # Extract out the event times
    T_events = sol.event_times

    # Depending on the mode we're in, we have three different variations
    # Mode 'f'
    A_f(t) = isnothing(Df) ? ForwardDiff.jacobian(y->sys.F(y), sol(t)) : Df(sol(t))
    # Mode 'g'
    A_g(t) = isnothing(Dg) ? ForwardDiff.jacobian(y->sys.G(y), sol(t)) : Dg(sol(t))
    # Mode 'k'. This version is dependent on ForwardDiff
    dh(y) = isnothing(guard_deriv) ? ForwardDiff.gradient(sys.h, y) : guard_deriv(y)
    λ(y) = dot(sys.G(y), dh(y)) / (dot(sys.G(y), dh(y)) - dot(sys.F(y), dh(y)))
    k_vf(y) = λ(y) * sys.F(y) + (1-λ(y)) * sys.G(y)
    A_k(t) = ForwardDiff.jacobian(y->k_vf(y), sol(t))
    
    # Initialize the variational matrix and time
    Φ = Matrix(I(size(A_f(sol.t[end]))[1]))
    t_past = sol.t[1]

    # Loop through events using their index to get pre/post states (these should be equal for Filippov)
    for k in 1:length(T_events)
        # The next event
        t_e = T_events[k]

        # The time timesegment
        times_segments = LinRange(t_past, t_e, t_s)
        dt = times_segments[2] - times_segments[1]

        # What mode are we in?
        idx = sol.event_indices[k]
        if idx == 1
            current_mode = sol.mode[idx]
        else
            current_mode = sol.mode[idx-1]
        end
        
        # Integrate the corresponding trajectory
        if current_mode == :f
            for i ∈ 1:t_s -1
                Φ = exp_step(A_f, Φ, dt, times_segments[i])
            end
        elseif current_mode == :g
            for i ∈ 1:t_s -1
                Φ = exp_step(A_g, Φ, dt, times_segments[i])
            end
        elseif current_mode == :k
            for i ∈ 1:t_s -1
                Φ = exp_step(A_k, Φ, dt, times_segments[i])
            end
        else
            error("Unknown mode type")
        end

        # Extract states/type at the event boundary
        x_hit = sol.x[idx]
        cross_type = sol.event_types[k]

        # Determine the correct differential
        dh_val = isnothing(guard_deriv) ? ForwardDiff.gradient(sys.h, x_hit) : guard_deriv(x_hit)
        if cross_type == :f_to_g
            Δ_augmented = I + 1/dot(dh_val, sys.F(x_hit))*(sys.G(x_hit)-sys.F(x_hit)) * dh_val'
        elseif cross_type == :f_to_k
            Δ_augmented = I + 1/dot(dh_val, sys.F(x_hit))*(k_vf(x_hit)-sys.F(x_hit)) * dh_val'
        elseif cross_type == :g_to_f
            Δ_augmented = I + 1/dot(dh_val, sys.G(x_hit))*(sys.F(x_hit)-sys.G(x_hit)) * dh_val'
        elseif cross_type == :g_to_k
            Δ_augmented = I + 1/dot(dh_val, sys.G(x_hit))*(k_vf(x_hit)-sys.G(x_hit)) * dh_val'
        else
            Δ_augmented = I
        end

        # Apply the transition and update time
        Φ = Δ_augmented * Φ
        t_past = t_e

    end

    # Determine the tail
    t_current = isempty(T_events) ? sol.t[1] : T_events[end]
    t_final = sol.t[end]
    times_segment = LinRange(t_current, t_final, t_s)
    dt = times_segment[2] - times_segment[1]
    final_mode = sol.mode[end]
    # Integrate
    if final_mode == :f
        for i ∈ 1:t_s -1
            Φ = exp_step(A_f, Φ, dt, times_segment[i])
        end
    elseif final_mode == :g
        for i ∈ 1:t_s -1
            Φ = exp_step(A_g, Φ, dt, times_segment[i])
        end
    elseif final_mode == :k
        for i ∈ 1:t_s -1
            Φ = exp_step(A_k, Φ, dt, times_segment[i])
        end
    else
        error("Unknown mode type")
    end

    return Φ
end

## The list of Lyapunov exponents
# Inputs: The system and initial condition
#         The rest of the inputs are performance parameters
# Outputs: The vector of Lyapunov exponents
"""
    LyapunovExponents(sys::Union{GeneralSystem, FilippovSystem, MechanicalSystem}, x0::AbstractArray; 
                solver::AbstractODESolver=RK45(),
                Df = nothing, Dg = nothing, res_deriv = nothing, 
                guard_deriv = nothing, t_s = 1000, 
                run_length::AbstractFloat=1.0, 
                run_iter::Int=Int(1e4), transient::AbstractFloat=0.0)

Compute the Lyapunov exponents for a trajectory from a General or Filippov system.

"""
function LyapunovExponents(sys::Union{GeneralSystem, FilippovSystem, MechanicalSystem}, x0::AbstractArray; 
    solver::AbstractODESolver=RK45(), tol = 1e-6, dt_initial=1e-3,
    Df = nothing, Dg = nothing, res_deriv = nothing, guard_deriv = nothing, 
    t_s = 1000, run_length::AbstractFloat=1.0, 
    run_iter::Int=Int(1e4), transient::AbstractFloat=0.0)

    # Discard the transient
    if transient > 0.0
        prob_transient = prob(sys, x0, (0.0, transient))
        sol_transient  = solve(prob_transient, solver, tol=tol, dt_initial=dt_initial)
        # If the solution cannot reach the end, we terminate and declare NaN exponents.
        if sol_transient.t[end] < transient
            return fill(NaN, length(x0))
        end
        x0 = sol_transient(transient)
    end
    # Loop through solutions and collect the exponents
    Φ = I(length(x0))
    λ = zeros(length(x0))
    current_time = 0.0
    for n ∈ 1:run_iter
        # Collect the STM for each run
        prob_iter = prob(sys, x0, (current_time, current_time+run_length))
        sol_iter  = solve(prob_iter, solver, tol=tol, dt_initial=dt_initial)
        # Exit early if solution doesn't exist
        if sol_iter.t[end] < current_time+run_length
            return fill(NaN, length(x0))
        end
        current_time = current_time + run_length
        x0 = sol_iter(sol_iter.t[end])
        Φ_iter = tangent_dynamics(sol_iter, sys; Df, Dg, res_deriv, guard_deriv, t_s)
        Φ = Φ_iter * Φ
        # The QR decomposition
        F = qr(Φ)
        Q, R = Matrix(F.Q), Matrix(F.R)
        D = diag(R)
        # Fix signs
        signs = ifelse.(D .>= 0, 1.0, -1.0)
        Q = Q * Diagonal(signs)
        R = Diagonal(signs) * R
        Φ = Q
        # Record exponents
        λ .+= log.(diag(R))
    end
    return λ ./ (run_length * run_iter)
end