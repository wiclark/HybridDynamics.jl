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
function tangent_dynamics(sol::GeneralSol, sys::GeneralSystem)
    # Extract out the event times
    T_events = sol.event_times
    # Create the tangent vector field
    A(t) = ForwardDiff.jacobian(y->sys.f(y,t), sol(t))
    # As well as the augmented differential
    # ****sol(impact_time) returns the state immedeately post-impact****
    #     I want a way to get the state immedeately pre-impact.
    #     The code below is wrong (as it's using the post-impact values)
    ϵ=1e-5 # <- I want this gone
    Δ_star(t) = ForwardDiff.jacobian(y->sys.Δ(y), sol(t-ϵ))
    f⁻(t) = sys.f(sol(t-ϵ), t)
    f⁺(t) = sys.f(sys.Δ(sol(t-ϵ)), t)
    dh⁻(t) = ForwardDiff.gradient(sys.h, sol(t-ϵ))
    Δ_augmented(t) = Δ_star(t)*(I-(f⁻(t)*dh⁻(t)')./dot(dh⁻(t), f⁻(t))) + (f⁺(t)*dh⁻(t)')./dot(dh⁻(t), f⁻(t))
    # We want to solve the temporally-triggered system
    #   Φ' = A(t)Φ,           t ̸∈ T_events
    #   Φ⁺ = Δ_augmented(t)Φ, t ∈ T_events
    # We want to return Φ(sol.t[end])
    # return A, Δ_augmented
    # The initial condition is the identity matrix
    Φ = Matrix( I(size(A(0))[1]) )
    # Let's just take 100 steps between each event
    t_past = sol.t[1] # <- The initial time
    pre_mature = false
    for t_e ∈ T_events
        if t_e - t_past < 1e-4 # <- Strange things happen if we allow faster switching
                               #    I think this has to do with ϵ above
            pre_mature = true
            println([t_past, t_e])
            break
        end
        times_segment = LinRange(t_past, t_e, 1_000)
        dt = times_segment[2] - times_segment[1]
        for i ∈ 1:999
            Φ = exp_step(A, Φ, dt, times_segment[i])
        end
        Φ = Δ_augmented(t_e) * Φ
        t_past = t_e
    end
    # Determine the tail. Notice that T_events may be empty
    if !pre_mature
        t_current = isempty(T_events) ? sol.t[1] : T_events[end]
        t_final = sol.t[end]
        times_segment = LinRange(t_current, t_final, 1_000)
        dt = times_segment[2] - times_segment[1]
        for i ∈ 1:999
            Φ = exp_step(A, Φ, dt, times_segment[i])
        end
    end
    return Φ, pre_mature
end

## The list of Lyapunov exponents
# Inputs: The system and initial condition
#         The rest of the inputs are performance parameters
# Outputs: The vector of Lyapunov exponents
function LyapunovExponents(sys::GeneralSystem, x0::AbstractArray; run_length::AbstractFloat=1.0, run_iter::Int=Int(1e4), transient::AbstractFloat=0.0)
    # Discard the transient
    if transient > 0.0
        prob_transient = prob(sys, x0, (0.0, transient))
        sol_transient  = solve(prob_transient)
        x0 = sol_transient(transient)
    end
    # Loop through solutions and collect the exponents
    Φ = I(length(x0))
    λ = zeros(length(x0))
    current_time = 0.0
    for n ∈ 1:run_iter
        # Collect the STM for each run
        prob_iter = prob(sys, x0, (current_time, current_time+run_length))
        sol_iter  = solve(prob_iter, RK45(), tol=1e-9)
        x0 = sol_iter(sol_iter.t[end])
        Φ_iter, pre_mature = tangent_dynamics(sol_iter, sys)
        if pre_mature
            @warn "Unable to complete trajectory"
        end
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

#=
#Function calcs the continuous time derivative for the augmented state. Runs the usual dx and dΦ.

#CURRENTLY THIS FUNCTION IS WORTHLESS BUT I WANT IT TO BE USEFUL LATER SO I WILL KEEP IT - DS
function variational_vector_field(f, U::AbstractMatrix, t)
    # Unpack state
    x = U[:, 1]
    Φ = U[:, 2:end]

    # Base dynamics: x' = f(x, p, t)
    dx = f(x, t)

    # Variational dynamics: Φ' = A(t)Φ
    A = ForwardDiff.jacobian(y -> f(y, t), x)
    dΦ = A * Φ

    # Return augmented derivative
    return hcat(dx, dΦ)
end

#Calcs the Δ_*^f at a boundary. We look how how much the vf mismatches before and after the jump (f⁺/f⁻), how much it scales things (DΔ⁻), and how the trajectory hits boundary (dh⁻).
#Then we create a matrix that takes all of this.
function compute_pushforward(f, Δ, h_guard, x⁻, t)
    n = length(x⁻)
    Id = I(n)

    # Eval field at boundaries (using p for parameters)
    f⁻ = f(x⁻, t)
    x⁺ = Δ(x⁻, t)
    f⁺ = f(x⁺, t)

    # Compute grads and jacob via ForwardDiff
    dh⁻ = ForwardDiff.gradient(h_guard, x⁻)
    DΔ⁻ = ForwardDiff.jacobian(y -> Δ(y, t), x⁻)

    # Check dh(x) * f(x) = 0
    denom = dot(dh⁻, f⁻)
    if abs(denom) < 1e-6
        @warn "Non-transversal crossing detected: Trajectory is tangent to guard surface."
    end

    # Outer prods
    term1 = Id - (f⁻ * dh⁻') ./ denom
    term2 = (f⁺ * dh⁻') ./ denom

    # Full pushforward Δᶠ_*
    Δ_star_f = DΔ⁻ * term1 + term2

    return Δ_star_f
end

#When our solver hits h(x)=0 this is called. We get the pushforward matrix and multiply by Φ⁻ to get new Φ⁺. 
function apply_variational_jump(U::AbstractMatrix, f, Δ, h_guard, t)
    x⁻ = U[:, 1]
    Φ⁻ = U[:, 2:end]

    # Compute the pushforward before state updates
    Δ_star_f = compute_pushforward(f, Δ, h_guard, x⁻, t)

    # Apply disc jump to base state x⁺ = Δ(x⁻)
    x⁺ = Δ(x⁻, t)

    # Apply pf mapping to fund matrix: Φ⁺ = Δ_*^f * Φ⁻
    Φ⁺ = Δ_star_f * Φ⁻

    # Update state vector in-place
    U[:, 1] .= x⁺
    U[:, 2:end] .= Φ⁺
    
    return U
end
=#