# A Filippov (discontinuous) dynamical system
struct FilippovSystem{F, G, H, N} <: AbstractHybridSystem
    F::F    # Function one, H(x) > 0
    G::G    # Function two, H(x) < 0
    h::H    # Guard
    N::N    # Normal to the guard, ∇H
end

# Default to auto diff to find ∇H
function FilippovSystem(F, G, H; N= (x-> ForwardDiff.gradient(H,x)))
    return FilippovSystem(F, G, H, N)
end

struct FilippovSol{T, X, DX, S, E, EI} <: AbstractHybridSolution
    t::T                # Time data
    x::X                # Position data
    dx::DX              # f(x) Derivative at each state x - only filled when dense_out = true
    s::S                # Time indices while sliding (this still needs added)
    event_times::E      # Times of where events take place
    event_indices::EI   # Indices of times where events take place
end

function FilippovSol(prob::prob{S, I, T}) where {S<:FilippovSystem, I, T}
    return FilippovSol([prob.tspan[1]],
        [prob.init],
        Vector{Vector{Float64}}(),
        Float64[],
        Float64[],
        Int64[])
end

function guard(sys::FilippovSystem, x::AbstractArray)
    x_phys = (x isa AbstractMatrix) ? x[:, 1] : x
    return isnothing(sys.h) ? nothing : sys.h(x_phys)
end

# INTERNAL
# Returns a vector field function at state `x` for a Filippov system
######
### WC: I do not believe that there is any reason for 'atol'
######
function filippov_vector_field(sys::FilippovSystem, x; Ftol=1e-7, atol=1e-7)
    F = sys.F
    G = sys.G
    h_val = sys.h(x)
    N_val = sys.N

    if h_val > Ftol
        return y -> F(y), false
    elseif h_val < -Ftol
        return y -> G(y), false
    end

    # Near the guard (Sliding mode check)
    a(y) = dot(N_val(y), F(y))
    b(y) = dot(N_val(y), G(y))

    if a(x) < -atol && b(x) > atol
        λ(y) = a(y) / (a(y) - b(y))
        return y -> (1 - λ(y)) * F(y) + λ(y) * G(y), true
    elseif a(x) < -atol && b(x) < -atol
        return y -> G(y), false
    elseif a(x) > atol && b(x) > atol
        return y -> F(y), false
    elseif a(x) > atol && b(x) < -atol
        return h_val >= 0.0 ? (y -> F(y), false) : (y -> G(y), false)
    else
        return y -> F(y), false
    end
end


function take_step_filippov!(solver, prob::prob{S,I,T}, Δt, tol, sol; 
    dense_out=true, stepper::AbstractODESolver=RK4(), 
    event_method::AbstractEventLocator=LinearLocator(), guard_direction=0, boundary_tol) where {S<:FilippovSystem, I, T}

    # Extract current sim state and time
    xₖ = sol.x[end]
    tₖ = sol.t[end]
    # Acces the system def from the problem
    sys = prob.sys

    # Define tolerances. guard_tol is for precise logic, boundary_layer for sliding stuff. I am sure we can add more fine tuning here but it works fine now
    guard_tol = max(tol * 10, 1e-7)
    boundary_layer = guard_tol * boundary_tol

    # Retrieve the active Filippov vector field and determine if we are sliding or not
    vf_fun, sliding_now = filippov_vector_field(sys, xₖ; Ftol=boundary_layer, atol=guard_tol)
    # Wrap the vector field for use in ODE solvers
    vf(x, t) = vf_fun(x)

    # Determine if we are currently on the switching surface (guard)
    h_val = guard(sys, xₖ)
    is_on_guard = !isnothing(h_val) && abs(h_val) <= guard_tol

    # When sliding or on guard, we mask the guard logic to prevent the ODE solvers from triggering events
    # due to minor numerical bugs around H(x) = 0
    if sliding_now || is_on_guard
        # Create a dummy system where the guard is effectively invisible.
        dummy_sys = FilippovSystem(sys.F, sys.G, x -> 1.0, sys.N)
        active_prob = HybridDynamics.prob(dummy_sys, prob.init, prob.tspan)
        #Use simple stepping for on guard states to prevent adaptive methods breaking things. 
        active_solver = (is_on_guard && !sliding_now) ? stepper : solver
    else
        # Otherwise, use standard system and solver. 
        active_prob = prob
        active_solver = solver
    end

    # Delegate the numerical integration to the active solver
    x_predict, eventtrigger, t_root, dt_used, dt_next = take_step(active_solver, active_prob, vf, xₖ, tₖ, Δt, tol, sol; guard_direction=guard_direction)

    # If a solver reports a very small step while sliding, it is likely caught
    # in a rejection cycle. We force a conservative step to keep things moving. 
    if (sliding_now || is_on_guard) && dt_used < 1e-10
        dt_used = max(Δt * 0.1, 1e-5) #Force larger step
        x_predict = xₖ .+ dt_used .* vf(xₖ, tₖ) #Estimate state via prjection
        eventtrigger = false    #Clear event trigger to break the loop
        t_root = NaN
    end

    # Drifting occurs when our solver pulls the state of the switching surface (H(x) = 0).
    # We use Newton projection step along the surface normal sys.N to snap the state
    # back onto the manifold, preventing the solver from leaking out the sliding mode. 
    if sliding_now && !eventtrigger
        h_err = guard(sys, x_predict)
        if !isnothing(h_err)
            n_vec = sys.N(x_predict)
            # Apply correction
            x_predict = x_predict .- (h_err / dot(n_vec, n_vec)) .* n_vec
        end
    end
    
    # Dampen the next step size if we are sliding to prioritize surface accuracy 
    if sliding_now || is_on_guard
        dt_next = min(dt_next, dt_used * 1.5, 0.01) 
    end

    # Handle event detection
    if eventtrigger
        t_star, x_star = locate_event(event_method, prob, solver, vf, xₖ, tₖ, dt_used, guard(sys, xₖ), tol, sol, stepper)
        push!(sol.event_times, t_star)
        push!(sol.t, t_star)
        push!(sol.x, x_star)
        push!(sol.event_indices, length(sol.t))
        if dense_out; push!(sol.dx, vf(x_star, t_star)); end
    else 
        push!(sol.t, tₖ + dt_used)
        push!(sol.x, x_predict)
        if dense_out; push!(sol.dx, vf(x_predict, tₖ + dt_used)); end
    end
    return x_predict, dt_used, dt_next, false, sliding_now
end

function solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45();
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial=0.01, dt_min = 1e-6, max_iter = 10^6,
               tol = 1e-6, boundary_tol = 10,
               stepper::AbstractODESolver=RK4(),
               guard_direction = 0, 
               ) where {S<:FilippovSystem, I, T}
    sys = prob.sys
    sol = FilippovSol(prob)
    _, t_end = prob.tspan
    Δt = dt_initial
    iter = 0
    sliding_prev = false

    while sol.t[end] < t_end
        iter += 1
        if iter > max_iter
            @info "Maximum Iterations $max_iter Reached."
            break
        end

        if t_end - sol.t[end] <= eps(t_end)
            @info "Time step below minimum threshold $dt_min. Terminating."
            break
        end

        Δt = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

        _, _, Δt, terminate, sliding_now = take_step_filippov!(solver, prob, Δt, tol, sol; dense_out=dense_out, stepper=stepper, event_method=event_method, guard_direction=guard_direction, boundary_tol)

        #handle sliding tracking 
        if sliding_now && !sliding_prev
            @info "Sliding mode entered at t = $(sol.t[end])"
        end

        if sliding_now
            push!(sol.s, sol.t[end])
        end
        sliding_prev = sliding_now
        if terminate
            break
        end
    end
    return sol
end