# A Filippov (discontinuous) dynamical system
struct FilippovSystem{F, G, H, N} <: AbstractHybridSystem
    F::F    # Function one, H(x) > 0
    G::G    # Function two, H(x) < 0
    h::H    # Guard
    N::N    # Normal to the guard, ∇H
end

# Default to auto diff to find ∇H
"""
    FilippovSystem(F, G, H; N= (x-> ForwardDiff.gradient(H,x)))

Construct a Filippov system of the form:
```math
	\\begin{cases}
        \\dot{x} = F(x), & H(x) > 0 \\
        \\dot{x} = G(x), & H(x) < 0
    \\end{cases}
```
"""
function FilippovSystem(F, G, H; N= (x-> ForwardDiff.gradient(H,x)))
    return FilippovSystem(F, G, H, N)
end

struct FilippovSol{T, X, DX, S, E, EI, ET, MO} <: AbstractHybridSolution
    t::T                # Time data
    x::X                # Position data
    dx::DX              # f(x) Derivative at each state x - only filled when dense_out = true
    s::S                # Time indices while sliding (this still needs added)
    event_times::E      # Times of where events take place
    event_indices::EI   # Indices of times where events take place
    event_types::ET     # The type of transition that occurred 
    mode::MO            # The mode the dynamics are in (a symbol)
end

function FilippovSol(prob::prob{S, I, T}) where {S<:FilippovSystem, I, T}
    return FilippovSol([prob.tspan[1]],
        [prob.init],
        Vector{Vector{Float64}}(),
        Float64[],
        Float64[],
        Int64[],
        Symbol[], 
        Symbol[])       
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
        return y -> F(y), :f
    elseif h_val < -Ftol
        return y -> G(y), :g
    end

    # Near the guard (Sliding mode check)
    a(y) = dot(N_val(y), F(y))
    b(y) = dot(N_val(y), G(y))
    
    # Evaluate at current state x for routing logic
    a_x = a(x)
    b_x = b(x)

    if a_x < -atol && b_x > atol
        λ(y) = a(y) / (a(y) - b(y))
        return y -> (1 - λ(y)) * F(y) + λ(y) * G(y), :k
    elseif a_x < atol && b_x < atol
        # Both point down (or one is tangent), enter H < 0 safely
        return y -> G(y), :g
    elseif a_x > -atol && b_x > -atol
        # Both point up (or one is tangent), enter H > 0 safely
        return y -> F(y), :f
    else
        # Repulsive boundary fallback
        return h_val >= 0.0 ? (y -> F(y), :f) : (y -> G(y), :g)
    end
end


function take_step_filippov!(solver, prob::prob{S,I,T}, Df, Δt, tol, sol; 
    dense_out=true, stepper::AbstractODESolver=RK4(), 
    event_method::AbstractEventLocator=LinearLocator(), guard_direction=0, boundary_tol, track_sliding) where {S<:FilippovSystem, I, T}

    # Extract current sim state and time
    xₖ = sol.x[end]
    tₖ = sol.t[end]
    # Acces the system def from the problem
    sys = prob.sys

    # Define tolerances. guard_tol is for precise logic, boundary_layer for sliding stuff. I am sure we can add more fine tuning here but it works fine now
    guard_tol = max(tol * 10, 1e-7)
    boundary_layer = guard_tol * boundary_tol

    # Retrieve the active Filippov vector field and determine if we are sliding or not
    vf_fun, current_mode = filippov_vector_field(sys, xₖ; Ftol=boundary_layer, atol=guard_tol)

    if current_mode == :k
        sliding_now = true
    else
        sliding_now = false
    end

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
        Δt = min(Δt, 0.01)
    else
        # Otherwise, use standard system and solver. 
        active_prob = prob
        active_solver = solver
    end

    # Delegate the numerical integration to the active solver
    x_predict, eventtrigger, t_root, dt_used, dt_next = take_step(active_solver, active_prob, vf, Df, xₖ, tₖ, Δt, tol, sol; guard_direction=guard_direction)

    # If a solver reports a very small step while sliding, it is likely caught
    # in a rejection cycle. We force a conservative step to keep things moving. 
    if (sliding_now || is_on_guard) && dt_used < 1e-10
        dt_used = max(Δt * 0.1, 1e-5) #Force larger step
        x_predict = xₖ .+ dt_used .* vf(xₖ, tₖ) #Estimate state via prjection
        eventtrigger = false    #Clear event trigger to break the loop
        t_root = NaN
    end

    # sliding exit detection
    sliding_exit_trigger = false
    exit_guard_fun = nothing

    # Drifting occurs when our solver pulls the state of the switching surface (H(x) = 0).
    # We use Newton projection step along the surface normal sys.N to snap the state
    # back onto the manifold, preventing the solver from leaking out the sliding mode. 
    if sliding_now && !eventtrigger
        # Iterate Newton projection to achieve machine-precision manifold adherence
        for _ in 1:10
            h_err = guard(sys, x_predict)
            (isnothing(h_err) || abs(h_err) < 1e-11) && break
            n_vec = sys.N(x_predict)
            denom = dot(n_vec, n_vec)
            if denom > 1e-14
                x_predict = x_predict .- (h_err / denom) .* n_vec
            else
                break
            end
        end

        # check if we crossed an exit boundary during this sliding step
        a_k = dot(sys.N(xₖ), sys.F(xₖ))
        b_k = dot(sys.N(xₖ), sys.G(xₖ))

        a_p = dot(sys.N(x_predict), sys.F(x_predict))
        b_p = dot(sys.N(x_predict), sys.G(x_predict))

        # exit implies a(x) goes to 0 or b(x) goes to 0
        exit_a = (a_k < -guard_tol) && (a_p >= -guard_tol)
        exit_b = (b_k > guard_tol) && (b_p <= guard_tol)

        sliding_exit_trigger = exit_a || exit_b
        if sliding_exit_trigger
            # set guard for the locator to whichever vector field decoupled
            exit_guard_fun = exit_a ? (x -> dot(sys.N(x), sys.F(x))) : (x -> dot(sys.N(x), sys.G(x)))
        end
    end
    
    # Dampen the next step size if we are sliding to prioritize surface accuracy 
    if sliding_now || is_on_guard
        dt_next = min(dt_next, dt_used * 1.5, 0.01) 
    end

    # Handle event detection
    t_next = tₖ + dt_used
    x_next = x_predict
    mode_next = current_mode
    is_event = false
    event_type = nothing

    if eventtrigger
        t_star, x_star = locate_event(event_method, prob, solver, vf, Df, xₖ, tₖ, dt_used, guard(sys, xₖ), tol, sol, stepper)

        _, mode_next = filippov_vector_field(sys, x_star; Ftol=boundary_layer, atol=guard_tol)
        event_type = Symbol(current_mode, :_to_, mode_next) # e.g. :f_to_g or :f_to_k
        
        t_next, x_next = t_star, x_star
        is_event = true

    elseif sliding_exit_trigger
        # create dummy system targeting the sliding exit condition 
        exit_sys = FilippovSystem(sys.F, sys.G, exit_guard_fun, sys.N)
        exit_prob = HybridDynamics.prob(exit_sys, prob.init, prob.tspan)

        # locate exit
        t_star, x_star = locate_event(event_method, exit_prob, solver, vf, Df, xₖ, tₖ, dt_used, exit_guard_fun(xₖ), tol, sol, stepper)

        # Iterative Newton re-projection to main switch surface
        for _ in 1:10
            h_err_star = guard(sys, x_star)
            (isnothing(h_err_star) || abs(h_err_star) < 1e-11) && break
            n_vec_star = sys.N(x_star)
            denom_star = dot(n_vec_star, n_vec_star)
            if denom_star > 1e-14
                x_star = x_star .- (h_err_star / denom_star) .* n_vec_star
            else
                break
            end
        end
        
        _, mode_next = filippov_vector_field(sys, x_star; Ftol=boundary_layer, atol=guard_tol)
        event_type = Symbol(:k_to_, mode_next)
        
        t_next, x_next = t_star, x_star
        is_event = true
        
        dt_used = t_star - tₖ
        x_predict = x_star
    end

    # Unified pushing
    push!(sol.t, t_next)
    push!(sol.x, x_next)
    push!(sol.mode, mode_next)
    
    if is_event
        push!(sol.event_times, t_next)
        push!(sol.event_indices, length(sol.t))
        push!(sol.event_types, event_type)
    end

    if dense_out
        push!(sol.dx, vf(x_next, t_next))
    end
    
    final_sliding_state = sliding_exit_trigger ? false : sliding_now
    return x_predict, dt_used, dt_next, false, final_sliding_state
end

"""
    solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45();
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial=0.01, dt_min = 1e-6, max_iter = 10^6,
               tol = 1e-6, boundary_tol = 10,
               stepper::AbstractODESolver=RK4(),
               guard_direction = 0, 
               ) where {S<:FilippovSystem, I, T}

Solve a Filippov system.

"""
function solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45();
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial=0.01, dt_min = 1e-6, max_iter = 10^6,
               tol = 1e-6, boundary_tol = 10,
               stepper::AbstractODESolver=RK4(),
               guard_direction = 0,
               track_sliding=:none,
               Df = nothing
               ) where {S<:FilippovSystem, I, T}
    sys = prob.sys
    sol = FilippovSol(prob)
    _, t_end = prob.tspan
    Δt = dt_initial
    iter = 0
    sliding_prev = false

    if !(track_sliding in (:both, :enter, :exit, :none))
        throw(ArgumentError("track_sliding must be :both, :enter, :exit, or :none"))
    end

    # Evaluate initial conditions for dense_out and check if starting on the guard
    x₀ = sol.x[end]
    t₀ = sol.t[end]
    
    # Boundary and guard tolerances
    guard_tol = max(tol * 10, 1e-7)
    boundary_layer = guard_tol * boundary_tol
    
    h_val = guard(sys, x₀)
    vf_fun, initial_mode = filippov_vector_field(sys, x₀; Ftol=boundary_layer, atol=guard_tol)

    push!(sol.mode, initial_mode)

    # Check if we start exactly on the switching surface
    if !isnothing(h_val) && abs(h_val) <= guard_tol
        # This really just means that we can start out sliding. This shouldn't be an issue.
        # @info "System started on the guard at t = $t₀."
        push!(sol.event_times, t₀)
        push!(sol.event_indices, length(sol.t))
        push!(sol.event_types, Symbol(:start_on_, initial_mode))
        
        if initial_mode == :k 
            if track_sliding in (:both, :enter)
                @info "Sliding mode entered at t = $t₀"
            end
            push!(sol.s, t₀)
            sliding_prev = true 
        end
    end

    # Initial derivative if dense_out is requested
    if dense_out
        push!(sol.dx, vf_fun(x₀))
    end

    while sol.t[end] < t_end
        iter += 1
        if iter > max_iter
            @info "Maximum Iterations $max_iter Reached."
            break
        end

        if t_end - sol.t[end] <= eps(t_end)
            sol.t[end] = t_end
            @info "Time step below minimum threshold $dt_min. Snapping final time to t = $t_end."
            break
        end

        Δt = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

        _, _, Δt, terminate, sliding_now = take_step_filippov!(solver, prob, Df, Δt, tol, sol; dense_out=dense_out, stepper=stepper, event_method=event_method, guard_direction=guard_direction, boundary_tol=boundary_tol, track_sliding=track_sliding)

        if sliding_now && !sliding_prev
            if track_sliding in (:both, :enter)
                @info "Sliding mode entered at t = $(sol.t[end])"
            end
        elseif !sliding_now && sliding_prev
            if track_sliding in (:both, :exit)
                @info "Sliding mode exited at t = $(sol.t[end])"
            end
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