# A general hybrid dynamical system
struct GeneralSystem <: AbstractHybridSystem
    f::Function     #Continuous Dynamics: (x,t) -> dx/dt
    h::Function     #Guard Surface: x -> real
    Δ::Function     #Reset map: x -> x⁺
    direction::Int  #Direction for the guard
end

#Constructor for above
"""
    GeneralSystem(f::Function, h::Function, Δ::Function; direction::Int=0)

Construct a general hybrid system of the form:
```math
    \\begin{cases}
        \\dot{x} = f(x, t), & h(x) \\ne 0 \\\\
        x^+ = \\Delta(x), & h(x) = 0
    \\end{cases}
```
"""
function GeneralSystem(f, h, Δ; direction::Int=0)
    return GeneralSystem(f, h, Δ, direction)
end

struct TimeTriggeredGeneral{F, V<:AbstractVector, D} <: AbstractHybridSystem
    f::F            # Continuous Dynamics
    event_times::V  # Times for jumps
    Δ::D            # Reset map x->x⁺
end

# Overload GeneralSystem constructor so the syntax Clark wanted works.
function GeneralSystem(f, times::AbstractVector, Δ)
    return TimeTriggeredGeneral(f, times, Δ)
end

# Dummy guard to make this work with our solvers with no headache - Maybe make it smarter later? - DS
function guard(sys::TimeTriggeredGeneral, x::AbstractArray)
    return 1.0 # Always positive, meaning it never crosses a 0-boundary
end

"""
    GeneralSol

Solution type returned by `solve` for a `GeneralSystem` or
`TimeTriggeredGeneral` system.

# Fields
- `t`: Vector of time values corresponding to the computed solution states.
- `x`: Vector of solution states.
- `dx`: Vector of derivatives evaluated at each solution state. Only populated when `dense_out=true`.
- `event_times`: Times at which reset events occur.
- `event_indices`: Indices into `t` and `x` corresponding to the post-reset states at each event.
"""
struct GeneralSol{T, X, DX} <: AbstractHybridSolution
    t::T                        #Time points of sim
    x::X                        #State traj: T is generic to support varying state types
    dx::DX 
    event_times::Vector{Float64} #Explicit storage of timestamps where resets occurred
    event_indices::Vector{Int}   #Map of event_times to indices in the x and t vectors
end

#Internal - to initialize the solution struct
function GeneralSol(prob::prob{S, I, T}) where {S<:Union{GeneralSystem, TimeTriggeredGeneral}, I, T}
    return GeneralSol([prob.tspan[1]], 
        [prob.init], 
        Vector{Vector{Float64}}(),
        Float64[],
        Int[])
end

#Internal
function guard(sys::GeneralSystem, x::AbstractArray)
    x_phys = (x isa AbstractMatrix) ? x[:, 1] : x
    val = sys.h(x_phys)
    # Return the most negative value (the most violated constraint)
    # If all are positive, it returns the smallest positive value.
    return val isa AbstractVector ? minimum(val) : val
end
#Internal
function apply_reset(sys::GeneralSystem, x::AbstractArray)
    return sys.Δ(x)
end

function take_step_general!(solver, prob::prob{S,I,T}, f, Df, Δt, tol, sol; 
    dense_out=true, stepper::AbstractODESolver=RK4(), 
    event_method::AbstractEventLocator, 
    guard_direction=default_guard_direction(prob.sys),
    event_max_iters, force_iters, guard_derivative, event_refinement_method) where {S<:GeneralSystem, I, T}

    xₖ = sol.x[end]
    tₖ = sol.t[end]

    sys = prob.sys

    x_predict, eventtrigger, _, dt_used, dt_next = take_step(solver, prob, f, Df, xₖ, tₖ, Δt, tol, sol; guard_direction=guard_direction, event_method=event_method)

    if eventtrigger
         
        t_star, x_star = locate_event(event_method, prob, stepper, f, Df, xₖ, tₖ, dt_used, guard(sys, xₖ), tol, sol, stepper; event_max_iters=event_max_iters, 
        force_iters=force_iters, guard_derivative=guard_derivative, event_refinement_method=event_refinement_method)

        if abs(guard(sys, x_star)) > 1e-3
            wrong_ness = guard(sys, x_star)
            @warn "Event Location is not located on the guard. h(x) = $wrong_ness at t = $t_star"
        end

        x⁺ = apply_reset(sys, x_star)

        if guard_direction == 0 && abs(guard(sys, x⁺)) > tol
            @warn "An Instantaneous reset is detected. Possible beating/blocking/Zeno."
            return x⁺, dt_used, dt_next, true
        end

        push!(sol.t, t_star, t_star)
        push!(sol.x, x_star, x⁺)
        
        push!(sol.event_times, t_star)
        push!(sol.event_indices, length(sol.t))

        if dense_out
            push!(sol.dx, f(x_star, t_star))
            push!(sol.dx, f(x⁺, t_star))
        end

        if length(sol.event_times) > 1 && sol.event_times[end]-sol.event_times[end-1] < 1e-3
            @warn "Fast switching detected. Possible beating/blocking/Zeno."
            return x⁺, dt_used, dt_next, true
        end
        return x⁺, dt_used, dt_next, false 
    else
        
        t_next = tₖ + dt_used

        if dense_out
            push!(sol.dx, f(x_predict, t_next)) 
        end

        push!(sol.t, t_next)
        push!(sol.x, x_predict)
    end
    return x_predict, dt_used, dt_next, false
end

"""
    solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45();
          event_method::AbstractEventLocator=LinearHermite(),
          dense_out=true,
          dt_initial=1e-3,
          dt_min=1e-6,
          max_iter=10^6,
          tol=1e-6,
          stepper::AbstractODESolver=RK4(),
          guard_direction=prob.sys.direction,
          Df=nothing,
          event_max_iters=100,
          force_iters=false,
          guard_derivative=nothing,
          event_refinement_method::AbstractEventRefinement=FalsePosition())
          where {S<:GeneralSystem, I, T}


The solver integrates the continuous dynamics until a guard crossing
is detected. The event location is then refined, and the reset map
is applied to obtain the post-event state.

# Arguments
- `prob`: Problem containing the `GeneralSystem`, initial condition, and time span.
- `solver`: ODE solver used for numerical integration. Defaults to `RK45()`.
- `event_method`: Event locator used to detect and locate guard crossings. Defaults to `LinearHermite()`.
- `dense_out`: If `true`, stores the derivative at each solution point in `sol.dx`. Defaults to `true`.
- `dt_initial`: Initial time step size. Defaults to `1e-3`.
- `dt_min`: Minimum allowable time step size. Defaults to `1e-6`.
- `max_iter`: Maximum number of integration iterations. Defaults to `10^6`.
- `tol`: Numerical tolerance used for integration and event detection. Defaults to `1e-6`.
- `stepper`: Fixed-step ODE solver used internally by the event locator to refine the event state. Defaults to `RK4()`.
- `guard_direction`: Direction of guard crossing to detect: `0` detects either direction, `1` detects rising crossings, and `-1` detects falling crossings. Defaults to `prob.sys.direction`.
- `Df`: Optional Jacobian or derivative information supplied to the numerical solver. Defaults to `nothing`.
- `event_max_iters`: Maximum number of iterations used by the event locator. Defaults to `100`.
- `force_iters`: If `true`, forces the event locator to perform the specified number of iterations rather than stopping early based on its convergence criteria. Defaults to `false`.
- `guard_derivative`: Optional function providing the derivative of the guard function, used by event refinement methods that require it. Defaults to `nothing`.
- `event_refinement_method`: Root-finding method used to refine the event location. Defaults to `FalsePosition()`.

# Returns
A `GeneralSol` containing the solution times, states, derivatives
(if `dense_out=true`), and detected events.
"""
function solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45();
               event_method::AbstractEventLocator=LinearHermite(),
               dense_out = true,
               dt_initial=1e-3, dt_min = 1e-6, max_iter = 10^6,
               tol = 1e-6,
               stepper::AbstractODESolver=RK4(),
               guard_direction = prob.sys.direction,
               Df = nothing,
               event_max_iters = 100,
               force_iters = false,
               guard_derivative = nothing,
               event_refinement_method::AbstractEventRefinement=FalsePosition()
               ) where {S<:GeneralSystem, I, T}

    sys = prob.sys
    f = sys.f
    sol = GeneralSol(prob)
    _, t_end = prob.tspan
    Δt = dt_initial
    iter = 0

    # Check if we start on the guard
    if abs(guard(sys, sol.x[end])) <= tol
        x₀ = sol.x[end]
        t₀ = sol.t[end]

        x⁺ = apply_reset(sys, x₀)

        push!(sol.t, t₀)
        push!(sol.x, x⁺)
        push!(sol.event_times, t₀)
        push!(sol.event_indices, length(sol.t))

        if dense_out
            push!(sol.dx, f(x₀, t₀))
            push!(sol.dx, f(x⁺, t₀))
        end

        @info "System started on the guard. Immediate jump applied at t = $t₀."
    else
        # Initial derivative if NOT on guard
        if dense_out
            push!(sol.dx, f(prob.init, prob.tspan[1]))
        end
    end

    while sol.t[end] < t_end
        
        iter += 1

        if iter > max_iter 
            @info "Maximum iterations $max_iter reached."
            break
        end

        if t_end - sol.t[end] < dt_min
            sol.t[end] = t_end  # Snap the final recorded time exactly to the end of tspan
            @info "Time to end of simulation below minimum time step. Snapping final time to t = $t_end"
            break 
        end

        dt_step = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

        _, _, Δt, terminate = take_step_general!(solver, prob, f, Df, dt_step, tol, sol; dense_out=dense_out, stepper=stepper, event_method=event_method, 
        guard_direction=guard_direction, event_max_iters=event_max_iters, force_iters=force_iters, guard_derivative=guard_derivative, event_refinement_method=event_refinement_method)

        if terminate
            break
        end

    end
    return sol
end

"""
    solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45();
          dense_out=true,
          dt_initial=1e-3,
          dt_min=1e-6,
          max_iter=10^6,
          tol=1e-6,
          Df=nothing)

The solver integrates the continuous dynamics until each specified
event time is reached. At each event time, the reset map is applied
to the state before integration continues.

Only event times inside the problem's time span are used.

# Arguments
- `prob`: Problem containing the `TimeTriggeredGeneral`, initial condition, and time span.
- `solver`: ODE solver used for numerical integration. Defaults to `RK45()`.
- `dense_out`: If `true`, stores the derivative at each solution point in `sol.dx`. Defaults to `true`.
- `dt_initial`: Initial time step size. Defaults to `1e-3`.
- `dt_min`: Minimum allowable time step size. Defaults to `1e-6`.
- `max_iter`: Maximum number of integration iterations. Defaults to `10^6`.
- `tol`: Numerical tolerance supplied to the numerical integration steps. Defaults to `1e-6`.
- `Df`: Optional Jacobian or derivative information supplied to the numerical solver. Defaults to `nothing`.

# Returns
A `GeneralSol` containing the solution times, states, derivatives
(if `dense_out=true`), and detected time-triggered events.
"""
function solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45();
               dense_out = true,
               dt_initial=1e-3, dt_min = 1e-6, max_iter = 10^6,
               tol = 1e-6,
               Df = nothing
               ) where {S<:TimeTriggeredGeneral, I, T}

    sys = prob.sys
    f = sys.f
    sol = GeneralSol(prob)
    t_start, t_end = prob.tspan
    Δt = dt_initial
    iter = 0

    # filter the specified times for events
    # orders and gets rid of any outside of tspan. Very useful!
    valid_times = filter(t -> t_start < t < t_end, sys.event_times)
    sort!(valid_times)
    targets = [valid_times; t_end]

    if dense_out
        push!(sol.dx, f(prob.init, t_start))
    end

    current_t = t_start
    current_x = prob.init

    # main loop
    for target_t in targets
        
        # Inner loop: continuous integration up to the next target
        while current_t < target_t - 1e-12
            iter += 1
            if iter > max_iter 
                @info "Maximum iterations $max_iter reached."
                return sol
            end

            # Do not overshoot the target time
            dt_step = min(Δt, target_t - current_t)

            x_predict, _, _, dt_used, dt_next = take_step(solver, prob, f, Df, current_x, current_t, dt_step, tol, sol; check=false)

            current_t += dt_used
            current_x = x_predict
            Δt = dt_next

            if dense_out
                push!(sol.dx, f(current_x, current_t))
            end

            push!(sol.t, current_t)
            push!(sol.x, current_x)
        end

        # Force exact alignment with the target time to prevent drift/infinite loops.
        # This became an issue when running certain tests. It may be gotten rid of later? - DS
        if current_t < target_t
            current_t = target_t
        end

        # Apply Jump
        if target_t != t_end
            x⁺ = sys.Δ(current_x)
            
            push!(sol.t, target_t)
            push!(sol.x, x⁺)
            push!(sol.event_times, target_t)
            push!(sol.event_indices, length(sol.t))

            if dense_out
                push!(sol.dx, f(current_x, target_t))
                push!(sol.dx, f(x⁺, target_t))
            end
            
            current_x = x⁺
        end
    end

    return sol
end