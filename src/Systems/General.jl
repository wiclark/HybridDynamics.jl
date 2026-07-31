# A general hybrid dynamical system
struct GeneralSystem <: AbstractHybridSystem
    f::Function     #Continuous Dynamics: (x,t) -> dx/dt
    h::Function     #Guard Surface: x -> real
    Δ::Function     #Reset map: x -> x⁺
    direction::Int  #Direction for the guard
end

#Constructor for above
"""
    GeneralSystem(f, h, Δ; direction::Int=0)

Construct a general hybrid system of the form:


"""
function GeneralSystem(f, h, Δ; direction::Int=0)
    return GeneralSystem(f, h, Δ, direction)
end

struct GeneralSol{T, X, DX} <: AbstractHybridSolution
    t::T                        #Time points of sim
    x::X                        #State traj: T is generic to support varying state types
    dx::DX 
    event_times::Vector{Float64} #Explicit storage of timestamps where resets occurred
    event_indices::Vector{Int}   #Map of event_times to indices in the x and t vectors
end

#Internal - to initialize the solution struct
function GeneralSol(prob::prob{S, I, T}) where {S<:GeneralSystem, I, T}
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

function take_step_general!(solver, prob::prob{S,I,T}, f, Δt, tol, sol; dense_out=true, stepper::AbstractODESolver=RK4(), event_method::AbstractEventLocator=LinearLocator(), guard_direction=default_guard_direction(prob.sys)) where {S<:GeneralSystem, I, T}

    xₖ = sol.x[end]
    tₖ = sol.t[end]

    sys = prob.sys

    x_predict, eventtrigger, _, dt_used, dt_next = take_step(solver, prob, f, xₖ, tₖ, Δt, tol, sol; guard_direction=guard_direction)

    if eventtrigger

        t_star, x_star = locate_event(event_method, prob, solver, f, xₖ, tₖ, Δt, guard(sys, xₖ), tol, sol, stepper)

        if abs(guard(sys, x_star)) > 1e-3
            @warn "Event Location is not located on the guard."
        end

        x⁺ = apply_reset(sys, x_star)

        if guard_direction == 0 && abs(guard(sys, x⁺)) > tol
            @warn "An Instantaneous reset is detected. Possible beating/blocking/Zeno."
            return x⁺, dt_used, dt_next, true
        end

        push!(sol.event_times, t_star)

        push!(sol.t, t_star)
        push!(sol.x, x_star)

        push!(sol.event_indices, length(sol.t))

        push!(sol.t, t_star)
        push!(sol.x, x⁺)

        if dense_out
            push!(sol.dx, f(x_star, t_star))
            push!(sol.dx, f(x⁺, t_star))
        end

        if length(sol.event_times) > 1 && sol.event_times[end]-sol.event_times[end-1] < 1e-3
            @warn "Fast switching detected. Possible beating/blocking/Zeno."
            return x⁺, dt_used, dt_next, true
        end
    else
        push!(sol.t, tₖ+dt_used)
        push!(sol.x, x_predict)

        if dense_out
            push!(sol.dx, f(x_predict, tₖ+dt_used))
        end
    end
    return x_predict, dt_used, dt_next, false
end

"""
    solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45(); kwargs...) where {S<:GeneralSystem, I, T}

Solve a general hybrid system.

# Arguments
## Required:
- 'prob': The problem definition containing the system dynamics 'sys', initial state, and time span.
- 'solver': The numerical integration method used for continous steps. Defaults to RK45().

## Optional:

### Simulation and Step Controls:
* 'dt_initial' (Float64, default '0.01'): The starting time step for the continuous solver.
* 'dt_min' (Float64, default '1e-6'): The absolute minimum allowable time step. If the solver tries to go below this, the simulation terminates.
* 'max_iter' (Int, default '10^6'): The maximum number of continuous integration steps allowed before forcing a timeout.
- 'tol' (Float64, default '1e-6'): The baseline numerical tolerance used across the solver. Acts as the foundational unit for multipliers below.

### Event Handling
* 'event_method' (AbstractEventLocator, default 'LinearLocator()'): The algorithm used to pinpoint the exact time and state of a guard crossing.
* 'stepper' (AbstractODESolver. default 'RK4()'): The secondary ODE solver used internally by the event detection locator to pinpoint the impact state.
* 'guard\\_direction' (Int, default 'default\\_guard\\_direction(prob.sys)'): 0 -> detects crossings in either direction. 1 -> detects increasing crossings. -1 -> detects decreating crossings.

"""
function solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45();
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial=0.01, dt_min = 1e-6, max_iter = 10^6,
               tol = 1e-6,
               stepper::AbstractODESolver=RK4(),
               guard_direction = prob.sys.direction
               ) where {S<:GeneralSystem, I, T}

    sys = prob.sys
    f = sys.f
    sol = GeneralSol(prob)
    _, t_end = prob.tspan
    Δt = dt_initial
    iter = 0

    # Initial derivative, pray it's not on the guard (for now)
    if dense_out
        push!(sol.dx, f(prob.init, prob.tspan[1]))
    end

    while sol.t[end] < t_end
        
        iter += 1

        if iter > max_iter 
            @info "Maximum iterations $max_iter reached."
            break
        end

        if t_end - sol.t[end] <= eps(t_end)
            @info "Time step below minimum threshold $dt_min. Terminating"
            break 
        end

        Δt = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

        _, _, Δt, terminate = take_step_general!(solver, prob, f, Δt, tol, sol; dense_out=dense_out, stepper=stepper, event_method=event_method, guard_direction=guard_direction)

        if terminate
            break
        end

    end
    return sol
end