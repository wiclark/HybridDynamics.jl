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
function filippov_vector_field(sys::FilippovSystem, x; Ftol=1e-4, atol=1e-8)
    F = sys.F
    G = sys.G
    h_val = sys.h(x)
    N_val = sys.N

    # Away from the guard
    if h_val > Ftol
        return y -> F(y), false
    elseif h_val < -Ftol
        return y -> G(y), false
    end

    # Near the guard (Sliding mode check)
    a(y) = dot(N_val(y), F(y))
    b(y) = dot(N_val(y), G(y))

    if (a(x) < -atol && b(x) > atol) || (a(x) > atol && b(x) < -atol)
        λ(y) = a(y) / (a(y) - b(y))
        return y -> (1 - λ(y)) * F(y) + λ(y) * G(y), true

    elseif a(x) > atol && b(x) > atol
        return y -> F(y), false

    elseif a(x) < -atol && b(x) < -atol
        return y -> G(y), false
    else
        @warn "Failed vector field determination at a=$(a(x)), b=$(b(x)). Defaulting to F."
        return y -> F(y), false
    end
end


function take_step_filippov!(solver, prob::prob{S,I,T}, Δt, tol, sol; 
    dense_out=true, stepper::AbstractODESolver=RK4(), 
    event_method::AbstractEventLocator=LinearLocator(), guard_direction=0) where {S<:FilippovSystem, I, T}

    xₖ = sol.x[end]
    tₖ = sol.t[end]
    sys = prob.sys

    vf_fun, sliding_now = filippov_vector_field(sys, xₖ)
    vf(x, t) = vf_fun(x)

    # Check if we are currently sitting ON the guard
    h_val = guard(sys, xₖ)
    is_on_guard = !isnothing(h_val) && abs(h_val) <= 1e-4

    # If sliding OR starting the step directly on the guard, 
    # pass the dummy system to prevent false zero-crossing detections.
    if sliding_now || is_on_guard
        dummy_sys = FilippovSystem(sys.F, sys.G, x -> 1.0, sys.N)
        active_prob = HybridDynamics.prob(dummy_sys, prob.init, prob.tspan)
    else
        active_prob = prob
    end

    x_predict, eventtrigger, _, dt_used, dt_next = take_step(solver, active_prob, vf, xₖ, tₖ, Δt, tol, sol; guard_direction=guard_direction)

    if eventtrigger
        t_star, x_star = locate_event(event_method, prob, solver, vf, xₖ, tₖ, Δt, guard(sys, xₖ), tol, sol, stepper)

        if abs(guard(sys, x_star)) > 1e-3
            @warn "Event Location is not located on the guard."
        end

        push!(sol.event_times, t_star)
        push!(sol.t, t_star)
        push!(sol.x, x_star)
        push!(sol.event_indices, length(sol.t))
        
        if dense_out
            push!(sol.dx, vf(x_star, t_star))
        end

    else 
        push!(sol.t, tₖ + dt_used)
        push!(sol.x, x_predict)

        if dense_out
            push!(sol.dx, vf(x_predict, tₖ + dt_used))
        end
    end

    return x_predict, dt_used, dt_next, false, sliding_now
end

function solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45();
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true,
               dt_initial=0.01, dt_min = 1e-6, max_iter = 10^6,
               tol = 1e-6,
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

        _, _, Δt, terminate, sliding_now = take_step_filippov!(solver, prob, Δt, tol, sol; dense_out=dense_out, stepper=stepper, event_method=event_method, guard_direction=guard_direction)

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

#=
function guard(sys::FilippovSystem, x)
    if isnothing(sys.h)
        return nothing
    else
        return sys.h(x)
    end
    
end

# EXTERNAL
# Filippov-specific solve
function solve(prob::prob{S,I,T}, solver::AbstractODESolver=RK45();
    dense_out = true,
    dt_initial = 0.001, max_iter = 10^6, tol = 1e-6, 
    stepper::AbstractODESolver=RK4(), event_method::AbstractEventLocator=LinearLocator(),
    kwargs...) where {S<:FilippovSystem, I, T}
    
    sys = prob.sys
    sol = FilippovSol(prob)

    _, t_end = prob.tspan     # Extract start and end times for bounds

    Δt = dt_initial                 #Initialize current time step with user input
    iter = 0                        #Start iteration counter
    sliding_prev = false

    # Run sim until end of specified time span
    while sol.t[end] < t_end 
        
    # Safties
        # Stop if we hit the iteration limit to avoid memory doomsday
        iter += 1
        if iter > max_iter 
            @warn "Maximum Iteration Count ($max_iter) exceeded."
            break
        end

        # Terminate if the remaining time is below machine precision
        if t_end - sol.t[end] <= eps(t_end)
            break
        end

        #Truncate time step if we overshoot the final sim time
        dt_step = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

    # Actually solve now

        xₖ = sol.x[end] #Retrieve current state at start of step
        tₖ = sol.t[end] #Retrieve current time at start of step

        # Choose vector field for current step based on current position relative to the guard
        vf_fun, sliding_now = filippov_vector_field(sys, xₖ)

        ## This doesn't quite work due to tolerances in the above function
        if sliding_now && !sliding_prev
            @warn "Sliding mode entered at t = $(tₖ)"
        end

        sliding_prev = sliding_now
        if sliding_now
            push!(sol.s, tₖ)
        end

        vf(x,t) = vf_fun(x)
        x_predict, _, _, dt_used, dt_next = take_step(solver, prob, vf, xₖ, tₖ, dt_step, tol, sol)

        tₖ += dt_used
        xₖ = x_predict
        if dense_out
            push!(sol.dx, vf(xₖ, tₖ))
        end
        push!(sol.t, tₖ)
        push!(sol.x, xₖ)

        Δt = dt_next

    end

    return sol
end
=#