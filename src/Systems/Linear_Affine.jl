#EVERYTHING FOR LINEAR SYSTEMS 

#Linear System Structs: Goal to define the date for a linear system and its problem.
#We subtype it to AbstractHybridSystem/Problem for campatibility with solvers 

struct LinearSystem <: AbstractHybridSystem
    A::Matrix{Float64} #State Transition matrix (dx/dt = Ax)
    λ::Vector{Float64} #Normal Vector for the Guard Surface
    C::Matrix{Float64} #Reset map matrix (x⁺ = Cx)
    direction::Int
end

#External
#external constructor to help user see data types 
"""
    LinearSystem(A::AbstractMatrix, λ::AbstractVector, C::AbstractMatrix; direction::Int=0)

Construct a linear hybrid system of the form:
```math
	\\mathcal{LH} = 
	\\begin{cases}
		\\dot x = Ax, & λx \\neq 0, \\\\
		x^+ = Cx, & λx = 0.
	\\end{cases}
```
"""
function LinearSystem(A::AbstractMatrix, λ::AbstractVector, C::AbstractMatrix; direction::Int=0)
    return LinearSystem(Float64.(A), Float64.(λ), Float64.(C), direction)
end

# EVERYTHING FOR AFFINE SYSTEMS
#AFFINE SYSTEM STRUCTS
#similar to the Linear. Defining the concrete structure for affine systems.
struct AffineSystem <: AbstractHybridSystem
    A::Matrix{Float64}  #State Matrix
    b::Vector{Float64}  #Continuous affine vector (dx/dt = Ax + b)
    λ::Vector{Float64}  #Guard Normal vector  
    a::Float64          #Guard Offset const dot(λ, x) + a = 0
    C::Matrix{Float64}  #Reset matrix
    κ::Vector{Float64}  #Discrete affine vector const x⁺ = Cx + κ
    direction::Int
end

#External
# external constructor to help user see data types
"""
    AffineSystem(A::AbstractMatrix, b::AbstractVector, λ::AbstractVector, a::Real, C::AbstractMatrix, κ::AbstractVector; direction::Int=0)

Construct an affine hybrid system of the form:
```math
	\\mathcal{AH} = 
	\\begin{cases}
		\\dot x = Ax + b, & λx + a \\neq 0, \\\\
		x^+ = Cx + κ, & λx + a = 0.
	\\end{cases}
```
"""
function AffineSystem(A::AbstractMatrix, b::AbstractVector, λ::AbstractVector, a::Real, C::AbstractMatrix, κ::AbstractVector; direction::Int=0)
    return AffineSystem(Float64.(A), Float64.(b), Float64.(λ), Float64(a), Float64.(C), Float64.(κ), direction)
end

#SOLUTION STRUCTS AND HELPERS. 
"""
    LinearAffineSol

Solution type returned by `solve` for a `LinearSystem` or
`AffineSystem`.

# Fields
- `t`: Vector of time values corresponding to the computed solution states.
- `x`: Vector of solution states.
- `dx`: Vector of derivatives evaluated at each solution state. Only populated when `dense_out=true`.
- `event_times`: Times at which reset events occur.
- `event_indices`: Indices into `t` and `x` corresponding to the post-reset states at each event.
"""
struct LinearAffineSol{T, DX} <: AbstractHybridSolution
    t::Vector{Float64}
    x::Vector{T}
    dx::DX 
    event_times::Vector{Float64}
    event_indices::Vector{Int}
end 

function LinearAffineSol(prob::prob{S, I, T}) where {S<:Union{LinearSystem, AffineSystem}, I, T}
    return LinearAffineSol([prob.tspan[1]], [prob.init], Vector{Vector{Float64}}(), Float64[], Int[])
end

#Internal
#Calc guard function for linear systems
function guard(sys::LinearSystem, x::AbstractArray)
    #if x is a matrix (usually for Variational equation) we isolate the first column to calc the guard
    x_first_column = x isa AbstractMatrix ? x[:, 1] : x
    return sys.λ' * x_first_column
end
#Internal
#Calc guard for affine systems
function guard(sys::AffineSystem, x::AbstractArray)
    #if x is a matrix (usually for Variational equation) we isolate the first column to calc the guard
    x_first_column = x isa AbstractMatrix ? x[:, 1] : x
    return sys.λ' * x_first_column + sys.a 
end
#internal
function apply_reset(sys::LinearSystem, x::AbstractArray)
    return sys.C * x
end
#internal
function apply_reset(sys::AffineSystem, x::AbstractArray)
    #first linear transformation Cx
    x_new = sys.C * x
    # If matrix, κ only applies to the physical state (column 1)
    if x isa AbstractMatrix
        x_new[:, 1] .+= sys.κ
    else
        x_new .+= sys.κ
    end
    return x_new
end

#Solution Initialization
#Goal is to setup the empty memory containers before solver starts running. 
#pre-allocating the vectors with the exact starting conditions ensures that the solver loop is stable 

"""
    beating_and_blocking_sets(sys::Union{LinearSystem, AffineSystem})

Compute the sequence of constraint matrices and offset vectors
defining the beating sets and final blocking set of a linear or
affine hybrid system.

The beating sets describe the successive guard constraints generated
by repeated application of the inverse reset map. The computation
continues until the constraint matrix no longer increases in rank
or the maximum number of iterations is reached.

# Arguments
- `sys`: A `LinearSystem` or `AffineSystem` whose beating and blocking sets are to be computed.

# Returns
A named tuple containing:
- `beating_sets`: Vector of constraint matrices defining the successive beating sets. Each row represents a guard constraint.
- `beating_offsets`: Vector of offset vectors corresponding to the successive beating sets.
- `blocking_set`: Final constraint matrix defining the blocking set.
- `blocking_offsets`: Offset vector corresponding to the blocking set.
- `k_blocking`: Discrete iteration at which the constraint sequence reaches its final rank.

For linear systems, the offset vectors are zero. For affine systems,
the offsets account for the affine reset `κ`.
"""
function beating_and_blocking_sets(sys::Union{LinearSystem, AffineSystem})

    λ, C = sys.λ, sys.C
    n = length(λ)
    C_inv = inv(C)

    #Extracts affine properties to see if they exist. Tells us if we have linear or affine system
    a = hasproperty(sys, :a) ? sys.a : 0.0
    κ = hasproperty(sys, :κ) ? sys.κ : zeros(Float64, n)

    #Initialize containers for Σ_0 = {x | λᵀx = -a}
    O_current = reshape(λ, 1, :)
    b_current = [-a]

    O = [O_current]
    b_vecs = [b_current]

    prev_rank = rank(O_current)
    k_∞ = 0
    row = λ'

    for k in 1:n
        #Advance constraint row through inverse reset map
        row = row * C_inv
        new_row = reshape(row, 1, :)

        #offset update: b_k = b_{k-1} - row * κ
        new_offset = b_current[end] - (row * κ)[1]

        #next candidate matrix and offset 
        O_next = vcat(O_current, new_row)
        b_next = vcat(b_current, [new_offset])

        #Check for rank defi
        current_rank = rank(O_current)
        if current_rank == prev_rank
            k_∞ = k - 1
            break
        end

        #update loop and log constraints
        O_current = O_next
        b_current = b_next

        push!(O, O_current)
        push!(b_vecs, b_current)

        prev_rank = current_rank
        k_∞ = k
    end

        return (beating_sets = O, beating_offsets = b_vecs, blocking_set = O[end], blocking_offsets = b_vecs[end], k_blocking = k_∞)
end

#TRIVIALLY BLOCKING
"""
    is_trivially_blocking(sys::Union{LinearSystem, AffineSystem})

Determine whether the blocking set of a linear or affine hybrid system collapses to the origin, that is, whether `Σ_∞ = {0}`.

A system is considered trivially blocking when the final blocking constraint matrix has full rank and its corresponding offset vector
is approximately zero.

# Arguments
- `sys`: A `LinearSystem` or `AffineSystem` to be analyzed.

# Returns
A Boolean indicating whether the system is trivially blocking.
The offset vector is considered zero using an absolute tolerance
of `1e-10`.
"""
function is_trivially_blocking(sys::Union{LinearSystem, AffineSystem})
    analysis = beating_and_blocking_sets(sys)
    n = length(sys.λ)

    #check full rank AND the point is the oriign
    return rank(analysis.blocking_set) == n && isapprox(norm(analysis.blocking_offsets), 0.0, atol=1e-10)
end

function take_step_linear_affine!(solver, prob::prob{S, I, T}, f, Df, Δt, tol, sol; dense_out=true, stepper::AbstractODESolver=RK45(), 
    event_method::AbstractEventLocator, guard_direction=prob.sys.direction, event_before_flow=false,
    #Pathology
    last_jump_time, last_intervals, zeno_count,
    instant_jump_count, zeno_ratio,
    max_zeno_jumps, max_instant_jumps, max_buffer_size, min_zeno_history,
    zeno_floor_mult, zeno_time_threshold, zeno_reset_mult,
    beating_tol_mult, adaptive_tol_mult, adaptive_dt_mult,
    dt_initial, dt_min, event_max_iters, force_iters, guard_derivative, event_refinement_method) where {S<:Union{LinearSystem, AffineSystem}, I, T}

    xₖ = sol.x[end]
    tₖ = sol.t[end]

    sys = prob.sys

    #Check for event before taking a cont step. Just inverted the usual solve method.
    if event_before_flow
        #State is on guard surface
        if abs(guard(sys, xₖ)) <= tol

            #apply reset 
            x⁺ = apply_reset(sys, xₖ)
            jump_interval = tₖ - last_jump_time

            zeno_count, instant_jump_count, status = check_system_pathology(
                jump_interval, last_intervals, zeno_count,
                instant_jump_count, tₖ, tol, zeno_ratio,
                max_zeno_jumps, max_instant_jumps, max_buffer_size;
                min_zeno_history=min_zeno_history, zeno_floor_mult=zeno_floor_mult,
                zeno_time_threshold=zeno_time_threshold, zeno_reset_mult=zeno_reset_mult,
                beating_tol_mult=beating_tol_mult)

            if status == :terminate
                return xₖ, 0.0, Δt, true, last_jump_time, zeno_count, instant_jump_count
            end

            last_jump_time = tₖ

            push!(sol.t, tₖ, tₖ)
            push!(sol.x, xₖ, x⁺)
            push!(sol.event_times, tₖ)
            push!(sol.event_indices, length(sol.t))

            if dense_out
                push!(sol.dx, f(xₖ, tₖ))
                push!(sol.dx, f(x⁺, tₖ))
            end

            Δt_next = min(dt_initial, jump_interval * .5)

            return x⁺, 0.0, Δt_next, false, last_jump_time, zeno_count, instant_jump_count
        end
    end

    _, t_end = prob.tspan
    dt_step = (tₖ + Δt > t_end) ? (t_end - tₖ) : Δt

    #=
    if abs(guard(sys, xₖ)) < tol * adaptive_tol_mult
        dt_step = min(dt_step, dt_min * adaptive_dt_mult)
    end
    =#
    
    x_predict, eventtrigger, _, dt_used, dt_next = take_step(solver, prob, f, Df, xₖ, tₖ, dt_step, tol, sol; guard_direction=guard_direction, event_method=event_method)

    if eventtrigger
        t_star, x_star = locate_event(event_method, prob, solver, f, Df, xₖ, tₖ, dt_used, guard(sys, xₖ), tol, sol, stepper; event_max_iters=event_max_iters, force_iters=force_iters, guard_derivative=guard_derivative, event_refinement_method=event_refinement_method)

        jump_interval = t_star - last_jump_time

        zeno_count, instant_jump_count, status = check_system_pathology(
                jump_interval, last_intervals, zeno_count,
                instant_jump_count, t_star, tol, zeno_ratio,
                max_zeno_jumps, max_instant_jumps, max_buffer_size;
                min_zeno_history=min_zeno_history, zeno_floor_mult=zeno_floor_mult,
                zeno_time_threshold=zeno_time_threshold, zeno_reset_mult=zeno_reset_mult,
                beating_tol_mult=beating_tol_mult)

        if status == :terminate
            return xₖ, dt_used, dt_next, true, last_jump_time, zeno_count, instant_jump_count
        end

        last_jump_time = t_star

        x⁺ = apply_reset(sys, x_star)

        push!(sol.t, t_star, t_star)
        push!(sol.x, x_star, x⁺)

        push!(sol.event_times, t_star)
        push!(sol.event_indices, length(sol.t))

        if dense_out
            push!(sol.dx, f(x_star, t_star))
            push!(sol.dx, f(x⁺, t_star))
        end

        Δt = min(dt_initial, jump_interval * 0.5)

        return x⁺, dt_used, Δt, false, last_jump_time, zeno_count, instant_jump_count
    else
        t_next = tₖ + dt_used

        if dense_out
            push!(sol.dx, f(xₖ, tₖ))
        end

        push!(sol.t, t_next)
        push!(sol.x, x_predict)

        return x_predict, dt_used, dt_next, false, last_jump_time, zeno_count, instant_jump_count
    end
end

"""
    solve(prob::prob{S, I, T}, solver::AbstractODESolver=RK45();
          event_method::AbstractEventLocator=LinearHermite(),
          event_before_flow=false,
          dense_out=true,
          dt_initial=0.01,
          dt_min=1e-6,
          max_iter=10^6,
          tol=1e-6,
          zeno_ratio=0.90,
          max_zeno_jumps=1,
          stepper::AbstractODESolver=ModifiedTrap(),
          max_buffer_size=5,
          max_instant_jumps=5,
          guard_direction=prob.sys.direction,
          Df=(x,t) -> prob.sys.A,
          event_max_iters=100,
          force_iters=false,
          guard_derivative=nothing,
          event_refinement_method::AbstractEventRefinement=FalsePosition(),
          min_zeno_history=2,
          zeno_floor_mult=2.0,
          zeno_time_threshold=1e-2,
          zeno_reset_mult=100.0,
          beating_tol_mult=1.0,
          adaptive_tol_mult=100.0,
          adaptive_dt_mult=10.0)

The solver integrates the continuous dynamics until a guard crossing
is detected. The event location is then refined, and the appropriate
reset map is applied to obtain the post-event state.

Additional pathology detection monitors repeated and increasingly
frequent events to identify potential exeptional pathologies.

# Arguments
- `prob`: Problem containing the `LinearSystem` or `AffineSystem`, initial condition, and time span.
- `solver`: ODE solver used for numerical integration. Defaults to `RK45()`.
- `event_method`: Event locator used to detect and locate guard crossings. Defaults to `LinearHermite()`.
- `event_before_flow`: If `true`, checks for and applies a reset before taking a continuous integration step when the current state is already on the guard. If `false`, performs continuous integration before detecting events. Defaults to `false`.
- `dense_out`: If `true`, stores the derivative at each solution point in `sol.dx`. Defaults to `true`.
- `dt_initial`: Initial time step size. Defaults to `0.01`.
- `dt_min`: Minimum allowable remaining time before the end of the simulation. If the remaining time falls below this threshold, the final time is snapped to `prob.tspan[2]`. Defaults to `1e-6`.
- `max_iter`: Maximum number of integration iterations. Defaults to `10^6`.
- `tol`: Numerical tolerance used for integration, event detection, and pathology classification. Defaults to `1e-6`.
- `zeno_ratio`: Contraction ratio used to identify Zeno behavior from successive event intervals. Defaults to `0.90`.
- `max_zeno_jumps`: Maximum number of detected Zeno contractions allowed before terminating the simulation. Defaults to `1`.
- `stepper`: Secondary ODE solver used internally by the event locator to refine event locations. Defaults to `ModifiedTrap()`.
- `max_buffer_size`: Maximum number of previous event intervals retained for evaluating Zeno behavior. Defaults to `5`.
- `max_instant_jumps`: Maximum number of instantaneous jumps allowed before terminating the simulation. Defaults to `5`.
- `guard_direction`: Direction of guard crossing to detect: `0` detects either direction, `1` detects rising crossings, and `-1` detects falling crossings. Defaults to `prob.sys.direction`.
- `Df`: Optional Jacobian or derivative information supplied to the numerical solver. Defaults to the constant matrix `A` of the system, provided as `(x,t) -> prob.sys.A`.
- `event_max_iters`: Maximum number of iterations used by the event locator. Defaults to `100`.
- `force_iters`: If `true`, forces the event locator to perform the specified number of iterations rather than stopping early based on its convergence criteria. Defaults to `false`.
- `guard_derivative`: Optional function providing the derivative of the guard function, used by event refinement methods that require it. Defaults to `nothing`.
- `event_refinement_method`: Root-finding method used to refine the event location. Defaults to `FalsePosition()`.
- `min_zeno_history`: Minimum number of recorded event intervals required before attempting to calculate a Zeno contraction ratio. Defaults to `2`.
- `zeno_floor_mult`: Multiplier defining the numerical floor `tol * zeno_floor_mult` for event intervals. Defaults to `2.0`.
- `zeno_time_threshold`: Maximum event interval duration eligible for Zeno contraction classification. Defaults to `1e-2`.
- `zeno_reset_mult`: Multiplier defining the event interval `tol * zeno_reset_mult` above which the Zeno counter is reset. Defaults to `100.0`.
- `beating_tol_mult`: Multiplier defining the time window `tol * beating_tol_mult` used to classify events as instantaneous jumps or beating. Defaults to `1.0`.
- `adaptive_tol_mult`: Reserved multiplier for adaptive near-guard tolerance handling. Defaults to `100.0`.
- `adaptive_dt_mult`: Reserved multiplier for adaptive time-step handling near the guard. Defaults to `10.0`.

# Returns
A `LinearAffineSol` containing the solution times, states,
derivatives (if `dense_out=true`), and detected events.
"""
function solve(prob::prob{S, I, T},
               solver::AbstractODESolver=RK45();
               event_method::AbstractEventLocator=LinearHermite(),
               event_before_flow = false,
               dense_out = true, 
               dt_initial=0.01, dt_min = 1e-6, max_iter = 10^6, 
               tol = 1e-6, 
               zeno_ratio = 0.90, max_zeno_jumps = 1,
               stepper::AbstractODESolver=ModifiedTrap(),
               max_buffer_size=5,
               max_instant_jumps = 5,
               guard_direction = prob.sys.direction,
               Df = (x, t) -> prob.sys.A,
               event_max_iters = 100,
               force_iters = false,
               guard_derivative=nothing,
               event_refinement_method::AbstractEventRefinement=FalsePosition(),
               #Tunable pathology tolerance parameters
               min_zeno_history = 2,
               zeno_floor_mult = 2.0,
               zeno_time_threshold = 1e-2,
               zeno_reset_mult = 100.0,
               beating_tol_mult = 1.0,
               adaptive_tol_mult = 100.0,
               adaptive_dt_mult = 10.0) where {S<:Union{LinearSystem, AffineSystem}, I, T<:Tuple{Float64, Float64}}
    
    sys = prob.sys

    # Initialize solution based on linear/affine
    sol = LinearAffineSol(prob)

    function f(x, t)
        dx = sys.A * x

        if hasproperty(sys, :b)
            if x isa AbstractMatrix
                dx[:, 1] .+= sys.b
            else
                dx .+= sys.b
            end
        end
        return dx
    end

    # Extract start and end times
    t_start, t_end = prob.tspan

    # Initialize time step and iter counter
    Δt = dt_initial
    iter = 0

    # Trackers for beating, blocking, and Zeno logic
    instant_jump_count = 0
    zeno_count = 0
    last_jump_time = t_start      
    last_intervals = Float64[]     

    #Check if we start on the guard. 
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
        # Initial derivative if not on guard
        if dense_out
            push!(sol.dx, f(prob.init, prob.tspan[1]))
        end
    end

    # Run until end time or max iter
    while sol.t[end] < t_end
        iter += 1
        if iter > max_iter 
            @warn "Maximum Iteration Count ($max_iter) reached."
            break
        end

        # Terminate if time is below machine precision
        if t_end - sol.t[end] < dt_min
            sol.t[end] = t_end
            @info "Time to end of simulation below minimum time step. Snapping final time to t = $t_end."
            break
        end

        # Truncate time step if we overshoot the final time
        _, _, Δt, terminate, last_jump_time, zeno_count, instant_jump_count = take_step_linear_affine!(
            solver, prob, f, Df, Δt, tol, sol; dense_out=dense_out,
            stepper=stepper, event_method=event_method, guard_direction=guard_direction,
            event_before_flow,
            last_jump_time=last_jump_time, last_intervals=last_intervals, zeno_count=zeno_count,
            instant_jump_count=instant_jump_count, zeno_ratio=zeno_ratio, max_zeno_jumps=max_zeno_jumps,
            max_instant_jumps=max_instant_jumps, max_buffer_size=max_buffer_size,
            min_zeno_history=min_zeno_history, zeno_floor_mult=zeno_floor_mult, zeno_time_threshold=zeno_time_threshold,
            zeno_reset_mult=zeno_reset_mult, beating_tol_mult=beating_tol_mult,
            adaptive_tol_mult=adaptive_tol_mult, adaptive_dt_mult=adaptive_dt_mult,
            dt_initial=dt_initial, dt_min=dt_min, event_max_iters=event_max_iters, force_iters=force_iters, guard_derivative=guard_derivative, event_refinement_method=event_refinement_method)

        if terminate
            break
        end
    end
    return sol
end

