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
function AffineSystem(A::AbstractMatrix, b::AbstractVector, λ::AbstractVector, a::Real, C::AbstractMatrix, κ::AbstractVector; direction::Int=0)
    return AffineSystem(Float64.(A), Float64.(b), Float64.(λ), Float64(a), Float64.(C), Float64.(κ), direction)
end

#SOLUTION STRUCTS AND HELPERS. 
#Goal to provide standard date for simulation outputs. Keeps cont trajectories and discrete events organized. Currently affine and linear are the samem but kept separate to allow specific plotting later?
struct LinearAffineSol{T, DX} <: AbstractHybridSolution
    t::Vector{Float64}
    x::Vector{T}
    dx::DX 
    event_times::Vector{Float64}
    event_indices::Vector{Int}
end 

# CK: Is this ever used?
function LinearAffineSol(prob::prob{S, I, T}, t::AbstractVector, x::AbstractVector,
    event_times::AbstractVector, event_indices::AbstractVector) where {S <: Union{LinearSystem, AffineSystem}, I, T}
    return LinearAffineSol(Float64.(t), Vector{I}(x), Vector{Vector{Float64}}(), Float64.(event_times), Vector{Int}(event_indices))
end

#Exact Linear Flow (matrix exponential)
#The goal is to instead of using numerical approximations for linear systems we can get exact solutions. 
struct LinearFlow{TM<:AbstractMatrix{Float64}}
    A::TM #dx/dt = Ax
end

#Internal
#Constructor with numerical/structural checks
function LinearFlow(A::AbstractMatrix)
    #Make sure matrix is square
    m, n = size(A)
    if m != n
        throw(DimensionMismatch("State matrix A must be square, but got ($m × $n) matrix."))
    end

    return LinearFlow(Float64.(A))
end

function flow(flowmap::LinearFlow, τ::Real, x::AbstractVector)
    #check matrix and state vector match dim 
    n = size(flowmap.A, 1)
    if length(x) != n
        throw(DimensionMismatch("State vector length ($(length(x))) does not match system matrix dimension ($n)."))
    end

    #Computes exact state x(t+τ) = exp(A*τ) * x(t)
    #This is very memory instensive as of now. Fixing will come later
    return exp(flowmap.A .* τ) * x
end

#Solution Initialization
#Goal is to setup the empty memory containers before solver starts running. 
#pre-allocating the vectors with the exact starting conditions ensures that the solver loop is stable 

#internal
function LinearAffineSol(prob::prob{S, I, T}) where {S<:Union{LinearSystem, AffineSystem}, I, T}
    return LinearAffineSol([prob.tspan[1]], [prob.init], Vector{Vector{Float64}}(), Float64[], Int[])
end

"""
    beating_and_blocking_sets(sys::Union{LinearSystem, AffineSystem})

Compute the discrete sequence of constraint matrices and offset vectors defining the beating sets and the final blocking set for both linear and affine systems

Returns:

A unified "NamedTuple" *fancy*
- 'beating_sets::Vector{Matrix{Float64}}': Matrix constraints where rows represent normal vectors. Constraint matrix is the [λ^T ; λ^TC^{-1} ; ...]
- 'beating_offsets::Vector{Vector{Float64}}': Vector offsets shifting the hyperplanes from the origin
- 'blocking_set::Matrix{Float64}': Final blocking constraint matrix
- 'blocking_offsets::Vector{Float64}': Final offset for blocking set
- 'k_blocking::Int': Discrete iteration where the sequence achieved blocking.


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

    Determine if the system's blocking set collapses strictly to the origin (Σ_∞ = {0}).
    For both linear and affine systems, this requires the final constraint matrix to 
    have full rank and the final offset vector to be (basically) zero.

"""
function is_trivially_blocking(sys::Union{LinearSystem, AffineSystem})
    analysis = beating_and_blocking_sets(sys)
    n = length(sys.λ)

    #check full rank AND the point is the oriign
    return rank(analysis.blocking_set) == n && isapprox(norm(analysis.blocking_offsets))
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

function take_step_linear_affine!(solver, prob::prob{S, I, T}, f, Δt, tol, sol; dense_out=true, stepper::AbstractODESolver=RK45(), event_method::AbstractEventLocator=LinearLocator(), guard_direction=prob.sys.direction,
    #Pathology
    last_jump_time, last_intervals, zeno_count,
    instant_jump_count, zeno_ratio,
    max_zeno_jumps, max_instant_jumps, max_buffer_size, min_zeno_history,
    zeno_floor_mult, zeno_time_threshold, zeno_reset_mult,
    beating_tol_mult, adaptive_tol_mult, adaptive_dt_mult,
    dt_initial, dt_min) where {S<:Union{LinearSystem, AffineSystem}, I, T}

    xₖ = sol.x[end]
    tₖ = sol.t[end]

    sys = prob.sys

    _, t_end = prob.tspan
    dt_step = (tₖ + Δt > t_end) ? (t_end - tₖ) : Δt

    if abs(guard(sys, xₖ)) < tol * adaptive_tol_mult
        dt_step = min(dt_step, dt_min * adaptive_dt_mult)
    end

    x_predict, eventtrigger, t_root, dt_used, dt_next = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol; guard_direction=guard_direction)

    if eventtrigger
        t_root = guard(sys, xₖ)

        t_star, x_star = locate_event(event_method, prob, solver, f, xₖ, tₖ, dt_used, t_root, tol, sol, stepper)

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

        push!(sol.t, t_star)
        push!(sol.x, x_star)

        push!(sol.event_times, t_star)
        push!(sol.event_indices, length(sol.t))

        push!(sol.t, t_star)
        push!(sol.x, x⁺)

        if dense_out
            push!(sol.dx, f(x_star, t_star))
            push!(sol.dx, f(x⁺, t_star))
        end

        Δt = min(dt_initial, jump_interval * 0.5)

        return x⁺, dt_used, Δt, false, last_jump_time, zeno_count, instant_jump_count
    else
        t_next = tₖ + dt_used

        if dense_out
            push!(sol.dx, f(x_predict, t_next))
        end

        push!(sol.t, t_next)
        push!(sol.x, x_predict)

        return x_predict, dt_used, dt_next, false, last_jump_time, zeno_count, instant_jump_count
    end
end

"""
    solve(prob::prob{F, I, T}, solver::AbstractODESolver=RK45(); kwargs...) where {F<:Union{LinearSystem, AffineSystem}, I, T<:Tuple{Float64, Float64}}

# ARGUMENT KEY

## Required

* 'prob': The problem definition containing the system dynamics 'sys', initial state, and time span.
* 'solver': The numerical integration method used for continous steps. Defaults to RK45().

## Optional

### Simulation and Step Controls:
* 'dt_initial' (Float64, default '0.01'): The starting time step for the continuous solver.
* 'dt_min' (Float64, default '1e-6'): The absolute minimum allowable time step. If the solver tries to go below this, the simulation terminates.
* 'max_iter' (Int, default '10^6'): The maximum number of continuous integration steps allowed before forcing a timeout.
'tol' (Float64, default '1e-6'): The baseline numerical tolerance used across the solver. Acts as the foundational unit for multipliers below.

### Event Handling
* 'event_method' (AbstractEventLocator, default 'LinearLocator()'): The algorithm used to pinpoint the exact time and state of a guard crossing.
* 'stepper' (AbstractODESolver. default 'RK4()'): The secondary ODE solver used internally by the event detection locator to pinpoint the impact state.
* 'guard_direction' (Symbol, default 'default_guard_direction(prob.sys)'): Dictates which zero-crossing direction triggers an event (e.g., positive-to-negative).

### Pathology Tuning
* 'zeno_ratio' (Float64, default '.90'): The ratio threshold for Zeno detection. If consecutive jump intervals contract by this ratio (or faster) it triggers a Zeno classification. 
* 'max_zeno_jumps' (Int, default '5'): The maximum number of consecutuve Zeno contractions allowed before terminating the simulation. 
* 'max_instant_jumps' (Int, default '5'): The maximum number of instantaneous jumps allowed before classifying the system as blocked and terminating. 
* 'max_buffer_size' (Int, default '5'): The number of previous jump intervals stored in memory to evaluate Zeno contractions.

### Fine-Tuning Multipliers (scaled against 'tol')
* 'min_zeno_history' (Int, default '2'): The minimum number of recorded jumps required before the solver will attempt to calculate a Zeno contraction ratio. 
* 'zeno_floor_mult' (Float64, default '2.0'): Defines the numerical floor ('tol * zeno_floor_mult'). If a jump interval falls below this, it maintains a Zeno state to prevent
machine precision drops into beating blocks.
* 'zeno_time_threshold' (Float64, default '1e-2'): The absolute maximum duration a jump interval can be while still being eligible for Zeno contraction classification.
* 'zeno_reset_mult' (Float64, default '100.0'): If a jump interval exceeds 'tol * zeno_reset_mult', the system is deemed safe and the Zeno counter is decremented.
* 'beating_tol_mult' (Float64, default '1.0'): Defines the time window ('tol * beating_tol_mult'). Jumps occurring within this window are classifed as instantaneous jumps or 'beating'
* 'adaptive_tol_mult' (Float64, defualt '100.0'): The boundary distance multiplier ('tol * adaptive_tol_mult'). When the system enters this distance from the guard, it shrinks step size to increase resolution.
* 'sliding_tol_mult' (Float64, default '10.0'): If the post-impact guard value is within 'tol * sliding_tol_mult', the solver enters sliding mode to suppress immediate erroneous events. This helps us avoid strange chattering.  

"""
function solve(prob::prob{S, I, T}, 
               solver::AbstractODESolver=RK45(); 
               event_method::AbstractEventLocator=LinearLocator(),
               dense_out = true, 
               dt_initial=0.01, dt_min = 1e-6, max_iter = 10^6, 
               tol = 1e-6, 
               zeno_ratio = 0.90, max_zeno_jumps = 5,
               stepper::AbstractODESolver=ModifiedTrap(),
               max_buffer_size=5,
               max_instant_jumps = 5,
               guard_direction = prob.sys.direction,
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

    # Run until end time or max iter
    while sol.t[end] < t_end
        iter += 1
        if iter > max_iter 
            @warn "Maximum Iteration Count ($max_iter) reached."
            break
        end

        # Terminate if time is below machine precision
        if t_end - sol.t[end] < dt_min
            @info "Time to end of simulation below minimum time step. Ending simulation at t = $(sol.t[end])"
            break
        end

        # Truncate time step if we overshoot the final time
        _, _, Δt, terminate, last_jump_time, zeno_count, instant_jump_count = take_step_linear_affine!(
            solver, prob, f, Δt, tol, sol; dense_out=dense_out,
            stepper=stepper, event_method=event_method, guard_direction=guard_direction,
            last_jump_time=last_jump_time, last_intervals=last_intervals, zeno_count=zeno_count,
            instant_jump_count=instant_jump_count, zeno_ratio=zeno_ratio, max_zeno_jumps=max_zeno_jumps,
            max_instant_jumps=max_instant_jumps, max_buffer_size=max_buffer_size,
            min_zeno_history=min_zeno_history, zeno_floor_mult=zeno_floor_mult, zeno_time_threshold=zeno_time_threshold,
            zeno_reset_mult=zeno_reset_mult, beating_tol_mult=beating_tol_mult,
            adaptive_tol_mult=adaptive_tol_mult, adaptive_dt_mult=adaptive_dt_mult,
            dt_initial=dt_initial, dt_min=dt_min)

        if terminate
            break
        end
    end
    return sol
end

