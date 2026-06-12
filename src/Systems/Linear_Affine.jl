#EVERYTHING FOR LINEAR SYSTEMS 

#Linear System Structs: Goal to define the date for a linear system and its problem.
#We subtype it to AbstractHybridSystem/Problem for campatibility with solvers 

struct LinearSystem <: AbstractHybridSystem
    A::Matrix{Float64} #State Transition matrix (dx/dt = Ax)
    λ::Vector{Float64} #Normal Vector for the Guard Surface
    C::Matrix{Float64} #Reset map matrix (x⁺ = Cx)
end

#External
#Constructor to help user see data types 
function LinearSystem(A::AbstractMatrix, λ::AbstractVector, C::AbstractMatrix)
    return LinearSystem(Float64.(A), Float64.(λ), Float64.(C))
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
end

#External
# Constructor to help user see data types
function AffineSystem(A::AbstractMatrix, b::AbstractVector, λ::AbstractVector, a::Real, C::AbstractMatrix, κ::AbstractVector)
    return AffineSystem(Float64.(A), Float64.(b), Float64.(λ), Float64(a), Float64.(C), Float64.(κ))
end

#SOLUTION STRUCTS AND HELPERS. 
#Goal to provide standard date for simulation outputs. Keeps cont trajectories and discrete events organized. Currently affine and linear are the samem but kept separate to allow specific plotting later?
struct HybridSolution <: AbstractHybridSolution
    t::Vector{Float64}
    x::Vector{Vector{Float64}}
    jump_times::Vector{Float64}
    jump_indices::Vector{Int}
end 

function CreateSolution(prob::prob{S, I, T}, t::AbstractVector, x::AbstractVector,
    jump_times::AbstractVector, jump_indices::AbstractVector) where {S <: Union{LinearSystem, AffineSystem}, I, T}
    return HybridSolution(Float64.(t), Vector{Float64}.x, Float64.(jump_times), Int{jump_indices})
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

#Beating/Zeno check
#Goal is to prevent loops in Linear and Affine systems.
#Through Multiple dispath, this overrides the default fallback in Definitions.jl only when the solver is Linear or Affine. 

#Internal (could be external if we want?)
function check_beating_status(sys::Union{LinearSystem, AffineSystem}, instant_jumps, n, x_current, t_current, tol, trivial_tol_multiplier)
    if (instant_jumps) > (n-1) 
        if norm(x_current) < tol * trivial_tol_multiplier
            @info "Trivial Blocking: System settled at origin at t = $t_current"
            return :blocking_trivial 
        end
        x_next = apply_reset(sys, x_current)

        if norm(x_next) < norm(x_current) * (1-tol)
            @info "Contractive Beating: State shrinking toward origin at t = $t_current"
            return :contractive_beating
        elseif norm(x_next) > norm(x_current) * (1+tol)
            @warn "Expansive Blocking: State trapped and expanding on guard at t = $t_current"
            return :blocking_expansive
        else 
            @warn "Expansive Blocking: State trapped on guard at t = $t_current"
            return :blocking_non_trivial 
        end
    end
    return :continue
end

#Solution Initialization
#Goal is to setup the empty memory containers before solver starts running. 
#pre-allocating the vectors with the exact starting conditions ensures that the solver loop is stable 

#internal
function init_solution(prob::prob{S, I, T}) where {S<:Union{LinearSystem, AffineSystem}, I, T}
    return HybridSolution([prob.tspan[1]], [prob.init], Float64[], Int[])
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
function guard(sys::LinearSystem, x::AbstractVector)
    return sys.λ' * x
end
#Internal
function guard(sys::AffineSystem, x::AbstractVector)
    return sys.λ' * x + sys.a 
end

#internal
function apply_reset(sys::LinearSystem, x::AbstractVector)
    return sys.C * x
end
#internal
function apply_reset(sys::AffineSystem, x::AbstractVector)
    return sys.C * x + sys.κ
end

#Internal
function get_dimension(sys::LinearSystem)
    n = size(sys.A, 1)
    if size(sys.A, 2) != n || length(sys.λ) != n || size(sys.C, 1) != n || size(sys.C, 2) != n
        throw(DimensionMismatch("LinearSystem dimensions are inconsistent."))
    end
    return n
end
function get_dimension(sys::AffineSystem)
    n = size(sys.A, 1)
    # Check consistency of vectors and matrices
    if size(sys.A, 2) != n || length(sys.b) != n || length(sys.λ) != n || 
       size(sys.C, 1) != n || size(sys.C, 2) != n || length(sys.κ) != n
        throw(DimensionMismatch("AffineSystem dimensions are inconsistent."))
    end
    return n
end

#SOLVER
#Very External
function solve(prob::prob{F, I, T}, 
               solver::AbstractODESolver=RK45(); 
               event_method::AbstractEventLocator=LinearLocator(), 
               dt_initial=.01, dt_min = 1e-6, max_iter = 10^6, 
               tol = 1e-6, 
               zeno_ratio = 0.90, max_zeno_jumps = 5,
               stepper::AbstractODESolver=ModifiedTrap(),
               max_buffer_size=5,
               beating_warn_threshold=3,
               max_instant_jumps = 5) where {F<:Union{LinearSystem, AffineSystem}, I<:AbstractVector{Float64}, T<:Tuple{Float64, Float64}}
    sys = prob.sys

    f = hasproperty(sys, :b) ? ((x,t) -> sys.A * x + sys.b) : ((x,t) -> sys.A * x) 

    #Initialize soltution based on linear/affine
    sol = init_solution(prob)

    #extract start and end times
    t_start, t_end = prob.tspan

    #initialize time step and iter counter
    Δt = dt_initial
    iter = 0

    #trackers for beating blocking and Zeno logic
    instant_jump_count = 0
    zeno_count = 0
    last_jump_time = t_start      #calc the current interval
    last_intervals = Float64[]  #History of intervals     
    n = length(prob.init)

    #run until end time or max iter
    while sol.t[end] < t_end
        iter += 1
        if iter > max_iter 
            @warn "Maximum Iteration Count ($max_iter) reached."
            break
        end

        #terminate if time is below machine precision
        if t_end - sol.t[end] < dt_min
            @info "Time to end of simulation below minimum time step. Ending simulation at t = $(sol.t[end])"
            break
        end

        #truncate time step if we overshoot the final time
        dt_step = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

        #continuous integration
        xₖ = sol.x[end]
        tₖ = sol.t[end]

        h_now = guard(sys, xₖ)

        #attempt continuous step
        x_predict, eventtriggered, _, dt_used, dt_next = take_step(solver, prob, f, xₖ, tₖ, dt_step, tol, sol)

        is_exactly_on_guard = abs(h_now) <= tol

        #discrete event logic
        if eventtriggered || is_exactly_on_guard

            if is_exactly_on_guard
                t_star, x_star = tₖ, xₖ
            else
                t_star, x_star = locate_event(event_method, prob, solver, f, xₖ, tₖ, dt_used, h_now, tol, sol, stepper)
            end
            
            #PATHOLOGY CHECK
            jump_interval = t_star - last_jump_time
            zeno_count, instant_jump_count, status = check_system_pathology(
                jump_interval, last_intervals, 
                zeno_count, instant_jump_count,
                t_star, tol, zeno_ratio, max_zeno_jumps, max_instant_jumps,
                max_buffer_size
            )

            if status == :terminate
                break
            end

            
            last_jump_time = t_star

            #Apply Reset
            x⁺ = apply_reset(sys, x_star)


            #Tracking data
            push!(sol.t, t_star, t_star)
            push!(sol.x, x_star, x⁺)

            #record event data for analysis
            if hasproperty(sol, :jump_times)
                push!(sol.jump_times, t_star)
                push!(sol.jump_indices, length(sol.t))
            end

            #shrink min step size to avoid overshooting
            Δt = dt_initial

        else
            t_next = tₖ + dt_used
            push!(sol.t, t_next)
            push!(sol.x, x_predict)

            #reset step size
            Δt = dt_next

        end

        end
        return sol
    end 

