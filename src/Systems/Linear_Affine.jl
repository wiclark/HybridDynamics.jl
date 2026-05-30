"""
LINEAR AND AFFINE SYSTEMS TOOLKIT

TABLE OF CONTENTS
1) System structs and Problem Definitions
2) Interface Fullfillment (get_dimension, vector_field, guard, apply_reset)
3) Exact Flow (Matrix Exponential)
4) Simulation Utilities (Bisection and Event Check)
"""
#EVERYTHING FOR LINEAR SYSTEMS 

#Linear System Structs: Goal to define the date for a linear system and its problem.
#We subtype it to AbstractHybridSystem/Problem for campatibility with solvers 

struct LinearSystem <: AbstractHybridSystem
    A::Matrix{Float64} #State Transition matrix (dx/dt = Ax)
    λ::Vector{Float64} #Normal Vector for the Guard Surface
    C::Matrix{Float64} #Reset map matrix (x⁺ = Cx)
end
struct LinearProblem <: AbstractHybridProblem
    sys::LinearSystem               #Phyisical System defined above
    x₀::Vector{Float64}             #Initial state vector
    tspan::Tuple{Float64, Float64}  #(start time, end time)
end

#Linear Interface Functions (see Definitions file)
#Goal to satisfy the 4 methods defined in the Def file. The solver never looks inside sys.A or sys.λ to hopefully make faster
get_dimension(sys::LinearSystem) = size(sys.A, 1)       #return n from the n × n state matrix
vector_field(sys::LinearSystem) = (x, t) -> sys.A * x   #Continuous dynamics: dx/dt = Ax
guard(sys::LinearSystem, x) = dot(sys.λ, x)             #Event Surface: Evals to 0 on guard
apply_reset(sys::LinearSystem, x) = sys.C * x           #Discrete Dynamics: applies the jump matrix


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

struct AffineProblem <: AbstractHybridProblem
    sys::AffineSystem               #Physical system defined above
    x₀::Vector{Float64}             #initial condition
    tspan::Tuple{Float64, Float64}  #time span
end

# Implementations of the interface functions. Shows we can add more to the solve loop easily
get_dimension(sys::AffineSystem) = size(sys.A, 1)                   #return n from the n × n state matrix
vector_field(sys::AffineSystem) = (x, t) -> sys.A * x + sys.b       #Continuous dynamics: dx/dt = Ax + b
guard(sys::AffineSystem, x) = dot(sys.λ, x) + sys.a                 #Event Surface with the 'a' offset
apply_reset(sys::AffineSystem, x) = sys.C * x + sys.κ               #Discrete Dynamics: applies the jump matrix with 'κ' offset

# Constructor to help user see data types
function AffineSystem(A::AbstractMatrix, b::AbstractVector, λ::AbstractVector, a::Real, C::AbstractMatrix, κ::AbstractVector)
    return AffineSystem(Float64.(A), Float64.(b), Float64.(λ), Float64(a), Float64.(C), Float64.(κ))
end

#SOLUTION STRUCTS AND HELPERS. 
#Goal to provide standard date for simulation outputs. Keeps cont trajectories and discrete events organized. Currently affine and linear are the samem but kept separate to allow specific plotting later?
struct LinearSolution
    t::Vector{Float64}          #Time history
    x::Vector{Vector{Float64}}  #State History
    jump_times::Vector{Float64} #Exact timestamps where an event has occurred 
    jump_indices::Vector{Int}   #Indices in 'x' and 't' where jumps map to
end
struct AffineSolution
    t::Vector{Float64}          #Time history
    x::Vector{Vector{Float64}}  #State History
    jump_times::Vector{Float64} #Exact timestamps where an event has occurred
    jump_indices::Vector{Int}   #Indices in 'x' and 't' where jumps map to
end

#Constructor
function CreateSolution(prob::LinearProblem, t::AbstractVector, x::AbstractVector, jump_times::AbstractVector, jump_indices::AbstractVector)
    return LinearSolution(Float64.(t), Vector{Float64}.(x), Float64.(jump_times), Int.(jump_indices))
end
#Constructor
function CreateSolution(prob::AffineProblem, t::AbstractVector, x::AbstractVector, jump_times::AbstractVector, jump_indices::AbstractVector)
    return AffineSolution(Float64.(t), Vector{Float64}.(x), Float64.(jump_times), Int.(jump_indices))
end

#Exact Linear Flow (matrix exponential)
#The goal is to instead of using numerical approximations for linear systems we can get exact solutions. 
struct LinearFlow{TM<:AbstractMatrix, TV<:AbstractVector}
    V::TM       #Eigenvectors
    Λ::TV       #Eigenvalues
    V_inv::TM   #inverse of eigenvectors
end

#Constructor that performs heavy math ONCE
function LinearFlow(A)
    eig = eigen(A) #Perform eigen decomp

    #We factor A = V * Λ * V{-1} so we can compute cont evolution of the system.
    return LinearFlow(eig.vectors, eig.values, inv(eig.vectors))
end

#Actual step function for the exact solver
function flow(flowmap::LinearFlow, τ, x) 

    #projects the state into the eigendecomposition space
    y = flowmap.V_inv * x 

    #Computes the exact state at time τ.
    #Note we use real. because complex conj pairs will result in tiny numerical evils.
    return real.(flowmap.V * (exp.(flowmap.Λ .* τ) .* y)) 
end

#Beating/Zeno check
#Goal is to prevent loops in Linear and Affine systems.
#Through Multiple dispath, this overrides the default fallback in Definitions.jl only when the solver is Linear or Affine. 
function check_beating_status(sys::Union{LinearSystem, AffineSystem}, instant_jumps, n, x_current, t_current, tol)
    if (instant_jumps) > (n-1) 
        if norm(x_current) < tol * 10
            @info "Trivial Blocking: System settled at origin at t = $t_current"
            return :blocking_trivial 
        end
        x_next = apply_reset(sys, x_current)

        if norm(x_next) < norm(x_current) * (1-tol)
            @info "Contractive Beating: State shrinking toward origin at t = $t_current"
            return :continue #Zeno logic here eventually I think
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
function init_solution(prob::LinearProblem) 
    return LinearSolution([prob.tspan[1]], [prob.x₀], Float64[], Int[])
end
function init_solution(prob::AffineProblem)
    return AffineSolution([prob.tspan[1]], [prob.x₀], Float64[], Int[])
end

#--------------------------------------------
#BEATING AND BLOCKING SETS
function beating_and_blocking_sets(sys::AbstractHybridSystem)
    #Deconstructs the system structure to get the matrices and vectors we need.
    λ, C = sys.λ, sys.C
    
    #Sote state space dimension; we want it to stabilize within n steps. 
    n = length(λ)

    #compute inverse of C. We need it to compute Σ_k = Σ ∩ CΣ_k-1
    C_inv = inv(C)

    #intialize contianer to store 
    O = Matrix{Float64}[] 

    #Convert the gaurd vector λ  1 times n matrix to represent λ^T x = 0
    O_current = reshape(λ, 1, :)

    #store inital guard constraint: \Sigma_0 = {x | λ^T x = 0}
    push!(O, O_current)

    #Establish initial rank to track dimension of subspace dim = n-rank
    prev_rank = rank(O_current)

    k_∞ = 0 #tracks when blocking occurs 

    row = sys.λ'

    #We iterate n times. There is a theorem about this I think
    for k in 1:n 
        #Calc the kth constraint row. This maps the guard back through k jumps. 
        row = (row * C_inv)
        new_row = reshape(row, 1, :)

        #Concatenate new constraint to previous matrix to form next candidate
        O_next = vcat(O_current, new_row)

        #Compute rank of new matrix to see if the added constraint is linearly independent. If not, we have found the blocking set and can stop.
        current_rank = rank(O_next)

        #If rank not increased the new constraint is redundant and we are done
        if current_rank == prev_rank

            #sequence stabilized; previous index marks start of blocking set. 
            k_∞ = k - 1
            break
        end

        #Update constraint matrix with ind row
        O_current = O_next
        
        #log matrix O_k representing the kth beating set
        push!(O, O_current)
        
        #update rank baseline for next iteration
        prev_rank = current_rank

        #Tentatively set blocking index to current iteration. 
        k_∞ = k
    end

    #return a tuple for access to whole thing
    return (
        beating_sets = O,       #Full sequence of constraint matrices which are our beating sets. Each one is a matrix where the rows represent linear constraints that define the set.
        blocking_set = O[end],  #Final matrix where we have blocking
        k_blocking = k_∞        #Smallest integer k such that we have the blocking set.
    )
end
function beating_and_blocking_sets(sys::AffineSystem)
    λ, a, C, κ = sys.λ, sys.a, sys.C, sys.κ

    n = length(λ)
    C_inv = inv(C)

    O = Matrix{Float64}[] #this will store the sequence of beating sets as matrices where the rows are the linear constraints defining the set.
    b_vecs = Vector{Float64}[] #this will store the corresponding constant terms for the affine constraints.

    #Σ_0 is defined by λ^T x + a = 0, which we can write as λ^T x = -a. 
    O_current = reshape(λ, 1, :) #constraint matrix for Σ_0
    b_current = [-a] #constant term for Σ_0

    push!(O, O_current)
    push!(b_vecs, b_current)

    prev_rank = rank(O_current)
    k_∞ = 0 #tracks when blocking occurs

    for k in 1:n 
        #linear constraints for λ^T(C^{-1}^k)
        new_row = reshape(λ' * (C_inv^k), 1, :)

        #for x to be in Σ_k, then Cx+κ must be in Σ_k-1
        #This takes reset constant κ back through powers of C_inv
        new_offset = -a - sum(λ' * (C_inv^j) * κ for j in 1:k)

        #append new row of O matrix
        O_next = vcat(O_current, new_row)

        #append new offset to constant b vector 
        b_next = vcat(b_current, [new_offset])

        #Current rank of matrix to see if ind
        current_rank = rank(O_next)

        #Check for blocking
        if current_rank == prev_rank
            k_∞ = k - 1
            break
        end

        #Update current constraints with new ind affine hyperplane 
        O_current = O_next
        b_current = b_next

        #Store matrix and vector that define Σ_k = {x:O_k x = b_k}
        push!(O,O_current)
        push!(b_vecs, b_current)

        prev_rank = current_rank
        k_∞ = k
    end
    return (
        beating_sets = O,                   #Array matrices where each row is a normal vector to guard preimage. Geometrically this gives the orientation of each set. 
        beating_offsets = b_vecs,           #Array of vectors defining the offset of the constraints. This is how far the above is shifted from the origin. 
        blocking_set = O[end],              #Specific linear matrix defining blocking set
        blocking_offsets = b_vecs[end],     #Specific Offset vector of final blocking set
        k_blocking = k_∞                    #Discrete index where we have blocking. 
    )
end

#TRIVIALLY BLOCKING
function is_trivially_blocking(sys::LinearSystem)

    #The blocking set is the null space of the matrix O we formed in the analysis part. 
    #Rank nullity says Rank(blockingset) + nullity(blockingset) = n, where n is the dimension of the state space. Rank is lin ind rows, and nullity is dim of blocking set
    
    #the set is trivially blocking if Σ_∞ = {0}.
    #dim(Σ_∞) = 0 iff nullity(O[end]) = 0 \implies rank(O[end])+0= n \implies rank(O[end]) = n. So we just check if the final matrix has full rank.
    analysis = beating_and_blocking_sets(sys)

    #gets dimension of state space from length of λ, which is the same as number of columns in O[end]
    n = length(sys.λ)

    #if rank == n then nullity = 0. 
    #if nullity == 0 then the solution to O[end] x = 0 is only x = 0, so the blocking set is just the origin, which is trivially blocking.
    #we return true is the blocking set is trivial and false otherwise. 
    return rank(analysis.blocking_set) == n

end
function is_trivially_blocking(sys::AffineSystem)
    #Perform beating and blocking sets analysis to find constraint matrix and offset
    analysis = beating_and_blocking_sets(sys)

    #dim of state space
    n = length(sys.λ)
     
    return rank(analysis.blocking_set) == n && isapprox(norm(analysis.blocking_offsets), 0, atol=1e-12)
end

function basis_beating_and_blocking_sets(sys::Union{LinearSystem, AffineSystem})
    #run function we already have to get constraint matrices
    analysis = beating_and_blocking_sets(sys)

    #convert the results from above to explicit bases. The nullspace function will return a matrix where columns form a orthonormal basis
    explicit_beating = Matrix{Float64}[]
    for i in analysis.beating_sets
        #compute the null space if it is just the origin this is an empty matrix
        basis_matrix = nullspace(i)
        push!(explicit_beating, basis_matrix)
    end

    #convert final blocking set to its explicit basis 
    explicit_blocking = nullspace(analysis.blocking_set)

    #return tuple with tangible basis vectors. This will return a n x x matrix where x is the dimension of the set. 
    return (
        explicit_beating_sets = explicit_beating,
        explicit_blocking_set = explicit_blocking,
        k_blocking = analysis.k_blocking
    )
end


