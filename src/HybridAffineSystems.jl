module HybridDynamics

using ForwardDiff
using LinearAlgebra


export HybridAffineSystem, solve_hybrid_affine_system
export forward_euler_step, modified_trap_step, modified_midpoint_step, rk_23_step
export beating_and_blocking_sets_affine, is_trivially_blocking_affine, basis_beating_and_blocking_sets


struct HybridAffineSystem
    A::Matrix{Float64} #Flow matrix
    b::Vector{Float64} #Flow Constant
    λ::Vector{Float64} #Guard normal
    a::Float64         #Gaurd Constant
    C::Matrix{Float64} #Reset Matrix
    κ::Vector{Float64} #Reset Constant
end

struct HybridAffineProblem
    sys::HybridAffineSystem
    x₀::Vector{Float64}
    tspan::Tuple{Float64,Float64}
end


include("ODE_solvers.jl")


function solve_hybrid_affine_system(problem::HybridAffineProblem, dt_initial::Float64; step_method=forward_euler_step, is_adaptive=false)
    
    #Deconstruct the system structure.
    sys = problem.sys
    A, b, λ, a, C, κ = sys.A, sys.b, sys.λ, sys.a, sys.C, sys.κ

    #unpacks time span
    t_start, t_end = problem.tspan

    #Initialize time and output vectors starting with initial conditions.
    t = [t_start]
    x = [problem.x₀]

    #empty list to store timestamps where jumps occur
    jump_times = Float64[]
    jump_indices = Int[]

    #Defines the vector field for continuous part
    f(x, t) = A * x + b

    #Sets initial time step size. This can change if the solver we use is adaptive or needs to hit a specific time boundary.
    Δt = dt_initial

    #This main loop continuous as long as the most recent recorded time is less than final
    while t[end] < t_end

        #If the proposed step would overshoot and jump past where we want, we shrink the step size to land exactly at the end.
        if t[end] + Δt > t_end
            Δt = t_end - t[end]
        end

        #Grabs most recent state and time. This lets us always start from latest point, including after jumps.
        xₖ = x[end]
        tₖ = t[end]

        #Logic for solver type
        if is_adaptive
            #If adaptive, the step_method  returns new state, the actual step size is used and a suggestion for the next step size.
            x_predict, Δt_used, Δt_next = step_method(f, xₖ, Δt, tₖ, t_end)
        else 
            #Execute numerical integration step (e.g. forward Euler) to predict the state at next time.
            x_predict = step_method(f, xₖ, Δt, tₖ)
            Δt_used = Δt
            Δt_next = Δt
        end

        #calc timestamp for predicted state
        t_next = tₖ + Δt_used

        #evaluates the affine version at start of next step
        h_now = dot(λ, xₖ) + a

        #evaluates the affine version at end of next step
        h_next = dot(λ, x_predict) + a

        #Jump condition
        if h_now * h_next < 0.0 && t_next < t_end

            #calculates the time of the jump using linear interpolation.
            θ = -h_now / (h_next - h_now)

            #calulates the time where trajectory hits the guard
            t_star = tₖ + θ * Δt_used

            #interpolates the state vector at jump time, directly before jump
            x⁻ = xₖ + θ * (x_predict - xₖ) #pre jump state at t_star

            #logs the state at moment of impact 
            push!(x, x⁻)
            push!(t, t_star)

            #JUMP
            #applies C and \kappa to the impact state and logs in t_star. 
            x⁺ = C * x⁻ + κ #post jump state at t_star

            push!(x, x⁺)
            push!(t,t_star)

            #update jump times and indices
            push!(jump_times, t_star)
            push!(jump_indices, length(x))

            #resets step size. 
            Δt = dt_initial

            #check for zeno by looking at the last two jump intervals
            if length(jump_times) > 2
                last_dt = jump_times[end] - jump_times[end - 1] 
                
                #if intervals are shrinking rapidly reduce Δt to follow it
                Δt = min(dt_initial, last_dt * 0.5)
            end
        else 

            #if no jump, the pred state and time are just appended to output.
            push!(x, x_predict)
            push!(t, t_next)

            #update step size if adaptive
            Δt = Δt_next
        end
    end
    return t, x, jump_times, jump_indices
end

function beating_and_blocking_sets_affine(sys::HybridAffineSystem)
    A, λ, a, C, κ = sys.A, sys.λ, sys.a, sys.C, sys.κ #we dont need b here since it just shifts the flow and doesnt change the geometry of the sets.

    n = length(λ)
    C_inv = inv(C)

    O = Matrix{Float64}[] #this will store the sequence of beating sets as matrices where the rows are the linear constraints defining the set.
    b_vecs = Vector{Float64}[] #this will store the corresponding constant terms for the affine constraints.

    #Σ_0 is defined by λ^T x + a = 0, which we can write as λ^T x = -a. 
    O_current = reshape(λ', 1, :) #constraint matrix for Σ_0
    b_current = [-a] #constant term for Σ_0

    push!(O, O_current)
    push!(b_vecs, b_current)

    prev_rank = rank(O_current)
    k_∞ = 0 #tracks when blocking occurs

    for k in 1:n 
        #linear constraints for λ^T(C^{-1}^k)
        new_row = (λ' * (C_inv^k))

        #for x to be in Σ_k, then Cx+κ must be in Σ_k-1
        #This takes reset constant κ back through powers of C_inv
        new_offset = -a - sum(λ' * (C_inv^j) * κ for j in 1:k)

        #append new row of O matrix
        O_next = vcat(O_current, new_row)

        #append new offset to constant b vector 
        b_next = vcat(b_current, new_offset)

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
    

function is_trivially_blocking_affine(sys::HybridAffineSystem)
    #Perform beating and blocking sets analysis to find constraint matrix and offset
    analysis = beating_and_blocking_sets_affine(sys)

    #dim of state space
    n = length(sys.λ)
     
    #Extract final constraint matrix and affine offset
    O_∞ = analysis.blocking_set
    b_∞ = analysis.blocking_offsets

    #The blocking set is a single point
    is_point = (rank(O_∞) == n)

    #This point needs to be the origin. We check if the norm of the offset is essentially zero
    is_at_origin = isapprox(norm(b_∞), 0, atol=1e-12)

    #both must be true to be trivially blocking
    return is_point && is_at_origin
end

function basis_beating_and_blocking_sets(sys::HybridAffineSystem)
    #run function we already have to get constraint matrices
    analysis = beating_and_blocking_sets_affine(sys)

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



end