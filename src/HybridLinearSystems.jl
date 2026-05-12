# I want a hybrid linear system and hybrid affine structures


# I want the following:
# 1. Generate trajectories
# 2. Determine whether or not the system is trivially blocking
# 3. Find the beating/blocking sets
# 4. (For affine) determine whether or not the system is Zeno and find the Zeno time

module HybridDynamics

using ForwardDiff
using LinearAlgebra


export HybridLinearSystem, solve_hybrid_linear_system, solve_hybrid_linear_system_exp
export forward_euler_step, modified_trap_step, modified_midpoint_step, rk_23_step
export beating_and_blocking_sets, is_trivially_blocking, basis_beating_and_blocking_sets


struct HybridLinearSystem
    A::Matrix{Float64}
    λ::Vector{Float64}
    C::Matrix{Float64}
end

struct HybridLinearProblem
    sys::HybridLinearSystem
    x₀::Vector{Float64}
    tspan::Tuple{Float64,Float64}
end

include("ODE_solvers.jl")

function solve_hybrid_linear_system(problem::HybridLinearProblem, dt_initial::Float64; step_method=forward_euler_step, is_adaptive=false)
    
    #deconstruct the system structure.
    sys = problem.sys
    A, λ, C = sys.A, sys.λ, sys.C

    #unpacks the time span into start and end times.
    t_start, t_end = problem.tspan
    
    #Initialize the time and output vectors starting with the initial conditions.
    t = [t_start]
    x = [problem.x₀]

    #empty list to store timestamps where jumps occur
    jump_times = Float64[]
    jump_indices = Int[]

    #Defines the vector field for continuous part
    f(x, t) = A * x
    
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

        #Calcs timestamp for predicted state. 
        t_next = tₖ + Δt_used

        #Evaluates the switching function at start of next step
        h_now = dot(λ, xₖ)
        
        #Evalutes the switching function at end of next step
        h_next = dot(λ, x_predict)

        #JUMP CONDITION
        if h_now * h_next < 0 && t_next < t_end

            #Calcs linear interpolation. We determine how far through the time step the corssing actually occurred as I didnt do this at first. 
            θ = -h_now / (h_next - h_now)

            #This calcs the exact time t_star where the trajectory hit the gaurd. 
            t_star = tₖ + θ * Δt_used

            #Interpolates the state vector to find its value at t_star, directly before the jump (x^-)
            x⁻ = xₖ + θ * (x_predict - xₖ) # Pre-jump state at t_star

            #Logs the state right at moment of impact
            push!(x, x⁻) 
            push!(t, t_star)
            
            #JUMP 
            #Applies C to the impact state and logs in t_star. This is the jump
            x⁺ = C * x⁻ # Post-jump state at t_star

            push!(x, x⁺)
            push!(t, t_star)

            #Updates event logs with pinpointed time and current index in the array. 
            push!(jump_times, t_star)
            push!(jump_indices, length(x))
            
            #Resets step size. Since a jump might happen early in a step, we resume with full sized step from new jump position. 
            Δt = dt_initial
        else 

            #if no crossing detected, the pred state and time are just appended to output. 
            push!(x, x_predict)
            push!(t, t_next)

            #updates step size for next iteration (mainly for adaptive solvers)
            Δt = Δt_next
        end
    end
    #returns the time history, state history, and specific detais of discrete jump events detected along the way.
    return t, x, jump_times, jump_indices
end

function solve_hybrid_linear_system_exp(problem::HybridLinearProblem, dt_initial::Float64; tol=1e-10)
    sys = problem.sys
    A,λ,C = sys.A, sys.λ, sys.C

    t_start, t_end = problem.tspan

    t = [t_start]
    x = [problem.x₀]

    jump_times = Float64[]
    jump_indices = Int[]

    #set initial step size
    Δt = dt_initial

    #precompute matrix exp for w/ step size. 
    Φ = exp(A * dt_initial)

    while t[end] < t_end

        #shrink step size if we are going to overshoot t_end
        if t[end] + Δt > t_end
            Δt = t_end - t[end] 
            Φ_step = exp(A * Δt) #compute final exp for final step
        else 
            Φ_step = Φ
        end

        #get most recent state and time
        xₖ = x[end]
        tₖ = t[end]

        #continuous evo using x_pred = e^(A * Δt) * xₖ
        x_predict = Φ_step * xₖ

        t_next = tₖ + Δt

        #eval the switching function at start and predicted end
        h_now = dot(λ, xₖ)
        h_next = dot(λ, x_predict)

        #jump condition
        if h_now * h_next < 0 && t_next <= t_end

            #bisection method root finding (maybe do others?)
            τ_left = 0.0
            τ_right = Δt
            h_left = h_now

            #loop until interval is smaller than tolerance
            while(τ_right - τ_left) > tol
                τ_mid = (τ_left + τ_right) / 2.0

                #eval state at midpoint
                x_mid = exp(A * τ_mid) * xₖ
                h_mid = dot(λ, x_mid)

                #narrow bis window
                if h_left * h_mid <= 0
                    τ_right = τ_mid
                else 
                    τ_left = τ_mid
                    h_left = h_mid
                end
            end

            #exact time int from tk to jump
            τ_star = (τ_left + τ_right) / 2.0
            t_star = tₖ + τ_star

            #calc exact state directly before jump
            x⁻ = exp(A * τ_star) * xₖ

            #log state at moment of impact
            push!(x,x⁻)
            push!(t,t_star)

            #jump, apply C to impact to get post jump state. 
            x⁺ = C * x⁻

            push!(x,x⁺)
            push!(t,t_star)

            #update event 
            push!(jump_times,t_star)
            push!(jump_indices, length(x))

            #reset step size
            Δt = dt_initial

        else 
            #if no crossing, append exact state 
            push!(x,x_predict)
            push!(t,t_next)

            #reset step in case shrunk at bd
            Δt = dt_initial
        end
    end

    return t,x,jump_times,jump_indices
end

function beating_and_blocking_sets(sys::HybridLinearSystem)
    #Deconstructs the system structure to get the matrices and vectors we need.
    A, λ, C = sys.A, sys.λ, sys.C
    
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

    #We iterate n times. There is a theorem about this I think
    for k in 1:n 
        #Calc the kth constraint row. This maps the guard back through k jumps. 
        new_row = (λ' * (C_inv^k))

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

function is_trivially_blocking(sys::HybridLinearSystem)

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

function basis_beating_and_blocking_sets(sys::HybridLinearSystem)
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

end