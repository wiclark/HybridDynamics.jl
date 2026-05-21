
function CreateSystem(A::AbstractMatrix, λ::AbstractVector, C::AbstractMatrix)
    #I am using Abstract matrix here as I believe it will make it easier to avoid type errors when "users" use the functions to see how the structs work. It might be worth just making my structs use this but I am not sure if that matters. 
    #this same reason is why I do Float64.(A), etc as I think it will be best to avoid errors. The issue would be if smth like [1 0; 0 1] gets typed as it would be Matrix{Int64} and not Matrix{Float64} hence why I did it this way
    #I do this in all the struct functions, I cant imagine rewriting these comments is better than just writing that here for now. When I write the inline thingy I will do it there. 
    return HybridLinearSystem(Float64.(A), Float64.(λ), Float64.(C))
end
struct HybridLinearSystem
    A::Matrix{Float64} 
    λ::Vector{Float64} 
    C::Matrix{Float64}
end

function CreateSystem(A::AbstractMatrix, b::AbstractVector, λ::AbstractVector, a::Real, C::AbstractMatrix, κ::AbstractVector)
    #different from linear as it inputs more stuff. That should work for what we need. 
    return HybridAffineSystem(Float64.(A), Float64.(b), Float64.(λ), Float64(a), Float64.(C), Float64.(κ))
end
struct HybridAffineSystem
    A::Matrix{Float64} #Flow matrix
    b::Vector{Float64} #Flow Constant
    λ::Vector{Float64} #Guard normal
    a::Float64         #Gaurd Constant
    C::Matrix{Float64} #Reset Matrix
    κ::Vector{Float64} #Reset Constant
end

function CreateProblem(sys::HybridLinearSystem, x₀::AbstractVector, tspan::Tuple{Real, Real})
    return HybridLinearProblem(sys, Float64.(x₀), (Float64(tspan[1]), Float64(tspan[2])))
end
struct HybridLinearProblem
    sys::HybridLinearSystem 
    x₀::Vector{Float64}
    tspan::Tuple{Float64,Float64}
end

function CreateProblem(sys::HybridAffineSystem, x₀::AbstractVector, tspan::Tuple{Real, Real})
    return HybridAffineProblem(sys, Float64.(x₀), (Float64(tspan[1]), Float64(tspan[2])))
end
struct HybridAffineProblem
    sys::HybridAffineSystem
    x₀::Vector{Float64}
    tspan::Tuple{Float64,Float64}
end

function CreateSolution(prob::HybridLinearProblem, t::AbstractVector, x::AbstractVector, jump_times::AbstractVector, jump_indices::AbstractVector)
    return HybridLinearSolution(Float64.(t), Vector{Float64}.(x), Float64.(jump_times), Int.(jump_indices))
end
struct HybridLinearSolution
    t::Vector{Float64}          #array storing cont time data
    x::Vector{Vector{Float64}}  #array storing cont x data 
    jump_times::Vector{Float64} #how many jump
    jump_indices::Vector{Int}   #when jump
end

function CreateSolution(prob::HybridAffineProblem, t::AbstractVector, x::AbstractVector, jump_times::AbstractVector, jump_indices::AbstractVector)
    return HybridAffineSolution(Float64.(t), Vector{Float64}.(x), Float64.(jump_times), Int.(jump_indices))
end
struct HybridAffineSolution
    t::Vector{Float64}
    x::Vector{Vector{Float64}}
    jump_times::Vector{Float64}
    jump_indices::Vector{Int}
end

function crossed_guard(h_now, h_next; tol=1e-12)
    #we check if switch function returns negative sign meaning we crossed guard. Edge case if we are ON guard. 
    return (h_now * h_next < 0) || (abs(h_now) <= tol && h_next < -tol)
end

guard(sys::HybridLinearSystem, x) = dot(sys.λ, x) #eval distance/orientation of current state relative to guard
guard(sys::HybridAffineSystem, x) = dot(sys.λ, x) + sys.a #same as above but for affine

apply_reset(sys::HybridLinearSystem, x) = sys.C * x #disc transition when we jump
apply_reset(sys::HybridAffineSystem, x) = sys.C * x + sys.κ #above

function interpolate_state(x₀, x₁, θ)
    return x₀ + θ * (x₁ - x₀)
end

function locate_guard_crossing(xₖ, x_predict, h_now, h_next, tₖ, Δt)
    θ = -h_now / (h_next - h_now)
    t_star = tₖ + θ * Δt
    x_star = interpolate_state(xₖ, x_predict, θ)
    return t_star, x_star
end

struct LinearFlow{TM<:AbstractMatrix, TV<:AbstractVector}
    V::TM
    Λ::TV
    V_inv::TM
end

function LinearFlow(A)
    eig = eigen(A) #eigen decomp for A

    #getting A=VΛV^{-1} to compute cont evo of system exactly (mostly)
    return LinearFlow(eig.vectors, eig.values, inv(eig.vectors))
end

function flow(flowmap::LinearFlow, τ, x) #provides exact cont state any any tiny step tau for exact sol
    y = flowmap.V_inv * x 
    return real.(flowmap.V * (exp.(flowmap.Λ .* τ) .* y)) 
end

function find_crossing_bisection(sys::HybridLinearSystem, xₖ, Δt, h_now; tol=1e-12, max_iter=100)
    τ_l, τ_r = 0.0, Δt #sets left and right bounds of search window
    h_l = h_now
    for _ in 1:max_iter
        if (τ_r - τ_l) < tol break end
        τ_m = (τ_l + τ_r) / 2.0 #midpoint
        x_m = xₖ + τ_m * (sys.A * xₖ)  #what state is at mp
        signbit(h_l) != signbit(guard(sys, x_m)) ? τ_r = τ_m : (τ_l = τ_m; h_l = guard(sys, x_m)) #if sign of guard eval at left is different then midpoint, crossing happened in first half, otherwise its second half. then we update
    end
    return τ_l
end
function find_crossing_bisection(sys::HybridAffineSystem, xₖ, Δt, h_now; tol=1e-12, max_iter=100)
    τ_l, τ_r = 0.0, Δt
    h_l = h_now
    for _ in 1:max_iter
        if (τ_r - τ_l) < tol break end
        τ_m = (τ_l + τ_r) / 2.0
        x_m = xₖ + τ_m * (sys.A * xₖ + sys.b) 
        signbit(h_l) != signbit(guard(sys, x_m)) ? τ_r = τ_m : (τ_l = τ_m; h_l = guard(sys, x_m))
    end
    return τ_l
end
function find_crossing_bisection_exp(sys, flowmap, xₖ, Δt, h_now; tol=1e-12, max_iter=100)
    τ_l, τ_r = 0.0, Δt
    h_l = h_now
    for _ in 1:max_iter
        if (τ_r - τ_l) < tol break end
        τ_m = (τ_l + τ_r) / 2.0
        x_m = flow(flowmap, τ_m, xₖ)
        signbit(h_l) != signbit(guard(sys, x_m)) ? τ_r = τ_m : (τ_l = τ_m; h_l = guard(sys, x_m))
    end
    return τ_l
end


function check_beating_status(sys, instant_jumps, n, x_current, t_current, tol)
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
#--------------------------------------------
#SOLVERS
function solve_hybrid_system(problem::HybridLinearProblem, dt_initial::Float64; step_method=forward_euler_step, is_adaptive=false, max_iter=10^6, tol = 1e-12)
    
    #deconstruct the system structure.
    sys = problem.sys
    A = sys.A

    n = size(sys.A, 1) #dimension for beating check

    #unpacks the time span into start and end times.
    t_start, t_end = problem.tspan
    
    #Initialize the time and output vectors starting with the initial conditions.
    t = [t_start]
    x = [problem.x₀]

    sizehint!(t, max_iter)
    sizehint!(x, max_iter)

    #empty list to store timestamps where jumps occur
    jump_times = Float64[]
    jump_indices = Int[]

    iter = 0

    instant_jumps = 0
    last_jump_t = -Inf #track time of  last jump

    #Defines the vector field for continuous part
    f(x, t) = A * x
    
    #Sets initial time step size. This can change if the solver we use is adaptive or needs to hit a specific time boundary. 
    Δt = dt_initial

    #This main loop continuous as long as the most recent recorded time is less than final
    while t[end] < t_end

        #prevent problems running away. Oh the issues this caused...
        iter += 1
        if iter > max_iter
            @warn "Possible Zeno detected"
            break
        end

        if t_end - t[end] <= eps(t_end) * 10
            break
        end

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
        h_now = guard(sys,xₖ)
        
        #Evalutes the switching function at end of next step
        h_next = guard(sys, x_predict)

        #JUMP CONDITION
        if crossed_guard(h_now, h_next; tol=tol) && t_next < t_end

            τ_star = find_crossing_bisection(sys, xₖ, Δt_used, h_now; tol=tol)
            t_star = tₖ + τ_star
            x⁻ = xₖ + τ_star * f(xₖ, tₖ)

            if abs(t_star - last_jump_t) < tol
                instant_jumps += 1
            else 
                instant_jumps = 1 
            end
            last_jump_t = t_star

            if check_beating_status(sys, instant_jumps, n, x⁻, t_star, tol) != :continue
                break
            end

            #Logs the state right at moment of impact
            push!(x, x⁻) 
            push!(t, t_star)
            
            #JUMP 
            #Applies C to the impact state and logs in t_star. This is the jump
            x⁺ = apply_reset(sys, x⁻)

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
    return HybridLinearSolution(t,x,jump_times,jump_indices)
end
function solve_hybrid_system(problem::HybridAffineProblem, dt_initial::Float64; step_method=forward_euler_step, is_adaptive=false, max_iter=10^6, tol = 1e-12)
    
    #Deconstruct the system structure.
    sys = problem.sys
    n = size(sys.A, 1)

    #unpacks time span
    t_start, t_end = problem.tspan

    #Initialize time and output vectors starting with initial conditions.
    t = [t_start]
    x = [problem.x₀]

    #empty list to store timestamps where jumps occur
    jump_times = Float64[]
    jump_indices = Int[]

    instant_jumps = 0
    last_jump_t = -Inf
    iter = 0

    #Sets initial time step size. This can change if the solver we use is adaptive or needs to hit a specific time boundary.
    Δt = dt_initial
    
    #Defines the vector field for continuous part
    f(x, t) = sys.A * x + sys.b

    #This main loop continuous as long as the most recent recorded time is less than final
    while t[end] < t_end
        iter += 1
        if iter > max_iter 
            @warn "Max iterations reached possible Zeno"
            break
        end

        if t_end - t[end] <= eps(t_end) * 10
            break
        end

        Δt_step = (t[end] + Δt > t_end) ? (t_end - t[end]) : Δt

        #Grabs most recent state and time. This lets us always start from latest point, including after jumps.
        xₖ = x[end]
        tₖ = t[end]

        #Logic for solver type
        if is_adaptive
            #If adaptive, the step_method  returns new state, the actual step size is used and a suggestion for the next step size.
            x_predict, Δt_used, Δt_next = step_method(f, xₖ, Δt_step, tₖ, t_end)
        else 
            #Execute numerical integration step (e.g. forward Euler) to predict the state at next time.
            x_predict = step_method(f, xₖ, Δt_step, tₖ)
            Δt_used = Δt_step
            Δt_next = dt_initial
        end

        #calc timestamp for predicted state
        t_next = tₖ + Δt_used

        #evaluates the affine version at start of next step
        h_now = guard(sys, xₖ)

        #evaluates the affine version at end of next step
        h_next = guard(sys, x_predict)

        #Jump condition
        if crossed_guard(h_now, h_next; tol=tol) && t_next <= t_end

            τ_star = find_crossing_bisection(sys, xₖ, Δt_used, h_now; tol=tol)
            t_star = tₖ + τ_star
            x⁻ = xₖ + τ_star * f(xₖ, tₖ)

            if abs(t_star - last_jump_t) < tol
                instant_jumps += 1
            else 
                instant_jumps = 1
            end
            last_jump_t = t_star

            if check_beating_status(sys, instant_jumps, n, x⁻, t_star, tol) != :continue
                break
            end

            #logs the state at moment of impact 
            push!(x, x⁻)
            push!(t, t_star)

            #JUMP
            #applies C and \kappa to the impact state and logs in t_star. 
            x⁺ = apply_reset(sys,x⁻) #post jump state at t_star

            push!(x, x⁺)
            push!(t,t_star)

            #update jump times and indices
            push!(jump_times, t_star)
            push!(jump_indices, length(x))

            #resets step size. 
            Δt = dt_initial
        else 
            #if no jump, the pred state and time are just appended to output.
            push!(x, x_predict)
            push!(t, t_next)

            #update step size if adaptive
            Δt = Δt_next
        end
    end
    return HybridAffineSolution(t, x, jump_times, jump_indices)
end

function solve_hybrid_system_exp(problem::HybridLinearProblem, dt_initial::Float64; tol=1e-10, max_iter=10^6)
    sys = problem.sys
    A = sys.A
    n = size(sys.A, 1)
    t_start, t_end = problem.tspan

    instant_jumps = 0
    last_jump_t = -Inf

    t = [t_start]
    x = [problem.x₀]

    jump_times = Float64[]
    jump_indices = Int[]

    #set initial step size
    Δt = dt_initial

    flowmap = LinearFlow(A) #before we run we calc matrix to make calc easier as I had issues. 

    Φ = real.(flowmap.V * Diagonal(exp.(flowmap.Λ .* dt_initial)) * flowmap.V_inv)

    iter = 0

    while t[end] < t_end

        #new max iteration check to avoid horrible errors
        iter += 1
        if iter > max_iter
            @warn "Maximum Iteration Count ($max_iter) exceeded. "
            break
        end

        if t_end - t[end] <= eps(t_end) * 10
            break
        end

        #shrink step size if we are going to overshoot t_end
        Δt_step = (t[end] + Δt > t_end) ? (t_end - t[end]) : Δt
        
        #check if loop had to shrink time step if yes recalc a transition matrix for that small slice. Otherwise we move on
        if Δt_step != Δt
            Φ_step = real.(flowmap.V * Diagonal(exp.(flowmap.Λ .* Δt_step)) * flowmap.V_inv)
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
        h_now = guard(sys, xₖ)
        h_next = guard(sys, x_predict)

        #jump condition
        if crossed_guard(h_now, h_next; tol=tol) && t_next <= t_end

            #once guard cross detected use bisection to find impact and calc exact state using flow. 
            τ_star = find_crossing_bisection_exp(sys, flowmap, xₖ, Δt_step, h_now; tol=tol)
            t_star = tₖ + τ_star
            x⁻ = flow(flowmap, τ_star, xₖ)

            if abs(t_star - last_jump_t) < tol
                instant_jumps += 1
            else 
                instant_jumps = 1
            end
            last_jump_t = t_star

            #calc exact state directly before jump
            x⁻ = flow(flowmap, τ_star, xₖ)

            status = check_beating_status(sys, instant_jumps, n, x⁻, t_star, tol)
            if status != :continue
                break
            end

            #jump, apply C to impact to get post jump state. 
            x⁺ = apply_reset(sys, x⁻)
            append!(x, [x⁻, x⁺])
            append!(t, [t_star, t_star])

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

    return HybridLinearSolution(t,x,jump_times,jump_indices)
end
#--------------------------------------------
#BEATING AND BLOCKING SETS
function beating_and_blocking_sets(sys::HybridLinearSystem)
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
function beating_and_blocking_sets(sys::HybridAffineSystem)
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
function is_trivially_blocking(sys::HybridAffineSystem)
    #Perform beating and blocking sets analysis to find constraint matrix and offset
    analysis = beating_and_blocking_sets(sys)

    #dim of state space
    n = length(sys.λ)
     
    return rank(analysis.blocking_set) == n && isapprox(norm(analysis.blocking_offsets), 0, atol=1e-12)
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
function basis_beating_and_blocking_sets(sys::HybridAffineSystem)
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

