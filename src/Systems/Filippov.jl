
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
struct FilippovSol{T, X, DX, S} <: AbstractHybridSolution
    t::T    # Time data
    x::X    # Position data
    dx::DX  # f(x) Derivative at each state x - only filled when dense_out = true
    s::S    # Time indices while sliding (this still needs added)
end

# Constructor for solution struct
function FilippovSol(T, X, DX; S = NaN)
    return FilippovSol(T, X, DX, S)
end

# Initialize solution struct
function initsol(prob::prob{<:FilippovSystem, I, T}) where {I, T}
    return FilippovSol([prob.tspan[1]], [prob.init], Vector{Vector{Float64}}())
end

# INTERNAL
# Returns a vector field function at state `x` for a Filippov system
function filippov_vector_field(sys, x;
        Ftol=1e-4,
        atol=1e-4)

    F = sys.F
    G = sys.G
    H = sys.h
    N = sys.N

    h = H(x)

# Away from the guard
    if h > Ftol
        return F, false
    elseif h < -Ftol
        return G, false
    end

# Near the guard

    a = dot(N(x), F(x))
    b = dot(N(x), G(x))

    if a < -atol && b > atol
        λ = a/(a-b)
        return y -> (1-λ)*F(y) + λ*G(y), true

    elseif a > atol && b < -atol
        λ = a/(a-b)
        return y -> (1-λ)*F(y) + λ*G(y), true

    elseif a > atol && b > atol
        return F, false

    elseif a < -atol && b < -atol
        return G, false
    else
        error("Failed vector field determination")
    end
end

function guard(sys::FilippovSystem, x)
    if isnothing(sys.h)
        return nothing
    else
        return sys.h(x)
    end
end

# EXTERNAL
# Filippov-specific solve
function solve(prob::prob{S,I,T};
    solver::AbstractODESolver=RK4(),
    dense_out = true,
    dt_initial = 0.01, max_iter = 10^6, tol = 1e-6, 
    kwargs...) where {S<:FilippovSystem, I, T}
    
    sys = prob.sys
    sol = initsol(prob)

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
        # if sliding_now && !sliding_prev
        #     @warn "Sliding mode entered at t = $(tₖ)"
        # end

        sliding_prev = sliding_now

        vf(x,t) = vf_fun(x)
        x_predict, _, _ = take_step(solver, prob, vf, xₖ, tₖ, dt_step, tol, sol)

        tₖ += dt_step
        xₖ = x_predict
        if dense_out
            push!(sol.dx, vf(xₖ, tₖ))
        end
        push!(sol.t, tₖ)
        push!(sol.x, xₖ)

       

    end

    return sol
end