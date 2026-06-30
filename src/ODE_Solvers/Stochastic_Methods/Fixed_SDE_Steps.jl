# This is an implementation of the Euler–Maruyama method
struct EulerMaruyama <: STO end

# This is for the Itô SDE: dx = f(x, t)dt + g(x, t)dW, where g is a n×m matrix function where m is the dim on Brownian motion

# Single take_step for stochastic methods
# Fixed step methods
# Helper function to take the step using multiple dispatch
compute_step(::EulerMaruyama, f, g, x, Δt, t) = euler_maruyama_step(f, g, x, Δt, t)

# Note sol is not used, we do this to make using the function easier.
function take_step(solver::STO, prob::prob{S, I, T}, f, g, xₖ, tₖ, Δt, tol, sol; check=true, guard_direction = default_guard_direction(prob.sys)) where {S<:StochasticSystem, I, T}
    sys = prob.sys
    x_predict = compute_step(solver, f, g, xₖ, Δt, tₖ)
    x_mid     = compute_step(solver, f, g, xₖ, Δt/2.0, tₖ)

    if check
        # Evaluate Guards
        h_now  = guard(sys, xₖ)
        h_mid  = guard(sys, x_mid)
        h_next = guard(sys, x_predict)

        #Use cross guard check
        eventtrigger, t_root, _ = crossed_guard(sys, h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol, direction=guard_direction)

        return x_predict, eventtrigger, t_root, Δt, Δt
    else
        return x_predict, false, NaN, Δt, Δt
    end
end

function euler_maruyama_step(f::Function, g::Function, z::AbstractArray, Δt::AbstractFloat, t::AbstractFloat)
    # The state and diffusion term
    fₖ, gₖ = f(z, t), g(z, t)
    # Check sizes
    if length(fₖ) ≠ size(gₖ)[1]
        @warn "Dimensions of drift and diffusion are not compatable"
    end
    # The noise term
    noise = Normal(0, √(Δt))
    ΔW = rand(noise, size(gₖ)[2])
    # The update
    return z .+ Δt*fₖ .+ gₖ*ΔW
end