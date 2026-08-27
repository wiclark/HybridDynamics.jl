# This way of calculating the deriv is redundant, but it would take a lot of rewriting to get f out of take_step.
# To add dense output to a solve dispatch, add the following:
# if dense_out
#     push!(sol.dx, f(xₖ, tₖ))
# end

function (sol::StochasticSol)(t::Real)
    error("Interpolation attempted on stochastic system type; this is undefined.")
end

# Cubic Hermite interpolant for dense output
function (sol::AbstractHybridSolution)(t::Real)

    # Safties to make sure it usually works
    isempty(sol.dx) && error("Solution struct did not return dense output")

    if t < sol.t[1]
        if abs(t-sol.t[1]) > 1e-6
            @warn("BoundsError(sol, t). Out of bounds: constant extrapolation used.")
        end
        return sol.x[1]
    end
    t == sol.t[1] && return sol.x[1]
    t == sol.t[end] && return sol.x[end]
    if t > sol.t[end]
        if abs(t-sol.t[end]) > 1e-6
            @warn("BoundsError(sol, t). Out of bounds: constant extrapolation used.")
        end
        return sol.x[end]
    end
    
    # Find the index of interest
    j = searchsortedfirst(sol.t, t)
    i = j - 1
    
    t₁, t₂ = sol.t[i], sol.t[i+1]
    x₁, x₂ = sol.x[i], sol.x[i+1]
    f₁, f₂ = sol.dx[i], sol.dx[i+1]

    Δt = t₂ - t₁
    θ = (t - t₁)/Δt

    # This is faster than matrix operations
    A1 = 2θ^3 - 3θ^2 + 1
    A2 = θ^3 - 2θ^2 + θ
    A3 = -2θ^3 + 3θ^2
    A4 = θ^3 - θ^2

    return A1*x₁ + A2*Δt*f₁ + A3*x₂ + A4*Δt*f₂
end