######
# Cubic Hermite interpolant for dense output 
# CK: I don't know where to put this. This also still needs to be integrated into solve dispatches
function (sol::AbstractHybridSolution)(t::AbstractFloat)
    t_data = sol.t
    x_data = sol.x
    f_data = sol.dx

    if isempty(sol.dx)
        error("Solution struct did not return dense output")
    end

    # The first step is to make sure that t∈t_data
    if t < t_data[1] || t > t_data[end]
        @warn "Time is out of bounds"
        return NaN
    end

    # Next, determine the interval t lives in
    idx_first = searchsortedfirst(t_data, t) - 1
    idx_second = idx_first + 1
    # Gather all of the useful information
    t₁, t₂ = t_data[idx_first], t_data[idx_second]
    x₁, x₂ = x_data[idx_first], x_data[idx_second]
    f₁, f₂ = f_data[idx_first], f_data[idx_second]
    Δt = t₂ - t₁

    # Determine the coefficients
    Aⁱ = [2/Δt^3 -2/Δt^3  1/Δt^2 1/Δt^2;
         -3/Δt^2  3/Δt^2 -2/Δt  -1/Δt;
          0       0       1      0;
          1       0       0      0]
    # bⁱ = [x₁, x₂, f₁, f₂]
    # Update to work for vectors
    bⁱ = vcat(x₁', x₂', f₁', f₂')
    #α, β, γ, δ = Aⁱ * bⁱ
    C = Aⁱ * bⁱ
    α, β, γ, δ = C[1,:], C[2,:], C[3,:], C[4,:]
    ts = t - t₁

    return α*ts^3 + β*ts^2 + γ*ts + δ
end