#Function calcs the continuous time derivative for the augmented state. Runs the usual dx and dΦ.

#CURRENTLY THIS FUNCTION IS WORTHLESS BUT I WANT IT TO BE USEFUL LATER SO I WILL KEEP IT - DS
function variational_vector_field(f, U::AbstractMatrix, t)
    # Unpack state
    x = U[:, 1]
    Φ = U[:, 2:end]

    # Base dynamics: x' = f(x, p, t)
    dx = f(x, t)

    # Variational dynamics: Φ' = A(t)Φ
    A = ForwardDiff.jacobian(y -> f(y, t), x)
    dΦ = A * Φ

    # Return augmented derivative
    return hcat(dx, dΦ)
end

#Calcs the Δ_*^f at a boundary. We look how how much the vf mismatches before and after the jump (f⁺/f⁻), how much it scales things (DΔ⁻), and how the trajectory hits boundary (dh⁻).
#Then we create a matrix that takes all of this.
function compute_pushforward(f, Δ, h_guard, x⁻, t)
    n = length(x⁻)
    Id = I(n)

    # Eval field at boundaries (using p for parameters)
    f⁻ = f(x⁻, t)
    x⁺ = Δ(x⁻, t)
    f⁺ = f(x⁺, t)

    # Compute grads and jacob via ForwardDiff
    dh⁻ = ForwardDiff.gradient(h_guard, x⁻)
    DΔ⁻ = ForwardDiff.jacobian(y -> Δ(y, t), x⁻)

    # Check dh(x) * f(x) = 0
    denom = dot(dh⁻, f⁻)
    if abs(denom) < 1e-6
        @warn "Non-transversal crossing detected: Trajectory is tangent to guard surface."
    end

    # Outer prods
    term1 = Id - (f⁻ * dh⁻') ./ denom
    term2 = (f⁺ * dh⁻') ./ denom

    # Full pushforward Δᶠ_*
    Δ_star_f = DΔ⁻ * term1 + term2

    return Δ_star_f
end

#When our solver hits h(x)=0 this is called. We get the pushforward matrix and multiply by Φ⁻ to get new Φ⁺. 
function apply_variational_jump(U::AbstractMatrix, f, Δ, h_guard, t)
    x⁻ = U[:, 1]
    Φ⁻ = U[:, 2:end]

    # Compute the pushforward before state updates
    Δ_star_f = compute_pushforward(f, Δ, h_guard, x⁻, t)

    # Apply disc jump to base state x⁺ = Δ(x⁻)
    x⁺ = Δ(x⁻, t)

    # Apply pf mapping to fund matrix: Φ⁺ = Δ_*^f * Φ⁻
    Φ⁺ = Δ_star_f * Φ⁻

    # Update state vector in-place
    U[:, 1] .= x⁺
    U[:, 2:end] .= Φ⁺
    
    return U
end