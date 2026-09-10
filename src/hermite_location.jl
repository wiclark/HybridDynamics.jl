# Uses third-order Hermite interpolation to determine whether or not a crossing happened
# This should replace 'evaluate_crossing'
function evaluate_crossing(h_now, h_next, hp_now, hp_next, t_now, t_next, direction::Int; tol=1e-6)
    # The first check is simply the linear crossing - this guarantees a crossing by the MVT
    valid_linear(h1, h2) = 
        (direction == 0 && h1 * h2 < 0) ||
        (direction == -1 && h1 > 0 && h2 < 0) ||
        (direction == 1 && h1 < 0 && h2 > 0)
    # The next function attempts to find an event based off of Hermite interpolation
    function valid_hermite()
        Δt = t_next - t_now
        α = -1/Δt^2*(3*h_now-3*h_next+2*Δt*hp_now+Δt*hp_next)
        β = 1/Δt^3*(2*h_now-2*h_next+Δt*hp_now+Δt*hp_next)
        D = 4*α^2-12*β*hp_now
        if D < 0 # Negative discriminant <- No real roots
            return false, NaN
        else
            # Is there a (or two) root(s) within the interval?
            r1 = (-2*α+sqrt(D)) / (6*β)
            r2 = (-2*α-sqrt(D)) / (6*β)
            R1 = minimum([r1, r2])
            R2 = maximum([r1, r2])
            if 0 < R1 < Δt
                # Find the value at the first root
                # This expression uses the interpolation to find the crossing
                # I suspect that a 'smarter' version would be to use 'take_step'
                y1 = h_now + hp_now*r1 + α*r1^2 + β*r1^3
                return valid_linear(h_now, y1), R1
            elseif 0 < R2 < Δt
                # Find the value at the first root
                y1 = h_now + hp_now*r2 + α*r2^2 + β*r2^3
                return valid_linear(h_now, y1), R2
            else
                return false, NaN
            end
        end 
    end
    # Do we suspect a crossing?
    if valid_linear(h_now, h_next)
        t_root = t_now - h_now * (t_next - t_now) / (h_next - h_now)
        return true, t_root, NaN
    elseif valid_hermite()[1]
        return true, valid_hermite()[2], NaN
    else
        return false, NaN, NaN
    end
end