# A collection of miscellaneous ODE integrators
# For all that follows:
#  1) f::Function is the vector field
#  2) z::Vector is the current state
#  3) h::Float is the step size (if fixed)
#  4) t::Float is the current time

## Single step, fully explicit methods

# Forward Euler
function forward_euler_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat)
    return z .+ h*f(z, t)
end

# Modified (fully explicit) Trapezoid Rule
function modified_trap_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat)
    z_guess = z .+ h*f(z)
    return z .+ 1/2*h*( f(z, t) + f(z_guess, t+h) )
end

# Modified (fully explicit) Midpoint Rule
function modified_midpoint_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat)
    z_guess = z .+ h/2*f(z, t)
    return z .+ h*f(z_guess, t+h/2)
end

## Adaptive Runge-Kutta methods

# A helper function to determine the adapted step size
function updated_step(LTE::AbstractFloat, tol::AbstractFloat, h::AbstractFloat, n::Int)
    # The safety parameters
    facmax = 3.
    facmin = 1/3
    fac    = 0.9
    # The predicted multiplier
    ε = abs( tol / LTE ) ^ (1/n)
    # The updated step
    return h * minimum( [ facmax, maximum( [ facmin, fac*ε ] ) ] )
end

# Runge-Kutta 23
function rk_23_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat, tf::AbstractFloat)
    # As this is an adaptive step solver, h is the step size from the previous step
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    h = minimum([h, tf-t])
    # Set the tolerence
    tol = 1e-4

    # Loop through to find an acceptable step
    while true
        # Compute the two predictions and their difference
        k1 = f(z, t)
        k2 = f(z+h*k1, t+h)
        k3 = f(z+h/4*(k1+k2), t+h/2)
        z1_3 = z + h*(1/6*k1+1/6*k2+2/3*k3)
        z1_2 = z + h*(1/2*k1+1/2*k2)
        LTE = la.norm(z1_2 - z1_3)
        # Reject or accept?
        h_new = updated_step(LTE, tol, h, 3)
        if LTE < tol
            return z1_2, h, h_new
        else
            # We reject and repeat the loop with an updated step
            h = h_new
        end
        if h < 1e-12
            @warn "Step size has decreased below 1e-12"
        end
    end
end

# Runge-Kutta 45
function rk_45_step(f::Function, z::Vector, h::AbstractFloat, t::AbstractFloat, tf::AbstractFloat)
    # As this is an adaptive step solver, h is the step size from the pervious step
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    h = minimum([h, tf-t])
    # Set the tolerence
    tol = 1e-4

    # Loop through to find an acceptable step
    while true
        # Compute the two predictions and their difference
        k1 = f(z, t)
        k2 = f(z+h*1/5, t+h*1/5)
        k3 = f(z+h*(3/40*k1+9/40*k2), t+h*3/10)
        k4 = f(z+h*(44/45*k1-56/15*k2+32/9*k3), t+h*4/5)
        k5 = f(z+h*(19372/6561*k1-25360/2187*k2+64448/6561*k3-212/729*k4),t+h*8/9)
        k6 = f(z+h*(9017/3168*k1-355/33*k2+46732/5247*k3+49/176*k4-5105/18656*k5), t+h)
        k7 = f(z+h*(35/384*k1+0*k2+500/1113*k3+125/192*k4-2187/6784*k5+11/84*k6),t+h)
        # The two updates
        z1_4 = z + h*k7
        z1_5 = z + h*(5179/57600*k1 + 0*k2 + 7571/16695*k3 + 393/640*k4 - 92097/339200*k5 + 187/2100*k6 + 1/40*k7)
        LTE = la.norm(z1_4 - z1_5)
        # Reject or accept?
        h_new = updated_step(LTE, tol, h, 5)
        if LTE < tol
            return z1_4, h, h_new
        else
            h = h_new
        end
        if h < 1e-12
            @warn "Step size has decreased below 1e-12"
        end
    end
end

