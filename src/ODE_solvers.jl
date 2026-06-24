# A collection of miscellaneous ODE integrators
# For all that follows:
#  1) f::Function is the vector field
#  2) z::Vector is the current state
#  3) h::Float is the step size (if fixed)
#  4) t::Float is the current time

#Runge Kutta / Single Step Family Solvers
#Fixed Runge Kutta parent tag
abstract type RK <: AbstractODESolver end
#Non-adaptive solvers
abstract type FixedRK <: RK end
struct ForwardEuler <: FixedRK end
struct ModifiedTrap <: FixedRK end
struct ModifiedMidpoint <: FixedRK end
struct RichardsonExtrapolation <: FixedRK end
struct RK4 <: FixedRK end
struct ImplicitEuler <: FixedRK end

#Adaptive Runge Kutta parent tag
abstract type AdaptiveRK <: RK end
#Adaptive solvers
struct RK45 <: AdaptiveRK end
struct RK23 <: AdaptiveRK end

#Linear Multistep Method Family Solvers
abstract type LMM <: AbstractODESolver end

abstract type FixedLMM <: LMM end

struct AdamsBashforth2 <: FixedLMM end
struct AdamsBashforth3 <: FixedLMM end
struct BDF2 <: FixedLMM end

abstract type AdaptiveLMM <: LMM end

struct AdaptiveABM2 <: AdaptiveLMM end
struct AdaptiveABM3 <: AdaptiveLMM end

#Other solvers Idk where to put
#Exponential Solver
struct ExponentialSolver <: AbstractODESolver end

struct MagnusLeapfrog <: FixedRK end


## Single step, fully explicit methods

# Forward Euler
function forward_euler_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    return z .+ h*f(z, t)
end

# Modified (fully explicit) Trapezoid Rule
function modified_trap_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    z_guess = z .+ h*f(z,t)
    return z .+ 1/2*h*( f(z, t) + f(z_guess, t+h) )
end

# Modified (fully explicit) Midpoint Rule
function modified_midpoint_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    z_guess = z .+ h/2*f(z, t)
    return z .+ h*f(z_guess, t+h/2)
end

#Classic Runge-Kutta 4
function rk_4_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    k1 = f(z,t)
    k2 = f(z + h/2 * k1, t + h/2)
    k3 = f(z + h/2 * k2, t + h/2)
    k4 = f(z + h * k3, t + h)

    #Shouldnt need to check guard in these inner stages here as no adpative step size.

    return z + h/6 * (k1 + 2*k2 + 2*k3 + k4)
end

#Implicit Euler (Stiff solver)
function implicit_euler_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    t_new = t + h
    #Initial Guess via exp Euler
    z_guess = forward_euler_step(f, z, h, t)

    #Imp euler 
    return implicit_newton_solve(f, z_guess, z, 1.0, h, t_new)
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
function rk_23_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat, tf::AbstractFloat, sys, tol; adaptive=true)
    # As this is an adaptive step solver, h is the step size from the previous step
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    h = minimum([h, tf-t])

    #Intialization to make sure they exist
    z2 = z; z3 = z

    # Loop through to find an acceptable step
    while true
        h_now = guard(sys, z)
        # Compute the two predictions and their difference
        k1 = f(z, t)

        z2 = z + h*k1
        k2 = f(z2, t+h)
        h2    = guard(sys, z2)

        z3 = z + h/4*(k1+k2)
        k3 = f(z+h/4*(k1+k2), t+h/2)
        h3    = guard(sys, z3)

        z1_3 = z + h*(1/6*k1+1/6*k2+2/3*k3)
        z1_2 = z + h*(1/2*k1+1/2*k2)

        #used to find true midpoint for quad guard check. 
        #This prevents solver from spinning up a whole new adaptive loop so x_mid is truly in the mid which is what we want for quad matrix. 
        if !adaptive
            return z1_2, h, h 
        end

        LTE = norm(z1_2 - z1_3)
        # Reject or accept?
        h_new = updated_step(LTE, tol, h, 3)
        
        h_end = guard(sys, z1_2)

        #did any intermediate stage cross guard?
        stage_crossed = (h_now * h2 < 0) || (h_now * h3 < 0)
        #Did final state completely miss crossing?
        end_missed = (h_now * h_end > 0)

        #if we crossed inside step but missed at end, force rejection
        if stage_crossed && end_missed
            h = h / 2.0 #force smaller step
            continue
        end

        if LTE < tol
            return z1_2, h, h_new
        else
            # We reject and repeat the loop with an updated step
            h = h_new
        end
        if h < 1e-12
            @warn "Step size has decreased below 1e-12"
            return z1_2, h, h_new #force break to avoid looping forever
        end
    end
end

# Runge-Kutta 45
function rk_45_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat, tf::AbstractFloat, sys, tol; adaptive=true)
    # As this is an adaptive step solver, h is the step size from the pervious step
    # As the step size is not of fixed size, we specify the terminal time, tf, of the problem
    h = minimum([h, tf-t])

    #Initialization to make them exist
    z2 = z; z3 = z; z4 = z; z5 = z; z6 = z

    # Loop through to find an acceptable step
    while true
        h_now = guard(sys, z)
        # Compute the two predictions and their difference
        k1 = f(z, t)

        z2 = z + h*1/5*k1
        k2 = f(z+h*1/5*k1, t+h*1/5)
        h2 = guard(sys, z2)

        z3 = z + h*(3/40*k1 + 9/40*k2)
        k3 = f(z+h*(3/40*k1+9/40*k2), t+h*3/10)
        h3 = guard(sys, z3)

        z4 = z + h*(44/45*k1 - 56/15*k2 + 32/9*k3)
        k4 = f(z+h*(44/45*k1-56/15*k2+32/9*k3), t+h*4/5)
        h4 = guard(sys, z4)

        z5 = z + h*(19372/6561*k1 - 25360/2187*k2 + 64448/6561*k3 - 212/729*k4)
        k5 = f(z+h*(19372/6561*k1-25360/2187*k2+64448/6561*k3-212/729*k4),t+h*8/9)
        h5 = guard(sys, z5)

        z6 = z + h*(9017/3168*k1 - 355/33*k2 + 46732/5247*k3 + 49/176*k4 - 5105/18656*k5)
        k6 = f(z+h*(9017/3168*k1-355/33*k2+46732/5247*k3+49/176*k4-5105/18656*k5), t+h)
        h6 = guard(sys, z6)

        k7 = f(z + h*(35/384*k1 + 0*k2 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6), t+h)

        # The two updates
        z1_4 = z + h*k7
        z1_5 = z + h*(5179/57600*k1 + 0*k2 + 7571/16695*k3 + 393/640*k4 - 92097/339200*k5 + 187/2100*k6 + 1/40*k7)

        #used to find true midpoint for quad guard check. 
        #This prevents solver from spinning up a whole new adaptive loop so x_mid is truly in the mid which is what we want for quad matrix. 
        if !adaptive
            return z1_4, h, h
        end

        LTE = norm(z1_4 - z1_5)
        # Reject or accept?
        h_new = updated_step(LTE, tol, h, 5)

        h_end = guard(sys, z1_4)

        stage_crossed = (h_now * h2 < 0) || (h_now * h3 < 0) || (h_now * h4 < 0) || (h_now * h5 < 0) || (h_now * h6 < 0)
        end_missed = (h_now * h_end > 0)

        if stage_crossed && end_missed
            h = h / 2.0 #force smaller step
            continue
        end

        if LTE < tol
            return z1_4, h, h_new
        else
            h = h_new
        end
        if h < 1e-12
            @warn "Step size has decreased below 1e-12"
            return z1_4, h, h_new
        end
    end
end

#Richardson Extrapolation based on modifed midpoint from Wiki
function richardson_step(f::Function, z::AbstractArray, h::AbstractFloat, t::AbstractFloat)
    #Take one full step size h
    z1 = modified_midpoint_step(f,z,h,t)

    #Take two smaller steps of h/2
    h_half = h / 2
    z_half = modified_midpoint_step(f,z,h_half,t)
    z2 = modified_midpoint_step(f, z_half,h_half,t+h_half)

    #Extrapolate to get rid of lower order error
    return (4 .* z2 .- z1) ./ 3.0
end

#Extra solvers
#MagnusLeapfrog step
function magnus_leapfrog_step(f::Function, U::AbstractMatrix, h::AbstractFloat, t::AbstractFloat)
    #Extract state and fund matrix
    xₖ = U[:, 1]
    Φₖ = U[:, 2:end]

    #Midpoint Eval for base trajectory 
    x_mid = xₖ .+ (h / 2.0) .* f(xₖ, t)
    t_mid = t + h / 2.0

    #Advance base state fully
    x_next = xₖ .+ h .* f(x_mid, t_mid)

    #Compute Jacobian exactly at midpoint
    A_mid = ForwardDiff.jacobian(y -> f(y, t_mid), x_mid)

    #Lie-Group step for fund matrix 
    Φ_next = exp(h * A_mid) * Φₖ

    return hcat(x_next, Φ_next)
end
#====================================#
#TAKE STEP SOLVERS CONSOLIDATION

#Single take_step for RK methods
#Fixed step methods 
#Helper function to take the step using multiple dispatch.
compute_step(::ForwardEuler, f, x, Δt, t) = forward_euler_step(f, x, Δt, t)
compute_step(::ModifiedTrap, f, x, Δt, t) = modified_trap_step(f, x, Δt, t)
compute_step(::ModifiedMidpoint, f, x, Δt, t) = modified_midpoint_step(f, x, Δt, t)
compute_step(::RichardsonExtrapolation, f, x, Δt, t) = richardson_step(f, x, Δt, t)
compute_step(::RK4, f, x, Δt, t) = rk_4_step(f, x, Δt, t)
compute_step(::ImplicitEuler, f, x, Δt, t) = implicit_euler_step(f, x, Δt, t)

#Adaptive step methods
#Helper function to take the step via multiple dispatch
compute_step(::RK23, f, x, Δt, t, tf, sys, tol; adaptive=true) = rk_23_step(f, x, Δt, t, tf, sys, tol; adaptive=adaptive)
compute_step(::RK45, f, x, Δt, t, tf, sys, tol; adaptive=true) = rk_45_step(f, x, Δt, t, tf, sys, tol; adaptive=adaptive)

#Extra Solvers
compute_step(::MagnusLeapfrog, f, U::AbstractMatrix, Δt, t) = magnus_leapfrog_step(f, U, Δt, t)

#Note sol is not used, we do this to make using the function easier. We would need an if/else statement everytime we use this function without it
function take_step(solver::FixedRK, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver=ModifiedMidpoint(); check=true, guard_direction=default_guard_direction(prob.sys)) 
    sys = prob.sys
    x_predict = compute_step(solver, f, xₖ, Δt, tₖ)
    x_mid     = compute_step(solver, f, xₖ, Δt / 2.0, tₖ)

    if check
        #Evaluate Guards
        h_now  = guard(sys, xₖ)
        h_mid  = guard(sys, x_mid)
        h_next = guard(sys, x_predict)

        #Use cross guard check
        eventtrigger, t_root, _ = crossed_guard(sys, h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol, direction=guard_direction)

        return x_predict, eventtrigger, t_root, Δt, Δt
    else
        return x_predict, NaN, NaN, NaN, NaN
    end
end

function take_step(solver::AdaptiveRK, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver=ModifiedMidpoint(); guard_direction=default_guard_direction(prob.sys))
    sys = prob.sys
    tf = prob.tspan[2] #terminal time

    #Take adaptive step (passes tf to prevent overshooting)
    #Disable adaptive loop to get the true midpoint. Without it it would reset the loop and ruin or previous adaptivity
    x_predict, dt_used, dt_next = compute_step(solver, f, xₖ, Δt, tₖ, tf, sys, tol; adaptive=false)

    #Get midpoint for quad guard check
    #For a half-step the local ceiling is just midpoint of time 
    x_mid, _, _ = compute_step(solver, f, xₖ, dt_used / 2.0, tₖ, tf, sys, tol)

    #eval guards
    h_now  = guard(sys, xₖ)
    h_mid  = guard(sys, x_mid)
    h_next = guard(sys, x_predict)

    eventtrigger, t_root, _ = crossed_guard(sys, h_now, h_mid, h_next, tₖ, tₖ + dt_used / 2.0, tₖ + dt_used; tol=tol, direction=guard_direction)

    return x_predict, eventtrigger, t_root, dt_used, dt_next
end

#===========================#
#Helper function to tell engine how many history steps are needed for LMM
lmm_order(::AdamsBashforth2) = 2
lmm_order(::AdamsBashforth3) = 3

lmm_order(::AdaptiveABM2) = 2
lmm_order(::AdaptiveABM3) = 3

lmm_order(::BDF2) = 2

function compute_lmm_step(::AdamsBashforth2, f, xₖ, tₖ, Δt, x_history, t_history)
    #x_history[end] is x_{k-1}
    x_prev = x_history[end]
    t_prev = t_history[end]

    #Calc previous step size
    dt_previous = tₖ - t_prev

    fₖ = f(xₖ, tₖ)
    f_prev_val = f(x_prev, t_prev)

    #Variable step AB2 Formula
    α = Δt / dt_previous
    return xₖ .+ Δt .* ((1.0 + .5 * α) .* fₖ .- (.5 * α) .* f_prev_val)
end

function compute_lmm_step(::AdamsBashforth3, f, xₖ, tₖ, Δt, x_history, t_history)
    # x_history[end] is x_{k-1}, x_history[end-1] is x_{k-2}
    x_prev1 = x_history[end]
    t_prev1 = t_history[end]
    
    x_prev2 = x_history[end-1]
    t_prev2 = t_history[end-1]
    
    fₖ = f(xₖ, tₖ)
    f_prev1_val = f(x_prev1, t_prev1)
    f_prev2_val = f(x_prev2, t_prev2)
    
    # Currently using fixed-step coefficients. 
    # To upgrade to full variable-step AB3 later, you would calculate the 
    # ratios between (tₖ - t_prev1) and (t_prev1 - t_prev2) to dynamically adjust these weights.
    #This will be done eventually but for now this is fine. 
    return xₖ .+ Δt .* ( (23/12) .* fₖ .- (16/12) .* f_prev1_val .+ (5/12) .* f_prev2_val )
end

function compute_lmm_step(::BDF2, f, zₖ, tₖ, Δt, x_history, t_history)
    #retrieve state at previous time step.
    z_prev = x_history[end]
    t_new = tₖ + Δt

    #doing some algebra to isolate the implicit part of BDF2 and "c" is the right side
    c = (4.0 / 3.0) .* zₖ .- (1.0 / 3.0) .* z_prev
    #coeff scaling step size 
    α = 2.0 / 3.0

    #Gen initial guess using Exp Euler to help the Newton Root finder
    z_guess = zₖ .+ Δt .* f(zₖ, tₖ)
    #solve implicit system G(z) = 0 using Newtons
    return implicit_newton_solve(f, z_guess, c, α, Δt, t_new)
end

function compute_lmm_step(::AdaptiveABM2, f, xₖ, tₖ, Δt, x_history, t_history)
    #Extract history
    x_prev = x_history[end]
    t_prev = t_history[end]

    #Calc previous step size
    dt_prev = tₖ - t_prev
    α = Δt / dt_prev

    fₖ = f(xₖ, tₖ)
    f_prev = f(x_prev, t_prev)

    #Predictor for variable step AB2
    x_predict = xₖ .+ Δt .* ((1.0 + 0.5 * α) .* fₖ .- (0.5 * α) .* f_prev)

    #Eval vector field at pred state
    f_predict = f(x_predict, tₖ + Δt)

    #Corrector: Implicit AM2 via predicted vf
    x_correct = xₖ .+ Δt .* (0.5 .* f_predict .+ 0.5 .* fₖ)

    #Error est: Difference between corrector and predictor give local truncation error
    LTE = norm(x_correct .- x_predict)

    return x_correct, LTE
end

function compute_lmm_step(::AdaptiveABM3, f, xₖ, tₖ, Δt, x_history, t_history)
    #extract history
    #x_history[2] is x_{k-1}, x_history[1] is x_{k-2}
    x_prev1 = x_history[2]
    t_prev1 = t_history[2]

    x_prev2 = x_history[1]
    t_prev2 = t_history[1]

    #Eval vf at past known coords
    fₖ      = f(xₖ, tₖ)
    f_prev1 = f(x_prev1, t_prev1)
    f_prev2 = f(x_prev2, t_prev2)

    #Calc non uniform time grid ints
    hk = Δt                     #Current proposed step size t_{k+1} - t_k
    hk1 = tₖ - t_prev1           #Previous step size duration t_k - t_{k-1}
    hk2 = t_prev1 - t_prev2      #Two steps back duraction t_{k-1} - t_{k-2}

    #Predictor Step: Fully variable-step explicit AB3
    β₀ = hk * (hk^2 / 3.0 + hk * (2.0 * hk1 + hk2) / 2.0 + hk1 * (hk1 + hk2)) / (hk1 * (hk1 + hk2))
    β₁ = -hk^2 * (2.0 * hk + 3.0 * hk1 + 3.0 * hk2) / (6.0 * hk1 * hk2)
    β₂ = hk^2 * (2.0 * hk + 3.0 * hk1) / (6.0 * hk2 * (hk1 + hk2))

    x_predict = xₖ .+ (β₀ .* fₖ .+ β₁ .* f_prev1 .+ β₂ .* f_prev2)

    #Eval vf at predicted state
    f_predict = f(x_predict, tₖ + hk)

    #Corrector step: Adams Moulton 3
    c_β_p1 = hk * (hk / 3.0 + hk1 / 2.0) / (hk + hk1)
    c_β₀   = hk * (hk + 3.0 * hk1) / (6.0 * hk1)
    c_β_m1 = -hk^3 / (6.0 * hk1 * (hk + hk1))
    
    x_correct = xₖ .+ (c_β_p1 .* f_predict .+ c_β₀ .* fₖ .+ c_β_m1 .* f_prev1)

    #ERROR EST
    LTE = norm(x_correct .- x_predict)

    return x_correct, LTE
end

function take_step(solver::FixedLMM, prob::AbstractHybridProblem, f, xₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver = RK4();  guard_direction=default_guard_direction(prob.sys))
    sys = prob.sys
    k = lmm_order(solver)

    #Determine how many continuous steps we have since the last jump
    history_len = isempty(sol.jump_indices) ? length(sol.x) : (length(sol.x) - sol.jump_indices[end] + 1)

    h_now = guard(sys, xₖ)

    if history_len < k
        #Startup phase: Use single step predictor
        #using midpoint here to feed the quadratic guard check
        x_predict = compute_step(stepper, f, xₖ, Δt, tₖ)
        x_mid     = compute_step(stepper, f, xₖ, Δt / 2.0, tₖ)

        h_mid  = guard(sys, x_mid)
        h_next = guard(sys, x_predict)

        eventtrigger, t_root, _ = crossed_guard(sys, h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol, direction=guard_direction)
        if eventtrigger
            if (t_root - tₖ) < (1e-4 * Δt) # Add a small buffer
                eventtrigger = false
                t_root = tₖ + Δt # Reset t_root to end of step
            end
        end
        return x_predict, eventtrigger, t_root, Δt, Δt
    else
        #Multistep phase: We do have rich enough history. Extract past states
        #If sol.x[end] is xₖ, then sol.x[end-1] is x_{k-1}
        x_history = sol.x[end - k + 1 : end - 1]
        t_history = sol.t[end - k + 1 : end - 1]

        #pass history arrays forward
        x_predict = compute_lmm_step(solver, f, xₖ, tₖ, Δt, x_history, t_history)
        
        h_next = guard(sys, x_predict)

        #Clark Fix: No calc of midpoint. 
        #Look back to the prev step evaluation to build quad guard
        t_prev = sol.t[end-1]
        x_prev = sol.x[end-1]
        h_prev = guard(sys, x_prev)

        eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_next, t_prev, tₖ, tₖ + Δt; tol = tol, direction=guard_direction)
        if eventtrigger
            if (t_root - tₖ) < (1e-4 * Δt) # Add a small buffer
                eventtrigger = false
                t_root = tₖ + Δt # Reset t_root to end of step
            end
        end

        return x_predict, eventtrigger, t_root, Δt, Δt
    end
    
end


"""
    take_step(solver::AdaptiveLMM, prob::AbstractHybridProblem, f, xₖ, tₖ, h, tol, sol, stepper::RK=RK4())

Executes single variable-step Linear Multstep Method update for hybrid systems
How adaptive LMMs work:
Adaptivity requires us to know how much error we haeve when we make the current step so we can shrink or grow the step size h. LMMs do this using a Predictor/Corrector system. 
Predictor: Uses histoy points to project a rough guess for next state
Corrector: Takes that guess from predictor and refines it using the current dynamics, acting as a stabilizer.
Because predictor and corrector have known, and linked, error bounds, the difference between their outputs gives us an estimate of the Local Truncation Error (LTE).

If LTE is too high h is shrunk, and history is interpolated or reset, and we try to step again. If LTE is below our tolerance the step is accepted and solver calcs a slightly larger h for next step. 



How it works:
1) History Validation: The solver checks the length of continuous history since the last disc jump. If the history buffer is smaller than the order of the LMM 'k', 
the step is delegated to the single-step stepper (default RK4).

2)Predictor/Corrector and Error Est: Once history is established, the solver extracts the past `k` states and calls `compute_lmm_step`. This helper function executes the explicit prediction, the implicit correction, and extracts the isolated LTE using Milne's device.

3) Adaptivity Loop: Operates within 'while true' rejection loop. If LTE exceeds 'tol' the step size 'h' shrinks, and the step is recalculated. If accepted, it computes the optimal h_next for the subsequent step

WHY I DID THINGS:
We use Milne's device for error est because it is easy to compute. Since we already are performing an explicit prediction and an implicit correction, the difference serves as a solid estimate.

"""
#user can specify if they want RK4 here
function take_step(solver::AdaptiveLMM, prob::AbstractHybridProblem, f, xₖ, tₖ, h, tol, sol, stepper::RK=RK4(); guard_direction=default_guard_direction(prob.sys))
    sys = prob.sys
    k = lmm_order(solver)
    tf = prob.tspan[2]

    h = minimum([h, tf - tₖ])

    #determine how many cont steps we have since last jump
    history_len = isempty(sol.jump_indices) ? length(sol.x) : (length(sol.x) - sol.jump_indices[end] + 1)

    #Startup phase: IF we dont have history we use RK stepper
    if history_len < k
        return take_step(stepper, prob, f, xₖ, tₖ, h, tol, sol, stepper; guard_direction=guard_direction)
    end

    #Extract history for LMM
    x_history = sol.x[end - k + 1 : end - 1]
    t_history = sol.t[end - k + 1 : end - 1]

    h_now = guard(sys, xₖ)

    #Adaptive step loop
    #Wanted to try using a max_iter variation here but that breaks things (for some reason idk)
    while true
        #compute step and retrieve LTE 
        x_next, LTE = compute_lmm_step(solver, f, xₖ, tₖ, h, x_history, t_history)

        #Calc proposed next step size using helper from beginning 
        h_new = updated_step(LTE, tol, h, k)

        if LTE < tol
            #step accepted
            h_next = guard(sys, x_next)

            #Guard eval looking back to prev step for the quad check 
            idx = max(1, length(sol.x) - 1)
            t_prev = sol.t[idx]
            h_prev = guard(sys, sol.x[idx])

            eventtrigger, t_root, _ = crossed_guard(sys, h_prev, h_now, h_next, t_prev, tₖ, tₖ + h; tol=tol, direction=guard_direction)

            return x_next, eventtrigger, t_root, h, h_new
        else
            #Step rejected: shrink and try again 
            h = h_new
            if h < 1e-6
                @warn "LMM Step size has decreased below 1e-6"
                #Force break to avoid inf loops
                return x_next, false, NaN, h, h_new
            end
        end
    end
end

#Extra Solvers that dont fit into RK or LMM
#Exponential Goes here

#Supposedly useful as when doing variational equation stuff the magnus expansion preserves Lie group structure (I dont know enough about Lie groups to tell you what that means)
#so we maintain the properties we want like volume and determinants which is very good for variational equations and eventual Lyapunov Exponents. 
function take_step(solver::MagnusLeapfrog, prob::AbstractHybridProblem, f, Uₖ, tₖ, Δt, tol, sol, stepper::AbstractODESolver=ModifiedMidpoint(), guard_direction=default_guard_direction(prob.sys))
    sys = prob.sys

    #Compute augmented matrix predictions
    U_predict = compute_step(solver, f, Uₖ, Δt, tₖ)
    U_mid     = compute_step(solver, f, Uₖ, Δt / 2.0, tₖ)

    #Guard checks only care about phys state
    xₖ = Uₖ[:, 1]
    x_mid = U_mid[:, 1]
    x_predict = U_predict[:, 1]

    #Eval guards
    h_now = guard(sys, xₖ)
    h_mid = guard(sys, x_mid)
    h_next = guard(sys, x_predict)

    #Use check
    eventtrigger, t_root, _ = crossed_guard(sys, h_now, h_mid, h_next, tₖ, tₖ + Δt / 2.0, tₖ + Δt; tol=tol, direction=guard_direction)

    return U_predict, eventtrigger, t_root, Δt, Δt
end

#Newton Raphson solver to find roots for implicit integration steps. 
#Solves for z in G(z) = z - c - α*h*f(z, t_new) = 0
function implicit_newton_solve(f::Function, z_guess::Vector, c::Vector, α::AbstractFloat, h::AbstractFloat, t_new::AbstractFloat; max_iter=50, tol=1e-6)
    z_curr = copy(z_guess)

    for _ in 1:max_iter
        #Eval vector field at curent guess for z_{n+1}
        val = f(z_curr, t_new)
        #construct residual function G(z) = z - c - αhf(z). When G(z) = 0 we have a sol
        G = z_curr .- c .- (α*h) .* val

        #Convergence check
        if norm(G) < tol
            return z_curr
        end

        #Calc Jacobian of vector field at our current guess
        J_f = ForwardDiff.jacobian(y -> f(y, t_new), z_curr)
        #Calc Jacobian of the residual G(z). 
        J_G = I - (α*h) * J_f

        #Update based on Newtons method. 
        z_curr = z_curr .- J_G \ G
    end
    #Just in case
    @warn "Newton Raphson solver failed to converge at t = $t_new. Consider a smaller step size."
    #return NaN to signal the sim that this step is invalid. 
    return fill(NaN, size(z_curr))
end

#====================================#
#Event detection utility. 
#Assigns a default crossing direction to a specific system to reduce issues (this may be made better for the front end later but for now this works)
default_guard_direction(sys::MechanicalSystem) = :falling
default_guard_direction(sys::NonholonomicSystem) = :falling
#Gen system lack specific constraints so we monitor crossings in both directions
default_guard_direction(sys::GeneralSystem) = :both
default_guard_direction(sys) = :both

#Wrapper function to interface between the system state and core logic
function crossed_guard(sys, h_prev, h_now, h_next, t_prev, t_now, t_next; 
                       tol=1e-6, direction=default_guard_direction(sys))   
    # Calls evaluator with the directionality
    return evaluate_crossing(h_prev, h_now, h_next, t_prev, t_now, t_next, direction; tol=tol)
end

#Core engine: determines if/when the guard function 'h' changes sign. 
function evaluate_crossing(h_prev, h_now, h_next, t_prev, t_now, t_next, direction::Symbol; tol=1e-6)
    #Helper to validate a linear sign change based on the required direction
    valid_linear(h1, h2) = 
        (direction == :both && h1 * h2 < 0) ||
        (direction == :falling && h1 > 0 && h2 < 0) ||
        (direction == :rising && h1 < 0 && h2 > 0)

    #First check: Simple linear crossing detection. Between previous and current
    if valid_linear(h_prev, h_now)
        #Linear interpolation to find the root
        t_root = t_prev - h_prev * (t_now - t_prev) / (h_now - h_prev)
        return true, t_root, NaN
    #Second check: Linear beween current and next
    elseif valid_linear(h_now, h_next)
        t_root = t_now - h_now * (t_next - t_now) / (h_next - h_now)
        return true, t_root, NaN 
    end

    #Quad version
    try
        #Set up linear system to solve for parabola coeffs  
        A = [t_prev^2 t_prev 1; t_now^2 t_now 1; t_next^2 t_next 1]
        a, b, c  = A \ [h_prev, h_now, h_next]

        #Ensure it is actually a parabola then check disc
        if abs(a) > 1e-6
            discriminant = b^2 - 4*a*c
            if discriminant > 0
                sqrt_d = sqrt(discriminant)
                r1 = (-b - sqrt_d) / (2*a)
                r2 = (-b + sqrt_d) / (2*a)

                valid_roots = Float64[]
                eps_t = 1e-9 #Buffer to ensure root isnt some weird artifact (it does seem necessary)

                #Derivative of parabola eval the slope of root. 
                h_prime(t) = 2*a*t + b

                #Helper: Does quad curve cross in the correct direction?
                valid_quad(r) = 
                    (direction == :both) ||
                    (direction == :falling && h_prime(r) < 0) ||
                    (direction == :rising && h_prime(r) > 0)

                #Check bounds and direction
                if (t_prev + eps_t) <= r1 <= t_next && valid_quad(r1) push!(valid_roots, r1) end 
                if (t_prev + eps_t) <= r2 <= t_next && valid_quad(r2) push!(valid_roots, r2) end

                #If root exists, return earliest and the parabolas critical point
                if !isempty(valid_roots)
                    proposed_root = minimum(valid_roots)
                    critical_point = -b / (2*a)
                    return true, proposed_root, critical_point
                end
            end
        end
    catch
        #If linear system is singular or math fails, return failure to cross. 
        return false, NaN, NaN
    end
    return false, NaN, NaN
end

#Locator Dispatches
#Isolates the root finding mathematics inside each one. This gets rid of global helpers so when we add new locators its really easy

#-----------------------
#LOCATOR TAGS
#These are similar to the solver tags. But these define how the solver finds the exact crossing time once an event is detected. 

#Tag to use Linear Interpolation. Very fast but can be innacurate for higher order methods. 
struct LinearLocator <: AbstractEventLocator end

#Tag to use a bisection method serach. Can be very accurate but also very slow with complex systems
struct BisectionLocator <: AbstractEventLocator end

#Tag for quadratic event locator
struct QuadraticLocator <: AbstractEventLocator end

#Tag to use Newtons method for event locators
struct NewtonLocator <: AbstractEventLocator end

#Bisection Method (Iterative)

function locate_event(::BisectionLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys
    τ_l, τ_r = 0.0, Δt
    h_l = h_now
    
    for _ in 1:100 # max_iter
        if (τ_r - τ_l) < tol break end

        #test midpoint
        τ_m = (τ_l + τ_r) / 2.0

        x_m, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_m, tol, sol, stepper)
        h_m = guard(sys, x_m)
    
        if signbit(h_l) != signbit(h_m)
            τ_r = τ_m
        else
            τ_l = τ_m
            h_l = h_m
        end
    end
    
    t_star = tₖ + τ_l
    x_star, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_l, tol, sol, stepper)
    return t_star, x_star
end
#Linear Interpolation

function locate_event(::LinearLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys
    τ_l, τ_r = 0.0, Δt
    h_l = h_now

    if abs(h_now) < tol || h_now < 0
        return tₖ, xₖ
    end

    #Get right side of boundary
    x_r, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_r, tol, sol, stepper)
    h_r = guard(sys, x_r)

    τ_star = Δt
    x_star = x_r

    for _ in 1:100
        if abs(τ_r - τ_l) < tol || abs(h_r) < tol
            τ_star = τ_r
            x_star = x_r
            break
        end

        #Linear Interp
        τ_m = τ_r - h_r * (τ_r - τ_l) / (h_r - h_l)

        x_m, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_m, tol, sol, stepper)
        h_m = guard(sys, x_m)

        #keep root bracketed
        if signbit(h_l) != signbit(h_m)
            τ_r = τ_m
            h_r = h_m
        else 
            τ_l = τ_m
            h_l = h_m
        end
    end
    t_star = tₖ + τ_star
    return t_star, x_star
end

function locate_event(::QuadraticLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys

    #Get three points 
    h₀ = h_now

    #middle point. We just take a half step instead of going one before the start point. I think itll be more stable. 
    x₁, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, Δt / 2.0, tol, sol, stepper)
    h₁ = guard(sys, x₁)

    #endpoint
    x₂, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, Δt, tol, sol, stepper)
    h₂ = guard(sys, x₂)

    #compute parabola coeffs
    c = h₀
    b = (-3.0 * h₀ + 4.0 * h₁ - h₂) / Δt
    a = 2.0 * (h₀ - 2.0 * h₁ + h₂) / (Δt ^ 2)

    #root finding with Fallbacks
    if abs(a) < tol
        #curve is basically zero, fallback to linear interp
        θ = -h₀ / (h₂ - h₀)
        τ_star = θ * Δt
    else 
        discriminant = b^2 - 4.0 * a * c

        if discriminant < 0
            #fallback to linear interp
            θ = -h₀ / (h₂ - h₀)
            τ_star = θ * Δt
        else
            #stable quad root finding
            q = -.5 * (b + sign(b) * sqrt(discriminant))
            root_1 = q / a
            root_2 = c / q

            valid_1 = 0.0 <= root_1 < Δt
            valid_2 = 0.0 <= root_2 < Δt

            #select right root
            if valid_1 && valid_2
                #if parabola crosses twice we pick first one
                τ_star = min(root_1, root_2)
            elseif valid_1
                τ_star = root_1
            elseif valid_2
                τ_star = root_2
            else
                #roots drifted to narnia? 
                θ = -h₀ / (h₂ - h₀)
                τ_star = θ * Δt
            end
        end
    end
    t_star = tₖ + τ_star
    x_star, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_star, tol, sol, stepper)

    return t_star, x_star
end

function locate_event(::NewtonLocator, prob, solver::AbstractODESolver, f, xₖ, tₖ, Δt, h_now, tol, sol, stepper::RK = RK4())
    sys = prob.sys

    τ_prev = 0.0
    h_prev = h_now

    τ_curr = Δt
    x_curr, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_curr, tol, sol, stepper)
    h_curr = guard(sys, x_curr)

    for _ in 1:100
        if abs(h_curr) < tol || abs(τ_curr - τ_prev) < tol
            break
        end

        #Estimate derivative 
        dh_dτ = (h_curr - h_prev) / (τ_curr - τ_prev)

        if abs(dh_dτ) < tol
            @warn "Derivative less than tolerance: Newton method unstable. Terminating."
            break
        end

        #Newton step
        τ_next = τ_curr - h_curr / dh_dτ

        #Bound next step to prevent over shooting 
        τ_next = clamp(τ_next, 0.0, Δt)

        #update for next iteration
        τ_prev = τ_curr
        h_prev = h_curr

        τ_curr = τ_next
        x_curr, _, _, _, _ = take_step(solver, prob, f, xₖ, tₖ, τ_curr, tol, sol, stepper)
        h_curr = guard(sys, x_curr)
    end
    t_star = tₖ + τ_curr
    return t_star, x_curr
end


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
