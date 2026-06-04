
struct GeneralSystem <: AbstractHybridSystem
    n::Int          #Dimension of system
    f::Function     #Continuous Dynamics: (x,t) -> dx/dt
    h::Function     #Guard Surface: x-> real
    Δ::Function     #Reset map: x-> x⁺
end

struct GeneralProblem <: AbstractHybridProblem
    sys::GeneralSystem
    x₀::Vector{Float64}
    tspan::Tuple{Float64, Float64}
end

struct GeneralSolution 
    t::Vector{Float64}
    x::Vector{Vector{Float64}}
    jump_times::Vector{Float64}
    jump_indices::Vector{Int}
end

function GeneralProblem(f::Function, h::Function, Δ::Function, x₀::Vector{Float64}, tspan)
    n = length(x₀)
    sys = GeneralSystem(n, f, h, Δ)
    return GeneralProblem(sys, x₀, tspan)
end

function GeneralProblem(f::Function, x₀::Vector{Float64}, tspan)
    n = length(x₀)
    dummy_h(x) = 1.0
    dummy_Δ(x) = x
    sys = GeneralSystem(n, f, dummy_h, dummy_Δ)
    return GeneralProblem(sys, x₀, tspan)
end

function init_solution(prob::GeneralProblem)
    return GeneralSolution([prob.tspan[1]], [prob.x₀], Float64[], Int[])
end

function guard(sys::GeneralSystem, x::AbstractVector)
    return sys.h(x)
end
function apply_reset(sys::GeneralSystem, x::AbstractVector)
    return sys.Δ(x)
end

function solve(prob::GeneralProblem, solver::AbstractODESolver=ModifiedMidpoint(); event_method::AbstractEventLocator=BisectionLocator(), dt_initial=.01, dt_min = 1e-6, max_iter = 10^6, tol = 1e-6)
    sys = prob.sys
    f = sys.f
    sol = init_solution(prob)
    t_start, t_end = prob.tspan
    Δt = dt_initial 
    iter = 0

    while sol.t[end] < t_end 
        iter += 1
        if iter > max_iter 
            @info "Maximum iterations $max_iter reached."
            break
        end

        #terminaite if time is below machine precision
        if t_end - sol.t[end] < dt_min
            @info "Time step below minimum threshold $dt_min. Terminating."
            break
        end

        #truncate time step if we overshoot
        dt_step = (sol.t[end] + Δt > t_end) ? (t_end - sol.t[end]) : Δt

        #continuous integration
        xₖ = sol.x[end]
        tₖ = sol.t[end]

        #attempt cont step
        x_predict, eventtriggered, h_now, dt_next = take_step(solver, sys, f, xₖ, tₖ, dt_step, tol, sol)

        #discrete event logic
        if eventtriggered
            #pinpoint exact time and state even happened
            t_star, x_star = locate_event(event_method, sys, solver, f, xₖ, tₖ, dt_step, h_now, tol, sol)

            #apply reset
            x⁺ = sys.Δ(x_star)

            #track data: pre and post jump states
            push!(sol.t, t_star, t_star)
            push!(sol.x, x_star, x⁺)

            #record for event analysis
            if hasproperty(sol, :jump_times)
                push!(sol.jump_times, t_star)
                push!(sol.jump_indices, length(sol.t))
            end

            #shrink main step size to avoid overshooting
            Δt = dt_min

        else
            t_next = tₖ + dt_step
            push!(sol.t, t_next)
            push!(sol.x, x_predict)
            Δt = dt_initial
        end
    end
    return sol
end