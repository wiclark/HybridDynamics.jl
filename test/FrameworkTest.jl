using HybridDynamics
using Plots

use_affine = true #CHANGE THIS FOR LINEAR VS AFFINE TEST

A = [0.0 -10.0; 10.0 -.1]
λ = [1.0, 0.0]
C = [1.5 1.0; 0.0 0.5]
x0 = [1.0, 0.5]
tspan = (0.0, 10.0)

function get_rid_of_jump_lines(sol)
    x1 = Float64[]
    x2 = Float64[]

    for i in 1:length(sol.x)
        push!(x1, sol.x[i][1])
        push!(x2, sol.x[i][2])
        if i < length(sol.t) && sol.t[i] == sol.t[i+1]
            push!(x1, NaN)
            push!(x2, NaN)
        end
    end
    return x1, x2
end

if use_affine
    b = [2., 0.0]
    a = .2
    κ = [.1, -.2]

    sys = AffineSystem(A, b, λ, a, C, κ)
    problem = prob(sys, x0, tspan)
    plot_title = "Affine System Test"

    my_solver = AdaptiveABM3()
    event_method = LinearLocator()
    #my_stepper = ForwardEuler()

    sol_trap = solve(problem, my_solver; event_method)

    #println("You just used $(typeof(my_solver)) and $(typeof(my_stepper)) for this run.")
    println("You just used $(typeof(my_solver))")

    p = plot(title=plot_title, xlabel="x1", ylabel="x2", grid=true, aspect_ratio=:equal)

    vline!(p, [-a], label="Guard", color=:black, alpha=.3)

    x1_trap, x2_trap = get_rid_of_jump_lines(sol_trap)
    plot!(p, x1_trap, x2_trap, label="Modified Trapezoid Method", color=:red, linestyle=:dash)

    display(p)
else 
    sys = LinearSystem(A, λ, C)
    problem = prob(sys, x0, tspan)
    plot_title = "Linear and Exact Solution Test"
    
    # Both solvers run to compare the purely linear system
    sol_trap = solve(problem, AdamsBashforth3())
    sol_exp  = solve(problem, ExponentialSolver())

    p = plot(title=plot_title, xlabel="x1", ylabel="x2", aspect_ratio=:equal)
    
    # Linear guard is always at the origin for x1
    vline!(p, [0.0], label="Guard", color=:black, alpha=0.3)
    
    x1_trap, x2_trap = get_rid_of_jump_lines(sol_trap)
    plot!(p, x1_trap, x2_trap, label="ModifiedTrap", color=:red, linestyle=:dash)
    
    x1_exp, x2_exp = get_rid_of_jump_lines(sol_exp)
    plot!(p, x1_exp, x2_exp, label="Exponential", color=:blue, alpha=0.7)

    display(p)
end
