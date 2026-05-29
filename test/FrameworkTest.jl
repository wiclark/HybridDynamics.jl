using Pkg
Pkg.activate(".")
using HybridDynamics
using LinearAlgebra
using Plots


A = [0.0 -1.0; 1.0 1.0]
λ = [1.0, 0.0]
C = [0.5 1.0; 0.0 0.5]
sys = LinearSystem(A, λ, C)
prob = LinearProblem(sys, [1.0, 0.5], (0.0, 10.0))


sol_trap = solve(prob, ModifiedTrap(); event_method=LinearLocator())
sol_exp  = solve(prob, ExponentialSolver(); event_method=BisectionLocator())

p = plot(title="Solver Verification", xlabel="x1", ylabel="x2", aspect_ratio=:equal)
vline!(p, [0.0], label="Guard", color=:black, alpha=0.3)

# Helper function to break lines at jump points
function get_rid_of_jump_lines(sol)
    x1 = Float64[]
    x2 = Float64[]

    for i in 1:length(sol.x)
        push!(x1, sol.x[i][1])
        push!(x2, sol.x[i][2])

        # If next timestamp is identical (jump event), break the line
        if i < length(sol.t) && sol.t[i] == sol.t[i+1]
            push!(x1, NaN)
            push!(x2, NaN)
        end
    end
    return x1, x2
end

# Plot ModifiedTrap
x1_trap, x2_trap = get_rid_of_jump_lines(sol_trap)
plot!(p, x1_trap, x2_trap, label="ModifiedTrap", color=:red, linestyle=:dash)

# Plot Exponential
x1_exp, x2_exp = get_rid_of_jump_lines(sol_exp)
plot!(p, x1_exp, x2_exp, label="Exponential", color=:blue, alpha=0.7)

display(p)