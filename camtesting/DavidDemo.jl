using Pkg
Pkg.activate(".")
using HybridDynamics
using LinearAlgebra
using Plots

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


#Setup System
A = [0.0 -1.0; 1.0 0.0] 
λ = [1.0, 0.0]           
C = [0.5 1.0; 0.0 0.5]   
sys = HybridLinearSystem(A, λ, C)

x₀ = [1.0, 0.5]
tspan = (0.0, 10.0)
dt = 0.05
prob = HybridLinearProblem(sys, x₀, tspan)

# Time the Standard Solver
t_std = @elapsed sol_std = solve_hybrid_system(prob, dt; step_method=modified_midpoint_step, is_adaptive=false,max_iter=10^6, tol = 1e-12) #euler wierd?
println("Standard Solver: $(round(t_std, digits=6)) seconds")

# Time the Exponential Solver
t_exp = @elapsed sol_exp = solve_hybrid_system_exp(prob, dt)
println("Exponential Solver: $(round(t_exp, digits=6)) seconds")

#Plot Trajectory Comparison
p1 = plot(title="Solver Comparison", xlabel="x1", ylabel="x2", aspect_ratio=:equal)
vline!(p1, [0.0], color=:black, alpha=0.3, label="Guard")

#Plot Standard
x1_std, x2_std = get_rid_of_jump_lines(sol_std)
plot!(p1, x1_std, x2_std, label="ODE Solver", color=:red, linestyle=:dash)

#Plot Exponential
x1_exp, x2_exp = get_rid_of_jump_lines(sol_exp)
plot!(p1, x1_exp, x2_exp, label="Exponential", color=:blue, alpha=0.7)

#AFFINE EXAMPLE THAT I DIDNT HAVE TIME TO TEST BUT SHOULD WORK?
#b = [.2, 0.0]  #flow drift
#a = -.5        #Guard const
#κ = [-.1, .1]  #Reset const
#sys_aff = HybridAffineSystem(A,b,λ,a,C,κ)
#prob_aff = HybridAffineProblem(sys_aff, x₀, tspan)
#t_aff = @elapsed sol_aff = solve_hybrid_system(prob_aff, dt; step_method=modified_midpoint_step, is_adaptive=false, max_iter=10^6, tol = 1e-12)
#println("Affine Solver: $(round(t_aff, digits=6)) seconds")

#p2 = plot(title="Affine Trajectory", xlabel="x1", ylabel="x2", aspect_ratio=:equal)
#vline!(p2, [-a / λ[1]], color=:black, alpha=0.3, label="Affine Guard")
#x1_aff = [pt[1] for pt in sol_aff.x]
#x2_aff = [pt[2] for pt in sol_aff.x]
#plot!(p2, x1_aff, x2_aff, label="Affine Solver", color=:green)
