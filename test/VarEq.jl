using HybridDynamics
using Plots
using LinearAlgebra
using ForwardDiff

function f_ball_var(x, t)
    g = 9.81
    velocity = x[2]
    return [velocity, -g - velocity * abs(velocity)]
end

function h_floor_var(x)
    return x[1] 
end

function Δ_bounce_var(x, t=0.0)
    restitution = 0.8
    new_position = abs(x[1]) 
    new_velocity = restitution * abs(x[2])
    return [new_position, new_velocity]
end

n = 2

f_aug = (u, t) -> variational_vector_field(f_ball_var, u, t, n)

h_aug = (u) -> h_floor_var(u[1:n])

Δ_aug = (u) -> begin
    u_copy = copy(u) # Ensure we do not overwrite the current step state prematurely
    return apply_variational_jump(u_copy, n, f_ball_var, Δ_bounce_var, h_floor_var, 0.0)
end

sys_aug = GeneralSystem(f_aug, h_aug, Δ_aug)

x0 = [10.0, 0.0]                      # Initial state (position, velocity)
Φ0 = Matrix{Float64}(I, n, n)         # Initial fundamental solution matrix (Identity)
u0 = vcat(x0, vec(Φ0))                # Augmented state vector [x1, x2, Φ11, Φ21, Φ12, Φ22]

tspan = (0.0, 5.0)   
problem_aug = prob(sys_aug, u0, tspan)

sol = solve(problem_aug)

function get_rid_of_jump_lines_var(sol)
    t = Float64[]
    p11, p21, p12, p22 = Float64[], Float64[], Float64[], Float64[]

    for i in 1:length(sol.x)
        push!(t, sol.t[i])
        push!(p11, sol.x[i][3])
        push!(p21, sol.x[i][4])
        push!(p12, sol.x[i][5])
        push!(p22, sol.x[i][6])

        if i < length(sol.t) && sol.t[i] == sol.t[i+1]
            push!(t, NaN)
            push!(p11, NaN); push!(p21, NaN); push!(p12, NaN); push!(p22, NaN)
        end
    end
    return t, p11, p21, p12, p22
end

t_plot, p11, p21, p12, p22 = get_rid_of_jump_lines_var(sol)

# Plot the sensitivity values over time
p_display = plot(t_plot, [p11 p21 p12 p22], 
                 label=["Φ11 (dx/dx0)" "Φ21 (dv/dx0)" "Φ12 (dx/dv0)" "Φ22 (dv/dv0)"], 
                 linewidth=2.0, 
                 title="Variational Equation Testing", 
                 xlabel="Time", 
                 ylabel="Sensitivity Value",
                 grid=true)

display(p_display)