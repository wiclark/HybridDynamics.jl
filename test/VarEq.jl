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

f_aug = (U, t) -> variational_vector_field(f_ball_var, U, t)

h_aug = (U) -> h_floor_var(U[:, 1])

Δ_aug = (U) -> begin
    U_copy = copy(U)
    return apply_variational_jump(U_copy, f_ball_var, Δ_bounce_var, h_floor_var, 0.0)
end

sys_aug = GeneralSystem(f_aug, h_aug, Δ_aug)

n = 2
x0 = [10.0, 0.0]                      
Φ0 = Matrix{Float64}(I, n, n)         
U0 = hcat(x0, Φ0)                

tspan = (0.0, 15.0)   
problem_aug = prob(sys_aug, U0, tspan)

sol = solve(problem_aug, AdamsBashforth3(); 
            event_method=LinearLocator(), 
            dt_initial=0.01, 
            dt_min=1e-6,
            tol=1e-6)

t_plot, x_plot = split_jumps(sol)

p11 = [state[1, 2] for state in x_plot]
p21 = [state[2, 2] for state in x_plot]
p12 = [state[1, 3] for state in x_plot]
p22 = [state[2, 3] for state in x_plot]

p_display = plot(t_plot, [p11 p21 p12 p22], 
                 label=["Φ11 (dx/dx0)" "Φ21 (dv/dx0)" "Φ12 (dx/dv0)" "Φ22 (dv/dv0)"], 
                 linewidth=2.0, 
                 title="Variational Equation Testing", 
                 xlabel="Time", 
                 ylabel="Sensitivity",
                 grid=true)

display(p_display)