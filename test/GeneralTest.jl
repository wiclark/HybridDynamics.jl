
using HybridDynamics
using Plots

#Bouncing Ball example
function f_ball(x, t)
    g = 9.81
    drag_coeff = 0.1
    velocity = x[2]
    
    return [velocity, -g - drag_coeff * velocity * abs(velocity)]
end

# Event Surface: The Floor
function h_floor(x)
    return x[1] # Triggers exactly when position hits 0
end

# Reset Map: The bounce
function Δ_bounce(x)
    restitution = 0.8
    # Keep position, reverse velocity and lose 20% 
    return [x[1], -restitution * x[2]]
end

sys = GeneralSystem(f_ball, h_floor, Δ_bounce)

x0 = [10.0, 0.0]
tspan = (0.0, 5.0)   

problem = prob(sys, x0, tspan)

sol = solve(problem, ModifiedTrap())

function get_rid_of_jump_lines(sol)
    t = Float64[]
    x1 = Float64[]
    x2 = Float64[]

    for i in 1:length(sol.x)
        push!(t, sol.t[i])
        push!(x1, sol.x[i][1])
        push!(x2, sol.x[i][2])

        # If next timestamp is identical (jump event), break the line
        if i < length(sol.t) && sol.t[i] == sol.t[i+1]
            push!(t, NaN)
            push!(x1, NaN)
            push!(x2, NaN)
        end
    end
    return t, x1, x2
end

t_plot, x1_plot, x2_plot = get_rid_of_jump_lines(sol)

p = plot(t_plot, x1_plot, label="Position (m)", linewidth=2.5, color=:blue,title="Nonlinear Bouncing Ball", xlabel="Time (s)", ylabel="State Value",grid=true)

# Add velocity to the same plot
plot!(p, t_plot, x2_plot, label="Velocity (m/s)", linewidth=1.5, color=:orange)

display(p)
