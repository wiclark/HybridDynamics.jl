import HybridDynamics as HD
import Plots as plt

f(x, t) = -0.5 .* x #Decay

Δ(x) = x .+ 1.0 #boosts by one

event_times = [2.0, 4.5, 8.0, 2.5] #jump times 

x0 = [1.0]
tspan = (0.0, 10.0)

sys = HD.GeneralSystem(f, event_times, Δ)
prob = HD.prob(sys, x0, tspan)
sol = HD.solve(prob)

for i in sol.event_indices
    t_jump = sol.t[i]

    x_before = sol.x[i-1][1]
    x_after = sol.x[i][1]
    diff = x_after - x_before

    # just for testing to see if the jump worked at the right time. Should always be difference of 1 with 'event_times' in correct order. 
    println("At t = $t_jump | x_before = $(round(x_before, digits=4)), x_after = $(round(x_after, digits=4)) | Difference = $(round(diff, digits=4))")
end

# its a little easier to see without split jumps tbh but its here. 
t_split, x_split = HD.split_jumps(sol)
x_vals = [x[1] for x in x_split]

p = plt.plot(t_split, x_vals, 
    label="Continuous State", 
    xlabel="Time (t)", 
    ylabel="State (x)", 
    title="Time-Triggered Jumps Demo",
    linewidth=2,
    color=:blue,
    legend=:topright)

# highlight post jump points for clarity
event_t = sol.t[sol.event_indices]
event_x = [sol.x[i][1] for i in sol.event_indices]

plt.scatter!(p, event_t, event_x, 
    label="Post-Jump State (x⁺)", 
    color=:red, 
    markersize=5)
    
plt.display(p)

