## Linear and Affine Systems
!!! info "Linear and Affine Problems"
	A Linear Hybrid System consists of the data $(A, λ, C)$ where 
	```math
		\mathcal{LH} = 
		\begin{cases}
			\dot x = Ax, & λx \neq 0, \\
			x^+ = Cx, & λx = 0.
		\end{cases}
	```
	An Affine Hybrid System consists of the data $(A, b, λ, a, C, κ)$ where
	```math
		\mathcal{AH} = 
		\begin{cases}
			\dot x = Ax + b, λx + a \neq 0, \\
			x^+ = Cx + κ, λx + a = 0.
		\end{cases}
	```

	1. .$A$ is the state matrix. 
	2. .$b$ is the continuous affine vector ($\dot x = Ax + b$).
	3. .$λ$ is the guard normal vector.
	4. .$a$ is the guard affine vector ($λx + a = 0$).
	5. .$C$ is the reset matrix
	6. .$κ$ is the reset affine vector $(x^+ = Cx + κ)$.


```@example linaf
import HybridDynamics as HD
import Plots as plt
using LaTeXStrings

LA_type = "Linear" # Options: "Linear", "Affine"
LA_A = [-0.5 -10.0; 10.0 -0.5]
LA_λ = [1.0, 0.0]
LA_C = [0.8 0.0; 0.0 0.8]
LA_x0 = [1.0, 0.5]
LA_tspan = (0.0, 2.0)

if LA_type == "Affine"
	LA_b = [2.0, 0.0]
	LA_a = 0.2
	LA_κ = [0.1, -0.2]
	LA_sys = HD.AffineSystem(LA_A, LA_b, LA_λ, LA_a, LA_C, LA_κ)
	LA_guard = -LA_a
	LA_prob = HD.prob(LA_sys, LA_x0, LA_tspan)
	LA_sol = HD.solve(LA_prob, HD.RK45())
else
	LA_sys = HD.LinearSystem(LA_A, LA_λ, LA_C)
	LA_guard = 0.0
	LA_prob = HD.prob(LA_sys, LA_x0, LA_tspan)
	LA_sol = HD.solve(LA_prob, HD.RK4())
end
```

```@example linaf
LA_t_list, LA_x_list = HD.split_jumps(LA_sol)

p = plt.plot(title="$(LA_type) System Trajectory", xlabel="x1", ylabel="x2", 
				aspect_ratio=:equal, grid=true, legend=false)

plt.vline!(p, [LA_guard], label="Guard", color=:black, alpha=0.3, lw=2)

plt.plot!(p, getindex.(LA_x_list, 1), getindex.(LA_x_list, 2), 
			color=:red, linestyle=:dash, lw=1.5, label="Trajectory")
p
```