## A General Example

!!! info "A General Problem"
	A general system consists of the data $(f, h, \Delta)$ where
	```math
		\begin{cases}
			\dot{x} = f(x, t), & h(x) \ne 0, \\
			x^+ = \Delta(x), & h(x) = 0.
		\end{cases}
	```

```math
	\begin{cases}
		\ddot{x} = -g -\alpha\cdot \dot{x} \cdot|\dot{x}|, & x > 0 \\[1ex]
		\dot{x}^+ = -e\dot{x}^-, & x = 0
	\end{cases}
```

```@example gen
import HybridDynamics as HD
import Plots as plt
using LaTeXStrings
```

```@example gen
function f_ball(x, t)
	g = 9.81
	α = 0.1
	q, v = x
	return [v, -g-α*v*abs(v)]
end

h_ball(x) = x[1]
Δ_ball(x) = [abs(x[1]), -0.8*x[2]]

sysG = HD.GeneralSystem(f_ball, h_ball, Δ_ball; direction=-1)
probG = HD.prob(sysG, [10.0, 0.0], (0.0, 15.0))
solG = HD.solve(probG)
```
```@example gen
times, states = HD.split_jumps(solG)
plt.plot(times, getindex.(states, 1), label="Position", lw=2, lc=:blue)
plt.plot!(times, getindex.(states, 2), label="Velocity", lw=2, lc=:orange)
plt.plot!(title = "Bouncing ball", xlabel = "Time", grid = true)
```
