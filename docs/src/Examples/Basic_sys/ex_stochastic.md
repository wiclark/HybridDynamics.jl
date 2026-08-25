## A Stochastic Example

!!! info "Stochastic Systems"
	A stochastic system contains the data $(f, g, h, Δ; δ)$
	```math
		\begin{cases}
			dx = f(x,t)dt + g(x,t)dW, & h(x) \ne 0, \\
			x^+ = \Delta(x), & h(x) = 0
		\end{cases}
	```
	Notice that $g$ must be given as a matrix.

```math
	\begin{cases}
		dx = -(x-3)dt + 0.2 dW, & x < 2 \\
		x^+ = x-1, & x = 2
	\end{cases}
```

```@example stoch
import HybridDynamics as HD
import Plots as plt
using LaTeXStrings

f_st(x, t) = -(x .-3.)
g_st(x, t) = [0.2;;]
h_st(x) = x[1]-2
Δ_st(x) = x .- 1.0

sysST = HD.StochasticSystem(f_st, g_st, h_st, Δ_st)
probST = HD.prob(sysST, [0.5], (0.0, 3.0))
solST = HD.solve(probST)
```

```@example stoch
ts, ss = HD.split_jumps(solST)
plt.plot(ts, getindex.(ss, 1), label="", lw=2, lc=:blue)
plt.plot!(title = "Stochastic example", xlabel = "Time", grid = true)
```