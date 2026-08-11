### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 6dab3570-7290-11f1-b67c-277aa2a67d52
begin
	using Pkg
	Pkg.activate("..")
	Pkg.resolve()
	Pkg.instantiate()
	Pkg.precompile()
end

# ╔═╡ 2b7869d2-dcd8-4ccf-b088-605b0c261b0b
begin
	import HybridDynamics as HD
	import Plots as plt
	using LaTeXStrings
	using PlutoUI
end

# ╔═╡ 1934f0b1-2e7f-43ca-a0b8-885d49d00a53
md"""
# The bouncing ball
The purpose of this notebook is to offer a tutorial of HybridDynamics.jl through the bouncing ball. This is the hybrid system given by
```math
	\begin{cases}
		\ddot{x} = -g, & x > 0, \\
		\dot{x}^+ = -e\dot{x}, & x = 0.
	\end{cases}
```
"""

# ╔═╡ 9619965e-6551-48c1-900b-9f4e4716d79c
md"""
## A general hybrid system
!!! info "General Systems"
	A general system consists of the data $(f, h, Δ; \delta)$ where
	```math
		\begin{cases}
			\dot{x} = f(x, t), & h(x)\ne 0, \\
			x^+ = \Delta(x), & h(x) = 0.
		\end{cases}
	```
	The optional argument of $\delta \in \{-1, 0, 1\}$ dictates the crossing direction.
"""

# ╔═╡ 4445d15b-0813-4ab0-9184-c1809219ad08
@bind e Slider(0.0:0.01:1.0, show_value=true, default=0.8)

# ╔═╡ 514aa411-0315-4737-b253-75c1559b4392
# The data for the general problem formulation
begin
	# The dynamics
	function f_ball(x, t)
		g = 9.81
		q, v = x
		return [v, -g]
	end
	# The impact condition
	h_ball(x) = x[1]
	# The reset map
	Δ_ball(x) = [x[1], -e*x[2]]
end

# ╔═╡ b8befa98-6049-424b-8f4c-e7d868dea35d
# Set up and solve the problem
begin
	sysG = HD.GeneralSystem(f_ball, h_ball, Δ_ball; direction=-1)
	probG = HD.prob(sysG, [10.0, 0.0], (0.0, 15.0))
	solG = HD.solve(probG, HD.RK4())
end

# ╔═╡ 9f2d777b-8e8c-4ace-865e-4ec7c6226ed3
begin
	times, states = HD.split_jumps(solG)
	plt.plot(times, getindex.(states, 1), label="Position", lw=2, lc=:blue)
	plt.plot!(times, getindex.(states, 2), label="Velocity", lw=2, lc=:orange)
	plt.plot!(title = "Bouncing ball", xlabel = "Time", grid = true)
end

# ╔═╡ ec50258b-cd93-4064-9f2a-cb371d7d10f8
md"""
Notice that the solution is Zeno. Altough the problem had a specified final time of 15, the simulation stopped prematurely. 

By exploiding the mechanical nature of this problem, we can extend past the Zeno point.
"""

# ╔═╡ 6886ba9c-28c9-4955-860c-582b4541b732
md"""
## A mechanical hybrid system
!!! info "Mechanical Systems"
	A Mechanical system contains the data $(M, V, G, N, R, E)$
	1. .$M(q)$ is the mass matrix.
	2. .$V(q)$ is the potential energy.
	3. .$G(q)$ is the event location function.
	4. .$R(q) = dG(q)$ is the normal direction.
	5. .$E$ is the coefficient of restitution.
"""

# ╔═╡ 13da5896-c5ef-4829-b934-917a42495d01
begin
	# The mass matrix
	M(q) = 1.0
	# The potential energy
	V(q) = 9.81*q[1]
	# The guard location and its differential (as a vector)
	h(q) = q[1]
	∇h(q) = [1.0]
end

# ╔═╡ 88bd4e41-04e8-4fa2-9223-66b25b5bd6b0
begin
	sysM = HD.MechanicalSystem(M, V; guard=h, normal=∇h, e=e)
	probM = HD.prob(sysM, [10.0, 0.0], (0.0, 15.0))
	solM = HD.solve(probM, solver=HD.RK45())
end

# ╔═╡ 9dd93417-455f-4be9-a9de-f00df58845e3
begin
	tm, sm = HD.split_jumps(solM)
	plt.plot(tm, getindex.(sm, 1), label="Position", lw=2, lc=:blue)
	plt.plot!(tm, getindex.(sm, 2), label="Velocity", lw=2, lc=:orange)
	plt.plot!(title = "Bouncing ball", xlabel = "Time", grid = true)
end

# ╔═╡ 526471c6-3624-4633-9caf-b2ce6a4c9903
tz = solM.zeno[1]

# ╔═╡ 6db93a1d-8da0-4232-aee9-ec3d5c77d1ef
md"""
We can even determine the first Zeno time. This happens at $t_z =$ $tz.
"""

# ╔═╡ 20cc2ebd-76c7-468c-8133-e271a958fa10
md"""
## A stochastic hybrid system
!!! info "Stochastic Systems"
	A stochastic system contains the data $(f, g, h, Δ; δ)$
	```math
		\begin{cases}
			dx = f(x,t)dt + g(x,t)dW, & h(x) \ne 0, \\
			x^+ = \Delta(x), & h(x) = 0
		\end{cases}
	```
	Notice that $g$ must be given as a matrix.
"""

# ╔═╡ cb32d96e-c98f-467b-8dc6-f33b1c8748d6
g(x, t) = [0.0; 0.5;;]

# ╔═╡ a44b7da6-d928-4758-ae77-f6276d8d752f
begin
	sysST = HD.StochasticSystem(f_ball, g, h_ball, Δ_ball; direction=-1)
	probST = HD.prob(sysST, [10.0, 0.0], (0.0, 15.0))
	solST = HD.solve(probST; dt_initial=1e-3)
end

# ╔═╡ 1f0e7c83-86d4-4717-99dc-087c32329285
begin
	ts, ss = HD.split_jumps(solST)
	plt.plot(ts, getindex.(ss, 1), label="Position", lw=2, lc=:blue)
	plt.plot!(ts, getindex.(ss, 2), label="Velocity", lw=2, lc=:orange)
	plt.plot!(title = "Stochastic Bouncing ball", xlabel = "Time", grid = true)
end

# ╔═╡ 61cc9dff-a67e-4733-a8bd-84a2f5586bdc
md"""
Notice that stochastic systems do not contain the sliding mode logic that mechanical systems do.
"""

# ╔═╡ Cell order:
# ╠═6dab3570-7290-11f1-b67c-277aa2a67d52
# ╠═2b7869d2-dcd8-4ccf-b088-605b0c261b0b
# ╟─1934f0b1-2e7f-43ca-a0b8-885d49d00a53
# ╟─9619965e-6551-48c1-900b-9f4e4716d79c
# ╠═514aa411-0315-4737-b253-75c1559b4392
# ╠═4445d15b-0813-4ab0-9184-c1809219ad08
# ╠═b8befa98-6049-424b-8f4c-e7d868dea35d
# ╠═9f2d777b-8e8c-4ace-865e-4ec7c6226ed3
# ╟─ec50258b-cd93-4064-9f2a-cb371d7d10f8
# ╠═6886ba9c-28c9-4955-860c-582b4541b732
# ╠═13da5896-c5ef-4829-b934-917a42495d01
# ╠═88bd4e41-04e8-4fa2-9223-66b25b5bd6b0
# ╠═9dd93417-455f-4be9-a9de-f00df58845e3
# ╟─6db93a1d-8da0-4232-aee9-ec3d5c77d1ef
# ╠═526471c6-3624-4633-9caf-b2ce6a4c9903
# ╟─20cc2ebd-76c7-468c-8133-e271a958fa10
# ╠═cb32d96e-c98f-467b-8dc6-f33b1c8748d6
# ╠═a44b7da6-d928-4758-ae77-f6276d8d752f
# ╠═1f0e7c83-86d4-4717-99dc-087c32329285
# ╟─61cc9dff-a67e-4733-a8bd-84a2f5586bdc
