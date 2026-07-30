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

# ╔═╡ a0e0bdeb-f73b-4917-b054-f35b73880d9d
begin
	import HybridDynamics as HD
	import Plots as plt
	using LaTeXStrings
	using PlutoUI
end

# ╔═╡ d1c084d0-809a-11f1-9ac4-a5a2296001e4
md"""
## Life after Zeno examples

In the 2006 paper by Ames et. al., "Is there life after Zeno? Taking executions past the breaking (Zeno) point", three examples are used to demonstrate a method of continuing the solutions of mechanical hybrid systems past Zeno points. Those three examples are presented here for demonstration.

---
Reference:

Ames, Aaron & Zheng, Haiyang & Gregg, Robert & Sastry, Shankar. (2006). Is there life after Zeno? Taking executions past the breaking (Zeno) point. 2006. 6 pp.. 10.1109/ACC.2006.1656623
"""

# ╔═╡ 4b4b0f78-8c57-48a2-b5ad-6986b343138c
md"""
### Bouncing ball on a sinusoidal surface
"""

# ╔═╡ 14e8698f-5d20-4218-b2e3-6d8b5e2c4a46
md"""
In two dimensions, it's easier to see what's going on. Consider the Hamiltonian system,
```math
	H = p_x^2 + p_y^2 + y,
```
where impacts occur when $y = \sin(x)$. Setting this up as a mechanical system requires the mass matrix, the potential energy, and the event function.
"""

# ╔═╡ 9c66f98c-8f51-4c12-8ae9-e2cb562654f0
begin
	M1a(q) = [1.0 0.0; 0.0 1.0]
	V1a(q) = q[2]
	h1a(q) = q[2] - sin(q[1])
end;

# ╔═╡ f447e6eb-7826-4b41-9c05-ff554ceb0ddd
md"""
We set up the problem. The initial conditions are in momentum space and have the form
```math
	z_0 = [q_0, p_0].
```
The simulation will run over the interval $0\leq t\leq 10$ and be solved via the fixed stepper RK4.
"""

# ╔═╡ 36635588-2a13-434b-b613-767a17e9437c
Coefficient_of_restitution_1a = @bind r1a Slider(0.0:0.01:1.0, show_value=true)

# ╔═╡ 2d00bb46-6df2-4087-ad97-dfa2700524a7
begin
	init1a = [0.0, 2.0, 2.3, -1.5]
	sys1a = HD.MechanicalSystem(M1a, V1a; guard=h1a, e=r1a)
	prob1a = HD.prob(sys1a, init1a, (0.0, 10.0))
	sol1a = HD.solve(prob1a, HD.RK4())
end

# ╔═╡ ea3d9084-6b36-49c7-ab55-5628d0c0ba33
begin
	plt.gr()
	xm1a = getindex.(sol1a.x, 1)
	ym1a = getindex.(sol1a.x, 2)
	
	plt.plot(xm1a, ym1a, lw=2,
			 label="Ball Trajectory", xlabel = L"x", ylabel = L"y",
			 aspect_ratio=1, dpi = 300)

	plt.scatter!([init1a[1]], [init1a[2]], label = "Initial position")
	
	θ1a = LinRange(minimum(xm1a)-2, maximum(xm1a)+2, 100)
	plt.plot!(θ1a, sin.(θ1a), label="", lc=:black)
end

# ╔═╡ f31a4457-1d1a-42dc-b329-8793e042ac09
md"""
The previous example can be extended to three dimensions by:
```math
	H = p_x^2 + p_y^2 + p_z^2 + z,
```
where the event function is now $z = 0.5\sin(y)$.
"""

# ╔═╡ 932a2a6a-7391-4fb0-b6c7-dda85db5ef37
begin
	M1b(q) = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] #I'm guessing it's this
	V1b(q) = q[3]
	h1b(q) = q[3] - 0.5*sin(q[2])
	∇h1b(q) = [0., - 0.5*cos(q[2]), 1.]
end;

# ╔═╡ 679c859e-3be9-4774-aa7f-3eeff2301e52
md"""
Initializing a MechanicalSystem can also take the normal of the event function as an optional argument.
"""

# ╔═╡ b8972101-5a59-4f3c-98fb-fa2c4ec29913
Coefficient_of_restitution_1b = @bind r1b Slider(0.0:0.01:1.0, show_value=true)

# ╔═╡ ca9dcc3a-c73e-413d-a78b-24a4005f3434
begin
	init1b = [0.0, 2.0, 1,
			  0.3, -1.5, 1]
	sys1b = HD.MechanicalSystem(M1b, V1b; guard=h1b, normal = ∇h1b, e=r1b)
	prob1b = HD.prob(sys1b, init1b, (0.0, 10.0))
	sol1b = HD.solve(prob1b, HD.RK4())
end

# ╔═╡ a6cdc0c7-5ab8-4f78-a9c0-a0c673b5c44a
begin
	# Plot the trajectory on the sinusoidal surface
	plt.plotly()
	# Extract out the three position components
	xm1b = getindex.(sol1b.x, 1)
	ym1b = getindex.(sol1b.x, 2)
	zm1b = getindex.(sol1b.x, 3)
	
	plt.plot(xm1b, ym1b, zm1b, lw=2,
			 label="Ball Trajectory", xlabel = "x", ylabel = "y", zlabel = "z")
	# Plot the surface
	xg1b = LinRange(minimum(xm1b)-5, maximum(xm1b)+5, 50)
    yg1b = LinRange(minimum(ym1b)-5, maximum(ym1b)+5, 100)
    zg1b = [0.5*sin(y) for y in yg1b, x in xg1b]

    plt.surface!(xg1b, yg1b, zg1b,
		c =:black, colorbar = false, alpha = 0.8,
        label = "Guard"
    )

end

# ╔═╡ 52aed06f-308f-44fa-9bfb-d32182373a2a
md"""
This trajectory becomes Zeno at time $t_z$= $(sol1b.zeno[1]).
"""

# ╔═╡ 7f9ade56-583a-4849-befc-053103754cf2
sol1b.zeno[1]

# ╔═╡ ac42fc4f-d2a2-4626-a032-6eca1e6a0f50
begin
	plt.gr()
	hvals = [h1b(x[1:3]) for x in sol1b.x]

	plt.plot(sol1b.t, hvals, label=L"$h(x(t))$", xlabel=L"$t$", ylabel="Height")
	plt.hline!([0], label="Surface")
	plt.plot!(title = "Distance from event surface over time")
end

# ╔═╡ Cell order:
# ╟─d1c084d0-809a-11f1-9ac4-a5a2296001e4
# ╠═a0e0bdeb-f73b-4917-b054-f35b73880d9d
# ╟─4b4b0f78-8c57-48a2-b5ad-6986b343138c
# ╟─14e8698f-5d20-4218-b2e3-6d8b5e2c4a46
# ╠═9c66f98c-8f51-4c12-8ae9-e2cb562654f0
# ╟─f447e6eb-7826-4b41-9c05-ff554ceb0ddd
# ╠═2d00bb46-6df2-4087-ad97-dfa2700524a7
# ╟─36635588-2a13-434b-b613-767a17e9437c
# ╟─ea3d9084-6b36-49c7-ab55-5628d0c0ba33
# ╟─f31a4457-1d1a-42dc-b329-8793e042ac09
# ╠═932a2a6a-7391-4fb0-b6c7-dda85db5ef37
# ╟─679c859e-3be9-4774-aa7f-3eeff2301e52
# ╠═ca9dcc3a-c73e-413d-a78b-24a4005f3434
# ╟─b8972101-5a59-4f3c-98fb-fa2c4ec29913
# ╠═a6cdc0c7-5ab8-4f78-a9c0-a0c673b5c44a
# ╟─52aed06f-308f-44fa-9bfb-d32182373a2a
# ╠═7f9ade56-583a-4849-befc-053103754cf2
# ╠═ac42fc4f-d2a2-4626-a032-6eca1e6a0f50
