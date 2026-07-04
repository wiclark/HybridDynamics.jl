### A Pluto.jl notebook ###
# v0.20.24

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

# ╔═╡ 6df38dc0-6787-11f1-a77e-9f1ff29a6e5b
begin
	using Pkg
	Pkg.activate("..")
	Pkg.resolve()
	Pkg.instantiate()
	Pkg.precompile()
end

# ╔═╡ 5f758439-587b-4bbe-bd02-f2f9ce6cc8fb
begin
	import HybridDynamics as HD
	import Plots as plt
	using LaTeXStrings
	using PlutoUI
end

# ╔═╡ 30a03940-11a2-45e8-8f91-501370ddcee6
md"""
## A Filippov Example
```math
	\begin{cases}
		\dot{x} = 3, \ \dot{y} = ξ, & y > \sin x \\[1ex]
		\dot{x} = 0, \ \dot{y} = 1, & y < \sin x
	\end{cases}
```
"""

# ╔═╡ 2880f5f2-fb01-4696-aba3-7faf4714c9d3
md"""
!!! info "A Filippov Problem"
	A Filippov system contains the data $(F,G,H,N)$ where
	```math
		\begin{cases}
			\dot{x} = F(x), & H(x) > 0, \\
			\dot{x} = G(x), & H(x) < 0,
		\end{cases}
	```
	and $N = \nabla H$.
"""

# ╔═╡ e6bc3aba-7b99-4ca8-8657-12b8609cd427
@bind ξ Slider(-2:0.1:0, show_value=true)

# ╔═╡ 9da11b33-0def-4ac1-bd19-574ef1cceb98
begin
	F(x) = [3, ξ]
	G(x) = [0, 1]
	H(x) = x[2] - sin(x[1])
	N(x) = [-cos(x[1]), 1]
end

# ╔═╡ 0cf59ef8-d158-42dc-913a-7766ac4d2420
begin
	sysF = HD.FilippovSystem(F, G, H, N)
	probF = HD.prob(sysF, [0.0, 1.0], (0.0, 10.0))
	solF = HD.solve(probF, HD.RK4())
end

# ╔═╡ 1959511b-8aa8-49d0-893c-b9b801a2bda0
begin
	xf = getindex.(solF.x, 1)
	yf = getindex.(solF.x, 2)
	xh = range(minimum(xf)-0.5, maximum(xf)+0.5, length=1_000)
end

# ╔═╡ bbc202e7-fd8d-4eb1-9de2-93b09a73a425
begin
	plt.plot(xf, yf, lw=2, label="Filippov Trajectory")
	plt.plot!(xh, sin.(xh), lw=2, lc=:black, ls=:dash, label=L"H(x)=0")
	plt.plot!(title = "Filippov Trajectory", label=L"x", ylabel=L"y")
end

# ╔═╡ 7b9a3bb2-0b4f-48bf-9c7a-1824281ad707
md"""
## A Mechancal Example
```math
	H(x,y,p_x,p_y) = \frac{1}{2}(p_x^2+p_y^2) + y
```
```math
	h(x,y) = 1 - (x^2+y^2) \geq 0.
```
"""

# ╔═╡ 80034ca8-e254-4beb-8486-da9de404f503
md"""
!!! info "A Mechanical Problem"
	A Mechanical system contains the data $(M, V, G, N, R, E)$
	1. .$M(q)$ is the mass matrix.
	2. .$V(q)$ is the potential energy.
	3. .$G(q)$ is the event location function.
	4. .$R(q) = dG(q)$ is the normal direction.
	5. .$E$ is the coefficient of restitution.
"""

# ╔═╡ 49b430c8-ac20-4988-9c49-83a38111285f
begin
	M(q) = [1.0 0.0; 0.0 1.0]
	V(q) = q[2]
	h(q) = 1 - (q[1]^2+q[2]^2)
	∇h(q) = [-2*q[1], -2*q[2]]
end

# ╔═╡ a5cd41cd-b341-49c7-a052-bfd20f631b0c
@bind r Slider(0.0:0.01:1.0, show_value=true)

# ╔═╡ c1a00790-005e-40ae-9e4c-20f8aac7b787
begin
	sysM = HD.MechanicalSystem(M, V; guard=h, normal=∇h, e=r)
	probM = HD.prob(sysM, [-0.95, 0.0, 0.2, -1.5], (0.0, 10.0))
	solM = HD.solve(probM, solver=HD.RK45())
end

# ╔═╡ fb3eb7c9-3f93-4986-aa89-c2ca478429d3
begin
	xm = getindex.(solM.x, 1)
	ym = getindex.(solM.x, 2)
	plt.plot(xm, ym, label="Ball Trajectory", lw=2)
	θ = LinRange(0, 2π, 100)
	plt.plot!(cos.(θ), sin.(θ), label="", lc=:black, aspect_ratio=1)
end

# ╔═╡ 597d4dce-acf4-44a1-bc7a-9d7ec35de2f7
md"""
## A Nonholonomic Example
#### The Chaplygin sleigh
```math
	M(q) = \begin{bmatrix}
		m & 0 & -ma\sin\theta \\
		0 & m & ma\cos\theta \\
		-ma\sin\theta & ma\cos\theta & I+ma^2
	\end{bmatrix}
```
```math
	A(q) = \begin{bmatrix}
		-\sin\theta & \cos\theta & 0
	\end{bmatrix}
```
```math
	h(x,y) = 1 - (x^2+y^2) \geq 0.
```
"""

# ╔═╡ 8e0f34ac-35ab-4808-a1ed-53525361c688
@bind a Slider(-0.5:0.01:0.5, show_value=true, default=0.0)

# ╔═╡ 6d47f9d2-c604-46ba-af77-254478d19ab0
begin
	M_sleigh(q) = [1.0 0.0 -a*sin(q[3]);
			0.0 1.0 a*cos(q[3]);
			-a*sin(q[3]) a*cos(q[3]) 1+a^2]
	V_sleigh(q) = 0.0
	A_sleigh(q) = [-sin(q[3]) cos(q[3]) 0.0]
	h_sleigh(q) = 1 - (q[1]^2+q[2]^2)
	∇h_sleigh(q) = [-2*q[1], -2*q[2], 0.0]
end

# ╔═╡ da267e54-7c2b-4a31-b429-077b47912814
begin
	sysNH = HD.NonholonomicSystem(M_sleigh, V_sleigh; A = A_sleigh, guard=h_sleigh, normal=∇h_sleigh, e=1)
	probNH = HD.prob(sysNH, [0.0, 0.0, 0.0, 1.0, 0.0, 1.5], (0.0, 100.0))
	solNH = HD.solve(probNH, solver=HD.RK4())
end

# ╔═╡ d8370b2e-077d-488f-afb2-ceffa00fc255
begin
	xs = getindex.(solNH.x, 1)
	ys = getindex.(solNH.x, 2)
	plt.plot(xs, ys, label="Sleigh Trajectory", lw=2)
	plt.plot!(cos.(θ), sin.(θ), label="", lc=:black, aspect_ratio=1)
end

# ╔═╡ eea32534-7258-478d-b96e-aee9bc212011
md"""
## A General Example
```math
	\begin{cases}
		\ddot{x} = -g -\alpha\cdot \dot{x} \cdot|\dot{x}|, & x > 0 \\[1ex]
		\dot{x}^+ = -e\dot{x}^-, & x = 0
	\end{cases}
```
"""

# ╔═╡ 15b10a81-2909-4df8-8d1c-0fcd7af5d9ff
md"""
!!! info "A General Problem"
	A general system consists of the data $(f, h, \Delta)$ where
	```math
		\begin{cases}
			\dot{x} = f(x, t), & h(x) \ne 0, \\
			x^+ = \Delta(x), & h(x) = 0.
		\end{cases}
	```
"""

# ╔═╡ 33b07ebd-6275-4aa5-a855-7a6b29b97de9
function f_ball(x, t)
	g = 9.81
	α = 0.1
	q, v = x
	return [v, -g-α*v*abs(v)]
end

# ╔═╡ 6e9ae3eb-e2d2-42b2-bcfe-c42c6dd4be65
h_ball(x) = x[1]

# ╔═╡ 8702895d-c0ea-459b-814a-b67f4c193860
Δ_ball(x) = [abs(x[1]), -0.8*x[2]]

# ╔═╡ cb05396f-92a2-468c-af34-f64e25a93e16
begin
	sysG = HD.GeneralSystem(f_ball, h_ball, Δ_ball; direction=-1)
	probG = HD.prob(sysG, [10.0, 0.0], (0.0, 15.0))
	solG = HD.solve(probG)
end

# ╔═╡ 7ca75d5b-b8a6-472f-b5a7-f746b581e67f
begin
	times, states = HD.split_jumps(solG)
	plt.plot(times, getindex.(states, 1), label="Position", lw=2, lc=:blue)
	plt.plot!(times, getindex.(states, 2), label="Velocity", lw=2, lc=:orange)
	plt.plot!(title = "Bouncing ball", xlabel = "Time", grid = true)
end

# ╔═╡ 293d4ac4-bbb7-41ea-ae29-ae0844fd3992
solG.jump_times

# ╔═╡ efe951a9-31d1-49e1-bf19-cdcd7a6a5f0e
md"""
## A Stochastic Example
```math
	\begin{cases}
		dx = -(x-3)dt + 0.2 dW, & x < 2 \\
		x^+ = x-1, & x = 2
	\end{cases}
```
"""

# ╔═╡ 8eae9c66-e573-4458-bab5-ac92e61f261a
md"""
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

# ╔═╡ 7e2f7386-647e-4503-9343-2ce78cc2752f
begin
	f_st(x, t) = -(x .-3.)
	g_st(x, t) = [0.2;;]
	h_st(x) = x[1]-2
	Δ_st(x) = x .- 1.0
end

# ╔═╡ 0a709e7a-f8a4-4b72-aa15-eb3fb92f517e
begin
	sysST = HD.StochasticSystem(f_st, g_st, h_st, Δ_st)
	probST = HD.prob(sysST, [0.5], (0.0, 3.0))
	solST = HD.solve(probST)
end

# ╔═╡ b31969a0-fac6-4e69-a9f1-391750214770
begin
	ts, ss = HD.split_jumps(solST)
	plt.plot(ts, getindex.(ss, 1), label="", lw=2, lc=:blue)
	plt.plot!(title = "Stochastic Bouncing ball", xlabel = "Time", grid = true)
end

# ╔═╡ 8e899a64-9654-4e7e-ac0f-2376b89aa913
md"""
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
"""

# ╔═╡ 154b65b0-f68b-4ec7-a035-b5d3ea78fb26
begin
	LA_type = "Linear" # Options: "Linear", "Affine"
	LA_A = [0.0 -10.0; 10.0 -.1]
	LA_λ = [1.0, 0.0]
	LA_C = [1.5 1.0; 0.0 0.5]
	LA_x0 = [1.0, 0.5]
	LA_tspan = (0.0, 10.0)
end

# ╔═╡ 7259e7cc-e810-41c7-b57b-acedd5ed6691
begin
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
        LA_sol = HD.solve(LA_prob, HD.AdamsBashforth3())
    end
end

# ╔═╡ 1e618fc0-28ab-4c28-9be5-8b1582a69cba
LA_t_list, LA_x_list = HD.split_jumps(LA_sol)

# ╔═╡ a77c28dc-8d08-42b6-8e51-1949ea04bb4a
begin
    p = plt.plot(title="$(LA_type) System Trajectory", xlabel="x1", ylabel="x2", 
                 aspect_ratio=:equal, grid=true, legend=false)
    
    plt.vline!(p, [LA_guard], label="Guard", color=:black, alpha=0.3, lw=2)
    
    plt.plot!(p, getindex.(LA_x_list, 1), getindex.(LA_x_list, 2), 
              color=:red, linestyle=:dash, lw=1.5, label="Trajectory")
    p
end

# ╔═╡ Cell order:
# ╠═6df38dc0-6787-11f1-a77e-9f1ff29a6e5b
# ╠═5f758439-587b-4bbe-bd02-f2f9ce6cc8fb
# ╟─30a03940-11a2-45e8-8f91-501370ddcee6
# ╟─2880f5f2-fb01-4696-aba3-7faf4714c9d3
# ╠═9da11b33-0def-4ac1-bd19-574ef1cceb98
# ╠═0cf59ef8-d158-42dc-913a-7766ac4d2420
# ╠═1959511b-8aa8-49d0-893c-b9b801a2bda0
# ╠═e6bc3aba-7b99-4ca8-8657-12b8609cd427
# ╠═bbc202e7-fd8d-4eb1-9de2-93b09a73a425
# ╟─7b9a3bb2-0b4f-48bf-9c7a-1824281ad707
# ╠═80034ca8-e254-4beb-8486-da9de404f503
# ╠═49b430c8-ac20-4988-9c49-83a38111285f
# ╠═c1a00790-005e-40ae-9e4c-20f8aac7b787
# ╠═a5cd41cd-b341-49c7-a052-bfd20f631b0c
# ╠═fb3eb7c9-3f93-4986-aa89-c2ca478429d3
# ╟─597d4dce-acf4-44a1-bc7a-9d7ec35de2f7
# ╠═8e0f34ac-35ab-4808-a1ed-53525361c688
# ╠═6d47f9d2-c604-46ba-af77-254478d19ab0
# ╠═da267e54-7c2b-4a31-b429-077b47912814
# ╠═d8370b2e-077d-488f-afb2-ceffa00fc255
# ╠═eea32534-7258-478d-b96e-aee9bc212011
# ╟─15b10a81-2909-4df8-8d1c-0fcd7af5d9ff
# ╠═33b07ebd-6275-4aa5-a855-7a6b29b97de9
# ╠═6e9ae3eb-e2d2-42b2-bcfe-c42c6dd4be65
# ╠═8702895d-c0ea-459b-814a-b67f4c193860
# ╠═cb05396f-92a2-468c-af34-f64e25a93e16
# ╠═7ca75d5b-b8a6-472f-b5a7-f746b581e67f
# ╠═293d4ac4-bbb7-41ea-ae29-ae0844fd3992
# ╠═a62b426b-cb7f-44f1-9371-9f446a494a59
# ╠═67463097-3bb9-48a6-bb63-008eb7936427
# ╟─8e899a64-9654-4e7e-ac0f-2376b89aa913
# ╠═154b65b0-f68b-4ec7-a035-b5d3ea78fb26
# ╠═7259e7cc-e810-41c7-b57b-acedd5ed6691
# ╠═1e618fc0-28ab-4c28-9be5-8b1582a69cba
# ╠═a77c28dc-8d08-42b6-8e51-1949ea04bb4a