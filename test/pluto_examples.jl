### A Pluto.jl notebook ###
# v1.0.1

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

# ╔═╡ e405e5fb-703e-452a-a704-962a162f9751
plt.plot(solF.s)

# ╔═╡ 1959511b-8aa8-49d0-893c-b9b801a2bda0
begin
	xf = getindex.(solF.x, 1)
	yf = getindex.(solF.x, 2)
	xh = range(minimum(xf)-0.5, maximum(xf)+0.5, length=1_000)
end

# ╔═╡ 7198ae6b-a1d8-4525-94f3-0c6c37c2346b
plt.plot(H.(solF.x))

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

# ╔═╡ 999c672e-fdff-4225-b6c8-80a7096432c8
solM(1.0)

# ╔═╡ fb3eb7c9-3f93-4986-aa89-c2ca478429d3
begin
	xm = getindex.(solM.x, 1)
	ym = getindex.(solM.x, 2)
	plt.plot(xm, ym, label="Ball Trajectory", lw=2)
	θ = LinRange(0, 2π, 100)
	plt.plot!(cos.(θ), sin.(θ), label="", lc=:black, aspect_ratio=1)
end

# ╔═╡ 5a82346d-c04c-4438-85b8-f6d66e8f32db
plt.plot(solM.t, 1 .- (xm.^2 + ym.^2))

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

# ╔═╡ e4b47eb0-38f7-490f-ae6d-8f6d0ad278d8
size([-0.0 1.0 0.0])[1]

# ╔═╡ c181233d-d767-4cc6-a114-1f6ee16a2fec
sysNH

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

# ╔═╡ 72a895f8-655e-4a0e-8025-90458f0f2819
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

# ╔═╡ 5bb10f03-ebc5-4b86-a150-04c0eb12e404
begin
	t_plot, x1_plot, x2_plot = get_rid_of_jump_lines(solG)
	plt.plot(t_plot, x1_plot, label="Position", lw=2, lc=:blue)
	plt.plot!(t_plot, x2_plot, label="Velocity", lw=2, color=:orange)
	plt.plot!(title = "Nonlinear bouncing ball", xlabel = "Time", grid = true)
end

# ╔═╡ 293d4ac4-bbb7-41ea-ae29-ae0844fd3992
solG.jump_times

# ╔═╡ a62b426b-cb7f-44f1-9371-9f446a494a59
solG.x

# ╔═╡ 67463097-3bb9-48a6-bb63-008eb7936427
solG.t

# ╔═╡ Cell order:
# ╠═6df38dc0-6787-11f1-a77e-9f1ff29a6e5b
# ╠═5f758439-587b-4bbe-bd02-f2f9ce6cc8fb
# ╟─30a03940-11a2-45e8-8f91-501370ddcee6
# ╟─2880f5f2-fb01-4696-aba3-7faf4714c9d3
# ╠═9da11b33-0def-4ac1-bd19-574ef1cceb98
# ╠═0cf59ef8-d158-42dc-913a-7766ac4d2420
# ╠═7198ae6b-a1d8-4525-94f3-0c6c37c2346b
# ╠═e405e5fb-703e-452a-a704-962a162f9751
# ╠═1959511b-8aa8-49d0-893c-b9b801a2bda0
# ╠═e6bc3aba-7b99-4ca8-8657-12b8609cd427
# ╠═bbc202e7-fd8d-4eb1-9de2-93b09a73a425
# ╟─7b9a3bb2-0b4f-48bf-9c7a-1824281ad707
# ╟─80034ca8-e254-4beb-8486-da9de404f503
# ╠═49b430c8-ac20-4988-9c49-83a38111285f
# ╠═c1a00790-005e-40ae-9e4c-20f8aac7b787
# ╠═999c672e-fdff-4225-b6c8-80a7096432c8
# ╠═a5cd41cd-b341-49c7-a052-bfd20f631b0c
# ╠═fb3eb7c9-3f93-4986-aa89-c2ca478429d3
# ╠═5a82346d-c04c-4438-85b8-f6d66e8f32db
# ╟─597d4dce-acf4-44a1-bc7a-9d7ec35de2f7
# ╠═8e0f34ac-35ab-4808-a1ed-53525361c688
# ╠═6d47f9d2-c604-46ba-af77-254478d19ab0
# ╠═da267e54-7c2b-4a31-b429-077b47912814
# ╠═e4b47eb0-38f7-490f-ae6d-8f6d0ad278d8
# ╠═c181233d-d767-4cc6-a114-1f6ee16a2fec
# ╠═d8370b2e-077d-488f-afb2-ceffa00fc255
# ╟─eea32534-7258-478d-b96e-aee9bc212011
# ╟─15b10a81-2909-4df8-8d1c-0fcd7af5d9ff
# ╠═33b07ebd-6275-4aa5-a855-7a6b29b97de9
# ╠═6e9ae3eb-e2d2-42b2-bcfe-c42c6dd4be65
# ╠═8702895d-c0ea-459b-814a-b67f4c193860
# ╠═cb05396f-92a2-468c-af34-f64e25a93e16
# ╟─72a895f8-655e-4a0e-8025-90458f0f2819
# ╠═5bb10f03-ebc5-4b86-a150-04c0eb12e404
# ╠═293d4ac4-bbb7-41ea-ae29-ae0844fd3992
# ╠═a62b426b-cb7f-44f1-9371-9f446a494a59
# ╠═67463097-3bb9-48a6-bb63-008eb7936427
