### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 6df38dc0-6787-11f1-a77e-9f1ff29a6e5b
begin
	using Pkg
	Pkg.activate("..")
	Pkg.resolve()
	Pkg.instantiate()
	Pkg.precompile()
end

# ╔═╡ 2b9b16b3-91fd-4b17-bee2-863c0e5ef544
using LaTeXStrings

# ╔═╡ e2f23c5b-afad-4d71-b3f9-3af04448d873
import HybridDynamics as HD

# ╔═╡ 70f106e0-7137-4fb4-9b20-371adb3efaa7
import Plots as plt

# ╔═╡ 30a03940-11a2-45e8-8f91-501370ddcee6
md"""
## A Filippov Example
```math
	\begin{cases}
		\dot{x} = 3, \ \dot{y} = -1, & y > \sin x \\[1ex]
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

# ╔═╡ 9da11b33-0def-4ac1-bd19-574ef1cceb98
begin
	F(x) = [3, -1]
	G(x) = [0, 1]
	H(x) = x[2] - sin(x[1])
	N(x) = [-cos(x[1]), 1]
end

# ╔═╡ 0cf59ef8-d158-42dc-913a-7766ac4d2420
begin
	sysF = HD.FilippovSys(F, G, H, N)
	probF = HD.prob(sysF, [0.0, 1.0], (0.0, 10.0))
	solF = HD.solve(probF, HD.RK4(); dt_initial=0.01)
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
	plt.plot!(title = "Filippov Trajectory", label=L"x", ylabel=L"y", dpi=500)
end

# ╔═╡ 7b9a3bb2-0b4f-48bf-9c7a-1824281ad707
md"""
## A Mechancal Example
"""

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
Δ_ball(x) = [abs(x[1]), 0.8*abs(x[2])]

# ╔═╡ cb05396f-92a2-468c-af34-f64e25a93e16
begin
	sysG = HD.GeneralSystem(f_ball, h_ball, Δ_ball)
	probG = HD.prob(sysG, [10.0, 0.0], (0.0, 5.0))
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

# ╔═╡ Cell order:
# ╠═6df38dc0-6787-11f1-a77e-9f1ff29a6e5b
# ╠═e2f23c5b-afad-4d71-b3f9-3af04448d873
# ╠═70f106e0-7137-4fb4-9b20-371adb3efaa7
# ╠═2b9b16b3-91fd-4b17-bee2-863c0e5ef544
# ╟─30a03940-11a2-45e8-8f91-501370ddcee6
# ╟─2880f5f2-fb01-4696-aba3-7faf4714c9d3
# ╠═9da11b33-0def-4ac1-bd19-574ef1cceb98
# ╠═0cf59ef8-d158-42dc-913a-7766ac4d2420
# ╠═1959511b-8aa8-49d0-893c-b9b801a2bda0
# ╠═bbc202e7-fd8d-4eb1-9de2-93b09a73a425
# ╠═7b9a3bb2-0b4f-48bf-9c7a-1824281ad707
# ╟─eea32534-7258-478d-b96e-aee9bc212011
# ╟─15b10a81-2909-4df8-8d1c-0fcd7af5d9ff
# ╠═33b07ebd-6275-4aa5-a855-7a6b29b97de9
# ╠═6e9ae3eb-e2d2-42b2-bcfe-c42c6dd4be65
# ╠═8702895d-c0ea-459b-814a-b67f4c193860
# ╠═cb05396f-92a2-468c-af34-f64e25a93e16
# ╟─72a895f8-655e-4a0e-8025-90458f0f2819
# ╠═5bb10f03-ebc5-4b86-a150-04c0eb12e404
