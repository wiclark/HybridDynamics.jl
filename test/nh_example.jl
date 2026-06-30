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

# ╔═╡ e1377a80-7493-11f1-ae92-dbdd66f67c8f
begin
	using Pkg
	Pkg.activate("..")
	Pkg.resolve()
	Pkg.instantiate()
	Pkg.precompile()
end

# ╔═╡ d48fccac-2d2d-41c7-8825-02aa9f2a5f2c
begin
	import HybridDynamics as HD
	import Plots as plt
	using LaTeXStrings
	using PlutoUI
end

# ╔═╡ cc04b485-6953-4f08-9b5c-06af7c27420c
using LinearAlgebra

# ╔═╡ cf77b769-2182-4377-8382-6107cee9c574
@bind a Slider(-0.5:0.01:0.5, show_value=true, default=0.0)

# ╔═╡ 879fc63c-303c-4e63-bf93-d60cc00ccb07
begin
	M_sleigh(q) = [1.0 0.0 -a*sin(q[3]);
			0.0 1.0 a*cos(q[3]);
			-a*sin(q[3]) a*cos(q[3]) 1+a^2]
	V_sleigh(q) = 0.0
	A_sleigh(q) = [-sin(q[3]) cos(q[3]) 0.0]
	h_sleigh(q) = 1 - (q[1]^2+q[2]^2)
	∇h_sleigh(q) = [-2*q[1], -2*q[2], 0.0]
end

# ╔═╡ fb350fee-7dd6-428f-aa5c-1214ea7984cb
begin
	sysNH = HD.NonholonomicSystem(M_sleigh, V_sleigh; A = A_sleigh, guard=h_sleigh, normal=∇h_sleigh, e=1)
	probNH = HD.prob(sysNH, [0.0, 0.0, 0.0, 1.0, 0.0, 1.5], (0.0, 4.3))
	solNH = HD.solve(probNH, solver=HD.RK4(), dt_intitial=1e-3)
end

# ╔═╡ ea6fc5b7-3677-410c-9eed-1a82f202a4dd
begin
	xs = getindex.(solNH.x, 1)
	ys = getindex.(solNH.x, 2)
	plt.plot(xs, ys, label="Sleigh Trajectory", lw=2)
	θ = LinRange(0, 2π, 100)
	plt.plot!(cos.(θ), sin.(θ), label="", lc=:black, aspect_ratio=1)
end

# ╔═╡ d9323fb7-db50-4a2b-aa22-b4086ab25e46
length(solNH.dx)

# ╔═╡ 0031dd50-c979-48e5-9980-a1d187b19dd1
solNH(1.0)

# ╔═╡ 13020d7d-c0de-4712-ab69-4def7c052495
length(solNH.x)

# ╔═╡ aeeb9cc0-98c8-4f0f-ba10-cb57c659945d
plt.plot(solNH.t, xs)

# ╔═╡ a6f01395-7454-4683-aa4e-b928143a9e81
H(t) = 1/2*dot(solNH.x(t)[4:6], M_sleigh(solNH.x(t)[1:3]) \ solNH.x(t)[4:6])

# ╔═╡ Cell order:
# ╠═e1377a80-7493-11f1-ae92-dbdd66f67c8f
# ╠═d48fccac-2d2d-41c7-8825-02aa9f2a5f2c
# ╠═cf77b769-2182-4377-8382-6107cee9c574
# ╠═879fc63c-303c-4e63-bf93-d60cc00ccb07
# ╠═fb350fee-7dd6-428f-aa5c-1214ea7984cb
# ╠═ea6fc5b7-3677-410c-9eed-1a82f202a4dd
# ╠═d9323fb7-db50-4a2b-aa22-b4086ab25e46
# ╠═0031dd50-c979-48e5-9980-a1d187b19dd1
# ╠═13020d7d-c0de-4712-ab69-4def7c052495
# ╠═aeeb9cc0-98c8-4f0f-ba10-cb57c659945d
# ╠═a6f01395-7454-4683-aa4e-b928143a9e81
# ╠═cc04b485-6953-4f08-9b5c-06af7c27420c
