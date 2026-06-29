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

# ╔═╡ 2c1a0a30-73cd-11f1-a35e-b14c626e34d6
begin
	using Pkg
	Pkg.activate("..")
	Pkg.resolve()
	Pkg.instantiate()
	Pkg.precompile()
end

# ╔═╡ dbb560aa-7fe7-4401-be28-f0f0d71e918f
begin
	import HybridDynamics as HD
	import Plots as plt
	using LaTeXStrings
	using PlutoUI
end

# ╔═╡ d2510d9a-1fda-461c-b7c2-13fb3244e00e
begin
	M(q) = [1.0 0.0; 0.0 1.0]
	V(q) = q[2]
	h(q) = 1 - (q[1]^2+q[2]^2)
	∇h(q) = [-2*q[1], -2*q[2]]
end

# ╔═╡ b46ddac0-9f57-4be6-aed2-3ef78944f6b6
@bind r Slider(0.0:0.01:1.0, show_value=true, default=0.07)

# ╔═╡ 4eefa8c3-f595-4878-ae46-6240ef512c32
begin
	sysM = HD.MechanicalSystem(M, V; guard=h, normal=∇h, e=r)
	probM = HD.prob(sysM, [-0.95, 0.0, 0.2, -1.5], (0.0, 10.0))
	solM = HD.solve(probM, solver=HD.RK45())
end

# ╔═╡ 9838cc5d-ee2d-47c2-b581-44bff7e9d4d5
solM.t

# ╔═╡ ac752ece-ef2d-4059-85cb-3dc887936ed5
solM.x

# ╔═╡ 5f31c7a7-4e78-4e49-b7ca-113c8a271dd3
length(solM.t)

# ╔═╡ c494b496-2058-4c0f-9446-90390ffb2e85
begin
	xm = getindex.(solM.x, 1)
	ym = getindex.(solM.x, 2)
	plt.plot(xm, ym, label="Ball Trajectory", lw=2)
	θ = LinRange(0, 2π, 100)
	plt.plot!(cos.(θ), sin.(θ), label="", lc=:black, aspect_ratio=1)
end

# ╔═╡ 34a81564-93ad-497b-b119-2ea307ea2a82
solM.zeno[end]

# ╔═╡ 64ea3f18-825d-4e24-ab7a-c4b452af2b1c
solM.t[end]

# ╔═╡ Cell order:
# ╠═2c1a0a30-73cd-11f1-a35e-b14c626e34d6
# ╠═dbb560aa-7fe7-4401-be28-f0f0d71e918f
# ╠═d2510d9a-1fda-461c-b7c2-13fb3244e00e
# ╠═4eefa8c3-f595-4878-ae46-6240ef512c32
# ╠═b46ddac0-9f57-4be6-aed2-3ef78944f6b6
# ╠═9838cc5d-ee2d-47c2-b581-44bff7e9d4d5
# ╠═ac752ece-ef2d-4059-85cb-3dc887936ed5
# ╠═5f31c7a7-4e78-4e49-b7ca-113c8a271dd3
# ╠═c494b496-2058-4c0f-9446-90390ffb2e85
# ╠═34a81564-93ad-497b-b119-2ea307ea2a82
# ╠═64ea3f18-825d-4e24-ab7a-c4b452af2b1c
